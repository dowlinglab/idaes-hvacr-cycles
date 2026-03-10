"""
PLR Vapor Compression Cycle with IDAES HeatExchangerNTU Units

Author: Shilpa Narasimhan
Codex Support: OpenAI Codex
QA/Testing: Shilpa Narasimhan

Description:
    Defines a copy-only IDAES vapor-compression cycle model that preserves the
    existing PLR post-correction workflow while replacing heater-style
    evaporator and condenser units with HeatExchangerNTU units.

Context Breadcrumb:
    This is an additive copy that keeps PLR structure while replacing heater
    coils with IDAES HeatExchangerNTU units for evaporator and condenser.

Notes:
    - All QA and validation are the responsibility of Shilpa Narasimhan.
    - This file was generated with Codex assistance.
"""

import logging
from enum import Enum

import CoolProp.CoolProp as CP
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Objective,
    Param,
    Var,
    maximize,
    value,
    units as pyunits,
    exp,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import calculate_scaling_factors
from pyomo.opt import TerminationCondition
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    AmountBasis,
    StateVars,
)
from idaes.models.unit_models import Compressor, HeatExchangerNTU, PressureChanger
from idaes.models_extra.power_generation.properties.flue_gas_ideal import (
    FlueGasParameterBlock,
)


C_TO_K = 273.15


class Mode(Enum):
    """Supported refrigerant state variable modes."""

    PH = "PH"
    IMPROVED_TPX = "improved_TPx"


class SimpleVaporCompressionCyclePLRNTU:
    """Copy-only PLR cycle with single NTU evaporator and condenser blocks."""

    def __init__(
        self,
        fluid_name,
        compressor_efficiency=0.75,
        PLR=0.75,
        CD=0.13,
        mode=Mode.PH,
        UA_evap_W_per_K=12000.0,
        UA_cond_W_per_K=15000.0,
    ):
        assert 0.0 < compressor_efficiency < 1.0, "Compressor efficiency must be in (0,1)"
        assert 0.0 <= PLR <= 1.0, "PLR must be in [0,1]"
        assert 0.0 <= CD <= 1.0, "CD must be in [0,1]"

        self.fluid_name = fluid_name
        self.compressor_efficiency = compressor_efficiency
        self.plr = PLR
        self.cd = CD
        self.mode = mode

        self.model = ConcreteModel()
        self.model.fs = FlowsheetBlock(dynamic=False)

        state_vars = StateVars.PH if mode == Mode.PH else StateVars.TPX
        self.model.fs.ref_props = HelmholtzParameterBlock(
            pure_component=fluid_name,
            state_vars=state_vars,
            amount_basis=AmountBasis.MASS,
        )
        self.model.fs.air_props = FlueGasParameterBlock()

        self.model.fs.ua_evap = Param(
            initialize=UA_evap_W_per_K,
            mutable=True,
            units=pyunits.W / pyunits.K,
        )
        self.model.fs.ua_cond = Param(
            initialize=UA_cond_W_per_K,
            mutable=True,
            units=pyunits.W / pyunits.K,
        )

        self._define_flowsheet()
        self.optimization_converged = None

    @staticmethod
    def _compute_plf(plr: float, cd: float) -> float:
        """Compute PLF from AHRI-style linear degradation relation."""

        if not (0.0 <= plr <= 1.0):
            raise ValueError("PLR must be in [0, 1]")
        if not (0.0 <= cd <= 1.0):
            raise ValueError("CD must be in [0, 1]")
        return max(0.0, min(1.0, 1.0 - cd * (1.0 - plr)))

    def _define_flowsheet(self):
        self.logger = logging.getLogger(__name__)

        fs = self.model.fs
        fs.evaporator = HeatExchangerNTU(
            hot_side={"property_package": fs.air_props, "has_pressure_change": False},
            cold_side={"property_package": fs.ref_props, "has_pressure_change": False},
        )
        fs.compressor = Compressor(property_package=fs.ref_props)
        fs.condenser = HeatExchangerNTU(
            hot_side={"property_package": fs.ref_props, "has_pressure_change": False},
            cold_side={"property_package": fs.air_props, "has_pressure_change": False},
        )
        fs.expansion_valve = PressureChanger(
            property_package=fs.ref_props,
            thermodynamic_assumption="adiabatic",
            compressor=False,
        )

        fs.evaporator_to_compressor = Arc(
            source=fs.evaporator.cold_side_outlet,
            destination=fs.compressor.inlet,
        )
        fs.compressor_to_condenser = Arc(
            source=fs.compressor.outlet,
            destination=fs.condenser.hot_side_inlet,
        )
        fs.condenser_to_expansion_valve = Arc(
            source=fs.condenser.hot_side_outlet,
            destination=fs.expansion_valve.inlet,
        )
        fs.expansion_valve_to_evaporator = Arc(
            source=fs.expansion_valve.outlet,
            destination=fs.evaporator.cold_side_inlet,
        )

        from pyomo.environ import TransformationFactory

        TransformationFactory("network.expand_arcs").apply_to(self.model)

        # Closed-loop flow equality is redundant.
        for cname in ["flow_mass_equality", "flow_mol_equality"]:
            c = getattr(fs.evaporator_to_compressor_expanded, cname, None)
            if c is not None:
                c.deactivate()

        # Use explicit high/low pressure variables.
        for arc in [
            fs.evaporator_to_compressor_expanded,
            fs.compressor_to_condenser_expanded,
            fs.condenser_to_expansion_valve_expanded,
            fs.expansion_valve_to_evaporator_expanded,
        ]:
            if hasattr(arc, "pressure_equality"):
                arc.pressure_equality.deactivate()

        fs.P_low = Var(initialize=200e3, bounds=(50e3, 2000e3), units=pyunits.Pa)
        fs.P_high = Var(initialize=1200e3, bounds=(100e3, 5000e3), units=pyunits.Pa)

        fs.P_low_target = Param(initialize=200e3, mutable=True, units=pyunits.Pa)
        fs.P_high_target = Param(initialize=1200e3, mutable=True, units=pyunits.Pa)

        fs.P_low_target_constraint = Constraint(expr=fs.P_low == fs.P_low_target)
        fs.P_high_target_constraint = Constraint(expr=fs.P_high == fs.P_high_target)

        fs.P_low_evap_in = Constraint(expr=fs.evaporator.cold_side_inlet.pressure[0] == fs.P_low)
        fs.P_low_evap_out = Constraint(expr=fs.evaporator.cold_side_outlet.pressure[0] == fs.P_low)
        fs.P_low_comp_in = Constraint(expr=fs.compressor.inlet.pressure[0] == fs.P_low)
        fs.P_low_valve_out = Constraint(expr=fs.expansion_valve.outlet.pressure[0] == fs.P_low)

        fs.P_high_comp_out = Constraint(expr=fs.compressor.outlet.pressure[0] == fs.P_high)
        fs.P_high_cond_in = Constraint(expr=fs.condenser.hot_side_inlet.pressure[0] == fs.P_high)
        fs.P_high_cond_out = Constraint(expr=fs.condenser.hot_side_outlet.pressure[0] == fs.P_high)
        fs.P_high_valve_in = Constraint(expr=fs.expansion_valve.inlet.pressure[0] == fs.P_high)
        # Redundant with HX no-pressure-drop equations plus inlet pressure constraints.
        fs.P_low_evap_out.deactivate()
        fs.P_high_cond_out.deactivate()

        fs.max_pressure_ratio = Param(initialize=20.0, mutable=True, units=pyunits.dimensionless)
        fs.pressure_ratio_constraint = Constraint(
            expr=fs.compressor.outlet.pressure[0]
            <= fs.max_pressure_ratio * fs.compressor.inlet.pressure[0]
        )

        # Refrigerant-side superheat/subcool and ambient approach hooks.
        fs.superheating = Param(initialize=3.0, mutable=True, units=pyunits.K)
        fs.subcooling = Param(initialize=3.0, mutable=True, units=pyunits.K)
        fs.ambient_T = Param(initialize=C_TO_K + 20.0, mutable=True, units=pyunits.K)
        fs.approach_T = Param(initialize=10.0, mutable=True, units=pyunits.K)

        fs.superheating_constraint = Constraint(
            expr=fs.evaporator.cold_side.properties_out[0].temperature
            >= fs.evaporator.cold_side.properties_out[0].temperature_sat + fs.superheating
        )
        fs.subcooling_constraint = Constraint(
            expr=fs.condenser.hot_side.properties_out[0].temperature
            <= fs.condenser.hot_side.properties_out[0].temperature_sat - fs.subcooling
        )
        fs.approach_constraint = Constraint(
            expr=fs.condenser.hot_side.properties_out[0].temperature_sat
            == fs.ambient_T + fs.approach_T
        )

        fs.evap_Tmin = Param(initialize=C_TO_K - 60.0, mutable=True, units=pyunits.K)
        fs.evap_Tmax = Param(initialize=C_TO_K - 10.0, mutable=True, units=pyunits.K)
        fs.cond_Tmin = Param(initialize=C_TO_K + 15.0, mutable=True, units=pyunits.K)
        fs.cond_Tmax = Param(initialize=C_TO_K + 70.0, mutable=True, units=pyunits.K)

        fs.evap_T_lower = Constraint(
            expr=fs.evaporator.cold_side.properties_out[0].temperature >= fs.evap_Tmin
        )
        fs.evap_T_upper = Constraint(
            expr=fs.evaporator.cold_side.properties_out[0].temperature <= fs.evap_Tmax
        )
        fs.cond_T_lower = Constraint(
            expr=fs.condenser.hot_side.properties_out[0].temperature >= fs.cond_Tmin
        )
        fs.cond_T_upper = Constraint(
            expr=fs.condenser.hot_side.properties_out[0].temperature <= fs.cond_Tmax
        )

        # Vapor-only compressor outlet guard.
        fs.vapor_constraint = Constraint(
            expr=fs.compressor.control_volume.properties_out[0].temperature
            >= fs.compressor.control_volume.properties_out[0].temperature_sat
        )

        # Phase-change dominant simplification for stability: epsilon = 1 - exp(-NTU).
        fs.evaporator_eps_relation = Constraint(
            expr=fs.evaporator.effectiveness[0] == (1.0 - exp(-fs.evaporator.NTU[0]))
        )
        fs.condenser_eps_relation = Constraint(
            expr=fs.condenser.effectiveness[0] == (1.0 - exp(-fs.condenser.NTU[0]))
        )

        fs.evaporator.area.fix(1.0)
        fs.condenser.area.fix(1.0)
        fs.evaporator.heat_transfer_coefficient[0].fix(value(fs.ua_evap))
        fs.condenser.heat_transfer_coefficient[0].fix(value(fs.ua_cond))
        fs.evaporator.effectiveness[0].setlb(1e-3)
        fs.evaporator.effectiveness[0].setub(0.999)
        fs.condenser.effectiveness[0].setlb(1e-3)
        fs.condenser.effectiveness[0].setub(0.999)
        fs.evaporator.effectiveness[0].fix(0.85)
        fs.condenser.effectiveness[0].fix(0.85)
        fs.evaporator_eps_relation.deactivate()
        fs.condenser_eps_relation.deactivate()

        fs.PLF = Param(initialize=1.0, mutable=True, units=pyunits.dimensionless)
        fs.plr = Param(initialize=self.plr, mutable=True, units=pyunits.dimensionless)
        fs.cd = Param(initialize=self.cd, mutable=True, units=pyunits.dimensionless)

        fs.cop = Var(initialize=3.0, bounds=(0.01, 100.0), units=pyunits.dimensionless)
        fs.compute_cop = Constraint(
            expr=fs.cop * fs.compressor.work_mechanical[0] == fs.evaporator.heat_duty[0]
        )
        fs.obj = Objective(expr=fs.cop, sense=maximize)
        fs.compute_cop.deactivate()
        fs.obj.deactivate()

        self._set_air_side_defaults()

    def _set_air_side_defaults(self):
        """Fix air-side inlets for both NTU blocks using dry-air composition."""

        fs = self.model.fs
        evap_in = fs.evaporator.hot_side_inlet
        cond_in = fs.condenser.cold_side_inlet

        # Nominal dry air composition in FlueGas package components.
        comp = {
            "H2O": 0.0,
            "CO2": 0.0004,
            "N2": 0.7900,
            "O2": 0.2096,
            "NO": 0.0,
            "SO2": 0.0,
        }

        evap_total = 70.0  # mol/s
        cond_total = 90.0  # mol/s

        for j, x in comp.items():
            evap_in.flow_mol_comp[0, j].fix(evap_total * x)
            cond_in.flow_mol_comp[0, j].fix(cond_total * x)

        evap_in.pressure[0].fix(101325.0)
        cond_in.pressure[0].fix(101325.0)
        evap_in.temperature[0].fix(C_TO_K - 20.0)
        cond_in.temperature[0].fix(C_TO_K + 20.0)

    def specify_initial_conditions(self, low_side_temperature=-20.0, high_side_temperature=30.0):
        """Set CoolProp-based refrigerant guesses for start-up initialization."""

        TL = low_side_temperature + C_TO_K
        TH = high_side_temperature + C_TO_K
        superheat = 3.0
        subcool = 3.0

        p_low = CP.PropsSI("P", "T", TL, "Q", 0, self.fluid_name)
        p_high = CP.PropsSI("P", "T", TH, "Q", 0, self.fluid_name)

        h4 = CP.PropsSI("H", "T", TL, "Q", 0.2, self.fluid_name)
        h1 = CP.PropsSI("H", "T", TL + superheat, "Q", 1, self.fluid_name)
        h2 = CP.PropsSI("H", "T", TH + superheat, "Q", 1, self.fluid_name)
        h3 = CP.PropsSI("H", "T", TH - subcool, "Q", 0, self.fluid_name)

        self._init = {
            "p_low": p_low,
            "p_high": p_high,
            "h1": h1,
            "h2": h2,
            "h3": h3,
            "h4": h4,
        }

    def initialize(self, verbose=False):
        """Load seed values for state variables without unit-level initializers."""

        fs = self.model.fs
        init = getattr(self, "_init", None)
        if init is None:
            self.specify_initial_conditions(low_side_temperature=-20.0, high_side_temperature=30.0)
            init = self._init

        fs.evaporator.cold_side_inlet.flow_mass[0].set_value(1.0)
        fs.evaporator.cold_side_inlet.pressure[0].set_value(init["p_low"])
        fs.evaporator.cold_side_inlet.enth_mass[0].set_value(init["h4"])
        fs.evaporator.cold_side_outlet.enth_mass[0].set_value(init["h1"])

        fs.compressor.inlet.pressure[0].set_value(init["p_low"])
        fs.compressor.inlet.enth_mass[0].set_value(init["h1"])
        fs.compressor.outlet.pressure[0].set_value(init["p_high"])
        fs.compressor.efficiency_isentropic[0].fix(self.compressor_efficiency)

        fs.condenser.hot_side_inlet.pressure[0].set_value(init["p_high"])
        fs.condenser.hot_side_inlet.enth_mass[0].set_value(init["h2"])
        fs.condenser.hot_side_outlet.enth_mass[0].set_value(init["h3"])

        fs.expansion_valve.inlet.pressure[0].set_value(init["p_high"])
        fs.expansion_valve.outlet.pressure[0].set_value(init["p_low"])

        # Keep epsilon fixed during unit initialization; activate correlation later.
        fs.evaporator_eps_relation.deactivate()
        fs.condenser_eps_relation.deactivate()
        fs.evaporator.effectiveness[0].fix(0.85)
        fs.condenser.effectiveness[0].fix(0.85)

        # Keep stream propagation lightweight; full solve occurs in optimize_COP.
        try:
            propagate_state(fs.evaporator_to_compressor)
            propagate_state(fs.compressor_to_condenser)
            propagate_state(fs.condenser_to_expansion_valve)
            propagate_state(fs.expansion_valve_to_evaporator)
        except Exception:
            pass

        if verbose:
            self.logger.info("Seed values loaded for NTU copy model.")

    def set_specifications(
        self,
        low_side_pressure=(200, 500),
        high_side_pressure=(1000, 3000),
        evaporator_temperature=(-20, 0),
        cold_storage_setpoint=None,
        evap_offset_bounds=(-35.0, 0.0),
        compressor_temperature=None,
        condenser_temperature=(30, 50),
        expansion_valve_temperature=None,
        subcooling=3,
        superheating=3,
        max_pressure_ratio=4,
        ambient_temperature=None,
        condenser_approach=None,
        evap_sat_temperature=None,
        debug_disable_arc_pressure_eq=False,
        plr=None,
        cd=None,
        UA_evap_total=None,
        UA_cond_total=None,
    ):
        """Apply operating specs while keeping PLR as post-correction only."""

        fs = self.model.fs

        if plr is not None:
            assert 0.0 <= plr <= 1.0
            self.plr = plr
        if cd is not None:
            assert 0.0 <= cd <= 1.0
            self.cd = cd

        fs.plr.set_value(self.plr)
        fs.cd.set_value(self.cd)
        fs.PLF.set_value(self._compute_plf(self.plr, self.cd))

        if UA_evap_total is not None:
            fs.ua_evap.set_value(UA_evap_total)
            fs.evaporator.heat_transfer_coefficient[0].fix(float(UA_evap_total))
        if UA_cond_total is not None:
            fs.ua_cond.set_value(UA_cond_total)
            fs.condenser.heat_transfer_coefficient[0].fix(float(UA_cond_total))

        # Unfix initialization anchors.
        fs.evaporator.cold_side_inlet.flow_mass[0].unfix()
        fs.evaporator.cold_side_inlet.pressure[0].unfix()
        fs.evaporator.cold_side_inlet.enth_mass[0].unfix()
        fs.evaporator.cold_side_outlet.enth_mass[0].unfix()

        fs.compressor.inlet.pressure[0].unfix()
        fs.compressor.inlet.enth_mass[0].unfix()
        fs.compressor.outlet.pressure[0].unfix()

        fs.condenser.hot_side_inlet.pressure[0].unfix()
        fs.condenser.hot_side_inlet.enth_mass[0].unfix()
        fs.condenser.hot_side_outlet.enth_mass[0].unfix()

        fs.expansion_valve.inlet.pressure[0].unfix()
        fs.expansion_valve.outlet.pressure[0].unfix()

        # Keep reference flow fixed (PLR is post-correction only).
        fs.evaporator.cold_side_inlet.flow_mass[0].fix(1.0)

        # Set pressure bounds and default targets.
        p_low_min, p_low_max = low_side_pressure
        p_high_min, p_high_max = high_side_pressure
        fs.P_low.setlb(p_low_min * 1e3)
        fs.P_low.setub(p_low_max * 1e3)
        fs.P_high.setlb(p_high_min * 1e3)
        fs.P_high.setub(p_high_max * 1e3)

        fs.P_low_target.set_value((p_low_min + p_low_max) * 0.5 * 1e3)

        # If condenser approach is active, let P_high be implied by Tsat relation.
        if ambient_temperature is not None and condenser_approach is not None:
            fs.ambient_T.set_value(ambient_temperature + C_TO_K)
            fs.approach_T.set_value(condenser_approach)
            fs.approach_constraint.activate()
            fs.P_high_target_constraint.deactivate()
            fs.condenser.cold_side_inlet.temperature[0].fix(ambient_temperature + C_TO_K)
        else:
            fs.approach_constraint.deactivate()
            fs.P_high_target.set_value((p_high_min + p_high_max) * 0.5 * 1e3)
            fs.P_high_target_constraint.activate()
        # Tie evaporator air inlet to cold-storage setpoint when provided.
        if cold_storage_setpoint is not None:
            fs.evaporator.hot_side_inlet.temperature[0].fix(cold_storage_setpoint + C_TO_K)
        else:
            fs.evaporator.hot_side_inlet.temperature[0].fix(C_TO_K - 20.0)
        if ambient_temperature is not None:
            fs.condenser.cold_side_inlet.temperature[0].fix(ambient_temperature + C_TO_K)

        fs.superheating.set_value(superheating)
        fs.subcooling.set_value(subcooling)
        fs.max_pressure_ratio.set_value(max_pressure_ratio)
        # Keep NTU copy numerically stable: temperature inequalities only.
        fs.superheating_constraint.deactivate()
        fs.subcooling_constraint.deactivate()
        fs.vapor_constraint.deactivate()

        if cold_storage_setpoint is not None:
            evap_lo = cold_storage_setpoint + evap_offset_bounds[0]
            evap_hi = cold_storage_setpoint + evap_offset_bounds[1]
        else:
            evap_lo, evap_hi = evaporator_temperature
        fs.evap_Tmin.set_value(evap_lo + C_TO_K)
        fs.evap_Tmax.set_value(evap_hi + C_TO_K)
        fs.cond_Tmin.set_value(condenser_temperature[0] + C_TO_K)
        fs.cond_Tmax.set_value(condenser_temperature[1] + C_TO_K)

        if evap_sat_temperature is not None:
            # Soft hook: pin low pressure from requested saturation temperature.
            p_low_sat = CP.PropsSI("P", "T", evap_sat_temperature + C_TO_K, "Q", 1, self.fluid_name)
            fs.P_low_target.set_value(p_low_sat)
            fs.P_low_target_constraint.activate()
        else:
            fs.P_low_target_constraint.activate()

        fs.compressor.efficiency_isentropic[0].fix(self.compressor_efficiency)
        fs.evaporator.effectiveness[0].unfix()
        fs.condenser.effectiveness[0].unfix()
        fs.evaporator_eps_relation.activate()
        fs.condenser_eps_relation.activate()

        calculate_scaling_factors(self.model)

    def optimize_COP(self, verbose=False, initialize=True, optimize=False):
        """Solve cycle; COP_part remains PLF*COP_full post-correction."""

        solver = get_solver()
        solver.options = {
            "max_iter": 4000,
            "tol": 1e-6,
            "acceptable_tol": 1e-5,
        }

        if initialize:
            self.model.fs.compute_cop.deactivate()
            self.model.fs.obj.deactivate()
            try:
                solver.solve(self.model, tee=verbose)
            except Exception:
                pass

        if optimize:
            self.model.fs.compute_cop.activate()
            self.model.fs.obj.activate()
        else:
            self.model.fs.compute_cop.deactivate()
            self.model.fs.obj.deactivate()

        try:
            results = solver.solve(self.model, tee=verbose)
            if results.solver.termination_condition != TerminationCondition.optimal:
                results = solver.solve(self.model, tee=verbose)
        except Exception:
            self.optimization_converged = False
            return float("nan"), False

        if (
            self.model.fs.compressor.work_mechanical[0].value is not None
            and self.model.fs.compressor.work_mechanical[0].value != 0
        ):
            self.model.fs.cop.set_value(
                self.model.fs.evaporator.heat_duty[0].value
                / self.model.fs.compressor.work_mechanical[0].value
            )

        converged = results.solver.termination_condition == TerminationCondition.optimal
        self.optimization_converged = converged

        if not converged:
            try:
                DiagnosticsToolbox(
                    self.model, constraint_residual_tolerance=1e-6
                ).display_constraints_with_large_residuals()
            except Exception:
                pass

        cop_full = value(self.model.fs.cop)
        self._last_cop_full = cop_full
        self._last_plf = self._compute_plf(value(self.model.fs.plr), value(self.model.fs.cd))
        self._last_cop_part = self._last_plf * cop_full

        # Report actual SH/SC from solved states (computed outputs, no forcing).
        try:
            tev_out = value(self.model.fs.evaporator.cold_side.properties_out[0].temperature)
            tev_sat = value(self.model.fs.evaporator.cold_side.properties_out[0].temperature_sat)
            tcd_out = value(self.model.fs.condenser.hot_side.properties_out[0].temperature)
            tcd_sat = value(self.model.fs.condenser.hot_side.properties_out[0].temperature_sat)
            self._last_sh_actual = max(0.0, tev_out - tev_sat)
            self._last_sc_actual = max(0.0, tcd_sat - tcd_out)
        except Exception:
            self._last_sh_actual = float("nan")
            self._last_sc_actual = float("nan")

        return cop_full, converged

    def get_full_load_cop(self):
        """Return full-load COP from current solve state."""

        return value(self.model.fs.cop)

    def get_part_load_cop(self):
        """Return part-load corrected COP = PLF * COP_full."""

        return self._compute_plf(value(self.model.fs.plr), value(self.model.fs.cd)) * value(self.model.fs.cop)

    def get_actual_sh_sc(self):
        """Return last solved SH/SC in K."""

        return getattr(self, "_last_sh_actual", float("nan")), getattr(self, "_last_sc_actual", float("nan"))
