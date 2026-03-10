"""
Lumped-Capacitance Heat Exchanger Vapor Compression PLR Model (Copy)

Author: Shilpa Narasimhan
Codex Support: OpenAI Codex
QA/Testing: Shilpa Narasimhan

Description:
    Copy-only vapor-compression refrigeration cycle model that preserves the
    PLR post-correction workflow while replacing heater-style evaporator and
    condenser units with IDAES HeatExchangerLumpedCapacitance units.

Context Breadcrumb:
    This variant uses IDAES lumped-capacitance HX blocks to mimic the prior
    non-IDAES lumped-UA workflow while keeping compressor/expansion and PLR
    post-correction structure aligned with the PLR code path.

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
    TransformationFactory,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.scaling import set_scaling_factor
from pyomo.opt import TerminationCondition
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    AmountBasis,
    StateVars,
)
from idaes.models.unit_models import (
    Compressor,
    HeatExchangerLumpedCapacitance,
    PressureChanger,
    HeatExchangerFlowPattern,
)
from idaes.models_extra.power_generation.properties.flue_gas_ideal import (
    FlueGasParameterBlock,
)


C_TO_K = 273.15


class Mode(Enum):
    """Supported refrigerant state variable modes.

    Attributes:
        PH: Pressure-enthalpy state variables. Preferred for robust cycle solves.
        IMPROVED_TPX: Temperature-pressure-quality style variables.
    """

    PH = "PH"
    IMPROVED_TPX = "improved_TPx"


class SimpleVaporCompressionCyclePLRHXLC:
    """Copy-only PLR cycle using HeatExchangerLumpedCapacitance for both coils.

    Mathematical notes:
        - COP is computed from refrigerant-side evaporator enthalpy rise:
          COP = Q_evap / W_comp
          Q_evap = m_dot_ref * (h_evap_out - h_evap_in)
        - Part-load COP post-correction uses:
          PLF = 1 - CD * (1 - PLR)
          COP_part = PLF * COP_full

    The compressor and expansion valve formulations are kept aligned with the
    existing PLR model family; only coil unit models are swapped to lumped
    capacitance HX units with explicit hot/cold UA parameters.
    """

    def __init__(
        self,
        fluid_name,
        compressor_efficiency=0.75,
        PLR=0.75,
        CD=0.13,
        mode=Mode.PH,
        evap_area_m2=6.0,
        cond_area_m2=100.0,
        UA_evap_hot=500.0,
        UA_evap_cold=1000.0,
        UA_cond_hot=900.0,
        UA_cond_cold=900.0,
        wall_heat_capacity_evap=5.0e4,
        wall_heat_capacity_cond=2.0e5,
        wall_temperature_evap_init=C_TO_K - 20.0,
        wall_temperature_cond_init=C_TO_K + 35.0,
        evap_air_mol_s=31.0,
        cond_air_mol_s=356.0,
    ):
        """Construct the lumped-capacitance HX PLR cycle copy.

        Args:
            fluid_name: Refrigerant name compatible with IDAES Helmholtz package.
            compressor_efficiency: Isentropic efficiency in (0,1).
            PLR: Part-load ratio in [0,1].
            CD: Cycling degradation coefficient in [0,1].
            mode: Refrigerant state variable mode (`PH` recommended).
            evap_area_m2/cond_area_m2: HX areas used in UA closure equations.
            UA_evap_hot/UA_evap_cold: Evaporator hot/cold side UA values (W/K).
            UA_cond_hot/UA_cond_cold: Condenser hot/cold side UA values (W/K).
            wall_heat_capacity_*: Lumped wall heat capacities (J/K).
            wall_temperature_*_init: Initial wall temperature seeds (K).
            evap_air_mol_s/cond_air_mol_s: Air-side molar flow anchors.
        """
        assert 0.0 < compressor_efficiency < 1.0
        assert 0.0 <= PLR <= 1.0
        assert 0.0 <= CD <= 1.0

        self.fluid_name = fluid_name
        self.idaes_fluid_name, self.cp_fluid_name = self._map_fluid_names(fluid_name)
        self.compressor_efficiency = compressor_efficiency
        self.plr = PLR
        self.cd = CD
        self.mode = mode

        self.model = ConcreteModel()
        self.model.fs = FlowsheetBlock(dynamic=False)

        sv = StateVars.PH if mode == Mode.PH else StateVars.TPX
        self.model.fs.ref_props = HelmholtzParameterBlock(
            pure_component=self.idaes_fluid_name,
            state_vars=sv,
            amount_basis=AmountBasis.MASS,
        )
        self.model.fs.air_props = FlueGasParameterBlock()

        fs = self.model.fs

        # Lumped-capacitance HX parameter anchors.
        fs.evap_area = Param(initialize=evap_area_m2, mutable=True, units=pyunits.m**2)
        fs.ua_evap_hot = Param(initialize=UA_evap_hot, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_evap_cold = Param(initialize=UA_evap_cold, mutable=True, units=pyunits.W / pyunits.K)
        fs.cw_evap = Param(initialize=wall_heat_capacity_evap, mutable=True, units=pyunits.J / pyunits.K)
        fs.tw_evap_init = Param(initialize=wall_temperature_evap_init, mutable=True, units=pyunits.K)
        fs.evap_air_mol_s = Param(initialize=evap_air_mol_s, mutable=True, units=pyunits.mol / pyunits.s)

        fs.cond_area = Param(initialize=cond_area_m2, mutable=True, units=pyunits.m**2)
        fs.ua_cond_hot = Param(initialize=UA_cond_hot, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_cond_cold = Param(initialize=UA_cond_cold, mutable=True, units=pyunits.W / pyunits.K)
        fs.cw_cond = Param(initialize=wall_heat_capacity_cond, mutable=True, units=pyunits.J / pyunits.K)
        fs.tw_cond_init = Param(initialize=wall_temperature_cond_init, mutable=True, units=pyunits.K)
        fs.cond_air_mol_s = Param(initialize=cond_air_mol_s, mutable=True, units=pyunits.mol / pyunits.s)

        self.logger = logging.getLogger(__name__)
        self.optimization_converged = None

        self._define_flowsheet()

    @staticmethod
    def _map_fluid_names(fluid_name: str):
        """Map external fluid names to IDAES and CoolProp conventions."""
        key = fluid_name.strip().lower()
        if key in {"r1234ze(e)", "r1234zee", "r1234zee", "r1234zee"}:
            return "R1234ze", "R1234ze(E)"
        if key == "r1234ze":
            return "R1234ze", "R1234ze(E)"
        return fluid_name, fluid_name

    @staticmethod
    def _compute_plf(plr: float, cd: float) -> float:
        """Compute AHRI-style part-load factor.

        Uses:
            PLF = 1 - CD * (1 - PLR)

        Args:
            plr: Part-load ratio in [0, 1].
            cd: Cycling degradation coefficient in [0, 1].

        Returns:
            Dimensionless part-load factor clipped to [0, 1].
        """
        if not (0.0 <= plr <= 1.0):
            raise ValueError("PLR must be in [0, 1]")
        if not (0.0 <= cd <= 1.0):
            raise ValueError("CD must be in [0, 1]")
        return max(0.0, min(1.0, 1.0 - cd * (1.0 - plr)))

    def _define_flowsheet(self):
        """Build the steady-state lumped-capacitance HX cycle structure."""
        fs = self.model.fs

        fs.evaporator = HeatExchangerLumpedCapacitance(
            hot_side={
                "property_package": fs.air_props,
                "has_phase_equilibrium": False,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            cold_side={
                "property_package": fs.ref_props,
                "has_phase_equilibrium": False,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
            dynamic_heat_balance=False,
        )

        fs.compressor = Compressor(property_package=fs.ref_props)

        fs.condenser = HeatExchangerLumpedCapacitance(
            hot_side={
                "property_package": fs.ref_props,
                "has_phase_equilibrium": False,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            cold_side={
                "property_package": fs.air_props,
                "has_phase_equilibrium": False,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
            dynamic_heat_balance=False,
        )

        fs.expansion_valve = PressureChanger(
            property_package=fs.ref_props,
            thermodynamic_assumption="adiabatic",
            compressor=False,
        )

        # Refrigerant loop connectivity.
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

        TransformationFactory("network.expand_arcs").apply_to(self.model)

        # Closed-loop mass equality on one arc is redundant.
        for cname in ["flow_mass_equality", "flow_mol_equality"]:
            c = getattr(fs.evaporator_to_compressor_expanded, cname, None)
            if c is not None:
                c.deactivate()

        # Use explicit pressure-level constraints for robustness.
        for arc in [
            fs.evaporator_to_compressor_expanded,
            fs.compressor_to_condenser_expanded,
            fs.condenser_to_expansion_valve_expanded,
            fs.expansion_valve_to_evaporator_expanded,
        ]:
            if hasattr(arc, "pressure_equality"):
                arc.pressure_equality.deactivate()

        self._fix_hx_lc_parameters()
        self._fix_air_side_defaults()
        self._build_common_cycle_constraints()

    def _fix_hx_lc_parameters(self):
        """Fix area, UA side parameters, and wall terms for LC HX blocks."""
        fs = self.model.fs

        fs.evaporator.area.fix(value(fs.evap_area))
        fs.condenser.area.fix(value(fs.cond_area))
        fs.evaporator.heat_capacity_wall.set_value(value(fs.cw_evap))
        fs.condenser.heat_capacity_wall.set_value(value(fs.cw_cond))
        fs.evaporator.thermal_resistance_wall.set_value(0.0)
        fs.condenser.thermal_resistance_wall.set_value(0.0)
        fs.evaporator.thermal_fouling_hot_side.set_value(0.0)
        fs.evaporator.thermal_fouling_cold_side.set_value(0.0)
        fs.condenser.thermal_fouling_hot_side.set_value(0.0)
        fs.condenser.thermal_fouling_cold_side.set_value(0.0)

        for t in fs.time:
            fs.evaporator.ua_hot_side[t].fix(value(fs.ua_evap_hot))
            fs.evaporator.ua_cold_side[t].fix(value(fs.ua_evap_cold))
            fs.condenser.ua_hot_side[t].fix(value(fs.ua_cond_hot))
            fs.condenser.ua_cold_side[t].fix(value(fs.ua_cond_cold))
            fs.evaporator.temperature_wall[t].set_value(value(fs.tw_evap_init))
            fs.condenser.temperature_wall[t].set_value(value(fs.tw_cond_init))

    def _fix_air_side_defaults(self):
        """Fix default air-side composition and nominal inlets for steady solves."""
        fs = self.model.fs

        comp = {
            "H2O": 0.0,
            "CO2": 0.0004,
            "N2": 0.7900,
            "O2": 0.2096,
            "NO": 0.0,
            "SO2": 0.0,
        }

        evap_total = value(fs.evap_air_mol_s)
        cond_total = value(fs.cond_air_mol_s)

        for j, x in comp.items():
            fs.evaporator.hot_side_inlet.flow_mol_comp[0, j].fix(evap_total * x)
            fs.condenser.cold_side_inlet.flow_mol_comp[0, j].fix(cond_total * x)

        fs.evaporator.hot_side_inlet.pressure[0].fix(101325.0)
        fs.condenser.cold_side_inlet.pressure[0].fix(101325.0)

        # Requested boundary convention:
        # evaporator air at cold-storage setpoint, condenser air at ambient.
        fs.evaporator.hot_side_inlet.temperature[0].fix(C_TO_K - 20.0)
        fs.condenser.cold_side_inlet.temperature[0].fix(C_TO_K + 20.0)

    def _build_common_cycle_constraints(self):
        """Add shared cycle constraints and performance expressions.

        Includes:
            - Low/high pressure decision variables (`P_low`, `P_high`)
            - Compressor pressure-ratio cap
            - Optional condenser approach-temperature coupling
            - Evaporator/condenser temperature window constraints
            - Refrigeration COP relation

        COP definition:
            COP = Q_evap / W_comp
            Q_evap = m_dot_ref * (h_evap_out - h_evap_in)

        PLR/CD are stored as parameters and used only in post-correction
        reporting (`get_part_load_cop`) to avoid double counting inside the
        thermodynamic balances.
        """
        fs = self.model.fs

        fs.P_low = Var(initialize=200e3, bounds=(50e3, 2000e3), units=pyunits.Pa)
        fs.P_high = Var(initialize=1200e3, bounds=(100e3, 5000e3), units=pyunits.Pa)

        fs.P_low_target = Param(initialize=200e3, mutable=True, units=pyunits.Pa)
        fs.P_high_target = Param(initialize=1200e3, mutable=True, units=pyunits.Pa)

        fs.P_low_target_constraint = Constraint(expr=fs.P_low == fs.P_low_target)
        fs.P_high_target_constraint = Constraint(expr=fs.P_high == fs.P_high_target)

        fs.P_low_evap_in = Constraint(
            expr=fs.evaporator.cold_side_inlet.pressure[0]
            == fs.expansion_valve.outlet.pressure[0]
        )
        fs.P_low_comp_in = Constraint(expr=fs.compressor.inlet.pressure[0] == fs.P_low)
        fs.P_low_valve_out = Constraint(expr=fs.expansion_valve.outlet.pressure[0] == fs.P_low)

        fs.P_high_comp_out = Constraint(expr=fs.compressor.outlet.pressure[0] == fs.P_high)
        fs.P_high_cond_in = Constraint(expr=fs.condenser.hot_side_inlet.pressure[0] == fs.P_high)
        fs.P_high_valve_in = Constraint(expr=fs.expansion_valve.inlet.pressure[0] == fs.P_high)

        fs.max_pressure_ratio = Param(initialize=20.0, mutable=True, units=pyunits.dimensionless)
        fs.pressure_ratio_constraint = Constraint(
            expr=fs.compressor.outlet.pressure[0]
            <= fs.max_pressure_ratio * fs.compressor.inlet.pressure[0]
        )

        fs.ambient_T = Param(initialize=C_TO_K + 20.0, mutable=True, units=pyunits.K)
        fs.approach_T = Param(initialize=30.0, mutable=True, units=pyunits.K)
        fs.approach_constraint = Constraint(
            expr=fs.condenser.hot_side.properties_out[0].temperature_sat
            == fs.ambient_T + fs.approach_T
        )
        fs.approach_constraint.deactivate()

        fs.evap_Tmin = Param(initialize=C_TO_K - 60.0, mutable=True, units=pyunits.K)
        fs.evap_Tmax = Param(initialize=C_TO_K - 10.0, mutable=True, units=pyunits.K)
        fs.cond_Tmin = Param(initialize=C_TO_K + 15.0, mutable=True, units=pyunits.K)
        fs.cond_Tmax = Param(initialize=C_TO_K + 70.0, mutable=True, units=pyunits.K)

        fs.evap_T_lower = Constraint(expr=fs.evaporator.cold_side.properties_out[0].temperature >= fs.evap_Tmin)
        fs.evap_T_upper = Constraint(expr=fs.evaporator.cold_side.properties_out[0].temperature <= fs.evap_Tmax)
        fs.cond_T_lower = Constraint(expr=fs.condenser.hot_side.properties_out[0].temperature >= fs.cond_Tmin)
        fs.cond_T_upper = Constraint(expr=fs.condenser.hot_side.properties_out[0].temperature <= fs.cond_Tmax)
        fs.evap_T_lower.deactivate()
        fs.evap_T_upper.deactivate()
        fs.cond_T_lower.deactivate()
        fs.cond_T_upper.deactivate()

        fs.vapor_constraint = Constraint(
            expr=fs.compressor.control_volume.properties_out[0].temperature
            >= fs.compressor.control_volume.properties_out[0].temperature_sat
        )
        fs.vapor_constraint.deactivate()

        fs.plr = Param(initialize=self.plr, mutable=True, units=pyunits.dimensionless)
        fs.cd = Param(initialize=self.cd, mutable=True, units=pyunits.dimensionless)

        fs.cop = Var(initialize=3.0, bounds=(0.01, 100.0), units=pyunits.dimensionless)
        fs.compute_cop = Constraint(
            expr=fs.cop * fs.compressor.work_mechanical[0]
            == fs.evaporator.cold_side_inlet.flow_mass[0]
            * (fs.evaporator.cold_side_outlet.enth_mass[0] - fs.evaporator.cold_side_inlet.enth_mass[0])
        )
        fs.obj = Objective(expr=fs.cop, sense=maximize)
        fs.compute_cop.deactivate()
        fs.obj.deactivate()

    def specify_initial_conditions(self, low_side_temperature=-20.0, high_side_temperature=30.0):
        """Generate initial pressure/enthalpy guesses from saturation anchors.

        Args:
            low_side_temperature: Low-side saturation anchor in degC.
            high_side_temperature: High-side saturation anchor in degC.

        Notes:
            - `h4` is seeded as two-phase (`Q=0.2`) to help valve/evaporator start.
            - `h1` and `h2` are seeded as superheated vapor points.
            - `h3` is seeded as a subcooled liquid point.
        """
        TL = low_side_temperature + C_TO_K
        TH = high_side_temperature + C_TO_K
        superheat = 3.0
        subcool = 3.0

        p_low = CP.PropsSI("P", "T", TL, "Q", 0, self.cp_fluid_name)
        p_high = CP.PropsSI("P", "T", TH, "Q", 0, self.cp_fluid_name)

        h4 = CP.PropsSI("H", "T", TL, "Q", 0.2, self.cp_fluid_name)
        h1 = CP.PropsSI("H", "T", TL + superheat, "Q", 1, self.cp_fluid_name)
        h2 = CP.PropsSI("H", "T", TH + superheat, "Q", 1, self.cp_fluid_name)
        h3 = CP.PropsSI("H", "T", TH - subcool, "Q", 0, self.cp_fluid_name)

        self._init = {
            "p_low": p_low,
            "p_high": p_high,
            "h1": h1,
            "h2": h2,
            "h3": h3,
            "h4": h4,
        }

    def initialize(self, verbose=False):
        """Initialize units sequentially with propagated states.

        Strategy:
            1. Fix a consistent refrigerant loop state guess.
            2. Initialize evaporator with an initial duty guess.
            3. Propagate state to compressor; initialize compressor.
            4. Propagate state to condenser; initialize condenser with duty guess.
            5. Propagate to expansion valve; initialize valve.

        Args:
            verbose: If True, print a full flowsheet report after initialization.
        """
        fs = self.model.fs
        init = getattr(self, "_init", None)
        if init is None:
            self.specify_initial_conditions(-20.0, 30.0)
            init = self._init

        fs.evaporator.cold_side_inlet.flow_mass[0].fix(1.0)
        fs.evaporator.cold_side_inlet.pressure[0].fix(init["p_low"])
        fs.evaporator.cold_side_inlet.enth_mass[0].fix(init["h4"])
        fs.evaporator.cold_side_outlet.enth_mass[0].fix(init["h1"])

        fs.compressor.inlet.pressure[0].fix(init["p_low"])
        fs.compressor.inlet.enth_mass[0].fix(init["h1"])
        fs.compressor.outlet.pressure[0].fix(init["p_high"])
        fs.compressor.efficiency_isentropic[0].fix(self.compressor_efficiency)

        fs.condenser.hot_side_inlet.pressure[0].fix(init["p_high"])
        fs.condenser.hot_side_inlet.enth_mass[0].fix(init["h2"])
        fs.condenser.hot_side_outlet.enth_mass[0].fix(init["h3"])

        fs.expansion_valve.inlet.pressure[0].fix(init["p_high"])
        fs.expansion_valve.outlet.pressure[0].fix(init["p_low"])

        # Stepwise initializer calls help LC-HX startup.
        try:
            fs.evaporator.initialize(outlvl=logging.WARNING)
        except Exception:
            pass
        try:
            propagate_state(fs.evaporator_to_compressor)
        except Exception:
            pass

        try:
            fs.compressor.initialize(outlvl=logging.WARNING)
        except Exception:
            pass
        try:
            propagate_state(fs.compressor_to_condenser)
        except Exception:
            pass

        try:
            fs.condenser.initialize(outlvl=logging.WARNING)
        except Exception:
            pass
        try:
            propagate_state(fs.condenser_to_expansion_valve)
        except Exception:
            pass

        try:
            fs.expansion_valve.initialize(outlvl=logging.WARNING)
        except Exception:
            pass

        try:
            propagate_state(fs.expansion_valve_to_evaporator)
        except Exception:
            pass

        if verbose:
            fs.report()

    def set_specifications(
        self,
        low_side_pressure=(200, 500),
        high_side_pressure=(1000, 4000),
        evaporator_temperature=(-20, 0),
        cold_storage_setpoint=None,
        evap_offset_bounds=(-35.0, 0.0),
        condenser_temperature=(30, 50),
        subcooling=3,
        superheating=3,
        max_pressure_ratio=4,
        ambient_temperature=None,
        condenser_approach=None,
        evap_sat_temperature=None,
        debug_disable_arc_pressure_eq=False,
        plr=None,
        cd=None,
        UA_evap_hot=None,
        UA_evap_cold=None,
        UA_cond_hot=None,
        UA_cond_cold=None,
    ):
        """Apply operating specifications for steady LC-HX solve.

        Args mirror the existing PLR interface for compatibility with scripts.

        Key behavior:
            - Refrigerant reference flow is fixed to 1.0 (normalized basis).
            - PLR/CD are retained for post-correction only.
            - Condenser air inlet can be tied to ambient.
            - Evaporator air inlet can be tied to cold-storage setpoint.
            - Optional approach relation links condenser saturation and ambient:
              T_cond,sat = T_ambient + approach

        Temperature bounds:
            `cold_storage_setpoint` with `evap_offset_bounds` overrides explicit
            `evaporator_temperature` bounds.
        """
        fs = self.model.fs

        if plr is not None:
            assert 0.0 <= plr <= 1.0
            self.plr = plr
        if cd is not None:
            assert 0.0 <= cd <= 1.0
            self.cd = cd
        fs.plr.set_value(self.plr)
        fs.cd.set_value(self.cd)

        # Update UA settings when provided.
        if UA_evap_hot is not None:
            fs.ua_evap_hot.set_value(float(UA_evap_hot))
        if UA_evap_cold is not None:
            fs.ua_evap_cold.set_value(float(UA_evap_cold))
        if UA_cond_hot is not None:
            fs.ua_cond_hot.set_value(float(UA_cond_hot))
        if UA_cond_cold is not None:
            fs.ua_cond_cold.set_value(float(UA_cond_cold))
        for t in fs.time:
            fs.evaporator.ua_hot_side[t].fix(value(fs.ua_evap_hot))
            fs.evaporator.ua_cold_side[t].fix(value(fs.ua_evap_cold))
            fs.condenser.ua_hot_side[t].fix(value(fs.ua_cond_hot))
            fs.condenser.ua_cold_side[t].fix(value(fs.ua_cond_cold))

        # Release initialization anchors.
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

        # Reference flow remains fixed; PLR remains post-correction.
        fs.evaporator.cold_side_inlet.flow_mass[0].fix(1.0)

        # Pressure bounds/targets.
        p_low_min, p_low_max = low_side_pressure
        p_high_min, p_high_max = high_side_pressure

        fs.P_low.setlb(p_low_min * 1e3)
        fs.P_low.setub(p_low_max * 1e3)
        fs.P_high.setlb(p_high_min * 1e3)
        fs.P_high.setub(p_high_max * 1e3)

        fs.P_low_target.set_value((p_low_min + p_low_max) * 0.5 * 1e3)
        fs.P_high_target.set_value((p_high_min + p_high_max) * 0.5 * 1e3)
        fs.P_low_target_constraint.activate()
        fs.P_high_target_constraint.activate()

        # Ambient coupling: condenser air inlet tracks ambient.
        if ambient_temperature is not None:
            fs.condenser.cold_side_inlet.temperature[0].fix(ambient_temperature + C_TO_K)
            fs.ambient_T.set_value(ambient_temperature + C_TO_K)

        fs.approach_constraint.deactivate()

        # Evaporator air side at cold storage setpoint convention.
        if cold_storage_setpoint is not None:
            fs.evaporator.hot_side_inlet.temperature[0].fix(cold_storage_setpoint + C_TO_K)
            evap_lo = cold_storage_setpoint + evap_offset_bounds[0]
            evap_hi = cold_storage_setpoint + evap_offset_bounds[1]
        else:
            fs.evaporator.hot_side_inlet.temperature[0].fix(C_TO_K - 20.0)
            evap_lo, evap_hi = evaporator_temperature

        fs.max_pressure_ratio.set_value(max_pressure_ratio)
        fs.compressor.efficiency_isentropic[0].fix(self.compressor_efficiency)

        fs.evap_Tmin.set_value(evap_lo + C_TO_K)
        fs.evap_Tmax.set_value(evap_hi + C_TO_K)
        fs.cond_Tmin.set_value(condenser_temperature[0] + C_TO_K)
        fs.cond_Tmax.set_value(condenser_temperature[1] + C_TO_K)

        # Lumped-equivalent closure: midpoint sat temperatures determine pressures.
        if evap_sat_temperature is None:
            evap_sat_temperature = 0.5 * (evap_lo + evap_hi)
        t_cond_sat_c = 0.5 * (float(condenser_temperature[0]) + float(condenser_temperature[1]))
        self._t_evap_sat_c = float(evap_sat_temperature)
        self._t_cond_sat_c = float(t_cond_sat_c)
        p_low_sat = CP.PropsSI("P", "T", self._t_evap_sat_c + C_TO_K, "Q", 1, self.cp_fluid_name)
        p_high_sat = CP.PropsSI("P", "T", self._t_cond_sat_c + C_TO_K, "Q", 0, self.cp_fluid_name)
        fs.P_low_target.set_value(float(p_low_sat))
        fs.P_high_target.set_value(float(p_high_sat))
        fs.P_low_target_constraint.activate()
        fs.P_high_target_constraint.activate()

        calculate_scaling_factors(self.model)
        set_scaling_factor(fs.evaporator.area, 1e-1)
        set_scaling_factor(fs.condenser.area, 1e-2)
        set_scaling_factor(fs.evaporator.overall_heat_transfer_coefficient[0], 1e-3)
        set_scaling_factor(fs.condenser.overall_heat_transfer_coefficient[0], 1e-3)
        set_scaling_factor(fs.evaporator.hot_side.heat[0], 1e-4)
        set_scaling_factor(fs.evaporator.cold_side.heat[0], 1e-4)
        set_scaling_factor(fs.condenser.hot_side.heat[0], 1e-4)
        set_scaling_factor(fs.condenser.cold_side.heat[0], 1e-4)
        set_scaling_factor(fs.compressor.control_volume.work[0], 1e-4)
        set_scaling_factor(fs.expansion_valve.control_volume.work[0], 1e-4)

    def optimize_COP(self, verbose=False, initialize=True, optimize=False):
        """Solve the steady cycle and return full-load COP.

        Args:
            verbose: If True, pass solver output to stdout.
            initialize: If True, run a pre-solve feasibility step.
            optimize: If True, activate COP objective/constraint for optimization.

        Returns:
            Tuple `(cop_full, converged)` where:
            - `cop_full` is `nan` when not converged.
            - `converged` is True only for optimal termination.

        Notes:
            COP is recomputed from solved state variables after the final solve.
        """
        solver = get_solver()
        solver.options = {
            "max_iter": 1200,
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

        fs = self.model.fs
        if fs.compressor.work_mechanical[0].value not in (None, 0):
            q_evap = fs.evaporator.cold_side_inlet.flow_mass[0].value * (
                fs.evaporator.cold_side_outlet.enth_mass[0].value
                - fs.evaporator.cold_side_inlet.enth_mass[0].value
            )
            fs.cop.set_value(q_evap / fs.compressor.work_mechanical[0].value)

        converged = results.solver.termination_condition == TerminationCondition.optimal
        self.optimization_converged = converged

        if not converged:
            pass

        try:
            tev = value(fs.evaporator.cold_side.properties_out[0].temperature)
            tev_sat = value(fs.evaporator.cold_side.properties_out[0].temperature_sat)
            tcd = value(fs.condenser.hot_side.properties_out[0].temperature)
            tcd_sat = value(fs.condenser.hot_side.properties_out[0].temperature_sat)
            self._last_sh_actual = max(0.0, tev - tev_sat)
            self._last_sc_actual = max(0.0, tcd_sat - tcd)
        except Exception:
            self._last_sh_actual = float("nan")
            self._last_sc_actual = float("nan")

        return value(fs.cop) if converged else float("nan"), converged

    def get_full_load_cop(self):
        """Return the solved full-load COP (`Q_evap / W_comp`)."""
        return value(self.model.fs.cop)

    def get_part_load_cop(self):
        """Return part-load corrected COP (`PLF * COP_full`)."""
        plf = self._compute_plf(value(self.model.fs.plr), value(self.model.fs.cd))
        return plf * value(self.model.fs.cop)

    def get_actual_sh_sc(self):
        """Return last solved SH/SC in K."""
        return getattr(self, "_last_sh_actual", float("nan")), getattr(self, "_last_sc_actual", float("nan"))
