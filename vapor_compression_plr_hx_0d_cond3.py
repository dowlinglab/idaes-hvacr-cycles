"""
PLR + 3-Zone Condenser HeatExchanger Model (Copy)

Author: Shilpa Narasimhan
Codex Support: OpenAI Codex
QA/Testing: Shilpa Narasimhan

Description:
    Copy-only vapor-compression refrigeration cycle model that preserves the
    PLR post-correction workflow while replacing the single condenser with a
    three-block 0-D HeatExchanger train:
    desuperheater -> condenser (two-phase) -> subcooler.

Context Breadcrumb:
    This variant keeps the same compressor/expansion and PLR logic as the PLR
    baseline, but applies zone-like condenser physics with standard IDAES
    HeatExchanger (0-D) blocks to better emulate zoned non-IDAES HX behavior.

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

from idaes.core import EnergyBalanceType, FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from pyomo.opt import TerminationCondition
from idaes.models.properties.general_helmholtz import (
    AmountBasis,
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
)
from idaes.models.unit_models import (
    Compressor,
    HeatExchanger,
    HeatExchangerFlowPattern,
    PressureChanger,
)
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_smooth_callback
from idaes.models_extra.power_generation.properties.flue_gas_ideal import FlueGasParameterBlock
from cycle_dx_hx_zonedUA.config import CycleConfig as ZonedCycleConfig
from cycle_dx_hx_zonedUA.cycle_model import solve_cycle_point as solve_zoned_cycle_point


C_TO_K = 273.15


class Mode(Enum):
    """Supported refrigerant state variable modes."""

    PH = "PH"
    IMPROVED_TPX = "improved_TPx"


class SimpleVaporCompressionCyclePLRHX0DCond3:
    """PLR cycle with 0-D evaporator and 3-zone condenser train.

    Mathematical notes:
        COP_full = Q_evap / W_comp
        PLF = 1 - CD * (1 - PLR)
        COP_part = PLF * COP_full

    Condenser train:
        compressor -> desuperheater -> condenser_tp -> subcooler -> valve

    The middle condenser block uses Underwood delta-T callback for smoother
    pinch behavior during latent heat exchange.
    """

    def __init__(
        self,
        fluid_name,
        compressor_efficiency=0.75,
        PLR=0.75,
        CD=0.13,
        mode=Mode.PH,
        evap_area_m2=6.0,
        cond_ds_area_m2=30.0,
        cond_tp_area_m2=40.0,
        cond_sc_area_m2=30.0,
        UA_evap_hot=500.0,
        UA_evap_cold=1000.0,
        UA_cond_ds_hot=450.0,
        UA_cond_ds_cold=450.0,
        UA_cond_tp_hot=900.0,
        UA_cond_tp_cold=900.0,
        UA_cond_sc_hot=350.0,
        UA_cond_sc_cold=350.0,
        evap_air_mol_s=31.0,
        cond_air_mol_s=356.0,
    ):
        """Construct the 0-D HX PLR cycle copy with zoned condenser blocks."""
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

        fs = self.model.fs
        fs.ref_props = HelmholtzParameterBlock(
            pure_component=self.idaes_fluid_name,
            state_vars=StateVars.PH,
            amount_basis=AmountBasis.MASS,
            phase_presentation=PhaseType.LG,
        )
        fs.air_props = FlueGasParameterBlock()

        fs.evap_area = Param(initialize=evap_area_m2, mutable=True, units=pyunits.m**2)
        fs.ua_evap_hot = Param(initialize=UA_evap_hot, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_evap_cold = Param(initialize=UA_evap_cold, mutable=True, units=pyunits.W / pyunits.K)
        fs.evap_air_mol_s = Param(initialize=evap_air_mol_s, mutable=True, units=pyunits.mol / pyunits.s)

        fs.cond_ds_area = Param(initialize=cond_ds_area_m2, mutable=True, units=pyunits.m**2)
        fs.cond_tp_area = Param(initialize=cond_tp_area_m2, mutable=True, units=pyunits.m**2)
        fs.cond_sc_area = Param(initialize=cond_sc_area_m2, mutable=True, units=pyunits.m**2)
        fs.cond_total_area = Param(
            initialize=(cond_ds_area_m2 + cond_tp_area_m2 + cond_sc_area_m2),
            mutable=True,
            units=pyunits.m**2,
        )

        fs.ua_cond_ds_hot = Param(initialize=UA_cond_ds_hot, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_cond_ds_cold = Param(initialize=UA_cond_ds_cold, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_cond_tp_hot = Param(initialize=UA_cond_tp_hot, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_cond_tp_cold = Param(initialize=UA_cond_tp_cold, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_cond_sc_hot = Param(initialize=UA_cond_sc_hot, mutable=True, units=pyunits.W / pyunits.K)
        fs.ua_cond_sc_cold = Param(initialize=UA_cond_sc_cold, mutable=True, units=pyunits.W / pyunits.K)
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

        Formula:
            PLF = 1 - CD * (1 - PLR)
        """
        if not (0.0 <= plr <= 1.0):
            raise ValueError("PLR must be in [0, 1]")
        if not (0.0 <= cd <= 1.0):
            raise ValueError("CD must be in [0, 1]")
        return max(0.0, min(1.0, 1.0 - cd * (1.0 - plr)))

    @staticmethod
    def _ua_to_u(ua_hot, ua_cold, area):
        """Convert hot/cold film UA pair to overall U using resistance addition."""
        ua_tot = 1.0 / (1.0 / max(ua_hot, 1e-9) + 1.0 / max(ua_cold, 1e-9))
        return ua_tot / max(area, 1e-9)

    def _define_flowsheet(self):
        """Build steady-state cycle with a 3-zone 0-D condenser train."""
        fs = self.model.fs

        fs.evaporator = HeatExchanger(
            hot_side={
                "property_package": fs.air_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            cold_side={
                "property_package": fs.ref_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            flow_pattern=HeatExchangerFlowPattern.cocurrent,
            delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
        )

        fs.compressor = Compressor(property_package=fs.ref_props)

        fs.desuperheater = HeatExchanger(
            hot_side={
                "property_package": fs.ref_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            cold_side={
                "property_package": fs.air_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            flow_pattern=HeatExchangerFlowPattern.cocurrent,
            delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
        )

        fs.condenser = HeatExchanger(
            hot_side={
                "property_package": fs.ref_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            cold_side={
                "property_package": fs.air_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            flow_pattern=HeatExchangerFlowPattern.cocurrent,
            delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
        )

        fs.subcooler = HeatExchanger(
            hot_side={
                "property_package": fs.ref_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            cold_side={
                "property_package": fs.air_props,
                "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                "has_pressure_change": False,
            },
            flow_pattern=HeatExchangerFlowPattern.cocurrent,
            delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
        )

        fs.expansion_valve = PressureChanger(
            property_package=fs.ref_props,
            thermodynamic_assumption="adiabatic",
            compressor=False,
        )

        fs.evaporator_to_compressor = Arc(source=fs.evaporator.cold_side_outlet, destination=fs.compressor.inlet)
        fs.compressor_to_desuperheater = Arc(source=fs.compressor.outlet, destination=fs.desuperheater.hot_side_inlet)
        fs.desuperheater_to_condenser = Arc(source=fs.desuperheater.hot_side_outlet, destination=fs.condenser.hot_side_inlet)
        fs.condenser_to_subcooler = Arc(source=fs.condenser.hot_side_outlet, destination=fs.subcooler.hot_side_inlet)
        fs.subcooler_to_expansion_valve = Arc(source=fs.subcooler.hot_side_outlet, destination=fs.expansion_valve.inlet)
        fs.expansion_valve_to_evaporator = Arc(source=fs.expansion_valve.outlet, destination=fs.evaporator.cold_side_inlet)
        # Serial condenser-train air path: SC -> COND -> DS.
        fs.air_subcooler_to_condenser = Arc(
            source=fs.subcooler.cold_side_outlet,
            destination=fs.condenser.cold_side_inlet,
        )
        fs.air_condenser_to_desuperheater = Arc(
            source=fs.condenser.cold_side_outlet,
            destination=fs.desuperheater.cold_side_inlet,
        )

        TransformationFactory("network.expand_arcs").apply_to(self.model)

        for cname in ["flow_mass_equality", "flow_mol_equality"]:
            c = getattr(fs.evaporator_to_compressor_expanded, cname, None)
            if c is not None:
                c.deactivate()

        for arc in [
            fs.evaporator_to_compressor_expanded,
            fs.compressor_to_desuperheater_expanded,
            fs.desuperheater_to_condenser_expanded,
            fs.condenser_to_subcooler_expanded,
            fs.subcooler_to_expansion_valve_expanded,
            fs.expansion_valve_to_evaporator_expanded,
        ]:
            if hasattr(arc, "pressure_equality"):
                arc.pressure_equality.deactivate()

        self._fix_hx_0d_parameters()
        self._fix_air_side_defaults()
        self._build_common_cycle_constraints()

    def _fix_hx_0d_parameters(self):
        """Fix HX geometry and U values from configured UA values."""
        fs = self.model.fs

        fs.evaporator.area.fix(value(fs.evap_area))
        fs.desuperheater.area.setlb(1.0)
        fs.condenser.area.setlb(1.0)
        fs.subcooler.area.setlb(1.0)
        # Keep zone areas fixed for closure in the base solve path.
        fs.desuperheater.area.fix(value(fs.cond_ds_area))
        fs.condenser.area.fix(value(fs.cond_tp_area))
        fs.subcooler.area.fix(value(fs.cond_sc_area))
        if hasattr(fs, "cond_area_sum"):
            fs.cond_area_sum.deactivate()

        fs.evaporator.overall_heat_transfer_coefficient[0].fix(
            self._ua_to_u(value(fs.ua_evap_hot), value(fs.ua_evap_cold), value(fs.evap_area))
        )
        fs.desuperheater.overall_heat_transfer_coefficient[0].fix(
            self._ua_to_u(value(fs.ua_cond_ds_hot), value(fs.ua_cond_ds_cold), value(fs.cond_ds_area))
        )
        fs.condenser.overall_heat_transfer_coefficient[0].fix(
            self._ua_to_u(value(fs.ua_cond_tp_hot), value(fs.ua_cond_tp_cold), value(fs.cond_tp_area))
        )
        fs.subcooler.overall_heat_transfer_coefficient[0].fix(
            self._ua_to_u(value(fs.ua_cond_sc_hot), value(fs.ua_cond_sc_cold), value(fs.cond_sc_area))
        )

    def _fix_air_side_defaults(self):
        """Fix air-side flow composition, pressures, and nominal temperatures."""
        fs = self.model.fs

        comp = {"H2O": 0.0, "CO2": 0.0004, "N2": 0.79, "O2": 0.2096, "NO": 0.0, "SO2": 0.0}

        evap_total = value(fs.evap_air_mol_s)
        cond_total = value(fs.cond_air_mol_s)

        for j, x in comp.items():
            fs.evaporator.hot_side_inlet.flow_mol_comp[0, j].fix(evap_total * x)
            # Condenser-train air inlet is fixed only on first block (subcooler).
            fs.subcooler.cold_side_inlet.flow_mol_comp[0, j].fix(cond_total * x)

        fs.evaporator.hot_side_inlet.pressure[0].fix(101325.0)
        fs.subcooler.cold_side_inlet.pressure[0].fix(101325.0)
        fs.condenser.cold_side_inlet.pressure[0].unfix()
        fs.desuperheater.cold_side_inlet.pressure[0].unfix()

        fs.evaporator.hot_side_inlet.temperature[0].fix(C_TO_K - 20.0)
        fs.subcooler.cold_side_inlet.temperature[0].fix(C_TO_K + 20.0)
        fs.condenser.cold_side_inlet.temperature[0].unfix()
        fs.desuperheater.cold_side_inlet.temperature[0].unfix()

    def _build_common_cycle_constraints(self):
        """Add pressure-level closure, operating bounds, and COP relation."""
        fs = self.model.fs

        fs.P_low = Var(initialize=200e3, bounds=(50e3, 2000e3), units=pyunits.Pa)
        fs.P_high = Var(initialize=1200e3, bounds=(100e3, 5000e3), units=pyunits.Pa)

        fs.P_low_target = Param(initialize=200e3, mutable=True, units=pyunits.Pa)
        fs.P_high_target = Param(initialize=1200e3, mutable=True, units=pyunits.Pa)

        fs.P_low_target_constraint = Constraint(expr=fs.P_low == fs.P_low_target)
        fs.P_high_target_constraint = Constraint(expr=fs.P_high == fs.P_high_target)

        fs.P_low_evap_in = Constraint(expr=fs.evaporator.cold_side_inlet.pressure[0] == fs.expansion_valve.outlet.pressure[0])
        fs.P_low_comp_in = Constraint(expr=fs.compressor.inlet.pressure[0] == fs.P_low)
        fs.P_low_valve_out = Constraint(expr=fs.expansion_valve.outlet.pressure[0] == fs.P_low)

        fs.P_high_comp_out = Constraint(expr=fs.compressor.outlet.pressure[0] == fs.P_high)
        fs.P_high_ds_in = Constraint(expr=fs.desuperheater.hot_side_inlet.pressure[0] == fs.P_high)
        fs.P_high_cond_in = Constraint(expr=fs.condenser.hot_side_inlet.pressure[0] == fs.P_high)
        fs.P_high_sc_in = Constraint(expr=fs.subcooler.hot_side_inlet.pressure[0] == fs.P_high)
        fs.P_high_valve_in = Constraint(expr=fs.expansion_valve.inlet.pressure[0] == fs.P_high)
        fs.cond_train_p12 = Constraint(
            expr=fs.desuperheater.hot_side_inlet.pressure[0] == fs.condenser.hot_side_inlet.pressure[0]
        )
        fs.cond_train_p23 = Constraint(
            expr=fs.condenser.hot_side_inlet.pressure[0] == fs.subcooler.hot_side_inlet.pressure[0]
        )
        fs.cond_train_p01 = Constraint(
            expr=fs.compressor.outlet.pressure[0] == fs.desuperheater.hot_side_inlet.pressure[0]
        )
        fs.subcooler_hot_dp0 = Constraint(
            expr=fs.subcooler.hot_side_outlet.pressure[0] == fs.subcooler.hot_side_inlet.pressure[0]
        )

        fs.max_pressure_ratio = Param(initialize=20.0, mutable=True, units=pyunits.dimensionless)
        fs.pressure_ratio_constraint = Constraint(
            expr=fs.compressor.outlet.pressure[0] <= fs.max_pressure_ratio * fs.compressor.inlet.pressure[0]
        )

        fs.ambient_T = Param(initialize=C_TO_K + 20.0, mutable=True, units=pyunits.K)
        fs.approach_T = Param(initialize=10.0, mutable=True, units=pyunits.K)
        fs.approach_constraint = Constraint(
            expr=fs.subcooler.hot_side.properties_out[0].temperature_sat == fs.ambient_T + fs.approach_T
        )
        fs.approach_constraint.deactivate()

        # Small pinch guards to keep iterative states out of temperature-crossing regimes.
        fs.pinch_margin = Param(initialize=0.5, mutable=True, units=pyunits.K)
        fs.evap_pinch_guard = Constraint(
            expr=fs.evaporator.hot_side_inlet.temperature[0]
            - fs.evaporator.cold_side.properties_out[0].temperature
            >= fs.pinch_margin
        )
        fs.cond_pinch_guard = Constraint(
            expr=fs.desuperheater.hot_side.properties_out[0].temperature
            - fs.desuperheater.cold_side_inlet.temperature[0]
            >= fs.pinch_margin
        )

        fs.evap_Tmin = Param(initialize=C_TO_K - 60.0, mutable=True, units=pyunits.K)
        fs.evap_Tmax = Param(initialize=C_TO_K - 10.0, mutable=True, units=pyunits.K)
        fs.cond_Tmin = Param(initialize=C_TO_K + 15.0, mutable=True, units=pyunits.K)
        fs.cond_Tmax = Param(initialize=C_TO_K + 70.0, mutable=True, units=pyunits.K)

        fs.evap_T_lower = Constraint(expr=fs.evaporator.cold_side.properties_out[0].temperature >= fs.evap_Tmin)
        fs.evap_T_upper = Constraint(expr=fs.evaporator.cold_side.properties_out[0].temperature <= fs.evap_Tmax)
        fs.cond_T_lower = Constraint(expr=fs.subcooler.hot_side.properties_out[0].temperature >= fs.cond_Tmin)
        fs.cond_T_upper = Constraint(expr=fs.subcooler.hot_side.properties_out[0].temperature <= fs.cond_Tmax)
        fs.subcooler_temperature_drop = Constraint(
            expr=fs.subcooler.hot_side.properties_out[0].temperature
            <= fs.subcooler.hot_side.properties_in[0].temperature - 0.5 * pyunits.K
        )
        # Air-side monotonic heating guard through each condenser-train block.
        fs.ds_air_heating_guard = Constraint(
            expr=fs.desuperheater.cold_side.properties_out[0].temperature
            >= fs.desuperheater.cold_side.properties_in[0].temperature + 0.1 * pyunits.K
        )
        fs.cond_air_heating_guard = Constraint(
            expr=fs.condenser.cold_side.properties_out[0].temperature
            >= fs.condenser.cold_side.properties_in[0].temperature + 0.1 * pyunits.K
        )
        fs.sc_air_heating_guard = Constraint(
            expr=fs.subcooler.cold_side.properties_out[0].temperature
            >= fs.subcooler.cold_side.properties_in[0].temperature + 0.1 * pyunits.K
        )
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
        fs.cond_area_sum = Constraint(
            expr=fs.desuperheater.area + fs.condenser.area + fs.subcooler.area == fs.cond_total_area
        )

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
        """Generate pressure/enthalpy seeds from saturation anchors.

        The condenser train is seeded at three key hot-side states:
            h2  : compressor discharge superheated vapor
            hg  : saturated vapor at high pressure
            hf  : saturated liquid at high pressure
            h3  : subcooled liquid
        """
        self._init_low_side_temperature_c = float(low_side_temperature)
        self._init_high_side_temperature_c = float(high_side_temperature)

        TL = low_side_temperature + C_TO_K
        TH = high_side_temperature + C_TO_K
        superheat = 3.0
        subcool = 3.0

        p_low = CP.PropsSI("P", "T", TL, "Q", 0, self.cp_fluid_name)
        p_high = CP.PropsSI("P", "T", TH, "Q", 0, self.cp_fluid_name)

        h4 = CP.PropsSI("H", "T", TL, "Q", 0.2, self.cp_fluid_name)
        h1 = CP.PropsSI("H", "T", TL + superheat, "Q", 1, self.cp_fluid_name)
        h2 = CP.PropsSI("H", "T", TH + superheat, "Q", 1, self.cp_fluid_name)
        h_g = CP.PropsSI("H", "T", TH, "Q", 1, self.cp_fluid_name)
        h_f = CP.PropsSI("H", "T", TH, "Q", 0, self.cp_fluid_name)
        h3 = CP.PropsSI("H", "T", TH - subcool, "Q", 0, self.cp_fluid_name)

        self._init = {
            "p_low": p_low,
            "p_high": p_high,
            "h1": h1,
            "h2": h2,
            "h_g": h_g,
            "h_f": h_f,
            "h3": h3,
            "h4": h4,
        }

    def initialize(self, verbose=False):
        """Initialize blocks sequentially with propagated states."""
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

        fs.desuperheater.hot_side_inlet.pressure[0].fix(init["p_high"])
        fs.desuperheater.hot_side_inlet.enth_mass[0].fix(init["h2"])
        fs.desuperheater.hot_side_outlet.enth_mass[0].setlb(init["h_f"])
        fs.desuperheater.hot_side_outlet.enth_mass[0].setub(init["h_g"])
        fs.desuperheater.hot_side_outlet.enth_mass[0].set_value(init["h_g"])

        fs.condenser.hot_side_inlet.pressure[0].fix(init["p_high"])
        # Keep refrigerant seeds continuous across DS -> COND.
        fs.condenser.hot_side_inlet.enth_mass[0].fix(init["h_g"])
        fs.condenser.hot_side_outlet.enth_mass[0].setlb(init["h_f"])
        fs.condenser.hot_side_outlet.enth_mass[0].setub(init["h_f"] + 5000.0)
        fs.condenser.hot_side_outlet.enth_mass[0].set_value(init["h_f"])

        fs.subcooler.hot_side_inlet.pressure[0].fix(init["p_high"])
        # Keep refrigerant seeds continuous across COND -> SC.
        fs.subcooler.hot_side_inlet.enth_mass[0].fix(init["h_f"])
        fs.subcooler.hot_side_outlet.enth_mass[0].setub(init["h_f"])
        fs.subcooler.hot_side_outlet.enth_mass[0].set_value(init["h3"])

        fs.expansion_valve.inlet.pressure[0].fix(init["p_high"])
        fs.expansion_valve.outlet.pressure[0].fix(init["p_low"])

        # Air-side enthalpy seed to keep flue-gas states in an ambient-like region.
        for st in [
            fs.subcooler.cold_side.properties_in[0],
            fs.subcooler.cold_side.properties_out[0],
            fs.condenser.cold_side.properties_in[0],
            fs.condenser.cold_side.properties_out[0],
            fs.desuperheater.cold_side.properties_in[0],
            fs.desuperheater.cold_side.properties_out[0],
        ]:
            try:
                st.enth_mol.set_value(-250.0)
            except Exception:
                pass

        deactivated_arc_constraints = []
        for arc in [
            fs.evaporator_to_compressor_expanded,
            fs.compressor_to_desuperheater_expanded,
            fs.desuperheater_to_condenser_expanded,
            fs.condenser_to_subcooler_expanded,
            fs.subcooler_to_expansion_valve_expanded,
            fs.expansion_valve_to_evaporator_expanded,
        ]:
            for cname in [
                "flow_mass_equality",
                "flow_mol_equality",
                "enth_mass_equality",
                "enth_mol_equality",
                "pressure_equality",
            ]:
                c = getattr(arc, cname, None)
                if c is not None and c.active:
                    c.deactivate()
                    deactivated_arc_constraints.append(c)

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
            propagate_state(fs.compressor_to_desuperheater)
        except Exception:
            pass

        try:
            fs.desuperheater.initialize(outlvl=logging.WARNING)
        except Exception:
            pass
        try:
            propagate_state(fs.desuperheater_to_condenser)
        except Exception:
            pass

        try:
            fs.condenser.initialize(outlvl=logging.WARNING)
        except Exception:
            pass
        try:
            propagate_state(fs.condenser_to_subcooler)
        except Exception:
            pass

        try:
            fs.subcooler.initialize(outlvl=logging.WARNING)
        except Exception:
            pass
        try:
            propagate_state(fs.subcooler_to_expansion_valve)
        except Exception:
            pass

        # Solve condenser train as a disconnected line first.
        try:
            solver = get_solver()
            solver.options = {"max_iter": 300, "tol": 1e-6, "acceptable_tol": 1e-5}
            solver.solve(self.model, tee=False)
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

        # Then solve with evaporator/valve included, still with arcs open.
        try:
            solver = get_solver()
            solver.options = {"max_iter": 300, "tol": 1e-6, "acceptable_tol": 1e-5}
            solver.solve(self.model, tee=False)
        except Exception:
            pass

        # Restore full loop closure.
        for c in deactivated_arc_constraints:
            c.activate()

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
        UA_cond_ds_hot=None,
        UA_cond_ds_cold=None,
        UA_cond_tp_hot=None,
        UA_cond_tp_cold=None,
        UA_cond_sc_hot=None,
        UA_cond_sc_cold=None,
    ):
        """Apply operating specifications for the 3-zone condenser cycle.

        Rules used to mirror the non-IDAES lumped operating policy:
            - evaporator air inlet tracks cold-storage setpoint when provided.
            - condenser train uses one air stream:
              subcooler air inlet = ambient,
              then subcooler outlet feeds condenser inlet,
              and condenser outlet feeds desuperheater inlet.
            - evap/cond saturation pressures are closed from sat temperature
              midpoints (or explicit evap sat target when provided).
        """
        _ = (subcooling, superheating, condenser_approach, debug_disable_arc_pressure_eq)

        fs = self.model.fs

        if plr is not None:
            assert 0.0 <= plr <= 1.0
            self.plr = plr
        if cd is not None:
            assert 0.0 <= cd <= 1.0
            self.cd = cd
        fs.plr.set_value(self.plr)
        fs.cd.set_value(self.cd)

        if UA_evap_hot is not None:
            fs.ua_evap_hot.set_value(float(UA_evap_hot))
        if UA_evap_cold is not None:
            fs.ua_evap_cold.set_value(float(UA_evap_cold))

        # Optional shorthand applies same condenser-side UA to all three zones.
        if UA_cond_hot is not None:
            fs.ua_cond_ds_hot.set_value(float(UA_cond_hot))
            fs.ua_cond_tp_hot.set_value(float(UA_cond_hot))
            fs.ua_cond_sc_hot.set_value(float(UA_cond_hot))
        if UA_cond_cold is not None:
            fs.ua_cond_ds_cold.set_value(float(UA_cond_cold))
            fs.ua_cond_tp_cold.set_value(float(UA_cond_cold))
            fs.ua_cond_sc_cold.set_value(float(UA_cond_cold))

        if UA_cond_ds_hot is not None:
            fs.ua_cond_ds_hot.set_value(float(UA_cond_ds_hot))
        if UA_cond_ds_cold is not None:
            fs.ua_cond_ds_cold.set_value(float(UA_cond_ds_cold))
        if UA_cond_tp_hot is not None:
            fs.ua_cond_tp_hot.set_value(float(UA_cond_tp_hot))
        if UA_cond_tp_cold is not None:
            fs.ua_cond_tp_cold.set_value(float(UA_cond_tp_cold))
        if UA_cond_sc_hot is not None:
            fs.ua_cond_sc_hot.set_value(float(UA_cond_sc_hot))
        if UA_cond_sc_cold is not None:
            fs.ua_cond_sc_cold.set_value(float(UA_cond_sc_cold))

        self._fix_hx_0d_parameters()

        # Release initialization anchors.
        fs.evaporator.cold_side_inlet.flow_mass[0].unfix()
        fs.evaporator.cold_side_inlet.pressure[0].unfix()
        fs.evaporator.cold_side_inlet.enth_mass[0].unfix()
        fs.evaporator.cold_side_outlet.enth_mass[0].unfix()

        fs.compressor.inlet.pressure[0].unfix()
        fs.compressor.inlet.enth_mass[0].unfix()
        fs.compressor.outlet.pressure[0].unfix()

        for hx in [fs.desuperheater, fs.condenser, fs.subcooler]:
            hx.hot_side_inlet.pressure[0].unfix()
            hx.hot_side_inlet.enth_mass[0].unfix()
            hx.hot_side_outlet.enth_mass[0].unfix()

        fs.expansion_valve.inlet.pressure[0].unfix()
        fs.expansion_valve.outlet.pressure[0].unfix()

        # Reference flow remains fixed; PLR stays post-correction only.
        fs.evaporator.cold_side_inlet.flow_mass[0].fix(1.0)

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

        if ambient_temperature is not None:
            t_amb_k = ambient_temperature + C_TO_K
            fs.ambient_T.set_value(t_amb_k)
            # Only the first condenser-train air inlet is fixed to ambient.
            fs.subcooler.cold_side_inlet.temperature[0].fix(t_amb_k)
            fs.subcooler.cold_side_inlet.pressure[0].fix(101325.0)
            # Downstream inlets are arc-driven; keep them unfixed.
            fs.condenser.cold_side_inlet.temperature[0].unfix()
            fs.condenser.cold_side_inlet.pressure[0].unfix()
            fs.desuperheater.cold_side_inlet.temperature[0].unfix()
            fs.desuperheater.cold_side_inlet.pressure[0].unfix()
            fs.condenser.cold_side_inlet.temperature[0].set_value(t_amb_k + 0.1)
            fs.desuperheater.cold_side_inlet.temperature[0].set_value(t_amb_k + 0.2)

        fs.approach_constraint.deactivate()

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

        if evap_sat_temperature is None:
            evap_sat_temperature = 0.5 * (evap_lo + evap_hi)
        t_cond_sat_c = 0.5 * (float(condenser_temperature[0]) + float(condenser_temperature[1]))
        self._t_evap_sat_c = float(evap_sat_temperature)
        self._t_cond_sat_c = float(t_cond_sat_c)
        p_low_sat = CP.PropsSI("P", "T", self._t_evap_sat_c + C_TO_K, "Q", 1, self.cp_fluid_name)
        p_high_sat = CP.PropsSI("P", "T", self._t_cond_sat_c + C_TO_K, "Q", 0, self.cp_fluid_name)
        fs.P_low_target.set_value(float(p_low_sat))
        fs.P_high_target.set_value(float(p_high_sat))
        fs.compressor.outlet.pressure[0].fix(float(p_high_sat))
        fs.P_high_target_constraint.deactivate()
        fs.P_high_ds_in.deactivate()
        fs.P_high_cond_in.deactivate()
        fs.P_high_sc_in.deactivate()
        fs.subcooler_hot_dp0.activate()
        # Avoid redundant zero-dP equations on the subcooler hot side.
        try:
            fs.subcooler.hot_side.pressure_balance[0].deactivate()
        except Exception:
            pass

        # Hard anti-crossing bounds to avoid HX temperature crossing during iterations.
        try:
            evap_air_in = value(fs.evaporator.hot_side_inlet.temperature[0])
            fs.evaporator.cold_side.properties_out[0].temperature.setub(evap_air_in - 0.1)
        except Exception:
            pass
        for hx in [fs.desuperheater, fs.condenser, fs.subcooler]:
            try:
                air_in = value(hx.cold_side_inlet.temperature[0])
                hx.hot_side.properties_out[0].temperature.setlb(air_in + 0.1)
            except Exception:
                pass

        # Air outlet temperature ceiling to keep iterates in physical ambient-cooling regime.
        try:
            fs.evaporator.hot_side.properties_out[0].temperature.setub(350.0)
        except Exception:
            pass
        for hx in [fs.desuperheater, fs.condenser, fs.subcooler]:
            try:
                hx.cold_side.properties_out[0].temperature.setub(350.0)
            except Exception:
                pass

        # Initialize air enthalpy states to avoid None values in diagnostics/solver.
        for st in [
            fs.evaporator.hot_side.properties_in[0],
            fs.evaporator.hot_side.properties_out[0],
            fs.desuperheater.cold_side.properties_in[0],
            fs.desuperheater.cold_side.properties_out[0],
            fs.condenser.cold_side.properties_in[0],
            fs.condenser.cold_side.properties_out[0],
            fs.subcooler.cold_side.properties_in[0],
            fs.subcooler.cold_side.properties_out[0],
        ]:
            try:
                st.enth_mol.set_value(-250.0)
            except Exception:
                pass

        calculate_scaling_factors(self.model)
        for hx in [fs.evaporator, fs.desuperheater, fs.condenser, fs.subcooler]:
            set_scaling_factor(hx.area, 1e-1)
            set_scaling_factor(hx.overall_heat_transfer_coefficient[0], 1e-3)
            set_scaling_factor(hx.hot_side.heat[0], 1e-4)
            set_scaling_factor(hx.cold_side.heat[0], 1e-4)
        set_scaling_factor(fs.compressor.control_volume.work[0], 1e-4)
        set_scaling_factor(fs.expansion_valve.control_volume.work[0], 1e-4)

        # Enthalpy scaling in J/kg space to reduce residual magnitude.
        def _scale_port_enthalpy(port, sf):
            if hasattr(port, "enth_mass"):
                set_scaling_factor(port.enth_mass[0], sf)
            elif hasattr(port, "enth_mol"):
                set_scaling_factor(port.enth_mol[0], sf)

        hx_blocks = [fs.evaporator, fs.desuperheater, fs.condenser, fs.subcooler]
        for hx in hx_blocks:
            _scale_port_enthalpy(hx.hot_side_inlet, 1e-5)
            _scale_port_enthalpy(hx.hot_side_outlet, 1e-5)
            _scale_port_enthalpy(hx.cold_side_inlet, 1e-5)
            _scale_port_enthalpy(hx.cold_side_outlet, 1e-5)

        # Air-side molar enthalpy scaling for condenser-train state variables.
        for hx in [fs.desuperheater, fs.condenser, fs.subcooler]:
            try:
                set_scaling_factor(hx.cold_side.properties_in[0].enth_mol, 1e-3)
            except Exception:
                pass
            try:
                set_scaling_factor(hx.cold_side.properties_out[0].enth_mol, 1e-3)
            except Exception:
                pass
        _scale_port_enthalpy(fs.compressor.inlet, 1e-5)
        _scale_port_enthalpy(fs.compressor.outlet, 1e-5)
        _scale_port_enthalpy(fs.expansion_valve.inlet, 1e-5)
        _scale_port_enthalpy(fs.expansion_valve.outlet, 1e-5)

    def optimize_COP(self, verbose=False, initialize=True, optimize=False):
        """Solve the steady cycle and return full-load COP."""
        # Keep arc expansion synchronized with any manual specification edits.
        try:
            TransformationFactory("network.expand_arcs").apply_to(self.model)
        except Exception:
            pass

        # Ensure key refrigerant loop arc equalities are active before solve.
        try:
            eq = self.model.fs.compressor_to_desuperheater_expanded.enth_mass_equality
            if not eq.active:
                eq.activate()
        except Exception:
            pass

        solver = get_solver()
        solver.options = {"max_iter": 1200, "tol": 1e-6, "acceptable_tol": 1e-5}

        if initialize:
            self._seed_from_zoned_profile()
            self._full_profile_warmstart(verbose=verbose)
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
            self._print_debug_snapshot("solver_exception")
            return float("nan"), False

        fs = self.model.fs
        if fs.compressor.work_mechanical[0].value not in (None, 0):
            q_evap = fs.evaporator.cold_side_inlet.flow_mass[0].value * (
                fs.evaporator.cold_side_outlet.enth_mass[0].value - fs.evaporator.cold_side_inlet.enth_mass[0].value
            )
            fs.cop.set_value(q_evap / fs.compressor.work_mechanical[0].value)

        converged = results.solver.termination_condition == TerminationCondition.optimal
        self.optimization_converged = converged
        if not converged:
            self._print_debug_snapshot(str(results.solver.termination_condition))

        try:
            tev = value(fs.evaporator.cold_side.properties_out[0].temperature)
            tev_sat = value(fs.evaporator.cold_side.properties_out[0].temperature_sat)
            tcd = value(fs.subcooler.hot_side.properties_out[0].temperature)
            tcd_sat = value(fs.subcooler.hot_side.properties_out[0].temperature_sat)
            self._last_sh_actual = max(0.0, tev - tev_sat)
            self._last_sc_actual = max(0.0, tcd_sat - tcd)
        except Exception:
            self._last_sh_actual = float("nan")
            self._last_sc_actual = float("nan")

        return value(fs.cop) if converged else float("nan"), converged

    def _seed_from_zoned_profile(self):
        """Seed refrigerant state profile using non-IDAES zoned model outputs.

        Uses the same ambient and saturation-temperature policy currently stored
        in this model (`self._t_evap_sat_c`, `self._t_cond_sat_c`, and inlet
        air temperatures), then maps zoned-cycle `(P, h)` states into IDAES
        stream variable initial values without fixing them.
        """
        fs = self.model.fs
        def _ua_total(ua_hot, ua_cold):
            return 1.0 / (1.0 / max(ua_hot, 1.0e-9) + 1.0 / max(ua_cold, 1.0e-9))

        try:
            ua_evap = _ua_total(value(fs.ua_evap_hot), value(fs.ua_evap_cold))
            ua_ds = _ua_total(value(fs.ua_cond_ds_hot), value(fs.ua_cond_ds_cold))
            ua_tp = _ua_total(value(fs.ua_cond_tp_hot), value(fs.ua_cond_tp_cold))
            ua_sc = _ua_total(value(fs.ua_cond_sc_hot), value(fs.ua_cond_sc_cold))
            cfg = ZonedCycleConfig(
                eta_isentropic=self.compressor_efficiency,
                m_dot_ref=1.0,  # IDAES reference-flow basis in this model
                plr=self.plr,
                cd=self.cd,
                UA_evap_total=ua_evap,
                UA_cond_total=(ua_ds + ua_tp + ua_sc),
                UA_evap_tp=0.67 * ua_evap,
                UA_evap_sh=0.33 * ua_evap,
                UA_cond_ds=ua_ds,
                UA_cond_tp=ua_tp,
                UA_cond_sc=ua_sc,
            )
            rz = solve_zoned_cycle_point(
                fluid=self.cp_fluid_name,
                t_evap_sat_c=float(getattr(self, "_init_low_side_temperature_c", getattr(self, "_t_evap_sat_c", -29.0))),
                t_cond_sat_c=float(getattr(self, "_init_high_side_temperature_c", getattr(self, "_t_cond_sat_c", 29.0))),
                cfg=cfg,
                t_air_evap_in_c=float(value(fs.evaporator.hot_side_inlet.temperature[0]) - C_TO_K),
                t_air_cond_in_c=float(value(fs.subcooler.cold_side_inlet.temperature[0]) - C_TO_K),
            )
            h1 = float(rz.diagnostics.get("h1_out", rz.diagnostics.get("h1_in")))
            h2 = float(rz.diagnostics["h2"])
            h3 = float(rz.diagnostics["h3"])
            h4 = float(rz.diagnostics["h4"])
            p_low = float(rz.p_evap_pa)
            p_high = float(rz.p_cond_pa)
            m_ref = 1.0

            q_ds = float(rz.diagnostics.get("cond_zone_q_ds", 0.0))
            q_tp = float(rz.diagnostics.get("cond_zone_q_tp", 0.0))
            q_sc = float(rz.diagnostics.get("cond_zone_q_sc", 0.0))
            h_ds_out = h2 - q_ds / m_ref
            h_cond_out = h_ds_out - q_tp / m_ref
            h_sc_out = h_cond_out - q_sc / m_ref

            t_air_sc_in = float(value(fs.subcooler.cold_side_inlet.temperature[0]))
            # Keep downstream seeds in an ambient-range band to avoid 500 K excursions.
            t_air_cond_in = t_air_sc_in + 0.1
            t_air_ds_in = t_air_sc_in + 0.2
            t_air_ds_out = t_air_sc_in + 0.3

            # Refrigerant loop seed values.
            fs.evaporator.cold_side_inlet.pressure[0].set_value(p_low)
            fs.evaporator.cold_side_inlet.enth_mass[0].set_value(h4)
            fs.evaporator.cold_side_outlet.enth_mass[0].set_value(h1)
            fs.compressor.inlet.pressure[0].set_value(p_low)
            fs.compressor.inlet.enth_mass[0].set_value(h1)
            fs.compressor.outlet.pressure[0].set_value(p_high)
            fs.desuperheater.hot_side_inlet.pressure[0].set_value(p_high)
            fs.desuperheater.hot_side_inlet.enth_mass[0].set_value(h2)
            fs.desuperheater.hot_side_outlet.enth_mass[0].set_value(h_ds_out)
            fs.condenser.hot_side_inlet.pressure[0].set_value(p_high)
            fs.condenser.hot_side_inlet.enth_mass[0].set_value(h_ds_out)
            fs.condenser.hot_side_outlet.enth_mass[0].set_value(h_cond_out)
            fs.subcooler.hot_side_inlet.pressure[0].set_value(p_high)
            fs.subcooler.hot_side_inlet.enth_mass[0].set_value(h_cond_out)
            fs.subcooler.hot_side_outlet.enth_mass[0].set_value(h_sc_out)
            fs.expansion_valve.inlet.pressure[0].set_value(p_high)
            fs.expansion_valve.outlet.pressure[0].set_value(p_low)

            # Air-side SC -> COND -> DS seeds.
            fs.subcooler.cold_side_inlet.temperature[0].set_value(t_air_sc_in)
            fs.condenser.cold_side_inlet.temperature[0].set_value(t_air_cond_in)
            fs.desuperheater.cold_side_inlet.temperature[0].set_value(t_air_ds_in)
            try:
                fs.desuperheater.cold_side_outlet.temperature[0].set_value(t_air_ds_out)
            except Exception:
                pass
        except Exception:
            # Deterministic fallback seed when zoned bridge is unavailable.
            try:
                init = getattr(self, "_init", None)
                if init is None:
                    self.specify_initial_conditions(-20.0, 30.0)
                    init = self._init

                p_low = float(value(fs.P_low_target))
                p_high = float(value(fs.P_high_target))
                h1 = float(init["h1"])
                h2 = float(init["h2"])
                h3 = float(init["h3"])
                h4 = float(init["h4"])
                h_g = float(init["h_g"])
                h_f = float(init["h_f"])

                fs.evaporator.cold_side_inlet.pressure[0].set_value(p_low)
                fs.evaporator.cold_side_inlet.enth_mass[0].set_value(h4)
                fs.evaporator.cold_side_outlet.enth_mass[0].set_value(h1)
                fs.compressor.inlet.pressure[0].set_value(p_low)
                fs.compressor.inlet.enth_mass[0].set_value(h1)
                fs.compressor.outlet.pressure[0].set_value(p_high)
                fs.compressor.outlet.enth_mass[0].set_value(h2)

                fs.desuperheater.hot_side_inlet.pressure[0].set_value(p_high)
                fs.desuperheater.hot_side_inlet.enth_mass[0].set_value(h2)
                fs.desuperheater.hot_side_outlet.enth_mass[0].set_value(h_g)
                fs.condenser.hot_side_inlet.pressure[0].set_value(p_high)
                fs.condenser.hot_side_inlet.enth_mass[0].set_value(h_g)
                fs.condenser.hot_side_outlet.enth_mass[0].set_value(h_f)
                fs.subcooler.hot_side_inlet.pressure[0].set_value(p_high)
                fs.subcooler.hot_side_inlet.enth_mass[0].set_value(h_f)
                fs.subcooler.hot_side_outlet.enth_mass[0].set_value(h_f - 5000.0)

                fs.expansion_valve.inlet.pressure[0].set_value(p_high)
                fs.expansion_valve.outlet.pressure[0].set_value(p_low)
                fs.expansion_valve.outlet.enth_mass[0].set_value(h4)

                t_air_sc_in = float(value(fs.subcooler.cold_side_inlet.temperature[0]))
                fs.condenser.cold_side_inlet.temperature[0].set_value(t_air_sc_in + 0.1)
                fs.desuperheater.cold_side_inlet.temperature[0].set_value(t_air_sc_in + 0.2)
            except Exception:
                pass

    def _full_profile_warmstart(self, verbose=False):
        """Solve a relaxed fully-connected cycle to generate a profile warm start.

        This routine runs one full-cycle solve with selected hard guards
        temporarily deactivated, then restores those guards. The solved state
        is kept in model variable values and used as initialization for the
        subsequent main solve.
        """
        fs = self.model.fs

        # Ensure arc equalities are active on all expanded arcs.
        arc_names = [
            "evaporator_to_compressor_expanded",
            "compressor_to_desuperheater_expanded",
            "desuperheater_to_condenser_expanded",
            "condenser_to_subcooler_expanded",
            "subcooler_to_expansion_valve_expanded",
            "expansion_valve_to_evaporator_expanded",
            "air_subcooler_to_condenser_expanded",
            "air_condenser_to_desuperheater_expanded",
        ]
        for an in arc_names:
            arc = getattr(fs, an, None)
            if arc is None:
                continue
            # Activate one flow basis and one enthalpy basis per arc to avoid
            # mixed-basis over-closure during warm-start.
            flow_mass = getattr(arc, "flow_mass_equality", None)
            flow_mol = getattr(arc, "flow_mol_equality", None)
            enth_mass = getattr(arc, "enth_mass_equality", None)
            enth_mol = getattr(arc, "enth_mol_equality", None)

            if flow_mass is not None:
                if an == "evaporator_to_compressor_expanded":
                    if flow_mass.active:
                        flow_mass.deactivate()
                elif not flow_mass.active:
                    flow_mass.activate()
                if flow_mol is not None and flow_mol.active:
                    flow_mol.deactivate()
            elif flow_mol is not None and not flow_mol.active:
                flow_mol.activate()

            if enth_mass is not None:
                if not enth_mass.active:
                    enth_mass.activate()
                if enth_mol is not None and enth_mol.active:
                    enth_mol.deactivate()
            elif enth_mol is not None and not enth_mol.active:
                enth_mol.activate()

        # Relax the most restrictive guards for warm-start solve only.
        relaxed = []
        for cname in [
            "evap_pinch_guard",
            "cond_pinch_guard",
            "subcooler_temperature_drop",
            "ds_air_heating_guard",
            "cond_air_heating_guard",
            "sc_air_heating_guard",
        ]:
            c = getattr(fs, cname, None)
            if c is not None and c.active:
                c.deactivate()
                relaxed.append(c)

        solver = get_solver()
        solver.options = {"max_iter": 500, "tol": 1e-6, "acceptable_tol": 1e-5}
        try:
            solver.solve(self.model, tee=verbose)
        except Exception:
            pass
        finally:
            for c in relaxed:
                c.activate()

    def get_full_load_cop(self):
        """Return solved full-load COP."""
        return value(self.model.fs.cop)

    def get_part_load_cop(self):
        """Return post-corrected part-load COP."""
        plf = self._compute_plf(value(self.model.fs.plr), value(self.model.fs.cd))
        return plf * value(self.model.fs.cop)

    def get_actual_sh_sc(self):
        """Return last solved SH/SC in K."""
        return getattr(self, "_last_sh_actual", float("nan")), getattr(self, "_last_sc_actual", float("nan"))

    def _print_debug_snapshot(self, reason):
        """Print compact condenser-train diagnostics for failed solves."""
        fs = self.model.fs
        try:
            print(f"[cond3-debug] reason={reason}")
            print(
                "[cond3-debug] Q_hot(W): ds={:.3f}, cond={:.3f}, sc={:.3f}".format(
                    value(fs.desuperheater.hot_side.heat[0]),
                    value(fs.condenser.hot_side.heat[0]),
                    value(fs.subcooler.hot_side.heat[0]),
                )
            )
            print(
                "[cond3-debug] h_hot_out(J/kg): comp={:.3f}, ds={:.3f}, cond={:.3f}, sc={:.3f}".format(
                    value(fs.compressor.outlet.enth_mass[0]),
                    value(fs.desuperheater.hot_side_outlet.enth_mass[0]),
                    value(fs.condenser.hot_side_outlet.enth_mass[0]),
                    value(fs.subcooler.hot_side_outlet.enth_mass[0]),
                )
            )
            print(
                "[cond3-debug] area(m2): ds={:.3f}, cond={:.3f}, sc={:.3f}, total={:.3f}".format(
                    value(fs.desuperheater.area),
                    value(fs.condenser.area),
                    value(fs.subcooler.area),
                    value(fs.desuperheater.area + fs.condenser.area + fs.subcooler.area),
                )
            )
        except Exception:
            pass
