"""
Generic two-loop cascade refrigeration cycle.

Hot loop and cold loop are independent vapor-compression cycles
(any general_helmholtz-registered fluid), coupled only via the
cascade heat exchanger: hot loop's evaporator <-> cold loop's condenser.

Default fluids: hot=r134a, cold=co2 (matches earlier validated baseline).
Swap by passing a different `fluids` dict to the constructor.

Untested -- IDAES/Pyomo cannot be run in the assistant's sandbox.
"""

import logging
from enum import Enum

from pyomo.environ import (
    ConcreteModel, Param, Var, Constraint, TransformationFactory,
    value, maximize, units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util import DiagnosticsToolbox
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock, StateVars, AmountBasis,
)
from idaes.models.unit_models import Heater, Compressor, PressureChanger

import CoolProp.CoolProp as CP

class Mode(Enum):
    IMPROVED_TPX = "improved_TPx"
    PH = "PH"

C_to_K = 273.15

DEFAULT_FLUIDS = {
    "hot": {
        "name": "r134a",          # general_helmholtz pure_component string
        "coolprop_name": "R134a", # CoolProp fluid string (for init guesses)
        "efficiency": 0.75,
        "superheat_cap": 3.0,     # K, evaporator outlet superheat cap
        "subcool_cap": 3.0,       # K, condenser outlet subcool cap
        "Tmax_C": None,           # optional safety cap on condenser outlet T
    },
    "cold": {
        "name": "co2",
        "coolprop_name": "CO2",
        "efficiency": 0.75,
        "superheat_cap": 5.0,
        "subcool_cap": 5.0,
        "Tmax_C": 25.0,           # keep well under CO2's ~31C critical T
    },
}

class CascadeCycle:
    """Generic hot-loop/cold-loop cascade cycle.

    `fluids` is a dict with keys "hot" and "cold", each a dict with:
      name           -- general_helmholtz pure_component string
      coolprop_name  -- CoolProp fluid string, used only for init guesses
      efficiency     -- compressor isentropic efficiency
      superheat_cap  -- K, evaporator outlet superheat upper bound
      subcool_cap    -- K, condenser outlet subcool upper bound
      Tmax_C         -- optional, condenser outlet T safety cap (None = off)
    """

    def __init__(self, fluids=None, mode=Mode.IMPROVED_TPX):
        self.mode = mode
        self.fluids = fluids if fluids is not None else DEFAULT_FLUIDS
        self._init_data = {}

        logging.basicConfig(level=logging.WARNING)
        self.logger = logging.getLogger(__name__)

        self.model = ConcreteModel()
        self.model.fs = FlowsheetBlock(dynamic=False)

        sv = StateVars.PH if mode == Mode.PH else StateVars.TPX

        for role in ("hot", "cold"):
            block = HelmholtzParameterBlock(
                pure_component=self.fluids[role]["name"],
                state_vars=sv, amount_basis=AmountBasis.MASS)
            setattr(self.model.fs, f"{role}_properties", block)

        for role in ("hot", "cold"):
            self._build_loop(role, getattr(self.model.fs, f"{role}_properties"))

        self._add_cascade_coupling()
        self._add_system_cop()

        self.optimization_converged = None

    def _build_loop(self, role, properties):
        fs = self.model.fs

        evaporator = Heater(property_package=properties)
        compressor = Compressor(property_package=properties)
        condenser = Heater(property_package=properties)
        expansion_valve = PressureChanger(
            property_package=properties,
            thermodynamic_assumption="adiabatic", compressor=False)

        setattr(fs, f"{role}_evaporator", evaporator)
        setattr(fs, f"{role}_compressor", compressor)
        setattr(fs, f"{role}_condenser", condenser)
        setattr(fs, f"{role}_expansion_valve", expansion_valve)

        setattr(fs, f"{role}_evap_to_comp", Arc(source=evaporator.outlet, destination=compressor.inlet))
        setattr(fs, f"{role}_comp_to_cond", Arc(source=compressor.outlet, destination=condenser.inlet))
        setattr(fs, f"{role}_cond_to_valve", Arc(source=condenser.outlet, destination=expansion_valve.inlet))
        setattr(fs, f"{role}_valve_to_evap", Arc(source=expansion_valve.outlet, destination=evaporator.inlet))

        TransformationFactory("network.expand_arcs").apply_to(self.model)

        getattr(fs, f"{role}_evap_to_comp_expanded").flow_mass_equality.deactivate()

        evaporator.superheating_max = Param(initialize=0, units=pyunits.K, mutable=True)
        evaporator.T_sat_set = Param(initialize=C_to_K, units=pyunits.K, mutable=True)

        @evaporator.Constraint(doc="Superheat cap")
        def superheating_upper_constraint(b):
            return (b.control_volume.properties_out[0].temperature
                    <= b.control_volume.properties_out[0].temperature_sat + b.superheating_max)
        evaporator.superheating_upper_constraint.deactivate()

        @evaporator.Constraint(doc="Evaporator saturation temperature setpoint")
        def evap_sat_constraint(b):
            return b.control_volume.properties_out[0].temperature_sat == b.T_sat_set
        evaporator.evap_sat_constraint.deactivate()

        condenser.subcooling_max = Param(initialize=0, units=pyunits.K, mutable=True)
        condenser.ambient_T = Param(initialize=C_to_K, units=pyunits.K, mutable=True)
        condenser.approach_T = Param(initialize=0, units=pyunits.K, mutable=True)

        @condenser.Constraint(doc="Subcool cap")
        def subcooling_lower_constraint(b):
            return (b.control_volume.properties_out[0].temperature
                    >= b.control_volume.properties_out[0].temperature_sat - b.subcooling_max)
        condenser.subcooling_lower_constraint.deactivate()

        @condenser.Constraint(doc="Condenser sat temp = ambient + approach")
        def approach_constraint(b):
            return b.control_volume.properties_out[0].temperature_sat == b.ambient_T + b.approach_T
        condenser.approach_constraint.deactivate()

        @compressor.Constraint(doc="Compressor outlet must be vapor")
        def vapor_constraint(b):
            return (b.control_volume.properties_out[0].temperature
                    >= b.control_volume.properties_out[0].temperature_sat)
        compressor.vapor_constraint.deactivate()

        compressor.ratioP.setub(8.0)
        compressor.ratioP.setlb(1.1)

        for u in (evaporator, compressor, condenser, expansion_valve):
            u.Tmin = Param(initialize=C_to_K, mutable=True)
            u.Tmax = Param(initialize=C_to_K + 10, mutable=True)

            @u.Constraint(doc="Temperature lower bound")
            def T_lower_bound(b):
                return b.control_volume.properties_out[0].temperature >= u.Tmin
            u.T_lower_bound.deactivate()

            @u.Constraint(doc="Temperature upper bound")
            def T_upper_bound(b):
                return b.control_volume.properties_out[0].temperature <= u.Tmax
            u.T_upper_bound.deactivate()

        # Generic safety cap: any fluid whose condenser sits near its
        # critical point gets an explicit Tmax on the condenser outlet.
        Tmax_C = self.fluids[role].get("Tmax_C")
        if Tmax_C is not None:
            condenser.Tmax.set_value(Tmax_C + C_to_K)
            condenser.T_upper_bound.activate()

    def _add_cascade_coupling(self):
        fs = self.model.fs
        fs.cascade_approach_dT = Param(initialize=3.0, units=pyunits.K, mutable=True)

        @fs.Constraint(doc="Cascade HX energy balance (adiabatic)")
        def cascade_energy_balance(b):
            return b.hot_evaporator.heat_duty[0] == -b.cold_condenser.heat_duty[0]

        @fs.Constraint(doc="Cascade HX temperature approach")
        def cascade_approach_constraint(b):
            return (b.hot_evaporator.control_volume.properties_out[0].temperature_sat
                    == b.cold_condenser.control_volume.properties_out[0].temperature_sat
                    - b.cascade_approach_dT)

        fs.cascade_energy_balance.deactivate()
        fs.cascade_approach_constraint.deactivate()

    def _add_system_cop(self):
        fs = self.model.fs
        fs.cop_cascade = Var(initialize=1, units=pyunits.dimensionless, bounds=(0.01, 100))

        @fs.Constraint(doc="Cascade COP definition")
        def compute_cop_cascade(b):
            return (b.cop_cascade * (b.hot_compressor.work_mechanical[0]
                                      + b.cold_compressor.work_mechanical[0])
                    == b.cold_evaporator.heat_duty[0])

        @fs.Objective(doc="Maximize cascade COP", sense=maximize)
        def obj(b):
            return b.cop_cascade

        fs.compute_cop_cascade.deactivate()
        fs.obj.deactivate()

    def specify_initial_conditions(self, hot_ambient_C=35.0, cold_evap_C=-27.0):
        hot_evap_guess_C = 5.0
        self._specify_loop_initial_conditions("hot", hot_evap_guess_C, hot_ambient_C)
        cold_cond_guess_C = hot_evap_guess_C + 3.0
        self._specify_loop_initial_conditions("cold", cold_evap_C, cold_cond_guess_C)

    def _specify_loop_initial_conditions(self, role, low_side_C, high_side_C):
        coolprop_name = self.fluids[role]["coolprop_name"]
        low_K, high_K = low_side_C + C_to_K, high_side_C + C_to_K
        superheat = subcool = 3

        low_P = CP.PropsSI('P', 'T', low_K, 'Q', 0, coolprop_name)
        high_P = CP.PropsSI('P', 'T', high_K, 'Q', 0, coolprop_name)
        low_liquid_H = CP.PropsSI('H', 'T', low_K, 'Q', 0.2, coolprop_name)
        low_vapor_H = CP.PropsSI('H', 'T', low_K + superheat, 'Q', 1, coolprop_name)
        high_liquid_H = CP.PropsSI('H', 'T', high_K + superheat, 'Q', 0, coolprop_name)
        high_vapor_H = CP.PropsSI('H', 'T', high_K - subcool, 'Q', 1, coolprop_name)

        self._init_data[role] = {
            "h": [low_vapor_H, high_vapor_H, high_liquid_H, low_liquid_H],
            "p": [low_P, high_P, high_P, low_P],
            "T": [low_K + superheat, high_K + superheat, high_K - subcool, low_K],
        }

    def initialize(self, verbose=False):
        self._initialize_loop("hot", verbose)
        self._initialize_loop("cold", verbose)

    def _initialize_loop(self, role, verbose=False):
        fs = self.model.fs
        evaporator = getattr(fs, f"{role}_evaporator")
        compressor = getattr(fs, f"{role}_compressor")
        condenser = getattr(fs, f"{role}_condenser")
        expansion_valve = getattr(fs, f"{role}_expansion_valve")
        h, p, T = (self._init_data[role][k] for k in ("h", "p", "T"))

        evaporator.inlet.flow_mass[0].fix(1)
        evaporator.inlet.pressure[0].fix(p[-1])
        if self.mode == Mode.PH:
            evaporator.inlet.enth_mass[0].fix(h[-1])
            evaporator.outlet.enth_mass[0].fix(h[0])
        else:
            evaporator.inlet.temperature[0].fix(T[-1])
            evaporator.outlet.temperature[0].fix(T[0])
        evaporator.initialize(outlvl=logging.WARNING)
        propagate_state(getattr(fs, f"{role}_evap_to_comp"))

        compressor.inlet.pressure[0].fix(p[0])
        if self.mode == Mode.PH:
            compressor.inlet.enth_mass[0].fix(h[0])
        else:
            compressor.inlet.temperature[0].fix(T[0])
        compressor.outlet.pressure[0].fix(p[1])
        compressor.efficiency_isentropic[0].fix(self.fluids[role]["efficiency"])
        compressor.initialize(outlvl=logging.WARNING)
        propagate_state(getattr(fs, f"{role}_comp_to_cond"))

        condenser.inlet.pressure[0].fix(p[1])
        if self.mode == Mode.PH:
            condenser.inlet.enth_mass[0].fix(h[1])
            condenser.outlet.enth_mass[0].fix(h[2])
        else:
            condenser.inlet.temperature[0].fix(T[1])
            condenser.outlet.temperature[0].fix(T[2])
        condenser.initialize(outlvl=logging.WARNING)
        propagate_state(getattr(fs, f"{role}_cond_to_valve"))

        expansion_valve.inlet.pressure[0].fix(p[2])
        expansion_valve.outlet.pressure[0].fix(p[3])
        expansion_valve.initialize(outlvl=logging.WARNING)
        propagate_state(getattr(fs, f"{role}_valve_to_evap"))

        if verbose:
            print(f"{role} ({self.fluids[role]['name']}) loop initialized.")

    def set_specifications(self, hot_ambient_C=35.0, hot_condenser_approach_C=9.0,
                            cold_evap_C=-27.0, cascade_approach_dT=3.0):
        fs = self.model.fs

        for role in ("hot", "cold"):
            self._unfix_loop(role)

        fs.cold_evaporator.inlet.flow_mass[0].fix(1.0)

        fs.hot_condenser.ambient_T.set_value(hot_ambient_C + C_to_K)
        fs.hot_condenser.approach_T.set_value(hot_condenser_approach_C)
        fs.hot_condenser.approach_constraint.activate()

        fs.cold_evaporator.T_sat_set.set_value(cold_evap_C + C_to_K)
        fs.cold_evaporator.evap_sat_constraint.activate()

        for role in ("hot", "cold"):
            self._activate_superheat_subcool(role)

        fs.cascade_approach_dT.set_value(cascade_approach_dT)
        fs.cascade_energy_balance.activate()
        fs.cascade_approach_constraint.activate()

        calculate_scaling_factors(self.model)

    def _unfix_loop(self, role):
        fs = self.model.fs
        for name in ("evaporator", "compressor", "condenser", "expansion_valve"):
            unit = getattr(fs, f"{role}_{name}")
            unit.inlet.flow_mass[0].unfix()
            unit.outlet.flow_mass[0].unfix()
            if self.mode == Mode.PH:
                unit.inlet.enth_mass[0].unfix()
                unit.outlet.enth_mass[0].unfix()
            else:
                unit.inlet.temperature[0].unfix()
                unit.inlet.vapor_frac[0].unfix()
                unit.outlet.temperature[0].unfix()
                unit.outlet.vapor_frac[0].unfix()
                unit.inlet.vapor_frac[0].setlb(0)
                unit.inlet.vapor_frac[0].setub(1)
            unit.inlet.pressure[0].unfix()
            unit.outlet.pressure[0].unfix()
            unit.T_lower_bound.deactivate()
            # Keep the condenser T_upper_bound active only if this fluid
            # role has a Tmax_C configured (generic replacement for the
            # old "CO2 condenser only" special case).
            has_Tmax = self.fluids[role].get("Tmax_C") is not None
            if not (name == "condenser" and has_Tmax):
                unit.T_upper_bound.deactivate()

        evaporator = getattr(fs, f"{role}_evaporator")
        compressor = getattr(fs, f"{role}_compressor")
        condenser = getattr(fs, f"{role}_condenser")
        expansion_valve = getattr(fs, f"{role}_expansion_valve")

        if self.mode == Mode.IMPROVED_TPX:
            evaporator.outlet.vapor_frac[0].fix(1.0)
            evaporator.control_volume.properties_out[0.0].eq_complementarity.deactivate()

            compressor.outlet.vapor_frac[0].fix(1.0)
            compressor.control_volume.properties_out[0.0].eq_complementarity.deactivate()
            compressor.vapor_constraint.activate()

            condenser.outlet.vapor_frac[0].fix(0.0)
            condenser.control_volume.properties_out[0.0].eq_complementarity.deactivate()

            expansion_valve.control_volume.properties_out[0.0].eq_complementarity.deactivate()
            expansion_valve.control_volume.properties_out[0.0].eq_sat.activate()

    def _activate_superheat_subcool(self, role):
        fs = self.model.fs
        fluid_cfg = self.fluids[role]
        getattr(fs, f"{role}_evaporator").superheating_max.set_value(fluid_cfg["superheat_cap"])
        getattr(fs, f"{role}_evaporator").superheating_upper_constraint.activate()
        getattr(fs, f"{role}_condenser").subcooling_max.set_value(fluid_cfg["subcool_cap"])
        getattr(fs, f"{role}_condenser").subcooling_lower_constraint.activate()

    def optimize_COP(self, verbose=False, initialize=True, optimize=True):
        solver = get_solver()
        solver.options = {"max_iter": 1000, "tol": 1e-6, "linear_solver": "ma57"}
        fs = self.model.fs

        if initialize:
            fs.compute_cop_cascade.deactivate()
            fs.obj.deactivate()
            results = solver.solve(self.model, tee=verbose)
            if results.solver.termination_condition != "optimal":
                self.logger.error("Feasibility initialization failed")

        if optimize:
            fs.compute_cop_cascade.activate()
            fs.obj.activate()
        else:
            fs.compute_cop_cascade.deactivate()
            fs.obj.deactivate()

        results = solver.solve(self.model, tee=verbose)
        for _ in range(2):
            if results.solver.termination_condition != "optimal":
                results = solver.solve(self.model, tee=verbose)

        converged = results.solver.termination_condition == "optimal"
        self.optimization_converged = converged
        if not converged:
            DiagnosticsToolbox(self.model, constraint_residual_tolerance=1e-6).display_constraints_with_large_residuals()

        return (value(fs.cop_cascade), converged) if converged else (None, converged)


if __name__ == "__main__":
    # Default fluids: hot=R134a, cold=CO2
    cycle = CascadeCycle()
    cycle.specify_initial_conditions()
    cycle.initialize(verbose=True)
    cycle.set_specifications()
    cop, converged = cycle.optimize_COP(verbose=True)
    print("Cascade COP:", cop, "| Converged:", converged)
