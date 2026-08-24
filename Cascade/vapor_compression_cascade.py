"""

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 2026-08-20

vapor_compression_cascade.py

Two-loop cascade vapor-compression cycle, built from this project's own
single-stage `vapor_compression_plr_r1234yf.py` baseline. Same PLR/CD
part-load metadata, matplotlib diagram helpers, and debug
arc-pressure-equality override as that file -- duplicated per loop
(hot, cold) and coupled through a cascade heat exchanger:

    hot loop's evaporator  <-- (adiabatic energy balance + approach dT) -->  cold loop's condenser

Fluids are NOT hardcoded -- pass `fluids={"hot": {...}, "cold": {...}}`
to the constructor. Each entry:
    gh_component   -- general_helmholtz `pure_component` string
    coolprop_name  -- CoolProp fluid string (may differ in
                       capitalization from gh_component, e.g. R1234yf)
    custom_json    -- path to a general_helmholtz parameter JSON to
                       register before use, or None for a fluid IDAES
                       already ships (e.g. "r134a", "co2"). Mirrors this
                       project's original R1234yf registration mechanism.
    efficiency     -- compressor isentropic efficiency
    superheat_max  -- K, evaporator outlet superheat cap
    subcool_max    -- K, condenser outlet subcool cap
    Tmax_C         -- optional condenser-outlet safety cap (None = off;
                       used for CO2 to stay clear of its ~31C Tc)

Default fluids: hot=R134a, cold=CO2 (both IDAES built-ins, no
registration needed) -- the cascade baseline agreed in BREADCRUMB.md.

Mass flow: exactly ONE loop's mass flow is fixed (`flow_fixed_role`,
default "cold" @ 1 kg/s); the other loop's mass flow is left free so the
solver can satisfy the cascade energy balance. Fixing both loops'
flows independently over-constrains that balance -- see BREADCRUMB.md
Section 4, open question 4.
"""

import os
import shutil
import logging
from enum import Enum

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.initialization import propagate_state
from idaes.core.util import DiagnosticsToolbox
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
    get_parameter_path,
    set_parameter_path,
)
from idaes.models.properties.general_helmholtz.helmholtz_parameters import WriteParameters
from idaes.models.unit_models import (
    Heater, Compressor, PressureChanger,
)
from pyomo.environ import (
    ConcreteModel, value, Objective, maximize, 
    TransformationFactory, Param, Var, Constraint, units as pyunits,
)
from pyomo.network import Arc


class Mode(Enum):
    ORIGINAL_TPX = "original_TPx"
    IMPROVED_TPX = "improved_TPx"
    PH = "PH"


C_to_K = 273.15
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_FLUIDS = {
    "hot": {
        "gh_component": "r134a",
        "coolprop_name": "R134a",
        "custom_json": None,
        "efficiency": 0.75,
        "superheat_max": 3.0,
        "subcool_max": 3.0,
        "Tmax_C": None,
    },
    "cold": {
        "gh_component": "co2",
        "coolprop_name": "CO2",
        "custom_json": None,
        "efficiency": 0.75,
        "superheat_max": 5.0,
        "subcool_max": 5.0,
        # Tmax_C=None: the 25C cap tried here earlier was an undemonstrated
        # guess ("insurance", ~20-25C under CO2's 31C Tc) and empirically
        # was NOT the operative constraint anyway -- removing it entirely
        # produced the identical restoration failure (same ~18.7 constraint
        # violation) as keeping it. The real fix belongs on the HOT loop's
        # own evaporator (currently unbounded -- see BREADCRUMB.md), not
        # here. Leave this None until/unless a real derived bound is found.
        "Tmax_C": None,
    },
}


def _register_custom_fluid(cfg):
    """Registers a general_helmholtz parameter JSON for a fluid that
    isn't one of IDAES's built-ins (mirrors this project's original
    R1234yf registration mechanism). No-op if `custom_json` is None --
    built-in fluids (r134a, co2, ...) need nothing beyond
    `pure_component=...`."""
    json_path = cfg.get("custom_json")
    if not json_path:
        return

    # WriteParameters.write() writes relative to CWD, not the JSON's
    # directory -- force it to land in _THIS_DIR.
    _prev_cwd = os.getcwd()
    try:
        os.chdir(_THIS_DIR)
        WriteParameters(parameters=json_path).write(dry_run=False)
    finally:
        os.chdir(_prev_cwd)

    # set_parameter_path() is destructive (clear_component_registry() +
    # rescan-only-this-dir) -- merge the current registry with _THIS_DIR
    # so we don't deregister anything already registered (e.g. the other
    # loop's fluid, if it's also custom, or IDAES's own built-ins).
    current_path = get_parameter_path()
    merged_dir = os.path.join(_THIS_DIR, "_merged_helmholtz_params")
    os.makedirs(merged_dir, exist_ok=True)
    for src_dir in (current_path, _THIS_DIR):
        for fname in os.listdir(src_dir):
            src = os.path.join(src_dir, fname)
            if not os.path.isfile(src):
                continue
            dst = os.path.join(merged_dir, fname)
            if not os.path.exists(dst) or os.path.getmtime(src) > os.path.getmtime(dst):
                shutil.copy2(src, dst)
    set_parameter_path(merged_dir)


# ---------------------------------------------------------------------------
# Fluid envelope helpers -- option (c): bounds are DERIVED from each fluid's
# own saturation curve at set_specifications time, never hardcoded.
#
# History (see BREADCRUMB.md): this file originally inherited the single-stage
# R134a model's constants -- low side 50-2000 kPa, high side 100-5000 kPa,
# condenser outlet 30-50 C. Those are correct for R134a (88.5 kPa at -29C,
# 1455 kPa at 54C) and catastrophically wrong for CO2 (2649 kPa at -10C,
# 7214 kPa at 30C). The result was a model that was infeasible by construction
# -- constant residuals of 18.7158 K (subcooling) and 0.6487 MPa (pressure)
# that never moved with operating conditions. Never hardcode a pressure or
# temperature bound in this file again; derive it from the fluid.
#
# Why staying clear of a fluid's critical point is a MODELING VALIDITY
# requirement, not just a numerical-safety margin: this file implements a
# SUBCRITICAL vapor-compression cycle -- the condenser's job is to reject
# heat by actually condensing vapor into liquid, a real phase change, which
# only exists below the critical point. Above Tc there is no liquid/vapor
# distinction left to condense between; a real system operating that warm
# runs a TRANSCRITICAL cycle instead (a gas cooler in place of a condenser,
# different unit operations and equations entirely), which this file does
# not implement. So a case that pushes a loop's condenser near or above its
# fluid's Tc isn't just numerically fragile -- it's asking this cycle
# architecture a question it cannot physically answer. This matters most
# for CO2 (Tc ~31C, uncomfortably close to plausible cascade operating
# temperatures) but the principle applies to any fluid/loop in this file.
# ---------------------------------------------------------------------------

def _fluid_envelope(coolprop_name):
    """Physical two-phase envelope of a fluid, from CoolProp."""
    return {
        "Tcrit_K": CP.PropsSI("Tcrit", coolprop_name),
        "Pcrit_Pa": CP.PropsSI("Pcrit", coolprop_name),
        "Ttriple_K": CP.PropsSI("Ttriple", coolprop_name),
        "Ptriple_Pa": CP.PropsSI("ptriple", coolprop_name),
    }


def _psat_Pa(coolprop_name, T_C):
    """Saturation pressure (Pa) at a given temperature (degC)."""
    return CP.PropsSI("P", "T", T_C + C_to_K, "Q", 0, coolprop_name)


def _clamp_sat_T_C(coolprop_name, T_C, env=None, T_margin_K=1.0):
    """Pull a requested saturation temperature back inside this fluid's
    triple-point-to-critical-point range, leaving a buffer of T_margin_K
    at each end. Right at Tc, liquid and vapor properties genuinely
    converge to the same values -- the two phases stop being physically
    distinguishable. That convergence is exactly what makes it numerically
    hard for general_helmholtz to keep tracking them as two separate
    roots (liquid density, vapor density) instead of one merged root.
    The buffer sidesteps that region rather than fixing it."""
    env = env or _fluid_envelope(coolprop_name)
    T_min_C = env["Ttriple_K"] - C_to_K + T_margin_K
    T_max_C = env["Tcrit_K"] - C_to_K - T_margin_K
    return max(T_min_C, min(T_C, T_max_C))


def derive_pressure_bounds(coolprop_name, T_sat_lo_C, T_sat_hi_C,
                            margin_frac=0.30, P_crit_frac=0.95):
    """Derive (P_min, P_max) in Pa spanning a saturation-temperature range.

    Widens by `margin_frac` so the solver has room to move, then clamps to
    the fluid's physical envelope: never below the triple-point pressure,
    never above `P_crit_frac` of the critical pressure. That upper clamp is
    a modeling-validity boundary, not just solver insurance -- see the
    module comment above `_fluid_envelope` for why this file's cycle
    architecture is only meaningful subcritical.

    Returns (P_min_Pa, P_max_Pa).
    """
    env = _fluid_envelope(coolprop_name)
    T_lo = _clamp_sat_T_C(coolprop_name, min(T_sat_lo_C, T_sat_hi_C), env)
    T_hi = _clamp_sat_T_C(coolprop_name, max(T_sat_lo_C, T_sat_hi_C), env)

    P_min = _psat_Pa(coolprop_name, T_lo) * (1.0 - margin_frac)
    P_max = _psat_Pa(coolprop_name, T_hi) * (1.0 + margin_frac)

    P_min = max(P_min, env["Ptriple_Pa"])
    P_max = min(P_max, P_crit_frac * env["Pcrit_Pa"])
    if P_max <= P_min:  # degenerate (fluid can't span the request) -- widen
        P_max = min(P_min * (1.0 + margin_frac), P_crit_frac * env["Pcrit_Pa"])
    return P_min, P_max


# Defaults for the per-loop kwargs accepted by set_specifications(hot=..., cold=...).
# Mirrors the single-stage file's set_specifications signature, split per role.
#
# EVERY bound below defaults to None, meaning "derive it" (option c). Passing an
# explicit value overrides the derivation for that one bound. Do not restore
# hardcoded defaults here -- see the header comment above.
_LOOP_SPEC_DEFAULTS = {
    "low_side_pressure": None,             # (lo, hi) kPa -- None = derive
    "high_side_pressure": None,            # (lo, hi) kPa -- None = derive
    "evaporator_temperature": None,        # (lo, hi) degC -- None = no bound
    "compressor_temperature": None,        # (lo, hi) degC -- None = no bound
    "condenser_temperature": None,         # (lo, hi) degC -- None = no bound
    "expansion_valve_temperature": None,   # (lo, hi) degC -- None = no bound
    "subcooling": 3,                       # degC
    "superheating": 3,                     # degC
    "subcooling_max": None,                # defaults to fluids[role]["subcool_max"]
    "superheating_max": None,              # defaults to fluids[role]["superheat_max"]
    "max_pressure_ratio": 8,
    "ambient_temperature": None,           # degC -- ambient-facing condenser only
    "condenser_approach": None,            # degC -- ambient-facing condenser only
    "evap_sat_temperature": None,          # degC -- load-facing evaporator only
}


class CascadeCycle:
    """Two-loop cascade vapor-compression cycle, generic hot/cold fluids.

    hot loop:  rejects heat to ambient (condenser), draws heat from the
               cascade HX (its evaporator).
    cold loop: draws heat from the cooled space (its evaporator), rejects
               heat into the cascade HX (its condenser).

    Cascade HX = 2 constraints reusing the existing Heater blocks
    (adiabatic energy balance + fixed approach dT between the two
    saturation temperatures), NOT a dedicated IDAES HeatExchanger unit
    model -- see BREADCRUMB.md Section 2c (real HX model deferred).
    """

    ROLES = ("hot", "cold")

    # Option (b): what each heat exchanger actually faces is an EXPLICIT
    # structural property of the cascade, not something inferred from which
    # kwargs a caller happens to pass.
    #
    # This matters: the original code decided whether to apply an
    # ambient-derived condenser temperature bound via
    #     use_condenser_temperature_bounds = not (ambient_temperature is not None
    #                                             and condenser_approach is not None)
    # The hot loop passes both kwargs -> bound skipped. The cold loop passes
    # only evap_sat_temperature -> bound APPLIED, silently pinning CO2's
    # cascade-facing condenser to an ambient-shaped 30-50 C. Making "faces"
    # explicit removes the possibility of that class of bug entirely.
    FACING = {
        "hot":  {"evaporator": "cascade", "condenser": "ambient"},
        "cold": {"evaporator": "load",    "condenser": "cascade"},
    }

    def __init__(self, fluids=None, PLR=0.5, CD=0.25, mode=Mode.IMPROVED_TPX,
                 flow_fixed_role="cold"):
        self.mode = mode
        self.fluids = fluids if fluids is not None else DEFAULT_FLUIDS
        assert flow_fixed_role in self.ROLES # hot/cold
        self.flow_fixed_role = flow_fixed_role

        for role in self.ROLES:
            eff = self.fluids[role]["efficiency"]
            assert 0 < eff < 1, f"{role} compressor efficiency must be in (0,1)"
        assert 0.0 <= PLR <= 1.0, "PLR must be in [0,1]"
        assert 0.0 <= CD <= 1.0, "CD must be in [0,1]"
        self.plr = PLR
        self.cd = CD

        logging.basicConfig(level=logging.WARNING)
        self.logger = logging.getLogger(__name__)

        for role in self.ROLES:
            _register_custom_fluid(self.fluids[role])

        self.model = ConcreteModel()
        self.model.fs = FlowsheetBlock(dynamic=False)

        sv = StateVars.PH if mode == Mode.PH else StateVars.TPX
        for role in self.ROLES:
            block = HelmholtzParameterBlock(
                pure_component=self.fluids[role]["gh_component"],
                state_vars=sv, amount_basis=AmountBasis.MASS)
            setattr(self.model.fs, f"{role}_properties", block)

        self.optimization_converged = None
        self._h_init = {}
        self._p_init = {}
        self._T_init = {}
        self._unit_operations = {}

        self._define_flowsheet()

    # @staticmethod: this is pure math on the two numbers passed in -- it
    # never reads or writes anything on `self` (no cascade-specific state,
    # no fluids, no solved model), so it doesn't need to know which
    # CascadeCycle instance (if any) is calling it. Removing @staticmethod
    # would NOT be a no-op: Python would still auto-pass the instance as an
    # extra hidden first argument on any `self._compute_plf(...)` call,
    # which this 2-argument signature has no slot for -- an immediate
    # "too many arguments" error, not a behavior change. To drop
    # @staticmethod safely you'd also need to add an unused `self` param.
    @staticmethod
    def _compute_plf(plr: float, cd: float) -> float:
        """Compute part-load factor using PLF = 1 - CD * (1 - PLR)."""
        if not (0.0 <= plr <= 1.0):
            raise ValueError("PLR must be in [0, 1]")
        if not (0.0 <= cd <= 1.0):
            raise ValueError("CD must be in [0, 1]")
        return max(0.0, min(1.0, 1.0 - cd * (1.0 - plr)))

    def _fluid_name(self, role):
        return self.fluids[role]["coolprop_name"]

    # ------------------------------------------------------------------
    # Flowsheet construction
    # ------------------------------------------------------------------

    def _define_flowsheet(self):
        fs = self.model.fs

        for role in self.ROLES:
            self._build_loop(role)

        self._add_cascade_coupling()

        # PLR/CD metadata only; does not alter thermodynamic equations
        # (system-level, applied to the cascade COP as a whole -- see
        # BREADCRUMB.md Section 4, open question 3).
        fs.plr = Param(initialize=self.plr, units=pyunits.dimensionless, mutable=True)
        fs.cd = Param(initialize=self.cd, units=pyunits.dimensionless, mutable=True)

        # `cop` is NOT computed after the fact -- it's a real decision Var,
        # solved SIMULTANEOUSLY with every pressure/temperature/flow in the
        # model (equation-oriented solve, not sequential calculation). It
        # has to be a Var, not a post-hoc `duty/work` calculation, because
        # it's the thing `obj` below maximizes -- the solver can only push
        # a real variable's value up during its search, not a number that
        # only exists after the search is already over. `compute_cop`
        # isn't "the formula that produces cop", it's a leash: it forces
        # whatever value `cop` takes at every step of the search to stay
        # consistent with the real physics (duty/work). At the final
        # solution `cop` does equal duty/work, same answer as computing it
        # by hand afterward -- it just got there by being solved alongside
        # everything else, not calculated as a separate last step.
        #
        # The constraint is also written as `cop * work == duty` rather
        # than `cop == duty / work` to avoid an actual division inside the
        # model -- if the solver ever visits a point where `work` is zero
        # or near-zero while searching, `duty/work` blows up numerically;
        # the multiplied-out form never has that problem for the same math.
        fs.cop = Var(initialize=1, units=pyunits.dimensionless, bounds=(0.01, 100))

        @fs.Constraint(doc="Cascade COP: cold-loop evaporator load / total compressor work")
        def compute_cop(b):
            return (b.cop * (b.hot_compressor.work_mechanical[0]
                              + b.cold_compressor.work_mechanical[0])
                    == b.cold_evaporator.heat_duty[0])

        @fs.Objective(doc="Maximize cascade COP", sense=maximize)
        def obj(b):
            return b.cop

        fs.compute_cop.deactivate()
        fs.obj.deactivate()

    def _build_loop(self, role):
        fs = self.model.fs
        properties = getattr(fs, f"{role}_properties")

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

        setattr(fs, f"{role}_evaporator_to_compressor",
                Arc(source=evaporator.outlet, destination=compressor.inlet))
        setattr(fs, f"{role}_compressor_to_condenser",
                Arc(source=compressor.outlet, destination=condenser.inlet))
        setattr(fs, f"{role}_condenser_to_expansion_valve",
                Arc(source=condenser.outlet, destination=expansion_valve.inlet))
        setattr(fs, f"{role}_expansion_valve_to_evaporator",
                Arc(source=expansion_valve.outlet, destination=evaporator.inlet))

        TransformationFactory('network.expand_arcs').apply_to(self.model)

        # Closed loop -- deactivate the redundant flow equality on one arc.
        #
        # `expand_arcs` turned each Arc above into a real bundle of equations
        # (flow/pressure/enthalpy match) named "<arc>_expanded" -- `getattr`
        # fetches this role's bundle by the same computed-name trick `setattr`
        # used to build it. This loop is closed (evaporator -> compressor ->
        # condenser -> valve -> back to evaporator), with 4 arcs, each
        # normally enforcing its own flow match. But if evap=comp, comp=cond,
        # and cond=valve are all already enforced, valve=evap is automatically
        # implied -- stating it too is one equation more than there are
        # genuinely unknown flows, an over-constrained/redundant equation, not
        # new information. So exactly one of the four (this one, evaporator->
        # compressor) is switched off; the other three already guarantee it.
        getattr(fs, f"{role}_evaporator_to_compressor_expanded").flow_mass_equality.deactivate()

        # --- Evaporator: superheat floor + cap, saturation-temperature setpoint ---
        evaporator.superheating = Param(initialize=0, units=pyunits.K, mutable=True)
        evaporator.superheating_max = Param(initialize=8, units=pyunits.K, mutable=True)
        evaporator.T_sat_set = Param(initialize=C_to_K, units=pyunits.K, mutable=True)

        evaporator.superheating_constraint = Constraint(
            expr=evaporator.control_volume.properties_out[0].temperature
            >= evaporator.control_volume.properties_out[0].temperature_sat + evaporator.superheating)
        evaporator.superheating_constraint.deactivate()

        evaporator.superheating_upper_constraint = Constraint(
            expr=evaporator.control_volume.properties_out[0].temperature
            <= evaporator.control_volume.properties_out[0].temperature_sat + evaporator.superheating_max)
        evaporator.superheating_upper_constraint.deactivate()

        evaporator.evap_sat_constraint = Constraint(
            expr=evaporator.control_volume.properties_out[0].temperature_sat == evaporator.T_sat_set)
        evaporator.evap_sat_constraint.deactivate()

        # --- Condenser: subcool floor + cap, ambient-approach setpoint ---
        condenser.subcooling = Param(initialize=0, units=pyunits.K, mutable=True)
        condenser.subcooling_max = Param(initialize=8, units=pyunits.K, mutable=True)
        condenser.ambient_T = Param(initialize=C_to_K, units=pyunits.K, mutable=True)
        condenser.approach_T = Param(initialize=0, units=pyunits.K, mutable=True)

        condenser.subcooling_constraint = Constraint(
            expr=condenser.control_volume.properties_out[0].temperature
            <= condenser.control_volume.properties_out[0].temperature_sat - condenser.subcooling)
        condenser.subcooling_constraint.deactivate()

        condenser.subcooling_lower_constraint = Constraint(
            expr=condenser.control_volume.properties_out[0].temperature
            >= condenser.control_volume.properties_out[0].temperature_sat - condenser.subcooling_max)
        condenser.subcooling_lower_constraint.deactivate()

        condenser.approach_constraint = Constraint(
            expr=condenser.control_volume.properties_out[0].temperature_sat
            == condenser.ambient_T + condenser.approach_T)
        condenser.approach_constraint.deactivate()

        # --- Compressor: outlet must be vapor ---
        compressor.vapor_constraint = Constraint(
            expr=compressor.control_volume.properties_out[0].temperature
            >= compressor.control_volume.properties_out[0].temperature_sat)
        compressor.vapor_constraint.deactivate()

        # --- Debug: explicit pressure levels (used when arc pressure equalities are disabled) ---
        P_low = Var(initialize=200e3, units=pyunits.Pa, bounds=(50e3, 2000e3))
        P_high = Var(initialize=1200e3, units=pyunits.Pa, bounds=(100e3, 5000e3))
        P_low_target = Param(initialize=200e3, units=pyunits.Pa, mutable=True)
        P_high_target = Param(initialize=1200e3, units=pyunits.Pa, mutable=True)
        setattr(fs, f"{role}_P_low", P_low)
        setattr(fs, f"{role}_P_high", P_high)
        setattr(fs, f"{role}_P_low_target", P_low_target)
        setattr(fs, f"{role}_P_high_target", P_high_target)

        debug_constraints = {
            "P_low_target_constraint": Constraint(expr=P_low == P_low_target),
            "P_high_target_constraint": Constraint(expr=P_high == P_high_target),
            "P_low_evap_in": Constraint(expr=evaporator.inlet.pressure[0] == P_low),
            "P_low_evap_out": Constraint(expr=evaporator.outlet.pressure[0] == P_low),
            "P_low_comp_in": Constraint(expr=compressor.inlet.pressure[0] == P_low),
            "P_low_valve_out": Constraint(expr=expansion_valve.outlet.pressure[0] == P_low),
            "P_high_comp_out": Constraint(expr=compressor.outlet.pressure[0] == P_high),
            "P_high_cond_in": Constraint(expr=condenser.inlet.pressure[0] == P_high),
            "P_high_cond_out": Constraint(expr=condenser.outlet.pressure[0] == P_high),
            "P_high_valve_in": Constraint(expr=expansion_valve.inlet.pressure[0] == P_high),
        }
        for name, con in debug_constraints.items():
            setattr(fs, f"{role}_{name}", con)
            con.deactivate()

        # --- Temperature bounds on all four units (needed for PH state vars) ---
        for u in (evaporator, compressor, condenser, expansion_valve):
            u.Tmin = Param(initialize=C_to_K, mutable=True)
            u.Tmax = Param(initialize=C_to_K + 10, mutable=True)
            u.T_lower_bound = Constraint(expr=u.control_volume.properties_out[0].temperature >= u.Tmin)
            u.T_lower_bound.deactivate()
            u.T_upper_bound = Constraint(expr=u.control_volume.properties_out[0].temperature <= u.Tmax)
            u.T_upper_bound.deactivate()

        # Generic safety cap: if this fluid's condenser sits near its
        # critical point (e.g. CO2), activate a hard Tmax on the outlet.
        Tmax_C = self.fluids[role].get("Tmax_C")
        if Tmax_C is not None:
            condenser.Tmax.set_value(Tmax_C + C_to_K)
            condenser.T_upper_bound.activate()

        self._unit_operations[role] = [evaporator, compressor, condenser, expansion_valve]

    def _add_cascade_coupling(self):
        fs = self.model.fs
        fs.cascade_approach_dT = Param(initialize=3.0, units=pyunits.K, mutable=True)

        fs.cascade_energy_balance = Constraint(
            doc="Cascade HX energy balance (adiabatic)",
            expr=fs.hot_evaporator.heat_duty[0] == -fs.cold_condenser.heat_duty[0])

        fs.cascade_approach_constraint = Constraint(
            doc="Cascade HX temperature approach",
            expr=fs.hot_evaporator.control_volume.properties_out[0].temperature_sat
            == fs.cold_condenser.control_volume.properties_out[0].temperature_sat - fs.cascade_approach_dT)

        fs.cascade_energy_balance.deactivate()
        fs.cascade_approach_constraint.deactivate()

    # ------------------------------------------------------------------
    # Diagrams
    # ------------------------------------------------------------------

    def draw_thermodynamic_diagrams(self, role):
        properties = getattr(self.model.fs, f"{role}_properties")
        properties.hp_diagram()
        plt.show()
        properties.pt_diagram()
        plt.show()
        properties.ts_diagram()
        plt.show()

    # ------------------------------------------------------------------
    # Initial conditions / initialization
    # ------------------------------------------------------------------

    # High-level entry point: pick rough starting-guess temperatures for
    # BOTH loops, before any real solving happens. Hot loop's evaporator
    # (cascade-facing side) gets an arbitrary 5C guess; its condenser gets
    # whatever ambient is passed in, since the hot loop rejects straight to
    # ambient. Cold loop's evaporator gets the real target delivery temp
    # (not arbitrary -- this is the actual setpoint you want); its condenser
    # gets a rough 3K-above-hot-evap guess, echoing the real approach-dT
    # idea just as a seed -- the real cascade coupling constraint takes
    # over once set_specifications runs for real.
    def specify_initial_conditions(self, hot_ambient_C=35.0, cold_evap_C=-27.0, plot=False):
        """Compute CoolProp-based initial guesses for both loops.

        hot_evap_guess_C is an arbitrary intermediate temperature for the
        hot loop's evaporator (= the cascade HX on the hot side); the
        cold loop's condenser guess is offset from it by a nominal 3K
        approach (overwritten for real once set_specifications runs).
        """
        hot_evap_guess_C = 5.0
        self._specify_loop_initial_conditions("hot", hot_evap_guess_C, hot_ambient_C, plot=plot)

        cold_cond_guess_C = hot_evap_guess_C + 3.0
        self._specify_loop_initial_conditions("cold", cold_evap_C, cold_cond_guess_C, plot=plot)

    # Pure CoolProp math, no Pyomo yet: given just a "low side" and "high
    # side" guessed temperature, compute rough starting values for all four
    # state points around ONE loop (evaporator/compressor/condenser/valve
    # outlets) -- a fixed 3K superheat/subcool assumption is used here only
    # as a seed, separate from whatever real superheat/subcool gets set
    # later in set_specifications. Results are stashed in _h_init/_p_init/
    # _T_init, consumed next by initialize()/_initialize_loop().
    #
    # NOTE -- possible label mismatch, not confirmed to matter: the inline
    # comments on the two middle variables look swapped relative to where
    # they land in h_init below. `high_side_liquid_H` is commented
    # "compressor outlet" but ends up in the CONDENSER's array slot;
    # `high_side_vapor_H` is commented "condenser outlet" but ends up in
    # the COMPRESSOR's slot. Since these are only rough starting guesses
    # for the solver (never enforced values), this likely doesn't break
    # anything -- but worth checking against whatever file this was copied
    # from, in case unit order changed during the copy and these
    # names/comments didn't get updated to match.
    def _specify_loop_initial_conditions(self, role, low_side_temperature_C, high_side_temperature_C, plot=False):
        fluid_name = self._fluid_name(role)
        low_K = low_side_temperature_C + C_to_K
        high_K = high_side_temperature_C + C_to_K
        superheat = subcool = 3

        low_side_pressure = CP.PropsSI('P', 'T', low_K, 'Q', 0, fluid_name)
        high_side_pressure = CP.PropsSI('P', 'T', high_K, 'Q', 0, fluid_name)

        low_side_liquid_H = CP.PropsSI('H', 'T', low_K, 'Q', 0.2, fluid_name)          # expansion valve outlet
        low_side_vapor_H = CP.PropsSI('H', 'T', low_K + superheat, 'Q', 1, fluid_name)  # evaporator outlet
        high_side_liquid_H = CP.PropsSI('H', 'T', high_K + superheat, 'Q', 0, fluid_name)  # compressor outlet
        high_side_vapor_H = CP.PropsSI('H', 'T', high_K - subcool, 'Q', 1, fluid_name)     # condenser outlet

        # Order matches self._unit_operations[role]: [evaporator, compressor, condenser, expansion_valve]
        h_init = np.array([low_side_vapor_H, high_side_vapor_H, high_side_liquid_H, low_side_liquid_H])
        p_init = np.array([low_side_pressure, high_side_pressure, high_side_pressure, low_side_pressure])
        T_init = np.array([low_K + superheat, high_K + superheat, high_K - subcool, low_K])

        self._h_init[role] = h_init
        self._p_init[role] = p_init
        self._T_init[role] = T_init

        if plot:
            properties = getattr(self.model.fs, f"{role}_properties")
            properties.hp_diagram()
            plt.plot(h_init / 1000, p_init / 1000, 'ko')
            plt.title(f"{role} loop ({fluid_name})")
            plt.show()

    def initialize(self, verbose=False):
        for role in self.ROLES:
            self._initialize_loop(role, verbose=verbose)

    def _initialize_loop(self, role, verbose=False):
        fs = self.model.fs
        evaporator, compressor, condenser, expansion_valve = self._unit_operations[role]
        h_init, p_init, T_init = self._h_init[role], self._p_init[role], self._T_init[role]
        p_scale = 1

        ## Evaporator
        evaporator.inlet.flow_mass[0].fix(1)
        evaporator.inlet.pressure[0].fix(p_init[-1] * p_scale)
        if self.mode == Mode.PH:
            evaporator.inlet.enth_mass[0].fix(h_init[-1])
            evaporator.outlet.enth_mass[0].fix(h_init[0])
        else:
            evaporator.inlet.temperature[0].fix(T_init[-1])
            evaporator.outlet.temperature[0].fix(T_init[0])
        self.logger.info(f"Initializing {role} evaporator...")
        evaporator.initialize(outlvl=logging.WARNING)
        if verbose:
            evaporator.report()
        propagate_state(getattr(fs, f"{role}_evaporator_to_compressor"))

        ## Compressor
        compressor.inlet.pressure[0].fix(p_init[0] * p_scale)
        if self.mode == Mode.PH:
            compressor.inlet.enth_mass[0].fix(h_init[0])
        else:
            compressor.inlet.temperature[0].fix(T_init[0])
        compressor.outlet.pressure[0].fix(p_init[1] * p_scale)
        compressor.efficiency_isentropic[0].fix(self.fluids[role]["efficiency"])
        self.logger.info(f"Initializing {role} compressor...")
        compressor.initialize(outlvl=logging.WARNING)
        if verbose:
            compressor.report()
        propagate_state(getattr(fs, f"{role}_compressor_to_condenser"))

        ## Condenser
        condenser.inlet.pressure[0].fix(p_init[1] * p_scale)
        if self.mode == Mode.PH:
            condenser.inlet.enth_mass[0].fix(h_init[1])
            condenser.outlet.enth_mass[0].fix(h_init[2])
        else:
            condenser.inlet.temperature[0].fix(T_init[1])
            condenser.outlet.temperature[0].fix(T_init[2])
        self.logger.info(f"Initializing {role} condenser...")
        condenser.initialize(outlvl=logging.WARNING)
        if verbose:
            condenser.report()
        propagate_state(getattr(fs, f"{role}_condenser_to_expansion_valve"))

        ## Expansion valve
        expansion_valve.inlet.pressure[0].fix(p_init[2] * p_scale)
        expansion_valve.outlet.pressure[0].fix(p_init[3] * p_scale)
        self.logger.info(f"Initializing {role} expansion valve...")
        expansion_valve.initialize(outlvl=logging.WARNING)
        if verbose:
            expansion_valve.report()
        propagate_state(getattr(fs, f"{role}_expansion_valve_to_evaporator"))

        if verbose:
            print(f"\nFinished initializing {role} ({self._fluid_name(role)}) loop.")

    # ------------------------------------------------------------------
    # Specifications
    # ------------------------------------------------------------------

    def _derive_loop_sat_ranges(self, role, spec, cascade_T_window_C, cascade_approach_dT):
        """Return ((T_low_lo, T_low_hi), (T_high_lo, T_high_hi)) in degC -- the
        saturation-temperature ranges this loop's evaporator (low side) and
        condenser (high side) are expected to operate across.

        Which side is pinned by a known setpoint and which floats across the
        cascade window is decided by FACING (option b), explicitly -- never
        inferred from which kwargs happen to be present.
        """
        facing = self.FACING[role]
        w_lo, w_hi = min(cascade_T_window_C), max(cascade_T_window_C)

        # --- low side (evaporator) ---
        if facing["evaporator"] == "load":
            # Pinned by the cold-side delivery setpoint.
            T = spec.get("evap_sat_temperature")
            if T is None:
                T = self._T_init[role][0] - C_to_K
            low = (T, T)
        else:
            # Cascade-facing evaporator sits one approach BELOW the window.
            low = (w_lo - cascade_approach_dT, w_hi - cascade_approach_dT)

        # --- high side (condenser) ---
        if facing["condenser"] == "ambient":
            # Pinned by approach_constraint: T_sat = ambient + approach.
            amb, app = spec.get("ambient_temperature"), spec.get("condenser_approach")
            T = (amb + app) if (amb is not None and app is not None) \
                else (self._T_init[role][2] - C_to_K)
            high = (T, T)
        else:
            # Cascade-facing condenser sits IN the window.
            high = (w_lo, w_hi)

        return low, high

    def set_specifications(self, hot=None, cold=None, cascade_approach_dT=3.0,
                            cascade_T_window_C=(-15.0, 20.0), bound_margin_frac=0.30,
                            plr=None, cd=None, debug_disable_arc_pressure_eq=False):
        """hot / cold are dicts overriding `_LOOP_SPEC_DEFAULTS` for that
        loop (same keys as the single-stage file's set_specifications
        kwargs, minus flow/pressure-debug plumbing which is handled
        globally here). See module docstring for the fluids dict and
        `flow_fixed_role` for which loop's mass flow stays fixed."""
        fs = self.model.fs

        assert cascade_approach_dT > 0, "Cascade approach dT must be > 0"

        if plr is not None:
            assert 0.0 <= plr <= 1.0, "PLR must be in [0, 1]"
            self.plr = plr
        if cd is not None:
            assert 0.0 <= cd <= 1.0, "CD must be in [0, 1]"
            self.cd = cd
        fs.plr.set_value(self.plr)
        fs.cd.set_value(self.cd)

        specs = {"hot": hot or {}, "cold": cold or {}}
        self._derived_bounds = {}
        for role in self.ROLES:
            merged = dict(_LOOP_SPEC_DEFAULTS)
            merged.update(specs[role])
            if merged["superheating_max"] is None:
                merged["superheating_max"] = self.fluids[role]["superheat_max"]
            if merged["subcooling_max"] is None:
                merged["subcooling_max"] = self.fluids[role]["subcool_max"]

            # --- option (c): derive pressure bounds from THIS fluid's own
            # saturation curve over the range this loop actually operates in.
            coolprop_name = self.fluids[role]["coolprop_name"]
            low_range, high_range = self._derive_loop_sat_ranges(
                role, merged, cascade_T_window_C, cascade_approach_dT)

            if merged["low_side_pressure"] is None:
                p_min, p_max = derive_pressure_bounds(
                    coolprop_name, low_range[0], low_range[1], margin_frac=bound_margin_frac)
                merged["low_side_pressure"] = (p_min / 1000.0, p_max / 1000.0)
            if merged["high_side_pressure"] is None:
                p_min, p_max = derive_pressure_bounds(
                    coolprop_name, high_range[0], high_range[1], margin_frac=bound_margin_frac)
                merged["high_side_pressure"] = (p_min / 1000.0, p_max / 1000.0)

            self._derived_bounds[role] = {
                "fluid": coolprop_name,
                "evaporator_faces": self.FACING[role]["evaporator"],
                "condenser_faces": self.FACING[role]["condenser"],
                "low_sat_range_C": low_range,
                "high_sat_range_C": high_range,
                "low_side_pressure_kPa": merged["low_side_pressure"],
                "high_side_pressure_kPa": merged["high_side_pressure"],
            }

            self._set_loop_specifications(role, merged, debug_disable_arc_pressure_eq)

        fs.cascade_approach_dT.set_value(cascade_approach_dT)
        fs.cascade_energy_balance.activate()
        fs.cascade_approach_constraint.activate()

        calculate_scaling_factors(self.model)

    def _set_loop_specifications(self, role, spec, debug_disable_arc_pressure_eq):
        fs = self.model.fs
        evaporator, compressor, condenser, expansion_valve = self._unit_operations[role]

        low_side_pressure = spec["low_side_pressure"]
        high_side_pressure = spec["high_side_pressure"]
        evaporator_temperature = spec["evaporator_temperature"]
        compressor_temperature = spec["compressor_temperature"]
        condenser_temperature = spec["condenser_temperature"]
        expansion_valve_temperature = spec["expansion_valve_temperature"]
        subcooling = spec["subcooling"]
        superheating = spec["superheating"]
        subcooling_max = spec["subcooling_max"]
        superheating_max = spec["superheating_max"]
        max_pressure_ratio = spec["max_pressure_ratio"]
        ambient_temperature = spec["ambient_temperature"]
        condenser_approach = spec["condenser_approach"]
        evap_sat_temperature = spec["evap_sat_temperature"]

        assert superheating >= 0, "Superheating must be >= 0"
        assert subcooling >= 0, "Subcooling must be >= 0"
        assert max_pressure_ratio > 1.2, "Maximum pressure ratio must be > 1.2"

        ## Unfix everything from initialization
        for unit in self._unit_operations[role]:
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
            unit.inlet.pressure[0].unfix()
            unit.outlet.pressure[0].unfix()
            unit.T_lower_bound.deactivate()
            unit.T_upper_bound.deactivate()
            if self.mode in (Mode.ORIGINAL_TPX, Mode.IMPROVED_TPX):
                unit.inlet.vapor_frac[0].setlb(0)
                unit.inlet.vapor_frac[0].setub(1)

            # Clear any temperature/pressure bounds left over from a previous
            # set_specifications call on this same model (the homotopy/warm-start
            # path reuses one model across setpoints). Without this, a bound set
            # once would silently persist even after the caller stopped asking
            # for it -- a quiet way to reintroduce exactly the 30-50 C bug.
            for port in (unit.inlet, unit.outlet):
                port.pressure[0].setlb(None)
                port.pressure[0].setub(None)
                if self.mode != Mode.PH:
                    port.temperature[0].setlb(None)
                    port.temperature[0].setub(None)

        if self.mode == Mode.IMPROVED_TPX:
            compressor.vapor_constraint.activate()
        elif self.mode == Mode.ORIGINAL_TPX:
            compressor.vapor_constraint.deactivate()

        def check_input(bounds):
            return bounds is not None and len(bounds) == 2

        evap_sat_active = evap_sat_temperature is not None
        cond_sat_active = (ambient_temperature is not None) and (condenser_approach is not None)

        # Pressure bounds are ALWAYS the per-fluid derived (or explicitly
        # supplied) values from set_specifications -- there is deliberately no
        # hardcoded "wide safety" branch here any more. The old one clamped
        # every fluid to R134a's envelope (50-2000 / 100-5000 kPa) and made the
        # CO2 loop infeasible by construction. See module header.
        if check_input(low_side_pressure):
            low_side_pressure_min, low_side_pressure_max = low_side_pressure[0] * 1000, low_side_pressure[1] * 1000
        else:
            low_side_pressure_min = low_side_pressure_max = None
        if check_input(high_side_pressure):
            high_side_pressure_min, high_side_pressure_max = high_side_pressure[0] * 1000, high_side_pressure[1] * 1000
        else:
            high_side_pressure_min = high_side_pressure_max = None

        # Fix mass flow only for the designated role; leave the other free
        # so the solver can satisfy the cascade energy balance.
        if role == self.flow_fixed_role:
            evaporator.inlet.flow_mass[0].fix(1)

        ## Evaporator
        if low_side_pressure_min:
            evaporator.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            evaporator.inlet.pressure[0].setub(low_side_pressure_max)

        if check_input(evaporator_temperature):
            if evaporator_temperature[0] is not None:
                if self.mode == Mode.PH:
                    evaporator.Tmin.set_value(evaporator_temperature[0] + C_to_K)
                    evaporator.T_lower_bound.activate()
                else:
                    evaporator.outlet.temperature[0].setlb(evaporator_temperature[0] + C_to_K)
            if evaporator_temperature[1] is not None:
                if self.mode == Mode.PH:
                    evaporator.Tmax.set_value(evaporator_temperature[1] + C_to_K)
                    evaporator.T_upper_bound.activate()
                else:
                    evaporator.outlet.temperature[0].setub(evaporator_temperature[1] + C_to_K)

        if self.mode == Mode.ORIGINAL_TPX:
            evaporator.outlet.vapor_frac[0].setlb(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            evaporator.outlet.vapor_frac[0].fix(1.0)
            evaporator.control_volume.properties_out[0.0].eq_complementarity.deactivate()

        if superheating > 0.1:
            evaporator.superheating.set_value(superheating)
            evaporator.superheating_constraint.activate()
            evaporator.superheating_max.set_value(superheating_max)
            evaporator.superheating_upper_constraint.activate()
        else:
            evaporator.superheating_constraint.deactivate()
            evaporator.superheating_upper_constraint.deactivate()

        if evap_sat_temperature is not None:
            evaporator.T_sat_set.set_value(evap_sat_temperature + C_to_K)
            evaporator.evap_sat_constraint.activate()
        else:
            evaporator.evap_sat_constraint.deactivate()

        ## Compressor
        if low_side_pressure_min:
            compressor.inlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            compressor.inlet.pressure[0].setub(low_side_pressure_max)
        if high_side_pressure_min:
            compressor.outlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            compressor.outlet.pressure[0].setub(high_side_pressure_max)

        if check_input(compressor_temperature):
            if compressor_temperature[0] is not None:
                if self.mode == Mode.PH:
                    compressor.Tmin.set_value(compressor_temperature[0] + C_to_K)
                    compressor.T_lower_bound.activate()
                else:
                    compressor.outlet.temperature[0].setlb(compressor_temperature[0] + C_to_K)
            if compressor_temperature[1] is not None:
                if self.mode == Mode.PH:
                    compressor.Tmax.set_value(compressor_temperature[1] + C_to_K)
                    compressor.T_upper_bound.activate()
                else:
                    compressor.outlet.temperature[0].setub(compressor_temperature[1] + C_to_K)

        if self.mode == Mode.ORIGINAL_TPX:
            compressor.outlet.vapor_frac[0].setlb(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            compressor.outlet.vapor_frac[0].fix(1.0)
            compressor.control_volume.properties_out[0.0].eq_complementarity.deactivate()
            compressor.vapor_constraint.activate()

        compressor.ratioP.setub(max_pressure_ratio)
        compressor.ratioP.setlb(1.1)
        compressor.efficiency_isentropic[0].fix(self.fluids[role]["efficiency"])

        ## Condenser
        if high_side_pressure_min:
            condenser.inlet.pressure[0].setlb(high_side_pressure_min)
            condenser.outlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            condenser.inlet.pressure[0].setub(high_side_pressure_max)
            condenser.outlet.pressure[0].setub(high_side_pressure_max)

        # Option (b): a condenser outlet temperature bound is applied ONLY when
        # the caller explicitly asks for one. It is never synthesized from an
        # ambient-shaped default, and never gated on which other kwargs happen
        # to be set -- that inference is exactly what pinned the cascade-facing
        # CO2 condenser to 30-50 C. FACING records what each side really faces;
        # the physics constraints (approach_constraint for ambient-facing,
        # cascade_approach_constraint for cascade-facing) do the actual work.
        if check_input(condenser_temperature):
            if condenser_temperature[0] is not None:
                if self.mode == Mode.PH:
                    condenser.Tmin.set_value(condenser_temperature[0] + C_to_K)
                    condenser.T_lower_bound.activate()
                else:
                    condenser.outlet.temperature[0].setlb(condenser_temperature[0] + C_to_K)
            if condenser_temperature[1] is not None:
                if self.mode == Mode.PH:
                    condenser.Tmax.set_value(condenser_temperature[1] + C_to_K)
                    condenser.T_upper_bound.activate()
                else:
                    condenser.outlet.temperature[0].setub(condenser_temperature[1] + C_to_K)
        elif self.mode == Mode.PH:
            condenser.T_lower_bound.deactivate()
            condenser.T_upper_bound.deactivate()

        # Re-apply the generic Tmax safety cap (e.g. CO2), independent of
        # whatever condenser_temperature bounds were/weren't set above.
        Tmax_C = self.fluids[role].get("Tmax_C")
        if Tmax_C is not None:
            condenser.Tmax.set_value(Tmax_C + C_to_K)
            condenser.T_upper_bound.activate()

        if self.mode == Mode.ORIGINAL_TPX:
            condenser.outlet.vapor_frac[0].setub(0.01)
        elif self.mode == Mode.IMPROVED_TPX:
            condenser.outlet.vapor_frac[0].fix(0.0)
            condenser.control_volume.properties_out[0.0].eq_complementarity.deactivate()

        if subcooling > 0.1:
            condenser.subcooling.set_value(subcooling)
            condenser.subcooling_constraint.activate()
            condenser.subcooling_max.set_value(subcooling_max)
            condenser.subcooling_lower_constraint.activate()
        else:
            condenser.subcooling_constraint.deactivate()
            condenser.subcooling_lower_constraint.deactivate()

        if cond_sat_active:
            condenser.ambient_T.set_value(ambient_temperature + C_to_K)
            condenser.approach_T.set_value(condenser_approach)
            condenser.approach_constraint.activate()
        else:
            condenser.approach_constraint.deactivate()

        ## Expansion valve -- debug pressure-level override
        P_low_target = getattr(fs, f"{role}_P_low_target")
        P_high_target = getattr(fs, f"{role}_P_high_target")
        P_low = getattr(fs, f"{role}_P_low")
        P_high = getattr(fs, f"{role}_P_high")
        debug_names = ["P_low_target_constraint", "P_high_target_constraint",
                        "P_low_evap_in", "P_low_evap_out", "P_low_comp_in", "P_low_valve_out",
                        "P_high_comp_out", "P_high_cond_in", "P_high_cond_out", "P_high_valve_in"]
        arcs = [
            getattr(fs, f"{role}_evaporator_to_compressor_expanded"),
            getattr(fs, f"{role}_compressor_to_condenser_expanded"),
            getattr(fs, f"{role}_condenser_to_expansion_valve_expanded"),
            getattr(fs, f"{role}_expansion_valve_to_evaporator_expanded"),
        ]

        if debug_disable_arc_pressure_eq:
            if evap_sat_temperature is not None:
                P_low_val = CP.PropsSI('P', 'T', evap_sat_temperature + C_to_K, 'Q', 1, self._fluid_name(role))
                P_low_target.set_value(P_low_val)
                P_low.set_value(P_low_val)
            if cond_sat_active:
                T_cond_sat = ambient_temperature + condenser_approach + C_to_K
                P_high_val = CP.PropsSI('P', 'T', T_cond_sat, 'Q', 1, self._fluid_name(role))
                P_high_target.set_value(P_high_val)
                P_high.set_value(P_high_val)

            for arc in arcs:
                if hasattr(arc, "pressure_equality"):
                    arc.pressure_equality.deactivate()

            for name in debug_names:
                getattr(fs, f"{role}_{name}").activate()
            getattr(fs, f"{role}_P_low_evap_out").deactivate()
            getattr(fs, f"{role}_P_high_cond_out").deactivate()
        else:
            for arc in arcs:
                if hasattr(arc, "pressure_equality"):
                    arc.pressure_equality.activate()
            for name in debug_names:
                getattr(fs, f"{role}_{name}").deactivate()

        if low_side_pressure_min:
            expansion_valve.outlet.pressure[0].setlb(low_side_pressure_min)
        if low_side_pressure_max:
            expansion_valve.outlet.pressure[0].setub(low_side_pressure_max)
        if high_side_pressure_min:
            expansion_valve.inlet.pressure[0].setlb(high_side_pressure_min)
        if high_side_pressure_max:
            expansion_valve.inlet.pressure[0].setub(high_side_pressure_max)

        if check_input(expansion_valve_temperature):
            if expansion_valve_temperature[0] is not None:
                if self.mode == Mode.PH:
                    expansion_valve.Tmin.set_value(expansion_valve_temperature[0] + C_to_K)
                    expansion_valve.T_lower_bound.activate()
                else:
                    expansion_valve.outlet.temperature[0].setlb(expansion_valve_temperature[0] + C_to_K)
            if expansion_valve_temperature[1] is not None:
                if self.mode == Mode.PH:
                    expansion_valve.Tmax.set_value(expansion_valve_temperature[1] + C_to_K)
                    expansion_valve.T_upper_bound.activate()
                else:
                    expansion_valve.outlet.temperature[0].setub(expansion_valve_temperature[1] + C_to_K)

        if self.mode == Mode.ORIGINAL_TPX:
            expansion_valve.outlet.vapor_frac[0].setlb(0.01)
            expansion_valve.outlet.vapor_frac[0].setub(0.99)
        elif self.mode == Mode.IMPROVED_TPX:
            expansion_valve.control_volume.properties_out[0.0].eq_complementarity.deactivate()
            expansion_valve.control_volume.properties_out[0.0].eq_sat.activate()

    # ------------------------------------------------------------------
    # Solve
    # ------------------------------------------------------------------

    # @staticmethod: same reasoning as _compute_plf above -- everything this
    # needs (`solver`, `model`) is passed in explicitly, it never touches
    # `self`, so it doesn't need to know which CascadeCycle (if any) is
    # calling it. See the note above _compute_plf for what removing
    # @staticmethod would actually require/break.
    @staticmethod
    def _safe_solve(solver, model, tee=False):
        """Solve, never raise. Returns (termination_condition_str, ok).

        Pyomo raises ValueError("Cannot load a SolverResults object with bad
        status: error") when IPOPT exits with an error status, which would
        otherwise abort a whole ambient sweep on one bad point.
        """
        try:
            results = solver.solve(model, tee=tee)
            tc = str(results.solver.termination_condition)
            return tc, tc == "optimal"
        except Exception as exc:  # noqa: BLE001 -- deliberately broad
            return f"error: {type(exc).__name__}", False

    def optimize_COP(self, verbose=False, initialize=True, optimize=True):
        solver = get_solver()
        solver.options = {'max_iter': 1000, 'tol': 1e-6, 'linear_solver': 'ma57'}
        fs = self.model.fs

        if initialize:
            # Feasibility pre-solve (objective off) purely as a warm start.
            #
            # This is ADVISORY ONLY and must never abort the run. With the
            # objective deactivated the problem has several degrees of freedom
            # and no direction to pin them down, so IPOPT can wander into
            # "Restoration phase is called at point that is almost feasible"
            # and exit with an error status -- even when the point is feasible
            # to ~1e-13 and the subsequent objective solve converges cleanly.
            # Observed exactly this for R134a/CO2 at evap <= -28 C: the
            # pre-solve errored while the real solve reached optimal. So:
            # try it, keep whatever warm start it produced, and carry on.
            self.logger.info("Warm-start pre-solve (no objective)...")
            fs.compute_cop.deactivate()
            fs.obj.deactivate()
            tc, ok = self._safe_solve(solver, self.model, tee=verbose)
            if ok:
                self.logger.info("Warm-start pre-solve successful")
            else:
                self.logger.info(
                    "Warm-start pre-solve did not converge (%s) -- continuing to "
                    "the objective solve, which is the one that matters", tc)
            if verbose:
                fs.report()
            try:
                work_total = (value(fs.hot_compressor.work_mechanical[0])
                              + value(fs.cold_compressor.work_mechanical[0]))
                if work_total:
                    fs.cop.set_value(value(fs.cold_evaporator.heat_duty[0]) / work_total)
            except Exception:  # noqa: BLE001 -- warm start only
                pass

        self.logger.info("Setting up the optimization problem...")
        if optimize:
            fs.compute_cop.activate()
            fs.obj.activate()
        else:
            fs.compute_cop.deactivate()
            fs.obj.deactivate()

        tc, ok = self._safe_solve(solver, self.model, tee=verbose)

        if not optimize:
            try:
                work_total = (value(fs.hot_compressor.work_mechanical[0])
                              + value(fs.cold_compressor.work_mechanical[0]))
                if work_total:
                    fs.cop.set_value(value(fs.cold_evaporator.heat_duty[0]) / work_total)
            except Exception:  # noqa: BLE001
                pass

        # Retry from the same (failed) state, not a fresh restart -- IPOPT's
        # solve path isn't perfectly deterministic, so a repeat attempt from
        # the exact same point occasionally converges when the first didn't.
        # Stops early the moment one attempt succeeds (`if not ok` short-
        # circuits once `ok` is True).
        for _ in range(10):
            if not ok:
                tc, ok = self._safe_solve(solver, self.model, tee=verbose)

        if ok:
            self.logger.info("Optimization successful")
            self.logger.info("COP: {:.2f}".format(value(fs.cop)))
            optimization_converged = True
        else:
            self.logger.error("Optimization failed (%s)", tc)
            try:
                diag = DiagnosticsToolbox(self.model, constraint_residual_tolerance=1e-6)
                diag.display_constraints_with_large_residuals()
            except Exception:  # noqa: BLE001 -- diagnostics are best-effort
                self.logger.error("Could not run diagnostics on the failed model")
            optimization_converged = False

        if verbose:
            fs.report()

        self.optimization_converged = optimization_converged
        cop_full = value(fs.cop)
        plf = self._compute_plf(value(fs.plr), value(fs.cd))
        self._last_cop_full = cop_full
        self._last_plf = plf
        self._last_cop_part = plf * cop_full
        return value(fs.cop), optimization_converged

    def get_full_load_cop(self):
        return value(self.model.fs.cop)

    def get_part_load_cop(self):
        plf = self._compute_plf(value(self.model.fs.plr), value(self.model.fs.cd))
        return plf * value(self.model.fs.cop)

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def report_solution(self, role=None):
        """role=None reports both loops; pass "hot" or "cold" for one."""
        print("Optimized cascade COP:", round(value(self.model.fs.cop), 3))
        print("Hot-loop duty (evap):", value(self.model.fs.hot_evaporator.heat_duty[0]))
        print("Cold-loop duty (cond):", value(self.model.fs.cold_condenser.heat_duty[0]))

        roles = self.ROLES if role is None else (role,)
        for r in roles:
            self._report_loop_solution(r)

    def _report_loop_solution(self, role):
        properties = getattr(self.model.fs, f"{role}_properties")
        unit_operations = self._unit_operations[role]
        n = len(self._h_init[role])

        h_sol = np.zeros(n)
        p_sol = np.zeros(n)
        T_sol = np.zeros(n)
        S_sol = np.zeros(n)
        for i, unit in enumerate(unit_operations):
            h_sol[i] = unit.control_volume.properties_out[0].enth_mass()
            p_sol[i] = unit.outlet.pressure[0].value
            T_sol[i] = unit.control_volume.properties_out[0].temperature()
            S_sol[i] = unit.control_volume.properties_out[0].entr_mass()

        def add_warning():
            if self.optimization_converged is False:
                xlim = plt.gca().get_xlim()
                ylim = plt.gca().get_ylim()
                x = xlim[1] - (xlim[1] - xlim[0]) * 0.1
                y = ylim[0] + (ylim[1] - ylim[0]) * 0.1
                plt.text(x, y, "Warning: did not converge", color="red", fontsize=12,
                          bbox=dict(facecolor='white', alpha=0.8), va='bottom', ha='right')

        properties.hp_diagram()
        plt.plot(h_sol / 1000, p_sol / 1000, 'ko')
        plt.title(f"{role} loop ({self._fluid_name(role)})")
        add_warning()
        plt.show()

        properties.pt_diagram()
        plt.plot(T_sol, p_sol / 1000, 'ko')
        plt.title(f"{role} loop ({self._fluid_name(role)})")
        add_warning()
        plt.show()

        properties.ts_diagram()
        plt.plot(S_sol / 1000, T_sol, 'ko')
        plt.title(f"{role} loop ({self._fluid_name(role)})")
        add_warning()
        plt.show()

        for unit in unit_operations:
            unit.report()


if __name__ == "__main__":
    # Default fluids: hot=R134a, cold=CO2
    cycle = CascadeCycle()
    cycle.specify_initial_conditions()
    cycle.initialize(verbose=True)
    cycle.set_specifications()
    cop, converged = cycle.optimize_COP(verbose=True)
    print("Cascade COP:", cop, "| Converged:", converged)
