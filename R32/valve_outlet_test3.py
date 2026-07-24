"""
valve_outlet_test3.py -- corrected fix for the expansion valve's OUTLET,
replacing the approach in valve_outlet_test2.py.

Why test2's approach was wrong: it independently, fully FIXED the outlet's
flow/composition/temperature (borrowing T and quality from the evaporator's
already-solved inlet), then called expansion_valve.initialize() afterward.
That's fine for a Heater (evaporator/condenser), because those units have a
free "heat_duty" variable with no constraint on it -- it just absorbs
whatever energy difference exists between an independently-chosen inlet and
outlet. The expansion valve has no such release valve: its adiabatic
PressureChanger assumption adds a HARD constraint, work[t] == 0, and models
no heat term at all, so the control-volume energy balance reduces to a
mandatory equality: inlet enthalpy == outlet enthalpy. By fixing the outlet
independently (via values borrowed from a different unit, the evaporator),
we gave the model no leftover freedom to satisfy that equality -- hence the
new "Too few degrees of freedom (rethrown)!" error once the full unit-level
expansion_valve.initialize() ran afterward and tried to enforce every
constraint (including the energy balance) with both ends already fully
pinned down.

Corrected approach: don't fix flow_mol/mole_frac_comp/temperature ourselves
at all -- let the unit's own initialize() handle those the normal way, same
as it always does. Only apply BOUNDS (not fixes) to temperature and
phase_frac["Vap"], fencing the solver away from the wrong branch while still
leaving the real isenthalpic energy-balance equation free to determine the
exact outlet state. This mirrors exactly what already worked for the valve's
INLET (a bound, not a fix, on phase_frac) -- just extended to also bound T,
since the outlet is genuinely two-phase and further from a safe default
guess than the inlet was.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode
import logging as pylog

vc = SimpleVaporCompressionCycle(
    "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX
)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)

# Run through the confirmed-working chain, failing (as expected) at the
# valve's outlet. We don't care about that failed attempt's leftover state --
# we're about to re-specify bounds ourselves and retry cleanly.
try:
    vc.initialize(verbose=False)
    print("vc.initialize() completed without raising (unexpected)\n")
except InitializationError as e:
    print(f"vc.initialize() failed as expected: {e}\n")

m = vc.model
evap_in = m.fs.evaporator.control_volume.properties_in[0]
valve_out = m.fs.expansion_valve.control_volume.properties_out[0]

T_target = value(evap_in.temperature)
quality_target = value(evap_in.phase_frac["Vap"])
print(f"Target (from evaporator inlet, same physical point): "
      f"T={T_target:.3f} K, phase_frac[Vap]={quality_target:.6f}")

# The earlier (failed) expansion_valve.initialize() attempt, run as part of
# vc.initialize() above, will have reactivated this outlet's sum_mole_frac_out
# as a side effect of its own internal bootstrapping (same behavior documented
# elsewhere in vapor_compression_cubic.py). It needs to stay ACTIVE this time,
# since we are NOT manually fixing mole_frac_comp ourselves -- reactivate/
# leave it active so composition is genuinely determined normally.
if not valve_out.sum_mole_frac_out.active:
    valve_out.sum_mole_frac_out.activate()

# Clean slate: the earlier (failed) expansion_valve.initialize() attempt
# inside vc.initialize() above ran the OLD recipe still baked into
# vapor_compression_cubic.py (fix flow_mol, fix mole_frac_comp, fix
# temperature at the end). Those fixes persist on this state block even
# though that attempt failed -- unfix them explicitly so we start from a
# genuinely clean, all-free state before applying our own bounds. Without
# this, DOF is wrongly negative before we've done anything ourselves.
valve_out.flow_mol.unfix()
valve_out.mole_frac_comp["R32"].unfix()
valve_out.temperature.unfix()
valve_out.phase_frac["Vap"].unfix()

print(f"Valve outlet DOF after clean-slate unfix: "
      f"{degrees_of_freedom(valve_out)} (expect 5: flow, T, P(fixed already), "
      f"comp, phase_frac -- P already fixed so this counts the other 4 free)")

# --- The corrected fix: BOUNDS only, no fixing of T or phase_frac. ---
valve_out.temperature.setlb(T_target - 5.0)
valve_out.temperature.setub(T_target + 5.0)
valve_out.temperature.set_value(T_target)  # warm start, not a fix

valve_out.phase_frac["Vap"].setlb(0.0)
valve_out.phase_frac["Vap"].setub(0.20)
valve_out.phase_frac["Vap"].set_value(quality_target)  # warm start, not a fix

print(f"\nValve outlet DOF before expansion_valve.initialize(): "
      f"{degrees_of_freedom(valve_out)} (informational only -- flow/composition "
      f"not yet fixed, that happens inside initialize() itself)")

print("\nRetrying expansion_valve.initialize() with T and phase_frac BOUNDED "
      "(not fixed)...")
try:
    m.fs.expansion_valve.initialize(
        outlvl=pylog.WARNING,
        state_args={
            "flow_mol": 1.0,
            "temperature": T_target,
            "pressure": value(m.fs.expansion_valve.outlet.pressure[0]),
            "mole_frac_comp": {"R32": 1.0},
        },
        optarg={"tol": 1e-4, "constr_viol_tol": 1e-4, "acceptable_tol": 1e-3},
    )
    print("expansion_valve.initialize() SUCCEEDED with bounds-only fix")
except Exception as e:
    print(f"expansion_valve.initialize() still failed: {type(e).__name__}: {e}")

# Remove the temporary bounds now that (hopefully) a real converged point has
# been found -- same pattern as the inlet's bound, removed once its job (as a
# fence during initialization only) is done.
valve_out.temperature.setlb(None)
valve_out.temperature.setub(None)
valve_out.phase_frac["Vap"].setlb(None)
valve_out.phase_frac["Vap"].setub(None)

print("\n=== Valve outlet after the corrected recipe ===")
T_out = value(valve_out.temperature)
tbub_out = value(valve_out.temperature_bubble["Vap", "Liq"])
vap_frac_out = value(valve_out.phase_frac["Vap"])
sum_frac = value(valve_out.phase_frac["Vap"]) + value(valve_out.phase_frac["Liq"])
print(f"  T = {T_out:.3f} K, tbub = {tbub_out:.3f} K")
print(f"  phase_frac[Vap] = {vap_frac_out:.6f}, sum of phase_frac = {sum_frac:.6f}")
print(f"  entr_mol = {value(valve_out.entr_mol):.4f} J/mol/K")

# Sanity check: does this state's enthalpy actually match the valve inlet's
# enthalpy (the real physical requirement for an adiabatic throttle)?
valve_in = m.fs.expansion_valve.control_volume.properties_in[0]
h_in = value(valve_in.enth_mol)
h_out = value(valve_out.enth_mol)
print(f"\n  valve inlet  enth_mol = {h_in:.4f} J/mol")
print(f"  valve outlet enth_mol = {h_out:.4f} J/mol")
print(f"  difference = {abs(h_in - h_out):.6f} J/mol "
      f"({'OK -- isenthalpic satisfied' if abs(h_in - h_out) < 1.0 else 'MISMATCH'})")

print(f"\nComparing against evaporator inlet target: "
      f"T diff = {abs(T_out - T_target):.4f} K, "
      f"quality diff = {abs(vap_frac_out - quality_target):.6f}")

print("\nPropagating valve outlet -> evaporator inlet...")
T_evap_before = value(evap_in.temperature)
vap_frac_evap_before = value(evap_in.phase_frac["Vap"])
propagate_state(m.fs.expansion_valve_to_evaporator)
T_evap_after = value(evap_in.temperature)
vap_frac_evap_after = value(evap_in.phase_frac["Vap"])
print(f"  Evaporator inlet T: {T_evap_before:.3f} -> {T_evap_after:.3f} K")
print(f"  Evaporator inlet phase_frac[Vap]: {vap_frac_evap_before:.6f} -> "
      f"{vap_frac_evap_after:.6f}")
