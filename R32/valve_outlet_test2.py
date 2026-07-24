"""
valve_outlet_test2.py -- apply the evaporator's FULL validated two-phase
recipe (free T with a tight bound, warm-start, fix phase_frac as a numerical
seed, solve directly, then revert) to the expansion valve's OUTLET, since
it's physically the same TYPE of state as the evaporator's inlet (genuinely
two-phase, sitting on the dome) -- not the same type of problem as the
valve's INLET (genuinely single-phase, just needed a bound to avoid a wrong
branch).

Also checks the follow-up question: after propagate_state() copies the
valve outlet's values into the evaporator's inlet, does the evaporator
inlet's phase_frac (which propagate_state() should NOT touch -- it only
copies flow/T/P/composition) actually stay unchanged, or does something
unexpected happen?

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.initialization import propagate_state
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode
import logging as pylog

vc = SimpleVaporCompressionCycle(
    "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX
)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)

# Run through the confirmed-working chain, failing (as expected) at the
# valve's outlet -- we don't care about that failed attempt's leftover
# state, we're about to explicitly re-specify everything ourselves.
try:
    vc.initialize(verbose=False)
    print("vc.initialize() completed without raising (unexpected)\n")
except InitializationError as e:
    print(f"vc.initialize() failed as expected: {e}\n")

m = vc.model
evap_in = m.fs.evaporator.control_volume.properties_in[0]
valve_out = m.fs.expansion_valve.control_volume.properties_out[0]

# Record the evaporator inlet's CURRENT state before we touch anything else,
# so we can check afterward whether propagate_state() disturbed it.
T_evap_before = value(evap_in.temperature)
vap_frac_evap_before = value(evap_in.phase_frac["Vap"])
print(f"Evaporator inlet BEFORE touching the valve outlet: "
      f"T={T_evap_before:.3f} K, phase_frac[Vap]={vap_frac_evap_before:.6f}")

# --- Evaporator-style two-phase recipe, applied to the valve's OUTLET ---
T_target = T_evap_before  # same physical point -- reuse the known-good value
quality_target = vap_frac_evap_before  # same physical point -- reuse the known-good quality

# The earlier (failed) expansion_valve.initialize() call, run as part of
# vc.initialize() just above, already reactivated this outlet's
# sum_mole_frac_out constraint as a side effect of its own internal
# bootstrapping (the exact same behavior documented in vapor_compression_
# cubic.py for why this deactivation needs reasserting right before the
# real solve). Redo it here, or fixing mole_frac_comp ourselves next would
# double-state the same fact and cause the same DOF=-1 bug we already fixed
# once this session (Task #24).
valve_out.sum_mole_frac_out.deactivate()

valve_out.flow_mol.fix(1.0)
valve_out.mole_frac_comp["R32"].fix(1.0)
# pressure is already fixed (line 615 in vapor_compression_cubic.py)

valve_out.temperature.unfix()
valve_out.temperature.setlb(T_target - 5.0)
valve_out.temperature.setub(T_target + 5.0)
valve_out.temperature.set_value(T_target)
valve_out.phase_frac["Vap"].fix(quality_target)

dof = degrees_of_freedom(valve_out)
print(f"\nValve outlet DOF before direct solve: {dof}")
assert dof == 0, f"expected DOF=0, got {dof}"

res = get_solver().solve(valve_out)
print(f"Valve outlet direct 2-phase solve: {res.solver.termination_condition}")

# Revert to the standard fixed-T pattern, same as the evaporator's own recipe
valve_out.phase_frac["Vap"].unfix()
valve_out.temperature.setlb(None)
valve_out.temperature.setub(None)
valve_out.temperature.fix(value(valve_out.temperature))

T_out = value(valve_out.temperature)
tbub_out = value(valve_out.temperature_bubble["Vap", "Liq"])
vap_frac_out = value(valve_out.phase_frac["Vap"])
print(f"\n=== Valve outlet after the recipe ===")
print(f"  T = {T_out:.3f} K, tbub = {tbub_out:.3f} K")
print(f"  phase_frac[Vap] = {vap_frac_out:.6f}")
print(f"  entr_mol = {value(valve_out.entr_mol):.4f} J/mol/K")

# --- Now check the follow-up question: does propagating into the
# evaporator inlet disturb its phase_frac? ---
print(f"\nPropagating valve outlet -> evaporator inlet...")
propagate_state(m.fs.expansion_valve_to_evaporator)

T_evap_after = value(evap_in.temperature)
vap_frac_evap_after = value(evap_in.phase_frac["Vap"])
print(f"\n=== Evaporator inlet AFTER propagation ===")
print(f"  T:            before={T_evap_before:.3f} K -> after={T_evap_after:.3f} K "
      f"(changed by {T_evap_after - T_evap_before:+.3f} K)")
print(f"  phase_frac[Vap]: before={vap_frac_evap_before:.6f} -> after={vap_frac_evap_after:.6f} "
      f"(changed by {vap_frac_evap_after - vap_frac_evap_before:+.6f})")

print("\n=== Verdict ===")
if abs(vap_frac_evap_after - vap_frac_evap_before) < 1e-9:
    print("  CONFIRMED: phase_frac on the evaporator inlet is untouched by")
    print("  propagate_state(), as expected (it only copies flow/T/P/composition).")
else:
    print("  UNEXPECTED: phase_frac changed. Something other than what we")
    print("  assumed is touching it -- worth investigating further.")
if abs(T_evap_after - T_evap_before) < 5.0:
    print("  Temperature shift is small/reasonable, not a big disruptive jump.")
else:
    print("  Temperature shift is LARGE -- worth checking if this is a problem.")
