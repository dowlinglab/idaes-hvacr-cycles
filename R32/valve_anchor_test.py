"""
valve_anchor_test.py -- test whether anchoring phase_frac["Vap"] near the
confirmed-correct value (1e-5, "almost all liquid") during the expansion
valve's inlet initialization fixes the wrong-branch convergence found in
valve_branch_debug.py (guess was correct, but the solve still drifted to
phase_frac[Vap] ~ 1.0, entropy mismatched from the condenser outlet feeding
it).

Same fix-then-release pattern already validated for the evaporator: fix the
phase variable at a known-good value to anchor the solver on the correct
branch, initialize, then unfix it afterward so downstream code sees a
normally free variable again.

Pass/fail check: does the valve inlet's entr_mol now match the condenser
outlet's entr_mol (the two are supposed to be the exact same physical
state)?

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from idaes.core.util.exceptions import InitializationError
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode
import logging as pylog

vc = SimpleVaporCompressionCycle(
    "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX
)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)

# Let this run through evaporator/compressor/condenser (all confirmed working)
# and fail at the valve, as we already know it does -- we don't care about
# that first failed attempt, we're about to override its leftover state.
try:
    vc.initialize(verbose=False)
    print("vc.initialize() completed without raising (unexpected at this point)\n")
except InitializationError as e:
    print(f"vc.initialize() failed as expected: {e}\n")

m = vc.model
cond_out = m.fs.condenser.control_volume.properties_out[0]
valve_in = m.fs.expansion_valve.control_volume.properties_in[0]

print("=== Before the anchored retry ===")
print(f"  condenser outlet entr_mol = {value(cond_out.entr_mol):.4f}")
print(f"  valve inlet   entr_mol = {value(valve_in.entr_mol):.4f}  "
      f"(phase_frac[Vap] = {value(valve_in.phase_frac['Vap']):.6f}, WRONG BRANCH)")

# --- The fix under test: BOUND (not fix) phase_frac near the correct value.
# Fixing it outright over-specifies the state (T, P, flow, composition are
# ALREADY fully fixed here -- confirmed via the first attempt's DOF=-1 error).
# A bound restricts the range without consuming a degree of freedom, fencing
# the solver away from the wrong (~1.0) branch without forcing an exact,
# thermodynamically-inappropriate value (this state isn't on the dome, so
# there's no real "correct nonzero quality" to fix to in the first place).
valve_in.phase_frac["Vap"].setub(0.01)

print("\nRetrying expansion_valve.initialize() with phase_frac bounded (not fixed)...")
try:
    m.fs.expansion_valve.initialize(
        outlvl=pylog.WARNING,
        optarg={"tol": 1e-4, "constr_viol_tol": 1e-4, "acceptable_tol": 1e-3},
    )
    print("expansion_valve.initialize() SUCCEEDED with phase_frac bounded")
except Exception as e:
    print(f"expansion_valve.initialize() still failed: {type(e).__name__}: {e}")

# Remove the temporary bound now that (hopefully) a real converged point has
# been found -- it shouldn't clip the real operating range once the full
# cycle is solved.
valve_in.phase_frac["Vap"].setub(None)

print("\n=== After the anchored retry ===")
entr_valve = value(valve_in.entr_mol)
entr_cond = value(cond_out.entr_mol)
print(f"  condenser outlet entr_mol = {entr_cond:.4f}")
print(f"  valve inlet   entr_mol = {entr_valve:.4f}")
print(f"  valve inlet   phase_frac[Vap] = {value(valve_in.phase_frac['Vap']):.6f}")
print(f"  difference = {abs(entr_valve - entr_cond):.6f}")

print("\n=== Verdict ===")
if abs(entr_valve - entr_cond) < 0.01:
    print("  MATCH -- the anchored retry converged to the same physical state")
    print("  as the condenser outlet. Fix confirmed.")
else:
    print("  STILL MISMATCHED -- anchoring phase_frac alone was not enough.")
    print("  Need to dig further (e.g. check if T also needs anchoring/bounding).")
