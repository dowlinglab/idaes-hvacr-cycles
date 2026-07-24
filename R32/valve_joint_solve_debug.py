"""
valve_joint_solve_debug.py -- stop guessing bounds, get Ipopt to show us
exactly which constraint is violated.

Two hypotheses have now been tried and ruled out for the "locally infeasible"
failure inside expansion_valve.initialize()'s joint (whole-unit) solve:
  1. phase_frac bound too tight (0.20) -- widening to 0.40 made things WORSE
     (diverged to nonsense: tdew defaulting to 450K, phase_frac summing to
     0.01), not better. Reverted.
  2. sum_mole_frac_out being reactivated -- turned out to already be
     deactivated (from __init__) by the time this code runs for the first
     time in a real, single-pass vc.initialize() call. Reverting this made
     no difference, confirming it was never live.

Instead of a third guess, this script reproduces the exact sequence
vapor_compression_cubic.py now does (valve_in fixed correctly, valve_out
bounded on T and phase_frac, NOT fixed) but replaces the opaque
expansion_valve.initialize() call with our OWN direct, visible steps:
  1. Do the same local warm-start init on properties_out only (mirrors what
     init_adiabatic() does internally).
  2. Solve the WHOLE unit (control volume + both state blocks) ourselves
     with tee=True, so Ipopt's full iteration log prints directly --
     showing constraint violations, not just a final "infeasible" verdict.
  3. Run IDAES's own DiagnosticsToolbox on the unit afterward to explicitly
     list which constraints have large residuals / are infeasible.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from idaes.core.util.exceptions import InitializationError
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode
import logging as pylog

vc = SimpleVaporCompressionCycle(
    "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX
)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)

# Run through the real initialize() sequence, which now fails inside
# expansion_valve.initialize()'s own joint solve. We don't care about that
# failure -- what matters is the state it leaves behind: valve_in correctly
# fixed, valve_out bounded (not fixed) on T and phase_frac.
try:
    vc.initialize(verbose=False)
    print("vc.initialize() completed without raising (unexpected)\n")
except InitializationError as e:
    print(f"vc.initialize() failed as expected: {e}\n")

m = vc.model
valve = m.fs.expansion_valve
valve_in = valve.control_volume.properties_in[0]
valve_out = valve.control_volume.properties_out[0]

print("=== State entering this script ===")
print(f"  valve_in : T={value(valve_in.temperature):.3f} K, "
      f"P={value(valve_in.pressure):.1f} Pa, "
      f"enth_mol={value(valve_in.enth_mol):.4f} J/mol, "
      f"fixed(T)={valve_in.temperature.fixed}, "
      f"fixed(flow)={valve_in.flow_mol.fixed}, "
      f"fixed(x[R32])={valve_in.mole_frac_comp['R32'].fixed}")
print(f"  valve_out: T={value(valve_out.temperature):.3f} K "
      f"(bounds={valve_out.temperature.lb},{valve_out.temperature.ub}), "
      f"P={value(valve_out.pressure):.1f} Pa, "
      f"enth_mol={value(valve_out.enth_mol):.4f} J/mol, "
      f"phase_frac[Vap]={value(valve_out.phase_frac['Vap']):.6f} "
      f"(bounds={valve_out.phase_frac['Vap'].lb},{valve_out.phase_frac['Vap'].ub}), "
      f"fixed(T)={valve_out.temperature.fixed}, "
      f"fixed(flow)={valve_out.flow_mol.fixed}, "
      f"fixed(x[R32])={valve_out.mole_frac_comp['R32'].fixed}, "
      f"sum_mole_frac_out.active={valve_out.sum_mole_frac_out.active}")

print(f"\n  Isenthalpic gap right now: "
      f"{abs(value(valve_in.enth_mol) - value(valve_out.enth_mol)):.4f} J/mol "
      f"(this is what the joint solve needs to close to ~0)")

print(f"\n  DOF(expansion_valve unit) = {degrees_of_freedom(valve)}")

print("\n=== Re-solving the WHOLE unit directly, with tee=True ===")
print("(Watch for 'infeasible' constraint names / large residuals in the "
      "Ipopt log below)\n")
res = get_solver().solve(valve, tee=True)
print(f"\nTermination condition: {res.solver.termination_condition}")

print("\n=== DiagnosticsToolbox: constraints with large residuals ===")
dt = DiagnosticsToolbox(valve)
try:
    dt.display_constraints_with_large_residuals()
except Exception as e:
    print(f"  (display_constraints_with_large_residuals raised: {e})")

print("\n=== DiagnosticsToolbox: variables at or near their bounds ===")
try:
    dt.display_variables_at_or_outside_bounds()
except Exception as e:
    print(f"  (display_variables_at_or_outside_bounds raised: {e})")

print("\n=== After the direct solve ===")
print(f"  valve_out: T={value(valve_out.temperature):.3f} K, "
      f"phase_frac[Vap]={value(valve_out.phase_frac['Vap']):.6f}, "
      f"enth_mol={value(valve_out.enth_mol):.4f} J/mol")
print(f"  isenthalpic gap: "
      f"{abs(value(valve_in.enth_mol) - value(valve_out.enth_mol)):.4f} J/mol")
