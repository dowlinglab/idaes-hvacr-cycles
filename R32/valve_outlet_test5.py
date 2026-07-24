"""
valve_outlet_test5.py -- replace reliance on expansion_valve.initialize()'s
built-in joint (whole-unit) solve with our OWN direct, explicit isenthalpic
solve of just the outlet.

Diagnosed via valve_joint_solve_debug.py (tee=True + DiagnosticsToolbox):
the built-in joint solve (`opt.solve(blk)` inside IDAES's init_adiabatic())
re-solves the ENTIRE unit -- inlet and outlet together. Even though the
inlet's 4 canonical state vars (flow/T/P/comp) are fixed, its DERIVED flash
variables (phase-specific flows, log mole fractions) are still free, and a
badly-guessed outlet warm start was dragging them into a bad restoration
spiral (50 Ipopt iterations, mostly restoration-phase) ending in "locally
infeasible" -- confirmed by a huge residual specifically on the INLET's own
component_flow_balances constraint, not anything on the outlet. Two rounds
of bound-tweaking (phase_frac ceiling 0.20 -> 0.40, sum_mole_frac_out
activation) were both ruled out as the cause.

Fix under test here: don't touch the inlet at all. Solve the outlet's own
state block directly, in isolation, with the ACTUAL physical requirement
(inlet enthalpy == outlet enthalpy) as an explicit constraint -- using the
valve's own already-correctly-solved inlet enthalpy as the target, not a
value borrowed from the evaporator. This is a clean, well-posed 1-variable
root-find (T free, phase_frac genuinely free/internally-determined by the
block's own flash equations, one added equation) that never touches the
inlet's internals at all -- should avoid the joint-solve fragility entirely.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value, Constraint
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.initialization import propagate_state
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode
import logging as pylog

vc = SimpleVaporCompressionCycle(
    "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX
)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)

# Run through the real (currently still-failing) sequence -- it fails inside
# expansion_valve.initialize()'s joint solve, but properties_in.initialize()
# (called internally, before the joint solve) already succeeds every time,
# leaving valve_in correctly solved. We just need that part.
try:
    vc.initialize(verbose=False)
    print("vc.initialize() completed without raising (unexpected)\n")
except InitializationError as e:
    print(f"vc.initialize() failed as expected: {e}\n")

m = vc.model
valve = m.fs.expansion_valve
valve_in = valve.control_volume.properties_in[0]
valve_out = valve.control_volume.properties_out[0]

h_target = value(valve_in.enth_mol)
flow_target = value(valve_in.flow_mol)
comp_target = value(valve_in.mole_frac_comp["R32"])
print(f"Valve inlet (already correctly solved): T={value(valve_in.temperature):.3f} K, "
      f"enth_mol={h_target:.4f} J/mol, flow_mol={flow_target:.6f}, "
      f"mole_frac_comp[R32]={comp_target:.6f}")

# --- Direct, isolated solve of the outlet, bypassing the joint unit solve ---
valve_out.sum_mole_frac_out.deactivate()  # Task #24 pattern -- always off,
# regardless of current status (deactivating an already-inactive constraint
# is a harmless no-op; the earlier failed attempt likely reactivated it).

valve_out.flow_mol.unfix()
valve_out.mole_frac_comp["R32"].unfix()
valve_out.temperature.unfix()

valve_out.flow_mol.fix(flow_target)
valve_out.mole_frac_comp["R32"].fix(comp_target)
# pressure already fixed elsewhere (expansion_valve.outlet.pressure)

T_guess = value(m.fs.evaporator.control_volume.properties_in[0].temperature)
valve_out.temperature.setlb(T_guess - 10.0)
valve_out.temperature.setub(T_guess + 10.0)
valve_out.temperature.set_value(T_guess)

# The actual physical requirement for an adiabatic throttle: enthalpy is
# conserved. Add it as an explicit, temporary constraint using the valve's
# OWN actual inlet enthalpy (not a value borrowed from another unit).
valve_out.isenthalpic_seed_con = Constraint(expr=valve_out.enth_mol == h_target)

dof = degrees_of_freedom(valve_out)
print(f"\nValve outlet DOF before direct isenthalpic solve: {dof}")
assert dof == 0, f"expected DOF=0, got {dof}"

# Relaxed tolerance (same fix already confirmed for the condenser earlier
# tonight): the diagnostic above showed only ONE tiny residual (1.2e-05,
# barely over the 1e-05 threshold) with enthalpy matching EXACTLY -- this is
# a near-converged point tripping Ipopt's strict default tolerance, not a
# wrong branch.
solver = get_solver(solver_options={"tol": 1e-4, "constr_viol_tol": 1e-4,
                                     "acceptable_tol": 1e-3})
res = solver.solve(valve_out, tee=False)
print(f"Direct isenthalpic solve: {res.solver.termination_condition}")

T_out = value(valve_out.temperature)
tbub_out = value(valve_out.temperature_bubble["Vap", "Liq"])
vap_frac_out = value(valve_out.phase_frac["Vap"])
h_out = value(valve_out.enth_mol)
print(f"\n=== Valve outlet after direct isenthalpic solve ===")
print(f"  T = {T_out:.3f} K, tbub = {tbub_out:.3f} K")
print(f"  phase_frac[Vap] = {vap_frac_out:.6f}")
print(f"  enth_mol = {h_out:.4f} J/mol  (target was {h_target:.4f})")
print(f"  isenthalpic gap = {abs(h_out - h_target):.6f} J/mol")

print("\n=== DiagnosticsToolbox: constraints with large residuals "
      "(BEFORE cleanup) ===")
dt = DiagnosticsToolbox(valve_out)
try:
    dt.display_constraints_with_large_residuals()
except Exception as e:
    print(f"  (raised: {e})")

print("\n=== DiagnosticsToolbox: variables at or outside bounds "
      "(BEFORE cleanup) ===")
try:
    dt.display_variables_at_or_outside_bounds()
except Exception as e:
    print(f"  (raised: {e})")

# Clean up: remove the temporary constraint, revert T to the standard fixed
# pattern (matches how every other unit's state ends up before the real
# coupled solve; set_specifications() unfixes things again later as needed).
valve_out.del_component(valve_out.isenthalpic_seed_con)
valve_out.temperature.setlb(None)
valve_out.temperature.setub(None)
valve_out.temperature.fix(T_out)

print("\nPropagating valve outlet -> evaporator inlet (sanity check only)...")
evap_in = m.fs.evaporator.control_volume.properties_in[0]
T_evap_before = value(evap_in.temperature)
propagate_state(m.fs.expansion_valve_to_evaporator)
T_evap_after = value(evap_in.temperature)
print(f"  Evaporator inlet T: {T_evap_before:.3f} -> {T_evap_after:.3f} K "
      f"(should be unchanged -- already fixed)")

print("\n=== Verdict ===")
if abs(h_out - h_target) < 0.01 and res.solver.termination_condition.name == "optimal":
    print("  CONFIRMED: direct isenthalpic solve of the outlet alone works "
          "cleanly, bypassing the fragile joint unit-level solve entirely.")
else:
    print("  Still not clean -- needs further investigation.")
