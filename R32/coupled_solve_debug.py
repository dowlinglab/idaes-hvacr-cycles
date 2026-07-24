"""
coupled_solve_debug.py -- diagnose why the real coupled solve
(optimize_COP(initialize=True, optimize=False)) fails at the actual
target spec (superheating=0, subcooling=0, per Shridhar 2016's ideal
cycle -- confirmed via PHASE3_NOTES.md and the source paper, NOT
something we can loosen).

vc.initialize() now succeeds end-to-end (valve outlet fix confirmed). But
the real coupled solve still fails ("infeasible" / "maxIterations"), with
the DiagnosticsToolbox flagging one large residual on the expansion_valve_
to_evaporator arc's mole_frac_comp_equality (~0.0125) plus several tiny
(~1e-5/1e-6) near-dome residuals elsewhere.

Hypothesis under test: Mode.IMPROVED_TPX pins phase_frac["Vap"] EXACTLY to
1.0 at the evaporator outlet and 0.0 at the condenser outlet (confirmed via
code inspection, lines ~879 and ~987). Physically, phase_frac=1 or 0 is
STILL a point ON the saturation dome (zero liquid / zero vapor, but still
a two-phase boundary point) -- the same Gibbs-phase-rule degeneracy (T,P
not independent) that caused every other bug this session. Unlike the
valve's outlet or the evaporator's own inlet, NEITHER of these two units'
outlets gets a tight T bound + warm-start anywhere in set_specifications()
-- so if Newton wanders to a trivial-root branch here (same failure mode
diagnosed multiple times already: hL=hV collapse, or T landing far from
Tsat), there's nothing fencing it in.

This script checks: right after set_specifications() (before solving),
are evaporator-outlet and condenser-outlet temperatures bounded/fixed at
all? Then attempts the solve and runs DiagnosticsToolbox on the full model
to see exactly which constraints/variables are the problem.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.solvers import get_solver
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

vc = SimpleVaporCompressionCycle("R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)
vc.initialize(verbose=False)
print("vc.initialize() succeeded.\n")

vc.set_specifications(
    ambient_temperature=20, condenser_approach=9, evap_sat_temperature=-29,
    superheating=0, subcooling=0, max_pressure_ratio=10,
)

m = vc.model
evap_out = m.fs.evaporator.control_volume.properties_out[0]
cond_out = m.fs.condenser.control_volume.properties_out[0]

print("=== State right after set_specifications(), BEFORE solving ===")
for name, blk in [("evaporator OUTLET", evap_out), ("condenser OUTLET", cond_out)]:
    T = blk.temperature
    pf = blk.phase_frac["Vap"]
    print(f"  {name}: T={value(T):.3f} K, fixed={T.fixed}, "
          f"lb={T.lb}, ub={T.ub}, phase_frac[Vap]={value(pf):.6f}, "
          f"fixed={pf.fixed}")

print(f"\nDOF(full model) = {degrees_of_freedom(m)}")

# Check the persistent, unchanging arc residual directly: is it a fixed-vs-
# fixed contradiction (comp never unfixed by set_specifications())?
valve_out = m.fs.expansion_valve.control_volume.properties_out[0]
evap_in = m.fs.evaporator.control_volume.properties_in[0]
vo_x = valve_out.mole_frac_comp["R32"]
ei_x = evap_in.mole_frac_comp["R32"]
print(f"\nvalve_out.mole_frac_comp['R32']: value={value(vo_x):.8f}, fixed={vo_x.fixed}")
print(f"evap_in.mole_frac_comp['R32']:   value={value(ei_x):.8f}, fixed={ei_x.fixed}")
print(f"raw difference: {abs(value(vo_x) - value(ei_x)):.8f}")

print("\n=== Attempting the coupled solve ===")
# Relaxed tolerance: the primal residuals are already near-zero at this
# point (confirmed in earlier runs) -- many log_mole_frac_tbub/tdew
# variables sit EXACTLY at 0.0 against a (None, 0) bound, which is
# expected/correct for a pure fluid (log(1.0)=0) but numerically degenerate
# right at the solution, tripping Ipopt's strict default dual-feasibility
# check even when the primal solution is fine. Same fix already confirmed
# for the condenser and valve outlet.
solver = get_solver(solver_options={"tol": 1e-4, "constr_viol_tol": 1e-4,
                                     "acceptable_tol": 1e-3})
res = solver.solve(m, tee=False)
print(f"Termination condition: {res.solver.termination_condition}\n")

print("=== State AFTER the (failed) solve ===")
for name, blk in [("evaporator OUTLET", evap_out), ("condenser OUTLET", cond_out)]:
    T = blk.temperature
    tbub = blk.temperature_bubble["Vap", "Liq"]
    pf = blk.phase_frac["Vap"]
    print(f"  {name}: T={value(T):.3f} K, tbub={value(tbub):.3f} K, "
          f"diff={value(T)-value(tbub):+.3f} K, phase_frac[Vap]={value(pf):.6f}")

print("\n=== DiagnosticsToolbox: constraints with large residuals ===")
dt = DiagnosticsToolbox(m)
try:
    dt.display_constraints_with_large_residuals()
except Exception as e:
    print(f"  (raised: {e})")

print("\n=== DiagnosticsToolbox: variables at or outside bounds ===")
try:
    dt.display_variables_at_or_outside_bounds()
except Exception as e:
    print(f"  (raised: {e})")
