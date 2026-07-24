"""
phase4_spgp_sweep.py -- Phase 4, SPGP parameter set only: run the ideal
vapor-compression cycle (cubic PR EoS) across the full ambient sweep
(10/15/20/25 C) and write phase4_cop_SPGP.csv.

SPGP is the expected problem child (flagged since Phase 0): its own fitted
Pitzer acentric factor is NEGATIVE (omega=-0.2741), which makes the PR
alpha-function's kappa parameter negative too --

    kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2 ~= -0.068

-- unphysical for a real substance (kappa should always be > 0; a negative
value means the alpha(T) temperature-dependence term can behave non-
monotonically). This is a candidate STRUCTURAL problem with SPGP's own
critical-property fit, not just a numerics/tolerance issue like the ones
already fixed for NIST/GCGP this session.

Per project direction ("we will need to get SPGP to converge... they want
to know how bad the error is"), this script makes a genuine third attempt
per cell, beyond the standard two (normal spec, then the
debug_disable_arc_pressure_eq fallback): a much more permissive Ipopt
tolerance, specifically relaxing acceptable_dual_inf_tol (which defaults to
~1e10 / effectively disabled already for "acceptable" exits -- see the
2026-07-23 breadcrumb entry on why "acceptable" can mask a wrong branch;
here we deliberately want the OPPOSITE risk profile: get ANY number out,
even a low-confidence one, and quantify it, rather than reporting nothing).

Every cell's max_residual (largest raw constraint violation, computed
directly rather than trusting DiagnosticsToolbox's report formatting) is
recorded whether or not that cell "converged" -- this is the direct,
quantitative answer to "how bad is it" for any cell that doesn't reach a
clean optimal, including ones where the third attempt still doesn't
satisfy Ipopt's own convergence criteria but the underlying state is at
least inspectable.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
import os
from phase4_common import run_sweep, write_csv

METHOD = "SPGP"

# Third-attempt-only override: much more permissive than the standard
# relaxed options already baked into optimize_COP() (tol=1e-4,
# constr_viol_tol=1e-4, acceptable_tol=1e-3). This is a deliberate,
# labeled last resort for SPGP ONLY -- any cell that only converges under
# this override is flagged (fallback/third_attempt columns in the CSV) so
# it's never confused with a clean NIST/GCGP-style result.
THIRD_ATTEMPT_OPTIONS = {
    "tol": 1e-3,
    "constr_viol_tol": 1e-3,
    "acceptable_tol": 1e-2,
    "acceptable_dual_inf_tol": 1e6,
    "max_iter": 3000,
}

rows = run_sweep(METHOD, extra_solver_options=THIRD_ATTEMPT_OPTIONS)
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        f"phase4_cop_{METHOD}.csv")
write_csv(rows, out_path)

print(f"\n{'='*70}\n  SPGP NOTE\n{'='*70}")
print("Any COP reported here that needed the third-attempt override "
      "(third_attempt=True in the CSV) is LOW CONFIDENCE -- it means "
      "standard/relaxed Ipopt tolerances were not enough and we deliberately "
      "loosened acceptable_dual_inf_tol to force an exit. Check max_residual "
      "in the CSV before trusting any such value; a large max_residual means "
      "the reported 'COP' is not a real converged operating point, just "
      "whatever Ipopt was sitting on when it gave up.")
