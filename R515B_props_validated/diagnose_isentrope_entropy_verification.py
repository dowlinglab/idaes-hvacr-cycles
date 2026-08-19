"""
Independent post-hoc verification: does every point our isentrope-walking
code emits actually have the entropy it claims to have?

Why this check, why now
------------------------
The R1234yf sister project (separate pure-fluid property package, same lab)
documented a bug in its own isentrope solver's BREADCRUMB.md: very close to
the critical point, `sol.success=True` from a least_squares/root call was
being trusted as proof of a correct answer, when the solved state was
actually off by ~100 J/(kg K) from the intended target entropy. The failure
mode: near the critical point, dP/drho -> 0 and the residual surface can
have a shallow valley the solver reports as "converged" even though it
hasn't actually zeroed the entropy residual -- `sol.success` only means
the underlying numerical method terminated normally, not that it found the
right root.

This project's own isentrope-walking code (compute_isentrope_liquid_side,
compute_isentrope_vapor_side, and the extension helper) calls
`scipy.optimize.root(_isentrope_2eq_residual, ..., method="hybr")` at every
step and gates acceptance ONLY on `sol.success` (see mixture_isentrope_
validation.py lines ~2733-2745, 2872-2873, 2972-2976) -- structurally the
SAME pattern R1234yf's Bug #3 exploited. By contrast, this project's own
bubble/dew VLE solver (solve_bubble_at_t / solve_dew_at_t, lines ~1198-1217,
~1390-1398) already does NOT trust sol.success alone -- it re-evaluates the
converged state and explicitly checks r_p<=1e-6 and r_mu<=1e-6 before
accepting. So there's a real asymmetry in our own codebase: one solver
already has the discipline R1234yf's bug taught them, the other doesn't.

This script checks, empirically, whether that gap is just theoretical or
has actually produced silently-wrong points in our real output -- exactly
the way R1234yf's own team diagnosed their bug (checking a specific
suspect entropy directly rather than assuming).

Method
------
Runs the REAL pipeline (run_true_vle_envelope + compute_isentrope_liquid_
side + compute_isentrope_vapor_side) at the project's standard parameters
(w1=0.911, Tmin=255.4K, Tmax=349.8K, n=80) -- IDENTICAL to what
mixture_isentrope_validation.py's own CLI does, so this checks the actual
production output, not a synthetic case.

Each solver call only returns (T_K, P_Pa, h_Jmol) per point -- rho is
solved internally but not exposed. So for each returned (T,P) point, this
script independently re-solves for rho via a SEPARATE damped-Newton search
(same method already vetted in validate_against_kang2024.py / validate_
pure_components_against_kang2024.py for this exact EOS), seeded by walking
along each isentrope's own point sequence (continuation from the previous
point, exactly mirroring how the original solve was seeded) rather than a
single global bracket -- deliberately avoiding the wide-brentq root-
collision failure mode already documented in this project for exactly this
kind of near-critical density solve.

Having independently recovered rho at each point, it computes the ACTUAL
entropy there via _mix_entropy_direct and compares to the s_target the
point was supposed to represent. Any point whose independently-recovered
entropy deviates from s_target by more than a strict tolerance is flagged
-- this is the direct, non-circular analogue of R1234yf's "checked directly
for s=1475... off by ~100 J/(kg K) even though sol.success said converged."

Read-only against the model: only imports existing functions, does not
modify mixture_isentrope_validation.py.
"""

import sys
from pathlib import Path

import numpy as np

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    mix_state,
    _dp_drho_fd,
    _mix_entropy_direct,
    run_true_vle_envelope,
    compute_isentrope_liquid_side,
    compute_isentrope_vapor_side,
    ISENTROPE_VALUES_BTU_LBMR,
    BTU_LBMR_TO_JKGK,
)

W1 = 0.911
FLUID1, FLUID2 = "r1234ze", "r227ea"
T_MIN, T_MAX, N = 255.4, 349.8, 80

# Independent verification tolerance: relative entropy mismatch above this
# is flagged as a silent solver failure (R1234yf's bug was ~100 J/(kg K)
# out of an ~1500-2000 J/(kg K) scale, i.e. several percent -- this
# tolerance is deliberately far stricter, matching the model's own r_s
# nominal solve tolerance of 1e-10 relative, with generous margin for
# independent-solve numerical noise).
ENTROPY_MISMATCH_REL_TOL = 1.0e-4


def solve_rho_at_TP(d1, d2, z1, t_k, p_target_pa, rho_seed_mol,
                     max_iter=80, tol_rel=1e-11):
    """Damped-Newton molar-density solve at fixed (T,P), seeded from a
    trusted nearby value (continuation), same method already vetted in
    validate_against_kang2024.py. Returns None on non-convergence rather
    than raising, so the caller can flag it as a broken continuation
    rather than crash the whole sweep."""
    rho_mol = rho_seed_mol
    for _ in range(max_iter):
        st = mix_state(d1, d2, t_k, rho_mol, z1)
        resid = st.p_pa - p_target_pa
        if abs(resid) < tol_rel * max(1.0, abs(p_target_pa)):
            return rho_mol
        dp_drho = _dp_drho_fd(d1, d2, t_k, rho_mol, z1)
        if dp_drho == 0 or not np.isfinite(dp_drho):
            return None
        step = -resid / dp_drho
        step = np.clip(step, -0.2 * rho_mol, 0.2 * rho_mol)
        rho_mol += step
        if rho_mol <= 0:
            return None
    return None


def verify_branch(label, d1, d2, z1, isentrope_dict, seed_mode):
    """seed_mode: 'liquid' (start dense, ~1000 kg/m3-ish in molar units via
    a generic dense guess) or 'vapor' (start ideal-gas-dilute)."""
    n_checked = 0
    n_flagged = 0
    n_seed_fail = 0
    flagged_detail = []
    for s_target, points in isentrope_dict.items():
        if not points:
            continue
        # Seed the FIRST point of this isentrope's sequence.
        t0, p0 = points[0]["T_K"], points[0]["P_Pa"]
        if seed_mode == "liquid":
            rho_seed = 10000.0  # mol/m^3, generic dense-liquid-ish starting guess
        else:
            rho_seed = p0 / (8.314 * t0)  # ideal-gas guess, mol/m^3
        rho_prev = None
        for i, pt in enumerate(points):
            t_k, p_pa = pt["T_K"], pt["P_Pa"]
            seed = rho_prev if rho_prev is not None else rho_seed
            rho = solve_rho_at_TP(d1, d2, z1, t_k, p_pa, seed)
            if rho is None and rho_prev is not None:
                # continuation lost -- retry once from the generic seed
                # before giving up on this point (branch-jump guard)
                rho = solve_rho_at_TP(d1, d2, z1, t_k, p_pa, rho_seed)
            if rho is None:
                n_seed_fail += 1
                rho_prev = None
                continue
            s_actual = _mix_entropy_direct(d1, d2, t_k, rho, z1)
            rel_err = abs(s_actual - s_target) / max(1.0, abs(s_target))
            n_checked += 1
            if rel_err > ENTROPY_MISMATCH_REL_TOL:
                n_flagged += 1
                flagged_detail.append({
                    "branch": label, "s_target_Jmol_K": s_target,
                    "point_index": i, "T_K": t_k, "P_Pa": p_pa,
                    "s_actual_Jmol_K": s_actual, "rel_err": rel_err,
                })
            rho_prev = rho
    return n_checked, n_flagged, n_seed_fail, flagged_detail


def main():
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1, mw2 = mw_from_json(d1), mw_from_json(d2)

    t_vals = np.linspace(T_MIN, T_MAX, N)
    print(f"Running real pipeline: w1={W1}, Tmin={T_MIN}, Tmax={T_MAX}, n={N} ...")
    bubble_rows, dew_rows, z1, crit_point = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals)
    mw_mix = z1 * mw1 + (1.0 - z1) * mw2

    s_values_jmolK = [s_btu * BTU_LBMR_TO_JKGK * mw_mix for s_btu in ISENTROPE_VALUES_BTU_LBMR]

    print("Computing isentropes (liquid side + vapor side), exactly as the real CLI does...")
    isentropes_liquid = compute_isentrope_liquid_side(d1, d2, z1, bubble_rows, s_values_jmolK, crit_point=crit_point)
    isentropes_vapor = compute_isentrope_vapor_side(d1, d2, z1, dew_rows, s_values_jmolK, crit_point=crit_point, bubble_rows=bubble_rows)

    print(f"\nLiquid side: {len(isentropes_liquid)} isentropes with a subcooled-liquid branch")
    print(f"Vapor side:  {len(isentropes_vapor)} isentropes with a superheated-vapor branch")

    print("\nIndependently re-solving rho at every emitted point and re-checking entropy...")
    nc_l, nf_l, nsf_l, det_l = verify_branch("liquid", d1, d2, z1, isentropes_liquid, "liquid")
    nc_v, nf_v, nsf_v, det_v = verify_branch("vapor", d1, d2, z1, isentropes_vapor, "vapor")

    print(f"\n=== LIQUID SIDE === checked={nc_l} flagged(>{ENTROPY_MISMATCH_REL_TOL:.0e} rel err)={nf_l} independent-seed-failures={nsf_l}")
    print(f"=== VAPOR SIDE  === checked={nc_v} flagged(>{ENTROPY_MISMATCH_REL_TOL:.0e} rel err)={nf_v} independent-seed-failures={nsf_v}")

    all_flagged = det_l + det_v
    if all_flagged:
        print(f"\n*** {len(all_flagged)} SILENT ENTROPY MISMATCHES FOUND (sol.success said converged, but the point's actual entropy doesn't match its target) ***")
        for d in sorted(all_flagged, key=lambda x: -x["rel_err"])[:20]:
            s_btu = d["s_target_Jmol_K"] / (mw_mix * BTU_LBMR_TO_JKGK)
            print(f"  [{d['branch']}] s_target={s_btu:.2f} Btu/lb-R (pt #{d['point_index']}): "
                  f"T={d['T_K']:.2f}K P={d['P_Pa']/6894.76:.1f}psia "
                  f"s_actual={d['s_actual_Jmol_K']:.4f} vs s_target={d['s_target_Jmol_K']:.4f} J/mol-K "
                  f"rel_err={d['rel_err']:.2e}")
    else:
        print("\n*** NO SILENT ENTROPY MISMATCHES FOUND. Every checked point's independently-recovered "
              "entropy matches its target within tolerance. sol.success was reliable across this real "
              "dataset (unlike R1234yf's Bug #3), though this checks CURRENT parameters only and the "
              "lack of an explicit post-hoc check remains a structural gap for future edits/parameter "
              "changes to rely on. ***")

    if nsf_l + nsf_v > 0:
        print(f"\nNote: {nsf_l + nsf_v} points could not be independently re-solved for rho at all "
              "(damped-Newton continuation lost) -- these are NOT necessarily wrong, just not "
              "independently checkable by this script's method; flagged separately from confirmed "
              "entropy mismatches.")


if __name__ == "__main__":
    sys.exit(main())
