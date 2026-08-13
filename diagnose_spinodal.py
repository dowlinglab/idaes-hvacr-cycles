#!/usr/bin/env python3
"""
Author: Shilpa Narasimhan
Technical support: Claude AI
QA/Testing Responsibility: Shilpa Narasimhan
Creation date: 2026-08-13
Purpose of file: DIAGNOSTIC ONLY (not part of the production CLI). Prints
the raw (dP/drho)_T values across a (T, rho) grid spanning the search
window used by solve_mixture_critical_point in
mixture_dome_validation_pseudo_pure.py, to visualize the actual shape of
the criticality surface. Two consecutive fix attempts to that solver
(window re-centering, then a 4-5x resolution increase) both failed with
"no outer sign-change bracket found" -- this script exists to determine
whether (dP/drho)_T ever actually goes negative anywhere in the intended
search region at all (a real unstable "hump," just poorly resolved) or
never does (a deeper issue: either the EOS surface itself doesn't show
the expected van-der-Waals-loop shape here, consistent with Bell & Jager
(2017)'s explicit warning that "multi-fluid model" L1=0 contours "can be
not smooth" for some mixture combinations, or finite-difference noise
from the same Z-cancellation sensitivity documented earlier in this
project's PROJECT_CONTEXT.md is corrupting the first derivative itself).
Dependencies: numpy, mixture_dome_validation_pseudo_pure, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Usage: python diagnose_spinodal.py
"""

import numpy as np

from linear_model_codex import (
    BELL_2023_R1234ZE_R227EA,
    bell2023_Tred_vred,
    load_idaes_helmholtz_json,
    mw_from_json,
)
from mixture_dome_validation_pseudo_pure import (
    _dp_drho_fd,
    _pressure_rho_derivatives_fd,
    w1_to_x1,
)

W1 = 0.911
FLUID1 = "r1234ze"
FLUID2 = "r227ea"

# Same window solve_mixture_critical_point derived from the last real
# converged bubble/dew run (rho_v_last=2864.2, rho_l_last=5345.8 at
# T=380K, 0.85x/1.15x margins) -- reproduced here as literals so this
# script stays self-contained and doesn't need a fresh sweep to run.
RHO_LO = 2435.5
RHO_HI = 6147.7
T_VALUES = [376.0, 378.0, 379.0, 380.0, 380.5, 381.0, 381.5, 382.0, 382.5, 383.0, 384.0]
N_RHO = 60


def main() -> None:
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = w1_to_x1(W1, mw1, mw2)

    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred_mix, vred_mix = bell2023_Tred_vred(z1, 1.0 - z1, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    print(f"z1={z1:.10f}  Tred_mix={tred_mix:.4f} K  rho_red_mix={1.0 / vred_mix:.2f} mol/m3")
    print(f"Scanning rho in [{RHO_LO}, {RHO_HI}] mol/m3, {N_RHO} points, at each of {len(T_VALUES)} T values.\n")

    rho_values = np.linspace(RHO_LO, RHO_HI, N_RHO)

    for t_k in T_VALUES:
        vals = []
        for rho in rho_values:
            try:
                dpdrho = _dp_drho_fd(d1, d2, t_k, float(rho), z1)
            except Exception:
                dpdrho = float("nan")
            vals.append(dpdrho)
        vals = np.array(vals, dtype=float)
        signs = "".join("+" if v > 0 else ("-" if v < 0 else ("0" if v == 0 else "X")) for v in vals)
        finite = vals[np.isfinite(vals)]
        if len(finite):
            min_v = float(finite.min())
            min_idx = int(np.nanargmin(vals))
            rho_at_min = float(rho_values[min_idx])
        else:
            min_v = float("nan")
            rho_at_min = float("nan")
        n_finite = int(len(finite))
        print(
            f"T={t_k:7.2f}K  min(dP/drho)={min_v: .6e} at rho={rho_at_min:8.2f}  "
            f"finite_pts={n_finite}/{N_RHO}  signs={signs}"
        )

        # If the minimum is negative (a genuine dip found), also report
        # d2P/drho2 there as a sanity cross-check against what Stage 2
        # of the real solver would see.
        if np.isfinite(min_v) and min_v < 0.0:
            try:
                _dp, d2p = _pressure_rho_derivatives_fd(d1, d2, t_k, rho_at_min, z1)
                print(f"           -> d2P/drho2 at that point = {d2p: .6e}")
            except Exception as e:
                print(f"           -> d2P/drho2 evaluation failed: {e}")

    print("\nDone. Look for any line containing a '-' in the sign pattern -- that means a real dip WAS found at that T.")
    print("If EVERY line is all '+' (or 'X'), the EOS never predicts a negative dP/drho anywhere in this window/T-range.")


if __name__ == "__main__":
    main()
