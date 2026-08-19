"""
Validation script for `ancillary_initial_guess.py`'s composite pure-fluid
ancillary saturation-density guess.

Purpose
-------
`mixture_ancillary_saturation_guess_molm3` is NOT a validated saturation
model -- it is seed-quality machinery intended only to warm-start the new
IDAES package's cold-start `initialize()` routine (see the module docstring
in `ancillary_initial_guess.py` for the full rationale). This script
characterizes empirically how close that seed lands relative to the
oracle's own real, converged 3-equation VLE bubble/dew solve
(`solve_bubble_at_t`/`solve_dew_at_t` in `mixture_fully_validated.py`),
across a representative sub-critical temperature range.

This is READ-ONLY reference/regression use of the oracle (permitted under
MASTER TASK spec rule 9's carve-out) -- it is invoked here ONLY for
comparison, never imported by `ancillary_initial_guess.py` itself, and
never modified.

Usage
-----
    cd R515B_idaes_package
    python3 validate_ancillary_guess.py
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from mixture_fully_validated import (  # noqa: E402
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    solve_bubble_at_t,
    solve_dew_at_t,
    solve_mixture_critical_point,
)
from linear_model_codex import bell2023_Tred_vred, BELL_2023_R1234ZE_R227EA  # noqa: E402

from ancillary_initial_guess import mixture_ancillary_saturation_guess_molm3  # noqa: E402

FLUID1, FLUID2 = "r1234ze", "r227ea"
W1 = 0.911
OUT_JSON = HERE / "ancillary_guess_validation_results.json"


def main():
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1, mw2 = mw_from_json(d1), mw_from_json(d2)
    x1 = w1_to_x1(W1, mw1, mw2)

    rhoc1 = float(d1["basic"]["rhoc"]) / mw1
    rhoc2 = float(d2["basic"]["rhoc"]) / mw2

    # Mixture critical temperature, for setting a safe upper bound on the
    # test-temperature grid (ancillary equations, pure or composite, are
    # only meaningful well below Tc).
    tc1_g, tc2_g = float(d1["basic"]["Tc"]), float(d2["basic"]["Tc"])
    vc1_g, vc2_g = 1.0 / rhoc1, 1.0 / rhoc2
    tred_mix_g, vred_mix_g = bell2023_Tred_vred(
        x1, 1.0 - x1, tc1_g, tc2_g, vc1_g, vc2_g, BELL_2023_R1234ZE_R227EA
    )
    crit = solve_mixture_critical_point(d1, d2, x1, t_guess=tred_mix_g, rho_guess=1.0 / vred_mix_g)
    Tc_mix = crit["T_K"]
    print(f"Mixture critical point (oracle solve): Tc={Tc_mix:.4f} K, "
          f"Pc={crit['P_Pa']/1e6:.4f} MPa, rhoc={crit['rho_molm3']:.4f} mol/m^3")

    # Representative sub-critical grid, comfortably clear of the critical
    # region on the high end (ancillary equations get progressively worse
    # approaching Tc by construction -- this is expected and documented,
    # not a bug).
    t_lo = 250.0
    t_hi = Tc_mix - 15.0
    t_grid = np.linspace(t_lo, t_hi, 12)

    rows = []
    rho_l_seed_prev = None
    rho_v_seed_prev = None
    for T in t_grid:
        # Ancillary composite guess (fast, non-iterative).
        rho_l_guess, rho_v_guess = mixture_ancillary_saturation_guess_molm3(d1, d2, x1, float(T))

        # Oracle's own real converged VLE solve, seeded either by
        # continuation from the previous point or (first point only) the
        # oracle's own 80%/1% pure-critical-density heuristic -- exactly
        # the same seeding recipe used in establish_reference_tolerances.py,
        # so this is a fair "real accepted state" reference, not a
        # hand-picked favorable point.
        if rho_l_seed_prev is None:
            rho_l_seed0 = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
            rho_v_seed0 = 0.01 * rho_l_seed0
        else:
            rho_l_seed0, rho_v_seed0 = rho_l_seed_prev, rho_v_seed_prev

        bubble = solve_bubble_at_t(d1, d2, float(T), x1, rho_l_seed0, rho_v_seed0, x1)
        dew = solve_dew_at_t(d1, d2, float(T), x1, rho_l_seed0, rho_v_seed0, x1)
        rho_l_ref = bubble["rho_l_molm3"]
        rho_v_ref = dew["rho_v_molm3"]
        rho_l_seed_prev, rho_v_seed_prev = rho_l_ref, rho_v_ref

        rel_err_l = abs(rho_l_guess - rho_l_ref) / rho_l_ref
        rel_err_v = abs(rho_v_guess - rho_v_ref) / rho_v_ref

        rows.append({
            "T_K": float(T),
            "rho_l_guess_molm3": float(rho_l_guess),
            "rho_l_ref_molm3": float(rho_l_ref),
            "rel_err_l": float(rel_err_l),
            "rho_v_guess_molm3": float(rho_v_guess),
            "rho_v_ref_molm3": float(rho_v_ref),
            "rel_err_v": float(rel_err_v),
            "bubble_converged": bool(bubble.get("converged", bubble.get("success", True))),
            "dew_converged": bool(dew.get("converged", dew.get("success", True))),
        })
        print(f"T={T:7.2f} K  rho_l guess/ref = {rho_l_guess:9.2f}/{rho_l_ref:9.2f} "
              f"(rel_err={rel_err_l:.3%})   rho_v guess/ref = {rho_v_guess:8.3f}/{rho_v_ref:8.3f} "
              f"(rel_err={rel_err_v:.3%})")

    rel_errs_l = np.array([r["rel_err_l"] for r in rows])
    rel_errs_v = np.array([r["rel_err_v"] for r in rows])
    summary = {
        "x1": x1, "w1": W1, "Tc_mix_K": Tc_mix,
        "n_states": len(rows),
        "rel_err_l_max": float(rel_errs_l.max()),
        "rel_err_l_mean": float(rel_errs_l.mean()),
        "rel_err_v_max": float(rel_errs_v.max()),
        "rel_err_v_mean": float(rel_errs_v.mean()),
    }
    print("\n=== Ancillary composite guess vs. oracle real VLE solve: summary ===")
    print(json.dumps(summary, indent=2))

    with open(OUT_JSON, "w") as f:
        json.dump({"summary": summary, "rows": rows}, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")


if __name__ == "__main__":
    sys.exit(main())
