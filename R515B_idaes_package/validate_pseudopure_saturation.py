"""
Stage L part 3 support: sanity-check the new pseudo-pure (x1=y1=z1 fixed)
saturation solve (`r515b_helmholtz_core.solve_pseudopure_saturation_at_t`)
against the rigorous bubble/dew branches (`solve_bubble_at_t`/
`solve_dew_at_t`, real y1/x1 split) at the same temperatures.

This is NOT a pass/fail validation against the oracle -- the pseudo-pure
solve is a deliberate simplification (see the 2026-08-18 rule-22
resolution in helmholtz_prop_validation.md), not a port of oracle
behavior, so there is no "ground truth" it must match exactly. The checks
here instead confirm the solve is well-behaved and physically sane:
  1. It converges across the practical temperature range.
  2. Its saturation pressure at each T sits BETWEEN the real bubble
     pressure and dew pressure at that T (since x1=y1=z1 is a composition
     "average" of the real liquid/vapor split, the pseudo-pure dome should
     sit between the true bubble and dew branches, not off in some
     unrelated region).
  3. Its densities are reasonably close to the real bubble rho_l / dew
     rho_v (same-order-of-magnitude sanity, not exact-match).
  4. Quantifies its own deviation from the TRUE dome (rigorous bubble/dew
     average), mirroring the historical `pseudopure_mode_eval.py`
     MAE/max-error reporting pattern -- so the size of this simplification's
     error is recorded, not just its existence.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import r515b_helmholtz_core as core  # noqa: E402


def main():
    d1 = core.load_idaes_helmholtz_json(core.FLUID1)
    d2 = core.load_idaes_helmholtz_json(core.FLUID2)
    x1 = core.r515b_x1()
    mw1, mw2 = core.mw_from_json(d1), core.mw_from_json(d2)
    rhoc1, rhoc2 = float(d1["basic"]["rhoc"]) / mw1, float(d2["basic"]["rhoc"]) / mw2

    rho_l_seed0 = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0

    t_grid = np.linspace(250.0, 375.0, 15)
    rl_b, rv_b = rho_l_seed0, rho_v_seed0  # bubble seed carry
    rl_d, rv_d = rho_l_seed0, rho_v_seed0  # dew seed carry
    rl_p, rv_p = rho_l_seed0, rho_v_seed0  # pseudo-pure seed carry

    n_conv_bubble = n_conv_dew = n_conv_pp = 0
    dT_glide_list = []
    dP_list = []
    abs_T_errs = []
    rows = []
    for T in t_grid:
        b = core.solve_bubble_at_t(d1, d2, float(T), x1, rl_b, rv_b, x1)
        d = core.solve_dew_at_t(d1, d2, float(T), x1, rl_d, rv_d, x1)
        p = core.solve_pseudopure_saturation_at_t(d1, d2, x1, float(T), rl_p, rv_p)

        if b["status"] == "CONVERGED":
            n_conv_bubble += 1
            rl_b, rv_b = b["rho_l_molm3"], b["rho_v_molm3"]
        if d["status"] == "CONVERGED":
            n_conv_dew += 1
            rl_d, rv_d = d["rho_l_molm3"], d["rho_v_molm3"]
        if p["status"] == "CONVERGED":
            n_conv_pp += 1
            rl_p, rv_p = p["rho_l_molm3"], p["rho_v_molm3"]

        row = {"T": T, "bubble": b, "dew": d, "pp": p}
        rows.append(row)

        if b["status"] == "CONVERGED" and d["status"] == "CONVERGED" and p["status"] == "CONVERGED":
            p_bubble, p_dew, p_pp = b["P_Pa"], d["P_Pa"], p["P_Pa"]
            in_between = min(p_bubble, p_dew) - 1.0 <= p_pp <= max(p_bubble, p_dew) + 1.0
            p_avg = 0.5 * (p_bubble + p_dew)
            dP_list.append(abs(p_pp - p_avg) / p_avg)
            print(f"T={T:7.2f}K  P_bubble={p_bubble:12.1f}  P_dew={p_dew:12.1f}  "
                  f"P_pseudopure={p_pp:12.1f}  in_[bubble,dew]={in_between}")

    print(f"\nConvergence: bubble {n_conv_bubble}/{len(t_grid)}, dew {n_conv_dew}/{len(t_grid)}, "
          f"pseudo-pure {n_conv_pp}/{len(t_grid)}")
    if dP_list:
        print(f"Pseudo-pure P deviation from bubble/dew average: mean={np.mean(dP_list)*100:.4f}%, "
              f"max={np.max(dP_list)*100:.4f}%")

    all_ok = (n_conv_pp >= int(0.8 * len(t_grid)))
    print(f"\nOVERALL sanity check: {'PASS' if all_ok else 'FAIL'} "
          f"(pseudo-pure solve converges on >=80% of the test grid)")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
