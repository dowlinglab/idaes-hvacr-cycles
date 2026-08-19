"""
Stage K (part 2) validation: isentrope machinery in the new independent
module `r515b_helmholtz_core.py` vs. `mixture_isentrope_validation.py`
(the insurance-patched sibling of the base oracle -- see this module's
own docstring for why THIS file, not the base `mixture_fully_validated.py`,
is the reference for isentropes specifically).

Uses a shared bubble/dew sweep + critical point (computed once via the
reference file) so both sides start from bit-identical input rows.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from mixture_isentrope_validation import (  # noqa: E402
    load_idaes_helmholtz_json as ref_load_json,
    solve_bubble_at_t as ref_bubble,
    solve_dew_at_t as ref_dew,
    solve_mixture_critical_point as ref_crit,
    compute_isentropes_two_phase as ref_isen_2p,
    compute_isentrope_liquid_side as ref_isen_liq,
    compute_isentrope_vapor_side as ref_isen_vap,
    ISENTROPE_VALUES_BTU_LBMR as REF_ISEN_BTU,
    BTU_LBMR_TO_JKGK as REF_BTU_CONV,
)
from linear_model_codex import bell2023_Tred_vred, BELL_2023_R1234ZE_R227EA  # noqa: E402
import r515b_helmholtz_core as core  # noqa: E402

TOL_REL = 1e-8


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def main():
    d1 = ref_load_json(core.FLUID1)
    d2 = ref_load_json(core.FLUID2)
    x1 = core.r515b_x1()
    mw1, mw2 = core.mw_from_json(d1), core.mw_from_json(d2)
    mw_mix = x1 * mw1 + (1.0 - x1) * mw2

    rhoc1, rhoc2 = float(d1["basic"]["rhoc"]) / mw1, float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed0 = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0

    t_grid = np.linspace(255.0, 349.0, 20)
    bubble_rows, dew_rows = [], []
    rl, rv = rho_l_seed0, rho_v_seed0
    for T in t_grid:
        b = ref_bubble(d1, d2, float(T), x1, rl, rv, x1)
        d = ref_dew(d1, d2, float(T), x1, rl, rv, x1)
        bubble_rows.append(b)
        dew_rows.append(d)
        if b["status"] == "CONVERGED":
            rl, rv = b["rho_l_molm3"], b["rho_v_molm3"]

    # Critical point, real-data-informed (Stage J's own validated pattern:
    # seed from a real bubble solve at T=350K, explicit t_scan_lo=350K).
    b350 = ref_bubble(d1, d2, 350.0, x1, rho_l_seed0, rho_v_seed0, x1)
    rho_guess_c = 0.5 * (b350["rho_l_molm3"] + b350["rho_v_molm3"])
    crit = ref_crit(d1, d2, x1, t_guess=372.0, rho_guess=rho_guess_c, t_scan_lo=350.0, t_scan_hi=385.0)
    print(f"crit_point converged={crit['converged']} T={crit['T_K']:.4f}K")

    # Exact conversion the oracle's own _cli() uses (mixture_isentrope_
    # validation.py line 3962): Btu/(lbm-R) -> J/(kg-K) -> J/(mol-K) via
    # mw_mix [kg/mol]. An earlier draft of this script divided by an extra
    # 1000 (assuming mw_mix was in g/mol) -- that was WRONG and made every
    # s_target ~1000x too small, hence 0 rows/points everywhere. Fixed here.
    s_values = [bv * REF_BTU_CONV * mw_mix for bv in REF_ISEN_BTU]

    all_pass = True

    # --- Two-phase isentropes ---
    ref_2p = ref_isen_2p(d1, d2, bubble_rows, dew_rows, s_values)
    new_2p = core.compute_isentropes_two_phase(d1, d2, bubble_rows, dew_rows, s_values)
    p2p_pass = True
    n_rows_2p = 0
    for s in s_values:
        rr, nr = ref_2p[s], new_2p[s]
        if len(rr) != len(nr):
            p2p_pass = False
            continue
        for a, b in zip(rr, nr):
            n_rows_2p += 1
            if rel(b["h_Jmol"], a["h_Jmol"]) > TOL_REL or rel(b["P_Pa"], a["P_Pa"]) > TOL_REL:
                p2p_pass = False
    all_pass = all_pass and p2p_pass
    print(f"Two-phase isentropes: {'PASS' if p2p_pass else 'FAIL'} ({n_rows_2p} rows)")

    # --- Liquid-side isentrope extension ---
    ref_liq = ref_isen_liq(d1, d2, x1, bubble_rows, s_values, crit_point=crit)
    new_liq = core.compute_isentrope_liquid_side(d1, d2, x1, bubble_rows, s_values, crit_point=crit)
    liq_pass = set(ref_liq.keys()) == set(new_liq.keys())
    n_liq_pts = 0
    max_liq_err = 0.0
    for s in ref_liq:
        rr, nr = ref_liq[s], new_liq[s]
        if len(rr) != len(nr):
            liq_pass = False
            continue
        for a, b in zip(rr, nr):
            n_liq_pts += 1
            e = max(rel(b["P_Pa"], a["P_Pa"]), rel(b["h_Jmol"], a["h_Jmol"]), rel(b["T_K"], a["T_K"]))
            max_liq_err = max(max_liq_err, e)
            if e > TOL_REL:
                liq_pass = False
    all_pass = all_pass and liq_pass
    print(f"Liquid-side isentrope ext: {'PASS' if liq_pass else 'FAIL'} ({len(ref_liq)} isentropes, "
          f"{n_liq_pts} points, max_err={max_liq_err:.3g})")

    # --- Vapor-side isentrope extension ---
    ref_vap = ref_isen_vap(d1, d2, x1, dew_rows, s_values, crit_point=crit, bubble_rows=bubble_rows)
    new_vap = core.compute_isentrope_vapor_side(d1, d2, x1, dew_rows, s_values, crit_point=crit, bubble_rows=bubble_rows)
    vap_pass = set(ref_vap.keys()) == set(new_vap.keys())
    n_vap_pts = 0
    max_vap_err = 0.0
    for s in ref_vap:
        rr, nr = ref_vap[s], new_vap[s]
        if len(rr) != len(nr):
            vap_pass = False
            print(f"  [len mismatch] s={s} ref_len={len(rr)} new_len={len(nr)}")
            continue
        for a, b in zip(rr, nr):
            n_vap_pts += 1
            e = max(rel(b["P_Pa"], a["P_Pa"]), rel(b["h_Jmol"], a["h_Jmol"]), rel(b["T_K"], a["T_K"]))
            max_vap_err = max(max_vap_err, e)
            if e > TOL_REL:
                vap_pass = False
    all_pass = all_pass and vap_pass
    print(f"Vapor-side isentrope ext: {'PASS' if vap_pass else 'FAIL'} ({len(ref_vap)} isentropes, "
          f"{n_vap_pts} points, max_err={max_vap_err:.3g})")

    print(f"\nOVERALL Stage K (part 2, isentropes): {'PASS' if all_pass else 'FAIL'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
