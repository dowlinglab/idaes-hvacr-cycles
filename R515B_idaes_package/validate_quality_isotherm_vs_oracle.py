"""
Stage K (part 1) validation: compute_quality_lines / compute_isotherms_*
in the new independent module `r515b_helmholtz_core.py` vs. the oracle.

Uses a small shared bubble/dew sweep (computed once via the oracle, since
both compute_quality_lines and compute_isotherms_two_phase take
already-solved bubble/dew rows as input -- not a new EOS call) so both
sides start from bit-identical input rows.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from mixture_fully_validated import (  # noqa: E402
    load_idaes_helmholtz_json as oracle_load_json,
    solve_bubble_at_t as oracle_bubble,
    solve_dew_at_t as oracle_dew,
    compute_quality_lines as oracle_quality_lines,
    compute_isotherms_two_phase as oracle_iso_2p,
    compute_isotherms_liquid_side as oracle_iso_liq,
    compute_isotherms_vapor_side as oracle_iso_vap,
    compute_isotherms_supercritical as oracle_iso_sc,
)
import r515b_helmholtz_core as core  # noqa: E402

TOL_REL = 1e-8


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def main():
    d1 = oracle_load_json(core.FLUID1)
    d2 = oracle_load_json(core.FLUID2)
    x1 = core.r515b_x1()

    mw1, mw2 = core.mw_from_json(d1), core.mw_from_json(d2)
    rhoc1, rhoc2 = float(d1["basic"]["rhoc"]) / mw1, float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed0 = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0

    # Shared bubble/dew sweep (oracle only -- these ARE the "already-solved
    # rows" input both compute_quality_lines and compute_isotherms_two_phase
    # consume; both sides get the SAME rows, isolating the test to whether
    # the new module's post-processing/root-finding matches, not whether
    # its VLE solve matches -- that's already validated in Stage I).
    t_grid = np.linspace(260.0, 350.0, 10)
    bubble_rows, dew_rows = [], []
    rl, rv = rho_l_seed0, rho_v_seed0
    for T in t_grid:
        b = oracle_bubble(d1, d2, float(T), x1, rl, rv, x1)
        d = oracle_dew(d1, d2, float(T), x1, rl, rv, x1)
        bubble_rows.append(b)
        dew_rows.append(d)
        rl, rv = b["rho_l_molm3"], b["rho_v_molm3"]

    all_pass = True

    # --- Quality lines ---
    ref_q = oracle_quality_lines(bubble_rows, dew_rows)
    new_q = core.compute_quality_lines(bubble_rows, dew_rows)
    q_pass = True
    for q in core.QUALITY_LINE_VALUES:
        ref_rows, new_rows = ref_q[q], new_q[q]
        if len(ref_rows) != len(new_rows):
            q_pass = False
            continue
        for rr, nr in zip(ref_rows, new_rows):
            if rel(nr["h_Jmol"], rr["h_Jmol"]) > TOL_REL or rel(nr["P_Pa"], rr["P_Pa"]) > TOL_REL:
                q_pass = False
    all_pass = all_pass and q_pass
    print(f"Quality lines: {'PASS' if q_pass else 'FAIL'} ({sum(len(v) for v in ref_q.values())} total rows checked)")

    # --- Two-phase isotherms ---
    ref_iso2p = oracle_iso_2p(d1, d2, x1, bubble_rows, dew_rows)
    new_iso2p = core.compute_isotherms_two_phase(d1, d2, x1, bubble_rows, dew_rows)
    iso2p_pass = set(ref_iso2p.keys()) == set(new_iso2p.keys())
    for tf in ref_iso2p:
        r, n = ref_iso2p[tf], new_iso2p[tf]
        if (rel(n["P_Pa"], r["P_Pa"]) > TOL_REL or rel(n["h_l_Jmol"], r["h_l_Jmol"]) > TOL_REL or
                rel(n["h_v_Jmol"], r["h_v_Jmol"]) > TOL_REL):
            iso2p_pass = False
    all_pass = all_pass and iso2p_pass
    print(f"Two-phase isotherms: {'PASS' if iso2p_pass else 'FAIL'} ({len(ref_iso2p)} T-values)")

    # --- Liquid-side / vapor-side / supercritical isotherm extensions ---
    ref_liq = oracle_iso_liq(d1, d2, x1, ref_iso2p)
    new_liq = core.compute_isotherms_liquid_side(d1, d2, x1, new_iso2p)
    liq_pass = set(ref_liq.keys()) == set(new_liq.keys())
    max_liq_err = 0.0
    for tf in ref_liq:
        for rp, npnt in zip(ref_liq[tf], new_liq[tf]):
            e = max(rel(npnt["P_Pa"], rp["P_Pa"]), rel(npnt["h_Jmol"], rp["h_Jmol"]))
            max_liq_err = max(max_liq_err, e)
            if e > TOL_REL:
                liq_pass = False
    all_pass = all_pass and liq_pass
    print(f"Liquid-side isotherm ext: {'PASS' if liq_pass else 'FAIL'} ({len(ref_liq)} T-values, max_err={max_liq_err:.3g})")

    ref_vap = oracle_iso_vap(d1, d2, x1, ref_iso2p)
    new_vap = core.compute_isotherms_vapor_side(d1, d2, x1, new_iso2p)
    vap_pass = set(ref_vap.keys()) == set(new_vap.keys())
    max_vap_err = 0.0
    for tf in ref_vap:
        for rp, npnt in zip(ref_vap[tf], new_vap[tf]):
            e = max(rel(npnt["P_Pa"], rp["P_Pa"]), rel(npnt["h_Jmol"], rp["h_Jmol"]))
            max_vap_err = max(max_vap_err, e)
            if e > TOL_REL:
                vap_pass = False
    all_pass = all_pass and vap_pass
    print(f"Vapor-side isotherm ext: {'PASS' if vap_pass else 'FAIL'} ({len(ref_vap)} T-values, max_err={max_vap_err:.3g})")

    ref_sc = oracle_iso_sc(d1, d2, x1)
    new_sc = core.compute_isotherms_supercritical(d1, d2, x1)
    sc_pass = set(ref_sc.keys()) == set(new_sc.keys())
    max_sc_err = 0.0
    for tf in ref_sc:
        for rp, npnt in zip(ref_sc[tf], new_sc[tf]):
            e = max(rel(npnt["P_Pa"], rp["P_Pa"]), rel(npnt["h_Jmol"], rp["h_Jmol"]))
            max_sc_err = max(max_sc_err, e)
            if e > TOL_REL:
                sc_pass = False
    all_pass = all_pass and sc_pass
    print(f"Supercritical isotherms: {'PASS' if sc_pass else 'FAIL'} ({len(ref_sc)} T-values, max_err={max_sc_err:.3g})")

    print(f"\nOVERALL Stage K (part 1, quality+isotherms): {'PASS' if all_pass else 'FAIL'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
