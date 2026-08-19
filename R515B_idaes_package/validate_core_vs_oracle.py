"""
Stage F/G/H re-validation of the NEW INDEPENDENT production module
`r515b_helmholtz_core.py` (not `linear_model_codex.compute_table1_properties`)
against the oracle, at bit-identical input states (computed once, reused for
both sides, to avoid any input-precision mismatch masquerading as a model
discrepancy -- P is extremely sensitive to rho on the deep-liquid branch,
see mix_state's own docstring note on Z~0.2 cancellation).

Read-only reference use of the oracle for comparison only (spec rule 9).
`r515b_helmholtz_core.py` itself never imports the oracle.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
ORACLE_DIR = HERE.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))

from mixture_fully_validated import (  # noqa: E402
    load_idaes_helmholtz_json as oracle_load_json,
    mix_state as oracle_mix_state,
    _mix_entropy_direct as oracle_entropy,
    chemical_potentials_analytic as oracle_mu,
    solve_bubble_at_t as oracle_bubble,
    solve_dew_at_t as oracle_dew,
    solve_mixture_critical_point as oracle_crit,
)
from linear_model_codex import bell2023_Tred_vred, BELL_2023_R1234ZE_R227EA  # noqa: E402
import r515b_helmholtz_core as core  # noqa: E402

TOL_REL = 1e-8


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def main():
    d1 = oracle_load_json(core.FLUID1)
    d2 = oracle_load_json(core.FLUID2)
    x1 = core.r515b_x1()
    print(f"x1={x1:.12f}")

    T_SAT = 300.0
    rhoc1 = float(d1["basic"]["rhoc"]) / core.mw_from_json(d1)
    rhoc2 = float(d2["basic"]["rhoc"]) / core.mw_from_json(d2)
    rho_l_seed0 = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0

    bubble_ref = oracle_bubble(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1)
    dew_ref = oracle_dew(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1)
    RHO_LIQ = 1.10 * bubble_ref["rho_l_molm3"]
    RHO_VAP = 0.50 * dew_ref["rho_v_molm3"]

    all_pass = True
    print("\n--- Direct-state re-check: r515b_helmholtz_core.py vs oracle (bit-identical inputs) ---")
    for label, T, rho in [
        ("liquid_direct_subcooled", T_SAT, RHO_LIQ),
        ("vapor_direct_superheated", T_SAT, RHO_VAP),
        ("supercritical_direct", 389.5, 3000.0),
    ]:
        st_ref = oracle_mix_state(d1, d2, T, rho, x1)
        st_new = core.mix_state(d1, d2, T, rho, x1)
        s_ref = oracle_entropy(d1, d2, T, rho, x1)
        s_new = core.mix_entropy_direct(d1, d2, T, rho, x1)
        mu1_ref, mu2_ref = oracle_mu(d1, d2, T, rho, x1)
        mu1_new, mu2_new = core.chemical_potentials_analytic(d1, d2, T, rho, x1)

        p_rel = rel(st_new.p_pa, st_ref.p_pa)
        h_rel = rel(st_new.h_jmol, st_ref.h_jmol)
        s_rel = rel(s_new, s_ref)
        mu1_rel = rel(mu1_new, mu1_ref)
        mu2_rel = rel(mu2_new, mu2_ref)
        ok = max(p_rel, h_rel, s_rel, mu1_rel, mu2_rel) <= TOL_REL
        all_pass = all_pass and ok
        print(f"[{label}] T={T:.3f} rho={rho:.6f}  status={'PASS' if ok else 'FAIL'}")
        print(f"    P rel_err={p_rel:.3g}  h rel_err={h_rel:.3g}  s(rebased) rel_err={s_rel:.3g}  "
              f"mu1 rel_err={mu1_rel:.3g}  mu2 rel_err={mu2_rel:.3g}")

    print(f"\nDirect-state OVERALL: {'PASS' if all_pass else 'FAIL'} (tolerance {TOL_REL:.0e})")

    # Stage I: bubble/dew VLE solve comparison across a small T sweep.
    print("\n--- Stage I: solve_bubble_at_t / solve_dew_at_t, new module vs oracle ---")
    TOL_REL_VLE = 1e-4  # tier-2 (iterative-equilibrium) frozen tolerance, Section 1
    vle_pass = True
    rho_l_prev, rho_v_prev = rho_l_seed0, rho_v_seed0
    for T in [270.0, 285.0, 300.0, 315.0, 330.0, 345.0]:
        b_ref = oracle_bubble(d1, d2, T, x1, rho_l_prev, rho_v_prev, x1)
        b_new = core.solve_bubble_at_t(d1, d2, T, x1, rho_l_prev, rho_v_prev, x1)
        d_ref = oracle_dew(d1, d2, T, x1, rho_l_prev, rho_v_prev, x1)
        d_new = core.solve_dew_at_t(d1, d2, T, x1, rho_l_prev, rho_v_prev, x1)

        b_ok = (b_ref["status"] == "CONVERGED" and b_new["status"] == "CONVERGED" and
                rel(b_new["P_Pa"], b_ref["P_Pa"]) <= TOL_REL_VLE and
                rel(b_new["rho_l_molm3"], b_ref["rho_l_molm3"]) <= TOL_REL_VLE and
                rel(b_new["rho_v_molm3"], b_ref["rho_v_molm3"]) <= TOL_REL_VLE and
                rel(b_new["y1_vap"], b_ref["y1_vap"]) <= TOL_REL_VLE)
        d_ok = (d_ref["status"] == "CONVERGED" and d_new["status"] == "CONVERGED" and
                rel(d_new["P_Pa"], d_ref["P_Pa"]) <= TOL_REL_VLE and
                rel(d_new["rho_l_molm3"], d_ref["rho_l_molm3"]) <= TOL_REL_VLE and
                rel(d_new["rho_v_molm3"], d_ref["rho_v_molm3"]) <= TOL_REL_VLE and
                rel(d_new["x1_liq"], d_ref["x1_liq"]) <= TOL_REL_VLE)
        vle_pass = vle_pass and b_ok and d_ok
        print(f"T={T:.1f}K  bubble: ref_status={b_ref['status']} new_status={b_new['status']} "
              f"P_rel={rel(b_new['P_Pa'], b_ref['P_Pa']):.3g} rho_l_rel={rel(b_new['rho_l_molm3'], b_ref['rho_l_molm3']):.3g} "
              f"rho_v_rel={rel(b_new['rho_v_molm3'], b_ref['rho_v_molm3']):.3g} y1_rel={rel(b_new['y1_vap'], b_ref['y1_vap']):.3g} "
              f"[{'PASS' if b_ok else 'FAIL'}]")
        print(f"           dew:    ref_status={d_ref['status']} new_status={d_new['status']} "
              f"P_rel={rel(d_new['P_Pa'], d_ref['P_Pa']):.3g} rho_l_rel={rel(d_new['rho_l_molm3'], d_ref['rho_l_molm3']):.3g} "
              f"rho_v_rel={rel(d_new['rho_v_molm3'], d_ref['rho_v_molm3']):.3g} x1_rel={rel(d_new['x1_liq'], d_ref['x1_liq']):.3g} "
              f"[{'PASS' if d_ok else 'FAIL'}]")

        # continuation seeding for next T, from the ORACLE's converged state
        # (matching run_true_vle_envelope's own continuation pattern)
        rho_l_prev, rho_v_prev = b_ref["rho_l_molm3"], b_ref["rho_v_molm3"]

    print(f"\nStage I OVERALL: {'PASS' if vle_pass else 'FAIL'} (tolerance {TOL_REL_VLE:.0e} relative, tier-2)")

    # Stage J: mixture critical-point solve comparison.
    print("\n--- Stage J: solve_mixture_critical_point, new module vs oracle ---")
    crit_pass = True

    # Case A: crude Bell-reducing-point guess (no real sweep data available) --
    # documented finding: this guess does NOT actually converge (returns the
    # fallback dict, T_K/rho_molm3 = the guess itself, converged=False) for
    # BOTH the oracle and the new module. This is NOT a bug -- it's a
    # pre-existing characteristic of this crude-guess call path, faithfully
    # reproduced. IMPORTANT CORRECTION vs. earlier Stage C usage: values
    # like "Tc_mix=381.5224K" quoted earlier this session (establish_
    # reference_tolerances.py, validate_ancillary_guess.py) came from this
    # NON-CONVERGED fallback path -- they are the Bell(2023) mixing-rule
    # reducing-point APPROXIMATION (Tred_mix), not a genuinely solved
    # critical point. Harmless where used (only as a safe upper bound for a
    # temperature grid), but must not be mistaken for a validated Tc.
    tc1_g, tc2_g = float(d1["basic"]["Tc"]), float(d2["basic"]["Tc"])
    vc1_g, vc2_g = 1.0 / rhoc1, 1.0 / rhoc2
    tred_g, vred_g = bell2023_Tred_vred(x1, 1.0 - x1, tc1_g, tc2_g, vc1_g, vc2_g, BELL_2023_R1234ZE_R227EA)
    ref_a = oracle_crit(d1, d2, x1, t_guess=tred_g, rho_guess=1.0 / vred_g)
    new_a = core.solve_mixture_critical_point(d1, d2, x1, t_guess=tred_g, rho_guess=1.0 / vred_g)
    case_a_ok = (ref_a["converged"] == new_a["converged"] == False and  # noqa: E712
                 rel(new_a["T_K"], ref_a["T_K"]) <= TOL_REL and rel(new_a["rho_molm3"], ref_a["rho_molm3"]) <= TOL_REL)
    crit_pass = crit_pass and case_a_ok
    print(f"[crude-guess, expected non-convergence] ref: T={ref_a['T_K']:.4f}K converged={ref_a['converged']}  "
          f"new: T={new_a['T_K']:.4f}K converged={new_a['converged']}  [{'PASS' if case_a_ok else 'FAIL'}]")
    print(f"    (Tred_mix approximation reported here: {tred_g:.4f} K -- NOT the true mixture Tc)")

    # Case B: real-data-informed guess/window (mirrors run_true_vle_envelope's
    # actual usage pattern when sweep data exists) -- genuinely converges.
    rho_l_seed0_ = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
    rho_v_seed0_ = 0.01 * rho_l_seed0_
    b350 = oracle_bubble(d1, d2, 350.0, x1, rho_l_seed0_, rho_v_seed0_, x1)
    rho_guess_b = 0.5 * (b350["rho_l_molm3"] + b350["rho_v_molm3"])
    ref_b = oracle_crit(d1, d2, x1, t_guess=372.0, rho_guess=rho_guess_b, t_scan_lo=350.0, t_scan_hi=385.0)
    new_b = core.solve_mixture_critical_point(d1, d2, x1, t_guess=372.0, rho_guess=rho_guess_b, t_scan_lo=350.0, t_scan_hi=385.0)
    case_b_ok = (ref_b["converged"] and new_b["converged"] and
                 rel(new_b["T_K"], ref_b["T_K"]) <= TOL_REL and
                 rel(new_b["rho_molm3"], ref_b["rho_molm3"]) <= TOL_REL and
                 rel(new_b["P_Pa"], ref_b["P_Pa"]) <= TOL_REL and
                 rel(new_b["h_Jmol"], ref_b["h_Jmol"]) <= TOL_REL)
    crit_pass = crit_pass and case_b_ok
    print(f"[real-data-informed guess, genuine convergence] ref: T={ref_b['T_K']:.6f}K rho={ref_b['rho_molm3']:.6f} "
          f"P={ref_b['P_Pa']:.4f}Pa h={ref_b['h_Jmol']:.4f}J/mol converged={ref_b['converged']}")
    print(f"                                                 new: T={new_b['T_K']:.6f}K rho={new_b['rho_molm3']:.6f} "
          f"P={new_b['P_Pa']:.4f}Pa h={new_b['h_Jmol']:.4f}J/mol converged={new_b['converged']}  [{'PASS' if case_b_ok else 'FAIL'}]")
    print("    THIS is the true validated mixture critical point: "
          f"Tc={ref_b['T_K']:.4f}K, rhoc={ref_b['rho_molm3']:.4f} mol/m^3, Pc={ref_b['P_Pa']/1e6:.6f} MPa")

    print(f"\nStage J OVERALL: {'PASS' if crit_pass else 'FAIL'} (tolerance {TOL_REL:.0e} relative)")

    return 0 if (all_pass and vle_pass and crit_pass) else 1


if __name__ == "__main__":
    sys.exit(main())
