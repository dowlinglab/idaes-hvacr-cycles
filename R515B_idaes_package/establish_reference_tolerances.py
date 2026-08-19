"""
Stage C (MASTER TASK spec, rule 41/52): establish and freeze regression
tolerances by measuring the oracle's own reproducibility.

Purpose
-------
Before any new R-515B property-package code is judged against
`mixture_fully_validated.py` (the designated oracle -- confirmed 2026-08-17
to be functionally equivalent to the actively-worked `mixture_isentrope_
validation.py`, see helmholtz_prop_validation.md Section 0 and
PROJECT_CONTEXT.md's 2026-08-17 "MASTER TASK kickoff" entry for the full
diff record), this script evaluates the oracle REPEATEDLY at deterministic
representative states and measures how much its own output varies run to
run. Pure floating-point numpy/scipy evaluation with no random seeding is
expected to be exactly bit-reproducible; this script confirms that rather
than assuming it, and uses the measured spread (if any) to set frozen
regression tolerances with a safety margin on top.

This is READ-ONLY use of the oracle for reference/regression purposes only
(permitted under MASTER TASK spec rule 9's explicit carve-out: "Only new
validation/regression code may invoke mixture_fully_validated.py for
reference comparisons"). `mixture_fully_validated.py` itself is never
modified, and the eventual production property package will NOT import
this module or any oracle file at runtime.

Representative states (deterministic, chosen to span the domain the spec
requires -- rule 41/52):
  1. Liquid state       : direct mix_state() at fixed (T, rho_mol) deep in the
                           subcooled-liquid single-phase region.
  2. Vapor state         : direct mix_state() at fixed (T, rho_mol) deep in the
                           superheated-vapor single-phase region.
  3. Saturation state     : solve_bubble_at_t() at a representative sub-critical T.
  4. Two-phase-derived    : quality=0.5 lever-rule state built from the bubble/dew
                           rows at the same T as (3), using the oracle's own
                           compute_quality_lines().
  5. Supercritical state  : direct mix_state() at fixed (T, rho_mol) above the
                           oracle's own solved mixture critical temperature.
  6. Near-critical state  : solve_mixture_critical_point() itself (the most
                           numerically delicate calculation in the whole file).

Usage
-----
    cd R515B_idaes_package
    python3 establish_reference_tolerances.py

Outputs a summary to stdout and a JSON file
(reference_repeatability_results.json) with full numeric detail, both of
which are cited from helmholtz_prop_validation.md Section 1.
"""

import json
import sys
from pathlib import Path

import numpy as np

# The oracle lives in a sibling directory and itself does `from
# linear_model_codex import (...)` assuming linear_model_codex.py is on
# sys.path (it sits right next to the oracle) -- so that sibling directory,
# not this one, must be on sys.path for the oracle's own import to resolve.
ORACLE_DIR = Path(__file__).parent.parent / "R515B_props_validated"
sys.path.insert(0, str(ORACLE_DIR))

from mixture_fully_validated import (  # noqa: E402  (see sys.path note above)
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    mix_state,
    solve_bubble_at_t,
    solve_dew_at_t,
    solve_mixture_critical_point,
    compute_quality_lines,
    run_true_vle_envelope,
)

FLUID1, FLUID2 = "r1234ze", "r227ea"
W1 = 0.911
N_REPEATS = 5  # number of repeated evaluations per state, for repeatability measurement

HERE = Path(__file__).parent
OUT_JSON = HERE / "reference_repeatability_results.json"


def _repeatability(label, fn, extract_fields):
    """Run fn() N_REPEATS times, extract named numeric fields from each
    result via extract_fields (a dict of name -> callable(result)->float),
    and report min/max/mean/spread per field."""
    runs = [fn() for _ in range(N_REPEATS)]
    report = {"label": label, "n_runs": N_REPEATS, "fields": {}}
    for name, getter in extract_fields.items():
        vals = np.array([getter(r) for r in runs], dtype=float)
        spread = float(vals.max() - vals.min())
        rel_spread = spread / max(1.0, abs(float(vals.mean())))
        report["fields"][name] = {
            "mean": float(vals.mean()),
            "min": float(vals.min()),
            "max": float(vals.max()),
            "abs_spread": spread,
            "rel_spread": rel_spread,
        }
    return report


def main():
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1, mw2 = mw_from_json(d1), mw_from_json(d2)
    x1 = w1_to_x1(W1, mw1, mw2)
    print(f"x1 (mole fraction R-1234ze(E))) = {x1:.10f} at w1={W1}")

    results = []

    # 3. Saturation state (bubble/dew solve at representative sub-critical T).
    # Computed FIRST so states 1/2 (subcooled liquid / superheated vapor) can
    # be anchored to these actual solved saturation densities at the SAME T,
    # rather than an arbitrary (T,rho) guess -- an early version of this
    # script picked an unanchored (T=280K, rho=9500 mol/m^3) "liquid" point
    # and got p_pa = -14.25 MPa: not a bug, just an unphysical/metastable
    # state outside the real liquid branch at that T for this EOS (this
    # mixture's liquid density at 280K is much higher than 9500 mol/m^3
    # would suggest is "deep liquid" -- picking density relative to the
    # actual bubble curve avoids ever landing in this kind of invalid
    # region again).
    # solve_bubble_at_t/solve_dew_at_t require explicit initial density/
    # composition guesses (no defaults) -- reuse the oracle's OWN first-point
    # seeding recipe from run_true_vle_envelope (rough 80%/1% pure-critical-
    # density-weighted guess, y1/x1 guess = z1) rather than inventing a
    # different seed, so this repeatability check exercises the exact same
    # numerical path a real sweep uses.
    T_SAT = 300.0
    rhoc1 = float(d1["basic"]["rhoc"]) / mw1
    rhoc2 = float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed0 = 0.8 * (x1 * rhoc1 + (1.0 - x1) * rhoc2)
    rho_v_seed0 = 0.01 * rho_l_seed0
    results.append(_repeatability(
        "bubble_saturation",
        lambda: solve_bubble_at_t(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1),
        {
            "P_Pa": lambda r: r["P_Pa"],
            "rho_l_molm3": lambda r: r["rho_l_molm3"],
            "rho_v_molm3": lambda r: r["rho_v_molm3"],
            "y1_vap": lambda r: r["y1_vap"],
        },
    ))
    results.append(_repeatability(
        "dew_saturation",
        lambda: solve_dew_at_t(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1),
        {
            "P_Pa": lambda r: r["P_Pa"],
            "rho_l_molm3": lambda r: r["rho_l_molm3"],
            "rho_v_molm3": lambda r: r["rho_v_molm3"],
            "x1_liq": lambda r: r["x1_liq"],
        },
    ))

    # 1/2. Subcooled-liquid / superheated-vapor states, anchored to the
    # ACTUAL solved saturation densities at T_SAT above (rather than an
    # unanchored guess -- see the note above state 3): 10% denser than the
    # real saturated-liquid density for the liquid state, half the real
    # saturated-vapor density for the vapor state, both at the same T_SAT so
    # they're genuinely single-phase (liquid side / vapor side respectively)
    # states for this mixture, not accidentally inside or past the dome.
    bubble_ref = solve_bubble_at_t(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1)
    dew_ref = solve_dew_at_t(d1, d2, T_SAT, x1, rho_l_seed0, rho_v_seed0, x1)
    RHO_LIQ = 1.10 * bubble_ref["rho_l_molm3"]
    RHO_VAP = 0.50 * dew_ref["rho_v_molm3"]
    results.append(_repeatability(
        "liquid_direct_subcooled",
        lambda: mix_state(d1, d2, T_SAT, RHO_LIQ, x1),
        {"p_pa": lambda r: r.p_pa, "h_jmol": lambda r: r.h_jmol, "g_jmol": lambda r: r.g_jmol},
    ))
    results.append(_repeatability(
        "vapor_direct_superheated",
        lambda: mix_state(d1, d2, T_SAT, RHO_VAP, x1),
        {"p_pa": lambda r: r.p_pa, "h_jmol": lambda r: r.h_jmol, "g_jmol": lambda r: r.g_jmol},
    ))

    # 4. Two-phase-derived quality state (quality=0.5 lever rule at T_SAT,
    # built from a small bubble/dew sweep via the oracle's own
    # run_true_vle_envelope + compute_quality_lines -- exercises the same
    # code path the oracle itself uses for interior two-phase states).
    def quality_state_fn():
        t_vals = np.array([T_SAT - 1.0, T_SAT, T_SAT + 1.0])
        bubble_rows, dew_rows, z1, crit_point = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals)
        qlines = compute_quality_lines(bubble_rows, dew_rows)
        # qlines: dict of quality -> list of {"T_K","P_Pa","h_Jmol",...} rows;
        # pick the row nearest T_SAT at quality 0.5.
        rows = qlines[0.5]
        idx = int(np.argmin([abs(r["T_K"] - T_SAT) for r in rows]))
        return rows[idx]
    results.append(_repeatability(
        "quality_0p5_at_Tsat",
        quality_state_fn,
        {"P_Pa": lambda r: r["P_Pa"], "h_Jmol": lambda r: r["h_Jmol"], "T_K": lambda r: r["T_K"]},
    ))

    # 5/6. Supercritical state + critical-point solve. solve_mixture_critical_
    # point requires explicit t_guess/rho_guess (no defaults) -- reuse the
    # oracle's OWN fallback guess recipe (Bell-mixing-rule reducing point,
    # rho_guess=1/vred_mix) from run_true_vle_envelope's own "no usable
    # sweep data" branch, since a standalone call here has no sweep to seed
    # from either.
    tc1_g, tc2_g = float(d1["basic"]["Tc"]), float(d2["basic"]["Tc"])
    vc1_g, vc2_g = 1.0 / rhoc1, 1.0 / rhoc2
    from linear_model_codex import bell2023_Tred_vred, BELL_2023_R1234ZE_R227EA
    tred_mix_g, vred_mix_g = bell2023_Tred_vred(x1, 1.0 - x1, tc1_g, tc2_g, vc1_g, vc2_g, BELL_2023_R1234ZE_R227EA)

    results.append(_repeatability(
        "critical_point_solve",
        lambda: solve_mixture_critical_point(d1, d2, x1, t_guess=tred_mix_g, rho_guess=1.0 / vred_mix_g),
        {"T_K": lambda r: r["T_K"], "P_Pa": lambda r: r["P_Pa"], "rho_molm3": lambda r: r["rho_molm3"]},
    ))
    crit = solve_mixture_critical_point(d1, d2, x1, t_guess=tred_mix_g, rho_guess=1.0 / vred_mix_g)
    Tc = crit["T_K"]
    T_SC = Tc + 8.0
    RHO_SC = 3000.0
    results.append(_repeatability(
        "supercritical_direct",
        lambda: mix_state(d1, d2, T_SC, RHO_SC, x1),
        {"p_pa": lambda r: r.p_pa, "h_jmol": lambda r: r.h_jmol, "g_jmol": lambda r: r.g_jmol},
    ))

    print("\n=== Reference repeatability (oracle: mixture_fully_validated.py) ===")
    for r in results:
        print(f"\n[{r['label']}] ({r['n_runs']} runs)")
        for name, stats in r["fields"].items():
            print(f"  {name}: mean={stats['mean']:.10g} spread_abs={stats['abs_spread']:.3g} "
                  f"spread_rel={stats['rel_spread']:.3g}")

    with open(OUT_JSON, "w") as f:
        json.dump({"x1": x1, "w1": W1, "n_repeats": N_REPEATS, "results": results}, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")

    # Determine whether ANY measured spread is non-zero -- if the oracle is
    # exactly bit-reproducible (expected for pure deterministic floating-point
    # evaluation with no RNG), frozen tolerances are then set by a fixed,
    # conservative floor rather than by the (zero) measured spread itself.
    max_rel_spread = max(
        stats["rel_spread"]
        for r in results
        for stats in r["fields"].values()
    )
    print(f"\nMax relative spread observed across all fields/states: {max_rel_spread:.3g}")
    if max_rel_spread == 0.0:
        print("Oracle is exactly bit-reproducible across all tested states (as expected for "
              "deterministic floating-point code with no RNG). Frozen tolerances will be set "
              "by conservative fixed floors, not by measured spread.")
    else:
        print("WARNING: non-zero repeatability spread detected -- inspect "
              "reference_repeatability_results.json before freezing tolerances.")


if __name__ == "__main__":
    sys.exit(main())
