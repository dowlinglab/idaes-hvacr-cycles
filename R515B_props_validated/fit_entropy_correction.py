"""
Fit a T-dependent entropy correction curve against 10 (s, P, H) points
hand-digitized from the Honeywell Solstice N15 (R-515B) TDS p-h chart.

Purpose
-------
_mix_entropy_direct() in mixture_isentrope_validation.py currently applies a
single FLAT offset (ENTROPY_REFERENCE_OFFSET_JMOLK = -1.4595 J/(mol*K)),
calibrated at exactly one point: Honeywell's stated reference state
(T=273.15K, sat. liq. at 0C). That flat offset is exactly right AT that one
point but is not expected to stay right elsewhere, especially near the
critical point where the model's composition-extrapolation bias (Bell 2023
departure function fit for x1=0.33-0.68 vs R-515B's real x1~0.9385) is
largest -- this is the confirmed root cause of the isentrope-vs-Honeywell
mismatch (2026-08-14 breadcrumb).

This script builds a genuine T-dependent correction:
1. For each of the 10 Honeywell (s_label, P, H) points, solve the model for
   the (T, rho) state whose OWN P and H match the point's P and H exactly
   (a 2-equation/2-unknown solve, same pattern as the existing isentrope
   P/s solve in mixture_isentrope_validation.py, just swapping the second
   target from s to h). This pins down "which state is Honeywell actually
   talking about" using only P and H -- quantities the model gets right
   without needing any entropy correction at all.
2. At that solved state, read off the model's own RAW (pre-offset) entropy
   and compare it against Honeywell's stated s_label for that point. The
   residual (s_honeywell - s_model_raw) is the correction the model is
   missing AT THAT STATE's temperature.
3. Fit residual(T) with both a line (a + b*T) and a quadratic
   (a + b*T + c*T^2) across all 10 points and compare RMS error, to check
   whether the extra curvature is actually justified by the data or would
   just be fitting noise (visual cue from the full reference chart: the
   isentropes themselves look close to straight lines, which is weak prior
   evidence -- but not proof -- that the residual is smooth/low-order too).

This script is READ-ONLY / diagnostic -- it does not modify
mixture_isentrope_validation.py. Once the fit is confirmed, its output
coefficients get hand-wired into a replacement for the flat
ENTROPY_REFERENCE_OFFSET_JMOLK constant there.

Usage
-----
Run from R515B_props_validated/ (same directory as
mixture_isentrope_validation.py and its own copy of linear_model_codex.py):
    python3 fit_entropy_correction.py
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import root

from mixture_isentrope_validation import (
    BTU_LBMR_TO_JKGK,
    BTU_LBM_TO_KJ_KG,
    ENTROPY_REFERENCE_OFFSET_JMOLK,
    PSI_TO_PA,
    load_idaes_helmholtz_json,
    mix_state,
    mw_from_json,
    w1_to_x1,
    _mix_entropy_direct,
)

FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911

# Honeywell's stated critical pressure is 507 psig (PHYSICAL PROPERTIES
# table convention) -> 521.696 psia absolute (507 + 14.696), confirmed in
# the 2026-08-14 pressure-unit breadcrumb entry. Point #4 below uses this
# literal value ("Pc" on the chart, not a numeric placeholder).
PC_PSIA = 507.0 + 14.696

# All 15 hand-digitized (s, P, H) points confirmed as of 2026-08-14
# (kept name POINTS_10 for the original set for continuity with earlier
# breadcrumb entries; POINTS_LIQUID_BATCH is the second batch of 6 points,
# with the #11/#16 duplicate collapsed to one row after the user confirmed
# both were the same point, "0.22.450,76" read as (0.22, 450, 76), with the
# second (0.22, 450, 86) a parallax misread of the same point).
# Columns: s [Btu/(lbm-R)], P [psia], H [Btu/lbm]
POINTS_10: List[Tuple[float, float, float]] = [
    (0.22, 12.0, 76.5),
    (0.43, 90.0, 194.0),
    (0.32, 1050.0, 132.0),
    (0.37, PC_PSIA, 168.0),
    (0.26, 60.0, 96.0),
    (0.28, 150.0, 106.0),
    (0.51, 15.0, 228.0),
    (0.49, 30.0, 222.0),
    (0.43, 60.0, 190.0),
    (0.37, 1200.0, 173.0),
]

# Second batch (2026-08-14, later): user-supplied, all low-H (76-88
# Btu/lbm) -- expected to land on the liquid/dense branch, given to
# reinforce that side of the fit. Duplicate (#11/#16, both (0.22,450,*))
# collapsed to a single row per user confirmation ("It is actually 76, I
# had a parallax before").
POINTS_LIQUID_BATCH: List[Tuple[float, float, float]] = [
    (0.22, 450.0, 76.0),
    (0.22, 75.0, 76.0),
    (0.22, 30.0, 76.0),
    (0.24, 1200.0, 88.0),
    (0.24, 1050.0, 88.0),  # corrected 2026-08-14: was misread as H=78, user confirmed H=88
]

ALL_POINTS: List[Tuple[float, float, float]] = POINTS_10 + POINTS_LIQUID_BATCH

# Pragmatic density-based branch split for fitting the entropy correction
# separately by region: NOT a rigorous phase check (that would need the
# actual mixture critical density from solve_mixture_critical_point), just
# a fast, defensible threshold near the mixture's known critical density
# ballpark to separate "dense/liquid-like branch" from "light/vapor-like or
# lean-supercritical branch" -- the two clusters the first 10-point fit
# already showed behaving differently.
BRANCH_DENSITY_THRESHOLD_KGM3 = 500.0

# Model's own solved mixture critical temperature (matches the comment
# already in this file at VAPOR_ONLY_ISOTHERM_VALUES_F: "solved mixture
# Tc=228.0F (382.045K, 0.04% off Honeywell)"). Used to draw the actual
# "liquid phase up to critical point" cutoff requested by the user --
# tighter/more correct than the density threshold alone, since a state can
# be dense (>500 kg/m3) while still being supercritical (T > Tc) at high
# enough P, e.g. point #10 below (T solves to ~409K/276F, well above Tc,
# despite rho=745 kg/m3 clearing the density threshold).
TC_MIX_K = 382.045


def _ph_match_residual(vars_, d1: Dict, d2: Dict, z1: float, p_target_pa: float, h_target_jmol: float) -> List[float]:
    """
    Residual for the 2-equation/2-unknown (P, H)-matching solve: given a
    target pressure and target enthalpy (both quantities the model already
    gets right, per the 2026-08-14 finding that only entropy needs
    correcting), solve jointly for (T, rho). Mirrors
    _isentrope_2eq_residual's normalized-residual pattern (dimensionless,
    so hybr's numerical Jacobian stays well-conditioned across the Pa-scale
    vs J/mol-scale mismatch).
    """
    t_k, rho = vars_
    if t_k <= 0 or rho <= 0:
        return [1.0e12, 1.0e12]
    st = mix_state(d1, d2, t_k, rho, z1)
    r_p = (st.p_pa - p_target_pa) / p_target_pa
    r_h = (st.h_jmol - h_target_jmol) / h_target_jmol
    return [r_p, r_h]


def solve_ph_state(d1: Dict, d2: Dict, z1: float, p_target_pa: float, h_target_jmol: float, mw_mix_kgmol: float) -> Tuple[float, float]:
    """
    Multi-seed scan for the (T, rho) state matching (p_target_pa,
    h_target_jmol). Unlike the isentrope extensions, these 10 points have no
    single natural anchor (no shared T-grid row to warm-start from) -- they
    span liquid, near-critical, and supercritical states at very different
    conditions -- so this tries a spread of trial temperatures with both a
    liquid-like and an ideal-gas-like density seed at each, and returns the
    first converged, physically sensible solution.

    Raises RuntimeError if no seed converges.
    """
    r_u = 8.31446261815324  # J/(mol*K), only used for the ideal-gas seed estimate
    trial_t_f = [-20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380]
    candidates: List[Tuple[float, float]] = []
    for t_f in trial_t_f:
        t_k = (t_f - 32.0) * 5.0 / 9.0 + 273.15
        # Liquid-like seed: ~1200 kg/m3 order of magnitude for this mixture
        # (matches the RHO_REF_KGM3=1258.4 reference state used elsewhere),
        # converted to mol/m3 via the actual mixture molar mass.
        rho_liquid_seed = 1200.0 / mw_mix_kgmol
        rho_vapor_seed = p_target_pa / (r_u * t_k)
        candidates.append((t_k, rho_vapor_seed))
        candidates.append((t_k, rho_liquid_seed))

    for t_k0, rho0 in candidates:
        try:
            sol = root(
                _ph_match_residual,
                x0=[t_k0, rho0],
                args=(d1, d2, z1, p_target_pa, h_target_jmol),
                method="hybr",
                tol=1.0e-11,
            )
        except Exception:
            continue
        if not sol.success:
            continue
        t_sol, rho_sol = float(sol.x[0]), float(sol.x[1])
        if t_sol <= 0 or rho_sol <= 0:
            continue
        # Confirm residuals are actually tight (hybr can report success on
        # a loose stall) before accepting.
        r_p, r_h = _ph_match_residual([t_sol, rho_sol], d1, d2, z1, p_target_pa, h_target_jmol)
        if abs(r_p) < 1.0e-6 and abs(r_h) < 1.0e-6:
            return t_sol, rho_sol

    raise RuntimeError(f"No seed converged for P={p_target_pa:.1f} Pa, H={h_target_jmol:.1f} J/mol")


def _fit_and_report(label: str, rows: List[Dict]) -> None:
    """Fit linear + quadratic residual(T) curves to one branch's rows and print a recommendation."""
    print(f"\n--- {label}: {len(rows)} point(s) ---")
    if len(rows) < 2:
        print("  Too few points to fit a line -- skipping.")
        return

    t_arr = np.array([r["t_k"] for r in rows])
    resid_arr = np.array([r["residual"] for r in rows])

    coeffs_lin = np.polyfit(t_arr, resid_arr, deg=1)
    resid_pred_lin = np.polyval(coeffs_lin, t_arr)
    rms_lin = float(np.sqrt(np.mean((resid_arr - resid_pred_lin) ** 2)))
    b_lin, a_lin = coeffs_lin
    print("  Linear fit:    residual(T[K]) = {:+.6f} + ({:+.8f})*T".format(a_lin, b_lin))
    print(f"    RMS error: {rms_lin:.4f} J/(mol*K)")

    if len(rows) < 4:
        print("  (fewer than 4 points -- skipping quadratic, would just interpolate/overfit)")
        return

    coeffs_quad = np.polyfit(t_arr, resid_arr, deg=2)
    resid_pred_quad = np.polyval(coeffs_quad, t_arr)
    rms_quad = float(np.sqrt(np.mean((resid_arr - resid_pred_quad) ** 2)))
    c_quad, b_quad, a_quad = coeffs_quad
    print("  Quadratic fit: residual(T[K]) = {:+.6f} + ({:+.8f})*T + ({:+.10f})*T^2".format(a_quad, b_quad, c_quad))
    print(f"    RMS error: {rms_quad:.4f} J/(mol*K)")

    improvement_pct = 100.0 * (rms_lin - rms_quad) / rms_lin if rms_lin > 0 else 0.0
    print(f"  Quadratic reduces RMS by {improvement_pct:.1f}% vs linear.")
    if improvement_pct < 15.0:
        print(f"  Recommendation: LINEAR -- quadratic's extra term isn't earning its keep on {len(rows)} points.")
    else:
        print("  Recommendation: QUADRATIC -- meaningful RMS reduction.")


def _fit_anchored_constrained(label: str, rows: List[Dict], t_anchor_k: float, resid_anchor: float) -> None:
    """
    Fit residual(T) = resid_anchor + b*(T-t_anchor) + c*(T-t_anchor)^2,
    i.e. a quadratic FORCED to pass exactly through (t_anchor, resid_anchor)
    by construction (b, c solved via ordinary least-squares on the
    remaining points' deviations from the anchor value; the anchor itself
    contributes no residual to that least-squares problem since it's
    satisfied exactly, not fit).

    This replaces the earlier "anchored" attempt, which just appended the
    exact reference-state row as an 11th ordinary data point to
    np.polyfit(deg=2) -- that MINIMIZES total squared error across all 11
    points but does NOT force exact pass-through any single one of them,
    including the anchor. Confirmed 2026-08-14: that version predicted
    -0.2804 J/(mol*K) at T=273.15K, a +1.1791 J/(mol*K) miss off the exact
    -1.4595 anchor value (equivalent to ~1% error in absolute entropy at
    Honeywell's own stated reference state) -- defeats the entire purpose
    of adding the anchor. This constrained version fixes that: by
    construction, plugging t_anchor_k back in always returns exactly
    resid_anchor, to floating-point precision.
    """
    print(f"\n--- {label}: {len(rows)} point(s) + 1 exact anchor (constrained) ---")
    t_arr = np.array([r["t_k"] for r in rows])
    resid_arr = np.array([r["residual"] for r in rows])
    dt = t_arr - t_anchor_k
    y = resid_arr - resid_anchor  # deviation from the anchor value; anchor itself has dt=0, y=0 by definition, not included as a row

    # Linear-in-dt fit (still forces value=resid_anchor at dt=0 since there's no intercept term): y = b*dt
    b_lin = float(np.sum(dt * y) / np.sum(dt * dt))
    resid_pred_lin = resid_anchor + b_lin * dt
    rms_lin = float(np.sqrt(np.mean((resid_arr - resid_pred_lin) ** 2)))
    print(f"  Linear (anchored):    residual(T[K]) = {resid_anchor:+.6f} + ({b_lin:+.8f})*(T-{t_anchor_k:.2f})")
    print(f"    RMS error (all {len(rows)} chart points): {rms_lin:.4f} J/(mol*K)")

    # Quadratic-in-dt fit: y = b*dt + c*dt^2 (no intercept -- forced 0 at dt=0)
    design = np.column_stack([dt, dt * dt])
    coeffs, *_ = np.linalg.lstsq(design, y, rcond=None)
    b_quad, c_quad = float(coeffs[0]), float(coeffs[1])
    resid_pred_quad = resid_anchor + b_quad * dt + c_quad * dt * dt
    rms_quad = float(np.sqrt(np.mean((resid_arr - resid_pred_quad) ** 2)))
    print(f"  Quadratic (anchored): residual(T[K]) = {resid_anchor:+.6f} + ({b_quad:+.8f})*(T-{t_anchor_k:.2f}) + ({c_quad:+.10f})*(T-{t_anchor_k:.2f})^2")
    print(f"    RMS error (all {len(rows)} chart points): {rms_quad:.4f} J/(mol*K)")

    improvement_pct = 100.0 * (rms_lin - rms_quad) / rms_lin if rms_lin > 0 else 0.0
    print(f"  Quadratic reduces RMS by {improvement_pct:.1f}% vs linear.")
    rec = "QUADRATIC" if improvement_pct >= 15.0 else "LINEAR"
    print(f"  Recommendation: {rec}")

    # Sanity check: confirm exact pass-through at the anchor (should be ~0.0, floating point only)
    check = (resid_anchor + b_quad * 0.0 + c_quad * 0.0 * 0.0) - resid_anchor
    print(f"  Sanity check -- value at anchor minus resid_anchor (should be 0): {check:.2e}")

    # Also report in EXPANDED a+b*T+c*T^2 form for direct use as a drop-in
    # replacement in _mix_entropy_direct (equivalent function, just not
    # written in the (T-t_anchor) form -- expand (T-t0)=u:
    # a0 + b0*u + c0*u^2, u=T-t0 => a0 - b0*t0 + c0*t0^2 + (b0-2*c0*t0)*T + c0*T^2
    a_exp = resid_anchor - b_quad * t_anchor_k + c_quad * t_anchor_k ** 2
    b_exp = b_quad - 2.0 * c_quad * t_anchor_k
    c_exp = c_quad
    print(f"  Expanded (a+b*T+c*T^2) form for direct wiring: a={a_exp:+.6f}, b={b_exp:+.8f}, c={c_exp:+.10f}")


def main() -> None:
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = w1_to_x1(W1, mw1, mw2)
    mw_mix = z1 * mw1 + (1.0 - z1) * mw2  # kg/mol

    print(f"z1 (R-1234ze(E) mole fraction) = {z1:.6f}, mw_mix = {mw_mix*1000:.4f} g/mol\n")
    print(f"Current flat offset: ENTROPY_REFERENCE_OFFSET_JMOLK = {ENTROPY_REFERENCE_OFFSET_JMOLK:+.4f} J/(mol*K)\n")
    print(f"Branch split: dense/liquid-like if solved mass density > {BRANCH_DENSITY_THRESHOLD_KGM3:.0f} kg/m3, else light/vapor-supercritical.\n")

    rows = []
    header = f"{'#':>2} {'s_label':>8} {'P_psia':>9} {'H_btulbm':>9} {'T_solved_F':>11} {'rho_kgm3':>9} {'branch':>10} {'s_honeywell':>12} {'s_model_raw':>12} {'residual':>10}"
    print(header)
    for i, (s_label, p_psia, h_btulbm) in enumerate(ALL_POINTS, start=1):
        p_pa = p_psia * PSI_TO_PA
        h_jmol = h_btulbm * BTU_LBM_TO_KJ_KG * 1000.0 * mw_mix
        s_honeywell_jmolK = s_label * BTU_LBMR_TO_JKGK * mw_mix

        try:
            t_k, rho_mol = solve_ph_state(d1, d2, z1, p_pa, h_jmol, mw_mix)
        except RuntimeError as exc:
            print(f"{i:>2} {s_label:>8.2f} {p_psia:>9.1f} {h_btulbm:>9.1f}   FAILED: {exc}")
            continue

        rho_kgm3 = rho_mol * mw_mix
        branch = "liquid" if rho_kgm3 > BRANCH_DENSITY_THRESHOLD_KGM3 else "vapor/SC"

        s_model_with_offset = _mix_entropy_direct(d1, d2, t_k, rho_mol, z1)
        s_model_raw = s_model_with_offset - ENTROPY_REFERENCE_OFFSET_JMOLK
        residual = s_honeywell_jmolK - s_model_raw

        t_f = (t_k - 273.15) * 9.0 / 5.0 + 32.0
        print(f"{i:>2} {s_label:>8.2f} {p_psia:>9.1f} {h_btulbm:>9.1f} {t_f:>11.2f} {rho_kgm3:>9.1f} {branch:>10} {s_honeywell_jmolK:>12.4f} {s_model_raw:>12.4f} {residual:>10.4f}")
        rows.append({"i": i, "t_k": t_k, "t_f": t_f, "residual": residual, "branch": branch})

    print(f"\n{len(rows)} of {len(ALL_POINTS)} points converged.")

    liquid_rows = [r for r in rows if r["branch"] == "liquid"]
    vapor_rows = [r for r in rows if r["branch"] == "vapor/SC"]

    _fit_and_report("LIQUID / DENSE branch (density-only split, reference)", liquid_rows)
    _fit_and_report("VAPOR / SUPERCRITICAL branch", vapor_rows)

    # Authoritative fit for wiring into mixture_isentrope_validation.py:
    # "liquid phase up to critical point" per user request 2026-08-14 --
    # dense branch AND T <= Tc_mix (excludes e.g. point #10, dense but
    # T~409K/276F, well above Tc=382.045K/228F -- that's supercritical
    # despite clearing the density threshold).
    liquid_up_to_tc_rows = [r for r in liquid_rows if r["t_k"] <= TC_MIX_K + 1.0]
    excluded = [r for r in liquid_rows if r not in liquid_up_to_tc_rows]
    if excluded:
        print(f"\nExcluded from the liquid-up-to-Tc fit (T > Tc={TC_MIX_K:.3f}K despite dense branch): "
              + ", ".join(f"#{r['i']} (T={r['t_f']:.1f}F)" for r in excluded))

    # Anchor the fit with the EXACT reference-state calibration point
    # (T=273.15K, sat. liq. at 0C) -- this one has NO chart-digitization
    # uncertainty at all (it comes directly from Honeywell's printed
    # footnote "s=1.00 kJ/kg-K" compared against the model's own exactly
    # computable raw entropy there, already validated in
    # check_entropy_reference.py), unlike the 10 hand-read (s,P,H) points.
    # Including it keeps the fitted curve consistent with the one anchor we
    # have zero uncertainty on, rather than letting a sparse quadratic fit
    # drift away from it.
    print(f"\n(Reference-state anchor: T=273.15K, residual={ENTROPY_REFERENCE_OFFSET_JMOLK:+.4f} J/(mol*K) exact, no chart-read uncertainty)")

    _fit_and_report("LIQUID branch, T<=Tc ONLY, density-branch points only (no ref anchor)", liquid_up_to_tc_rows)
    # NOTE 2026-08-14: the naive "append the anchor as an 11th np.polyfit
    # row" version previously reported here did NOT force exact pass-through
    # at the anchor (confirmed: predicted -0.2804 vs the exact -1.4595 target,
    # a +1.1791 J/(mol*K) miss) -- replaced with the constrained fit below,
    # which forces exact pass-through by construction. This is now the
    # authoritative fit for wiring in.
    _fit_anchored_constrained(
        "LIQUID branch, T<=Tc, CONSTRAINED through exact ref-state anchor (AUTHORITATIVE, for wiring in)",
        liquid_up_to_tc_rows, t_anchor_k=273.15, resid_anchor=ENTROPY_REFERENCE_OFFSET_JMOLK,
    )


if __name__ == "__main__":
    main()
