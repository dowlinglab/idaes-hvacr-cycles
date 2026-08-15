"""
Standalone Honeywell-density/pressure correction layer for the R-515B
(R-1234ze/R-227ea) mixture model.

Purpose
-------
The Bell (2023) R-1234ze(E)/R-227ea departure function was only ever fit and
validated for x1=0.33-0.68 (see PROJECT_CONTEXT.md, 2026-08-12 entries), while
the real Honeywell R-515B blend sits at x1~0.9385 -- a genuine composition
extrapolation, not a coding bug (every parameter/formula in the mixing/
departure layer was independently re-verified against the paper). This shows
up as a confirmed ~8.8% critical-density error and smaller-but-real errors in
saturated liquid/vapor density elsewhere (see PROJECT_CONTEXT.md, 2026-08-14
entries "Honeywell's stated critical properties" and "Critical-density error
connected to the known 2026-08-12 composition-extrapolation root cause").

Rather than continuing to chase this inside the Bell (2023) departure
function itself, this module builds a PURE POST-PROCESSING correction layer:
it digitizes Honeywell's own stated saturation density/pressure values (from
the TDS "PHYSICAL PROPERTIES" table) and fits a smooth residual curve
(Honeywell value - model value) vs. temperature, separately for liquid
density, vapor density, and saturation pressure. The underlying EOS/solver
(mix_state, solve_bubble_at_t, solve_dew_at_t, solve_mixture_critical_point,
etc.) is NOT modified anywhere by this file -- it is only ever called, never
edited, and this module's corrected values are reported as EXTRA fields
alongside (never replacing) the model's own raw output.

Data source: Honeywell Solstice N15 (R-515B) TDS, "PHYSICAL PROPERTIES"
table (screenshotted 2026-08-14). 5 stated points, all sharing the fixed
z1=0.9385-ish (w1=0.911 mass fraction) blend composition:
  - Boiling temperature @ 0 psig = -2.0 F        (defines the vapor P anchor)
  - Vapor Density at 0 psig Boiling Point = 0.367 lbm/ft3
  - Liquid Density at 32 F   = 78.56 lbm/ft3
  - Liquid Density at 77 F   = 73.65 lbm/ft3
  - Vapor Density at 77 F    = 1.69 lbm/ft3
  - Saturated Pressure at 77 F = 57.45 psig
  - Critical temperature = 228.0 F, Critical pressure = 507 psig,
    Critical density = 31.03 lbm/ft3

Status: with only 3 points per curve, a degree-2 (quadratic) fit passes
through all 3 exactly -- this is NOT yet a validated functional form, just
the simplest smooth placeholder (see PROJECT_CONTEXT.md discussion on why
linear is likely wrong -- classical/mean-field EOS critical exponent mismatch
vs. real fluids -- and why quadratic can't yet be distinguished from other
3-parameter forms with only 3 data points). Expected to be refit once more
Honeywell data points are digitized.

UPDATE 2026-08-14: the quadratic fit was found to overshoot BETWEEN anchors
(vapor-density residual dips to -2.44 kg/m3 around T=277K, well past either
neighboring data point, with zero supporting data -- a pure artifact of
forcing one parabola to also satisfy the huge +43.75 kg/m3 jump at the
critical point). Both `_fit_residual_curve` (quadratic, kept for comparison)
and `_fit_pchip_curve` (shape-preserving PCHIP, added to fix the overshoot)
are available; `main()` prints both side by side. PCHIP is the recommended
choice until more data points are digitized -- it passes through the same 3
known points exactly, but doesn't invent unsupported behavior in between.

Usage
-----
Run from R515B_props_validated/ (same directory as
mixture_isentrope_validation.py and its own copy of linear_model_codex.py):
    python3 honeywell_saturation_correction.py

Read-only with respect to mixture_isentrope_validation.py -- only imports
and calls existing functions, never modifies them.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

from mixture_isentrope_validation import (
    run_true_vle_envelope,
    solve_bubble_at_t,
    solve_dew_at_t,
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
)

# ---- Fixed blend definition, matches every other script in this folder ----
FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
TMIN, TMAX, N = 255.4, 380.0, 80  # normal sweep range, used only to seed the exact-T solves below

# ---- Honeywell's digitized anchor points (TDS "PHYSICAL PROPERTIES" table, 2026-08-14) ----
# Liquid branch: 32 F, 77 F, critical (228.0 F)
HONEYWELL_LIQUID_ANCHORS_TK = [273.15, 298.15, 382.0389]
HONEYWELL_LIQUID_ANCHORS_KGM3 = [1258.4105, 1179.7598, 497.0529]

# Vapor branch: -2.0 F (boiling pt @ 0 psig), 77 F, critical (228.0 F)
HONEYWELL_VAPOR_ANCHORS_TK = [254.2611, 298.15, 382.0389]
HONEYWELL_VAPOR_ANCHORS_KGM3 = [5.8788, 27.0712, 497.0529]

# Saturation pressure: -2.0 F (=0 psig by definition), 77 F (stated 57.45 psig), critical (507 psig)
HONEYWELL_PSAT_ANCHORS_TK = [254.2611, 298.15, 382.0389]
HONEYWELL_PSAT_ANCHORS_PA = [101325.4, 497429.2, 3596967.3]


def _nearest_row(rows, target_t_k):
    """Nearest CONVERGED row to target_t_k -- used only as a solver seed, not as the answer itself."""
    candidates = [r for r in rows if r["status"] == "CONVERGED"]
    if not candidates:
        raise RuntimeError("No converged rows available to seed from.")
    return min(candidates, key=lambda r: abs(r["T_K"] - target_t_k))


def solve_model_at_anchors():
    """
    Runs the normal sweep once (for seeds), then solves the model's OWN
    bubble/dew states at the EXACT Honeywell anchor temperatures (not
    whatever happens to land on the T-grid). Returns three dicts keyed by
    anchor T, each mapping to the model's own value at that T:
        rho_liq_model_kgm3[T], rho_vap_model_kgm3[T], p_model_pa[T]
    (p_model_pa is populated at whichever anchor T's the pressure anchors
    need -- 254.2611, 298.15, and the critical point.)
    """
    t_vals = np.linspace(TMIN, TMAX, N)
    print(f"Running dome sweep ({FLUID1}/{FLUID2}, w1={W1}, T={TMIN}-{TMAX} K, n={N}) for seeds...")
    bubble_rows, dew_rows, z1, crit_point = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals)

    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = z1 * mw1 + (1.0 - z1) * mw2  # kg/mol

    rho_liq_model_kgm3 = {}
    rho_vap_model_kgm3 = {}
    p_model_pa = {}
    tc_model = None  # the model's OWN solved Tc -- differs slightly from Honeywell's stated 382.0389 K,
                     # so this (not the Honeywell constant) is the correct x-coordinate to fit against,
                     # since the correction curve gets evaluated later in the model's own T-coordinate space

    print("\n=== Liquid-branch exact-T solves (solve_bubble_at_t) ===")
    for t_k in HONEYWELL_LIQUID_ANCHORS_TK:
        if abs(t_k - 382.0389) < 1e-3:
            continue  # critical point handled separately below
        seed = _nearest_row(bubble_rows, t_k)
        row = solve_bubble_at_t(d1, d2, t_k, z1, seed["rho_l_molm3"], seed["rho_v_molm3"], seed["y1_vap"])
        if row["status"] != "CONVERGED":
            print(f"  T={t_k:.4f} K: DID NOT CONVERGE ({row['status']}) -- seed was T={seed['T_K']:.2f} K")
            continue
        rho_kgm3 = row["rho_l_molm3"] * mw_mix
        rho_liq_model_kgm3[t_k] = rho_kgm3
        p_model_pa[t_k] = row["P_Pa"]
        print(f"  T={t_k:.4f} K: rho_l = {rho_kgm3:.4f} kg/m3   P = {row['P_Pa']:.1f} Pa   (seed T={seed['T_K']:.2f} K)")

    print("\n=== Vapor-branch exact-T solves (solve_dew_at_t) ===")
    for t_k in HONEYWELL_VAPOR_ANCHORS_TK:
        if abs(t_k - 382.0389) < 1e-3:
            continue  # critical point handled separately below
        seed = _nearest_row(dew_rows, t_k)
        row = solve_dew_at_t(d1, d2, t_k, z1, seed["rho_l_molm3"], seed["rho_v_molm3"], seed["x1_liq"])
        if row["status"] != "CONVERGED":
            print(f"  T={t_k:.4f} K: DID NOT CONVERGE ({row['status']}) -- seed was T={seed['T_K']:.2f} K")
            continue
        rho_kgm3 = row["rho_v_molm3"] * mw_mix
        rho_vap_model_kgm3[t_k] = rho_kgm3
        p_model_pa[t_k] = row["P_Pa"]  # will overwrite the 298.15 K liquid-side P if also solved there -- see cross-check print below
        print(f"  T={t_k:.4f} K: rho_v = {rho_kgm3:.4f} kg/m3   P = {row['P_Pa']:.1f} Pa   (seed T={seed['T_K']:.2f} K)")

    print("\n=== Critical point (already solved, from solve_mixture_critical_point) ===")
    if crit_point is not None and crit_point.get("converged", False):
        tc = crit_point["T_K"]
        rho_c_kgm3 = crit_point["rho_molm3"] * mw_mix
        pc = crit_point["P_Pa"]
        rho_liq_model_kgm3[tc] = rho_c_kgm3
        rho_vap_model_kgm3[tc] = rho_c_kgm3  # same point, shared by both branches
        p_model_pa[tc] = pc
        tc_model = tc
        print(f"  Tc={tc:.4f} K: rho_c = {rho_c_kgm3:.4f} kg/m3   Pc = {pc:.1f} Pa"
              f"   (Honeywell states Tc=382.0389 K -- {tc - 382.0389:+.4f} K different; "
              f"using the MODEL's own Tc={tc:.4f} K as this anchor's T-coordinate for fitting, "
              f"paired with Honeywell's critical density/pressure VALUES)")
    else:
        print("  WARNING: critical point did not converge -- critical anchor unavailable.")

    # Cross-check: at T=298.15K we solved P from BOTH a bubble and a dew call.
    # These should be close (R-515B is called an azeotropic blend -- minimal glide).
    bubble_seed_298 = _nearest_row(bubble_rows, 298.15)
    dew_seed_298 = _nearest_row(dew_rows, 298.15)
    row_b298 = solve_bubble_at_t(d1, d2, 298.15, z1, bubble_seed_298["rho_l_molm3"], bubble_seed_298["rho_v_molm3"], bubble_seed_298["y1_vap"])
    row_d298 = solve_dew_at_t(d1, d2, 298.15, z1, dew_seed_298["rho_l_molm3"], dew_seed_298["rho_v_molm3"], dew_seed_298["x1_liq"])
    if row_b298["status"] == "CONVERGED" and row_d298["status"] == "CONVERGED":
        print(f"\nCross-check at T=298.15K: bubble-solve P = {row_b298['P_Pa']:.1f} Pa, "
              f"dew-solve P = {row_d298['P_Pa']:.1f} Pa "
              f"(diff = {row_b298['P_Pa'] - row_d298['P_Pa']:+.2f} Pa, "
              f"{100*(row_b298['P_Pa']-row_d298['P_Pa'])/row_b298['P_Pa']:+.4f}% -- "
              f"should be small for a near-azeotropic blend)")
        p_model_pa[298.15] = row_b298["P_Pa"]  # use the bubble-solve value as the canonical 298.15K pressure anchor

    if tc_model is None:
        raise RuntimeError("Critical point did not converge -- cannot build the critical anchor for fitting.")

    return rho_liq_model_kgm3, rho_vap_model_kgm3, p_model_pa, mw_mix, tc_model


def _fit_residual_curve(anchor_t_list, honeywell_val_list, model_val_dict, degree=2):
    """
    residual(T) = honeywell - model, fit as a degree-`degree` polynomial
    through the anchor points. Returns (poly1d callable, list of raw
    residuals actually used in the fit) -- the raw residuals are returned
    for the sanity check that the fitted curve reproduces them exactly.
    """
    model_vals = [model_val_dict[t] for t in anchor_t_list]
    residuals = [h - m for h, m in zip(honeywell_val_list, model_vals)]
    coeffs = np.polyfit(anchor_t_list, residuals, degree)
    return np.poly1d(coeffs), residuals


def _fit_pchip_curve(anchor_t_list, honeywell_val_list, model_val_dict):
    """
    Shape-preserving (monotone) interpolation through the same anchor points
    used by _fit_residual_curve, via scipy's PCHIP (Piecewise Cubic Hermite
    Interpolating Polynomial). Unlike the raw quadratic fit above, PCHIP is
    constructed so the curve never overshoots past its neighboring data
    points -- added 2026-08-14 after finding the quadratic vapor-density fit
    dips to -2.44 kg/m3 around T=277K, well past either neighboring anchor
    (-0.22 at 254K, -0.59 at 298K) with zero supporting data (see
    PROJECT_CONTEXT.md, "did we fit correctly" discussion). Same anchor
    points, same residuals -- just a different, overshoot-safe interpolation
    method. Requires the anchor T's sorted ascending (handled internally).
    """
    model_vals = [model_val_dict[t] for t in anchor_t_list]
    residuals = [h - m for h, m in zip(honeywell_val_list, model_vals)]
    order = np.argsort(anchor_t_list)
    t_sorted = np.asarray(anchor_t_list)[order]
    r_sorted = np.asarray(residuals)[order]
    return PchipInterpolator(t_sorted, r_sorted), residuals


def honeywell_corrected_saturation_state(t_k, phase, rho_model_kgm3, p_model_pa,
                                          rho_liq_correction_fn, rho_vap_correction_fn, p_correction_fn):
    """
    Pure post-processing -- does NOT touch mix_state()/the EOS solve. Call
    this on top of an already-solved bubble/dew row to get a
    Honeywell-corrected (rho, P) pair for reporting/plotting, alongside
    (not replacing) the raw solver output that everything else internally
    still uses.

    phase: "liquid" or "vapor"

    NOTE: this correction is calibrated ONLY on the saturation curve (bubble/
    dew lines + critical point). Applying it to off-dome states (subcooled
    liquid, superheated vapor, isentrope/isotherm extensions) would be an
    extrapolation of THIS correction itself, not yet validated -- do not use
    this function outside bubble_rows/dew_rows/crit_point without revisiting
    this assumption first.
    """
    rho_corr_fn = rho_liq_correction_fn if phase == "liquid" else rho_vap_correction_fn
    return (rho_model_kgm3 + rho_corr_fn(t_k),
            p_model_pa + p_correction_fn(t_k))


def main():
    rho_liq_model_kgm3, rho_vap_model_kgm3, p_model_pa, mw_mix, tc_model = solve_model_at_anchors()

    print(f"\nmw_mix = {mw_mix*1000:.4f} g/mol\n")

    # Fit-anchor T's: identical to the HONEYWELL_*_TK constants EXCEPT the critical
    # entry, which uses the MODEL's own solved Tc (tc_model) instead of Honeywell's
    # stated 382.0389 K. The correction curve gets evaluated later against the
    # model's own T's, so the critical anchor's x-coordinate must be the model's Tc,
    # not Honeywell's -- these differ by ~0.145 K (see solve_model_at_anchors print).
    # The corresponding Honeywell VALUE at that anchor (critical density/pressure) is
    # unchanged -- only which T that value gets plotted/fit against shifts slightly.
    liquid_fit_T = HONEYWELL_LIQUID_ANCHORS_TK[:-1] + [tc_model]
    vapor_fit_T = HONEYWELL_VAPOR_ANCHORS_TK[:-1] + [tc_model]
    psat_fit_T = HONEYWELL_PSAT_ANCHORS_TK[:-1] + [tc_model]

    print("=== Fitting liquid-density residual curve ===")
    rho_liq_correction_fn, rho_liq_residuals = _fit_residual_curve(
        liquid_fit_T, HONEYWELL_LIQUID_ANCHORS_KGM3, rho_liq_model_kgm3)
    for t_k, hw, resid in zip(liquid_fit_T, HONEYWELL_LIQUID_ANCHORS_KGM3, rho_liq_residuals):
        model_val = rho_liq_model_kgm3[t_k]
        print(f"  T={t_k:.4f} K: Honeywell={hw:.4f} kg/m3   Model={model_val:.4f} kg/m3   "
              f"Residual={resid:+.4f} kg/m3 ({100*resid/hw:+.3f}%)")
    print(f"  Fitted quadratic coefficients (highest power first): {rho_liq_correction_fn.coefficients}")

    print("\n=== Fitting vapor-density residual curve ===")
    rho_vap_correction_fn, rho_vap_residuals = _fit_residual_curve(
        vapor_fit_T, HONEYWELL_VAPOR_ANCHORS_KGM3, rho_vap_model_kgm3)
    for t_k, hw, resid in zip(vapor_fit_T, HONEYWELL_VAPOR_ANCHORS_KGM3, rho_vap_residuals):
        model_val = rho_vap_model_kgm3[t_k]
        print(f"  T={t_k:.4f} K: Honeywell={hw:.4f} kg/m3   Model={model_val:.4f} kg/m3   "
              f"Residual={resid:+.4f} kg/m3 ({100*resid/hw:+.3f}%)")
    print(f"  Fitted quadratic coefficients (highest power first): {rho_vap_correction_fn.coefficients}")

    print("\n=== Fitting saturation-pressure residual curve ===")
    p_correction_fn, p_residuals = _fit_residual_curve(
        psat_fit_T, HONEYWELL_PSAT_ANCHORS_PA, p_model_pa)
    for t_k, hw, resid in zip(psat_fit_T, HONEYWELL_PSAT_ANCHORS_PA, p_residuals):
        model_val = p_model_pa[t_k]
        print(f"  T={t_k:.4f} K: Honeywell={hw:.1f} Pa   Model={model_val:.1f} Pa   "
              f"Residual={resid:+.1f} Pa ({100*resid/hw:+.3f}%)")
    print(f"  Fitted quadratic coefficients (highest power first): {p_correction_fn.coefficients}")

    print("\n=== Sanity check: fitted curves must reproduce the known residuals exactly ===")
    for t_k, resid in zip(liquid_fit_T, rho_liq_residuals):
        fitted = rho_liq_correction_fn(t_k)
        ok = "OK" if abs(fitted - resid) < 1e-6 else "MISMATCH"
        print(f"  liquid  T={t_k:.4f} K: fitted={fitted:+.6f}  raw={resid:+.6f}  {ok}")
    for t_k, resid in zip(vapor_fit_T, rho_vap_residuals):
        fitted = rho_vap_correction_fn(t_k)
        ok = "OK" if abs(fitted - resid) < 1e-6 else "MISMATCH"
        print(f"  vapor   T={t_k:.4f} K: fitted={fitted:+.6f}  raw={resid:+.6f}  {ok}")
    for t_k, resid in zip(psat_fit_T, p_residuals):
        fitted = p_correction_fn(t_k)
        ok = "OK" if abs(fitted - resid) < 1e-3 else "MISMATCH"
        print(f"  Psat    T={t_k:.4f} K: fitted={fitted:+.3f}  raw={resid:+.3f}  {ok}")

    print("\n=== Fitting PCHIP (shape-preserving) curves for comparison ===")
    rho_liq_pchip_fn, _ = _fit_pchip_curve(liquid_fit_T, HONEYWELL_LIQUID_ANCHORS_KGM3, rho_liq_model_kgm3)
    rho_vap_pchip_fn, _ = _fit_pchip_curve(vapor_fit_T, HONEYWELL_VAPOR_ANCHORS_KGM3, rho_vap_model_kgm3)
    p_pchip_fn, _ = _fit_pchip_curve(psat_fit_T, HONEYWELL_PSAT_ANCHORS_PA, p_model_pa)
    print("  (fitted -- both curves pass through the same 3 known points exactly by construction;")
    print("   the difference only shows up BETWEEN anchors, checked next)")

    print("\n=== Quadratic vs. PCHIP: checking the overshoot found earlier, between the two low-T anchors ===")
    print("  Vapor-density residual (this is the curve that showed the -2.44 kg/m3 dip):")
    lo_v, hi_v = sorted(vapor_fit_T)[0], sorted(vapor_fit_T)[1]  # 254.26 -> 298.15
    for t in np.linspace(lo_v, hi_v, 10):
        print(f"    T={t:7.2f} K: quadratic={rho_vap_correction_fn(t):+8.4f} kg/m3   PCHIP={rho_vap_pchip_fn(t):+8.4f} kg/m3")

    print("\n  Liquid-density residual (for comparison -- this curve was already well-behaved):")
    lo_l, hi_l = sorted(liquid_fit_T)[0], sorted(liquid_fit_T)[1]  # 273.15 -> 298.15
    for t in np.linspace(lo_l, hi_l, 6):
        print(f"    T={t:7.2f} K: quadratic={rho_liq_correction_fn(t):+8.4f} kg/m3   PCHIP={rho_liq_pchip_fn(t):+8.4f} kg/m3")

    print("\n  Saturation-pressure residual (for comparison -- this curve was already fairly well-behaved):")
    lo_p, hi_p = sorted(psat_fit_T)[0], sorted(psat_fit_T)[1]  # 254.26 -> 298.15
    for t in np.linspace(lo_p, hi_p, 6):
        print(f"    T={t:7.2f} K: quadratic={p_correction_fn(t):+10.2f} Pa   PCHIP={p_pchip_fn(t):+10.2f} Pa")

    print("\n=== Relative-magnitude comparison: is pressure really 'protected' vs. density? ===")
    max_rho_liq_pct = max(abs(100*r/hw) for r, hw in zip(rho_liq_residuals, HONEYWELL_LIQUID_ANCHORS_KGM3))
    max_rho_vap_pct = max(abs(100*r/hw) for r, hw in zip(rho_vap_residuals, HONEYWELL_VAPOR_ANCHORS_KGM3))
    max_p_pct = max(abs(100*r/hw) for r, hw in zip(p_residuals, HONEYWELL_PSAT_ANCHORS_PA))
    print(f"  Max |residual| as %% of Honeywell value: liquid density = {max_rho_liq_pct:.3f}%%, "
          f"vapor density = {max_rho_vap_pct:.3f}%%, pressure = {max_p_pct:.3f}%%")


if __name__ == "__main__":
    main()
