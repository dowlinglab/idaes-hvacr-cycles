"""
Check the achievable entropy range along the dew (saturated vapor) branch,
and compare it directly against the mixture's critical-point entropy.

Purpose
-------
diagnose_isentrope_vapor_037.py reported ZERO anchor candidates for the
0.37 Btu/lb-R vapor-side isentrope -- every one of the 80 converged dew
rows across T=255.4-380 K had s_v >= s_target. That's consistent with two
different explanations:

  (a) The mixture's critical entropy s_c is itself above 0.37 Btu/lb-R,
      so NO point on the sampled dew branch (whose entropy is bounded
      below by s_c, if s_v(T) decreases monotonically toward Tc) can ever
      dip below 0.37 -- this isentrope simply doesn't intersect the dew
      branch in the physical sense, and 0.37 not appearing is CORRECT,
      not a bug.

  (b) Something narrower/T-range-specific -- e.g. the sampled Tmax=380 K
      is still ~2 K short of the actual Tc (~382.04 K per earlier
      breadcrumb notes), so the branch's lowest achievable entropy at the
      sampled points hasn't yet dropped low enough, even though the TRUE
      critical entropy is below 0.37.

This script reports the min/max s_v actually achieved across the sampled
dew branch, plus the entropy at the solved critical point itself (via the
same solve_mixture_critical_point() call already used to close the dome),
so we can tell which of (a)/(b) is happening from real numbers instead of
guessing.

Usage
-----
Run from R515B_props_validated/ (same directory as
mixture_isentrope_validation.py and its own copy of linear_model_codex.py):
    python3 check_dew_entropy_range.py

Read-only diagnostic -- does not modify anything.
"""

from mixture_isentrope_validation import (
    run_true_vle_envelope,
    load_idaes_helmholtz_json,
    mw_from_json,
    _mix_entropy_direct,
    BTU_LBMR_TO_JKGK,
)
import numpy as np

FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
TMIN, TMAX, N = 255.4, 380.0, 80
S_TARGETS_BTU = [0.35, 0.37]  # compare the "checks out (slightly low)" one against the "missing" one

t_vals = np.linspace(TMIN, TMAX, N)
print(f"Running dome sweep ({FLUID1}/{FLUID2}, w1={W1}, T={TMIN}-{TMAX} K, n={N})...")
bubble_rows, dew_rows, z1, crit_point = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals)

d1 = load_idaes_helmholtz_json(FLUID1)
d2 = load_idaes_helmholtz_json(FLUID2)
mw1 = mw_from_json(d1)
mw2 = mw_from_json(d2)
mw_mix = z1 * mw1 + (1.0 - z1) * mw2

converged_dew = [r for r in dew_rows if r["status"] == "CONVERGED"]
s_v_list = [(r["T_K"], _mix_entropy_direct(d1, d2, r["T_K"], r["rho_v_molm3"], z1)) for r in converged_dew]

print(f"\n{len(converged_dew)}/{len(dew_rows)} dew rows converged.\n")
print("=== Dew-branch entropy range (sampled T=255.4-380 K) ===")
t_lo, s_lo = min(s_v_list, key=lambda p: p[0])
t_hi, s_hi = max(s_v_list, key=lambda p: p[0])
s_min_pt = min(s_v_list, key=lambda p: p[1])
s_max_pt = max(s_v_list, key=lambda p: p[1])
print(f"At T={t_lo:.2f} K (coldest sampled): s_v = {s_v_list[0][1]:.4f} J/(mol*K) = {s_v_list[0][1]/(BTU_LBMR_TO_JKGK*mw_mix):.4f} Btu/lb-R")
print(f"At T={t_hi:.2f} K (warmest sampled): s_v = {s_v_list[-1][1]:.4f} J/(mol*K) = {s_v_list[-1][1]/(BTU_LBMR_TO_JKGK*mw_mix):.4f} Btu/lb-R")
print(f"MIN s_v over whole branch: {s_min_pt[1]:.4f} J/(mol*K) = {s_min_pt[1]/(BTU_LBMR_TO_JKGK*mw_mix):.4f} Btu/lb-R  at T={s_min_pt[0]:.2f} K")
print(f"MAX s_v over whole branch: {s_max_pt[1]:.4f} J/(mol*K) = {s_max_pt[1]/(BTU_LBMR_TO_JKGK*mw_mix):.4f} Btu/lb-R  at T={s_max_pt[0]:.2f} K")

print("\n=== Critical point ===")
if crit_point is not None and crit_point.get("converged", False):
    tc = crit_point["T_K"]
    rhoc = crit_point["rho_molm3"]
    x1c = crit_point.get("x1", z1)
    s_c = _mix_entropy_direct(d1, d2, tc, rhoc, x1c)
    s_c_btu = s_c / (BTU_LBMR_TO_JKGK * mw_mix)
    print(f"Tc = {tc:.4f} K   Pc = {crit_point['P_Pa']:.1f} Pa   rho_c = {rhoc:.4f} mol/m3")
    print(f"s_c (entropy AT critical point) = {s_c:.4f} J/(mol*K) = {s_c_btu:.4f} Btu/lb-R")
    print(f"(sampled Tmax={TMAX} K is {tc-TMAX:.3f} K below Tc -- gap not yet sampled by the sweep)")
else:
    print("Critical point solve did not converge or was not computed -- no s_c available.")
    s_c_btu = None

print("\n=== Comparison against isentrope targets ===")
for s_btu in S_TARGETS_BTU:
    s_target = s_btu * BTU_LBMR_TO_JKGK * mw_mix
    print(f"{s_btu} Btu/lb-R = {s_target:.4f} J/(mol*K)  "
          f"{'BELOW' if s_target < s_min_pt[1] else 'within/above'} sampled dew-branch minimum ({s_min_pt[1]:.4f})"
          + (f"; {'BELOW' if s_c_btu is not None and s_btu < s_c_btu else 'ABOVE/AT'} s_c ({s_c_btu:.4f})" if s_c_btu is not None else ""))
