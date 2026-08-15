"""
Diagnose why the liquid-side isentrope at s=0.37 Btu/lb-R appears NOT to
begin at the critical point on the rendered plot, even though its target
entropy (182.00 J/(mol*K)) is close to the critical entropy and its
crit-point-fallback anchor lands very close to Tc in temperature and pressure.

Purpose
-------
User asked (2026-08-14): "0.37 should begin at critical and it doesn't."
This script traces through exactly why, using three checks:

1. Confirms compute_isentrope_liquid_side()'s crit-point-fallback anchor for
   s=0.37 (used because 182.00 J/(mol*K) exceeds the T-grid's own sampled
   bubble-branch maximum, ~179.4 J/(mol*K) at T=380K, per
   diagnose_vapor_isentrope_anchor_collision.py's earlier finding) against
   an independently-computed, much finer T-grid pushed to within ~0.04K of
   the true Tc -- to check whether the fallback's (T, P, h) is actually a
   good approximation of the TRUE bubble-line point at that entropy, or if
   it's the fallback itself that's inaccurate.
2. Compares that anchor point's own enthalpy against the critical point's
   own enthalpy directly -- this is the number that actually determines
   whether the line LOOKS like it starts at the dome's tip on a p-h plot
   (proximity in P and T alone is not enough; h is the plotted x-axis).
3. Flags a secondary, unrelated-but-adjacent finding: the critical point
   solve itself (solve_mixture_critical_point, invoked internally by
   run_true_vle_envelope) is sensitive to the T-grid it's given -- two
   different T-grids here produce two different T_c/s_c values, a
   discrepancy worth tracking separately.

What this reveals
------------------
The crit-point-fallback anchor IS numerically accurate (matches the
independently-interpolated true bubble-line point to ~0.02% in h) -- so the
fallback logic itself is not the bug. But that true point's own enthalpy is
still ~1300+ J/mol (~2.9%) below the critical point's enthalpy, despite
being within ~0.5% of critical in pressure and ~0.1% in temperature. This is
consistent with real near-critical divergence (dh/ds along the saturation
curve gets very steep approaching Tc) -- NOT obviously a computational bug.

Separately, mixture_isentrope_validation.py's compute_isentropes_two_phase()
-- the function that would draw each isentrope's actual path through the
INTERIOR of the two-phase dome (which is what would let a sub-critical-
entropy isentrope emerge from inside the dome rather than only ever
asymptotically approaching its edge from outside) -- is currently DISABLED
(see the "DISABLED 2026-08-14" comment right before the _cli() plotting
call). That is a separate, deliberate prior design decision (isentropes
believed to only belong in single-phase regions per the Honeywell reference
chart), not a bug in this script's checks -- but it is the reason 0.37 (and
any other low-entropy, near-critical isentrope) can only ever show a
liquid-side segment that approaches, but does not visually touch, the dome
tip. Whether that design decision is itself correct for isentropes this
close to critical is an OPEN QUESTION, flagged for the user -- not resolved
by this script.

Usage
-----
Run from R515B_props_validated/:
    python3 diagnose_037_critical_approach.py

Read-only -- does not modify anything.
"""

import numpy as np

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    run_true_vle_envelope,
    _mix_entropy_direct,
    mix_state,
    compute_isentrope_liquid_side,
    BTU_LBMR_TO_JKGK,
)

FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
S_TARGET_BTU = 0.37

d1 = load_idaes_helmholtz_json(FLUID1)
d2 = load_idaes_helmholtz_json(FLUID2)
mw1 = mw_from_json(d1)
mw2 = mw_from_json(d2)
z1 = w1_to_x1(W1, mw1, mw2)
mw_mix = z1 * mw1 + (1.0 - z1) * mw2
s_target = S_TARGET_BTU * BTU_LBMR_TO_JKGK * mw_mix
print(f"Target entropy for {S_TARGET_BTU} Btu/lb-R: {s_target:.3f} J/(mol*K)")

# ---------------------------------------------------------------------
# STEP 1: reproduce the actual fallback anchor via the normal CLI grid
# ---------------------------------------------------------------------
TMIN, TMAX, N = 255.4, 380.0, 80
t_vals_normal = np.linspace(TMIN, TMAX, N)
print(f"\n=== STEP 1: normal CLI grid (T={TMIN}-{TMAX}K, n={N}) ===")
bubble_rows, dew_rows, _, crit_point_normal = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals_normal)
liq_ext = compute_isentrope_liquid_side(d1, d2, z1, bubble_rows, [s_target], crit_point=crit_point_normal)
pts = liq_ext.get(s_target, [])
print(f"Liquid-side extension for s={S_TARGET_BTU}: {len(pts)} points")
first = pts[0]
t_c_n, rho_c_n = crit_point_normal["T_K"], crit_point_normal["rho_molm3"]
h_c_n = mix_state(d1, d2, t_c_n, rho_c_n, z1).h_jmol
print(f"Fallback first point: T={first['T_K']:.4f}K, P={first['P_Pa']/1e6:.4f}MPa, h={first['h_Jmol']:.2f} J/mol")
print(f"Critical point (this grid): T={t_c_n:.4f}K, P={crit_point_normal['P_Pa']/1e6:.4f}MPa, h={h_c_n:.2f} J/mol")
print(f"h gap (fallback point vs critical point): {h_c_n - first['h_Jmol']:.2f} J/mol "
      f"({100*(h_c_n - first['h_Jmol'])/h_c_n:.2f}% of h_c)")

# ---------------------------------------------------------------------
# STEP 2: independent, much finer T-grid pushed close to Tc, to check
# whether the fallback point matches the TRUE (interpolated) bubble line
# ---------------------------------------------------------------------
print("\n=== STEP 2: independent fine grid (T=370-381.85K, n=60) for cross-check ===")
t_vals_fine = np.linspace(370.0, 381.85, 60)
bubble_rows_fine, dew_rows_fine, _, crit_point_fine = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals_fine)
converged_fine = [r for r in bubble_rows_fine if r["status"] == "CONVERGED"]
rows = sorted(
    (r["T_K"], _mix_entropy_direct(d1, d2, r["T_K"], r["rho_l_molm3"], z1), r["P_Pa"], r["rho_l_molm3"])
    for r in converged_fine
)
ts = np.array([r[0] for r in rows])
ss = np.array([r[1] for r in rows])
ps = np.array([r[2] for r in rows])
rhos = np.array([r[3] for r in rows])
t_c_f, rho_c_f = crit_point_fine["T_K"], crit_point_fine["rho_molm3"]
h_c_f = mix_state(d1, d2, t_c_f, rho_c_f, z1).h_jmol
s_c_f = _mix_entropy_direct(d1, d2, t_c_f, rho_c_f, z1)
print(f"Critical point (fine grid): T={t_c_f:.4f}K, s_c={s_c_f:.4f} J/(mol*K), "
      f"P={crit_point_fine['P_Pa']/1e6:.4f}MPa, h={h_c_f:.2f} J/mol")

if ss[0] <= s_target <= ss[-1]:
    t_i = np.interp(s_target, ss, ts)
    p_i = np.interp(s_target, ss, ps)
    rho_i = np.interp(s_target, ss, rhos)
    h_i = mix_state(d1, d2, t_i, rho_i, z1).h_jmol
    print(f"Interpolated TRUE bubble-line point at s={s_target:.2f}: "
          f"T={t_i:.4f}K (Tc-T={t_c_f-t_i:.4f}K), P={p_i/1e6:.4f}MPa, h={h_i:.2f} J/mol")
    print(f"Compare to Step 1's fallback point: T={first['T_K']:.4f}K, "
          f"P={first['P_Pa']/1e6:.4f}MPa, h={first['h_Jmol']:.2f} J/mol")
    print(f"h difference (fallback vs true-interpolated): "
          f"{first['h_Jmol'] - h_i:.2f} J/mol ({100*abs(first['h_Jmol']-h_i)/h_i:.3f}% -- "
          f"small means the fallback IS a good approximation of the real bubble-line point)")
    print(f"h gap (true point vs critical point): {h_c_f - h_i:.2f} J/mol "
          f"({100*(h_c_f - h_i)/h_c_f:.2f}% of h_c -- this is the REAL, physical gap, "
          f"not a fallback artifact)")
else:
    print(f"s_target={s_target:.2f} outside fine-grid bubble range [{ss[0]:.2f}, {ss[-1]:.2f}]")

# ---------------------------------------------------------------------
# STEP 3: flag the critical-point-solve grid-sensitivity as a separate note
# ---------------------------------------------------------------------
print("\n=== STEP 3: critical point solve sensitivity to input T-grid (separate finding) ===")
print(f"Normal grid  (T up to {TMAX}K): T_c={t_c_n:.4f}K, P_c={crit_point_normal['P_Pa']/1e6:.4f}MPa")
print(f"Fine grid    (T up to 381.85K): T_c={t_c_f:.4f}K, P_c={crit_point_fine['P_Pa']/1e6:.4f}MPa")
print(f"Difference: {abs(t_c_n - t_c_f):.4f}K in T_c -- worth tracking separately, "
      f"not investigated further here.")
