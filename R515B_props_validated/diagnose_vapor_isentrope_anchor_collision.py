"""
Diagnose why the vapor-side isentropes above 0.37 Btu/lb-R (0.39, 0.41, 0.43,
0.45, 0.47, 0.49) render wrong on the Honeywell-units plot.

Purpose
-------
User reported (2026-08-14, via side-by-side screenshot comparison against the
real Honeywell Solstice N15 p-h chart) that all vapor-side isentropes above
0.37 Btu/lb-R look wrong. This script walks through the exact diagnostic
steps used to find the root cause, in order, so the result can be
independently reproduced and checked.

Diagnostic steps (run in this order, matching how the bug was actually found)
-------------------------------------------------------------------------
1. Confirm all 15 standard isentropes (ISENTROPE_VALUES_BTU_LBMR) fully
   converge at the CLI's normal resolution (n=80) -- ruling out a solver
   crash as the cause (they all do converge; the lines are drawn, just in
   the wrong place/shape).
2. Print the dew branch's own saturation entropy s_v(T) across the whole
   sampled temperature range, T-sorted, to see its actual shape.
3. For each vapor-side isentrope target (0.37 through 0.49 Btu/lb-R),
   reproduce the SAME anchor-selection rule used inside
   compute_isentrope_vapor_side() (mixture_isentrope_validation.py, "pick
   the converged dew row with s_v(T) closest to, but still below, the
   target"), and print which (T, s_v, P) row gets picked for each target,
   plus the entropy gap between that anchor and the target.

What this reveals
------------------
The dew branch's s_v(T) is NOT monotonic across the sampled range -- it
rises from the cold end, peaks at an INTERIOR temperature, then falls back
down approaching the critical point (a real, physical "hump" shape, not a
numerical artifact -- see printed values below). Its peak value is a hard
ceiling on what the "closest row below target" anchor rule can ever find.
Every isentrope target above that peak (0.41 through 0.49 Btu/lb-R, since
0.39's target happens to sit just below the peak) ends up picking the
EXACT SAME anchor row (the peak itself), regardless of how far above the
peak the target actually is. compute_isentrope_vapor_side then asks a
single 2-unknown Newton solve (scipy.optimize.root, "hybr") to jump
directly from that one shared seed to the target entropy, in ONE step,
before any of the normal incremental pressure-walk continuation even
starts. The size of that first jump grows from 5.5 J/(mol*K) (0.41) up to
44.8 J/(mol*K) (0.49) -- a large, increasingly extreme ask from an
identical starting point, which is the likely reason these five isentropes
end up in visually wrong/implausible positions rather than fanning out
smoothly like the lower-s ones.

Usage
-----
Run from R515B_props_validated/ (same directory as
mixture_isentrope_validation.py and its own copy of linear_model_codex.py):
    python3 diagnose_vapor_isentrope_anchor_collision.py

Read-only diagnostic -- does not modify anything.
"""

import numpy as np

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    run_true_vle_envelope,
    compute_isentrope_liquid_side,
    compute_isentrope_vapor_side,
    _mix_entropy_direct,
    BTU_LBMR_TO_JKGK,
    ISENTROPE_VALUES_BTU_LBMR,
)

FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
TMIN, TMAX, N = 255.4, 380.0, 80

d1 = load_idaes_helmholtz_json(FLUID1)
d2 = load_idaes_helmholtz_json(FLUID2)
mw1 = mw_from_json(d1)
mw2 = mw_from_json(d2)
z1 = w1_to_x1(W1, mw1, mw2)
mw_mix = z1 * mw1 + (1.0 - z1) * mw2

t_vals = np.linspace(TMIN, TMAX, N)
print(f"Running dome sweep ({FLUID1}/{FLUID2}, w1={W1}, T={TMIN}-{TMAX}K, n={N})...")
bubble_rows, dew_rows, z1_out, crit_point = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals)

# ---------------------------------------------------------------------
# STEP 1: confirm all 15 standard isentropes fully converge at n=80
# (rules out "solver crash" -- the lines ARE being drawn, just wrong)
# ---------------------------------------------------------------------
s_values_jmolK = [s * BTU_LBMR_TO_JKGK * mw_mix for s in ISENTROPE_VALUES_BTU_LBMR]
liq_ext = compute_isentrope_liquid_side(d1, d2, z1, bubble_rows, s_values_jmolK, crit_point=crit_point)
vap_ext = compute_isentrope_vapor_side(d1, d2, z1, dew_rows, s_values_jmolK, crit_point=crit_point)

print("\n=== STEP 1: point counts per isentrope (out of up to 40 each) ===")
print(f"{'s(Btu/lb-R)':>12} {'liq pts':>8} {'vap pts':>8}")
for s_btu, s_target in zip(ISENTROPE_VALUES_BTU_LBMR, s_values_jmolK):
    print(f"{s_btu:>12.2f} {len(liq_ext.get(s_target, [])):>8} {len(vap_ext.get(s_target, [])):>8}")
print("(all vapor-side counts above should read 40 -- confirms this is NOT a convergence failure)")

# ---------------------------------------------------------------------
# STEP 2: print the dew branch's own saturation entropy s_v(T) shape
# ---------------------------------------------------------------------
converged_dew = [r for r in dew_rows if r["status"] == "CONVERGED"]
s_v_list = sorted(
    (r["T_K"], _mix_entropy_direct(d1, d2, r["T_K"], r["rho_v_molm3"], z1), r["P_Pa"])
    for r in converged_dew
)
print(f"\n=== STEP 2: dew-branch s_v(T) shape ({len(s_v_list)} converged rows) ===")
print("First 3 rows (cold end):")
for t, s, p in s_v_list[:3]:
    print(f"  T={t:.2f}K ({(t-273.15)*9/5+32:.1f}F): s_v={s:.3f} J/(mol*K), P={p/6894.757:.1f} psia")
peak = max(s_v_list, key=lambda row: row[1])
print(f"PEAK (interior T, not at either end): T={peak[0]:.2f}K ({(peak[0]-273.15)*9/5+32:.1f}F), "
      f"s_v={peak[1]:.3f} J/(mol*K), P={peak[2]/6894.757:.1f} psia")
print("Last 3 rows (warm end, approaching Tc):")
for t, s, p in s_v_list[-3:]:
    print(f"  T={t:.2f}K ({(t-273.15)*9/5+32:.1f}F): s_v={s:.3f} J/(mol*K), P={p/6894.757:.1f} psia")
print("(s_v(T) rises from the cold end, peaks at an INTERIOR temperature, then falls -- non-monotonic)")

# ---------------------------------------------------------------------
# STEP 3: reproduce the anchor-selection rule for each target and show
# which row gets picked + the entropy gap it must bridge in one Newton step
# ---------------------------------------------------------------------
print("\n=== STEP 3: anchor row picked per isentrope target (same rule as compute_isentrope_vapor_side) ===")
picked_anchors = {}
for s_btu in [0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49]:
    s_target = s_btu * BTU_LBMR_TO_JKGK * mw_mix
    best_row, best_gap = None, None
    for t, s_v, p in s_v_list:
        if s_v >= s_target:
            continue  # not reachable via expansion from this row
        gap = s_target - s_v
        if best_gap is None or gap < best_gap:
            best_gap = gap
            best_row = (t, s_v, p)
    if best_row is None:
        print(f"  s={s_btu}: target={s_target:.2f} J/(mol*K) -- NO ordinary row qualifies (falls back to crit_point or is unreachable)")
        continue
    t, s_v, p = best_row
    picked_anchors[s_btu] = (t, s_v, p)
    print(f"  s={s_btu}: target={s_target:>7.2f} J/(mol*K)  ->  anchor T={t:.2f}K ({(t-273.15)*9/5+32:>6.1f}F), "
          f"anchor P={p/6894.757:>6.1f} psia, GAP TO BRIDGE IN ONE NEWTON STEP = {best_gap:>6.2f} J/(mol*K)")

# Highlight the collision explicitly
anchor_temps = {s: row[0] for s, row in picked_anchors.items() if s >= 0.41}
if len(set(anchor_temps.values())) == 1 and len(anchor_temps) > 1:
    shared_t = next(iter(anchor_temps.values()))
    print(f"\n*** CONFIRMED: isentropes {sorted(anchor_temps.keys())} ALL share the exact same anchor "
          f"(T={shared_t:.2f}K) despite having very different targets. ***")
    print("Each is asked to jump a different, increasingly large entropy gap from that SAME starting")
    print("point in a single 2-unknown Newton solve, before any gradual pressure-walk continuation")
    print("begins -- the likely cause of the wrong/implausible shapes seen on the plot.")
