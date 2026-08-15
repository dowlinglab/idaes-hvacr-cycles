"""
Plot the bubble line's own saturation entropy s_l(T) -- the liquid-side
counterpart to plot_dew_entropy_shape.py -- to check whether the bubble
branch has the same "hump, not a ramp" non-monotonic shape that was found
on the dew (vapor) branch, or whether it behaves differently.

Purpose
-------
plot_dew_entropy_shape.py showed the dew (saturated vapor) branch's own
entropy s_v(T) rises from the cold end, peaks at an interior temperature
(~153.4F), then falls back down approaching the critical point. This script
does the exact same computation and plot, but for the bubble (saturated
liquid) branch's entropy s_l(T), so the two can be compared directly.

What it does
------------
1. Runs the same true-VLE dome sweep used everywhere else in this repo
   (run_true_vle_envelope) over the same T range/resolution as the normal
   CLI (255.4-380.0 K, n=80).
2. For every CONVERGED bubble-branch row, computes s_l(T) directly via
   _mix_entropy_direct (same function the model itself uses -- no
   correction, no fitting, just the raw EOS-based entropy).
3. Plots s_l vs. T (bottom x-axis, K) and vs. T (top x-axis, deg F) on the
   same figure, marks both the min and max explicitly (since we don't know
   in advance whether this branch is monotonic or not), and overlays the 15
   standard isentrope target entropies (ISENTROPE_VALUES_BTU_LBMR) as
   horizontal reference lines, same as the dew-branch plot, for direct
   comparison.
4. Saves the figure to bubble_entropy_shape.png in the current directory.

Usage
-----
Run from R515B_props_validated/ (same directory as
mixture_isentrope_validation.py and its own copy of linear_model_codex.py):
    python3 plot_bubble_entropy_shape.py

Read-only -- does not modify the model or any other file.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    run_true_vle_envelope,
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

converged_bubble = [r for r in bubble_rows if r["status"] == "CONVERGED"]
rows = sorted(
    (r["T_K"], _mix_entropy_direct(d1, d2, r["T_K"], r["rho_l_molm3"], z1))
    for r in converged_bubble
)
t_k = np.array([r[0] for r in rows])
s_l = np.array([r[1] for r in rows])
t_f = (t_k - 273.15) * 9.0 / 5.0 + 32.0

print(f"Converged bubble rows: {len(rows)}")
print("First 3 rows (cold end):")
for t, s in rows[:3]:
    print(f"  T={t:.2f}K ({(t-273.15)*9/5+32:.1f}F): s_l={s:.3f} J/(mol*K)")
print("Last 3 rows (warm end, approaching Tc):")
for t, s in rows[-3:]:
    print(f"  T={t:.2f}K ({(t-273.15)*9/5+32:.1f}F): s_l={s:.3f} J/(mol*K)")

max_idx = int(np.argmax(s_l))
min_idx = int(np.argmin(s_l))
print(f"MAX: T={t_k[max_idx]:.2f}K ({t_f[max_idx]:.1f}F), s_l={s_l[max_idx]:.3f} J/(mol*K)")
print(f"MIN: T={t_k[min_idx]:.2f}K ({t_f[min_idx]:.1f}F), s_l={s_l[min_idx]:.3f} J/(mol*K)")

is_monotonic_inc = np.all(np.diff(s_l) >= -1e-9)
is_monotonic_dec = np.all(np.diff(s_l) <= 1e-9)
if is_monotonic_inc or is_monotonic_dec:
    print("Bubble branch s_l(T) is MONOTONIC across the sampled range (a ramp, not a hump).")
else:
    print("Bubble branch s_l(T) is NON-MONOTONIC across the sampled range (has an interior peak/valley).")

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(t_k, s_l, "o-", color="tab:green", markersize=3, linewidth=1.5,
         label="Bubble-branch saturation entropy $s_l(T)$ (raw EOS, no correction)")

# Mark max and min explicitly (whichever is the "interesting" extremum
# depends on what the curve actually looks like -- show both).
ax.plot(t_k[max_idx], s_l[max_idx], "*", color="tab:red", markersize=16, zorder=5)
ax.annotate(
    f"MAX\nT={t_k[max_idx]:.1f}K ({t_f[max_idx]:.1f}°F)\ns_l={s_l[max_idx]:.2f} J/(mol·K)",
    xy=(t_k[max_idx], s_l[max_idx]),
    xytext=(15, 10), textcoords="offset points",
    fontsize=9, color="tab:red",
    arrowprops=dict(arrowstyle="->", color="tab:red"),
)
ax.plot(t_k[min_idx], s_l[min_idx], "*", color="tab:purple", markersize=16, zorder=5)
ax.annotate(
    f"MIN\nT={t_k[min_idx]:.1f}K ({t_f[min_idx]:.1f}°F)\ns_l={s_l[min_idx]:.2f} J/(mol·K)",
    xy=(t_k[min_idx], s_l[min_idx]),
    xytext=(15, -35), textcoords="offset points",
    fontsize=9, color="tab:purple",
    arrowprops=dict(arrowstyle="->", color="tab:purple"),
)

# Overlay the 15 standard isentrope target entropies as horizontal lines,
# colored by whether they exceed the branch's max (would collide) or not.
for s_btu in ISENTROPE_VALUES_BTU_LBMR:
    s_target = s_btu * BTU_LBMR_TO_JKGK * mw_mix
    exceeds = s_target > s_l[max_idx] or s_target < s_l[min_idx]
    color = "tab:red" if exceeds else "0.75"
    lw = 1.0 if exceeds else 0.6
    ax.axhline(s_target, color=color, linestyle="--", linewidth=lw, alpha=0.6)
    ax.text(t_k[-1] + 0.5, s_target, f"{s_btu:.2f}", fontsize=7,
            color=color, va="center", ha="left")

ax.set_xlabel("Temperature, T (K)")
ax.set_ylabel("Saturation entropy, $s_l$ (J/(mol·K))")
ax.set_title(
    "Bubble-branch entropy shape (liquid side)\n"
    "compare against dew_entropy_shape.png (vapor side) for the hump-vs-ramp check"
)
ax.legend(loc="best", fontsize=8)
ax.grid(True, alpha=0.3)

ax2 = ax.secondary_xaxis("top", functions=(lambda t: (t - 273.15) * 9.0 / 5.0 + 32.0,
                                             lambda f: (f - 32.0) * 5.0 / 9.0 + 273.15))
ax2.set_xlabel("Temperature, T (°F)")

fig.tight_layout()
out_path = "bubble_entropy_shape.png"
fig.savefig(out_path, dpi=150)
print(f"Saved: {out_path}")
