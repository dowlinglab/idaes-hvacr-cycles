"""
Plot the dew line's own saturation entropy s_v(T) to SHOW the "hump, not a
ramp" shape that is the root cause behind the vapor-side isentropes >0.37
Btu/lb-R rendering wrong (see diagnose_vapor_isentrope_anchor_collision.py
and PROJECT_CONTEXT.md, 2026-08-14 entry).

Purpose
-------
Step 2 of that diagnostic claimed (in words) that s_v(T) along the dew
(saturated vapor) branch rises from the cold end, peaks at an INTERIOR
temperature (~153.4F), then falls back down approaching the critical point --
a hump shape, not a monotonic ramp. This script draws that curve directly so
it can be seen, not just read as printed numbers.

What it does
------------
1. Runs the same true-VLE dome sweep used everywhere else in this repo
   (run_true_vle_envelope) over the same T range/resolution as the normal
   CLI (255.4-380.0 K, n=80).
2. For every CONVERGED dew-branch row, computes s_v(T) directly via
   _mix_entropy_direct (same function the model itself uses -- no
   correction, no fitting, just the raw EOS-based entropy).
3. Plots s_v vs. T (bottom x-axis, K) and vs. T (top x-axis, deg F) on the
   same figure, marks the peak explicitly, and overlays the 15 standard
   isentrope target entropies (ISENTROPE_VALUES_BTU_LBMR) as horizontal
   reference lines so it's visually obvious which targets exceed the peak
   (and therefore collide onto the same anchor row).
4. Saves the figure to dew_entropy_shape.png in the current directory.

Usage
-----
Run from R515B_props_validated/ (same directory as
mixture_isentrope_validation.py and its own copy of linear_model_codex.py):
    python3 plot_dew_entropy_shape.py

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

converged_dew = [r for r in dew_rows if r["status"] == "CONVERGED"]
rows = sorted(
    (r["T_K"], _mix_entropy_direct(d1, d2, r["T_K"], r["rho_v_molm3"], z1))
    for r in converged_dew
)
t_k = np.array([r[0] for r in rows])
s_v = np.array([r[1] for r in rows])
t_f = (t_k - 273.15) * 9.0 / 5.0 + 32.0

peak_idx = int(np.argmax(s_v))
print(f"Peak: T={t_k[peak_idx]:.2f}K ({t_f[peak_idx]:.1f}F), s_v={s_v[peak_idx]:.3f} J/(mol*K)")

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(t_k, s_v, "o-", color="tab:blue", markersize=3, linewidth=1.5,
         label="Dew-branch saturation entropy $s_v(T)$ (raw EOS, no correction)")
ax.axvline(t_k[peak_idx], color="tab:red", linestyle=":", linewidth=1)
ax.plot(t_k[peak_idx], s_v[peak_idx], "*", color="tab:red", markersize=18, zorder=5)
ax.annotate(
    f"PEAK\nT={t_k[peak_idx]:.1f}K ({t_f[peak_idx]:.1f}°F)\ns_v={s_v[peak_idx]:.2f} J/(mol·K)",
    xy=(t_k[peak_idx], s_v[peak_idx]),
    xytext=(15, -35), textcoords="offset points",
    fontsize=9, color="tab:red",
    arrowprops=dict(arrowstyle="->", color="tab:red"),
)

# Overlay the 15 standard isentrope target entropies as horizontal lines,
# colored by whether they exceed the peak (collide) or not (well-matched).
for s_btu in ISENTROPE_VALUES_BTU_LBMR:
    s_target = s_btu * BTU_LBMR_TO_JKGK * mw_mix
    exceeds = s_target > s_v[peak_idx]
    color = "tab:red" if exceeds else "0.75"
    lw = 1.0 if exceeds else 0.6
    ax.axhline(s_target, color=color, linestyle="--", linewidth=lw, alpha=0.6)
    ax.text(t_k[-1] + 0.5, s_target, f"{s_btu:.2f}", fontsize=7,
            color=color, va="center", ha="left")

ax.set_xlabel("Temperature, T (K)")
ax.set_ylabel("Saturation entropy, $s_v$ (J/(mol·K))")
ax.set_title(
    "Dew-branch entropy is a HUMP, not a ramp\n"
    "(red dashed lines = isentrope targets that exceed the peak → all collide onto the same anchor row)"
)
ax.legend(loc="lower left", fontsize=8)
ax.grid(True, alpha=0.3)

# Secondary top axis in deg F for readability against the Honeywell chart
ax2 = ax.secondary_xaxis("top", functions=(lambda t: (t - 273.15) * 9.0 / 5.0 + 32.0,
                                             lambda f: (f - 32.0) * 5.0 / 9.0 + 273.15))
ax2.set_xlabel("Temperature, T (°F)")

fig.tight_layout()
out_path = "dew_entropy_shape.png"
fig.savefig(out_path, dpi=150)
print(f"Saved: {out_path}")
