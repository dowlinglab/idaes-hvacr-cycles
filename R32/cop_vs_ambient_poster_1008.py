"""
cop_vs_ambient_poster_1008.py -- the R-32 ambient-vs-COP sweep behind the
poster table, drawn as a figure.

SOURCE OF THE NUMBERS
    phase6_combined_report.csv, written by phase6_final_0924.py on
    2026-09-28 19:07. All four poster rows are in that file to six decimals:

        model               10      15      20      25
        Helmholtz         3.95    3.53    3.19    2.91
        PR / NIST         3.57    3.40    3.13    2.80
        PR / GCGP         3.51    3.26    3.05    2.89
        PR / SPGP         3.74    3.55    3.18    3.04

    NOTE: the poster's "SPGP" column is fp_single_0919 -- the Colon group's
    09/19 delivery, displayed under the SPGP label for slide continuity. Same
    relabelling convention as compare_cp_methods_0903data_SPGPlabel.py. It is
    NOT the 09/03 SPGP fit, which had a very different Tc.

CYCLE CONDITIONS (R-32)
    T_evap = -29 C, T_cond = T_amb + 9, superheat 0, subcool 0,
    compressor isentropic efficiency 0.9999, max pressure ratio 10.

This script READS the CSV rather than hard-coding the values, so re-running
phase6 and re-running this keeps the figure honest. If the CSV is missing it
falls back to the published values above and says so loudly.

    python cop_vs_ambient_poster_1008.py

Writes cop_vs_ambient_poster_1008.png.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = "phase6_combined_report.csv"

# (csv column prefix, poster label, colour, marker, linestyle)
SERIES = [
    ("Helmholtz",      "Helmholtz EoS / CoolProp", "k",          "o", "-"),
    ("NIST",           "PR EoS / NIST",            "tab:blue",   "s", "--"),
    ("GCGP",           "PR EoS / GCGP",            "tab:green",  "^", "--"),
    ("fp_single_0919", "PR EoS / SPGP",            "tab:red",    "D", "--"),
]

# Published fallback -- the poster's own numbers, used only if the CSV is gone.
FALLBACK = {
    "Helmholtz":      [3.948947, 3.526194, 3.190784, 2.914391],
    "NIST":           [3.567653, 3.397422, 3.131608, 2.800516],
    "GCGP":           [3.509379, 3.258122, 3.054898, 2.890954],
    "fp_single_0919": [3.739092, 3.549206, 3.184143, 3.041566],
}
AMB_FALLBACK = [10, 15, 20, 25]


def load():
    """Return (ambients, {prefix: [cop...]}, source_string)."""
    if not os.path.exists(CSV):
        print(f"!! {CSV} not found -- using the published poster values.")
        print("!! Re-run phase6_final_0924.py to regenerate it.")
        return np.array(AMB_FALLBACK, float), dict(FALLBACK), "published values"

    import csv as _csv
    amb, cols = [], {k: [] for k, *_ in SERIES}
    with open(CSV, newline="") as fh:
        for row in _csv.DictReader(fh):
            amb.append(float(row["ambient_C"]))
            for key, *_ in SERIES:
                v = row.get(f"{key}_cop", "")
                conv = row.get(f"{key}_converged", "").strip().lower() == "true"
                cols[key].append(float(v) if (v and conv) else np.nan)
    return np.array(amb), cols, CSV


amb, cop, src = load()
print(f"source: {src}\n")
print(f"{'T_amb [C]':>10}" + "".join(f"{lbl:>28}" for _, lbl, *_ in SERIES))
for i, a in enumerate(amb):
    line = f"{a:>10.0f}"
    for key, *_ in SERIES:
        v = cop[key][i]
        line += f"{v:>28.4f}" if np.isfinite(v) else f"{'--':>28}"
    print(line)

# Poster panel proportions -- same 6.4 x 5.0 as the T-s and p-h figures,
# so the three sit together without one looking stretched.
fig, ax = plt.subplots(figsize=(6.4, 5.0))

# Straight segments between the computed points. No interpolation: every
# line on this figure connects two values the model actually produced.
for key, label, colour, marker, ls in SERIES:
    y = np.asarray(cop[key], float)
    ok = np.isfinite(y)
    lw = 2.4 if key == "Helmholtz" else 1.3
    ms = 5 if key == "Helmholtz" else 4.5
    ax.plot(amb[ok], y[ok], ls, color=colour, marker=marker, ms=ms, lw=lw,
            mfc="none" if key != "Helmholtz" else colour, mew=1.1, label=label)

ax.set_xlabel("Ambient temperature [deg C]")
ax.set_ylabel("COP")
ax.set_xticks(amb)
# Pad the y range so the curves do not run into the frame. The data
# spans 2.80-3.95; 5 % headroom each side keeps the spread readable
# without exaggerating it by cropping tight.
_all = np.concatenate([np.asarray(cop[k], float) for k, *_ in SERIES])
_lo, _hi = np.nanmin(_all), np.nanmax(_all)
_pad = 0.05 * (_hi - _lo)
ax.set_ylim(_lo - _pad, _hi + _pad)
ax.set_xlim(amb.min() - 1.0, amb.max() + 1.0)
ax.grid(True, ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
fig.tight_layout()
fig.savefig("cop_vs_ambient_poster_1008.png", dpi=200)
plt.close(fig)
print("\nwrote cop_vs_ambient_poster_1008.png")
