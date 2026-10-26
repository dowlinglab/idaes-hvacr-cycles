"""
poster_figures_1008.py -- R-32 T-s and p-h diagrams in poster formatting.

Same five series, same colours and line weights as the EARTH poster figure:

    Linde (datasheet)       black, dashed, small round markers
    NIST                    thin blue
    GCGP                    thin green
    SPGP                    thin red      (= fp_single_0919, relabelled)
    Helmholtz (CoolProp)    thick black

Differences from compare_cp_methods_0929_Helmholtz_SPGPlabel.py, which is
where the curves come from:
    - no figure title (the poster supplies its own heading)
    - axis labels in the poster's wording: "T [deg C]", "s [kJ/kg/K]",
      "h [kJ/kg]", "P [bar]"
    - legend inside, upper left, small
    - sized for a poster panel rather than a report page

This script IMPORTS that module rather than re-deriving anything, so the
curves are identical to the ones already validated. Importing it also re-runs
its own figures and error table -- harmless, and the table is a useful check
that nothing has drifted.

    python poster_figures_1008.py

Writes ts_diagram_poster_1008.png and ph_diagram_poster_1008.png.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import compare_cp_methods_0929_Helmholtz_SPGPlabel as M

DPI = 200
FIGSIZE = (6.4, 5.0)          # poster panel proportions

# series name -> (colour, linewidth)
STYLE = {
    "NIST":  ("tab:blue",  1.3),
    "GCGP":  ("tab:green", 1.3),
    "SPGP":  ("tab:red",   1.3),
}
HELM_STYLE = ("k", 2.4)
# Linde as OPEN markers with no connecting line, drawn LAST so it sits on
# top. Filled black dots on a black dashed line were invisible: Linde's h and
# s are Tillner-Roth & Yokozeki EoS output, the same equation CoolProp uses,
# so the two curves lie exactly on top of each other. Open circles read as
# discrete reference points sitting ON the Helmholtz curve, which is also the
# honest picture of what the data is.
LINDE_STYLE = dict(color="k", ls="none", marker="o", ms=5.5,
                   mfc="white", mec="k", mew=1.2)


def dome_xy(d, kind):
    """Liquid branch then reversed vapour branch, for a closed dome."""
    if kind == "ts":
        return (np.r_[d["sf"], d["sg"][::-1]],
                np.r_[d["T_C"], d["T_C"][::-1]])
    return (np.r_[d["hf"], d["hg"][::-1]],
            np.r_[d["Pd"], d["Pd"][::-1]])


# =============================================================================
# T-s
# =============================================================================
fig, ax = plt.subplots(figsize=FIGSIZE)

for d in M.domes:
    colour, lw = STYLE[d["name"]]
    x, y = dome_xy(d, "ts")
    ax.plot(x, y, color=colour, lw=lw, label=d["name"])
colour, lw = HELM_STYLE
x, y = dome_xy(M.helm_dome, "ts")
ax.plot(x, y, color=colour, lw=lw, label="Helmholtz (CoolProp)")
ax.plot(np.r_[M.L_SF, M.L_SG[::-1]], np.r_[M.L_TC, M.L_TC[::-1]],
        label="Linde (datasheet)", zorder=10, **LINDE_STYLE)

ax.set_xlabel("s [kJ/kg/K]")
ax.set_ylabel("T [deg C]")
ax.grid(True, ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
fig.tight_layout()
fig.savefig("ts_diagram_poster_1008.png", dpi=DPI)
plt.close(fig)
print("wrote ts_diagram_poster_1008.png")


# =============================================================================
# p-h
# =============================================================================
fig, ax = plt.subplots(figsize=FIGSIZE)

for d in M.domes:
    colour, lw = STYLE[d["name"]]
    x, y = dome_xy(d, "ph")
    ax.semilogy(x, y, color=colour, lw=lw, label=d["name"])
colour, lw = HELM_STYLE
x, y = dome_xy(M.helm_dome, "ph")
ax.semilogy(x, y, color=colour, lw=lw, label="Helmholtz (CoolProp)")
ax.semilogy(np.r_[M.L_HF, M.L_HG[::-1]], np.r_[M.L_P, M.L_P[::-1]],
            label="Linde (datasheet)", zorder=10, **LINDE_STYLE)

ax.set_xlabel("h [kJ/kg]")
ax.set_ylabel("P [bar]")
ax.grid(True, which="both", ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="lower center", framealpha=0.9)
fig.tight_layout()
fig.savefig("ph_diagram_poster_1008.png", dpi=DPI)
plt.close(fig)
print("wrote ph_diagram_poster_1008.png")

print("\nNOTE: 'SPGP' here is fp_single_0919 -- the Colon group's 09/19")
print("      delivery under the poster's display label, same convention as")
print("      compare_cp_methods_0903data_SPGPlabel.py.")
