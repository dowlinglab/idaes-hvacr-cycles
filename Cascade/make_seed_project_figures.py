"""
make_seed_project_figures.py

Makes the two "Seed Project: Optimization Results" figures for the EARTH
DVCT deck (slide 5):

    seed_project_cop.png             COP vs ambient temperature
    seed_project_pressure_ratio.png  compressor pressure ratio vs ambient
                                     (R-134a comparison only)

Style: journal / MATLAB-like -- boxed legend inside the axes, every line
named in full, refrigerant names per ASHRAE 34 (R-134a, R-1234yf, R-515B).

Data -- both from the "matched settings" runs (evap -29 C, condenser
approach 9 C, compressor eta 0.80, SH/SC cap 8 K):
    Cascade/cascade_vs_single_stage.csv
        cascade COP + pressure ratios, Carnot COP, single-stage COP
        (written by compare_cascade_vs_single_stage.py)
    R1234yf/cop_vs_ambient_plr_benchmark_combined.csv
        single-stage R-134a pressure ratio
        (written by R1234yf/run_plr_benchmark_r1234yf_vs_r134a.py)

Run from anywhere:
    python Cascade/make_seed_project_figures.py
The PNGs are written next to this script.
"""

import os

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ------------------------------------------------------------------ paths
# Build every path from this file's own location, so the script works no
# matter which folder you run it from.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CASCADE_CSV = os.path.join(THIS_DIR, "cascade_vs_single_stage.csv")
SINGLE_CSV = os.path.join(THIS_DIR, "..", "R1234yf",
                          "cop_vs_ambient_plr_benchmark_combined.csv")

# ----------------------------------------------------------------- colors
# One color per working fluid, same on both figures.
# Cascade (blue) = the cascade system on the COP plot, and its CO2 stage on
# the pressure-ratio plot. The cascade's high-temperature stage runs on
# R-134a, so on the pressure-ratio plot it takes R-134a's purple.
CASCADE = "#2a78d6"   # blue
R134A = "#4a3aa7"     # purple
R1234YF = "#eb6834"   # orange
R515B = "#1baf7a"     # green
CAP_RED = "#e34948"
INK = "#0b0b0b"       # text
GRID = "#d8d7d2"
SPINE = "#9a9990"

LINE_WIDTH = 2.2
MARKER_SIZE = 8

# ------------------------------------------------------------------- data
casc = pd.read_csv(CASCADE_CSV).sort_values("ambient_C")
single = pd.read_csv(SINGLE_CSV).sort_values("ambient_C")


# ------------------------------------------------------- helper functions
def new_figure():
    """Blank figure with a boxed, gridded axis."""
    fig, ax = plt.subplots(figsize=(6.6, 4.9), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.grid(True, color=GRID, linewidth=0.7, zorder=0)
    for spine in ax.spines.values():          # full box around the plot
        spine.set_visible(True)
        spine.set_color(SPINE)
        spine.set_linewidth(1.0)
    ax.tick_params(colors="#3d3d3a", labelsize=12)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(5))
    ax.set_xlim(8, 47)
    ax.set_xlabel("Ambient temperature (deg C)", fontsize=14, color=INK)
    return fig, ax


def plot_line(ax, x, y, color, marker, label, linestyle="-"):
    """One data series: line + filled markers, all in one color."""
    ax.plot(x, y, linestyle=linestyle, linewidth=LINE_WIDTH,
            marker=marker, markersize=MARKER_SIZE,
            color=color, markerfacecolor=color, markeredgecolor=color,
            label=label, zorder=3)


def boxed_legend(ax, location):
    """MATLAB-style legend: square frame, inside the axes."""
    legend = ax.legend(loc=location, frameon=True, fancybox=False,
                       edgecolor="#3d3d3a", framealpha=1.0, fontsize=11,
                       handlelength=2.6, borderpad=0.6, labelspacing=0.45)
    legend.get_frame().set_linewidth(0.8)
    for text in legend.get_texts():
        text.set_color(INK)


def save(fig, filename):
    fig.tight_layout()
    path = os.path.join(THIS_DIR, filename)
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print("saved:", path)


# ------------------------------------------------------- figure 1: COP
fig, ax = new_figure()
T = casc["ambient_C"]

ax.plot(T, casc["cop_carnot"], "--", color=INK, linewidth=1.6,
        label="Carnot limit", zorder=2)
plot_line(ax, T, casc["cascade_cop_full"], CASCADE, "o", "Cascade R-134a/CO$_2$")
plot_line(ax, T, casc["r134a_cop_full"], R134A, "s", "Single-stage R-134a")
plot_line(ax, T, casc["r1234yf_cop_full"], R1234YF, "s", "Single-stage R-1234yf")
plot_line(ax, T, casc["r515b_cop_full"], R515B, "s", "Single-stage R-515B")

ax.set_ylabel("COP", fontsize=14, color=INK)
ax.set_ylim(1, 7.5)                 # headroom so the legend clears the data
boxed_legend(ax, "upper right")
save(fig, "seed_project_cop.png")

# -------------------------------------- figure 2: compressor pressure ratio
# R-134a comparison only: single-stage R-134a vs the R-134a/CO2 cascade.
fig, ax = new_figure()

plot_line(ax, single["ambient_C"], single["r134a_comp_ratioP"],
          R134A, "s", "Single-stage R-134a")
plot_line(ax, T, casc["cascade_ratioP_hot"],
          R134A, "o", "Cascade high-temperature stage (R-134a)")
plot_line(ax, T, casc["cascade_ratioP_cold"],
          CASCADE, "o", "Cascade low-temperature stage (CO$_2$)", linestyle="--")
ax.axhline(8.0, color=CAP_RED, linestyle=":", linewidth=1.8,
           label="Pressure ratio cap (8)", zorder=2)

ax.set_ylabel("Compressor pressure ratio", fontsize=14, color=INK)
ax.set_ylim(0, 20)
ax.yaxis.set_major_locator(mticker.MultipleLocator(2))
boxed_legend(ax, "upper left")
save(fig, "seed_project_pressure_ratio.png")
