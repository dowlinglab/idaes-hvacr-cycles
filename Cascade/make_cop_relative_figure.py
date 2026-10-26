"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/16/2026

What this does
    Single stage cycle AND CO2 cascade, both cases on ONE plot (no condensing limit and
    condensing >= 21 C). Two figures:
        1. optimal COP versus ambient temperature
        2. compressor pressure ratio versus ambient temperature
    R134a is the curve. R515B and R1234yf are points. No titles.
    The Carnot COP is the grey line in the COP figure.

    Line / marker style tells the case:
        no limit      solid line, filled points
        21 C limit    dashed line, open points

Inputs (R1234yf folder, written by run_single_stage_sweep.py)
    sweep_single_<fluid>_nolimit.csv
    sweep_single_<fluid>_21C.csv
    columns used: ambient_C, cop_full, cop_carnot, ratio

Inputs for the cascade (this folder, written by the cascade sweeps)
    sweep_cascade_<fluid>_nolimit.csv, sweep_cascade_<fluid>_21C.csv
    columns used: ambient_C, cop_full, cop_carnot, hot_ratio, cold_ratio
    (the R515B files call the ratios hot_comp_ratioP and cold_comp_ratioP)
    The cascade pressure ratio plotted is the HIGHER of the two compressors.

Outputs (png and pdf each, in the current folder)
    single_stage_cop_vs_ambient
    single_stage_pressure_ratio_vs_ambient
    cascade_cop_vs_ambient
    cascade_pressure_ratio_vs_ambient
    cop_bars_no_condensing_limit
    cop_bars_condensing_limit

Bar figures
    Grouped bars of COP at selected ambient temperatures, single stage next to
    cascade for each fluid. Colour = fluid. Solid bar = single stage,
    hatched open bar = cascade.

Run
    cd ~/idaes-hvacr-cycles/Cascade
    conda activate myidaesenv
    python make_cop_relative_figure.py
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# =====================================================================
# 1. SETTINGS  (change things here, nowhere else)
# =====================================================================

DATA_DIR = os.path.join("..", "R1234yf")   # folder with the sweep_single csv files
DATA_DIR_CASCADE = "."                      # folder with the sweep_cascade csv files

CASES = ["nolimit", "21C"]
FLUIDS = ["R134a", "R515B", "R1234yf"]

# colors (no orange, no green): black, blue, purple, grey
COLORS = {"R134a": "black", "R515B": "#0072B2", "R1234yf": "#AA3377"}
# names shown on the figures
FLUID_LABEL = {"R134a": "R-134a", "R515B": "R-515B", "R1234yf": "R-1234yf"}
CO2_LABEL = "R-744"          # CO2 in the cascade

COLOR_CARNOT = "0.30"          # darker than before so it is easy to see

# R134a is a curve; the other two are points
MARKERS = {"R515B": "o", "R1234yf": "^"}
MARKER_SIZE = {"R515B": 5.5, "R1234yf": 6.5}

# single-stage pressure ratio limit drawn on the pressure ratio figure
# (literature: about 8 to 9 for a single-stage reciprocating compressor)
SINGLE_STAGE_RATIO_LIMIT = 8.0
COLOR_RATIO_LIMIT = "red"

LIMIT_AMBIENT_C = 12.0   # below this ambient the 21 C limit holds condensing fixed

# extra empty room at the top of the y axis, so the legend does not sit on the data
# (1.0 = none; the COP curve is high only on the left, so it needs none)
HEADROOM_COP = 1.0
HEADROOM_RATIO = 1.65

# bar figures: which ambients to show, and room at the top for the legend
AMBIENTS_BAR = [-30, -20, -10, 0, 10, 20, 30, 40, 50]
HEADROOM_BAR = {"nolimit": 1.08, "21C": 1.7}
BAR_WIDTH = 0.12
PAIR_GAP = 1.25   # spacing between fluid pairs, in units of the pair width

plt.rcParams["font.size"] = 11
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["DejaVu Serif"]


# =====================================================================
# 2. LOAD THE SWEEP FILES
#    data[(fluid, case)] is a table with one row per ambient temperature
# =====================================================================

data = {}
for fluid in FLUIDS:
    for case in CASES:
        file_name = "sweep_single_" + fluid + "_" + case + ".csv"
        data[(fluid, case)] = pd.read_csv(os.path.join(DATA_DIR, file_name))

# cascade tables, same layout: cascade_data[(fluid, case)]
cascade_data = {}
for fluid in FLUIDS:
    for case in CASES:
        file_name = "sweep_cascade_" + fluid + "_" + case + ".csv"
        table = pd.read_csv(os.path.join(DATA_DIR_CASCADE, file_name))
        # the R515B files name the two compressor ratios differently
        table = table.rename(columns={"hot_comp_ratioP": "hot_ratio",
                                      "cold_comp_ratioP": "cold_ratio"})
        # one number per ambient: the higher of the two compressor ratios
        table["ratio"] = table[["hot_ratio", "cold_ratio"]].max(axis=1)
        cascade_data[(fluid, case)] = table


# =====================================================================
# 3. ONE FUNCTION DRAWS ONE FIGURE (y_column picks COP or pressure ratio)
# =====================================================================

def draw_figure(tables, y_column, y_label, out_name, with_carnot, headroom,
                ratio_limit=None, y_limits=None, label_tail=""):

    fig, ax = plt.subplots(figsize=(7, 5.2))

    for case in CASES:

        # solid / filled for no limit, dashed / open for the 21 C limit
        if case == "nolimit":
            line_style = "-"
            fill = True
        else:
            line_style = "--"
            fill = False

        # Carnot (grey), same line style as the case
        if with_carnot:
            table = tables[("R134a", case)]
            ax.plot(table["ambient_C"], table["cop_carnot"], line_style,
                    color=COLOR_CARNOT, linewidth=2.0)

        for fluid in FLUIDS:
            table = tables[(fluid, case)]
            T = table["ambient_C"]
            y = table[y_column]

            if fluid == "R134a":
                ax.plot(T, y, line_style, color=COLORS[fluid], linewidth=1.8)
            else:
                if fill:
                    face = COLORS[fluid]
                else:
                    face = "none"
                ax.plot(T, y, linestyle="none", marker=MARKERS[fluid],
                        markersize=MARKER_SIZE[fluid], markerfacecolor=face,
                        markeredgecolor=COLORS[fluid], markeredgewidth=1.4)

    # mark where the 21 C limit starts to act
    ax.axvline(LIMIT_AMBIENT_C, color="0.6", linestyle=":", linewidth=1)

    # single-stage pressure ratio limit (pressure ratio figure only)
    if ratio_limit is not None:
        ax.axhline(ratio_limit, color=COLOR_RATIO_LIMIT, linestyle="-.", linewidth=1.4)

    ax.set_xlabel(r"$T_{\mathrm{amb}}\ (^\circ\mathrm{C})$")
    ax.set_ylabel(y_label)
    ax.grid(True, color="0.8", linewidth=0.6)
    for side in ax.spines.values():
        side.set_visible(True)
        side.set_linewidth(1.0)

    # legend 1: which fluid
    fluid_handles = [
        Line2D([0], [0], color="black", linewidth=1.8, label=FLUID_LABEL["R134a"] + label_tail),
        Line2D([0], [0], linestyle="none", marker="o", markersize=5.5,
               color=COLORS["R515B"], label=FLUID_LABEL["R515B"] + label_tail),
        Line2D([0], [0], linestyle="none", marker="^", markersize=6.5,
               color=COLORS["R1234yf"], label=FLUID_LABEL["R1234yf"] + label_tail),
    ]
    if with_carnot:
        fluid_handles.append(Line2D([0], [0], color=COLOR_CARNOT, linewidth=2.0,
                                    label="Carnot (ideal)"))

    if ratio_limit is not None:
        fluid_handles.append(Line2D([0], [0], color=COLOR_RATIO_LIMIT, linestyle="-.",
                                    linewidth=1.4,
                                    label=r"$\mathrm{PR} = " + str(int(ratio_limit)) + "$"))

    # legend 2: which case
    case_handles = [
        Line2D([0], [0], color="0.3", linestyle="-", marker="o", markersize=5,
               markerfacecolor="0.3", label="Solid, filled: no condensing limit"),
        Line2D([0], [0], color="0.3", linestyle="--", marker="o", markersize=5,
               markerfacecolor="none", label=r"Dashed, open: condensing $\geq 21\,^\circ$C"),
        Line2D([0], [0], color="0.6", linestyle=":", linewidth=1.2,
               label=r"Dotted: limit starts at 12$\,^\circ$C ambient"),
    ]

    # one legend, top right, inside the plot
    all_handles = fluid_handles + case_handles
    ax.legend(handles=all_handles, loc="upper right",
              fontsize=9, frameon=True, edgecolor="black")

    # leave empty room at the top so the legend does not sit on the data
    y_low, y_high = ax.get_ylim()
    ax.set_ylim(y_low, y_low + (y_high - y_low) * headroom)

    # optional: use the y axis of another figure so the two can be compared
    if y_limits is not None:
        ax.set_ylim(y_limits)
    final_limits = ax.get_ylim()

    fig.tight_layout()
    fig.savefig(out_name + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out_name + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print("saved", out_name + ".png and .pdf")
    return final_limits


# =====================================================================
# 4. MAKE THE FIGURES
#    the cascade figures reuse the single stage y axes, so they compare directly
# =====================================================================

cop_limits = draw_figure(data, "cop_full", r"$\mathrm{COP}$",
                         "single_stage_cop_vs_ambient", with_carnot=True,
                         headroom=HEADROOM_COP)

ratio_limits = draw_figure(data, "ratio",
                           r"Pressure ratio $P_{\mathrm{high}}/P_{\mathrm{low}}$",
                           "single_stage_pressure_ratio_vs_ambient", with_carnot=False,
                           headroom=HEADROOM_RATIO, ratio_limit=SINGLE_STAGE_RATIO_LIMIT)

draw_figure(cascade_data, "cop_full", r"$\mathrm{COP}$",
            "cascade_cop_vs_ambient", with_carnot=True,
            headroom=HEADROOM_COP, y_limits=cop_limits, label_tail="/" + CO2_LABEL)

draw_figure(cascade_data, "ratio", r"Highest compressor pressure ratio",
            "cascade_pressure_ratio_vs_ambient", with_carnot=False,
            headroom=HEADROOM_RATIO, ratio_limit=SINGLE_STAGE_RATIO_LIMIT,
            y_limits=ratio_limits, label_tail="/" + CO2_LABEL)


# =====================================================================
# 5. BAR FIGURES: COP, single stage next to cascade, one figure per case
# =====================================================================

def draw_bar_figure(case, out_name):

    fig, ax = plt.subplots(figsize=(9, 5.2))

    # six bars per ambient: for each fluid, a pair (single stage, then cascade)

    for group, ambient in enumerate(AMBIENTS_BAR):
        for i, fluid in enumerate(FLUIDS):

            single = data[(fluid, case)].set_index("ambient_C")["cop_full"][ambient]
            cascade = cascade_data[(fluid, case)].set_index("ambient_C")["cop_full"][ambient]

            # centre of this fluid's pair, then the two bars either side of it
            pair_centre = group + (i - 1) * 2 * BAR_WIDTH * PAIR_GAP
            x_single = pair_centre - BAR_WIDTH / 2
            x_cascade = pair_centre + BAR_WIDTH / 2

            # single stage: solid fill
            ax.bar(x_single, single, width=BAR_WIDTH,
                   color=COLORS[fluid], edgecolor=COLORS[fluid], linewidth=1.2)
            # cascade: open bar with hatching
            ax.bar(x_cascade, cascade, width=BAR_WIDTH, facecolor="white",
                   edgecolor=COLORS[fluid], hatch="////", linewidth=1.2)

    ax.set_xticks(range(len(AMBIENTS_BAR)))
    ax.set_xticklabels([str(t).replace("-", "\u2212") for t in AMBIENTS_BAR])
    ax.set_xlabel(r"$T_{\mathrm{amb}}\ (^\circ\mathrm{C})$")
    ax.set_ylabel(r"$\mathrm{COP}$")
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color="0.8", linewidth=0.6)
    for side in ax.spines.values():
        side.set_visible(True)
        side.set_linewidth(1.0)

    # legend: colour = fluid, fill = cycle
    handles = [
        Patch(facecolor=COLORS["R134a"], edgecolor=COLORS["R134a"], label=FLUID_LABEL["R134a"]),
        Patch(facecolor=COLORS["R515B"], edgecolor=COLORS["R515B"], label=FLUID_LABEL["R515B"]),
        Patch(facecolor=COLORS["R1234yf"], edgecolor=COLORS["R1234yf"], label=FLUID_LABEL["R1234yf"]),
        Patch(facecolor="0.4", edgecolor="0.4", label="Single stage (solid)"),
        Patch(facecolor="white", edgecolor="0.4", hatch="////",
              label="Cascade with " + CO2_LABEL + " (hatched)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=9,
              frameon=True, edgecolor="black")

    y_low, y_high = ax.get_ylim()
    ax.set_ylim(0, y_high * HEADROOM_BAR[case])

    fig.tight_layout()
    fig.savefig(out_name + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out_name + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print("saved", out_name + ".png and .pdf")


draw_bar_figure("nolimit", "cop_bars_no_condensing_limit")
draw_bar_figure("21C", "cop_bars_condensing_limit")
