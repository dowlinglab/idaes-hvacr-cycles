"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/11/2026

What this does
    The two COP vs ambient figures (one per case) from the 12 sweep CSVs,
    plus one combined CSV per case (for Zhuorui).
        nolimit  condensing = ambient + 9 C
        21C      head pressure control, condensing >= 21 C

Each figure: six curves.
    Cascade (fluid/CO2)  solid line
    Single stage         dashed line
    One color + one marker per fluid (R134a, R1234yf, R515B)

Two versions per case: full-load COP, and part-load COP (PLF = 1 - CD(1 - PLR)
= 0.9675) with the Carnot COP as a black dotted line.

Only points that count are drawn: converged, certified, and passing every
physics check (physics_checks.py). A missing point leaves a gap.
Full-load and part-load versions (see above).

Inputs
    Cascade/sweep_cascade_<fluid>_<case>.csv   (run_cascade_sweep_hpc.py,
                                                run_cascade_sweep_r515b.py)
    R1234yf/sweep_single_<fluid>_<case>.csv    (run_single_stage_sweep.py)
Outputs (this folder)
    cop_vs_ambient_<case>.png / .pdf            full load
    cop_vs_ambient_<case>_part.png / .pdf       part load + Carnot
    cop_vs_ambient_all_<case>.csv / _part.csv

Run
    python make_cop_figures.py
"""

import os

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SINGLE_DIR = os.path.join(os.path.dirname(HERE), "R1234yf")

FLUIDS = {                       # color, marker, label
    "R134a": ("#4a3aa7", "o", "R-134a"),
    "R1234yf": ("#eb6834", "s", "R-1234yf"),
    "R515B": ("#1baf7a", "^", "R-515B"),
}
CASES = {"nolimit": "No condensing limit (condensing = ambient + 9 C)",
         "21C": "Head pressure control (condensing $\\geq$ 21 C)"}
INK, MUTED, GRID = "#0b0b0b", "#3d3d3a", "#d8d7d2"
EVAP_C, APPROACH_C, T_COND_MIN_C = -29.0, 9.0, 21.0     # for the Carnot line


def load(path):
    """Ambient + COP for the points that count. Handles both CSV styles."""
    if not os.path.exists(path):
        print("missing:", os.path.basename(path))
        return None
    df = pd.read_csv(path)
    if "counts" in df:
        good = df["counts"].astype(str) == "True"
    else:                                   # R515B cascade runner
        good = (df["converged"] == 1) & df["checks_failed"].fillna("").eq("")
    out = df.loc[good, ["ambient_C", "cop_full", "cop_part"]].sort_values("ambient_C")
    print(f"{os.path.basename(path):40s} {good.sum():2d}/{len(df)} points count")
    return out


def carnot(ambient_c, case):
    """Carnot COP = T_evap / (T_cond - T_evap), K. T_evap = -29 C; T_cond = ambient + 9,
    or max(ambient + 9, 21) with head pressure control. Same for every fluid and cycle."""
    t_evap = EVAP_C + 273.15
    t_cond = ambient_c + APPROACH_C + 273.15 if case == "nolimit" else max(ambient_c + APPROACH_C, T_COND_MIN_C) + 273.15
    return t_evap / (t_cond - t_evap)


def figure(case, which="full"):
    """which = "full": full-load COP. which = "part": part-load COP (PLF x COP) plus Carnot."""
    col_name = "cop_full" if which == "full" else "cop_part"
    fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=300)
    ax.grid(True, color=GRID, lw=0.7, zorder=0)
    for s in ax.spines.values():
        s.set_color("#9a9990")
    table = None
    for fluid, (color, marker, label) in FLUIDS.items():
        for cycle, ls, path in (
                ("Cascade", "-", os.path.join(HERE, f"sweep_cascade_{fluid}_{case}.csv")),
                ("Single stage", "--", os.path.join(SINGLE_DIR, f"sweep_single_{fluid}_{case}.csv"))):
            df = load(path)
            if df is None or df.empty:
                continue
            name = f"{cycle} {label}" + ("/CO$_2$" if cycle == "Cascade" else "")
            ax.plot(df.ambient_C, df[col_name], ls=ls, color=color, lw=2.0, marker=marker, ms=5.5,
                    mfc=color if cycle == "Cascade" else "white", mec=color, mew=1.4,
                    label=name, zorder=3)
            col = df[["ambient_C", col_name]].rename(columns={col_name: f"{cycle.lower().replace(' ', '_')}_{fluid}"})
            table = col if table is None else table.merge(col, on="ambient_C", how="outer")
    if which == "part":
        amb = list(range(-30, 55, 5))
        ax.plot(amb, [carnot(a, case) for a in amb], ls=":", color=INK, lw=2.0, label="Carnot", zorder=2)
        if table is not None:
            table["carnot"] = [carnot(a, case) for a in table["ambient_C"]]
    ax.set_xlim(-32, 52)
    ax.xaxis.set_major_locator(mt.MultipleLocator(10))
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Ambient temperature (C)", fontsize=13, color=INK)
    ax.set_ylabel("COP (full load)" if which == "full" else "COP (part load, PLR 0.75, CD 0.13)", fontsize=13, color=INK)
    ax.set_title(CASES[case], fontsize=12, color=INK)
    ax.tick_params(colors=MUTED, labelsize=11)
    # 21C: the flat part + kink are at the top, so the legend goes bottom left
    leg = ax.legend(loc="lower left" if case == "21C" else "upper right", fontsize=9, frameon=True, fancybox=False,
                    edgecolor=MUTED, framealpha=1, ncol=2)
    leg.get_frame().set_linewidth(0.8)
    fig.tight_layout()
    stem = f"cop_vs_ambient_{case}" + ("" if which == "full" else "_part")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{stem}.{ext}"), facecolor="white")
    plt.close(fig)
    if table is not None:
        table.sort_values("ambient_C").to_csv(
            os.path.join(HERE, f"cop_vs_ambient_all_{case}" + ("" if which == "full" else "_part") + ".csv"), index=False)
    print(f"saved {stem}.png/.pdf + table")


if __name__ == "__main__":
    for c in CASES:
        figure(c, "full")
        figure(c, "part")
