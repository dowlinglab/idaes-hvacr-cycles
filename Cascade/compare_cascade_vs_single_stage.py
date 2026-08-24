"""
compare_cascade_vs_single_stage.py

Compares the R134a/CO2 cascade against the existing single-stage
R134a / R1234yf / R515B benchmark, at MATCHED conditions:

    evap_sat = -29 C, condenser approach = 9 C, PLR = 0.75, CD = 0.13,
    ambient swept 10-45 C in 5 C steps

Single-stage numbers are read from the existing benchmark output,
`../R1234yf/cop_vs_ambient_plr_benchmark_combined.csv` (produced by
R1234yf/run_plr_benchmark_r1234yf_vs_r134a.py). The cascade is solved here.

Writes: cascade_vs_single_stage.csv, cascade_vs_single_stage.png
"""

import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from validate_cascade import (stage3_ambient_sweep, BENCHMARK_EVAP_SAT_C,
                              BENCHMARK_FLUIDS)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SINGLE_STAGE_CSV = os.path.join(
    _THIS_DIR, "..", "R1234yf", "cop_vs_ambient_plr_benchmark_combined.csv")


def load_single_stage(path=SINGLE_STAGE_CSV):
    if not os.path.exists(path):
        print(f"  (single-stage benchmark not found at {path} -- cascade only)")
        return {}
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            out[float(r["ambient_C"])] = {
                "r134a_cop": float(r["r134a_cop_full"]),
                "r134a_ratioP": float(r["r134a_comp_ratioP"]),
                "r1234yf_cop": float(r["r1234yf_cop_full"]),
                "r515b_cop": float(r["r515b_cop_full"]),
                "cop_carnot": float(r["cop_carnot"]),
            }
    return out


def main():
    # MUST use BENCHMARK_FLUIDS (eta 0.80, SH/SC cap 8K) -- matching the
    # single-stage benchmark's own settings. Using the cascade's as-designed
    # DEFAULT_FLUIDS here would compare eta 0.75 against eta 0.80 and make
    # the cascade look ~7% worse purely from the mismatch. See BREADCRUMB.md.
    cascade_rows = stage3_ambient_sweep(cold_evap_C=BENCHMARK_EVAP_SAT_C,
                                        fluids=BENCHMARK_FLUIDS,
                                        label="benchmark-matched")
    single = load_single_stage()

    merged = []
    for r in cascade_rows:
        if not r["converged"]:
            continue
        amb = r["ambient_C"]
        s = single.get(amb, {})
        row = {
            "ambient_C": amb,
            "cascade_cop_full": r["cop_full"],
            "cascade_cop_part": r["cop_part"],
            "cascade_ratioP_hot": r["hot_ratioP"],
            "cascade_ratioP_cold": r["cold_ratioP"],
            "cascade_ratioP_max": max(r["hot_ratioP"], r["cold_ratioP"]),
            "cascade_intermediate_T_C": r["intermediate_T_C"],
            "cop_carnot": r["cop_carnot"],
        }
        row.update({
            "r134a_cop_full": s.get("r134a_cop"),
            "r134a_ratioP": s.get("r134a_ratioP"),
            "r1234yf_cop_full": s.get("r1234yf_cop"),
            "r515b_cop_full": s.get("r515b_cop"),
        })
        if s.get("r134a_cop"):
            row["cascade_vs_r134a_pct"] = 100.0 * (r["cop_full"] / s["r134a_cop"] - 1.0)
        merged.append(row)

    out_csv = os.path.join(_THIS_DIR, "cascade_vs_single_stage.csv")
    if merged:
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(merged[0].keys()))
            w.writeheader()
            w.writerows(merged)
        print(f"\nWrote {out_csv}")

    print()
    print("  COP (full load), matched conditions")
    print(f"  {'amb':>4} {'cascade':>9} {'R134a':>9} {'R1234yf':>9} {'R515B':>9} "
          f"{'casc vs R134a':>14}")
    for m in merged:
        def f(x):
            return f"{x:9.4f}" if x is not None else f"{'--':>9}"
        d = m.get("cascade_vs_r134a_pct")
        dstr = f"{d:+.1f}%" if d is not None else "--"
        print(f"  {m['ambient_C']:4.0f} {f(m['cascade_cop_full'])} {f(m['r134a_cop_full'])} "
              f"{f(m['r1234yf_cop_full'])} {f(m['r515b_cop_full'])} {dstr:>14}")

    print()
    print("  Compressor pressure ratio -- the cascade's real advantage here")
    print(f"  {'amb':>4} {'1-stage R134a':>14} {'cascade hot':>12} {'cascade cold':>13} "
          f"{'cascade max':>12}")
    for m in merged:
        r1 = m.get("r134a_ratioP")
        r1str = f"{r1:.3f}" if r1 else "--"
        print(f"  {m['ambient_C']:4.0f} {r1str:>14} "
              f"{m['cascade_ratioP_hot']:12.3f} {m['cascade_ratioP_cold']:13.3f} "
              f"{m['cascade_ratioP_max']:12.3f}")

    # ---- plot ----
    amb = [m["ambient_C"] for m in merged]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.plot(amb, [m["cop_carnot"] for m in merged], "k--", lw=1, label="Carnot")
    ax1.plot(amb, [m["cascade_cop_full"] for m in merged], "o-", lw=2,
             label="Cascade R134a/CO2")
    if any(m.get("r134a_cop_full") for m in merged):
        ax1.plot(amb, [m["r134a_cop_full"] for m in merged], "s-", label="1-stage R134a")
        ax1.plot(amb, [m["r1234yf_cop_full"] for m in merged], "^-", label="1-stage R1234yf")
        ax1.plot(amb, [m["r515b_cop_full"] for m in merged], "v-", label="1-stage R515B")
    ax1.set_xlabel("Ambient temperature (C)")
    ax1.set_ylabel("COP (full load)")
    ax1.set_title("COP vs ambient -- MATCHED settings\n(evap -29 C, approach 9 C, eta 0.80, SH/SC cap 8 K)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    if any(m.get("r134a_ratioP") for m in merged):
        ax2.plot(amb, [m["r134a_ratioP"] for m in merged], "s-", label="1-stage R134a")
    ax2.plot(amb, [m["cascade_ratioP_hot"] for m in merged], "o-", label="Cascade hot (R134a)")
    ax2.plot(amb, [m["cascade_ratioP_cold"] for m in merged], "o--", label="Cascade cold (CO2)")
    ax2.axhline(8.0, color="r", ls=":", label="ratio cap = 8")
    ax2.set_xlabel("Ambient temperature (C)")
    ax2.set_ylabel("Compressor pressure ratio")
    ax2.set_title("Pressure ratio vs ambient\nCascade keeps both stages under the cap")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(_THIS_DIR, "cascade_vs_single_stage.png")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"\nWrote {out_png}")


if __name__ == "__main__":
    main()
