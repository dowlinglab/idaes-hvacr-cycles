#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/testing: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Overlay REFPROP-style isotherms on the existing R515B dome plot.
Dependencies: numpy, pandas, matplotlib, scripts/plot_refprop_style_isotherm.py
Context reference: PROJECT_CONTEXT.md

BREADCRUMB:
- Date: 2026-03-03
- Context: User requested isotherms shown with dome after removing old Layer-4 sweep.
- Scope: plotting-only composition of dome layers with REFPROP-style isotherm traces.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plot_refprop_style_isotherm import generate_refprop_style_isotherm


def main() -> None:
    """
    Build a dome plot with REFPROP-style isotherm overlays.
    """
    p = argparse.ArgumentParser(description="Overlay REFPROP-style isotherms on R515B dome.")
    p.add_argument("--temps-c", default="-19,-5,10,25,40,55,70,85")
    p.add_argument("--w1", type=float, default=0.911)
    p.add_argument("--Pmin-bar", type=float, default=1.0)
    p.add_argument("--Pmax-bar", type=float, default=100.0)
    p.add_argument("--nv", type=int, default=60)
    p.add_argument("--nl", type=int, default=60)
    p.add_argument("--eps-kpa", type=float, default=1.0)
    args = p.parse_args()

    stamp = datetime.now().strftime("%Y%m%d")
    out_dir = ROOT / "diagnostics/plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fig = out_dir / f"ph_dome_with_refprop_style_isotherms_{stamp}.png"
    out_csv = out_dir / f"ph_dome_with_refprop_style_isotherms_{stamp}.csv"

    d2 = pd.read_csv(ROOT / "diagnostics/plots/ph_layer2_true_vle_boundary_20260303.csv").sort_values("T_K")
    d3 = pd.read_csv(ROOT / "diagnostics/plots/ph_layer3_beta_grid_20260303.csv")
    d5 = pd.read_csv(ROOT / "diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.csv")

    fig, ax = plt.subplots(figsize=(10.8, 6.6), dpi=190)
    ax.set_yscale("log")
    ax.set_xlim(150.0, 500.0)
    ax.set_ylim(1.0, 100.0)
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.grid(True, which="both", alpha=0.30)

    ax.plot(d2["h_kJkg_bubble"], d2["P_bar_bubble"], color="#0b4f6c", lw=2.0, label="True VLE bubble (liq)")
    ax.plot(d2["h_kJkg_dew"], d2["P_bar_dew"], color="#c44536", lw=2.0, label="True VLE dew (vap)")

    beta0 = d3[np.isclose(d3["beta"], 0.0)][["pair_index", "P_pair_bar", "h_beta_kJkg"]].rename(
        columns={"h_beta_kJkg": "h_l_kJkg"}
    )
    beta1 = d3[np.isclose(d3["beta"], 1.0)][["pair_index", "h_beta_kJkg"]].rename(
        columns={"h_beta_kJkg": "h_v_kJkg"}
    )
    band = beta0.merge(beta1, on="pair_index").sort_values("P_pair_bar")
    ax.fill_betweenx(
        band["P_pair_bar"], band["h_l_kJkg"], band["h_v_kJkg"],
        color="#f4d35e", alpha=0.18, label="Two-phase beta-band"
    )

    d5_ok = d5[d5["status"] == "OK"].copy()
    ax.plot(d5_ok["h_l_kJkg_pseudopure"], d5_ok["P_bar"], linestyle=":", color="#6a4c93", lw=2.0, label="Pseudo-pure overlay")
    ax.plot(d5_ok["h_v_kJkg_pseudopure"], d5_ok["P_bar"], linestyle=":", color="#6a4c93", lw=2.0)

    rows = []
    temps = [float(x.strip()) for x in args.temps_c.split(",") if x.strip()]
    cmap = plt.get_cmap("viridis")
    for k, t_c in enumerate(temps):
        df_iso, _ = generate_refprop_style_isotherm(
            T_C=t_c,
            w1=args.w1,
            P_min_bar=args.Pmin_bar,
            P_max_bar=args.Pmax_bar,
            n_vapor=args.nv,
            n_liquid=args.nl,
            eps_kPa=args.eps_kpa,
        )
        color = cmap(k / max(1, len(temps) - 1))
        dv = df_iso[df_iso["phase_flag"] == "vapor"]
        dl = df_iso[df_iso["phase_flag"] == "liquid"]
        dc = df_iso[df_iso["phase_flag"] == "two_phase_connector"]
        if len(dv):
            ax.plot(dv["h_kJkg"], dv["P_bar"], color=color, lw=1.0, alpha=0.9)
        if len(dl):
            ax.plot(dl["h_kJkg"], dl["P_bar"], color=color, lw=1.0, alpha=0.9)
        if len(dc) == 2:
            ax.plot(dc["h_kJkg"], dc["P_bar"], color=color, lw=0.9, alpha=0.75, linestyle="--")
        ax.plot([], [], color=color, lw=1.5, label=f"Isotherm {t_c:.0f}C")
        for _, r in df_iso.iterrows():
            rows.append({"T_C": t_c, "P_bar": float(r["P_bar"]), "h_kJkg": float(r["h_kJkg"]), "phase_flag": str(r["phase_flag"])})

    ax.set_title("R515B P-h Dome with REFPROP-style Isotherms")
    ax.legend(loc="lower right", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")

    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"saved_figure={out_fig}")
    print(f"saved_csv={out_csv}")


if __name__ == "__main__":
    main()
