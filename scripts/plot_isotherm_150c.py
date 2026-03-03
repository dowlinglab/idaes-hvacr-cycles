#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot R515x-style isotherm at 150 degC for density sweep 0.1-0.5 g/cc.

Output plot axes are limited to:
- Pressure: 20-100 bar
- Enthalpy: 410-500 kJ/kg
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pressure_validated_model import (
    ChartReference,
    chart_offset,
    compute_props,
    default_interaction,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot 150C isotherm in requested P-h window.")
    parser.add_argument("--h-basis", choices=["raw", "chart"], default="raw", help="Enthalpy basis for plotting")
    parser.add_argument("--pref", type=float, default=None, help="Reference saturation pressure [kPa] at Tref for liquid-root rho_ref computation")
    parser.add_argument("--Tref", type=float, default=273.15, help="Reference temperature [K] for enthalpy shift")
    parser.add_argument("--rhoref", type=float, default=1258.4, help="Fallback reference density [kg/m^3] if --pref is not provided")
    parser.add_argument("--href", type=float, default=200.0, help="Reference enthalpy [kJ/kg] at (Tref, rhoref or solved rho_ref)")
    args = parser.parse_args()

    fluid1 = "r1234ze"
    fluid2 = "r227ea"
    w1 = 0.911
    w2 = 1.0 - w1
    T_K = 150.0 + 273.15
    interaction = default_interaction()

    ref = ChartReference(
        T_ref_K=args.Tref,
        rho_ref_kgm3=args.rhoref,
        p_ref_kPa=args.pref,
        h_ref_kJkg=args.href,
    )

    rho_gcc = np.linspace(0.1, 0.5, 300)
    rho_kgm3 = rho_gcc * 1000.0

    p_bar = []
    h_kJkg = []
    rho_gcc_ok = []

    # Compute reference offset once for chart basis (if used).
    h_off = 0.0
    rho_ref_used = ref.rho_ref_kgm3
    if args.h_basis == "chart":
        h_off, rho_ref_used = chart_offset(
            ref=ref,
            comp1_json=f"{fluid1}.json",
            comp2_json=f"{fluid2}.json",
            w1=w1,
            w2=w2,
            interaction=interaction,
        )

    for rg, rk in zip(rho_gcc, rho_kgm3):
        st = compute_props(
            T=T_K,
            rho=rk,
            comp1_json=f"{fluid1}.json",
            comp2_json=f"{fluid2}.json",
            w1=w1,
            w2=w2,
            interaction=interaction,
        )
        h = st.h if args.h_basis == "raw" else st.h + h_off
        p = st.p / 100.0
        if 20.0 <= p <= 100.0 and 410.0 <= h <= 500.0:
            p_bar.append(p)
            h_kJkg.append(h)
            rho_gcc_ok.append(rg)

    fig, ax = plt.subplots(figsize=(8, 6))

    if len(p_bar) > 0:
        sc = ax.scatter(
            h_kJkg,
            p_bar,
            c=rho_gcc_ok,
            cmap="viridis",
            s=20,
            edgecolors="none",
        )
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label("Density [g/cc]")
    else:
        ax.text(
            0.5,
            0.5,
            "No points in requested window",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )

    ax.set_title(f"Isotherm at 150 degC (R515x model, {args.h_basis}-basis h)")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_xlim(410, 500)
    ax.set_ylim(20, 100)
    ax.grid(True, alpha=0.3)

    out = ROOT / "isotherm_150C_P20to100bar_H410to500_rho01to05gcc.png"
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    print(f"Saved figure: {out}")
    print(f"Points in window: {len(p_bar)}")
    print(f"Enthalpy basis: {args.h_basis}")
    print(f"Reference density used [kg/m^3]: {rho_ref_used:.6g}")
    print(f"Applied enthalpy offset [kJ/kg]: {h_off:.6g}")


if __name__ == "__main__":
    main()
