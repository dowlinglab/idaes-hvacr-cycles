#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/testing: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Generate REFPROP-style mixture P-h isotherm segments at fixed T
using single-phase roots plus a saturation-endpoint connector.
Dependencies: numpy, pandas, matplotlib, pressure_validated_model
Context reference: PROJECT_CONTEXT.md

BREADCRUMB:
- Date: 2026-03-03
- Context: R515B isotherm layer scrapped; replaced with REFPROP-style isotherms.
- Scope: plotting-only algorithm; EOS, VLE equations, and mixture model unchanged.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pressure_validated_model import compute_props, default_interaction, solve_rho_mass_for_P


def interpolate_endpoints_at_T(boundary_csv: Path, T_C: float) -> Dict[str, float]:
    """
    Interpolate bubble/dew endpoint coordinates at fixed temperature.

    Thermodynamic basis
    -------------------
    The two-phase region is represented by saturation endpoints only; no interior
    two-phase density sweep is used. This follows REFPROP-style property-plot logic
    where metastable interior traces are optional and generally non-physical.

    References
    ----------
    - Lemmon, E. W., Huber, M. L., & McLinden, M. O. (2018). NIST REFPROP
      Documentation, Version 10.0.
    - Bell, I. H. et al. (2014). Pure and Pseudo-pure Fluid Thermophysical
      Property Evaluation and the Open-Source Thermophysical Property Library
      CoolProp. Industrial & Engineering Chemistry Research, 53(6), 2498-2508.

    Breadcrumb context
    ------------------
    Added during R515B isotherm-layer replacement (2026-03-03) to provide
    saturation endpoint extraction for REFPROP-style plotting.
    """
    df = pd.read_csv(boundary_csv).sort_values("T_K")
    T_K = T_C + 273.15
    x = df["T_K"].to_numpy(dtype=float)
    pb = np.interp(T_K, x, df["P_bar_bubble"].to_numpy(dtype=float)) * 100.0
    pdw = np.interp(T_K, x, df["P_bar_dew"].to_numpy(dtype=float)) * 100.0
    hb = np.interp(T_K, x, df["h_kJkg_bubble"].to_numpy(dtype=float))
    hd = np.interp(T_K, x, df["h_kJkg_dew"].to_numpy(dtype=float))
    return {
        "T_K": float(T_K),
        "P_bub_kPa": float(pb),
        "h_bub_kJkg": float(hb),
        "P_dew_kPa": float(pdw),
        "h_dew_kJkg": float(hd),
    }


def single_phase_branch_points(
    T_K: float,
    P_kPa_grid: np.ndarray,
    phase: str,
    w1: float,
    w2: float,
) -> List[Tuple[float, float, str]]:
    """
    Build single-phase isotherm points at fixed T from pressure-root solves.

    Thermodynamic basis
    -------------------
    Helmholtz EOS uses (T, rho) as natural variables. For plotting at fixed (T, P),
    rho is obtained from phase-selected roots:
    - low-density root for vapor branch
    - high-density root for liquid branch
    This keeps traces in the physically intended single-phase basins.

    References
    ----------
    - Span, R. (2000). Multiparameter Equations of State.
    - Lemmon, E. W., Huber, M. L., & McLinden, M. O. (2018). NIST REFPROP
      Documentation, Version 10.0.

    Breadcrumb context
    ------------------
    Added during R515B isotherm-layer replacement (2026-03-03) to keep only
    phase-selected single-phase points in plotted isotherms.
    """
    out: List[Tuple[float, float, str]] = []
    interaction = default_interaction()
    for P_kPa in P_kPa_grid:
        rho = solve_rho_mass_for_P(
            T_K=T_K,
            P_target_kPa=float(P_kPa),
            comp1_json="r1234ze.json",
            comp2_json="r227ea.json",
            w1=w1,
            w2=w2,
            interaction=interaction,
            phase_hint=phase,
        )
        if rho is None or not np.isfinite(rho):
            continue
        st = compute_props("r1234ze.json", "r227ea.json", T_K, float(rho), w1, w2, interaction)
        out.append((float(st.p * 1.0e-2), float(st.h), phase))
    return out


def generate_refprop_style_isotherm(
    T_C: float = 70.0,
    w1: float = 0.911,
    P_min_bar: float = 1.0,
    P_max_bar: float = 100.0,
    n_vapor: int = 80,
    n_liquid: int = 80,
    eps_kPa: float = 1.0,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Generate REFPROP-style mixture isotherm data: vapor + connector + liquid.

    Thermodynamic basis
    -------------------
    Plot only single-phase segments plus a saturation-endpoint connector.
    No two-phase interior (metastable) density sweep points are generated.
    For mixtures, connector joins (h_bub, P_bub) to (h_dew, P_dew) at the same T.

    References
    ----------
    - Lemmon, E. W., Huber, M. L., & McLinden, M. O. (2018). NIST REFPROP
      Documentation, Version 10.0.
    - REFPROP-docs (NIST GitHub), notes on metastable two-phase plotting behavior.
    - EES REFPROP interface plotting note: straight line between saturated
      liquid and vapor points by default in property plots.

    Breadcrumb context
    ------------------
    Added during R515B isotherm-layer replacement (2026-03-03) to eliminate
    interior metastable loops from default isotherm visualization.
    """
    boundary_csv = ROOT / "diagnostics/plots/ph_layer2_true_vle_boundary_20260303.csv"
    endp = interpolate_endpoints_at_T(boundary_csv, T_C)
    T_K = endp["T_K"]
    w2 = 1.0 - w1

    P_dew_bar = endp["P_dew_kPa"] * 1.0e-2
    P_bub_bar = endp["P_bub_kPa"] * 1.0e-2

    vapor_hi_bar = max(P_min_bar, P_dew_bar - eps_kPa * 1.0e-2)
    liquid_lo_bar = min(P_max_bar, P_bub_bar + eps_kPa * 1.0e-2)

    vapor_points: List[Tuple[float, float, str]] = []
    liquid_points: List[Tuple[float, float, str]] = []
    if vapor_hi_bar > P_min_bar:
        P_v_kPa = np.logspace(np.log10(P_min_bar * 100.0), np.log10(vapor_hi_bar * 100.0), int(n_vapor))
        vapor_points = single_phase_branch_points(T_K=T_K, P_kPa_grid=P_v_kPa, phase="vapor", w1=w1, w2=w2)
    if P_max_bar > liquid_lo_bar:
        P_l_kPa = np.logspace(np.log10(liquid_lo_bar * 100.0), np.log10(P_max_bar * 100.0), int(n_liquid))
        liquid_points = single_phase_branch_points(T_K=T_K, P_kPa_grid=P_l_kPa, phase="liquid", w1=w1, w2=w2)

    rows = []
    for p_bar, h_kJkg, phase in vapor_points:
        rows.append({"P_bar": p_bar, "h_kJkg": h_kJkg, "phase_flag": phase})
    rows.append({"P_bar": P_bub_bar, "h_kJkg": endp["h_bub_kJkg"], "phase_flag": "two_phase_connector"})
    rows.append({"P_bar": P_dew_bar, "h_kJkg": endp["h_dew_kJkg"], "phase_flag": "two_phase_connector"})
    for p_bar, h_kJkg, phase in liquid_points:
        rows.append({"P_bar": p_bar, "h_kJkg": h_kJkg, "phase_flag": phase})
    return pd.DataFrame(rows), endp


def main() -> None:
    """
    CLI for REFPROP-style mixture isotherm generation.

    Thermodynamic basis
    -------------------
    Single-phase branches are traced with phase-selected roots at fixed T and P,
    and the two-phase interior is represented only by an endpoint connector.
    This avoids plotting metastable interior loops while preserving EOS state calls.

    References
    ----------
    - NIST REFPROP Documentation (v10, 2018) and REFPROP-docs notes on
      metastable two-phase plotting semantics.

    Breadcrumb context
    ------------------
    Added during R515B isotherm-layer replacement (2026-03-03).
    """
    p = argparse.ArgumentParser(description="Generate REFPROP-style P-h isotherm.")
    p.add_argument("--T-C", type=float, default=70.0)
    p.add_argument("--w1", type=float, default=0.911)
    p.add_argument("--Pmin-bar", type=float, default=1.0)
    p.add_argument("--Pmax-bar", type=float, default=100.0)
    p.add_argument("--nv", type=int, default=80)
    p.add_argument("--nl", type=int, default=80)
    p.add_argument("--eps-kpa", type=float, default=1.0)
    args = p.parse_args()

    df, endp = generate_refprop_style_isotherm(
        T_C=args.T_C,
        w1=args.w1,
        P_min_bar=args.Pmin_bar,
        P_max_bar=args.Pmax_bar,
        n_vapor=args.nv,
        n_liquid=args.nl,
        eps_kPa=args.eps_kpa,
    )

    out_dir = ROOT / "diagnostics/isotherms_refprop_style"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_path = out_dir / f"T_{int(round(args.T_C))}C_REFPROP_style_isotherm.png"
    csv_path = out_dir / f"T_{int(round(args.T_C))}C_REFPROP_style_isotherm.csv"
    df.to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(8.4, 5.9), dpi=180)
    dv = df[df["phase_flag"] == "vapor"]
    dc = df[df["phase_flag"] == "two_phase_connector"]
    dl = df[df["phase_flag"] == "liquid"]
    if len(dv):
        ax.plot(dv["h_kJkg"], dv["P_bar"], lw=1.3, color="#2e8b57", label="vapor")
    if len(dc) == 2:
        ax.plot(dc["h_kJkg"], dc["P_bar"], lw=1.3, color="#f4a259", label="two_phase_connector")
    if len(dl):
        ax.plot(dl["h_kJkg"], dl["P_bar"], lw=1.3, color="#7f3c8d", label="liquid")
    ax.set_yscale("log")
    ax.set_xlim(150.0, 500.0)
    ax.set_ylim(1.0, 100.0)
    ax.set_xlabel("h [kJ/kg]")
    ax.set_ylabel("P [bar]")
    ax.set_title(f"T={args.T_C:.1f} C REFPROP-style isotherm")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")

    print(f"P_bub({args.T_C:.1f}C)={endp['P_bub_kPa']*1e-2:.6f} bar, h_bub={endp['h_bub_kJkg']:.6f} kJ/kg")
    print(f"P_dew({args.T_C:.1f}C)={endp['P_dew_kPa']*1e-2:.6f} bar, h_dew={endp['h_dew_kJkg']:.6f} kJ/kg")
    print(f"n_vapor_points={len(dv)} n_liquid_points={len(dl)}")
    print(f"saved_figure={fig_path}")
    print(f"saved_csv={csv_path}")


if __name__ == "__main__":
    main()
