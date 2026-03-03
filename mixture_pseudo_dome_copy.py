#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: Compute an approximate fixed-composition pseudo-saturation
P-h dome for a binary mixture using Helmholtz EOS helpers.
Dependencies: numpy, scipy, matplotlib, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-02
# Assumptions:
# - Fixed composition pseudo-dome approximation (not full mixture VLE).
# - Coexistence is approximated with P_l=P_v and g_mix,l=g_mix,v at fixed x.
# - Composition derivatives through reducing rules are neglected here.
# TODO: Replace with full mixture fugacity-equality flash formulation.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

from linear_model_codex import (
    BELL_2023_R1234ZE_R227EA,
    R_u,
    bell2023_Tred_vred,
    load_idaes_helmholtz_json,
    mixture_alpha0_alphar_derivs,
    mw_from_json,
)


@dataclass
class MixtureState:
    """Thermodynamic state for fixed-composition mixture pseudo-dome solve."""

    p_pa: float
    h_jmol: float
    g_jmol: float
    rho_mass: float
    rho_mol: float


def w1_to_x1(w1: float, mw1: float, mw2: float) -> float:
    """
    Convert mass fraction to mole fraction.

    Inputs
    ------
    w1 : float [kg/kg]
    mw1 : float [kg/mol]
    mw2 : float [kg/mol]

    Outputs
    -------
    x1 : float [mol/mol]
    """
    w2 = 1.0 - w1
    n1 = w1 / mw1
    n2 = w2 / mw2
    return float(n1 / (n1 + n2))


def _mixture_state_from_rhomol(
    d1: Dict,
    d2: Dict,
    x1: float,
    x2: float,
    t_k: float,
    rho_mol: float,
) -> MixtureState:
    """
    Evaluate fixed-composition mixture p/h/g at (T, rho_mol).

    Inputs
    ------
    d1, d2 : dict [unitless]
    x1, x2 : float [mol/mol]
    t_k : float [K]
    rho_mol : float [mol/m^3]

    Outputs
    -------
    MixtureState
      p_pa [Pa], h_jmol [J/mol], g_jmol [J/mol], rho_mass [kg/m^3], rho_mol [mol/m^3]
    """
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = x1 * mw1 + x2 * mw2

    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k
    delta = rho_mol * vred
    rho_red_mol = 1.0 / vred

    a0, a0_tau, ar, ar_tau, ar_del = mixture_alpha0_alphar_derivs(
        d1=d1,
        d2=d2,
        x1=x1,
        x2=x2,
        tau=tau,
        delta=delta,
        Tred=tred,
        rho_red_mol=rho_red_mol,
        pair_key="r1234ze|r227ea",
    )

    z = 1.0 + delta * ar_del
    p_pa = rho_mol * R_u * t_k * z
    h_rt = 1.0 + tau * (a0_tau + ar_tau) + delta * ar_del
    h_jmol = h_rt * R_u * t_k
    g_jmol = R_u * t_k * (1.0 + (a0 + ar) + delta * ar_del)
    rho_mass = rho_mol * mw_mix
    return MixtureState(p_pa=float(p_pa), h_jmol=float(h_jmol), g_jmol=float(g_jmol), rho_mass=float(rho_mass), rho_mol=float(rho_mol))


def _acceptance_metrics(st_l: MixtureState, st_v: MixtureState, t_k: float) -> Tuple[float, float, bool]:
    """Compute residual gates for pseudo-saturation acceptance."""
    r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
    r_mu = abs(st_l.g_jmol - st_v.g_jmol) / (R_u * t_k)
    rho_ok = st_l.rho_mol > st_v.rho_mol * (1.0 + 1e-8)
    return float(r_p), float(r_mu), bool(rho_ok)


def pseudo_saturation_point_at_t(
    fluid1: str,
    fluid2: str,
    w1: float,
    t_k: float,
    rho_l_guess: Optional[float] = None,
    rho_v_guess: Optional[float] = None,
) -> Dict:
    """
    Solve one fixed-composition pseudo-saturation point.

    Inputs
    ------
    fluid1, fluid2 : str [unitless]
      Binary pair names.
    w1 : float [kg/kg]
      Mass fraction of fluid1.
    t_k : float [K]
      Temperature.
    rho_l_guess, rho_v_guess : float | None [mol/m^3]
      Initial guesses.

    Outputs
    -------
    result : dict [mixed units]
      Includes status, rho_l/v [mol/m^3], p [Pa], h_l/v [J/mol], r_P, r_mu.

    Assumptions
    -----------
    - Fixed composition with no phase split in composition (pseudo approximation).

    Failure modes
    -------------
    - Returns non-CONVERGED status if solver fails strict gates.
    """
    d1 = load_idaes_helmholtz_json(fluid1)
    d2 = load_idaes_helmholtz_json(fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    x1 = w1_to_x1(w1, mw1, mw2)
    x2 = 1.0 - x1

    if rho_l_guess is None:
        rho_l_guess = 0.85 * (float(d1["basic"]["rhoc"]) / mw1 * x1 + float(d2["basic"]["rhoc"]) / mw2 * x2)
    if rho_v_guess is None:
        rho_v_guess = 0.02 * rho_l_guess

    rv0 = max(min(rho_v_guess, 0.95 * rho_l_guess), 1e-12)
    gap0 = max(rho_l_guess - rv0, 1e-9)
    y0 = np.array([np.log(rv0), np.log(gap0)], dtype=float)

    def residual(y):
        rho_v = float(np.exp(y[0]))
        rho_l = float(rho_v + np.exp(y[1]))
        st_l = _mixture_state_from_rhomol(d1, d2, x1, x2, t_k, rho_l)
        st_v = _mixture_state_from_rhomol(d1, d2, x1, x2, t_k, rho_v)
        r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
        r2 = (st_l.g_jmol - st_v.g_jmol) / (R_u * t_k)
        return np.array([r1, r2], dtype=float)

    sol = root(residual, y0, method="hybr", tol=1e-12)
    rho_v = float(np.exp(sol.x[0]))
    rho_l = float(rho_v + np.exp(sol.x[1]))
    st_l = _mixture_state_from_rhomol(d1, d2, x1, x2, t_k, rho_l)
    st_v = _mixture_state_from_rhomol(d1, d2, x1, x2, t_k, rho_v)
    r_p, r_mu, rho_ok = _acceptance_metrics(st_l, st_v, t_k)
    converged = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_ok)

    status = "CONVERGED" if converged else "DIVERGED"
    return {
        "status": status,
        "T_K": float(t_k),
        "rho_l_molm3": float(rho_l),
        "rho_v_molm3": float(rho_v),
        "rho_l_kgm3": float(st_l.rho_mass),
        "rho_v_kgm3": float(st_v.rho_mass),
        "p_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
        "h_l_Jmol": float(st_l.h_jmol),
        "h_v_Jmol": float(st_v.h_jmol),
        "r_P": float(r_p),
        "r_mu": float(r_mu),
        "iterations": int(sol.nfev),
        "notes": str(sol.message),
        "x1_molmol": float(x1),
        "w1_kgkg": float(w1),
    }


def compute_pseudo_dome(
    fluid1: str,
    fluid2: str,
    w1: float,
    t_vals: np.ndarray,
) -> List[Dict]:
    """Compute pseudo-dome rows with continuation in density guesses."""
    rows: List[Dict] = []
    rho_l_guess = None
    rho_v_guess = None
    for t_k in np.asarray(t_vals, dtype=float):
        row = pseudo_saturation_point_at_t(
            fluid1=fluid1,
            fluid2=fluid2,
            w1=w1,
            t_k=float(t_k),
            rho_l_guess=rho_l_guess,
            rho_v_guess=rho_v_guess,
        )
        rows.append(row)
        rho_l_guess = row["rho_l_molm3"]
        rho_v_guess = row["rho_v_molm3"]
    return rows


def save_run_csv(rows: List[Dict], path: str | Path) -> None:
    """Save pseudo-dome run rows to CSV."""
    fieldnames = [
        "T_K",
        "status",
        "rho_l_molm3",
        "rho_v_molm3",
        "rho_l_kgm3",
        "rho_v_kgm3",
        "p_Pa",
        "h_l_Jmol",
        "h_v_Jmol",
        "r_P",
        "r_mu",
        "iterations",
        "x1_molmol",
        "w1_kgkg",
        "notes",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def plot_ph(rows: List[Dict], out_clean: str | Path, out_fail: str | Path, fluid1: str, fluid2: str) -> None:
    """Plot clean and diagnostic pseudo-dome P-h curves."""
    d1 = load_idaes_helmholtz_json(fluid1)
    d2 = load_idaes_helmholtz_json(fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    w1 = float(rows[0]["w1_kgkg"])
    x1 = w1_to_x1(w1, mw1, mw2)
    mw_mix = x1 * mw1 + (1.0 - x1) * mw2

    conv = [r for r in rows if r["status"] == "CONVERGED"]
    fail = [r for r in rows if r["status"] != "CONVERGED"]

    def hk(row, key):
        return (float(row[key]) / mw_mix) * 1e-3

    if conv:
        p_bar = np.array([float(r["p_Pa"]) * 1e-5 for r in conv], dtype=float)
        hl = np.array([hk(r, "h_l_Jmol") for r in conv], dtype=float)
        hv = np.array([hk(r, "h_v_Jmol") for r in conv], dtype=float)
    else:
        p_bar = np.array([], dtype=float)
        hl = np.array([], dtype=float)
        hv = np.array([], dtype=float)

    # Clean plot.
    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if len(p_bar):
        ax.plot(hl, p_bar, lw=2.0, label="Pseudo sat. liquid")
        ax.plot(hv, p_bar, lw=2.0, label="Pseudo sat. vapor")
        ax.fill_betweenx(p_bar, hl, hv, alpha=0.15, label="Pseudo two-phase region")
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("R515B pseudo-dome (fixed composition)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_clean, bbox_inches="tight")

    # Plot with failed points.
    fig2, ax2 = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if len(p_bar):
        ax2.plot(hl, p_bar, lw=2.0, label="Pseudo sat. liquid (CONVERGED)")
        ax2.plot(hv, p_bar, lw=2.0, label="Pseudo sat. vapor (CONVERGED)")
    if fail:
        p_fail = np.array([float(r["p_Pa"]) * 1e-5 for r in fail], dtype=float)
        hl_fail = np.array([hk(r, "h_l_Jmol") for r in fail], dtype=float)
        hv_fail = np.array([hk(r, "h_v_Jmol") for r in fail], dtype=float)
        ax2.plot(hl_fail, p_fail, "rx", ms=5, label="Failed liquid estimate")
        ax2.plot(hv_fail, p_fail, "rx", ms=5, label="Failed vapor estimate")
    ax2.set_yscale("log")
    ax2.set_xlabel("Enthalpy [kJ/kg]")
    ax2.set_ylabel("Pressure [bar]")
    ax2.set_title("R515B pseudo-dome (with failures)")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(out_fail, bbox_inches="tight")


def _cli() -> None:
    """CLI entry point for approximate fixed-composition pseudo-dome run."""
    p = argparse.ArgumentParser(description="Compute fixed-composition pseudo-dome for a binary mixture.")
    p.add_argument("--fluid1", default="r1234ze", help="fluid1 stem")
    p.add_argument("--fluid2", default="r227ea", help="fluid2 stem")
    p.add_argument("--w1", required=True, type=float, help="mass fraction of fluid1 [kg/kg]")
    p.add_argument("--Tmin", required=True, type=float, help="minimum T [K]")
    p.add_argument("--Tmax", required=True, type=float, help="maximum T [K]")
    p.add_argument("--n", default=140, type=int, help="number of temperature points")
    p.add_argument("--out", default="verification/r515b_pseudo_dome.csv", help="output CSV")
    p.add_argument("--fig-clean", default="verification/r515b_pseudo_ph_dome_clean.png", help="clean plot")
    p.add_argument("--fig-fail", default="verification/r515b_pseudo_ph_dome_with_failures.png", help="plot with failures")
    p.add_argument("--metadata", default="verification/r515b_pseudo_dome_metadata.json", help="output metadata JSON")
    args = p.parse_args()

    t_vals = np.linspace(args.Tmin, args.Tmax, args.n)
    rows = compute_pseudo_dome(args.fluid1, args.fluid2, args.w1, t_vals)
    save_run_csv(rows, args.out)
    plot_ph(rows, args.fig_clean, args.fig_fail, args.fluid1, args.fluid2)

    n_converged = int(sum(r["status"] == "CONVERGED" for r in rows))
    meta = {
        "fluid1": args.fluid1,
        "fluid2": args.fluid2,
        "w1_kgkg": float(args.w1),
        "Tmin_K": float(args.Tmin),
        "Tmax_K": float(args.Tmax),
        "n_points": int(args.n),
        "n_converged": n_converged,
        "n_failed": int(len(rows) - n_converged),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "approximation": "fixed-composition pseudo-dome (not full mixture VLE)",
    }
    mpath = Path(args.metadata)
    mpath.parent.mkdir(parents=True, exist_ok=True)
    with mpath.open("w") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved CSV: {args.out}")
    print(f"Saved figure: {args.fig_clean}")
    print(f"Saved failure figure: {args.fig_fail}")
    print(f"Saved metadata: {args.metadata}")


if __name__ == "__main__":
    _cli()
