#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: Compute true binary VLE envelopes (bubble/dew) for a
fixed-overall-composition refrigerant blend using Helmholtz EOS mixture model.
Dependencies: numpy, scipy, matplotlib, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-02
# Assumptions:
# - Binary VLE solved from P equality and component chemical-potential equality.
# - Chemical potentials are evaluated numerically from total Helmholtz energy
#   A(T,V,n1,n2) using finite differences at fixed T,V,n_j.
# - Bubble curve fixes liquid composition x=z; dew curve fixes vapor composition y=z.
# - This module is isolated from compute_pressure_enthalpy workflow.
# TODO: Replace finite-difference chemical potentials with analytic expressions.
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


EPS_X = 1e-12


@dataclass
class MixState:
    """State container for one phase."""

    p_pa: float
    h_jmol: float
    g_jmol: float
    rho_mol: float
    rho_mass: float
    x1: float


def w1_to_x1(w1: float, mw1: float, mw2: float) -> float:
    """
    Purpose
    -------
    Convert component-1 mass fraction to mole fraction.

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


def _clip_x(x1: float) -> float:
    """Clip composition to open interval for log-stable transforms."""
    return float(min(max(x1, EPS_X), 1.0 - EPS_X))


def _mix_alpha_and_derivs(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> Tuple[float, float, float]:
    """
    Evaluate mixture alpha and reduced derivatives at fixed composition.

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]

    Outputs
    -------
    alpha_mix : float [unitless]
    alpha_tau_mix : float [unitless]
    ar_del_mix : float [unitless]
    """
    x1 = _clip_x(x1)
    x2 = 1.0 - x1

    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
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

    # Add ideal mixing entropy contribution in Helmholtz form.
    a0_mix = a0 + x1 * np.log(x1) + x2 * np.log(x2)
    return float(a0_mix + ar), float(a0_tau + ar_tau), float(ar_del)


def mix_state(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> MixState:
    """
    Compute pressure/enthalpy/Gibbs state at fixed composition.

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]

    Outputs
    -------
    MixState
      p [Pa], h [J/mol], g [J/mol], rho_mol [mol/m^3], rho_mass [kg/m^3], x1 [mol/mol]
    """
    x1 = _clip_x(x1)
    x2 = 1.0 - x1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = x1 * mw1 + x2 * mw2

    alpha, alpha_tau, ar_del = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1)
    # Recover tau by invert identity alpha_tau contribution requires reduced tau.
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k
    delta = rho_mol * vred

    z = 1.0 + delta * ar_del
    p_pa = rho_mol * R_u * t_k * z
    h_jmol = R_u * t_k * (1.0 + tau * alpha_tau + delta * ar_del)
    g_jmol = R_u * t_k * (1.0 + alpha + delta * ar_del)
    return MixState(
        p_pa=float(p_pa),
        h_jmol=float(h_jmol),
        g_jmol=float(g_jmol),
        rho_mol=float(rho_mol),
        rho_mass=float(rho_mol * mw_mix),
        x1=float(x1),
    )


def total_helmholtz_a(d1: Dict, d2: Dict, t_k: float, v_m3: float, n1_mol: float, n2_mol: float) -> float:
    """
    Evaluate total Helmholtz energy A(T,V,n1,n2).

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    v_m3 : float [m^3]
    n1_mol : float [mol]
    n2_mol : float [mol]

    Outputs
    -------
    A : float [J]
    """
    n1 = max(n1_mol, EPS_X)
    n2 = max(n2_mol, EPS_X)
    n = n1 + n2
    x1 = _clip_x(n1 / n)
    rho_mol = n / v_m3
    alpha, _, _ = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1)
    return float(n * R_u * t_k * alpha)


def chemical_potentials_fd(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> Tuple[float, float]:
    """
    Compute component chemical potentials by finite-difference of A(T,V,n).

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]

    Outputs
    -------
    mu1, mu2 : float [J/mol]
    """
    x1 = _clip_x(x1)
    x2 = 1.0 - x1
    n_tot = 1.0
    n1 = x1 * n_tot
    n2 = x2 * n_tot
    v_m3 = n_tot / rho_mol

    dn1 = max(1e-8, 1e-6 * n1)
    dn2 = max(1e-8, 1e-6 * n2)

    a_p1 = total_helmholtz_a(d1, d2, t_k, v_m3, n1 + dn1, n2)
    a_m1 = total_helmholtz_a(d1, d2, t_k, v_m3, max(EPS_X, n1 - dn1), n2)
    mu1 = (a_p1 - a_m1) / (2.0 * dn1)

    a_p2 = total_helmholtz_a(d1, d2, t_k, v_m3, n1, n2 + dn2)
    a_m2 = total_helmholtz_a(d1, d2, t_k, v_m3, n1, max(EPS_X, n2 - dn2))
    mu2 = (a_p2 - a_m2) / (2.0 * dn2)
    return float(mu1), float(mu2)


def _sigmoid(z: float) -> float:
    """Stable logistic transform."""
    if z >= 0:
        ez = np.exp(-z)
        return float(1.0 / (1.0 + ez))
    ez = np.exp(z)
    return float(ez / (1.0 + ez))


def solve_bubble_at_t(
    d1: Dict,
    d2: Dict,
    t_k: float,
    z1: float,
    rho_l0: float,
    rho_v0: float,
    y10: float,
) -> Dict:
    """
    Solve bubble state at fixed liquid composition x=z.

    Unknowns: rho_l, rho_v, y1.
    Equations: P_l=P_v, mu1_l=mu1_v, mu2_l=mu2_v.
    """
    z1 = _clip_x(z1)

    def res(u):
        rho_l = float(np.exp(u[0]))
        rho_v = float(np.exp(u[1]))
        y1 = _clip_x(_sigmoid(float(u[2])))
        st_l = mix_state(d1, d2, t_k, rho_l, z1)
        st_v = mix_state(d1, d2, t_k, rho_v, y1)
        mu1_l, mu2_l = chemical_potentials_fd(d1, d2, t_k, rho_l, z1)
        mu1_v, mu2_v = chemical_potentials_fd(d1, d2, t_k, rho_v, y1)
        r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
        r2 = (mu1_l - mu1_v) / (R_u * t_k)
        r3 = (mu2_l - mu2_v) / (R_u * t_k)
        return np.array([r1, r2, r3], dtype=float)

    u0 = np.array([np.log(max(rho_l0, 1e-9)), np.log(max(rho_v0, 1e-12)), np.log(y10 / (1.0 - y10))], dtype=float)
    sol = root(res, u0, method="hybr", tol=1e-10)

    rho_l = float(np.exp(sol.x[0]))
    rho_v = float(np.exp(sol.x[1]))
    y1 = _clip_x(_sigmoid(float(sol.x[2])))
    st_l = mix_state(d1, d2, t_k, rho_l, z1)
    st_v = mix_state(d1, d2, t_k, rho_v, y1)
    mu1_l, mu2_l = chemical_potentials_fd(d1, d2, t_k, rho_l, z1)
    mu1_v, mu2_v = chemical_potentials_fd(d1, d2, t_k, rho_v, y1)
    r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
    r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)
    ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * (1.0 + 1e-8))
    return {
        "status": "CONVERGED" if ok else "DIVERGED",
        "T_K": float(t_k),
        "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
        "rho_l_molm3": float(rho_l),
        "rho_v_molm3": float(rho_v),
        "x1_liq": float(z1),
        "y1_vap": float(y1),
        "h_l_Jmol": float(st_l.h_jmol),
        "h_v_Jmol": float(st_v.h_jmol),
        "r_P": float(r_p),
        "r_mu": float(r_mu),
        "iterations": int(sol.nfev),
        "notes": str(sol.message),
    }


def solve_dew_at_t(
    d1: Dict,
    d2: Dict,
    t_k: float,
    z1: float,
    rho_l0: float,
    rho_v0: float,
    x10: float,
) -> Dict:
    """
    Solve dew state at fixed vapor composition y=z.

    Unknowns: rho_l, rho_v, x1.
    Equations: P_l=P_v, mu1_l=mu1_v, mu2_l=mu2_v.
    """
    z1 = _clip_x(z1)

    def res(u):
        rho_l = float(np.exp(u[0]))
        rho_v = float(np.exp(u[1]))
        x1 = _clip_x(_sigmoid(float(u[2])))
        st_l = mix_state(d1, d2, t_k, rho_l, x1)
        st_v = mix_state(d1, d2, t_k, rho_v, z1)
        mu1_l, mu2_l = chemical_potentials_fd(d1, d2, t_k, rho_l, x1)
        mu1_v, mu2_v = chemical_potentials_fd(d1, d2, t_k, rho_v, z1)
        r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
        r2 = (mu1_l - mu1_v) / (R_u * t_k)
        r3 = (mu2_l - mu2_v) / (R_u * t_k)
        return np.array([r1, r2, r3], dtype=float)

    u0 = np.array([np.log(max(rho_l0, 1e-9)), np.log(max(rho_v0, 1e-12)), np.log(x10 / (1.0 - x10))], dtype=float)
    sol = root(res, u0, method="hybr", tol=1e-10)

    rho_l = float(np.exp(sol.x[0]))
    rho_v = float(np.exp(sol.x[1]))
    x1 = _clip_x(_sigmoid(float(sol.x[2])))
    st_l = mix_state(d1, d2, t_k, rho_l, x1)
    st_v = mix_state(d1, d2, t_k, rho_v, z1)
    mu1_l, mu2_l = chemical_potentials_fd(d1, d2, t_k, rho_l, x1)
    mu1_v, mu2_v = chemical_potentials_fd(d1, d2, t_k, rho_v, z1)
    r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
    r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)
    ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * (1.0 + 1e-8))
    return {
        "status": "CONVERGED" if ok else "DIVERGED",
        "T_K": float(t_k),
        "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
        "rho_l_molm3": float(rho_l),
        "rho_v_molm3": float(rho_v),
        "x1_liq": float(x1),
        "y1_vap": float(z1),
        "h_l_Jmol": float(st_l.h_jmol),
        "h_v_Jmol": float(st_v.h_jmol),
        "r_P": float(r_p),
        "r_mu": float(r_mu),
        "iterations": int(sol.nfev),
        "notes": str(sol.message),
    }


def run_true_vle_envelope(
    fluid1: str,
    fluid2: str,
    w1: float,
    t_vals: np.ndarray,
) -> Tuple[List[Dict], List[Dict], float]:
    """Run bubble/dew solves across temperature grid with continuation."""
    d1 = load_idaes_helmholtz_json(fluid1)
    d2 = load_idaes_helmholtz_json(fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = w1_to_x1(w1, mw1, mw2)

    rhoc1 = float(d1["basic"]["rhoc"]) / mw1
    rhoc2 = float(d2["basic"]["rhoc"]) / mw2
    rho_l_guess = 0.8 * (z1 * rhoc1 + (1.0 - z1) * rhoc2)
    rho_v_guess = 0.01 * rho_l_guess
    y1_guess = z1
    x1_guess = z1

    bubble_rows: List[Dict] = []
    dew_rows: List[Dict] = []
    for t_k in np.asarray(t_vals, dtype=float):
        b = solve_bubble_at_t(d1, d2, float(t_k), z1, rho_l_guess, rho_v_guess, y1_guess)
        bubble_rows.append(b)
        if b["status"] == "CONVERGED":
            rho_l_guess = b["rho_l_molm3"]
            rho_v_guess = b["rho_v_molm3"]
            y1_guess = b["y1_vap"]

        d = solve_dew_at_t(d1, d2, float(t_k), z1, rho_l_guess, rho_v_guess, x1_guess)
        dew_rows.append(d)
        if d["status"] == "CONVERGED":
            rho_l_guess = d["rho_l_molm3"]
            rho_v_guess = d["rho_v_molm3"]
            x1_guess = d["x1_liq"]
    return bubble_rows, dew_rows, z1


def save_csv(rows: List[Dict], path: str | Path) -> None:
    """Save VLE branch rows to CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def plot_envelope(
    bubble_rows: List[Dict],
    dew_rows: List[Dict],
    mw_mix: float,
    out_fig: str | Path,
) -> None:
    """Plot P-h envelope from converged bubble/dew states."""
    b = [r for r in bubble_rows if r["status"] == "CONVERGED"]
    d = [r for r in dew_rows if r["status"] == "CONVERGED"]

    hb = np.array([(r["h_l_Jmol"] / mw_mix) * 1e-3 for r in b], dtype=float)
    pb = np.array([r["P_Pa"] * 1e-5 for r in b], dtype=float)
    hd = np.array([(r["h_v_Jmol"] / mw_mix) * 1e-3 for r in d], dtype=float)
    pd = np.array([r["P_Pa"] * 1e-5 for r in d], dtype=float)

    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if len(b):
        ax.plot(hb, pb, lw=2.0, label="Bubble line (liq)")
    if len(d):
        ax.plot(hd, pd, lw=2.0, label="Dew line (vap)")
    if len(b) and len(d):
        m = min(len(b), len(d))
        ax.fill_betweenx(pb[:m], hb[:m], hd[:m], alpha=0.15, label="Two-phase region")
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("R515B true VLE envelope (mu-equality solve)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")


def _cli() -> None:
    """CLI for true VLE envelope run."""
    p = argparse.ArgumentParser(description="Compute binary true VLE envelope from chemical-potential equality.")
    p.add_argument("--fluid1", default="r1234ze")
    p.add_argument("--fluid2", default="r227ea")
    p.add_argument("--w1", required=True, type=float, help="mass fraction fluid1 [kg/kg]")
    p.add_argument("--Tmin", required=True, type=float, help="minimum temperature [K]")
    p.add_argument("--Tmax", required=True, type=float, help="maximum temperature [K]")
    p.add_argument("--n", default=80, type=int, help="temperature points")
    p.add_argument("--bubble-csv", default="verification/r515b_true_vle_bubble.csv")
    p.add_argument("--dew-csv", default="verification/r515b_true_vle_dew.csv")
    p.add_argument("--fig", default="verification/r515b_true_vle_envelope.png")
    p.add_argument("--metadata", default="verification/r515b_true_vle_metadata.json")
    args = p.parse_args()

    t_vals = np.linspace(args.Tmin, args.Tmax, args.n)
    bubble_rows, dew_rows, z1 = run_true_vle_envelope(args.fluid1, args.fluid2, args.w1, t_vals)
    save_csv(bubble_rows, args.bubble_csv)
    save_csv(dew_rows, args.dew_csv)

    d1 = load_idaes_helmholtz_json(args.fluid1)
    d2 = load_idaes_helmholtz_json(args.fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = z1 * mw1 + (1.0 - z1) * mw2
    plot_envelope(bubble_rows, dew_rows, mw_mix, args.fig)

    nb = int(sum(r["status"] == "CONVERGED" for r in bubble_rows))
    nd = int(sum(r["status"] == "CONVERGED" for r in dew_rows))
    meta = {
        "fluid1": args.fluid1,
        "fluid2": args.fluid2,
        "w1_kgkg": float(args.w1),
        "z1_molmol": float(z1),
        "Tmin_K": float(args.Tmin),
        "Tmax_K": float(args.Tmax),
        "n_points": int(args.n),
        "bubble_converged": nb,
        "bubble_failed": int(args.n - nb),
        "dew_converged": nd,
        "dew_failed": int(args.n - nd),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "P equality + mu1/mu2 equality, FD chemical potentials from A(T,V,n1,n2)",
    }
    mpath = Path(args.metadata)
    mpath.parent.mkdir(parents=True, exist_ok=True)
    with mpath.open("w") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved bubble CSV: {args.bubble_csv}")
    print(f"Saved dew CSV: {args.dew_csv}")
    print(f"Saved figure: {args.fig}")
    print(f"Saved metadata: {args.metadata}")


if __name__ == "__main__":
    _cli()
