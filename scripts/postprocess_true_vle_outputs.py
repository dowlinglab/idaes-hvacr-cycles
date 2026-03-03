#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: Post-process true-VLE bubble/dew CSV files into clean
converged outputs and failure diagnostics CSV/plots for P-h envelope review.
Dependencies: argparse, csv, pathlib, numpy, matplotlib, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from linear_model_codex import load_idaes_helmholtz_json, mw_from_json, x1_from_w1


def _read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))


def _to_float(v: str) -> float:
    return float(v.strip())


def _write_csv(path: Path, rows: List[Dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _enthalpy_kjkg(h_jmol: float, mw_mix: float) -> float:
    return (h_jmol / mw_mix) * 1e-3


def main() -> None:
    p = argparse.ArgumentParser(description="Post-process true-VLE bubble/dew outputs.")
    p.add_argument("--fluid1", default="r1234ze")
    p.add_argument("--fluid2", default="r227ea")
    p.add_argument("--w1", required=True, type=float, help="mass fraction fluid1 [kg/kg]")
    p.add_argument("--bubble-csv", required=True, help="input bubble CSV")
    p.add_argument("--dew-csv", required=True, help="input dew CSV")
    p.add_argument("--diag-csv", default="verification/r515b_true_vle_diagnostics.csv", help="output diagnostics CSV")
    p.add_argument("--conv-csv", default="verification/r515b_true_vle_converged_only.csv", help="output converged-only CSV")
    p.add_argument("--fig-clean", default="verification/r515b_true_vle_envelope_clean.png", help="output clean envelope plot")
    p.add_argument(
        "--fig-fail",
        default="verification/r515b_true_vle_envelope_with_failures.png",
        help="output diagnostics envelope plot with failed points",
    )
    args = p.parse_args()

    bubble_rows = _read_rows(Path(args.bubble_csv))
    dew_rows = _read_rows(Path(args.dew_csv))

    d1 = load_idaes_helmholtz_json(args.fluid1)
    d2 = load_idaes_helmholtz_json(args.fluid2)
    x1 = x1_from_w1(d1, d2, args.w1)
    x2 = 1.0 - x1
    mw_mix = x1 * mw_from_json(d1) + x2 * mw_from_json(d2)

    diagnostics = []
    converged = []

    def append_rows(rows: List[Dict[str, str]], branch: str, h_key: str) -> None:
        for r in rows:
            status = r["status"]
            t_k = _to_float(r["T_K"])
            p_pa = _to_float(r["P_Pa"])
            h_jmol = _to_float(r[h_key])
            rec = {
                "branch": branch,
                "status": status,
                "T_K": t_k,
                "P_Pa": p_pa,
                "P_bar": p_pa * 1e-5,
                "h_Jmol": h_jmol,
                "h_kJkg": _enthalpy_kjkg(h_jmol, mw_mix),
                "rho_l_molm3": _to_float(r["rho_l_molm3"]),
                "rho_v_molm3": _to_float(r["rho_v_molm3"]),
                "x1_liq": _to_float(r["x1_liq"]),
                "y1_vap": _to_float(r["y1_vap"]),
                "r_P": _to_float(r["r_P"]),
                "r_mu": _to_float(r["r_mu"]),
                "iterations": int(float(r["iterations"])),
                "notes": r["notes"].replace("\n", " ").strip(),
            }
            diagnostics.append(rec)
            if status == "CONVERGED":
                converged.append(rec)

    append_rows(bubble_rows, "bubble_liq", "h_l_Jmol")
    append_rows(dew_rows, "dew_vap", "h_v_Jmol")

    diag_fields = [
        "branch",
        "status",
        "T_K",
        "P_Pa",
        "P_bar",
        "h_Jmol",
        "h_kJkg",
        "rho_l_molm3",
        "rho_v_molm3",
        "x1_liq",
        "y1_vap",
        "r_P",
        "r_mu",
        "iterations",
        "notes",
    ]
    _write_csv(Path(args.diag_csv), diagnostics, diag_fields)
    _write_csv(Path(args.conv_csv), converged, diag_fields)

    # Plot clean envelope
    bconv = [r for r in diagnostics if r["branch"] == "bubble_liq" and r["status"] == "CONVERGED"]
    dconv = [r for r in diagnostics if r["branch"] == "dew_vap" and r["status"] == "CONVERGED"]
    bconv = sorted(bconv, key=lambda r: r["T_K"])
    dconv = sorted(dconv, key=lambda r: r["T_K"])

    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if bconv:
        hb = np.array([r["h_kJkg"] for r in bconv], dtype=float)
        pb = np.array([r["P_bar"] for r in bconv], dtype=float)
        ax.plot(hb, pb, lw=2.0, label="Bubble line (CONVERGED)")
    if dconv:
        hd = np.array([r["h_kJkg"] for r in dconv], dtype=float)
        pd = np.array([r["P_bar"] for r in dconv], dtype=float)
        ax.plot(hd, pd, lw=2.0, label="Dew line (CONVERGED)")
    if bconv and dconv:
        m = min(len(bconv), len(dconv))
        ax.fill_betweenx(
            np.array([r["P_bar"] for r in bconv[:m]], dtype=float),
            np.array([r["h_kJkg"] for r in bconv[:m]], dtype=float),
            np.array([r["h_kJkg"] for r in dconv[:m]], dtype=float),
            alpha=0.15,
            label="Two-phase region (CONVERGED)",
        )
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    Path(args.fig_clean).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.fig_clean, bbox_inches="tight")

    # Plot diagnostics with failed points
    fig2, ax2 = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if bconv:
        hb = np.array([r["h_kJkg"] for r in bconv], dtype=float)
        pb = np.array([r["P_bar"] for r in bconv], dtype=float)
        ax2.plot(hb, pb, lw=2.0, label="Bubble line (CONVERGED)")
    if dconv:
        hd = np.array([r["h_kJkg"] for r in dconv], dtype=float)
        pd = np.array([r["P_bar"] for r in dconv], dtype=float)
        ax2.plot(hd, pd, lw=2.0, label="Dew line (CONVERGED)")

    failed = [r for r in diagnostics if r["status"] != "CONVERGED"]
    if failed:
        hf = np.array([r["h_kJkg"] for r in failed], dtype=float)
        pf = np.array([r["P_bar"] for r in failed], dtype=float)
        ax2.plot(hf, pf, "rx", ms=5, label="Failed points")

    ax2.set_yscale("log")
    ax2.set_xlabel("Enthalpy [kJ/kg]")
    ax2.set_ylabel("Pressure [bar]")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    Path(args.fig_fail).parent.mkdir(parents=True, exist_ok=True)
    fig2.savefig(args.fig_fail, bbox_inches="tight")

    n_b = len([r for r in bubble_rows if r["status"] == "CONVERGED"])
    n_d = len([r for r in dew_rows if r["status"] == "CONVERGED"])
    print(f"bubble_converged={n_b}/{len(bubble_rows)}")
    print(f"dew_converged={n_d}/{len(dew_rows)}")
    print(f"Saved diagnostics CSV: {args.diag_csv}")
    print(f"Saved converged-only CSV: {args.conv_csv}")
    print(f"Saved clean figure: {args.fig_clean}")
    print(f"Saved diagnostics figure: {args.fig_fail}")


if __name__ == "__main__":
    main()
