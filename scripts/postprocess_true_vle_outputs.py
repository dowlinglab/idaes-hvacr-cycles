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

Version: v0.2.0
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
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


def detect_crossover_indices(h_liq_kjkg: np.ndarray, h_vap_kjkg: np.ndarray) -> List[int]:
    """
    Purpose
    -------
    Detect temperature-indexed crossover locations where h_liq > h_vap.

    Inputs
    ------
    h_liq_kjkg : np.ndarray [kJ/kg]
    h_vap_kjkg : np.ndarray [kJ/kg]

    Outputs
    -------
    List[int] [index]
        Indices i where h_liq[i] > h_vap[i].

    Assumptions
    -----------
    - Arrays are same length and correspond to aligned temperature rows.

    Failure modes
    -------------
    - Raises ValueError if lengths differ.

    References
    ----------
    - Zeotropic branch-order diagnostics for P-h visualization.

    Notes on numerical stability
    ----------------------------
    - Uses direct elementwise comparison without differencing.
    """
    if h_liq_kjkg.shape != h_vap_kjkg.shape:
        raise ValueError("h_liq_kjkg and h_vap_kjkg must have identical shape")
    return np.where(h_liq_kjkg > h_vap_kjkg)[0].astype(int).tolist()


def _segmented_line_data(x: np.ndarray, y: np.ndarray, crossover_indices: List[int]) -> tuple[np.ndarray, np.ndarray]:
    """Insert NaN line breaks at crossover indices while preserving scatter points."""
    xs = x.copy()
    ys = y.copy()
    for idx in crossover_indices:
        if 0 <= idx < len(xs):
            xs[idx] = np.nan
            ys[idx] = np.nan
    return xs, ys


def get_crossovers(dome_csv_path: str | Path) -> List[Dict[str, float]]:
    """
    Purpose
    -------
    Return crossover diagnostics from a postprocessed converged-only dome CSV.

    Inputs
    ------
    dome_csv_path : str | Path [filesystem path]

    Outputs
    -------
    List[Dict[str, float]] [mixed units]
        Crossover dictionaries with T_K, P_bubble_kPa, P_dew_kPa, h_liq_kJkg,
        h_vap_kJkg, x_liq, y_vap, index.

    Assumptions
    -----------
    - Input CSV is produced by this script and contains bubble/dew rows.

    Failure modes
    -------------
    - Raises ValueError when bubble/dew branches are misaligned in length.

    References
    ----------
    - Internal plotting diagnostics API.

    Notes on numerical stability
    ----------------------------
    - Not applicable.
    """
    rows = _read_rows(Path(dome_csv_path))
    b = sorted([r for r in rows if r["branch"] == "bubble_liq" and r["status"] == "CONVERGED"], key=lambda r: _to_float(r["T_K"]))
    d = sorted([r for r in rows if r["branch"] == "dew_vap" and r["status"] == "CONVERGED"], key=lambda r: _to_float(r["T_K"]))
    if len(b) != len(d):
        raise ValueError("Bubble/dew branches are not aligned for crossover extraction")

    h_liq = np.array([_to_float(r["h_kJkg"]) for r in b], dtype=float)
    h_vap = np.array([_to_float(r["h_kJkg"]) for r in d], dtype=float)
    idxs = detect_crossover_indices(h_liq, h_vap)

    out: List[Dict[str, float]] = []
    for i in idxs:
        out.append(
            {
                "index": int(i),
                "T_K": float(_to_float(b[i]["T_K"])),
                "P_bubble_kPa": float(_to_float(b[i]["P_Pa"]) * 1e-3),
                "P_dew_kPa": float(_to_float(d[i]["P_Pa"]) * 1e-3),
                "h_liq_kJkg": float(_to_float(b[i]["h_kJkg"])),
                "h_vap_kJkg": float(_to_float(d[i]["h_kJkg"])),
                "x_liq": float(_to_float(b[i]["x1_liq"])),
                "y_vap": float(_to_float(d[i]["y1_vap"])),
            }
        )
    return out


def _build_continuous_envelope_proxy(
    t_k: np.ndarray,
    p_b_bar: np.ndarray,
    p_d_bar: np.ndarray,
    h_b_kjkg: np.ndarray,
    h_d_kjkg: np.ndarray,
    beta_grid: np.ndarray,
) -> Dict[str, np.ndarray]:
    """
    Purpose
    -------
    Build a diagnostic continuous envelope proxy by beta-interpolation between
    bubble and dew branch points on the same temperature grid.

    Inputs
    ------
    t_k : np.ndarray [K]
    p_b_bar, p_d_bar : np.ndarray [bar]
    h_b_kjkg, h_d_kjkg : np.ndarray [kJ/kg]
    beta_grid : np.ndarray [unitless]

    Outputs
    -------
    Dict[str, np.ndarray] [mixed units]

    Assumptions
    -----------
    - This is a visualization-only proxy and not a TP-flash solve.

    Failure modes
    -------------
    - Raises ValueError for shape mismatch.

    References
    ----------
    - Mixture envelope diagnostics for reviewer interpretation.

    Notes on numerical stability
    ----------------------------
    - Uses linear interpolation only.
    """
    if not (t_k.shape == p_b_bar.shape == p_d_bar.shape == h_b_kjkg.shape == h_d_kjkg.shape):
        raise ValueError("Input arrays must be shape-aligned")

    beta = np.asarray(beta_grid, dtype=float)
    h_rows = []
    p_rows = []
    t_rows = []
    beta_rows = []
    for b in beta:
        h_rows.append((1.0 - b) * h_b_kjkg + b * h_d_kjkg)
        p_rows.append((1.0 - b) * p_b_bar + b * p_d_bar)
        t_rows.append(t_k)
        beta_rows.append(np.full_like(t_k, b, dtype=float))

    return {
        "h_kJkg": np.concatenate(h_rows),
        "p_bar": np.concatenate(p_rows),
        "T_K": np.concatenate(t_rows),
        "beta": np.concatenate(beta_rows),
    }


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
    p.add_argument("--fig-pt", default="verification/r515b_true_vle_pt_overlay.png", help="P-T overlay figure")
    p.add_argument("--fig-ht", default="verification/r515b_true_vle_ht_overlay.png", help="h-T overlay figure")
    p.add_argument(
        "--fig-ph-colored",
        default="verification/r515b_true_vle_ph_colored_by_comp.png",
        help="P-h plot with composition colorization",
    )
    p.add_argument(
        "--cross-csv",
        default="",
        help="optional crossover diagnostic CSV path (default auto-dated in verification/)",
    )
    p.add_argument("--strict", action="store_true", help="return non-zero if any crossover is detected")
    p.add_argument(
        "--continuous-envelope",
        action="store_true",
        help="add diagnostic beta-interpolated continuous envelope proxy to P-h colorized plot",
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

    bconv = [r for r in diagnostics if r["branch"] == "bubble_liq" and r["status"] == "CONVERGED"]
    dconv = [r for r in diagnostics if r["branch"] == "dew_vap" and r["status"] == "CONVERGED"]
    bconv = sorted(bconv, key=lambda r: r["T_K"])
    dconv = sorted(dconv, key=lambda r: r["T_K"])

    if len(bconv) != len(dconv):
        raise ValueError("Bubble and dew converged rows must be length-aligned for crossover diagnostics")

    t_b = np.array([r["T_K"] for r in bconv], dtype=float)
    h_b = np.array([r["h_kJkg"] for r in bconv], dtype=float)
    p_b = np.array([r["P_bar"] for r in bconv], dtype=float)
    x_b = np.array([r["x1_liq"] for r in bconv], dtype=float)

    t_d = np.array([r["T_K"] for r in dconv], dtype=float)
    h_d = np.array([r["h_kJkg"] for r in dconv], dtype=float)
    p_d = np.array([r["P_bar"] for r in dconv], dtype=float)
    y_d = np.array([r["y1_vap"] for r in dconv], dtype=float)

    if not np.allclose(t_b, t_d, rtol=0.0, atol=1e-10):
        raise ValueError("Bubble/dew T grids must match for paired diagnostics")

    crossover_indices = detect_crossover_indices(h_b, h_d)
    cross_count = len(crossover_indices)

    if args.cross_csv.strip():
        cross_csv_path = Path(args.cross_csv)
    else:
        stamp = datetime.now().strftime("%Y%m%d")
        cross_csv_path = Path(f"verification/r515b_dome_crossovers_{stamp}.csv")

    cross_rows = []
    for i in crossover_indices:
        cross_rows.append(
            {
                "index": int(i),
                "T[K]": float(t_b[i]),
                "P_bubble[kPa]": float(p_b[i] * 100.0),
                "P_dew[kPa]": float(p_d[i] * 100.0),
                "h_liq[kJ/kg]": float(h_b[i]),
                "h_vap[kJ/kg]": float(h_d[i]),
                "x_liq": float(x_b[i]),
                "y_vap": float(y_d[i]),
            }
        )

    cross_fields = ["index", "T[K]", "P_bubble[kPa]", "P_dew[kPa]", "h_liq[kJ/kg]", "h_vap[kJ/kg]", "x_liq", "y_vap"]
    _write_csv(cross_csv_path, cross_rows, cross_fields)

    hb_line, pb_line = _segmented_line_data(h_b, p_b, crossover_indices)
    hd_line, pd_line = _segmented_line_data(h_d, p_d, crossover_indices)

    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    ax.plot(hb_line, pb_line, lw=2.0, color="#1f77b4", label="Bubble line")
    ax.plot(hd_line, pd_line, lw=2.0, color="#ff7f0e", label="Dew line")
    ax.scatter(h_b, p_b, s=9, color="#1f77b4", alpha=0.6)
    ax.scatter(h_d, p_d, s=9, color="#ff7f0e", alpha=0.6)

    red_label_used = False
    for i in range(len(t_b)):
        if h_b[i] <= h_d[i]:
            ax.plot([h_b[i], h_d[i]], [p_b[i], p_d[i]], ls="--", lw=0.4, color="gray", alpha=0.28)
        else:
            lbl = "h ordering inverted" if not red_label_used else None
            ax.plot([h_b[i], h_d[i]], [p_b[i], p_d[i]], ls="--", lw=0.7, color="red", alpha=0.7, label=lbl)
            red_label_used = True

    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    ax.text(
        0.02,
        0.98,
        f"Crossovers detected: {cross_count} points\\nSee {cross_csv_path}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#999", "alpha": 0.9},
    )
    fig.tight_layout()
    Path(args.fig_clean).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.fig_clean, bbox_inches="tight")

    fig2, ax2 = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    ax2.plot(hb_line, pb_line, lw=2.0, color="#1f77b4", label="Bubble line")
    ax2.plot(hd_line, pd_line, lw=2.0, color="#ff7f0e", label="Dew line")
    ax2.scatter(h_b, p_b, s=10, color="#1f77b4", alpha=0.6)
    ax2.scatter(h_d, p_d, s=10, color="#ff7f0e", alpha=0.6)

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
    ax2.text(
        0.02,
        0.98,
        f"Crossovers detected: {cross_count} points\\nSee {cross_csv_path}",
        transform=ax2.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#999", "alpha": 0.9},
    )
    fig2.tight_layout()
    Path(args.fig_fail).parent.mkdir(parents=True, exist_ok=True)
    fig2.savefig(args.fig_fail, bbox_inches="tight")

    # P-T overlay
    fig3, ax3 = plt.subplots(figsize=(7.0, 5.0), dpi=160)
    ax3.plot(t_b, p_b, color="#1f77b4", lw=2.0, label="Bubble P(T)")
    ax3.plot(t_d, p_d, color="#ff7f0e", lw=2.0, label="Dew P(T)")
    ax3.set_xlabel("Temperature [K]")
    ax3.set_ylabel("Pressure [bar]")
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=8)
    fig3.tight_layout()
    Path(args.fig_pt).parent.mkdir(parents=True, exist_ok=True)
    fig3.savefig(args.fig_pt, bbox_inches="tight")

    # h-T overlay with crossover highlight
    fig4, ax4 = plt.subplots(figsize=(7.0, 5.0), dpi=160)
    ax4.plot(t_b, h_b, color="#1f77b4", lw=2.0, label="h_liq(T)")
    ax4.plot(t_d, h_d, color="#ff7f0e", lw=2.0, label="h_vap(T)")
    if crossover_indices:
        tc = t_b[np.array(crossover_indices, dtype=int)]
        ax4.scatter(tc, h_b[np.array(crossover_indices, dtype=int)], color="red", s=18, label="Crossover points")
    ax4.set_xlabel("Temperature [K]")
    ax4.set_ylabel("Enthalpy [kJ/kg]")
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=8)
    fig4.tight_layout()
    Path(args.fig_ht).parent.mkdir(parents=True, exist_ok=True)
    fig4.savefig(args.fig_ht, bbox_inches="tight")

    # P-h colored by composition
    fig5, ax5 = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    sc1 = ax5.scatter(h_b, p_b, c=x_b, cmap="viridis", s=18, label="Bubble (color=x_liq)")
    sc2 = ax5.scatter(h_d, p_d, c=y_d, cmap="plasma", s=18, marker="s", label="Dew (color=y_vap)")
    ax5.plot(hb_line, pb_line, lw=1.2, color="#1f77b4", alpha=0.7)
    ax5.plot(hd_line, pd_line, lw=1.2, color="#ff7f0e", alpha=0.7)
    ax5.set_yscale("log")
    ax5.set_xlabel("Enthalpy [kJ/kg]")
    ax5.set_ylabel("Pressure [bar]")
    ax5.grid(True, which="both", alpha=0.3)
    ax5.legend(fontsize=8)
    cbar1 = fig5.colorbar(sc1, ax=ax5, pad=0.01)
    cbar1.set_label("x_liq [mol/mol]")
    cbar2 = fig5.colorbar(sc2, ax=ax5, pad=0.08)
    cbar2.set_label("y_vap [mol/mol]")

    if args.continuous_envelope:
        beta_grid = np.linspace(0.0, 1.0, 11)
        proxy = _build_continuous_envelope_proxy(t_b, p_b, p_d, h_b, h_d, beta_grid)
        ax5.scatter(
            proxy["h_kJkg"],
            proxy["p_bar"],
            c=proxy["beta"],
            cmap="cividis",
            s=6,
            alpha=0.35,
            label="Continuous envelope proxy (beta interpolation)",
        )
        ax5.text(
            0.02,
            0.02,
            "continuous-envelope uses interpolation proxy\\n(not TP-flash)",
            transform=ax5.transAxes,
            fontsize=7,
            ha="left",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "#999", "alpha": 0.85},
        )
    fig5.tight_layout()
    Path(args.fig_ph_colored).parent.mkdir(parents=True, exist_ok=True)
    fig5.savefig(args.fig_ph_colored, bbox_inches="tight")

    n_b = len([r for r in bubble_rows if r["status"] == "CONVERGED"])
    n_d = len([r for r in dew_rows if r["status"] == "CONVERGED"])
    print(f"bubble_converged={n_b}/{len(bubble_rows)}")
    print(f"dew_converged={n_d}/{len(dew_rows)}")
    print(f"crossovers={cross_count}")
    print(f"Saved diagnostics CSV: {args.diag_csv}")
    print(f"Saved converged-only CSV: {args.conv_csv}")
    print(f"Saved clean figure: {args.fig_clean}")
    print(f"Saved diagnostics figure: {args.fig_fail}")
    print(f"Saved P-T overlay: {args.fig_pt}")
    print(f"Saved h-T overlay: {args.fig_ht}")
    print(f"Saved P-h colored plot: {args.fig_ph_colored}")
    print(f"Saved crossover CSV: {cross_csv_path}")

    if cross_count > 0:
        print(
            f"WARNING: Detected {cross_count} crossover points where h_liq > h_vap. "
            f"Bubble and dew branches are not conjugate pairs; P-h dome may cross. "
            f"See {cross_csv_path} and PROJECT_CONTEXT.md."
        )
        if args.strict:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
