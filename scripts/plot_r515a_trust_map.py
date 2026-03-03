#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-03
Purpose of file: Build a temperature trust map for R515A p-H envelope results
using solver convergence diagnostics and available external validation.
Dependencies: csv, pathlib, numpy, matplotlib
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-03
# Assumptions:
# - Uses existing diagnostics from verification/r515a_true_vle_diagnostics_20260303.csv.
# - Uses available external validation at T=277.6 K from
#   verification/r515a_reference_suite_validation.csv.
# - Trust classes are heuristic and are intended for review triage, not publication.
# TODO: Replace heuristic trust classes with multi-point external validation thresholds.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_diag(path: Path) -> dict[float, dict[str, dict[str, float | str]]]:
    """
    Purpose
    -------
    Load branch diagnostics grouped by temperature.

    Inputs
    ------
    path : Path [file path]

    Outputs
    -------
    dict [K -> {branch -> row fields}]

    Assumptions
    -----------
    CSV schema follows scripts/postprocess_true_vle_outputs.py output.

    Failure modes
    -------------
    Raises on missing file or malformed numeric values.

    References
    ----------
    Project diagnostic CSV schema.

    Numerical stability notes
    -------------------------
    Not applicable.
    """
    out: dict[float, dict[str, dict[str, float | str]]] = {}
    with path.open() as f:
        for r in csv.DictReader(f):
            t = float(r["T_K"])
            branch = str(r["branch"])
            out.setdefault(t, {})[branch] = {
                "status": str(r["status"]),
                "P_bar": float(r["P_bar"]),
                "h_kJkg": float(r["h_kJkg"]),
                "r_P": float(r["r_P"]),
                "r_mu": float(r["r_mu"]),
            }
    return out


def read_validation_point(path: Path) -> tuple[float, float, float]:
    """
    Purpose
    -------
    Extract the R515A external validation temperature and relative errors.

    Inputs
    ------
    path : Path [file path]

    Outputs
    -------
    tuple [T_K, pressure_error_abs_pct, latent_heat_error_abs_pct]

    Assumptions
    -----------
    CSV contains R515A bubble/dew pressure_kPa and hfg_kJkg entries.

    Failure modes
    -------------
    Raises ValueError if required rows are not found.

    References
    ----------
    verification/r515a_reference_suite_validation.csv

    Numerical stability notes
    -------------------------
    Not applicable.
    """
    rows = list(csv.DictReader(path.open()))
    p_rows = [r for r in rows if r["case"].startswith("R515A_") and r["property"] == "pressure_kPa"]
    h_rows = [r for r in rows if r["case"].startswith("R515A_") and r["property"] == "hfg_kJkg"]
    if not p_rows or not h_rows:
        raise ValueError("Missing R515A reference-suite rows")

    p_err = float(np.mean([abs(float(r["rel_err"])) for r in p_rows]) * 100.0)
    h_err = float(np.mean([abs(float(r["rel_err"])) for r in h_rows]) * 100.0)
    return 277.6, p_err, h_err


def classify_trust(
    bubble_status: str,
    dew_status: str,
    bubble_rp: float,
    bubble_rmu: float,
    dew_rp: float,
    dew_rmu: float,
    ext_checked: bool,
    ext_ok: bool,
) -> str:
    """
    Purpose
    -------
    Assign trust class for one temperature point.

    Inputs
    ------
    bubble_status, dew_status : str [unitless]
    bubble_rp, bubble_rmu, dew_rp, dew_rmu : float [dimensionless]
    ext_checked : bool [unitless]
    ext_ok : bool [unitless]

    Outputs
    -------
    str [unitless]
      One of {'GREEN','YELLOW','RED'}.

    Assumptions
    -----------
    Residuals are normalized as in solver diagnostics.

    Failure modes
    -------------
    None; deterministic classification.

    References
    ----------
    Project acceptance gates and user-requested trust map categories.

    Numerical stability notes
    -------------------------
    Threshold comparisons only.
    """
    both_conv = bubble_status == "CONVERGED" and dew_status == "CONVERGED"
    one_conv = (bubble_status == "CONVERGED") ^ (dew_status == "CONVERGED")

    strict = max(bubble_rp, bubble_rmu, dew_rp, dew_rmu) <= 1.0e-8

    if both_conv and strict and (not ext_checked or ext_ok):
        return "GREEN"
    if both_conv:
        return "YELLOW"
    if one_conv:
        return "YELLOW"
    return "RED"


def main() -> None:
    """
    Purpose
    -------
    Create trust-map figure and trust-label CSV for R515A envelope review.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    Writes:
      verification/r515a_trust_map_20260303.png
      verification/r515a_trust_map_20260303.csv

    Assumptions
    -----------
    Required diagnostic and validation CSV files already exist.

    Failure modes
    -------------
    Raises if inputs are missing.

    References
    ----------
    Project diagnostics and validation outputs.

    Numerical stability notes
    -------------------------
    Not applicable.
    """
    diag_path = Path("verification/r515a_true_vle_diagnostics_20260303.csv")
    val_path = Path("verification/r515a_reference_suite_validation.csv")
    out_fig = Path("verification/r515a_trust_map_20260303.png")
    out_csv = Path("verification/r515a_trust_map_20260303.csv")

    diag = read_diag(diag_path)
    t_ref, p_err_pct, h_err_pct = read_validation_point(val_path)

    ext_ok = (p_err_pct <= 5.0) and (h_err_pct <= 5.0)

    t_vals = sorted(diag.keys())
    trust_rows = []
    for t_k in t_vals:
        b = diag[t_k].get("bubble_liq")
        d = diag[t_k].get("dew_vap")
        if b is None or d is None:
            continue
        ext_checked = abs(t_k - t_ref) <= 0.75
        trust = classify_trust(
            bubble_status=str(b["status"]),
            dew_status=str(d["status"]),
            bubble_rp=float(b["r_P"]),
            bubble_rmu=float(b["r_mu"]),
            dew_rp=float(d["r_P"]),
            dew_rmu=float(d["r_mu"]),
            ext_checked=ext_checked,
            ext_ok=ext_ok,
        )
        trust_rows.append(
            {
                "T_K": t_k,
                "T_C": t_k - 273.15,
                "bubble_status": b["status"],
                "dew_status": d["status"],
                "bubble_rP": b["r_P"],
                "bubble_rmu": b["r_mu"],
                "dew_rP": d["r_P"],
                "dew_rmu": d["r_mu"],
                "trust_class": trust,
                "external_check_used": ext_checked,
                "external_check_passed": ext_ok if ext_checked else "",
            }
        )

    # Write trust CSV.
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(trust_rows[0].keys()))
        w.writeheader()
        w.writerows(trust_rows)

    # Build the trust map plot.
    t_c = np.array([float(r["T_C"]) for r in trust_rows], dtype=float)
    y = np.ones_like(t_c)
    cls = [str(r["trust_class"]) for r in trust_rows]
    colors = ["#2ca02c" if c == "GREEN" else "#ffbf00" if c == "YELLOW" else "#d62728" for c in cls]

    fig, ax = plt.subplots(figsize=(9.0, 3.5), dpi=180)
    ax.scatter(t_c, y, c=colors, s=48, edgecolors="black", linewidths=0.3)

    # draw reference point marker
    ax.axvline(t_ref - 273.15, color="#1f77b4", linestyle="--", linewidth=1.2, label="External validation point (277.6 K)")

    ax.set_ylim(0.7, 1.3)
    ax.set_yticks([])
    ax.set_xlabel("Temperature [C]")
    ax.set_title("R515A p-H Trust Map (Convergence + Available External Validation)")
    ax.grid(True, axis="x", alpha=0.25)

    # Legend proxies
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label="GREEN: both branches converged + strict residuals", markerfacecolor="#2ca02c", markeredgecolor="black", markersize=7),
        Line2D([0], [0], marker="o", color="w", label="YELLOW: partial or looser convergence", markerfacecolor="#ffbf00", markeredgecolor="black", markersize=7),
        Line2D([0], [0], marker="o", color="w", label="RED: neither branch converged", markerfacecolor="#d62728", markeredgecolor="black", markersize=7),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=7)

    ax.text(
        0.99,
        0.03,
        f"External point mean abs errors at 277.6 K: P={p_err_pct:.2f}%, h_fg={h_err_pct:.2f}%",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color="#1f77b4",
    )

    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")
    print(f"Saved: {out_fig}")
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
