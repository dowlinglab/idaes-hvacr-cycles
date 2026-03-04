#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: OpenAI Codex (version: GPT-5)
QA/testing: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Layer-by-layer P-h diagram plotting driver scaffold for R515B.
Dependencies: argparse, json, pathlib, datetime, matplotlib
Context reference: PROJECT_CONTEXT.md

Context breadcrumbs
-------------------
Date: 2026-03-03
Purpose: STEP 0 scaffold only (no thermodynamic calculations)
Artifacts produced:
- diagnostics/plots/ph_layer0_scaffold_<date>.png
- diagnostics/plots/ph_layer_status_<date>.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def run_step0(out_dir: Path, stamp: str) -> tuple[Path, Path]:
    """
    Purpose
    -------
    Generate STEP 0 scaffold plot and layer status JSON.

    Inputs
    ------
    out_dir : Path [filesystem path]
    stamp : str [YYYYMMDD]

    Outputs
    -------
    scaffold_path : Path [filesystem path]
    status_path : Path [filesystem path]

    Assumptions
    -----------
    - Output directory is writable.

    Failure modes
    -------------
    - Raises on file/plot write errors.

    References
    ----------
    - Protocol: Layer-by-layer P-h Diagram Build, STEP 0.

    Notes on numerical stability
    ----------------------------
    - Not applicable (no numerical solve in STEP 0).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.0, 6.0), dpi=170)
    ax.set_yscale("log")
    ax.set_xlabel("Specific enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("R515B P-h Diagram Layers (STEP 0 scaffold)")
    ax.grid(True, which="both", alpha=0.28)

    scaffold_path = out_dir / f"ph_layer0_scaffold_{stamp}.png"
    fig.tight_layout()
    fig.savefig(scaffold_path, bbox_inches="tight")

    status = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "completed_layers": [0],
        "inputs_used": [],
        "note": "STEP 0 scaffold only; no data layers plotted.",
    }
    status_path = out_dir / f"ph_layer_status_{stamp}.json"
    with status_path.open("w") as f:
        json.dump(status, f, indent=2)

    return scaffold_path, status_path


def run_step6(out_dir: Path, stamp: str, status_path: Path) -> tuple[Path, Path, Path]:
    """
    Purpose
    -------
    Build STEP 6 Honeywell-style final P-h figure using existing layer CSV artifacts.

    Inputs
    ------
    out_dir : Path [filesystem path]
    stamp : str [YYYYMMDD]
    status_path : Path [filesystem path]
        Existing layer status JSON to update.

    Outputs
    -------
    png_path : Path [filesystem path]
    pdf_path : Path [filesystem path]
    manifest_path : Path [filesystem path]

    Assumptions
    -----------
    - Pseudo-pure saturation CSV artifact exists in diagnostics/plots.

    Failure modes
    -------------
    - Raises if required input CSVs are missing or malformed.

    References
    ----------
    - Layer-by-layer P-h Diagram Build Protocol, STEP 6.

    Notes on numerical stability
    ----------------------------
    - Plot-only step; no thermodynamic solve is performed.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    layer5_csv = out_dir / f"ph_layer5_pseudopure_overlay_{stamp}.csv"

    d5 = pd.read_csv(layer5_csv).sort_values("P_bar")
    d5_ok = d5[d5["status"] == "OK"].copy()

    fig, ax = plt.subplots(figsize=(10.8, 6.6), dpi=190)
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_xlim(150.0, 500.0)
    ax.set_ylim(1.0, 35.0)
    ax.grid(True, which="both", alpha=0.30)

    # Pseudo-pure saturation boundaries (Honeywell-style chart representation).
    ax.fill_betweenx(
        d5_ok["P_bar"], d5_ok["h_l_kJkg_pseudopure"], d5_ok["h_v_kJkg_pseudopure"],
        color="#f4d35e", alpha=0.18, label="Two-phase region"
    )
    ax.plot(d5_ok["h_l_kJkg_pseudopure"], d5_ok["P_bar"], color="#0b4f6c", lw=2.0, label="Sat. liquid boundary")
    ax.plot(d5_ok["h_v_kJkg_pseudopure"], d5_ok["P_bar"], color="#c44536", lw=2.0, label="Sat. vapor boundary")

    # Quality lines x=0.1..0.9 on pseudo-pure saturation envelope.
    for q in np.arange(0.1, 1.0, 0.1):
        hq = (1.0 - q) * d5_ok["h_l_kJkg_pseudopure"].to_numpy(dtype=float) + q * d5_ok["h_v_kJkg_pseudopure"].to_numpy(dtype=float)
        ax.plot(hq, d5_ok["P_bar"], color="#6a4c93", lw=0.7, alpha=0.45)
    ax.plot([], [], color="#6a4c93", lw=1.0, alpha=0.7, label="Quality lines x=0.1..0.9")

    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("R515B P-h Diagram (Pseudo-pure Honeywell-style composite)")
    fig.tight_layout()

    png_path = out_dir / f"ph_layer6_final_honeywell_style_{stamp}.png"
    pdf_path = out_dir / f"ph_layer6_final_honeywell_style_{stamp}.pdf"
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")

    git_hash = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_hash_short": git_hash,
        "inputs": [
            str(layer5_csv),
        ],
        "outputs": [
            str(png_path),
            str(pdf_path),
        ],
        "axis_limits": {"h_kJkg": [150.0, 500.0], "p_bar": [1.0, 35.0]},
        "log_pressure_axis": True,
    }
    manifest_path = out_dir / f"ph_layer6_manifest_{stamp}.json"
    with manifest_path.open("w") as f:
        json.dump(manifest, f, indent=2)

    status = {}
    if status_path.exists():
        with status_path.open("r") as f:
            status = json.load(f)
    completed = set(status.get("completed_layers", []))
    completed.add(6)
    status["completed_layers"] = sorted(completed)
    status["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
    status["layer6_summary"] = {
        "outputs": [str(png_path), str(pdf_path), str(manifest_path)],
        "legend": [
            "Sat. liquid boundary",
            "Sat. vapor boundary",
            "Two-phase region",
            "Quality lines x=0.1..0.9",
        ],
    }
    with status_path.open("w") as f:
        json.dump(status, f, indent=2)

    return png_path, pdf_path, manifest_path


def main() -> None:
    """
    Purpose
    -------
    CLI entrypoint for P-h layering scaffold generation.

    Inputs
    ------
    --out-dir : str [filesystem path]
    --stamp : str [YYYYMMDD], optional

    Outputs
    -------
    None [unitless]

    Assumptions
    -----------
    - Called for STEP 0 only in current protocol stage.

    Failure modes
    -------------
    - Propagates exceptions from plotting/I/O.

    References
    ----------
    - Layer-by-layer protocol STEP 0.

    Notes on numerical stability
    ----------------------------
    - Not applicable.
    """
    parser = argparse.ArgumentParser(description="Layered P-h workflow plot driver.")
    parser.add_argument("--out-dir", default="diagnostics/plots")
    parser.add_argument("--stamp", default=datetime.now().strftime("%Y%m%d"))
    parser.add_argument("--step", choices=["0", "6"], default="0")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    stamp = str(args.stamp)
    if args.step == "0":
        scaffold_path, status_path = run_step0(out_dir, stamp)
        print(f"scaffold_png={scaffold_path}")
        print(f"status_json={status_path}")
        return

    status_path = out_dir / f"ph_layer_status_{stamp}.json"
    png_path, pdf_path, manifest_path = run_step6(out_dir, stamp, status_path)
    print(f"layer6_png={png_path}")
    print(f"layer6_pdf={pdf_path}")
    print(f"layer6_manifest={manifest_path}")
    print(f"status_json={status_path}")


if __name__ == "__main__":
    main()
