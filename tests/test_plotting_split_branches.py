# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-03
Purpose of file: Unit tests for crossover detection and split-branch plotting diagnostics.
Dependencies: pytest, numpy, csv, pathlib, scripts.postprocess_true_vle_outputs
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import numpy as np

from scripts.postprocess_true_vle_outputs import detect_crossover_indices, get_crossovers


# === SECTION: Synthetic Crossover Fixture ===
# Rationale: Build a minimal bubble/dew dataset with a known h-order inversion.
def _write_branch_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "status",
        "T_K",
        "P_Pa",
        "rho_l_molm3",
        "rho_v_molm3",
        "x1_liq",
        "y1_vap",
        "h_l_Jmol",
        "h_v_Jmol",
        "r_P",
        "r_mu",
        "iterations",
        "notes",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def test_detect_crossover_indices_synthetic() -> None:
    h_liq = np.array([100.0, 210.0, 180.0], dtype=float)
    h_vap = np.array([300.0, 205.0, 220.0], dtype=float)
    idx = detect_crossover_indices(h_liq, h_vap)
    assert idx == [1]


def test_postprocess_writes_crossover_csv_and_query(tmp_path: Path) -> None:
    bubble_csv = tmp_path / "bubble.csv"
    dew_csv = tmp_path / "dew.csv"
    conv_csv = tmp_path / "conv.csv"
    diag_csv = tmp_path / "diag.csv"
    cross_csv = tmp_path / "cross.csv"
    fig_clean = tmp_path / "clean.png"
    fig_fail = tmp_path / "fail.png"
    fig_pt = tmp_path / "pt.png"
    fig_ht = tmp_path / "ht.png"
    fig_ph_col = tmp_path / "ph_col.png"

    bubble_rows = [
        {
            "status": "CONVERGED",
            "T_K": "300.0",
            "P_Pa": "500000",
            "rho_l_molm3": "10000",
            "rho_v_molm3": "100",
            "x1_liq": "0.6",
            "y1_vap": "0.8",
            "h_l_Jmol": "12000",
            "h_v_Jmol": "25000",
            "r_P": "1e-12",
            "r_mu": "1e-12",
            "iterations": "20",
            "notes": "ok",
        },
        {
            "status": "CONVERGED",
            "T_K": "301.0",
            "P_Pa": "550000",
            "rho_l_molm3": "9800",
            "rho_v_molm3": "120",
            "x1_liq": "0.61",
            "y1_vap": "0.79",
            "h_l_Jmol": "26000",
            "h_v_Jmol": "24000",
            "r_P": "1e-12",
            "r_mu": "1e-12",
            "iterations": "22",
            "notes": "ok",
        },
    ]

    dew_rows = [
        {
            "status": "CONVERGED",
            "T_K": "300.0",
            "P_Pa": "510000",
            "rho_l_molm3": "9500",
            "rho_v_molm3": "110",
            "x1_liq": "0.55",
            "y1_vap": "0.82",
            "h_l_Jmol": "13000",
            "h_v_Jmol": "27000",
            "r_P": "1e-12",
            "r_mu": "1e-12",
            "iterations": "21",
            "notes": "ok",
        },
        {
            "status": "CONVERGED",
            "T_K": "301.0",
            "P_Pa": "560000",
            "rho_l_molm3": "9400",
            "rho_v_molm3": "130",
            "x1_liq": "0.56",
            "y1_vap": "0.81",
            "h_l_Jmol": "14000",
            "h_v_Jmol": "20000",
            "r_P": "1e-12",
            "r_mu": "1e-12",
            "iterations": "20",
            "notes": "ok",
        },
    ]

    _write_branch_csv(bubble_csv, bubble_rows)
    _write_branch_csv(dew_csv, dew_rows)

    subprocess.run(
        [
            "python",
            "scripts/postprocess_true_vle_outputs.py",
            "--fluid1",
            "r1234ze",
            "--fluid2",
            "r227ea",
            "--w1",
            "0.911",
            "--bubble-csv",
            str(bubble_csv),
            "--dew-csv",
            str(dew_csv),
            "--diag-csv",
            str(diag_csv),
            "--conv-csv",
            str(conv_csv),
            "--cross-csv",
            str(cross_csv),
            "--fig-clean",
            str(fig_clean),
            "--fig-fail",
            str(fig_fail),
            "--fig-pt",
            str(fig_pt),
            "--fig-ht",
            str(fig_ht),
            "--fig-ph-colored",
            str(fig_ph_col),
        ],
        check=True,
    )

    assert cross_csv.exists()
    cross = get_crossovers(conv_csv)
    assert len(cross) == 1
    assert abs(cross[0]["T_K"] - 301.0) < 1e-12
    assert fig_clean.exists()
    assert fig_fail.exists()
