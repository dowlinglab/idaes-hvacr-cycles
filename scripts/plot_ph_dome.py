#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: CLI wrapper to compute and plot pure-fluid P-h saturation dome.
Dependencies: numpy, helmholtz_saturation
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from helmholtz_saturation import compute_saturation_dome, plot_ph_dome_from_saturation_data, save_dome_csv


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute pure-fluid saturation dome and save CSV/plot.")
    parser.add_argument("--fluid", required=True, help="fluid stem, e.g. r1234ze")
    parser.add_argument("--Tmin", required=True, type=float, help="minimum temperature [K]")
    parser.add_argument("--Tmax", required=True, type=float, help="maximum temperature [K], below Tc")
    parser.add_argument("--n", default=200, type=int, help="number of temperature points")
    parser.add_argument("--tol", default=1e-10, type=float, help="Newton residual tolerance")
    parser.add_argument("--maxiter", default=50, type=int, help="max Newton iterations per point")
    parser.add_argument("--csv", default="dome.csv", help="output CSV path")
    parser.add_argument("--fig", default="ph_dome.pdf", help="output figure path")
    parser.add_argument("--metadata", default="verification/saturation_run_metadata.json", help="output metadata JSON")
    args = parser.parse_args()

    T_vals = np.linspace(args.Tmin, args.Tmax, args.n)
    dome = compute_saturation_dome(args.fluid, T_vals, tol=args.tol, maxiter=args.maxiter)

    save_dome_csv(dome, args.csv)
    plot_ph_dome_from_saturation_data(dome["h_l_kJkg"], dome["h_v_kJkg"], dome["p_kPa"], args.fig)

    meta = {
        "fluid": args.fluid,
        "Tmin_K": float(args.Tmin),
        "Tmax_K": float(args.Tmax),
        "n_points_requested": int(args.n),
        "n_points_solved": int(len(dome["T_K"])),
        "tol": float(args.tol),
        "maxiter": int(args.maxiter),
        "critical_reached": bool(dome["critical_reached"][0]),
        "Tc_K": float(dome["Tc_K"][0]),
        "outputs": {"csv": args.csv, "fig": args.fig},
    }
    meta_path = Path(args.metadata)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    with meta_path.open("w") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved CSV: {args.csv}")
    print(f"Saved figure: {args.fig}")
    print(f"Saved metadata: {args.metadata}")


if __name__ == "__main__":
    main()
