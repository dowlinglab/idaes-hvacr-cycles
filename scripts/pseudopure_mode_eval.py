#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Step-4 diagnostics for near-azeotropic pseudo-pure mode using
fixed-composition (x=y=z) pressure-density roots and P-T/P-h comparison plots.
Dependencies: numpy, pandas, scipy, matplotlib, mixture_true_vle_copy, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

FROZEN REFERENCE NOTE
---------------------
This file is intentionally separate from mixture_vle_true_reference.py.
No edits are made to the frozen true-VLE reference implementation.

# BREADCRUMB:
# Date: 2026-03-03
# Rationale: Step-4 Case A (near-azeotropic) diagnostics requested by protocol.
# Assumptions:
# - Use fixed composition z from Honeywell mass fraction conversion.
# - Use midpoint temperature T_mid=(T_bubble+T_dew)/2 from fixed-P glide test.
# - Solve two density roots rho_l/rho_v for P(T_mid, rho, z)=P_target.
# TODO: Replace T_mid heuristic with direct constrained pseudo-pure flash if approved.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import mixture_true_vle_copy as m
from density_bracketing import (
    bracket_vapor_liquid_density_roots,
    default_density_bounds_from_critical,
    format_bounds_log,
)
from linear_model_codex import load_idaes_helmholtz_json, mw_from_json


def run_pseudopure_case_a(
    glide_csv: Path,
    out_csv: Path,
    out_plot: Path,
    out_summary_json: Path,
    fluid1: str = "r1234ze",
    fluid2: str = "r227ea",
    w1: float = 0.911,
) -> Dict[str, float]:
    """
    Purpose
    -------
    Build pseudo-pure saturation diagnostics for near-azeotropic mode.

    Inputs
    ------
    glide_csv : Path [file path]
        Output from fixed-pressure glide test.
    out_csv : Path [file path]
        Pseudo-pure state table destination.
    out_plot : Path [file path]
        Side-by-side comparison figure destination.
    out_summary_json : Path [file path]
        Summary metrics output.
    fluid1, fluid2 : str [unitless]
    w1 : float [kg/kg]
        Honeywell mass fraction of component 1.

    Outputs
    -------
    metrics : dict [mixed units]

    Assumptions
    -----------
    - Near-azeotropic mode represented by fixed composition x=y=z.
    - Midpoint temperature from glide test is acceptable diagnostic proxy.

    Failure modes
    -------------
    - Rows can fail if density roots are not bracketed/found.

    References
    ----------
    - Lever-rule pseudo-pure charting approach for near-azeotropic diagnostics.

    Notes on numerical stability
    ----------------------------
    - Uses bracketing + Brent root solves for low/high density roots.
    """
    gdf = pd.read_csv(glide_csv).sort_values("P_kPa")
    d1 = load_idaes_helmholtz_json(fluid1)
    d2 = load_idaes_helmholtz_json(fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = (w1 / mw1) / ((w1 / mw1) + ((1.0 - w1) / mw2))
    z2 = 1.0 - z1
    mw_mix = z1 * mw1 + z2 * mw2
    rhoc1_mass = float(d1["basic"]["rhoc"])
    rhoc2_mass = float(d2["basic"]["rhoc"])
    bounds = default_density_bounds_from_critical(
        (rhoc1_mass, rhoc2_mass),
        rho_min_mass_kgm3=1.0e-4,
        rho_split_mass_kgm3=50.0,
        rho_max_cap_mass_kgm3=2000.0,
        rhoc_scale_factor=3.0,
    )

    def _mass_to_molar(rho_mass_kgm3: float) -> float:
        return float(rho_mass_kgm3 / mw_mix)

    h_jump_max_kjkg = 20.0
    rows: List[Dict[str, float | str]] = []
    candidate_rows: List[Dict[str, float | str]] = []
    date_tag = datetime.now().strftime("%Y%m%d")
    rho_l_prev: float | None = None
    h_l_prev_kjkg: float | None = None

    def _all_sign_change_brackets(grid_molm3: np.ndarray) -> List[Tuple[float, float]]:
        vals = np.array([f_rho(float(x)) for x in grid_molm3], dtype=float)
        out: List[Tuple[float, float]] = []
        for i in range(len(grid_molm3) - 1):
            v1 = vals[i]
            v2 = vals[i + 1]
            if not (np.isfinite(v1) and np.isfinite(v2)):
                continue
            if v1 == 0.0:
                out.append((float(grid_molm3[i]), float(grid_molm3[i])))
            elif v1 * v2 < 0.0:
                out.append((float(grid_molm3[i]), float(grid_molm3[i + 1])))
        return out

    for _, r in gdf.iterrows():
        P_kpa = float(r["P_kPa"])
        Tb = float(r["T_bubble_K"])
        Td = float(r["T_dew_K"])
        if not (np.isfinite(Tb) and np.isfinite(Td)):
            rows.append({
                "P_kPa": P_kpa,
                "T_mid_K": np.nan,
                "rho_l_molm3": np.nan,
                "rho_v_molm3": np.nan,
                "h_l_kJkg": np.nan,
                "h_v_kJkg": np.nan,
                "max_abs_ln_f_ratio": np.nan,
                "status": "FAILED_NO_T",
            })
            continue

        Tm = 0.5 * (Tb + Td)

        def f_rho(rho_molm3: float) -> float:
            return m.mix_state(d1, d2, Tm, rho_molm3, z1).p_pa * 1e-3 - P_kpa

        vap_lo = _mass_to_molar(bounds.rho_min_mass_kgm3)
        vap_hi = _mass_to_molar(bounds.rho_split_mass_kgm3)
        liq_lo = _mass_to_molar(bounds.rho_split_mass_kgm3)
        liq_hi = _mass_to_molar(bounds.rho_max_mass_kgm3)
        vapor_grid = np.logspace(np.log10(vap_lo), np.log10(vap_hi), 260)
        liquid_grid = np.logspace(np.log10(liq_lo), np.log10(liq_hi), 320)
        full_grid = np.logspace(np.log10(vap_lo), np.log10(liq_hi), 1400)
        vapor_brackets = _all_sign_change_brackets(vapor_grid)
        liquid_brackets = _all_sign_change_brackets(liquid_grid)
        full_brackets = _all_sign_change_brackets(full_grid)

        for phase_window, brackets in (("vapor", vapor_brackets), ("liquid", liquid_brackets), ("full", full_brackets)):
            for root_rank, (blo, bhi) in enumerate(brackets):
                rho_root = np.nan
                h_root = np.nan
                root_status = "OK"
                try:
                    rho_root = float(brentq(f_rho, blo, bhi, xtol=1e-10, rtol=1e-10, maxiter=200))
                    st_root = m.mix_state(d1, d2, Tm, rho_root, z1)
                    h_root = float((st_root.h_jmol / mw_mix) * 1e-3)
                except Exception:
                    root_status = "FAILED_ROOT"
                candidate_rows.append({
                    "P_kPa": P_kpa,
                    "T_mid_K": Tm,
                    "phase_window": phase_window,
                    "root_rank": int(root_rank),
                    "rho_bracket_lo": float(blo),
                    "rho_bracket_hi": float(bhi),
                    "rho_root": float(rho_root),
                    "h_root_kJkg": float(h_root),
                    "root_status": root_status,
                })

        all_ok = [
            row for row in candidate_rows
            if row["P_kPa"] == P_kpa and row["phase_window"] == "full" and row["root_status"] == "OK"
        ]
        all_ok_sorted = sorted(all_ok, key=lambda rr: float(rr["rho_root"]))
        vapor_ok_sorted = all_ok_sorted[:1]
        liquid_ok_sorted = all_ok_sorted[1:]

        br = bracket_vapor_liquid_density_roots(
            func_of_rho_molm3=f_rho,
            mw_mix_kgmol=mw_mix,
            bounds=bounds,
            vapor_points=260,
            liquid_points=320,
            fallback_points=1400,
        )
        low_br = br.vapor_bracket_molm3
        high_br = br.liquid_bracket_molm3

        if low_br is None or high_br is None:
            bounds_log = format_bounds_log(bounds)
            print(
                f"[density-bracket-fail] P={P_kpa} kPa, T_mid={Tm} K, "
                f"note={br.note}, bounds={bounds_log}"
            )
            rows.append({
                "P_kPa": P_kpa,
                "T_mid_K": Tm,
                "rho_l_molm3": np.nan,
                "rho_v_molm3": np.nan,
                "h_l_kJkg": np.nan,
                "h_v_kJkg": np.nan,
                "max_abs_ln_f_ratio": np.nan,
                "bracket_note": br.note,
                "bracket_bounds": bounds_log,
                "status": "FAILED_NO_RHO_BRACKET",
            })
            continue

        try:
            rv = float(brentq(f_rho, low_br[0], low_br[1], xtol=1e-10, rtol=1e-10, maxiter=200))
            rl = float(brentq(f_rho, high_br[0], high_br[1], xtol=1e-10, rtol=1e-10, maxiter=200))
            liquid_select_note = "LIQUID_STANDARD_BRACKET_ROOT"
            if len(liquid_ok_sorted) > 0:
                if rho_l_prev is None:
                    # First accepted point: keep standard window-bracket convention.
                    bracket_lo, bracket_hi = high_br
                    matched = [
                        rr for rr in liquid_ok_sorted
                        if abs(float(rr["rho_root"]) - rl) <= 1e-8 * max(1.0, abs(rl))
                        and abs(float(rr["rho_bracket_lo"]) - bracket_lo) <= 1e-12 * max(1.0, abs(bracket_lo))
                        and abs(float(rr["rho_bracket_hi"]) - bracket_hi) <= 1e-12 * max(1.0, abs(bracket_hi))
                    ]
                    if len(matched) == 0:
                        rl = float(liquid_ok_sorted[-1]["rho_root"])
                        liquid_select_note = "LIQUID_FIRSTPOINT_FALLBACK_HIGHEST_DENSITY"
                else:
                    nearest = min(liquid_ok_sorted, key=lambda rr: abs(float(rr["rho_root"]) - rho_l_prev))
                    rl = float(nearest["rho_root"])
                    liquid_select_note = "LIQUID_CONTINUITY_NEAREST_RHO_PREV"
            else:
                rows.append({
                    "P_kPa": P_kpa,
                    "T_mid_K": Tm,
                    "rho_l_molm3": np.nan,
                    "rho_v_molm3": np.nan,
                    "h_l_kJkg": np.nan,
                    "h_v_kJkg": np.nan,
                    "max_abs_ln_f_ratio": np.nan,
                    "bracket_note": br.note,
                    "bracket_bounds": format_bounds_log(bounds),
                    "liquid_select_note": "FAILED_NO_LIQUID_CANDIDATE_FOR_SELECTION",
                    "vapor_select_note": "",
                    "h_jump_warn": "",
                    "status": "FAILED_NO_LIQUID_CANDIDATE",
                })
                continue

            # Keep vapor selection anchored to standard window bracket if present.
            vapor_select_note = "VAPOR_STANDARD_BRACKET_ROOT"
            if len(vapor_ok_sorted) > 0:
                bracket_lo, bracket_hi = low_br
                matched = [
                    rr for rr in vapor_ok_sorted
                    if abs(float(rr["rho_root"]) - rv) <= 1e-8 * max(1.0, abs(rv))
                    and abs(float(rr["rho_bracket_lo"]) - bracket_lo) <= 1e-12 * max(1.0, abs(bracket_lo))
                    and abs(float(rr["rho_bracket_hi"]) - bracket_hi) <= 1e-12 * max(1.0, abs(bracket_hi))
                ]
                if len(matched) == 0:
                    rv = float(vapor_ok_sorted[0]["rho_root"])
                    vapor_select_note = "VAPOR_FIRST_AVAILABLE_IN_WINDOW"

            st_v = m.mix_state(d1, d2, Tm, rv, z1)
            st_l = m.mix_state(d1, d2, Tm, rl, z1)
            h_l_kjkg = float((st_l.h_jmol / mw_mix) * 1e-3)
            h_v_kjkg = float((st_v.h_jmol / mw_mix) * 1e-3)
            mu1_l, mu2_l = m.chemical_potentials_analytic(d1, d2, Tm, rl, z1)
            mu1_v, mu2_v = m.chemical_potentials_analytic(d1, d2, Tm, rv, z1)
            lnfr1 = (mu1_l - mu1_v) / (m.R_u * Tm)
            lnfr2 = (mu2_l - mu2_v) / (m.R_u * Tm)
            jump_warn = ""
            if h_l_prev_kjkg is not None and abs(h_l_kjkg - h_l_prev_kjkg) > h_jump_max_kjkg:
                jump_warn = f"HJUMP_WARN_{abs(h_l_kjkg - h_l_prev_kjkg):.6f}_kJkg"
                print(
                    f"[h-jump-warning] P={P_kpa} kPa, "
                    f"h_l_prev={h_l_prev_kjkg:.6f} kJ/kg, h_l={h_l_kjkg:.6f} kJ/kg, "
                    f"delta={h_l_kjkg-h_l_prev_kjkg:.6f} kJ/kg, candidates_csv=diagnostics/pseudopure_root_candidates_{date_tag}.csv"
                )
            rows.append({
                "P_kPa": P_kpa,
                "T_mid_K": Tm,
                "rho_l_molm3": rl,
                "rho_v_molm3": rv,
                "h_l_kJkg": h_l_kjkg,
                "h_v_kJkg": h_v_kjkg,
                "max_abs_ln_f_ratio": max(abs(lnfr1), abs(lnfr2)),
                "bracket_note": br.note,
                "bracket_bounds": format_bounds_log(bounds),
                "liquid_select_note": liquid_select_note,
                "vapor_select_note": vapor_select_note,
                "h_jump_warn": jump_warn,
                "status": "OK",
            })
            rho_l_prev = rl
            h_l_prev_kjkg = h_l_kjkg
        except Exception:
            bounds_log = format_bounds_log(bounds)
            print(
                f"[density-root-fail] P={P_kpa} kPa, T_mid={Tm} K, "
                f"note={br.note}, bounds={bounds_log}"
            )
            rows.append({
                "P_kPa": P_kpa,
                "T_mid_K": Tm,
                "rho_l_molm3": np.nan,
                "rho_v_molm3": np.nan,
                "h_l_kJkg": np.nan,
                "h_v_kJkg": np.nan,
                "max_abs_ln_f_ratio": np.nan,
                "bracket_note": br.note,
                "bracket_bounds": bounds_log,
                "status": "FAILED_ROOT",
            })

    odf = pd.DataFrame(rows).sort_values("P_kPa")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    odf.to_csv(out_csv, index=False)
    candidates_csv = out_csv.parent / f"pseudopure_root_candidates_{date_tag}.csv"
    pd.DataFrame(candidate_rows).to_csv(candidates_csv, index=False)

    ok = odf[odf["status"] == "OK"].copy()
    ref_pt = pd.read_csv("verification/r515b_honeywell_pt_comparison_tuned_solver_20260303.csv")[["T_C", "P_chart_kPa"]]
    ref_pt = ref_pt.rename(columns={"P_chart_kPa": "P_kPa"})
    merged = ref_pt.merge(odf[["P_kPa", "T_mid_K", "status"]], on="P_kPa", how="left")
    merged["dT_mid_minus_honeywell_K"] = merged["T_mid_K"] - (merged["T_C"] + 273.15)

    # Side-by-side figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.4), dpi=170)

    ax1.plot(merged["T_C"] + 273.15, merged["P_kPa"], "o", ms=4, label="Honeywell PT")
    ok_pt = merged[merged["status"] == "OK"]
    ax1.plot(ok_pt["T_mid_K"], ok_pt["P_kPa"], "s", ms=3.8, label="Pseudo-pure T_mid(P)")
    ax1.set_xlabel("Temperature [K]")
    ax1.set_ylabel("Pressure [kPa]")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    if len(ok):
        ax2.plot(ok["h_l_kJkg"], ok["P_kPa"] * 1e-2, lw=2.0, label="Pseudo-pure sat. liquid")
        ax2.plot(ok["h_v_kJkg"], ok["P_kPa"] * 1e-2, lw=2.0, label="Pseudo-pure sat. vapor")
        ax2.fill_betweenx(ok["P_kPa"] * 1e-2, ok["h_l_kJkg"], ok["h_v_kJkg"], alpha=0.15)
    ax2.set_yscale("log")
    ax2.set_xlabel("Enthalpy [kJ/kg]")
    ax2.set_ylabel("Pressure [bar]")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)

    fig.suptitle("Step-4 Case A: Pseudo-pure/near-azeotropic diagnostic mode", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, bbox_inches="tight")

    mae_t = float(np.nanmean(np.abs(merged["dT_mid_minus_honeywell_K"])))
    max_t = float(np.nanmax(np.abs(merged["dT_mid_minus_honeywell_K"])))
    mape_p = float(np.nanmean(np.abs((ok["P_kPa"] - ok["P_kPa"]) / ok["P_kPa"])) * 100.0) if len(ok) else np.nan
    # p metric above is trivially zero because pressure is fixed by construction in this mode.

    metrics = {
        "n_total": int(len(odf)),
        "n_ok": int((odf["status"] == "OK").sum()),
        "n_failed": int((odf["status"] != "OK").sum()),
        "mae_Tmid_vs_honeywell_K": mae_t,
        "max_abs_Tmid_vs_honeywell_K": max_t,
        "pseudo_mode_pressure_error_metric_pct": mape_p,
        "csv": str(out_csv),
        "candidates_csv": str(candidates_csv),
        "plot": str(out_plot),
    }

    with out_summary_json.open("w") as f:
        json.dump(metrics, f, indent=2)

    # One requested diagnostic: P(rho) sweep at 2400 kPa using T_mid from glide row.
    sweep_row = gdf[np.isclose(gdf["P_kPa"], 2400.0)]
    if len(sweep_row):
        Tm_2400 = float(0.5 * (float(sweep_row.iloc[0]["T_bubble_K"]) + float(sweep_row.iloc[0]["T_dew_K"])))
        rho_mass_sweep = np.logspace(np.log10(bounds.rho_min_mass_kgm3), np.log10(bounds.rho_max_mass_kgm3), 500)
        sweep = []
        for rho_mass in rho_mass_sweep:
            rho_mol = _mass_to_molar(float(rho_mass))
            try:
                p_kpa = m.mix_state(d1, d2, Tm_2400, rho_mol, z1).p_pa * 1e-3
            except Exception:
                p_kpa = np.nan
            sweep.append({"rho_mass_kgm3": float(rho_mass), "rho_molm3": float(rho_mol), "P_kPa": float(p_kpa)})
        sweep_df = pd.DataFrame(sweep)
        sweep_csv = out_csv.parent / f"pseudopure_prho_sweep_2400kPa_{datetime.now().strftime('%Y%m%d')}.csv"
        sweep_df.to_csv(sweep_csv, index=False)
        metrics["prho_sweep_2400kPa_csv"] = str(sweep_csv)
        metrics["prho_sweep_T_mid_K"] = Tm_2400
        with out_summary_json.open("w") as f:
            json.dump(metrics, f, indent=2)

    return metrics


def main() -> None:
    """
    Purpose
    -------
    CLI entrypoint for Step-4 pseudo-pure diagnostics.

    Inputs
    ------
    --glide-csv : str [path]
    --out-csv : str [path]
    --out-plot : str [path]
    --out-summary : str [path]

    Outputs
    -------
    None [unitless]
      Writes diagnostics CSV, plot, and summary JSON.

    Assumptions
    -----------
    - Step-1 glide CSV exists and is valid.

    Failure modes
    -------------
    - Raises on file I/O or solver errors that escape row-level handling.

    References
    ----------
    - Step-4 protocol for near-azeotropic pseudo-pure charting mode.

    Notes on numerical stability
    ----------------------------
    - Row-level failures are captured with status flags in output CSV.
    """
    parser = argparse.ArgumentParser(description="Run pseudo-pure Case-A diagnostics from fixed-P glide data.")
    parser.add_argument("--glide-csv", default="diagnostics/honeywell_glide_test_20260303.csv")
    parser.add_argument("--out-csv", default=f"diagnostics/r515b_pseudopure_dome_{datetime.now().strftime('%Y%m%d')}.csv")
    parser.add_argument("--out-plot", default=f"diagnostics/r515b_pseudopure_vs_honeywell_{datetime.now().strftime('%Y%m%d')}.png")
    parser.add_argument("--out-summary", default=f"diagnostics/r515b_pseudopure_summary_{datetime.now().strftime('%Y%m%d')}.json")
    args = parser.parse_args()

    metrics = run_pseudopure_case_a(Path(args.glide_csv), Path(args.out_csv), Path(args.out_plot), Path(args.out_summary))
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
