#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/testing: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Audit 70 C pseudo-pure isotherm root selection by logging all
sign-change brackets, selected bracket/root, and fallback usage for vapor and
liquid branches.
Dependencies: numpy, pandas, pressure_validated_model, plot_pseudopure_isodiagram
Context reference: PROJECT_CONTEXT.md

Version tag: v0.1.0

# === SECTION: Root-Selection Audit ===
# Rationale: Diagnose branch loss/jumps by exposing bracket choices and fallback
# behavior without changing EOS or solver physics.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pressure_validated_model import default_interaction, solve_rho_mass_for_P
from scripts.plot_pseudopure_isodiagram import (
    _find_sign_change_bracket_near,
    _solve_rho_pressure_with_bracket,
    _mw_mix_kgmol,
    _dPdrho_num_kpa_per_kgm3,
    interp_pseudopure_sat_at_T,
    load_pseudopure_saturation_table,
    compute_props,
)


@dataclass
class AuditConfig:
    """
    Purpose
    -------
    Collect numerical settings for the 70 C root-selection audit.

    Inputs
    ------
    t_c : float [degC]
        Isotherm temperature.
    w1 : float [kg/kg]
        Mass fraction of component 1.
    pmin_bar : float [bar]
        Lower pressure bound for plotted branch sampling.
    pmax_bar : float [bar]
        Upper pressure bound for plotted branch sampling.
    eps_kpa : float [kPa]
        Offset from saturation pressure for branch endpoints.
    nv : int [-]
        Number of vapor pressure points.
    nl : int [-]
        Number of liquid pressure points.
    rho_floor : float [kg/m^3]
        Minimum density scanned in full bracket search.
    rho_cap : float [kg/m^3]
        Maximum density scanned in full bracket search.
    nscan : int [-]
        Number of scan points for global sign-change detection.

    Outputs
    -------
    AuditConfig [dataclass]
        Immutable-style settings bundle for the diagnostic run.

    Assumptions
    -----------
    - Values are positive and suitable for log-space pressure/density grids.

    Failure modes
    -------------
    - Invalid (non-positive) bounds can cause log-space construction failure.

    References
    ----------
    - Numerical root bracketing diagnostics for Helmholtz EOS pressure roots.

    Notes on numerical stability
    ----------------------------
    - Higher `nscan` improves interval detection resolution at extra cost.
    """

    t_c: float = 70.0
    w1: float = 0.911
    pmin_bar: float = 1.0
    pmax_bar: float = 100.0
    eps_kpa: float = 1.0
    nv: int = 80
    nl: int = 80
    rho_floor: float = 1e-6
    rho_cap: float = 2000.0
    nscan: int = 3000


def _pressure_residual_kpa(t_k: float, rho_mass: float, p_target_kpa: float, w1: float, w2: float) -> float:
    """
    Purpose
    -------
    Evaluate pressure residual F(rho)=P(T,rho)-P_target for audit logging.

    Inputs
    ------
    t_k : float [K]
    rho_mass : float [kg/m^3]
    p_target_kpa : float [kPa]
    w1 : float [kg/kg]
    w2 : float [kg/kg]

    Outputs
    -------
    float [kPa]
        Pressure residual.

    Assumptions
    -----------
    - State evaluation in `compute_props` succeeds for the queried density.

    Failure modes
    -------------
    - Can return NaN/raise if EOS path fails at queried state.

    References
    ----------
    - Root condition for branch tracing: P(T,rho)=P_target.

    Notes on numerical stability
    ----------------------------
    - Residual sign is used for bracket detection; NaN values are ignored.
    """
    st = compute_props("r1234ze.json", "r227ea.json", t_k, float(rho_mass), w1, w2, default_interaction())
    return float(st.p - p_target_kpa)


def _all_sign_change_intervals(
    t_k: float,
    p_target_kpa: float,
    w1: float,
    w2: float,
    rho_floor: float,
    rho_cap: float,
    nscan: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Purpose
    -------
    Find all sign-change density intervals for F(rho)=P(T,rho)-P_target.

    Inputs
    ------
    t_k : float [K]
    p_target_kpa : float [kPa]
    w1 : float [kg/kg]
    w2 : float [kg/kg]
    rho_floor : float [kg/m^3]
    rho_cap : float [kg/m^3]
    nscan : int [-]

    Outputs
    -------
    rho_grid : np.ndarray [kg/m^3]
    idx_pairs : np.ndarray [index pairs]
        Array of index i where interval [i, i+1] contains a sign change.

    Assumptions
    -----------
    - Full scan range includes all physically relevant roots for diagnosis.

    Failure modes
    -------------
    - Missed intervals if scan is too coarse for steep transitions.

    References
    ----------
    - Standard bracketing by sign change in one-dimensional root finding.

    Notes on numerical stability
    ----------------------------
    - Log-space scan mitigates undersampling at low density.
    """
    rho_grid = np.logspace(np.log10(max(rho_floor, 1e-12)), np.log10(rho_cap), int(nscan))
    vals = np.empty_like(rho_grid)
    vals[:] = np.nan
    for i, rr in enumerate(rho_grid):
        try:
            vals[i] = _pressure_residual_kpa(t_k, float(rr), p_target_kpa, w1, w2)
        except Exception:
            vals[i] = np.nan
    signs = np.sign(vals)
    idx_list: List[int] = []
    for i in range(len(rho_grid) - 1):
        s1 = signs[i]
        s2 = signs[i + 1]
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and s1 * s2 <= 0.0:
            idx_list.append(i)
    return rho_grid, np.asarray(idx_list, dtype=int)


def _run_branch_audit(
    cfg: AuditConfig,
    branch: str,
    t_k: float,
    p_grid_kpa: np.ndarray,
    rho_anchor: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Purpose
    -------
    Audit one branch (vapor or liquid) along a fixed pressure grid.

    Inputs
    ------
    cfg : AuditConfig [-]
    branch : str [-]
        `vapor` or `liquid`.
    t_k : float [K]
    p_grid_kpa : np.ndarray [kPa]
    rho_anchor : float [kg/m^3]
        Initial continuation density seed.

    Outputs
    -------
    summary_df : pd.DataFrame
        One row per pressure step with selected root/bracket/fallback markers.
    intervals_df : pd.DataFrame
        One row per detected sign-change interval.

    Assumptions
    -----------
    - Uses the same helper functions as the plotting script for comparability.

    Failure modes
    -------------
    - If all solves fail, selected root columns remain NaN.

    References
    ----------
    - Audit mirrors the current isotherm tracing implementation.

    Notes on numerical stability
    ----------------------------
    - Continuation seed is updated only after a successful accepted root.
    """
    w1 = cfg.w1
    w2 = 1.0 - w1
    rho_prev = float(rho_anchor)
    summary_rows: List[Dict] = []
    interval_rows: List[Dict] = []

    for k, p_target in enumerate(p_grid_kpa):
        rho_grid, idx_pairs = _all_sign_change_intervals(
            t_k, float(p_target), w1, w2, cfg.rho_floor, cfg.rho_cap, cfg.nscan
        )
        for rank, i0 in enumerate(idx_pairs, start=1):
            lo = float(rho_grid[i0])
            hi = float(rho_grid[i0 + 1])
            interval_rows.append(
                {
                    "branch": branch,
                    "step_index": k,
                    "P_target_kPa": float(p_target),
                    "interval_rank": rank,
                    "rho_lo_kgm3": lo,
                    "rho_hi_kgm3": hi,
                }
            )

        selected_source = "primary_near_bracket"
        selected_lo = np.nan
        selected_hi = np.nan
        selected_rho = np.nan
        selected_p = np.nan
        selected_h = np.nan
        selected_dPdrho = np.nan
        fallback_used = False
        error_note = ""

        try:
            if branch == "vapor":
                near_floor = cfg.rho_floor
                near_max = max(rho_prev, cfg.rho_floor * 1.01)
            else:
                near_floor = rho_anchor
                near_max = cfg.rho_cap
            lo, hi = _find_sign_change_bracket_near(
                t_k,
                float(p_target),
                w1,
                w2,
                float(rho_prev),
                float(near_floor),
                float(near_max),
            )
            selected_lo = float(lo)
            selected_hi = float(hi)
            selected_rho = _solve_rho_pressure_with_bracket(t_k, float(p_target), w1, w2, lo, hi)
        except Exception as exc:
            fallback_used = True
            selected_source = "fallback_solve_rho_mass_for_P"
            error_note = str(exc)
            rho_fb = solve_rho_mass_for_P(
                T_K=t_k,
                P_target_kPa=float(p_target),
                comp1_json="r1234ze.json",
                comp2_json="r227ea.json",
                w1=w1,
                w2=w2,
                interaction=default_interaction(),
                phase_hint=branch,
            )
            if rho_fb is not None and np.isfinite(rho_fb):
                selected_rho = float(rho_fb)

        if np.isfinite(selected_rho):
            st = compute_props("r1234ze.json", "r227ea.json", t_k, float(selected_rho), w1, w2, default_interaction())
            selected_p = float(st.p)
            selected_h = float(st.h)
            selected_dPdrho = float(_dPdrho_num_kpa_per_kgm3(t_k, float(selected_rho), w1, w2))
            rho_prev = float(selected_rho)

        selected_rank = np.nan
        if np.isfinite(selected_lo):
            for rank, i0 in enumerate(idx_pairs, start=1):
                lo_i = float(rho_grid[i0])
                hi_i = float(rho_grid[i0 + 1])
                if abs(lo_i - selected_lo) <= max(1e-12, 1e-6 * max(1.0, abs(selected_lo))) and abs(
                    hi_i - selected_hi
                ) <= max(1e-12, 1e-6 * max(1.0, abs(selected_hi))):
                    selected_rank = float(rank)
                    break

        summary_rows.append(
            {
                "branch": branch,
                "step_index": k,
                "P_target_kPa": float(p_target),
                "n_intervals": int(len(idx_pairs)),
                "selected_source": selected_source,
                "selected_interval_rank": selected_rank,
                "selected_bracket_lo_kgm3": selected_lo,
                "selected_bracket_hi_kgm3": selected_hi,
                "selected_rho_kgm3": selected_rho,
                "selected_p_kPa": selected_p,
                "selected_h_kJkg": selected_h,
                "selected_dPdrho_kPa_per_kgm3": selected_dPdrho,
                "fallback_used": bool(fallback_used),
                "error_note": error_note,
            }
        )

    return pd.DataFrame(summary_rows), pd.DataFrame(interval_rows)


def main() -> None:
    """
    Purpose
    -------
    Execute strict 70 C root-selection audit and export summary/interval CSVs.

    Inputs
    ------
    None (uses internal audit defaults for current debug protocol).

    Outputs
    -------
    summary CSV : diagnostics/root_audit_70C_summary_<date>.csv
    intervals CSV : diagnostics/root_audit_70C_intervals_<date>.csv
    stdout summary : counts of multi-root and fallback events

    Assumptions
    -----------
    - Uses current pseudo saturation CSV as the reference for Psat interpolation.

    Failure modes
    -------------
    - Missing saturation CSV will raise file-not-found error.

    References
    ----------
    - Root-selection diagnostics requested for 70 C isotherm branch debugging.

    Notes on numerical stability
    ----------------------------
    - No solver settings are changed; this is read-only diagnostics.
    """
    cfg = AuditConfig()
    sat_df = load_pseudopure_saturation_table(ROOT / "diagnostics/r515b_pseudopure_dome_20260303.csv")
    sat = interp_pseudopure_sat_at_T(sat_df, cfg.t_c)
    t_k = sat["T_K"]
    p_sat_kpa = sat["P_sat_kPa"]
    w2 = 1.0 - cfg.w1
    mw_mix = _mw_mix_kgmol(cfg.w1, w2)
    t_axis = sat_df["T_mid_K"].to_numpy(dtype=float)
    rho_f_molm3 = float(np.interp(t_k, t_axis, sat_df["rho_l_molm3"].to_numpy(dtype=float)))
    rho_g_molm3 = float(np.interp(t_k, t_axis, sat_df["rho_v_molm3"].to_numpy(dtype=float)))
    rho_f_mass = rho_f_molm3 * mw_mix
    rho_g_mass = rho_g_molm3 * mw_mix

    p_v_kpa = np.logspace(
        np.log10(cfg.pmin_bar * 100.0),
        np.log10(max(cfg.pmin_bar * 100.0, (p_sat_kpa - cfg.eps_kpa))),
        cfg.nv,
    )
    p_l_kpa = np.logspace(
        np.log10(min(cfg.pmax_bar * 100.0, (p_sat_kpa + cfg.eps_kpa))),
        np.log10(cfg.pmax_bar * 100.0),
        cfg.nl,
    )

    vapor_summary, vapor_intervals = _run_branch_audit(cfg, "vapor", t_k, p_v_kpa, rho_g_mass)
    liquid_summary, liquid_intervals = _run_branch_audit(cfg, "liquid", t_k, p_l_kpa, rho_f_mass)

    summary = pd.concat([vapor_summary, liquid_summary], ignore_index=True)
    intervals = pd.concat([vapor_intervals, liquid_intervals], ignore_index=True)

    stamp = datetime.now().strftime("%Y%m%d")
    out_summary = ROOT / f"diagnostics/root_audit_70C_summary_{stamp}.csv"
    out_intervals = ROOT / f"diagnostics/root_audit_70C_intervals_{stamp}.csv"
    summary.to_csv(out_summary, index=False)
    intervals.to_csv(out_intervals, index=False)

    multi = int((summary["n_intervals"] > 1).sum())
    fb = int(summary["fallback_used"].sum())
    liq_fb = int(liquid_summary["fallback_used"].sum())
    vap_fb = int(vapor_summary["fallback_used"].sum())
    print(f"T_C={cfg.t_c:.2f}, Psat_kPa={p_sat_kpa:.9f}")
    print(f"summary_csv={out_summary}")
    print(f"intervals_csv={out_intervals}")
    print(f"steps_total={len(summary)}, multi_interval_steps={multi}, fallback_steps={fb}")
    print(f"fallback_vapor={vap_fb}, fallback_liquid={liq_fb}")
    print("liquid_multi_interval_top5:")
    print(
        liquid_summary.loc[liquid_summary["n_intervals"] > 1, ["step_index", "P_target_kPa", "n_intervals", "selected_source", "selected_interval_rank"]]
        .head(5)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
