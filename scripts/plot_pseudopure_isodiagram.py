#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/testing: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Build pseudo-pure P-h dome overlays and pseudo-pure isotherm
segments for R515B chart-style visualization.
Dependencies: numpy, pandas, matplotlib, pressure_validated_model
Context reference: PROJECT_CONTEXT.md

Pseudo-pure iso-diagram alignment decision (2026-03-03):
Diagram uses pseudo-pure saturation representation to match Honeywell datasheet
style for R-515B (Solstice N15).
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pressure_validated_model import compute_props, default_interaction, solve_rho_mass_for_P, load_json


def load_pseudopure_saturation_table(csv_path: Path) -> pd.DataFrame:
    """
    Load pseudo-pure saturation data used for chart-style dome construction.

    Thermodynamic basis
    -------------------
    Diagram uses pseudo-pure saturation endpoints (h_f, h_g, P_sat) as the
    saturation boundary representation for R515B chart alignment.

    Why no two-phase interior
    -------------------------
    Interior two-phase values are not traced from rho-sweeps in this plot mode.
    Two-phase is represented by saturation endpoint connectors only.

    References
    ----------
    - Honeywell Solstice N15 TDS (page-2 P-h chart style; azeotropic/zero-glide framing).
    - Lemmon, Huber, McLinden (2018), NIST REFPROP Documentation.

    Breadcrumb context
    ------------------
    Pseudo-pure iso-diagram alignment decision (2026-03-03).
    """
    df = pd.read_csv(csv_path).copy()
    if "status" in df.columns:
        df = df[df["status"] == "OK"].copy()
    df = df.sort_values("T_mid_K")
    need = {"T_mid_K", "P_kPa", "h_l_kJkg", "h_v_kJkg"}
    missing = need.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in pseudo saturation table: {sorted(missing)}")
    return df


def interp_pseudopure_sat_at_T(sat_df: pd.DataFrame, T_C: float) -> Dict[str, float]:
    """
    Interpolate pseudo-pure saturation endpoint state at a target temperature.

    Thermodynamic basis
    -------------------
    Uses pseudo-pure saturation curve P_sat(T), h_f(T), h_g(T) for chart-level
    plotting, rather than bubble/dew split endpoints.

    Why no two-phase interior
    -------------------------
    Connector is drawn at constant P_sat(T) from h_f(T) to h_g(T), avoiding
    interior metastable tracing.

    References
    ----------
    - Honeywell Solstice N15 TDS chart representation (page 2).
    - Span (2000), Helmholtz-EOS state framework in (T, rho).

    Breadcrumb context
    ------------------
    Pseudo-pure iso-diagram alignment decision (2026-03-03).
    """
    T_K = T_C + 273.15
    x = sat_df["T_mid_K"].to_numpy(dtype=float)
    return {
        "T_K": float(T_K),
        "P_sat_kPa": float(np.interp(T_K, x, sat_df["P_kPa"].to_numpy(dtype=float))),
        "h_f_kJkg": float(np.interp(T_K, x, sat_df["h_l_kJkg"].to_numpy(dtype=float))),
        "h_g_kJkg": float(np.interp(T_K, x, sat_df["h_v_kJkg"].to_numpy(dtype=float))),
        "rho_f_molm3": float(np.interp(T_K, x, sat_df["rho_l_molm3"].to_numpy(dtype=float))),
    }


def single_phase_points_at_TP_grid(
    T_K: float,
    P_kPa_grid: np.ndarray,
    phase_hint: str,
    w1: float,
    w2: float,
) -> List[Tuple[float, float, str]]:
    """
    Compute single-phase branch points at fixed T over pressure grid.

    Thermodynamic basis
    -------------------
    EOS states are solved as (T, P) with phase-selected rho roots:
    low-density root for vapor-like branch and high-density root for liquid-like branch.

    Why no two-phase interior
    -------------------------
    Pressure grids are defined strictly outside P_sat(T) +/- eps, so only
    single-phase points are produced.

    References
    ----------
    - Span (2000), Multiparameter Equations of State.
    - Lemmon et al. (2018), NIST REFPROP Documentation.

    Breadcrumb context
    ------------------
    Pseudo-pure iso-diagram alignment decision (2026-03-03).
    """
    interaction = default_interaction()
    out: List[Tuple[float, float, str]] = []
    for P_kPa in P_kPa_grid:
        rho = solve_rho_mass_for_P(
            T_K=T_K,
            P_target_kPa=float(P_kPa),
            comp1_json="r1234ze.json",
            comp2_json="r227ea.json",
            w1=w1,
            w2=w2,
            interaction=interaction,
            phase_hint=phase_hint,
        )
        if rho is None or not np.isfinite(rho):
            continue
        st = compute_props("r1234ze.json", "r227ea.json", T_K, float(rho), w1, w2, interaction)
        out.append((float(st.p * 1.0e-2), float(st.h), phase_hint))
    return out


def _mw_mix_kgmol(w1: float, w2: float) -> float:
    d1 = load_json("r1234ze.json")
    d2 = load_json("r227ea.json")
    mw1 = float(d1["basic"]["MW"]) * 1e-3
    mw2 = float(d2["basic"]["MW"]) * 1e-3
    return 1.0 / (w1 / mw1 + w2 / mw2)


def _dPdrho_num_kpa_per_kgm3(T_K: float, rho: float, w1: float, w2: float, eps_rel: float = 1e-6) -> float:
    interaction = default_interaction()
    dr = max(1e-8, eps_rel * max(float(rho), 1.0))
    p_plus = compute_props("r1234ze.json", "r227ea.json", T_K, rho + dr, w1, w2, interaction).p
    p_minus = compute_props("r1234ze.json", "r227ea.json", T_K, max(rho - dr, 1e-9), w1, w2, interaction).p
    return float((p_plus - p_minus) / (2.0 * dr))


def _solve_rho_pressure_with_bracket(
    T_K: float,
    P_target_kPa: float,
    w1: float,
    w2: float,
    lo: float,
    hi: float,
) -> float:
    interaction = default_interaction()

    def f(r):
        return compute_props("r1234ze.json", "r227ea.json", T_K, float(r), w1, w2, interaction).p - P_target_kPa

    f_lo = f(lo)
    f_hi = f(hi)
    if not np.isfinite(f_lo) or not np.isfinite(f_hi) or f_lo * f_hi > 0.0:
        raise ValueError("invalid bracket")
    return float(brentq(f, lo, hi, maxiter=300, xtol=1e-10, rtol=1e-10))


def _find_sign_change_bracket_near(
    T_K: float,
    P_target_kPa: float,
    w1: float,
    w2: float,
    rho_center: float,
    rho_floor: float,
    rho_max: float,
) -> Tuple[float, float]:
    interaction = default_interaction()

    def f(r):
        return compute_props("r1234ze.json", "r227ea.json", T_K, float(r), w1, w2, interaction).p - P_target_kPa

    factors = [0.98, 0.95, 0.92, 0.90, 0.88, 0.85]
    ups = [1.02, 1.05, 1.08, 1.12, 1.20, 1.35, 1.60, 2.00]
    for fl in factors:
        for fu in ups:
            lo = max(rho_floor, rho_center * fl)
            hi = min(rho_max, rho_center * fu)
            if hi <= lo:
                continue
            try:
                f_lo = f(lo)
                f_hi = f(hi)
                if np.isfinite(f_lo) and np.isfinite(f_hi) and f_lo * f_hi <= 0.0:
                    return float(lo), float(hi)
            except Exception:
                continue

    # last resort: scan [rho_floor, rho_max] and pick first sign change above rho_floor
    grid = np.logspace(np.log10(max(rho_floor, 1e-6)), np.log10(rho_max), 1200)
    vals = []
    for r in grid:
        try:
            vals.append(f(float(r)))
        except Exception:
            vals.append(np.nan)
    vals = np.asarray(vals, dtype=float)
    for i in range(len(grid) - 1):
        v1 = vals[i]
        v2 = vals[i + 1]
        if np.isfinite(v1) and np.isfinite(v2) and v1 * v2 <= 0.0:
            return float(grid[i]), float(grid[i + 1])
    raise ValueError("no sign-change bracket found")


def build_pseudopure_isotherm_segments(
    sat_df: pd.DataFrame,
    T_C: float,
    w1: float = 0.911,
    P_plot_min_bar: float = 1.0,
    P_plot_max_bar: float = 35.0,
    eps_kPa: float = 1.0,
    n_vapor: int = 80,
    n_liquid: int = 80,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Build pseudo-pure isotherm segments: vapor, connector, liquid.

    Thermodynamic basis
    -------------------
    Uses pseudo-pure saturation at fixed T:
      - vapor branch for P <= P_sat - eps
      - horizontal connector at P_sat from h_f to h_g
      - liquid branch for P >= P_sat + eps

    Why no two-phase interior
    -------------------------
    Interior two-phase region is represented by the saturation connector only.
    No interior rho-sweep points are used.

    References
    ----------
    - Honeywell Solstice N15 TDS (page-2 chart style target).
    - Lemmon et al. (2018), NIST REFPROP Documentation.
    - EES REFPROP interface note: straight saturated-endpoint connector.

    Breadcrumb context
    ------------------
    Pseudo-pure iso-diagram alignment decision (2026-03-03).
    """
    sat = interp_pseudopure_sat_at_T(sat_df, T_C)
    T_K = sat["T_K"]
    w2 = 1.0 - w1
    mw_mix = _mw_mix_kgmol(w1, w2)
    rho_f_mass = sat["rho_f_molm3"] * mw_mix
    P_sat_bar = sat["P_sat_kPa"] * 1.0e-2

    vapor_hi_bar = max(P_plot_min_bar, P_sat_bar - eps_kPa * 1.0e-2)
    liquid_lo_bar = min(P_plot_max_bar, P_sat_bar + eps_kPa * 1.0e-2)

    vapor_points: List[Tuple[float, float, str]] = []
    if vapor_hi_bar > P_plot_min_bar:
        P_v_kPa = np.logspace(np.log10(P_plot_min_bar * 100.0), np.log10(vapor_hi_bar * 100.0), int(n_vapor))
        vapor_points = single_phase_points_at_TP_grid(T_K, P_v_kPa, "vapor", w1, w2)

    liquid_points: List[Tuple[float, float, str]] = []
    if P_plot_max_bar > liquid_lo_bar:
        P_l_kPa = np.logspace(np.log10(liquid_lo_bar * 100.0), np.log10(P_plot_max_bar * 100.0), int(n_liquid))
        # Continuation-based liquid branch tracing anchored at rho_f(T).
        rho_prev = float(rho_f_mass)
        first_diag = None
        for i, Pk in enumerate(P_l_kPa):
            # First point: bracket around rho_f(T).
            if i == 0:
                try:
                    lo, hi = _find_sign_change_bracket_near(T_K, float(Pk), w1, w2, rho_prev, rho_prev, 2000.0)
                    rho = _solve_rho_pressure_with_bracket(T_K, float(Pk), w1, w2, lo, hi)
                    if rho <= 0.9 * rho_f_mass:
                        # mandatory re-solve tighter near rho_f(T)
                        lo, hi = _find_sign_change_bracket_near(T_K, float(Pk), w1, w2, rho_f_mass, rho_f_mass, 2000.0)
                        rho = _solve_rho_pressure_with_bracket(T_K, float(Pk), w1, w2, lo, hi)
                except Exception:
                    rho = solve_rho_mass_for_P(
                        T_K=T_K,
                        P_target_kPa=float(Pk),
                        comp1_json="r1234ze.json",
                        comp2_json="r227ea.json",
                        w1=w1,
                        w2=w2,
                        interaction=default_interaction(),
                        phase_hint="liquid",
                    )
                    if rho is None or not np.isfinite(rho):
                        continue
            else:
                try:
                    lo, hi = _find_sign_change_bracket_near(T_K, float(Pk), w1, w2, rho_prev, rho_f_mass, 2000.0)
                    rho = _solve_rho_pressure_with_bracket(T_K, float(Pk), w1, w2, lo, hi)
                except Exception:
                    rho = solve_rho_mass_for_P(
                        T_K=T_K,
                        P_target_kPa=float(Pk),
                        comp1_json="r1234ze.json",
                        comp2_json="r227ea.json",
                        w1=w1,
                        w2=w2,
                        interaction=default_interaction(),
                        phase_hint="liquid",
                    )
                    if rho is None or not np.isfinite(rho):
                        continue

            # Hard guard against vapor root capture for liquid branch.
            if rho < rho_f_mass:
                try:
                    lo, hi = _find_sign_change_bracket_near(T_K, float(Pk), w1, w2, max(rho_prev, rho_f_mass), rho_f_mass, 2000.0)
                    rho = _solve_rho_pressure_with_bracket(T_K, float(Pk), w1, w2, lo, hi)
                except Exception:
                    rho = solve_rho_mass_for_P(
                        T_K=T_K,
                        P_target_kPa=float(Pk),
                        comp1_json="r1234ze.json",
                        comp2_json="r227ea.json",
                        w1=w1,
                        w2=w2,
                        interaction=default_interaction(),
                        phase_hint="liquid",
                    )
                    if rho is None or not np.isfinite(rho):
                        continue
                if rho < rho_f_mass:
                    continue

            dPdrho = _dPdrho_num_kpa_per_kgm3(T_K, rho, w1, w2)
            if not np.isfinite(dPdrho) or dPdrho <= 0.0:
                continue

            st = compute_props("r1234ze.json", "r227ea.json", T_K, float(rho), w1, w2, default_interaction())
            liquid_points.append((float(st.p * 1.0e-2), float(st.h), "liquid"))
            rho_prev = float(rho)
            if i == 0:
                first_diag = (rho_f_mass, rho, float(st.h), dPdrho)

        if first_diag is not None:
            sat["rho_f_mass_kgm3"] = first_diag[0]
            sat["rho_first_liquid_kgm3"] = first_diag[1]
            sat["h_first_liquid_kJkg"] = first_diag[2]
            sat["dPdrho_first_liquid"] = first_diag[3]

    rows = []
    for p_bar, h_kJkg, _ in vapor_points:
        rows.append({"P_bar": p_bar, "h_kJkg": h_kJkg, "phase_flag": "vapor"})
    rows.append({"P_bar": P_sat_bar, "h_kJkg": sat["h_f_kJkg"], "phase_flag": "two_phase_connector"})
    rows.append({"P_bar": P_sat_bar, "h_kJkg": sat["h_g_kJkg"], "phase_flag": "two_phase_connector"})
    for p_bar, h_kJkg, _ in liquid_points:
        rows.append({"P_bar": p_bar, "h_kJkg": h_kJkg, "phase_flag": "liquid"})
    return pd.DataFrame(rows), sat


def main() -> None:
    """
    Generate pseudo-pure dome with pseudo-pure isotherm overlays and export artifacts.

    Thermodynamic basis
    -------------------
    Diagram uses pseudo-pure saturation representation to match Honeywell
    datasheet style for R-515B.

    Why no two-phase interior
    -------------------------
    Isotherm two-phase part is plotted as a constant-pressure connector at P_sat(T).

    References
    ----------
    - Honeywell Solstice N15 TDS (page-2 P-h chart target style).
    - Lemmon et al. (2018), NIST REFPROP Documentation.

    Breadcrumb context
    ------------------
    Pseudo-pure iso-diagram alignment decision (2026-03-03).
    """
    p = argparse.ArgumentParser(description="Pseudo-pure dome + isotherm overlay plotter.")
    p.add_argument("--sat-csv", default="diagnostics/r515b_pseudopure_dome_20260303.csv")
    p.add_argument("--temps-c", default="70")
    p.add_argument("--w1", type=float, default=0.911)
    p.add_argument("--Pmin-bar", type=float, default=1.0)
    p.add_argument("--Pmax-bar", type=float, default=35.0)
    p.add_argument("--eps-kpa", type=float, default=1.0)
    p.add_argument("--nv", type=int, default=80)
    p.add_argument("--nl", type=int, default=80)
    args = p.parse_args()

    sat_df = load_pseudopure_saturation_table(ROOT / args.sat_csv)
    temps = [float(x.strip()) for x in args.temps_c.split(",") if x.strip()]
    stamp = datetime.now().strftime("%Y%m%d")

    out_dir = ROOT / "diagnostics/pseudopure_iso"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fig = out_dir / "ph_dome_pseudopure_iso_overlays.png"
    out_iso_csv = out_dir / f"T_{int(round(temps[0]))}C_pseudopure_isotherm.csv"

    fig, ax = plt.subplots(figsize=(10.8, 6.6), dpi=190)
    ax.set_yscale("log")
    ax.set_xlim(150.0, 500.0)
    ax.set_ylim(args.Pmin_bar, args.Pmax_bar)
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.grid(True, which="both", alpha=0.30)

    # Pseudo-pure dome boundaries + fill
    pbar = sat_df["P_kPa"].to_numpy(dtype=float) * 1.0e-2
    hf = sat_df["h_l_kJkg"].to_numpy(dtype=float)
    hg = sat_df["h_v_kJkg"].to_numpy(dtype=float)
    ax.fill_betweenx(pbar, hf, hg, color="#f4d35e", alpha=0.18, label="Two-phase region")
    ax.plot(hf, pbar, color="#0b4f6c", lw=2.0, label="Sat. liquid boundary")
    ax.plot(hg, pbar, color="#c44536", lw=2.0, label="Sat. vapor boundary")

    # quality lines
    for q in np.arange(0.1, 1.0, 0.1):
        hq = (1.0 - q) * hf + q * hg
        ax.plot(hq, pbar, color="#6a4c93", lw=0.7, alpha=0.45)
    ax.plot([], [], color="#6a4c93", lw=1.0, alpha=0.7, label="Quality lines x=0.1..0.9")

    all_iso_rows = []
    for i, t_c in enumerate(temps):
        iso_df, sat = build_pseudopure_isotherm_segments(
            sat_df=sat_df,
            T_C=t_c,
            w1=args.w1,
            P_plot_min_bar=args.Pmin_bar,
            P_plot_max_bar=args.Pmax_bar,
            eps_kPa=args.eps_kpa,
            n_vapor=args.nv,
            n_liquid=args.nl,
        )
        color = plt.get_cmap("viridis")(i / max(1, len(temps) - 1))
        dv = iso_df[iso_df["phase_flag"] == "vapor"]
        dc = iso_df[iso_df["phase_flag"] == "two_phase_connector"]
        dl = iso_df[iso_df["phase_flag"] == "liquid"]
        if len(dv):
            ax.plot(dv["h_kJkg"], dv["P_bar"], color=color, lw=1.2)
        if len(dc) == 2:
            ax.plot(dc["h_kJkg"], dc["P_bar"], color=color, lw=1.2, linestyle="--")
        if len(dl):
            ax.plot(dl["h_kJkg"], dl["P_bar"], color=color, lw=1.2)
        ax.plot([], [], color=color, lw=1.5, label=f"Isotherm {t_c:.0f}C")

        for _, r in iso_df.iterrows():
            all_iso_rows.append({"T_C": t_c, "P_bar": float(r["P_bar"]), "h_kJkg": float(r["h_kJkg"]), "phase_flag": str(r["phase_flag"])})

        if abs(t_c - 70.0) < 1e-9:
            print(f"P_sat(70C)={sat['P_sat_kPa']*1e-2:.9f} bar")
            print(f"h_f(70C)={sat['h_f_kJkg']:.9f} kJ/kg, h_g(70C)={sat['h_g_kJkg']:.9f} kJ/kg")
            if "rho_f_mass_kgm3" in sat:
                print(f"rho_f(70C)={sat['rho_f_mass_kgm3']:.9f} kg/m^3")
                print(f"rho_first_liquid(70C)={sat['rho_first_liquid_kgm3']:.9f} kg/m^3")
                print(f"h_first_liquid(70C)={sat['h_first_liquid_kJkg']:.9f} kJ/kg")
                print(f"dPdrho_first_liquid(70C)={sat['dPdrho_first_liquid']:.9f} kPa/(kg/m^3)")
            if len(dv):
                print(f"vapor_Pmax_bar={dv['P_bar'].max():.9f} (target P_sat-eps={(sat['P_sat_kPa']-args.eps_kpa)*1e-2:.9f})")
            if len(dl):
                print(f"liquid_Pmin_bar={dl['P_bar'].min():.9f} (target P_sat+eps={(sat['P_sat_kPa']+args.eps_kpa)*1e-2:.9f})")
            if len(dc) == 2:
                print(f"connector_pressure_bar={dc['P_bar'].iloc[0]:.9f} (exact P_sat)")

        if i == 0:
            iso_df.to_csv(out_iso_csv, index=False)

    ax.set_title("R515B P-h Dome with Pseudo-pure Isotherm Overlays")
    ax.legend(loc="lower right", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")

    pd.DataFrame(all_iso_rows).to_csv(out_dir / f"pseudopure_isotherm_overlays_{stamp}.csv", index=False)
    print(f"saved_figure={out_fig}")
    print(f"saved_isotherm_csv={out_iso_csv}")


if __name__ == "__main__":
    main()
