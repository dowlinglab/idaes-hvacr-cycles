#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Third model wrapper:
- Keep existing pressure path from linear_model_codex.
- Report both raw EOS enthalpy and Honeywell chart-basis enthalpy.

Default Honeywell anchor:
    T_ref = 273.15 K
    rho_ref = 1258.4 kg/m^3
    h_chart_ref = 200.0 kJ/kg
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

try:
    from scipy.optimize import brentq as _brentq
except Exception:
    _brentq = None

from linear_model_codex import (
    compute_pressure_enthalpy,
)


@dataclass(frozen=True)
class ChartReference:
    T_ref_K: float = 273.15
    rho_ref_kgm3: float = 1258.4
    p_ref_kPa: Optional[float] = None
    h_chart_ref_kJkg: float = 200.0


def raw_ph_from_w1(fluid1: str, fluid2: str, w1: float, T_K: float, rho_kgm3: float) -> Tuple[float, float]:
    """Return raw model p [kPa], h [kJ/kg] for (T, rho, w1)."""
    return compute_pressure_enthalpy(fluid1, fluid2, w1, T_K, rho_kgm3)


def _solve_rho_mass_for_pressure(
    fluid1: str,
    fluid2: str,
    w1: float,
    T_K: float,
    p_target_kPa: float,
    phase_hint: str = "liquid",
    rho_min: float = 1e-3,
    rho_max: float = 2e3,
    ngrid: int = 800,
) -> Optional[float]:
    """
    Solve rho_mass [kg/m^3] for p(T, rho)=p_target_kPa.
    phase_hint: "liquid" (high-density root) or "vapor" (low-density root).
    """
    rhos = np.logspace(np.log10(rho_min), np.log10(rho_max), int(ngrid))
    vals = np.array([raw_ph_from_w1(fluid1, fluid2, w1, T_K, r)[0] - p_target_kPa for r in rhos])

    brackets = []
    for i in range(len(rhos) - 1):
        f1 = vals[i]
        f2 = vals[i + 1]
        if np.isnan(f1) or np.isnan(f2):
            continue
        if f1 == 0.0:
            brackets.append((rhos[i], rhos[i]))
        elif f1 * f2 < 0.0:
            brackets.append((rhos[i], rhos[i + 1]))
    if not brackets:
        return None

    roots = []
    for a, b in brackets:
        if a == b:
            roots.append(float(a))
            continue
        if _brentq is not None:
            root = _brentq(lambda rr: raw_ph_from_w1(fluid1, fluid2, w1, T_K, rr)[0] - p_target_kPa, a, b, maxiter=250)
        else:
            fa = raw_ph_from_w1(fluid1, fluid2, w1, T_K, a)[0] - p_target_kPa
            fb = raw_ph_from_w1(fluid1, fluid2, w1, T_K, b)[0] - p_target_kPa
            for _ in range(250):
                m = 0.5 * (a + b)
                fm = raw_ph_from_w1(fluid1, fluid2, w1, T_K, m)[0] - p_target_kPa
                if fa * fm <= 0.0:
                    b, fb = m, fm
                else:
                    a, fa = m, fm
            root = 0.5 * (a + b)
        roots.append(float(root))

    roots = sorted(set(roots))
    if phase_hint == "vapor":
        return roots[0]
    return roots[-1]


def _reference_density(fluid1: str, fluid2: str, w1: float, ref: ChartReference) -> float:
    """
    Reference density policy:
    - If p_ref_kPa is provided, compute liquid-root density at (T_ref_K, p_ref_kPa).
    - Otherwise use rho_ref_kgm3 as provided.
    """
    if ref.p_ref_kPa is not None:
        rho_liq = _solve_rho_mass_for_pressure(
            fluid1=fluid1,
            fluid2=fluid2,
            w1=w1,
            T_K=ref.T_ref_K,
            p_target_kPa=ref.p_ref_kPa,
            phase_hint="liquid",
        )
        if rho_liq is None:
            raise RuntimeError("Could not solve liquid reference density from Tref/Pref.")
        return rho_liq
    return ref.rho_ref_kgm3


def honeywell_h_offset(fluid1: str, fluid2: str, w1: float, ref: ChartReference) -> Tuple[float, float]:
    """
    Compute h offset so chart-basis enthalpy matches reference:
        h_chart = h_raw + h_offset

    Returns:
        h_offset_kJkg, rho_ref_used_kgm3
    """
    rho_ref_used = _reference_density(fluid1, fluid2, w1, ref)
    _p_ref_kPa, h_raw_ref = raw_ph_from_w1(
        fluid1=fluid1,
        fluid2=fluid2,
        w1=w1,
        T_K=ref.T_ref_K,
        rho_kgm3=rho_ref_used,
    )
    return ref.h_chart_ref_kJkg - h_raw_ref, rho_ref_used


def compute_ph_with_chart_basis(
    fluid1: str,
    fluid2: str,
    w1: float,
    T_K: float,
    rho_kgm3: float,
    ref: ChartReference = ChartReference(),
) -> Tuple[float, float, float, float, float]:
    """
    Return:
        p_kPa, h_raw_kJkg, h_chart_kJkg, h_offset_kJkg, rho_ref_used_kgm3
    """
    p_kPa, h_raw_kJkg = raw_ph_from_w1(fluid1, fluid2, w1, T_K, rho_kgm3)
    h_offset_kJkg, rho_ref_used = honeywell_h_offset(fluid1, fluid2, w1, ref)
    h_chart_kJkg = h_raw_kJkg + h_offset_kJkg
    return p_kPa, h_raw_kJkg, h_chart_kJkg, h_offset_kJkg, rho_ref_used


def _cli() -> None:
    parser = argparse.ArgumentParser(
        description="Compute p and both raw/chart-basis h for R515x-style points."
    )
    parser.add_argument("--fluid1", default="r1234ze", help="Fluid 1 stem/name")
    parser.add_argument("--fluid2", default="r227ea", help="Fluid 2 stem/name")
    parser.add_argument("--w1", type=float, default=0.911, help="Mass fraction of fluid1")
    parser.add_argument("--T", type=float, required=True, help="Temperature [K]")
    parser.add_argument("--rho", type=float, required=True, help="Mass density [kg/m^3]")
    parser.add_argument("--Tref", type=float, default=273.15, help="Chart reference temperature [K]")
    parser.add_argument("--rhoref", type=float, default=1258.4, help="Chart reference density [kg/m^3]")
    parser.add_argument("--pref", type=float, default=None, help="Chart reference saturation pressure [kPa] at Tref; if provided, rhoref is solved (liquid root)")
    parser.add_argument("--href", type=float, default=200.0, help="Chart reference enthalpy [kJ/kg]")
    args = parser.parse_args()

    ref = ChartReference(
        T_ref_K=args.Tref,
        rho_ref_kgm3=args.rhoref,
        p_ref_kPa=args.pref,
        h_chart_ref_kJkg=args.href,
    )
    p_kPa, h_raw, h_chart, h_off, rho_ref_used = compute_ph_with_chart_basis(
        fluid1=args.fluid1,
        fluid2=args.fluid2,
        w1=args.w1,
        T_K=args.T,
        rho_kgm3=args.rho,
        ref=ref,
    )
    print(f"p_kPa={p_kPa:.9g}")
    print(f"h_raw_kJkg={h_raw:.9g}")
    print(f"h_chart_kJkg={h_chart:.9g}")
    print(f"h_offset_kJkg={h_off:.9g}")
    print(f"rho_ref_used_kgm3={rho_ref_used:.9g}")


if __name__ == "__main__":
    _cli()
