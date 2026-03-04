#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Shared density bracketing utility with explicit mass/molar units
for vapor and liquid root searches.
Dependencies: dataclasses, numpy
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class DensityBracketBounds:
    """
    Purpose
    -------
    Store mass-density bracketing bounds and adaptation settings.

    Inputs
    ------
    rho_min_mass_kgm3 : float [kg/m^3]
    rho_split_mass_kgm3 : float [kg/m^3]
    rho_max_mass_kgm3 : float [kg/m^3]
    rho_max_cap_mass_kgm3 : float [kg/m^3]
    rhoc_scale_factor : float [unitless]

    Outputs
    -------
    DensityBracketBounds [dataclass]

    Assumptions
    -----------
    - rho_min < rho_split < rho_max.

    Failure modes
    -------------
    - Invalid ordering raises ValueError via validate().

    References
    ----------
    - Common engineering bracketing practice for two-phase root scans.

    Notes on numerical stability
    ----------------------------
    - Fixed bounds reduce root-search degeneracy and missed sign changes.
    """

    rho_min_mass_kgm3: float
    rho_split_mass_kgm3: float
    rho_max_mass_kgm3: float
    rho_max_cap_mass_kgm3: float
    rhoc_scale_factor: float

    def validate(self) -> None:
        if not (self.rho_min_mass_kgm3 > 0.0):
            raise ValueError("rho_min_mass_kgm3 must be > 0")
        if not (self.rho_min_mass_kgm3 < self.rho_split_mass_kgm3 < self.rho_max_mass_kgm3):
            raise ValueError("Require rho_min < rho_split < rho_max")


@dataclass(frozen=True)
class DensityBracketResult:
    """
    Purpose
    -------
    Return vapor/liquid brackets and diagnostics for logging.

    Inputs
    ------
    vapor_bracket_molm3 : tuple[float, float] | None [mol/m^3]
    liquid_bracket_molm3 : tuple[float, float] | None [mol/m^3]
    used_fallback : bool [unitless]
    bounds : DensityBracketBounds
    note : str [unitless]

    Outputs
    -------
    DensityBracketResult [dataclass]

    Assumptions
    -----------
    - Root function is continuous over local intervals.

    Failure modes
    -------------
    - Brackets may be None if no sign changes are detected.

    References
    ----------
    - Brent-type bracketing prerequisites.

    Notes on numerical stability
    ----------------------------
    - Captures fallback usage and bounds for reproducibility.
    """

    vapor_bracket_molm3: Tuple[float, float] | None
    liquid_bracket_molm3: Tuple[float, float] | None
    used_fallback: bool
    bounds: DensityBracketBounds
    note: str


def default_density_bounds_from_critical(
    rhoc_values_mass_kgm3: Tuple[float, ...] | list[float],
    rho_min_mass_kgm3: float = 1.0e-4,
    rho_split_mass_kgm3: float = 50.0,
    rho_max_cap_mass_kgm3: float = 2000.0,
    rhoc_scale_factor: float = 3.0,
) -> DensityBracketBounds:
    """
    Purpose
    -------
    Build adaptive default mass-density bounds from critical-density scale.

    Inputs
    ------
    rhoc_values_mass_kgm3 : tuple/list [kg/m^3]
    rho_min_mass_kgm3 : float [kg/m^3]
    rho_split_mass_kgm3 : float [kg/m^3]
    rho_max_cap_mass_kgm3 : float [kg/m^3]
    rhoc_scale_factor : float [unitless]

    Outputs
    -------
    DensityBracketBounds

    Assumptions
    -----------
    - At least one critical density is provided.

    Failure modes
    -------------
    - Raises ValueError if input list is empty or bounds are invalid.

    References
    ----------
    - Adaptive upper-bound heuristic based on critical density scale.

    Notes on numerical stability
    ----------------------------
    - Using min(cap, scale*max(rhoc)) prevents runaway high-density searches.
    """
    if len(rhoc_values_mass_kgm3) == 0:
        raise ValueError("rhoc_values_mass_kgm3 must be non-empty")
    rhoc_max = float(np.max(np.asarray(rhoc_values_mass_kgm3, dtype=float)))
    rho_max_mass_kgm3 = float(min(rho_max_cap_mass_kgm3, rhoc_scale_factor * rhoc_max))
    # Keep liquid window valid even if adaptive cap is low.
    if rho_max_mass_kgm3 <= rho_split_mass_kgm3:
        rho_max_mass_kgm3 = float(max(rho_split_mass_kgm3 * 1.2, rho_split_mass_kgm3 + 1.0))
    out = DensityBracketBounds(
        rho_min_mass_kgm3=float(rho_min_mass_kgm3),
        rho_split_mass_kgm3=float(rho_split_mass_kgm3),
        rho_max_mass_kgm3=float(rho_max_mass_kgm3),
        rho_max_cap_mass_kgm3=float(rho_max_cap_mass_kgm3),
        rhoc_scale_factor=float(rhoc_scale_factor),
    )
    out.validate()
    return out


def mass_to_molar_density(rho_mass_kgm3: float, mw_mix_kgmol: float) -> float:
    """
    Purpose
    -------
    Convert mass density to molar density.

    Inputs
    ------
    rho_mass_kgm3 : float [kg/m^3]
    mw_mix_kgmol : float [kg/mol]

    Outputs
    -------
    rho_molm3 : float [mol/m^3]

    Assumptions
    -----------
    - Positive mixture molecular weight.

    Failure modes
    -------------
    - Raises ValueError for non-positive molecular weight.

    References
    ----------
    - rho_mol = rho_mass / MW.

    Notes on numerical stability
    ----------------------------
    - Linear mapping; numerically stable for finite values.
    """
    if mw_mix_kgmol <= 0.0:
        raise ValueError("mw_mix_kgmol must be > 0")
    return float(rho_mass_kgm3 / mw_mix_kgmol)


def _find_sign_change_bracket(func: Callable[[float], float], grid: np.ndarray) -> Tuple[float, float] | None:
    vals = np.array([func(float(x)) for x in grid], dtype=float)
    for i in range(len(grid) - 1):
        v1 = vals[i]
        v2 = vals[i + 1]
        if not np.isfinite(v1) or not np.isfinite(v2):
            continue
        if v1 == 0.0:
            return float(grid[i]), float(grid[i])
        if v1 * v2 < 0.0:
            return float(grid[i]), float(grid[i + 1])
    return None


def _find_all_sign_change_brackets(func: Callable[[float], float], grid: np.ndarray) -> list[Tuple[float, float]]:
    vals = np.array([func(float(x)) for x in grid], dtype=float)
    out: list[Tuple[float, float]] = []
    for i in range(len(grid) - 1):
        v1 = vals[i]
        v2 = vals[i + 1]
        if not np.isfinite(v1) or not np.isfinite(v2):
            continue
        if v1 == 0.0:
            out.append((float(grid[i]), float(grid[i])))
        elif v1 * v2 < 0.0:
            out.append((float(grid[i]), float(grid[i + 1])))
    return out


def format_bounds_log(bounds: DensityBracketBounds) -> str:
    """
    Purpose
    -------
    Format bounds for consistent failure logging.

    Inputs
    ------
    bounds : DensityBracketBounds

    Outputs
    -------
    text : str [unitless]

    Assumptions
    -----------
    - bounds already validated.

    Failure modes
    -------------
    - None.

    References
    ----------
    - Internal reproducibility logging convention.

    Notes on numerical stability
    ----------------------------
    - Not applicable.
    """
    return (
        f"rho_min={bounds.rho_min_mass_kgm3} kg/m^3, "
        f"rho_split={bounds.rho_split_mass_kgm3} kg/m^3, "
        f"rho_max={bounds.rho_max_mass_kgm3} kg/m^3, "
        f"rho_max_cap={bounds.rho_max_cap_mass_kgm3} kg/m^3, "
        f"rhoc_scale={bounds.rhoc_scale_factor}"
    )


def bracket_vapor_liquid_density_roots(
    func_of_rho_molm3: Callable[[float], float],
    mw_mix_kgmol: float,
    bounds: DensityBracketBounds,
    vapor_points: int = 260,
    liquid_points: int = 320,
    fallback_points: int = 1400,
) -> DensityBracketResult:
    """
    Purpose
    -------
    Compute vapor/liquid density brackets using shared two-window strategy.

    Inputs
    ------
    func_of_rho_molm3 : callable [mol/m^3 -> residual]
    mw_mix_kgmol : float [kg/mol]
    bounds : DensityBracketBounds
    vapor_points : int [count]
    liquid_points : int [count]
    fallback_points : int [count]

    Outputs
    -------
    DensityBracketResult

    Assumptions
    -----------
    - Two-phase state should exhibit two sign changes over full density range.

    Failure modes
    -------------
    - Returns None brackets if no sign-change pair is found.

    References
    ----------
    - Two-window bracketing for vapor/liquid roots in saturation-like solves.

    Notes on numerical stability
    ----------------------------
    - Uses log-spaced scans plus adaptive fallback first/last sign changes.
    """
    bounds.validate()

    vap_lo = mass_to_molar_density(bounds.rho_min_mass_kgm3, mw_mix_kgmol)
    vap_hi = mass_to_molar_density(bounds.rho_split_mass_kgm3, mw_mix_kgmol)
    liq_lo = mass_to_molar_density(bounds.rho_split_mass_kgm3, mw_mix_kgmol)
    liq_hi = mass_to_molar_density(bounds.rho_max_mass_kgm3, mw_mix_kgmol)

    low_grid = np.logspace(np.log10(vap_lo), np.log10(vap_hi), int(vapor_points))
    high_grid = np.logspace(np.log10(liq_lo), np.log10(liq_hi), int(liquid_points))

    low_br = _find_sign_change_bracket(func_of_rho_molm3, low_grid)
    high_br = _find_sign_change_bracket(func_of_rho_molm3, high_grid)

    used_fallback = False
    note = "WINDOW_OK"

    if low_br is None or high_br is None:
        used_fallback = True
        full_grid = np.logspace(np.log10(vap_lo), np.log10(liq_hi), int(fallback_points))
        pairs = _find_all_sign_change_brackets(func_of_rho_molm3, full_grid)
        if len(pairs) >= 2:
            low_br = pairs[0]
            high_br = pairs[-1]
            note = "FALLBACK_FIRST_LAST_SIGN_CHANGE"
        else:
            note = "FAILED_NO_SIGN_CHANGE_PAIR"

    return DensityBracketResult(
        vapor_bracket_molm3=low_br,
        liquid_bracket_molm3=high_br,
        used_fallback=used_fallback,
        bounds=bounds,
        note=note,
    )
