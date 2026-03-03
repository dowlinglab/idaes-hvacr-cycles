#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa
Technical Support: AI assistant
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: Provide a clean API and CLI wrapper for querying pressure and
enthalpy at a specified thermodynamic point (T, rho_mass, composition) using
the strict implementation in linear_model_codex.py. External composition input
uses mass fraction and is converted to internal mole fraction.
Dependencies: argparse, typing, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.1
"""

from __future__ import annotations

import argparse
from typing import Tuple

from linear_model_codex import (
    compute_pressure_enthalpy,
)


# === SECTION: Error Classification ===
# Rationale: Surface strict derivative-mapping gaps with an actionable message.
def _format_strictness_error_message(exc: NotImplementedError) -> str:
    """
    Purpose
    -------
    Build a clear user-facing error message for strict derivative mapping gaps.

    Inputs
    ------
    exc : NotImplementedError [unitless]
        Exception raised from strict mixture derivative mapping path.

    Outputs
    -------
    msg : str [unitless]
        Human-readable message identifying the current blocker.

    Assumptions
    -----------
    - The strict path in `linear_model_codex.py` raises `NotImplementedError`
      until chain-rule derivative mapping is implemented.

    Failure modes
    -------------
    - If `exc` has an empty message, this function still returns a generic
      blocker message.

    References
    ----------
    - PROJECT_CONTEXT.md

    Notes on numerical stability
    ----------------------------
    - Not applicable; string formatting only.
    """
    base = str(exc).strip()
    prefix = (
        "Blocked by strict derivative mapping in "
        "mixture_alpha0_alphar_derivs: "
    )
    if base:
        return f"{prefix}{base}"
    return (
        prefix
        + "required chain-rule transform for pure-fluid-to-mixture derivatives "
        "is not implemented."
    )


# === SECTION: Public API ===
# Rationale: Expose a single typed entry point returning engineering units.
def get_ph_kpa_kjkg(
    fluid1: str,
    fluid2: str,
    w1: float,
    T_K: float,
    rho_mass_kgm3: float,
) -> Tuple[float, float]:
    """
    Purpose
    -------
    Compute pressure and specific enthalpy for a binary mixture state point.

    Inputs
    ------
    fluid1 : str [unitless]
        Name/stem of fluid 1 JSON parameter file.
    fluid2 : str [unitless]
        Name/stem of fluid 2 JSON parameter file.
    w1 : float [kg/kg]
        Mass fraction of fluid 1 (must satisfy 0 < w1 < 1).
    T_K : float [K]
        Mixture temperature.
    rho_mass_kgm3 : float [kg/m^3]
        Mixture mass density.

    Outputs
    -------
    p_kPa : float [kPa]
        Mixture pressure.
    h_kJkg : float [kJ/kg]
        Mixture specific enthalpy.

    Assumptions
    -----------
    - Underlying strict model in `linear_model_codex.py` is used unchanged.
    - Point definition is `(T_K, rho_mass_kgm3, w1_mass)` with `w2_mass = 1 - w1`.
    - Internal model API accepts mass fraction `w1`.

    Failure modes
    -------------
    - Raises ValueError for invalid mass-fraction bounds.
    - Raises FileNotFoundError/KeyError if fluid JSON data are unavailable or
      malformed.
    - Raises NotImplementedError when strict derivative mapping in
      `mixture_alpha0_alphar_derivs` is missing.

    References
    ----------
    - linear_model_codex.compute_pressure_enthalpy
    - PROJECT_CONTEXT.md

    Notes on numerical stability
    ----------------------------
    - Delegates all numerical behavior to `compute_pressure_enthalpy`.
    """
    try:
        if not (0.0 < w1 < 1.0):
            raise ValueError("w1 must be in (0,1) as a mass fraction [kg/kg].")

        return compute_pressure_enthalpy(fluid1, fluid2, w1, T_K, rho_mass_kgm3)
    except NotImplementedError as exc:
        msg = _format_strictness_error_message(exc)
        raise NotImplementedError(msg) from exc


# === SECTION: CLI ===
# Rationale: Enable direct terminal queries for one state point.
def _cli() -> None:
    """
    Purpose
    -------
    Parse command-line inputs and print p,h in engineering units.

    Inputs
    ------
    CLI args [unitless]:
    --fluid1, --fluid2, --w1 [kg/kg], --T [K], --rho [kg/m^3]

    Outputs
    -------
    None [unitless]
        Writes results or a clear blocking error message to stdout/stderr.

    Assumptions
    -----------
    - User provides one fully specified thermodynamic point.

    Failure modes
    -------------
    - Exits with code 2 on argument parsing failure.
    - Exits with code 1 on strict derivative mapping NotImplementedError.

    References
    ----------
    - argparse module documentation
    - PROJECT_CONTEXT.md

    Notes on numerical stability
    ----------------------------
    - Not applicable; orchestration only.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Query pressure [kPa] and enthalpy [kJ/kg] at a point defined by "
            "(T [K], rho [kg/m^3], w1 [kg/kg])."
        )
    )
    parser.add_argument("--fluid1", required=True, help="Fluid 1 name/stem")
    parser.add_argument("--fluid2", required=True, help="Fluid 2 name/stem")
    parser.add_argument("--w1", required=True, type=float, help="w1 mass fraction [kg/kg]")
    parser.add_argument("--T", required=True, type=float, help="Temperature [K]")
    parser.add_argument("--rho", required=True, type=float, help="Mass density [kg/m^3]")
    args = parser.parse_args()

    try:
        p_kpa, h_kjkg = get_ph_kpa_kjkg(
            fluid1=args.fluid1,
            fluid2=args.fluid2,
            w1=args.w1,
            T_K=args.T,
            rho_mass_kgm3=args.rho,
        )
    except NotImplementedError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    print(f"p_kPa = {p_kpa:.9g}")
    print(f"h_kJkg = {h_kjkg:.9g}")


if __name__ == "__main__":
    _cli()
