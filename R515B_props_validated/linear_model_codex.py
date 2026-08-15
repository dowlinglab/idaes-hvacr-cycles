#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mixture_helmholtz.py

Author: Shilpa Narasimhan

Acknowledgment:
Code refactoring / scaffolding assistance provided by OpenAI Codex.
Final design, validation, and scientific responsibility belong to the human author.

Purpose
-------
Compute pressure and enthalpy for a binary refrigerant mixture given:
- two refrigerant names (IDAES Helmholtz JSON pure-fluid parameter files),
- mole fraction x1 (x2 = 1-x1),
- temperature T [K],
- mass density rho_mass [kg/m^3].

Outputs:
- pressure p [kPa]
- specific enthalpy h [kJ/kg]

Breadcrumbs
-----------
Last updated date: 2026-03-02
Validation status: runnable for single-point calls, but not yet benchmark-validated end-to-end.
Known limitations: composition-derivative coupling terms are not included in the current derivative treatment.
Next steps:
1. Add integration tests that validate p,h against known benchmark state points.
2. Add/verify composition-derivative coupling terms when differentiating beyond fixed-composition assumptions.
3. Preserve strict semantics; do not add silent approximations in this module.

IMPORTANT
---------
Do NOT paste citation placeholders like ":contentReference[oaicite:...]" into code.
They are not Python and will cause SyntaxError.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, List

import numpy as np

from idaes.models.properties.general_helmholtz import get_parameter_path


# -----------------------------
# Units / constants (SI basis)
# -----------------------------
R_u = 8.314462618  # J/mol/K
GMOL_TO_KG_PER_MOL = 1e-3   # g/mol -> kg/mol
PA_TO_KPA = 1e-3


# -----------------------------
# JSON loading
# -----------------------------
def load_idaes_helmholtz_json(fluid: str, parameter_path: str | Path | None = None) -> Dict:
    """
    Purpose
    -------
    Load an IDAES Helmholtz EOS JSON parameter file for one pure fluid.

    Inputs
    ------
    fluid : str [unitless]
        Refrigerant name or filename stem. If it does not end with ".json",
        ".json" is appended. Example: "r1234ze" -> "r1234ze.json".
    parameter_path : str | Path | None [unitless]
        Directory containing IDAES JSON files. If None, uses
        `idaes.models.properties.general_helmholtz.get_parameter_path()`.

    Outputs
    -------
    data : dict [unitless]
        Parsed JSON structure used by downstream EOS evaluators.

    Assumptions
    -----------
    - JSON schema follows IDAES Helmholtz conventions.
    - The target file exists under `parameter_path`.

    Failure modes
    -------------
    - FileNotFoundError if the JSON file is missing.
    - json.JSONDecodeError if file content is invalid JSON.
    - KeyError in downstream consumers if required schema keys are absent.

    JSON schema expectations
    ------------------------
    Required root keys include:
    - `basic` with `MW` [g/mol], `Tc` [K], `rhoc` [kg/m^3]
    - `eos` with ideal/residual coefficient dictionaries and type selectors.

    References
    ----------
    - IDAES general Helmholtz parameter file conventions.

    Notes on numerical stability
    ----------------------------
    - Not applicable; I/O and parsing only.
    """
    if parameter_path is None:
        parameter_path = get_parameter_path()
    parameter_path = Path(parameter_path)

    fname = fluid if fluid.endswith(".json") else f"{fluid}.json"
    fpath = parameter_path / fname
    with open(fpath, "r") as f:
        return json.load(f)


def mw_from_json(data: Dict) -> float:
    """
    Purpose
    -------
    Convert molecular weight from JSON units to SI molar mass units.

    Inputs
    ------
    data : dict [unitless]
        Parsed IDAES fluid JSON dictionary containing `basic.MW`.

    Outputs
    -------
    MW : float [kg/mol]
        Molar mass derived from `basic.MW` [g/mol].

    Assumptions
    -----------
    - `data["basic"]["MW"]` exists and is numeric in g/mol.

    Failure modes
    -------------
    - KeyError if `basic.MW` is missing.
    - ValueError/TypeError if `basic.MW` cannot be cast to float.

    References
    ----------
    - SI conversion: 1 g/mol = 1e-3 kg/mol.

    Notes on numerical stability
    ----------------------------
    - Stable linear unit conversion.
    """
    return float(data["basic"]["MW"]) * GMOL_TO_KG_PER_MOL


def x1_from_w1(data1: Dict, data2: Dict, w1: float) -> float:
    """
    Purpose
    -------
    Convert component-1 mass fraction to component-1 mole fraction.

    Inputs
    ------
    data1 : dict [unitless]
        Parsed IDAES JSON dictionary for fluid 1 (must contain `basic.MW`).
    data2 : dict [unitless]
        Parsed IDAES JSON dictionary for fluid 2 (must contain `basic.MW`).
    w1 : float [kg/kg]
        Mass fraction of fluid 1.

    Outputs
    -------
    x1 : float [mol/mol]
        Mole fraction of fluid 1.

    Assumptions
    -----------
    - `0 < w1 < 1`.
    - Molecular weights in JSON are valid and positive.

    Failure modes
    -------------
    - ValueError if `w1` is outside (0, 1).
    - KeyError/TypeError/ValueError for malformed JSON MW entries.

    References
    ----------
    - Standard basis conversion between mass and mole fractions.

    Notes on numerical stability
    ----------------------------
    - Stable for physical compositions strictly within (0, 1).
    """
    if not (0.0 < w1 < 1.0):
        raise ValueError("w1 must be in (0,1) as a mass fraction [kg/kg].")
    mw1 = mw_from_json(data1)
    mw2 = mw_from_json(data2)
    w2 = 1.0 - w1
    n1 = w1 / mw1
    n2 = w2 / mw2
    return float(n1 / (n1 + n2))


# -----------------------------
# Pure-fluid reduced variables
# -----------------------------
def pure_tau_delta_from_T_rhomol(data: Dict, T: float, rho_mol: float) -> Tuple[float, float]:
    """
    Purpose
    -------
    Compute pure-fluid reduced Helmholtz variables from (T, rho_mol).

    Inputs
    ------
    data : dict [unitless]
        Parsed IDAES fluid JSON dictionary with critical properties.
    T : float [K]
        Absolute temperature.
    rho_mol : float [mol/m^3]
        Molar density.

    Outputs
    -------
    tau : float [unitless]
        Reduced inverse temperature, Tc/T.
    delta : float [unitless]
        Reduced density, rho_mol/rhoc_mol.

    Assumptions
    -----------
    - JSON contains `basic.Tc`, `basic.rhoc`, and `basic.MW`.
    - T > 0 K and rhoc_mol > 0.

    Failure modes
    -------------
    - KeyError for missing critical-property keys.
    - ZeroDivisionError if T == 0 or MW == 0.
    - ValueError/TypeError for non-numeric values.

    JSON schema expectations
    ------------------------
    - `basic.Tc` [K]
    - `basic.rhoc` [kg/m^3]
    - `basic.MW` [g/mol]

    References
    ----------
    - Reduced-variable definitions used in Helmholtz EOS formulations.

    Notes on numerical stability
    ----------------------------
    - Conditioning degrades as T approaches 0 K.
    """
    Tc = float(data["basic"]["Tc"])
    rhoc_mass = float(data["basic"]["rhoc"])
    MW = mw_from_json(data)
    rhoc_mol = rhoc_mass / MW
    tau = Tc / T
    delta = rho_mol / rhoc_mol
    return tau, delta


# -----------------------------
# IDAES ideal part (alpha^0) evaluation with derivatives
# -----------------------------
def alpha0_idaes_with_derivs(eos: Dict, tau: float, delta: float) -> Tuple[float, float]:
    """
    Purpose
    -------
    Evaluate ideal Helmholtz term alpha^0 and d(alpha^0)/d(tau) from IDAES EOS
    coefficients.

    Inputs
    ------
    eos : dict [unitless]
        EOS sub-dictionary containing ideal-term coefficients/metadata:
        `n0`, `g0`, `phi_ideal_type`, `last_term_ideal`.
    tau : float [unitless]
        Reduced inverse temperature.
    delta : float [unitless]
        Reduced density.

    Outputs
    -------
    alpha0 : float [unitless]
        Ideal Helmholtz free-energy contribution.
    alpha0_tau : float [unitless]
        First derivative d(alpha0)/d(tau).

    Assumptions
    -----------
    - `tau > 0` and `delta > 0` for logarithm/exponential terms.
    - `phi_ideal_type` values map to implemented branches {1, 2, 3}.

    Failure modes
    -------------
    - KeyError for missing EOS coefficient keys.
    - FloatingPointError/RuntimeWarning domain issues for invalid tau/delta.
    - ValueError if `phi_ideal_type` is unsupported.

    JSON schema expectations
    ------------------------
    - Keys in `eos`: `n0`, `g0`, `phi_ideal_type`, `last_term_ideal`.

    References
    ----------
    - IDAES Helmholtz EOS ideal-term parameterization.

    Notes on numerical stability
    ----------------------------
    - Logarithms/exponentials may lose precision at extreme tau.
    """
    tau = max(float(tau), 1e-12)
    delta = max(float(delta), 1e-300)

    n0 = eos["n0"]
    g0 = eos["g0"]
    phi = int(eos["phi_ideal_type"])
    last = eos["last_term_ideal"]

    alpha0 = np.log(delta) + float(n0["1"]) + float(n0["2"]) * tau + float(n0["3"]) * np.log(tau)
    alpha0_tau = float(n0["2"]) + float(n0["3"]) * (1.0 / tau)

    def d_ln_1_minus_exp_neg(a: float, tau_: float) -> float:
        e = np.exp(-a * tau_)
        return (a * e) / (1.0 - e)

    if phi == 1:
        h = int(last)
        for k in range(4, h + 1):
            nk = float(n0[str(k)])
            ak = float(g0[str(k)])
            alpha0 += nk * np.log(1.0 - np.exp(-ak * tau))
            alpha0_tau += nk * d_ln_1_minus_exp_neg(ak, tau)

    elif phi == 2:
        h1 = int(last[0])
        h2 = int(last[1])
        for k in range(4, h1 + 1):
            nk = float(n0[str(k)])
            ak = float(g0[str(k)])
            alpha0 += nk * (tau ** ak)
            alpha0_tau += nk * ak * (tau ** (ak - 1.0))
        for k in range(h1 + 1, h2 + 1):
            nk = float(n0[str(k)])
            ak = float(g0[str(k)])
            alpha0 += nk * np.log(1.0 - np.exp(-ak * tau))
            alpha0_tau += nk * d_ln_1_minus_exp_neg(ak, tau)

    elif phi == 3:
        h = int(last)
        for k in range(4, h + 1):
            nk = float(n0[str(k)])
            ak = float(g0[str(k)])
            alpha0 += nk * (tau ** ak)
            alpha0_tau += nk * ak * (tau ** (ak - 1.0))

    else:
        raise ValueError(f"Unexpected phi_ideal_type={phi}")

    return float(alpha0), float(alpha0_tau)


# -----------------------------
# IDAES residual part (alpha^r) evaluation with derivatives
# -----------------------------
def alphar_idaes_with_derivs(eos: Dict, tau: float, delta: float) -> Tuple[float, float, float]:
    """
    Purpose
    -------
    Evaluate residual Helmholtz term alpha^r and first derivatives with respect
    to tau and delta from IDAES EOS coefficients.

    Inputs
    ------
    eos : dict [unitless]
        EOS sub-dictionary containing residual-term coefficients/metadata:
        `phi_residual_type`, `last_term_residual`, and coefficient maps
        (`n`, `d`, `t`, `c`, `a`, `b`, `e`, `g` as needed by branch).
    tau : float [unitless]
        Reduced inverse temperature.
    delta : float [unitless]
        Reduced density.

    Outputs
    -------
    alphar : float [unitless]
        Residual Helmholtz free-energy contribution.
    alphar_tau : float [unitless]
        First derivative d(alphar)/d(tau).
    alphar_del : float [unitless]
        First derivative d(alphar)/d(delta).

    Assumptions
    -----------
    - `tau > 0` and `delta > 0`.
    - `phi_residual_type` values map to implemented branches {1, 2, 3, 4}.

    Failure modes
    -------------
    - KeyError for missing coefficient keys.
    - ValueError if `phi_residual_type` is unsupported.
    - Floating-point overflow/underflow for extreme tau/delta.

    JSON schema expectations
    ------------------------
    - Keys in `eos`: `phi_residual_type`, `last_term_residual`, `n`, `d`, `t`.
    - Additional keys required by branch: `c`, `a`, `b`, `e`, `g`.

    References
    ----------
    - IDAES Helmholtz EOS residual-term parameterization.

    Notes on numerical stability
    ----------------------------
    - Exponential damping helps high-delta behavior but can underflow at large
      exponents; this is expected.
    """
    tau = max(float(tau), 1e-12)
    delta = max(float(delta), 1e-300)

    phi = int(eos["phi_residual_type"])
    hlist = eos["last_term_residual"]

    n = eos["n"]
    d = eos["d"]
    t = eos["t"]
    c = eos["c"]
    a = eos["a"]
    b = eos["b"]
    e = eos["e"]
    g = eos["g"]

    alphar = 0.0
    alphar_tau = 0.0
    alphar_del = 0.0

    def add_term(poly_coeff: float, di: float, ti: float, extra: float, extra_tau: float, extra_del: float):
        nonlocal alphar, alphar_tau, alphar_del
        base = poly_coeff * (delta ** di) * (tau ** ti) * extra
        alphar += base
        alphar_tau += poly_coeff * (delta ** di) * (tau ** ti) * extra_tau + base * (ti / tau)
        alphar_del += poly_coeff * (delta ** di) * (tau ** ti) * extra_del + base * (di / delta)

    if phi == 1:
        h1 = int(hlist[0]); h2 = int(hlist[1])
        for i in range(1, h1 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            add_term(ni, di, ti, extra=1.0, extra_tau=0.0, extra_del=0.0)
        for i in range(h1 + 1, h2 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            ci = float(c[str(i)])
            expv = np.exp(-(delta ** ci))
            extra_del = expv * (-(ci * (delta ** (ci - 1.0))))
            add_term(ni, di, ti, extra=expv, extra_tau=0.0, extra_del=extra_del)

    elif phi == 2:
        h1 = int(hlist[0]); h2 = int(hlist[1]); h3 = int(hlist[2])
        for i in range(1, h1 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            add_term(ni, di, ti, extra=1.0, extra_tau=0.0, extra_del=0.0)
        for i in range(h1 + 1, h2 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)]); ci = float(c[str(i)])
            expv = np.exp(-(delta ** ci))
            extra_del = expv * (-(ci * (delta ** (ci - 1.0))))
            add_term(ni, di, ti, extra=expv, extra_tau=0.0, extra_del=extra_del)
        for i in range(h2 + 1, h3 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            ai = float(a[str(i)]); bi = float(b[str(i)])
            ei = float(e[str(i)]); gi = float(g[str(i)])
            expv = np.exp(-ai * ((delta - ei) ** 2) - bi * ((tau - gi) ** 2))
            extra_tau = expv * (-(2.0 * bi * (tau - gi)))
            extra_del = expv * (-(2.0 * ai * (delta - ei)))
            add_term(ni, di, ti, extra=expv, extra_tau=extra_tau, extra_del=extra_del)

    elif phi == 3:
        h1 = int(hlist[0]); h2 = int(hlist[1]); h3 = int(hlist[2])
        for i in range(1, h1 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            add_term(ni, di, ti, extra=1.0, extra_tau=0.0, extra_del=0.0)
        for i in range(h1 + 1, h2 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)]); ci = float(c[str(i)])
            expv = np.exp(-(delta ** ci))
            extra_del = expv * (-(ci * (delta ** (ci - 1.0))))
            add_term(ni, di, ti, extra=expv, extra_tau=0.0, extra_del=extra_del)
        for i in range(h2 + 1, h3 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            ci = float(c[str(i)]); bi = float(b[str(i)])
            expd = np.exp(-(delta ** ci))
            expt = np.exp(-(tau ** bi))
            extra = expd * expt
            extra_tau = extra * (-(bi * (tau ** (bi - 1.0))))
            extra_del = extra * (-(ci * (delta ** (ci - 1.0))))
            add_term(ni, di, ti, extra=extra, extra_tau=extra_tau, extra_del=extra_del)

    elif phi == 4:
        h0 = int(hlist[0])
        for i in range(1, h0 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            add_term(ni, di, ti, extra=1.0, extra_tau=0.0, extra_del=0.0)

        m = len(hlist) - 1
        for j in range(1, m + 1):
            hjm1 = int(hlist[j - 1])
            hj = int(hlist[j])
            expv = np.exp(-(delta ** j))
            sumv = 0.0
            sumv_tau = 0.0
            sumv_del = 0.0
            for i in range(hjm1 + 1, hj + 1):
                ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
                term = ni * (delta ** di) * (tau ** ti)
                sumv += term
                sumv_tau += term * (ti / tau)
                sumv_del += term * (di / delta)

            alphar += expv * sumv
            alphar_tau += expv * sumv_tau
            expv_del = expv * (-(j * (delta ** (j - 1.0))))
            alphar_del += expv_del * sumv + expv * sumv_del

    else:
        raise ValueError(f"Unexpected phi_residual_type={phi}")

    return float(alphar), float(alphar_tau), float(alphar_del)


# -----------------------------
# Bell (2023) mixture reducing functions + departure (R1234ze(E)/227ea ONLY)
# -----------------------------
@dataclass(frozen=True)
class Bell2023PairParams:
    """
    Purpose
    -------
    Hold Bell (2023) binary interaction parameters for reducing functions.

    Inputs
    ------
    beta_T : float [unitless]
    beta_v : float [unitless]
    gamma_T : float [unitless]
    gamma_v : float [unitless]

    Outputs
    -------
    Bell2023PairParams [dataclass]
      Immutable parameter bundle for reducing-rule evaluations.

    Assumptions
    -----------
    - Parameters are fitted for a specific binary pair and not transferable.

    Failure modes
    -------------
    - No runtime validation; incorrect values propagate to reducing calculations.

    References
    ----------
    - Bell (2023), J. Phys. Chem. Ref. Data 52, 013101.

    Notes on numerical stability
    ----------------------------
    - Not applicable; this class stores scalar constants.
    """
    beta_T: float
    beta_v: float
    gamma_T: float
    gamma_v: float


BELL_2023_R1234ZE_R227EA = Bell2023PairParams(
    beta_T=1.001247,
    beta_v=0.99290,
    gamma_T=0.989180,
    gamma_v=1.001581,
)

# Table 7 coefficients for R-1234ze(E)/227ea
BELL_2023_DEP_COEFFS: List[Tuple[float, float, float, float]] = [
    (-0.057178, 1.290298, 1.0, 1.0),
    ( 0.031318, 0.038796, 2.0, 1.0),
    (-0.027496, 2.640532, 3.0, 1.0),
]


def bell2023_Tred_vred(
    x1: float, x2: float, Tc1: float, Tc2: float, vc1: float, vc2: float,
    params: Bell2023PairParams
) -> Tuple[float, float]:
    """
    Purpose
    -------
    Evaluate Bell (2023) binary reducing functions for reduced temperature and
    reduced molar volume.

    Inputs
    ------
    x1 : float [mol/mol]
        Mole fraction of component 1.
    x2 : float [mol/mol]
        Mole fraction of component 2.
    Tc1 : float [K]
        Critical temperature of component 1.
    Tc2 : float [K]
        Critical temperature of component 2.
    vc1 : float [m^3/mol]
        Critical molar volume of component 1.
    vc2 : float [m^3/mol]
        Critical molar volume of component 2.
    params : Bell2023PairParams [unitless]
        Binary interaction parameters beta/gamma for T and v reducing rules.

    Outputs
    -------
    T_red : float [K]
        Mixture reducing temperature.
    v_red : float [m^3/mol]
        Mixture reducing molar volume.

    Assumptions
    -----------
    - x1 + x2 = 1 for physical composition.
    - Pair parameters correspond to the intended binary system.

    Failure modes
    -------------
    - Division by zero if composition-dependent denominators vanish.
    - Non-physical inputs can produce non-physical reducing values.

    References
    ----------
    - Bell (2023) reducing-function correlations for refrigerant mixtures.

    Notes on numerical stability
    ----------------------------
    - Algebraic operations are well-conditioned for typical compositions away
      from degenerate denominators.
    """
    beta_T, beta_v, gamma_T, gamma_v = params.beta_T, params.beta_v, params.gamma_T, params.gamma_v

    theta_T = (x1 + x2) / ((beta_T ** 2) * x1 + x2)
    theta_v = (x1 + x2) / ((beta_v ** 2) * x1 + x2)

    Tc12 = beta_T * gamma_T * np.sqrt(Tc1 * Tc2)
    vc12 = beta_v * gamma_v * ((vc1 ** (1.0/3.0) + vc2 ** (1.0/3.0)) ** 3) / 8.0

    Tred = (x1**2) * Tc1 + (x2**2) * Tc2 + 2.0 * x1 * x2 * theta_T * Tc12
    vred = (x1**2) * vc1 + (x2**2) * vc2 + 2.0 * x1 * x2 * theta_v * vc12
    return float(Tred), float(vred)


def bell2023_departure_alphar(x1: float, x2: float, tau: float, delta: float) -> Tuple[float, float, float]:
    """
    Purpose
    -------
    Evaluate Bell (2023) departure residual Helmholtz term and derivatives for
    the R-1234ze(E)/227ea binary pair.

    Inputs
    ------
    x1 : float [mol/mol]
        Mole fraction of component 1.
    x2 : float [mol/mol]
        Mole fraction of component 2.
    tau : float [unitless]
        Mixture reduced inverse temperature.
    delta : float [unitless]
        Mixture reduced density.

    Outputs
    -------
    dalphar : float [unitless]
        Departure alpha^r contribution.
    dalphar_tau : float [unitless]
        d(dalphar)/d(tau).
    dalphar_del : float [unitless]
        d(dalphar)/d(delta).

    Assumptions
    -----------
    - Coefficients in `BELL_2023_DEP_COEFFS` are valid for this pair only.

    Failure modes
    -------------
    - Overflow/underflow in exponential terms at extreme delta.

    References
    ----------
    - Bell (2023), Table 7 coefficients for R-1234ze(E)/227ea.

    Notes on numerical stability
    ----------------------------
    - Exponential damping can underflow harmlessly for very large delta.
    """
    pref = x1 * x2
    val = 0.0
    val_tau = 0.0
    val_del = 0.0

    for nk, tk, dk, lk in BELL_2023_DEP_COEFFS:
        expv = np.exp(-(delta ** lk))
        term = nk * (tau ** tk) * (delta ** dk) * expv
        val += term
        val_tau += nk * tk * (tau ** (tk - 1.0)) * (delta ** dk) * expv
        exp_del = expv * (-(lk * (delta ** (lk - 1.0))))
        val_del += nk * (tau ** tk) * ((dk * (delta ** (dk - 1.0)) * expv) + (delta ** dk) * exp_del)

    return float(pref * val), float(pref * val_tau), float(pref * val_del)


def bell2023_departure_base(tau: float, delta: float) -> float:
    """
    Evaluate base binary departure sum without composition prefactor x1*x2.

    Inputs
    ------
    tau : float [unitless]
        Mixture reduced inverse temperature.
    delta : float [unitless]
        Mixture reduced density.

    Outputs
    -------
    base_val : float [unitless]
        Sum_k n_k * tau^t_k * delta^d_k * exp(-delta^l_k).
    """
    base = 0.0
    for nk, tk, dk, lk in BELL_2023_DEP_COEFFS:
        base += nk * (tau ** tk) * (delta ** dk) * np.exp(-(delta ** lk))
    return float(base)


def mixture_alpha0_alphar_derivs(
    d1: Dict,
    d2: Dict,
    x1: float,
    x2: float,
    tau: float,
    delta: float,
    Tred: float,
    rho_red_mol: float,
    pair_key: str,
) -> Tuple[float, float, float, float, float]:
    """
    Purpose
    -------
    Evaluate mixture ideal/residual Helmholtz terms and first derivatives with
    respect to mixture reduced variables (tau, delta) at fixed composition.

    Inputs
    ------
    d1 : dict [unitless]
        Parsed IDAES JSON for fluid 1.
    d2 : dict [unitless]
        Parsed IDAES JSON for fluid 2.
    x1 : float [mol/mol]
        Mole fraction of fluid 1.
    x2 : float [mol/mol]
        Mole fraction of fluid 2.
    tau : float [unitless]
        Mixture reduced inverse temperature.
    delta : float [unitless]
        Mixture reduced density.
    Tred : float [K]
        Mixture reducing temperature.
    rho_red_mol : float [mol/m^3]
        Mixture reducing molar density.
    pair_key : str [unitless]
        Fluid-pair identifier, currently informational.

    Outputs
    -------
    a0_mix : float [unitless]
        Mixture ideal Helmholtz term.
    a0_tau_mix : float [unitless]
        d(a0_mix)/d(tau) at fixed composition.
    ar_mix : float [unitless]
        Mixture residual Helmholtz term including departure.
    ar_tau_mix : float [unitless]
        d(ar_mix)/d(tau) at fixed composition.
    ar_del_mix : float [unitless]
        d(ar_mix)/d(delta) at fixed composition.

    Assumptions
    -----------
    - Differentiation is performed at fixed composition.
    - Reducing functions are treated as composition-only for this derivative
      path (no explicit composition-coupling derivative terms).
    - Pure-fluid derivative evaluators return derivatives with respect to each
      pure fluid's own reduced variables (tau_i, delta_i).

    Failure modes
    -------------
    - KeyError/ValueError/TypeError for malformed/missing JSON schema values.
    - ZeroDivisionError if reducing or critical scales are invalid.

    Implemented vs not implemented
    ------------------------------
    - Implemented: fixed-composition chain-rule mapping from
      (tau_i, delta_i)-derivatives to (tau, delta)-derivatives and Bell
      departure contribution in residual term.
    - Not implemented: composition-derivative coupling terms beyond fixed
      composition treatment.

    References
    ----------
    - Chain rule mapping used:
      tau_i = (Tc_i / Tred) * tau, delta_i = (rho_red / rho_c_i) * delta
      d(alpha)/d(tau) = d(alpha)/d(tau_i) * (Tc_i / Tred)
      d(alpha)/d(delta) = d(alpha)/d(delta_i) * (rho_red / rho_c_i)

    Notes on numerical stability
    ----------------------------
    - Conditioning follows the underlying pure-fluid derivative evaluators.
    """
    Tc1 = float(d1["basic"]["Tc"])
    Tc2 = float(d2["basic"]["Tc"])
    MW1 = mw_from_json(d1)
    MW2 = mw_from_json(d2)
    rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2

    c1 = Tc1 / Tred
    c2 = Tc2 / Tred
    k1 = rho_red_mol / rhoc1_mol
    k2 = rho_red_mol / rhoc2_mol

    tau1 = c1 * tau
    tau2 = c2 * tau
    delta1 = k1 * delta
    delta2 = k2 * delta

    a01, a01_tau_i = alpha0_idaes_with_derivs(d1["eos"], tau1, delta1)
    a02, a02_tau_i = alpha0_idaes_with_derivs(d2["eos"], tau2, delta2)
    ar1, ar1_tau_i, ar1_del_i = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)
    ar2, ar2_tau_i, ar2_del_i = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)

    dar, dar_tau, dar_del = bell2023_departure_alphar(x1, x2, tau, delta)

    a0_mix = x1 * a01 + x2 * a02
    a0_tau_mix = x1 * (a01_tau_i * c1) + x2 * (a02_tau_i * c2)

    ar_mix = x1 * ar1 + x2 * ar2 + dar
    ar_tau_mix = x1 * (ar1_tau_i * c1) + x2 * (ar2_tau_i * c2) + dar_tau
    ar_del_mix = x1 * (ar1_del_i * k1) + x2 * (ar2_del_i * k2) + dar_del

    _ = pair_key
    return (
        float(a0_mix),
        float(a0_tau_mix),
        float(ar_mix),
        float(ar_tau_mix),
        float(ar_del_mix),
    )


def compute_pressure_enthalpy(
    fluid1: str, fluid2: str,
    w1: float,
    T: float,
    rho_mass: float
) -> Tuple[float, float]:
    """
    Purpose
    -------
    Compute pressure and specific enthalpy for a binary mixture at a point
    defined by temperature, mass density, and composition.

    Inputs
    ------
    fluid1 : str [unitless]
        Fluid 1 JSON stem/name.
    fluid2 : str [unitless]
        Fluid 2 JSON stem/name.
    w1 : float [kg/kg]
        Mass fraction of fluid 1; must satisfy 0 < w1 < 1.
    T : float [K]
        Mixture temperature.
    rho_mass : float [kg/m^3]
        Mixture mass density.

    Outputs
    -------
    p_kPa : float [kPa]
        Mixture pressure.
    h_kJkg : float [kJ/kg]
        Mixture mass-specific enthalpy.

    Assumptions
    -----------
    - Bell (2023) reducing/departure parameterization is being used for
      R-1234ze(E)/227ea constants currently hardcoded in this module.
    - Input composition is binary on mass basis and converted internally to
      mole basis before reducing/departure evaluations.
    - Derivative mapping uses fixed-composition chain-rule transforms.

    Failure modes
    -------------
    - ValueError for out-of-range w1.
    - File/JSON schema errors when fluid files are missing or malformed.
    - May raise numeric/domain errors from EOS evaluations at non-physical
      reduced states.

    JSON schema expectations
    ------------------------
    For each fluid JSON:
    - `basic`: `MW`, `Tc`, `rhoc`
    - `eos`: ideal/residual parameter maps required by helper evaluators

    Implemented vs not implemented
    ------------------------------
    - Implemented: pure-fluid alpha0/alphar evaluators, fixed-composition
      chain-rule derivative mapping, Bell reducing rules, and p/h assembly.
    - Not implemented: explicit composition-coupling derivative extensions.

    References
    ----------
    - IDAES Helmholtz JSON fluid parameterization.
    - Bell (2023) reducing/departure framework for this binary.
    - Lemmon and Tillner-Roth mixture Helmholtz formulation:
      h/(RT) = h0/(RT) + tau*(d alpha_mix^r / d tau) + delta*(d alpha_mix^r / d delta),
      with h0/(RT) = 1 + sum_i x_i*tau_i*(d alpha_i^0 / d tau_i).

    Notes on numerical stability
    ----------------------------
    - Exponential/log operations in EOS terms can be sensitive at extreme
      reduced states.
    """
    d1 = load_idaes_helmholtz_json(fluid1)
    d2 = load_idaes_helmholtz_json(fluid2)
    x1 = x1_from_w1(d1, d2, w1)
    x2 = 1.0 - x1

    MW1 = mw_from_json(d1)
    MW2 = mw_from_json(d2)
    MWmix = x1 * MW1 + x2 * MW2  # kg/mol
    rho_mol = rho_mass / MWmix   # mol/m^3

    Tc1 = float(d1["basic"]["Tc"]); Tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    Tred, vred = bell2023_Tred_vred(x1, x2, Tc1, Tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = Tred / T
    delta = rho_mol * vred

    rho_red_mol = 1.0 / vred
    c1 = Tc1 / Tred
    c2 = Tc2 / Tred
    k1 = rho_red_mol / rhoc1_mol
    k2 = rho_red_mol / rhoc2_mol
    tau1 = c1 * tau
    tau2 = c2 * tau
    delta1 = k1 * delta
    delta2 = k2 * delta

    # === SECTION: Pure/Departure Derivatives for Lemmon-style Mixture h ===
    # Rationale: Assemble h/(RT) explicitly from component tau_i derivatives.
    _a01, a01_tau_i = alpha0_idaes_with_derivs(d1["eos"], tau1, delta1)
    _a02, a02_tau_i = alpha0_idaes_with_derivs(d2["eos"], tau2, delta2)
    _ar1, ar1_tau_i, ar1_del_i = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)
    _ar2, ar2_tau_i, ar2_del_i = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)
    _dar, dar_tau, dar_del = bell2023_departure_alphar(x1, x2, tau, delta)

    ar_del_mix = x1 * (ar1_del_i * k1) + x2 * (ar2_del_i * k2) + dar_del
    ar_tau_mix = x1 * (ar1_tau_i * c1) + x2 * (ar2_tau_i * c2) + dar_tau

    Z = 1.0 + delta * ar_del_mix
    p_Pa = rho_mol * R_u * T * Z
    p_kPa = p_Pa * PA_TO_KPA

    # Strict Lemmon-style enthalpy assembly:
    # h0/(RT) = 1 + sum_i x_i * tau_i * (d alpha0_i / d tau_i)
    H0_over_RT = 1.0 + x1 * tau1 * a01_tau_i + x2 * tau2 * a02_tau_i
    # h/(RT) = h0/(RT) + tau*d(alpha_mix^r)/d(tau) + delta*d(alpha_mix^r)/d(delta)
    H_over_RT = H0_over_RT + tau * ar_tau_mix + delta * ar_del_mix
    h_molar_Jmol = H_over_RT * R_u * T
    h_kJkg = (h_molar_Jmol / MWmix) * 1e-3

    return float(p_kPa), float(h_kJkg)


# === SECTION: Table-1 Mixture Thermodynamic Relations (Lemmon/Tillner-Roth) ===
# Rationale: Provide a full property evaluator including fugacity/mu_i-oriented terms.
def _mixture_reduced_state(
    d1: Dict,
    d2: Dict,
    x1: float,
    T: float,
    rho_mol: float,
) -> Tuple[float, float, float, float, float, float, float, float]:
    """
    Compute mixture reduced variables and component reduced mappings.

    Inputs
    ------
    d1, d2 : dict [unitless]
        Parsed IDAES JSON dictionaries for components 1 and 2.
    x1 : float [mol/mol]
        Mole fraction of component 1.
    T : float [K]
        Temperature.
    rho_mol : float [mol/m^3]
        Mixture molar density.

    Outputs
    -------
    tau, delta : float [unitless]
        Mixture reduced variables.
    tau1, tau2 : float [unitless]
        Component reduced inverse temperatures.
    delta1, delta2 : float [unitless]
        Component reduced densities.
    Tred : float [K]
        Mixture reducing temperature.
    vred : float [m^3/mol]
        Mixture reducing molar volume.
    """
    x2 = 1.0 - x1
    MW1 = mw_from_json(d1)
    MW2 = mw_from_json(d2)
    Tc1 = float(d1["basic"]["Tc"])
    Tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    Tred, vred = bell2023_Tred_vred(x1, x2, Tc1, Tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = Tred / T
    delta = rho_mol * vred
    rho_red_mol = 1.0 / vred

    c1 = Tc1 / Tred
    c2 = Tc2 / Tred
    k1 = rho_red_mol / rhoc1_mol
    k2 = rho_red_mol / rhoc2_mol
    tau1 = c1 * tau
    tau2 = c2 * tau
    delta1 = k1 * delta
    delta2 = k2 * delta
    return tau, delta, tau1, tau2, delta1, delta2, Tred, vred


def _mixture_alpha_eval(
    d1: Dict,
    d2: Dict,
    x1: float,
    T: float,
    rho_mol: float,
) -> Dict[str, float]:
    """
    Evaluate mixture alpha terms and first derivatives at one state.

    Inputs
    ------
    d1, d2 : dict [unitless]
        Parsed IDAES JSON dictionaries.
    x1 : float [mol/mol]
        Mole fraction of component 1.
    T : float [K]
        Temperature.
    rho_mol : float [mol/m^3]
        Molar density.

    Outputs
    -------
    out : dict [mixed]
        alpha/derivative terms required for Table-1 property formulas.
    """
    x2 = 1.0 - x1
    (
        tau,
        delta,
        tau1,
        tau2,
        delta1,
        delta2,
        _Tred,
        _vred,
    ) = _mixture_reduced_state(d1, d2, x1, T, rho_mol)

    a01, a01_tau_i = alpha0_idaes_with_derivs(d1["eos"], tau1, delta1)
    a02, a02_tau_i = alpha0_idaes_with_derivs(d2["eos"], tau2, delta2)
    ar1, ar1_tau_i, ar1_del_i = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)
    ar2, ar2_tau_i, ar2_del_i = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)
    dar, dar_tau, dar_del = bell2023_departure_alphar(x1, x2, tau, delta)

    # Ideal mixture entropy-of-mixing contribution in alpha^0.
    a0_mix = x1 * a01 + x2 * a02 + x1 * np.log(x1) + x2 * np.log(x2)
    # Table-1 consistent ideal tau contribution.
    h0_over_rt = 1.0 + x1 * tau1 * a01_tau_i + x2 * tau2 * a02_tau_i

    # Chain-rule mapped residual derivatives wrt mixture tau/delta.
    MW1 = mw_from_json(d1)
    MW2 = mw_from_json(d2)
    rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2
    rho_red_mol = 1.0 / _vred
    c1 = float(d1["basic"]["Tc"]) / _Tred
    c2 = float(d2["basic"]["Tc"]) / _Tred
    k1 = rho_red_mol / rhoc1_mol
    k2 = rho_red_mol / rhoc2_mol
    ar_tau_mix = x1 * (ar1_tau_i * c1) + x2 * (ar2_tau_i * c2) + dar_tau
    ar_del_mix = x1 * (ar1_del_i * k1) + x2 * (ar2_del_i * k2) + dar_del
    ar_mix = x1 * ar1 + x2 * ar2 + dar

    return {
        "tau": float(tau),
        "delta": float(delta),
        "tau1": float(tau1),
        "tau2": float(tau2),
        "a0_mix": float(a0_mix),
        "ar_mix": float(ar_mix),
        "h0_over_rt": float(h0_over_rt),
        "a01": float(a01),
        "a02": float(a02),
        "a01_tau_i": float(a01_tau_i),
        "a02_tau_i": float(a02_tau_i),
        "ar1_tau_i": float(ar1_tau_i),
        "ar2_tau_i": float(ar2_tau_i),
        "ar_tau_mix": float(ar_tau_mix),
        "ar_del_mix": float(ar_del_mix),
    }


def _bell2023_reducing_derivs_binary(
    x1: float,
    Tc1: float,
    Tc2: float,
    vc1: float,
    vc2: float,
    params: Bell2023PairParams,
) -> Tuple[float, float]:
    """
    Derivatives of binary Bell reducing functions wrt x1 for x2=1-x1.

    Inputs
    ------
    x1 : float [mol/mol]
        Mole fraction of component 1.
    Tc1, Tc2 : float [K]
        Component critical temperatures.
    vc1, vc2 : float [m^3/mol]
        Component critical molar volumes.
    params : Bell2023PairParams [unitless]
        Beta/gamma pair parameters.

    Outputs
    -------
    dTred_dx1 : float [K]
    dvred_dx1 : float [m^3/mol]
    """
    x2 = 1.0 - x1
    beta_T, beta_v, gamma_T, gamma_v = params.beta_T, params.beta_v, params.gamma_T, params.gamma_v

    DT = (beta_T ** 2) * x1 + x2
    DV = (beta_v ** 2) * x1 + x2
    theta_T = 1.0 / DT
    theta_v = 1.0 / DV
    dtheta_T = -((beta_T ** 2) - 1.0) / (DT * DT)
    dtheta_v = -((beta_v ** 2) - 1.0) / (DV * DV)

    Tc12 = beta_T * gamma_T * np.sqrt(Tc1 * Tc2)
    vc12 = beta_v * gamma_v * ((vc1 ** (1.0 / 3.0) + vc2 ** (1.0 / 3.0)) ** 3) / 8.0

    dxx_theta_T = (x2 - x1) * theta_T + x1 * x2 * dtheta_T
    dxx_theta_v = (x2 - x1) * theta_v + x1 * x2 * dtheta_v

    dTred_dx1 = 2.0 * x1 * Tc1 - 2.0 * x2 * Tc2 + 2.0 * Tc12 * dxx_theta_T
    dvred_dx1 = 2.0 * x1 * vc1 - 2.0 * x2 * vc2 + 2.0 * vc12 * dxx_theta_v
    return float(dTred_dx1), float(dvred_dx1)


def _fd_2d(func, tau: float, delta: float, rel: float = 1e-6) -> Tuple[float, float, float]:
    """
    Numerical 2D derivatives wrt (tau, delta): f_tau, f_deldel, f_tautau, f_taudel.

    Inputs
    ------
    func : callable [unitless]
        Function f(tau, delta) -> float.
    tau : float [unitless]
    delta : float [unitless]
    rel : float [unitless]
        Relative finite-difference step.

    Outputs
    -------
    f_tautau, f_deldel, f_taudel : float [unitless]
        Second derivatives.

    Assumptions
    -----------
    - Function is smooth near the target state.

    Failure modes
    -------------
    - Numerical noise if state is near singular surfaces.

    Notes on numerical stability
    ----------------------------
    - Finite-difference approximation; not exact analytic derivatives.
    """
    ht = max(1e-8, rel * abs(tau))
    hd = max(1e-8, rel * abs(delta))

    f00 = func(tau, delta)
    ftp = func(tau + ht, delta)
    ftm = func(max(1e-12, tau - ht), delta)
    fdp = func(tau, delta + hd)
    fdm = func(tau, max(1e-12, delta - hd))
    fpp = func(tau + ht, delta + hd)
    fpm = func(tau + ht, max(1e-12, delta - hd))
    fmp = func(max(1e-12, tau - ht), delta + hd)
    fmm = func(max(1e-12, tau - ht), max(1e-12, delta - hd))

    f_tautau = (ftp - 2.0 * f00 + ftm) / (ht * ht)
    f_deldel = (fdp - 2.0 * f00 + fdm) / (hd * hd)
    f_taudel = (fpp - fpm - fmp + fmm) / (4.0 * ht * hd)
    return float(f_tautau), float(f_deldel), float(f_taudel)


def _nares_from_n(
    d1: Dict,
    d2: Dict,
    T: float,
    V: float,
    n1: float,
    n2: float,
) -> float:
    """
    Evaluate n * alpha_mix^r for numerical fugacity derivative.

    Inputs
    ------
    d1, d2 : dict [unitless]
    T : float [K]
    V : float [m^3]
    n1, n2 : float [mol]

    Outputs
    -------
    n_ar : float [unitless * mol]
        n * alpha_mix^r.
    """
    n1 = max(1e-16, n1)
    n2 = max(1e-16, n2)
    n = n1 + n2
    x1 = n1 / n
    rho_mol = n / V
    vals = _mixture_alpha_eval(d1, d2, x1, T, rho_mol)
    return float(n * vals["ar_mix"])


def compute_table1_properties(
    fluid1: str,
    fluid2: str,
    w1: float,
    T: float,
    rho_mass: float,
    fd_rel: float = 1e-6,
) -> Dict[str, float]:
    """
    Compute Table-1 style mixture thermodynamic properties.

    Inputs
    ------
    fluid1, fluid2 : str [unitless]
        Fluid JSON stems.
    w1 : float [kg/kg]
        Mass fraction of component 1.
    T : float [K]
    rho_mass : float [kg/m^3]
    fd_rel : float [unitless]
        Relative finite-difference step for second/composition derivatives.

    Outputs
    -------
    props : dict [mixed units]
        Includes Z, p, h, u, s, cv, cp, speed_of_sound, fugacity terms.

    Assumptions
    -----------
    - Binary R1234ze(E)/R227ea Bell reducing/departure model.
    - Second derivatives and composition partial derivatives use finite
      differences in this implementation.

    Failure modes
    -------------
    - Raises on invalid composition/domain.
    - Numerical noise possible near critical/two-phase singular regions.
    """
    d1 = load_idaes_helmholtz_json(fluid1)
    d2 = load_idaes_helmholtz_json(fluid2)
    x1 = x1_from_w1(d1, d2, w1)
    x2 = 1.0 - x1

    MW1 = mw_from_json(d1)
    MW2 = mw_from_json(d2)
    MWmix = x1 * MW1 + x2 * MW2
    rho_mol = rho_mass / MWmix

    vals = _mixture_alpha_eval(d1, d2, x1, T, rho_mol)
    tau = vals["tau"]
    delta = vals["delta"]
    ar = vals["ar_mix"]
    ar_tau = vals["ar_tau_mix"]
    ar_del = vals["ar_del_mix"]
    h0_over_rt = vals["h0_over_rt"]
    a0_mix = vals["a0_mix"]

    def ar_func(tau_v: float, delta_v: float) -> float:
        # Reconstruct rho_mol from delta, keeping composition fixed.
        # delta = rho_mol * vred, with vred at fixed composition.
        _tau, _delta, _tau1, _tau2, _d1, _d2, Tred, vred = _mixture_reduced_state(d1, d2, x1, T, rho_mol)
        _ = _tau, _delta, _tau1, _tau2, _d1, _d2
        rho_from_delta = delta_v / vred
        # Scale T to honor tau variation: tau = Tred / T_eval -> T_eval = Tred/tau_v.
        T_eval = Tred / tau_v
        return _mixture_alpha_eval(d1, d2, x1, T_eval, rho_from_delta)["ar_mix"]

    ar_tautau, ar_deldel, ar_taudel = _fd_2d(ar_func, tau, delta, rel=fd_rel)

    # Table-1 real-gas relations.
    Z = 1.0 + delta * ar_del
    p_Pa = rho_mol * R_u * T * Z

    h_over_rt = h0_over_rt + tau * ar_tau + delta * ar_del
    u0_over_rt = h0_over_rt - 1.0
    u_over_rt = u0_over_rt + tau * ar_tau
    s0_over_r = h0_over_rt - a0_mix - 1.0
    s_over_r = s0_over_r + tau * ar_tau - ar

    # Ideal cv per Lemmon/Table-1 style:
    # cv0/R = -sum_i x_i * tau_i^2 * (d2 alpha_i^0 / d tau_i^2)
    _tau, _delta, tau1, tau2, delta1, delta2, _Tred, _vred = _mixture_reduced_state(d1, d2, x1, T, rho_mol)
    _ = _tau, _delta, _Tred, _vred

    def _a0_tautau_i(eos: Dict, tau_i: float, delta_i: float) -> float:
        h = max(1e-8, fd_rel * abs(tau_i))
        _, a_tau_p = alpha0_idaes_with_derivs(eos, tau_i + h, delta_i)
        _, a_tau_m = alpha0_idaes_with_derivs(eos, max(1e-12, tau_i - h), delta_i)
        return (a_tau_p - a_tau_m) / (2.0 * h)

    a01_tautau_i = _a0_tautau_i(d1["eos"], tau1, delta1)
    a02_tautau_i = _a0_tautau_i(d2["eos"], tau2, delta2)
    cv0_over_r = -(x1 * (tau1 * tau1) * a01_tautau_i + x2 * (tau2 * tau2) * a02_tautau_i)
    cv_over_r = cv0_over_r - (tau * tau) * ar_tautau
    cp_over_r = cv_over_r + ((1.0 + delta * ar_del - delta * tau * ar_taudel) ** 2) / (
        1.0 + 2.0 * delta * ar_del + (delta * delta) * ar_deldel
    )
    speed2_over_rt = (cp_over_r / cv_over_r) * (
        1.0 + 2.0 * delta * ar_del + (delta * delta) * ar_deldel
    )
    speed_sound = float(np.sqrt(max(0.0, speed2_over_rt * R_u * T / MWmix)))

    # Fugacity terms from analytic composition derivative of n*alpha^r.
    Tc1 = float(d1["basic"]["Tc"])
    Tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    Tred, vred = bell2023_Tred_vred(x1, x2, Tc1, Tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    dTred_dx1, dvred_dx1 = _bell2023_reducing_derivs_binary(
        x1, Tc1, Tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA
    )
    dtau_dx1 = dTred_dx1 / T
    ddelta_dx1 = rho_mol * dvred_dx1

    _tau_s, _delta_s, tau1_s, tau2_s, delta1_s, delta2_s, _Tred_s, _vred_s = _mixture_reduced_state(
        d1, d2, x1, T, rho_mol
    )
    _ = _tau_s, _delta_s, _Tred_s, _vred_s
    ar1, _, _ = alphar_idaes_with_derivs(d1["eos"], tau1_s, delta1_s)
    ar2, _, _ = alphar_idaes_with_derivs(d2["eos"], tau2_s, delta2_s)
    _, dep_tau, dep_del = bell2023_departure_alphar(x1, x2, tau, delta)
    dep_base = bell2023_departure_base(tau, delta)
    ddep_dx1 = (x2 - x1) * dep_base + dep_tau * dtau_dx1 + dep_del * ddelta_dx1

    # ar_x at fixed T and rho.
    dar_dx1 = (ar1 - ar2) + ddep_dx1
    # ar_rho at fixed T and x.
    dar_drho = ar_del * vred
    # n*ar derivatives at constant T,V,n_j for binary.
    d_na_dn1 = ar + rho_mol * dar_drho + x2 * dar_dx1
    d_na_dn2 = ar + rho_mol * dar_drho - x1 * dar_dx1

    f1_Pa = x1 * rho_mol * R_u * T * np.exp(d_na_dn1)
    f2_Pa = x2 * rho_mol * R_u * T * np.exp(d_na_dn2)
    phi1 = f1_Pa / max(1e-300, x1 * p_Pa)
    phi2 = f2_Pa / max(1e-300, x2 * p_Pa)

    # Chemical potential from fugacity.
    mu1_Jmol = R_u * T * np.log(max(1e-300, f1_Pa))
    mu2_Jmol = R_u * T * np.log(max(1e-300, f2_Pa))

    return {
        "x1_molmol": float(x1),
        "x2_molmol": float(x2),
        "rho_mol_molm3": float(rho_mol),
        "Z": float(Z),
        "p_Pa": float(p_Pa),
        "p_kPa": float(p_Pa * PA_TO_KPA),
        "h_over_rt": float(h_over_rt),
        "h_molar_Jmol": float(h_over_rt * R_u * T),
        "h_kJkg": float((h_over_rt * R_u * T / MWmix) * 1e-3),
        "u_over_rt": float(u_over_rt),
        "u_molar_Jmol": float(u_over_rt * R_u * T),
        "s_over_r": float(s_over_r),
        "s_molar_JmolK": float(s_over_r * R_u),
        "cv_over_r": float(cv_over_r),
        "cv_molar_JmolK": float(cv_over_r * R_u),
        "cp_over_r": float(cp_over_r),
        "cp_molar_JmolK": float(cp_over_r * R_u),
        "speed_sound_ms": float(speed_sound),
        "f1_Pa": float(f1_Pa),
        "f2_Pa": float(f2_Pa),
        "phi1": float(phi1),
        "phi2": float(phi2),
        "mu1_Jmol": float(mu1_Jmol),
        "mu2_Jmol": float(mu2_Jmol),
        "d_nares_dn1": float(d_na_dn1),
        "d_nares_dn2": float(d_na_dn2),
        "fd_rel": float(fd_rel),
    }


def _cli():
    """
    Purpose
    -------
    Command-line interface for one-off p,h computation calls.

    Inputs
    ------
    CLI args [unitless]:
    --fluid1, --fluid2, --w1 [kg/kg], --T [K], --rho [kg/m^3]

    Outputs
    -------
    None [unitless]
        Prints p_kPa and h_kJkg to stdout.

    Assumptions
    -----------
    - Arguments form a valid single-point query for the strict model.

    Failure modes
    -------------
    - Argument parsing errors terminate with exit code 2.
    - Propagates runtime exceptions from `compute_pressure_enthalpy`.

    References
    ----------
    - `argparse` standard library.

    Notes on numerical stability
    ----------------------------
    - Not applicable; CLI orchestration only.
    """
    p = argparse.ArgumentParser(description="Binary mixture Helmholtz p,h from IDAES JSON + Bell(2023) for R1234ze(E)/227ea.")
    p.add_argument("--fluid1", required=True, help="e.g. r1234ze")
    p.add_argument("--fluid2", required=True, help="e.g. r227ea")
    p.add_argument("--w1", required=True, type=float, help="mass fraction of fluid1 [kg/kg]")
    p.add_argument("--T", required=True, type=float, help="temperature [K]")
    p.add_argument("--rho", required=True, type=float, help="mass density [kg/m^3]")
    args = p.parse_args()

    pkPa, hkJkg = compute_pressure_enthalpy(args.fluid1, args.fluid2, args.w1, args.T, args.rho)
    print(f"p_kPa = {pkPa:.6g}")
    print(f"h_kJkg = {hkJkg:.6g}")


if __name__ == "__main__":
    _cli()
