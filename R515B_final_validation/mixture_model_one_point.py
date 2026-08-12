#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mixture_helmholtz.py

Author: Shilpa Narasimhan

Support: Claude AI, Codex

Date (Re)-created: 08/12/2026

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
    Hold Bell (2023) binary interaction parameters for the asymmetric
    quadratic reducing functions Tred(x1), vred(x1) (PDF Eqs. 9-10). These
    four numbers are the ONLY fitted "mixing rule" parameters in the
    reducing-function step of the model; the departure-function
    coefficients (BELL_2023_DEP_COEFFS) are a separate, independently
    fitted set (paper Table 7).

    Inputs
    ------
    beta_T : float [unitless]
        Asymmetry parameter for the temperature reducing function
        (PDF Eq. 9). Enters via theta_T = 1/(beta_T^2*x1 + x2) and the
        cross critical temperature Tc12 = beta_T*gamma_T*sqrt(Tc1*Tc2)
        (PDF Eq. 5/7; paper Table 2, column "beta_T,12").
    beta_v : float [unitless]
        Asymmetry parameter for the volume reducing function (PDF Eq. 10).
        Enters via theta_v = 1/(beta_v^2*x1 + x2) and the cross critical
        volume vc12 (PDF Eq. 6/8; paper Table 2, column "beta_v,12").
    gamma_T : float [unitless]
        Scale parameter for the cross critical temperature Tc12
        (paper Table 2, column "gamma_T,12").
    gamma_v : float [unitless]
        Scale parameter for the cross critical volume vc12
        (paper Table 2, column "gamma_v,12").

    Outputs
    -------
    Bell2023PairParams [dataclass]
      Immutable parameter bundle for reducing-rule evaluations.

    Assumptions
    -----------
    - Parameters are fitted for a specific binary pair and not transferable
      (these four values are specific to R-1234ze(E)/R-227ea and must not
      be reused for a different fluid pair).

    Failure modes
    -------------
    - No runtime validation; incorrect values propagate silently to
      reducing calculations (see the CONFIRMED beta_v transcription-error
      note on BELL_2023_R1234ZE_R227EA below).

    References
    ----------
    - Bell (2023), J. Phys. Chem. Ref. Data 52(1), 013101,
      DOI 10.1063/5.0135368, Table 2.
    - Helmholtz_R515B_mixutre_math.pdf (user's reference doc), Eqs. 5-10.

    Notes on numerical stability
    ----------------------------
    - Not applicable; this class stores scalar constants.
    """
    beta_T: float   ## PDF Eq. 5/7; paper Table 2 column "beta_T,12" -- temperature-reducing asymmetry parameter
    beta_v: float   ## PDF Eq. 6/8; paper Table 2 column "beta_v,12" -- volume-reducing asymmetry parameter
    gamma_T: float  ## paper Table 2 column "gamma_T,12" -- scale factor on the cross critical temperature Tc12
    gamma_v: float  ## paper Table 2 column "gamma_v,12" -- scale factor on the cross critical volume vc12


# Verified against Bell (2023), J. Phys. Chem. Ref. Data 52(1), 013101
# (DOI 10.1063/5.0135368), Table 2, row R-1234ze(E)/R-227ea, via direct
# `pdftotext -layout` extraction of the user-uploaded paper PDF
# (013101_1_online_1.pdf) on 2026-08-12. Paper Table 2 prints:
#   beta_T,12 = 1.001247   gamma_T,12 = 0.989180
#   beta_v,12 = 0.999290   gamma_v,12 = 1.001581
#
# CONFIRMED TRANSCRIPTION ERROR: beta_v below is 0.99290, but the paper
# prints 0.999290 (missing a "9"; ~0.65% relative difference). This is a
# real numeric bug, not a display/rounding artifact -- confirmed via
# deterministic pdftotext extraction, not an AI read of the table. Not
# fixed here since it is a science/parameter change and this file is
# being hand-edited by the user (see PROJECT_CONTEXT.md breadcrumb dated
# 2026-08-12, "CONFIRMED: beta_v Transcription Error in
# BELL_2023_R1234ZE_R227EA", for full verification detail).
BELL_2023_R1234ZE_R227EA = Bell2023PairParams(
    beta_T=1.001247,
    beta_v=0.999290,   ## BUG (confirmed): paper Table 2 says 0.999290, not 0.99290 -- see note above
    gamma_T=0.989180,
    gamma_v=1.001581,
)

# Bell (2023) Table 7 departure-function coefficients for R-1234ze(E)/227ea:
# each row is (n_k, t_k, d_k, l_k) feeding
# alpha^r_dep = x1*x2 * sum_k n_k * tau^t_k * delta^d_k * exp(-delta^l_k)
# (PDF Eq. 14/16). Verified correct against the paper via pdftotext on
# 2026-08-12 -- no discrepancy found in these three rows.
BELL_2023_DEP_COEFFS: List[Tuple[float, float, float, float]] = [
    (-0.057178, 1.290298, 1.0, 1.0),  ## paper Table 7, row 1: n_1, t_1, d_1, l_1
    ( 0.031318, 0.038796, 2.0, 1.0),  ## paper Table 7, row 2: n_2, t_2, d_2, l_2
    (-0.027496, 2.640532, 3.0, 1.0),  ## paper Table 7, row 3: n_3, t_3, d_3, l_3
]


def bell2023_Tred_vred(
    x1: float, x2: float, Tc1: float, Tc2: float, vc1: float, vc2: float,
    params: Bell2023PairParams
) -> Tuple[float, float]:
    """
    Purpose
    -------
    Evaluate the Bell (2023) binary reducing functions (PDF Eqs. 9-10,
    after Lemmon & Jacobsen 2004's asymmetric quadratic mixing-rule form,
    later reused by GERG-2008). Tred and vred are the mixture-level
    reducing temperature/volume that tau = Tred/T and delta = rho_mol*vred
    are built from -- NOT the true mixture critical point.

    IMPORTANT: Tred, vred are mixing-rule CONSTRUCTS, not real physical
    critical properties of the mixture at the given composition. They only
    coincide with the true pure-fluid critical temperature/volume at the
    two composition endpoints x1=0 and x1=1 (where the cross/asymmetry
    terms vanish and Tred->Tc2 or Tc1, vred->vc2 or vc1 by construction).
    For any 0<x1<1, the actual mixture critical point (where it exists at
    all, for a partially-miscible or single-phase binary) is a separate
    computed quantity, not Tred/vred themselves.

    Inputs
    ------
    x1 : float [mol/mol]
        Mole fraction of component 1.
    x2 : float [mol/mol]
        Mole fraction of component 2 (expected x2 = 1 - x1).
    Tc1 : float [K]
        Real pure-fluid critical temperature of component 1.
    Tc2 : float [K]
        Real pure-fluid critical temperature of component 2.
    vc1 : float [m^3/mol]
        Real pure-fluid critical molar volume of component 1.
    vc2 : float [m^3/mol]
        Real pure-fluid critical molar volume of component 2.
    params : Bell2023PairParams [unitless]
        Fitted binary interaction parameters beta_T/beta_v/gamma_T/gamma_v
        for this fluid pair (paper Table 2).

    Outputs
    -------
    T_red : float [K]
        Mixture REDUCING temperature (PDF Eq. 9) -- used only to form
        tau = Tred/T; not the mixture's true critical temperature.
    v_red : float [m^3/mol]
        Mixture REDUCING molar volume (PDF Eq. 10) -- used only to form
        delta = rho_mol*vred; not the mixture's true critical volume.

    Assumptions
    -----------
    - x1 + x2 = 1 for physical composition.
    - Pair parameters correspond to the intended binary system (here,
      R-1234ze(E)/R-227ea via BELL_2023_R1234ZE_R227EA).

    Failure modes
    -------------
    - Division by zero if composition-dependent denominators
      (beta_T^2*x1 + x2) or (beta_v^2*x1 + x2) vanish -- not expected for
      physically reasonable beta_T, beta_v (both close to 1) over
      0<=x1<=1.
    - Non-physical inputs (negative Tc/vc, x1 outside [0,1]) can produce
      non-physical reducing values with no explicit guard here.

    References
    ----------
    - Bell (2023), J. Phys. Chem. Ref. Data 52(1), 013101,
      DOI 10.1063/5.0135368, Eqs. for Tred/vred (paper Sec. II) and
      Table 2 (fitted beta/gamma parameters).
    - Lemmon & Jacobsen (2004) for the origin of this asymmetric quadratic
      reducing-function functional form, later adopted by Kunz & Wagner
      (GERG-2008) and by Bell (2023).
    - Helmholtz_R515B_mixutre_math.pdf (user's reference doc), Eqs. 9-10.

    Notes on numerical stability
    ----------------------------
    - Algebraic operations are well-conditioned for typical compositions
      away from degenerate denominators.
    """
    beta_T, beta_v, gamma_T, gamma_v = params.beta_T, params.beta_v, params.gamma_T, params.gamma_v

    ## theta_T, theta_v: composition-dependent ASYMMETRY WEIGHTS on the
    ## cross (1-2) term in Tred/vred (PDF Eqs. 9-10). Because beta_T,
    ## beta_v != 1 in general, theta_T/theta_v are NOT symmetric in
    ## x1<->x2 -- this is what makes the reducing function "asymmetric
    ## quadratic" rather than a simple mole-fraction-squared mixing rule.
    theta_T = (x1 + x2) / ((beta_T ** 2) * x1 + x2)
    theta_v = (x1 + x2) / ((beta_v ** 2) * x1 + x2)

    ## Tc12, vc12: fitted CROSS (1-2) critical temperature/volume (PDF
    ## Eqs. 5-8; paper Table 2) -- not real physical properties of any
    ## single fluid, but fitted combining-rule outputs specific to this
    ## fluid pair, scaled by gamma_T/gamma_v off a geometric-mean (Tc) or
    ## Lorentz-Berthelot-style (vc) combining rule.
    Tc12 = beta_T * gamma_T * np.sqrt(Tc1 * Tc2)
    vc12 = beta_v * gamma_v * ((vc1 ** (1.0/3.0) + vc2 ** (1.0/3.0)) ** 3) / 8.0

    ## Tred, vred (PDF Eqs. 9-10): mixture REDUCING temperature/volume --
    ## quadratic in (x1, x2) with a cross term weighted by theta_T/theta_v.
    ## These equal Tc1/vc1 or Tc2/vc2 only at the pure-component limits
    ## x1=1 or x1=0; for 0<x1<1 they are mixing-rule constructs, NOT the
    ## mixture's actual critical temperature/volume (see docstring above).
    Tred = (x1**2) * Tc1 + (x2**2) * Tc2 + 2.0 * x1 * x2 * theta_T * Tc12
    vred = (x1**2) * vc1 + (x2**2) * vc2 + 2.0 * x1 * x2 * theta_v * vc12
    return float(Tred), float(vred)


def bell2023_departure_alphar(x1: float, x2: float, tau: float, delta: float) -> Tuple[float, float, float]:
    """
    Purpose
    -------
    Evaluate the Bell (2023) departure residual Helmholtz term and its
    first tau/delta derivatives for the R-1234ze(E)/227ea binary pair
    (PDF Eq. 14/16). This is the x1*x2-weighted correction term ADDED to
    the mole-fraction-weighted pure-fluid alphar sum to form the full
    mixture alphar: ar_mix = x1*ar1 + x2*ar2 + dalphar.

    Inputs
    ------
    x1 : float [mol/mol]
        Mole fraction of component 1.
    x2 : float [mol/mol]
        Mole fraction of component 2 (expected x2 = 1 - x1).
    tau : float [unitless]
        Mixture reduced inverse temperature (= Tred/T).
    delta : float [unitless]
        Mixture reduced density (= rho_mol*vred).

    Outputs
    -------
    dalphar : float [unitless]
        Departure alpha^r contribution = x1*x2*sum_k(...) (PDF Eq. 14/16).
    dalphar_tau : float [unitless]
        d(dalphar)/d(tau), holding x1, x2, delta fixed.
    dalphar_del : float [unitless]
        d(dalphar)/d(delta), holding x1, x2, tau fixed.

    Assumptions
    -----------
    - Coefficients in `BELL_2023_DEP_COEFFS` (n_k, t_k, d_k, l_k triples,
      paper Table 7) are valid for this pair only; do not reuse for a
      different fluid pair.

    Failure modes
    -------------
    - Overflow/underflow in exponential terms at extreme delta.

    References
    ----------
    - Bell (2023), J. Phys. Chem. Ref. Data 52(1), 013101,
      DOI 10.1063/5.0135368, Table 7 coefficients for R-1234ze(E)/227ea.
    - Helmholtz_R515B_mixutre_math.pdf (user's reference doc), Eq. 14/16.

    Notes on numerical stability
    ----------------------------
    - Exponential damping can underflow harmlessly for very large delta.
    """
    pref = x1 * x2   ## composition prefactor x1*x2 on the whole departure term (PDF Eq. 14/16) -- vanishes at both pure-component limits
    val = 0.0        ## accumulator for dalphar itself: sum_k n_k*tau^t_k*delta^d_k*exp(-delta^l_k)
    val_tau = 0.0    ## accumulator for d(val)/d(tau) (before applying the x1*x2 prefactor)
    val_del = 0.0    ## accumulator for d(val)/d(delta) (before applying the x1*x2 prefactor)

    for nk, tk, dk, lk in BELL_2023_DEP_COEFFS:  ## one (n_k, t_k, d_k, l_k) tuple per Table-7 row (3 rows for this pair)
        expv = np.exp(-(delta ** lk))  ## exp(-delta^l_k) damping factor for this term
        term = nk * (tau ** tk) * (delta ** dk) * expv  ## this term's contribution to val: n_k*tau^t_k*delta^d_k*exp(-delta^l_k)
        val += term
        val_tau += nk * tk * (tau ** (tk - 1.0)) * (delta ** dk) * expv  ## d(term)/d(tau) via power rule on tau^t_k
        exp_del = expv * (-(lk * (delta ** (lk - 1.0))))  ## d/d(delta)[exp(-delta^l_k)] via chain rule
        val_del += nk * (tau ** tk) * ((dk * (delta ** (dk - 1.0)) * expv) + (delta ** dk) * exp_del)  ## d(term)/d(delta) via product rule on delta^d_k * exp(-delta^l_k)

    return float(pref * val), float(pref * val_tau), float(pref * val_del)


def bell2023_departure_base(tau: float, delta: float) -> float:
    """
    Evaluate the base binary departure sum WITHOUT the composition
    prefactor x1*x2 -- i.e. this returns exactly `val` from
    `bell2023_departure_alphar` (same PDF Eq. 14/16 sum), computed
    independently here for callers that need the un-weighted sum on its
    own (e.g. composition-derivative formulas that need to differentiate
    the x1*x2 prefactor and this sum separately via the product rule).

    Inputs
    ------
    tau : float [unitless]
        Mixture reduced inverse temperature.
    delta : float [unitless]
        Mixture reduced density.

    Outputs
    -------
    base_val : float [unitless]
        Sum_k n_k * tau^t_k * delta^d_k * exp(-delta^l_k), using the same
        BELL_2023_DEP_COEFFS (paper Table 7) as bell2023_departure_alphar.
    """
    base = 0.0
    for nk, tk, dk, lk in BELL_2023_DEP_COEFFS:
        base += nk * (tau ** tk) * (delta ** dk) * np.exp(-(delta ** lk))
    return float(base)


# -----------------------------------------------------------------------
# DEAD CODE -- commented out 2026-08-12. Confirmed no call sites anywhere
# in this file (compute_pressure_enthalpy re-derives the same a0/ar
# assembly inline; _mixture_alpha_eval is the evaluator actually used by
# compute_table1_properties). Other files in this repo (mixture_pseudo_dome.py,
# mixture_true_vle_copy.py, mixture_vle_true_reference.py, mixture_pseudo_dome_copy.py,
# scripts/validate_r515a_reference_suite.py, scripts/audit_residual_helmholtz.py)
# import a same-named function from linear_model_codex.py, NOT from this
# file -- so commenting this out here does not affect any other module.
# Kept below (rather than deleted) as a possible refactor target -- see
# PROJECT_CONTEXT.md breadcrumb dated 2026-08-12 for the original
# "CURRENTLY UNUSED" finding this commenting-out implements.
# -----------------------------------------------------------------------
# def mixture_alpha0_alphar_derivs(
    # d1: Dict,
    # d2: Dict,
    # x1: float,
    # x2: float,
    # tau: float,
    # delta: float,
    # Tred: float,
    # rho_red_mol: float,
    # pair_key: str,
# ) -> Tuple[float, float, float, float, float]:
    # """
    # Purpose
    # -------
    # Evaluate mixture ideal/residual Helmholtz terms and first derivatives
    # with respect to mixture reduced variables (tau, delta) at fixed
    # composition, using the same Bell (2023) chain-rule mapping as
    # `_mixture_alpha_eval` (PDF Eqs. 12-16).
#
    # CURRENTLY UNUSED (checked 2026-08-12): this function has no call sites
    # elsewhere in this file. `compute_pressure_enthalpy` re-derives the same
    # a0/ar mixture assembly inline (its own tau1/tau2/delta1/delta2/c1/c2/
    # k1/k2/a0_mix/ar_mix/etc. block) rather than calling this function, and
    # `_mixture_alpha_eval` is the evaluator actually used by
    # `compute_table1_properties`. This function is effectively dead code as
    # of the current file state -- kept here (not deleted) since it is
    # equivalent in intent to `_mixture_alpha_eval` and may be a refactor
    # target, but any caller wiring it in should double check it matches
    # `_mixture_alpha_eval`'s behavior (notably: this function does NOT add
    # the ideal entropy-of-mixing term x1*ln(x1)+x2*ln(x2) to a0_mix, unlike
    # `_mixture_alpha_eval`'s a0_mix -- see Outputs below).
#
    # Inputs
    # ------
    # d1 : dict [unitless]
        # Parsed IDAES JSON for fluid 1.
    # d2 : dict [unitless]
        # Parsed IDAES JSON for fluid 2.
    # x1 : float [mol/mol]
        # Mole fraction of fluid 1.
    # x2 : float [mol/mol]
        # Mole fraction of fluid 2 (expected x2 = 1 - x1).
    # tau : float [unitless]
        # Mixture reduced inverse temperature (= Tred/T).
    # delta : float [unitless]
        # Mixture reduced density (= rho_mol*vred).
    # Tred : float [K]
        # Mixture reducing temperature (PDF Eq. 9), passed in by the caller
        # rather than recomputed here.
    # rho_red_mol : float [mol/m^3]
        # Mixture reducing molar density (= 1/vred).
    # pair_key : str [unitless]
        # Fluid-pair identifier, currently informational only (unused in the
        # body below beyond the `_ = pair_key` no-op).
#
    # Outputs
    # -------
    # a0_mix : float [unitless]
        # Mixture ideal Helmholtz term = x1*a01 + x2*a02 ONLY -- unlike
        # `_mixture_alpha_eval`'s a0_mix, this does NOT include the ideal
        # entropy-of-mixing term x1*ln(x1) + x2*ln(x2) (PDF Eq. 16). Callers
        # expecting the full Eq.-16 a0_mix must add that term themselves.
    # a0_tau_mix : float [unitless]
        # d(a0_mix)/d(tau) at fixed composition.
    # ar_mix : float [unitless]
        # Mixture residual Helmholtz term including departure (PDF Eq. 14).
    # ar_tau_mix : float [unitless]
        # d(ar_mix)/d(tau) at fixed composition.
    # ar_del_mix : float [unitless]
        # d(ar_mix)/d(delta) at fixed composition.
#
    # Assumptions
    # -----------
    # - Differentiation is performed at fixed composition.
    # - Reducing functions are treated as composition-only for this derivative
      # path (no explicit composition-coupling derivative terms).
    # - Pure-fluid derivative evaluators return derivatives with respect to each
      # pure fluid's own reduced variables (tau_i, delta_i).
#
    # Failure modes
    # -------------
    # - KeyError/ValueError/TypeError for malformed/missing JSON schema values.
    # - ZeroDivisionError if reducing or critical scales are invalid.
#
    # Implemented vs not implemented
    # ------------------------------
    # - Implemented: fixed-composition chain-rule mapping from
      # (tau_i, delta_i)-derivatives to (tau, delta)-derivatives and Bell
      # departure contribution in residual term.
    # - Not implemented: composition-derivative coupling terms beyond fixed
      # composition treatment; the ideal entropy-of-mixing term in a0_mix
      # (see Outputs note above).
#
    # References
    # ----------
    # - Chain rule mapping used:
      # tau_i = (Tc_i / Tred) * tau, delta_i = (rho_red / rho_c_i) * delta
      # d(alpha)/d(tau) = d(alpha)/d(tau_i) * (Tc_i / Tred)
      # d(alpha)/d(delta) = d(alpha)/d(delta_i) * (rho_red / rho_c_i)
    # - Helmholtz_R515B_mixutre_math.pdf (user's reference doc), Eqs. 12-16.
#
    # Notes on numerical stability
    # ----------------------------
    # - Conditioning follows the underlying pure-fluid derivative evaluators.
    # """
    # Tc1 = float(d1["basic"]["Tc"])                        ## fluid-1 critical temperature [K]
    # Tc2 = float(d2["basic"]["Tc"])                        ## fluid-2 critical temperature [K]
    # MW1 = mw_from_json(d1)                                ## fluid-1 molar mass [kg/mol]
    # MW2 = mw_from_json(d2)                                ## fluid-2 molar mass [kg/mol]
    # rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1           ## fluid-1 critical density [mol/m^3]
    # rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2           ## fluid-2 critical density [mol/m^3]
#
    # c1 = Tc1 / Tred   ## temperature-side chain-rule scale factor, c1 = Tc1/Tred (PDF Eqs. 12-13)
    # c2 = Tc2 / Tred   ## c2 = Tc2/Tred
    # k1 = rho_red_mol / rhoc1_mol   ## density-side chain-rule scale factor, k1 = rho_red/rhoc1
    # k2 = rho_red_mol / rhoc2_mol   ## k2 = rho_red/rhoc2
#
    # tau1 = c1 * tau       ## fluid-1's own reduced inverse temperature (= Tc1/T)
    # tau2 = c2 * tau       ## fluid-2's own reduced inverse temperature (= Tc2/T)
    # delta1 = k1 * delta   ## fluid-1's own reduced density
    # delta2 = k2 * delta   ## fluid-2's own reduced density
#
    # ## Getting Helmholtz parameters from IDAES json files for individual fluids
    # a01, a01_tau_i = alpha0_idaes_with_derivs(d1["eos"], tau1, delta1)              ## fluid-1 alpha0 and d(alpha0)/d(tau1)
    # a02, a02_tau_i = alpha0_idaes_with_derivs(d2["eos"], tau2, delta2)              ## fluid-2 alpha0 and d(alpha0)/d(tau2)
    # ar1, ar1_tau_i, ar1_del_i = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)   ## fluid-1 alphar, d/d(tau1), d/d(delta1)
    # ar2, ar2_tau_i, ar2_del_i = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)   ## fluid-2 alphar, d/d(tau2), d/d(delta2)
#
    # dar, dar_tau, dar_del = bell2023_departure_alphar(x1, x2, tau, delta)  ## Bell (2023) departure alphar and its d/d(tau), d/d(delta)
#
    # a0_mix = x1 * a01 + x2 * a02                              ## mole-fraction-weighted pure alpha0 sum ONLY -- no ideal entropy-of-mixing term here (see Outputs note)
    # a0_tau_mix = x1 * (a01_tau_i * c1) + x2 * (a02_tau_i * c2) ## mixture d(a0)/d(tau): rescale each pure fluid's d/d(tau_i) by c_i, weight by mole fraction
#
    # ar_mix = x1 * ar1 + x2 * ar2 + dar                                          ## mixture alphar (PDF Eq. 14): mole-fraction-weighted pure alphar sum + departure term
    # ar_tau_mix = x1 * (ar1_tau_i * c1) + x2 * (ar2_tau_i * c2) + dar_tau        ## mixture d(alphar)/d(tau), same chain-rule mapping + departure d/d(tau)
    # ar_del_mix = x1 * (ar1_del_i * k1) + x2 * (ar2_del_i * k2) + dar_del        ## mixture d(alphar)/d(delta), same chain-rule mapping + departure d/d(delta)
#
    # _ = pair_key  ## pair_key is currently informational only; not used in this derivative path
    # return (
        # float(a0_mix),
        # float(a0_tau_mix),
        # float(ar_mix),
        # float(ar_tau_mix),
        # float(ar_del_mix),
    # )


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

    ## DEBUG: print every suspected quantity in the delta / ar_del_mix chain,
    ## in the ORDER they're computed above, so a bad value is easy to spot.
    print("x1 =", x1, " x2 =", x2)
    print("MWmix =", MWmix, " rho_mol =", rho_mol)
    print("Tc1 =", Tc1, " Tc2 =", Tc2, " vc1 =", vc1, " vc2 =", vc2)
    print("Tred =", Tred, " vred =", vred)
    print("tau =", tau, " delta =", delta)
    print("k1 =", k1, " k2 =", k2, " c1 =", c1, " c2 =", c2)
    print("tau1 =", tau1, " tau2 =", tau2, " delta1 =", delta1, " delta2 =", delta2)
    print("ar1_del_i =", ar1_del_i, " ar2_del_i =", ar2_del_i)
    print("dar_del =", dar_del)
    print("ar_del_mix =", ar_del_mix)

    Z = 1.0 + delta * ar_del_mix
    print("Z =", Z)
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
    Compute mixture reduced variables (tau, delta) and their per-component
    mappings (tau1/delta1, tau2/delta2) at one (x1, T, rho_mol) state, via
    the Bell (2023) reducing functions (PDF Eqs. 4, 9-13). This is the
    shared "get everyone onto reduced coordinates" helper used by both
    `_mixture_alpha_eval` and (independently/inline) `compute_pressure_enthalpy`.

    Purpose
    -------
    - Evaluate the mixture reducing temperature/volume Tred, vred via
      `bell2023_Tred_vred` (PDF Eqs. 9-10).
    - Form the mixture's own reduced state tau = Tred/T, delta = rho_mol*vred
      (PDF Eq. 4).
    - Map that mixture reduced state onto each pure fluid's own reduced
      coordinates tau_i, delta_i via the chain-rule scale factors
      c_i = Tc_i/Tred, k_i = rho_red/rhoc_i (PDF Eqs. 12-13), so that each
      pure fluid's IDAES Helmholtz evaluator can be called at ITS OWN
      correct reduced state rather than at the mixture's.

    Inputs
    ------
    d1, d2 : dict [unitless]
        Parsed IDAES JSON dictionaries for components 1 and 2
        (R-1234ze(E) and R-227ea).
    x1 : float [mol/mol]
        Mole fraction of component 1. x2 = 1 - x1.
    T : float [K]
        Mixture temperature.
    rho_mol : float [mol/m^3]
        Mixture molar density.

    Outputs
    -------
    tau : float [unitless]
        Mixture reduced inverse temperature, tau = Tred/T (PDF Eq. 4).
    delta : float [unitless]
        Mixture reduced density, delta = rho_mol*vred (PDF Eq. 4).
    tau1, tau2 : float [unitless]
        Per-component reduced inverse temperatures, tau_i = (Tc_i/Tred)*tau
        = Tc_i/T (PDF Eq. 12).
    delta1, delta2 : float [unitless]
        Per-component reduced densities, delta_i = (rho_red/rhoc_i)*delta
        (PDF Eq. 13).
    Tred : float [K]
        Mixture REDUCING temperature (PDF Eq. 9) -- a mixing-rule
        construct, NOT the mixture's true critical temperature except at
        the pure-component limits x1 in {0, 1} (see bell2023_Tred_vred
        docstring for the full explanation).
    vred : float [m^3/mol]
        Mixture REDUCING molar volume (PDF Eq. 10) -- same caveat as Tred.

    Assumptions
    -----------
    - x1 strictly in (0, 1) for a physical binary mixture (x1 in {0, 1}
      degenerates to a pure fluid, for which this reducing-function
      machinery is unnecessary but not explicitly guarded against here).
    - BELL_2023_R1234ZE_R227EA parameters are appropriate for the (d1, d2)
      pair actually passed in; this function does not verify that d1/d2
      correspond to R-1234ze(E)/R-227ea.

    Failure modes
    -------------
    - T == 0 K -> division by zero forming tau.
    - vred == 0 -> division by zero forming rho_red_mol (not expected for
      physical vc1, vc2, but not explicitly guarded against).

    References
    ----------
    - Helmholtz_R515B_mixutre_math.pdf (user's reference doc), Eqs. 4,
      9-13.
    - Bell (2023), J. Phys. Chem. Ref. Data 52(1), 013101, for the
      underlying reducing-function form (see bell2023_Tred_vred).

    Notes on numerical stability
    ----------------------------
    - Well-conditioned for physical T, rho_mol, x1 away from the
      degenerate limits noted above.
    """
    x2 = 1.0 - x1
    MW1 = mw_from_json(d1)  ## fluid-1 molar mass [kg/mol]
    MW2 = mw_from_json(d2)  ## fluid-2 molar mass [kg/mol]
    Tc1 = float(d1["basic"]["Tc"])  ## fluid-1 critical temperature [K] (real pure-fluid property)
    Tc2 = float(d2["basic"]["Tc"])  ## fluid-2 critical temperature [K]
    rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1  ## fluid-1 critical density [mol/m^3]
    rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2  ## fluid-2 critical density [mol/m^3]
    vc1 = 1.0 / rhoc1_mol  ## fluid-1 critical molar volume [m^3/mol], input to Bell reducing functions
    vc2 = 1.0 / rhoc2_mol  ## fluid-2 critical molar volume [m^3/mol]

    ## Tred, vred: mixture REDUCING temperature/volume (PDF Eqs. 9-10) --
    ## mixing-rule constructs, not the mixture's real critical point (see
    ## docstring above). Everything below derives from these two values.
    Tred, vred = bell2023_Tred_vred(x1, x2, Tc1, Tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = Tred / T           ## mixture reduced inverse temperature (PDF Eq. 4)
    delta = rho_mol * vred   ## mixture reduced density (PDF Eq. 4)
    print("delta = rho_mol*vred =", delta)
    print("vred = ", vred)
    rho_red_mol = 1.0 / vred  ## mixture reducing molar DENSITY (reciprocal of reducing volume)

    ## c1, c2: temperature-side chain-rule scale factors (PDF Eq. 12-13),
    ## c_i = Tc_i/Tred. Used both to build tau_i and, elsewhere, to rescale
    ## d(alpha_i)/d(tau_i) into d(alpha_i)/d(tau).
    c1 = Tc1 / Tred
    c2 = Tc2 / Tred
    ## k1, k2: density-side chain-rule scale factors, k_i = rho_red/rhoc_i.
    k1 = rho_red_mol / rhoc1_mol
    k2 = rho_red_mol / rhoc2_mol
    tau1 = c1 * tau       ## fluid-1's own reduced inverse temperature (= Tc1/T)
    tau2 = c2 * tau       ## fluid-2's own reduced inverse temperature (= Tc2/T)
    delta1 = k1 * delta   ## fluid-1's own reduced density
    delta2 = k2 * delta   ## fluid-2's own reduced density
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

    Purpose
    -------
    - Convert the mixture-level (x1, T, rho_mol) state into per-component
      reduced states (tau_i, delta_i) via `_mixture_reduced_state`.
    - Call each pure fluid's IDAES ideal-gas (alpha0) and residual
      (alphar) Helmholtz evaluators at its own reduced state.
    - Call the Bell (2023) departure function at the mixture reduced
      state.
    - Assemble the mole-fraction-weighted sums (PDF Eq. 14/16) that define
      the mixture's alpha0 and alphar, plus their first tau/delta
      derivatives (needed for p, h, and as building blocks for the second
      derivatives computed elsewhere by finite differences).

    Inputs
    ------
    d1, d2 : dict [unitless]
        Parsed IDAES Helmholtz JSON dictionaries for fluid 1
        (R-1234ze(E)) and fluid 2 (R-227ea) respectively.
    x1 : float [mol/mol]
        Mole fraction of component 1 (fluid 1). x2 = 1 - x1.
    T : float [K]
        Mixture temperature.
    rho_mol : float [mol/m^3]
        Mixture molar density.

    Outputs
    -------
    out : dict [mixed]
        "tau", "delta" : float [unitless]
            Mixture reduced inverse temperature / reduced density
            (PDF Eq. 4), evaluated at the Bell (2023) Tred/vred.
        "tau1", "tau2" : float [unitless]
            Per-component reduced inverse temperatures (= Tc_i / T),
            passed through from `_mixture_reduced_state` for reuse by
            callers (e.g. fugacity/chemical-potential formulas that need
            tau_i directly).
        "a0_mix" : float [unitless]
            Mixture ideal-gas alpha0, INCLUDING the ideal entropy-of-mixing
            term x1*ln(x1) + x2*ln(x2) (PDF Eq. 16). This is alpha0 only —
            NOT an entropy. No entropy (s/R) is computed anywhere in this
            function; s_over_r is assembled separately downstream in
            `compute_table1_properties` from a0/ar and their tau
            derivatives.
        "ar_mix" : float [unitless]
            Mixture residual alpha^r = x1*ar1 + x2*ar2 + dar (PDF Eq. 14),
            i.e. the mole-fraction-weighted pure-fluid residual terms plus
            the Bell (2023) departure-function correction `dar`.
        "h0_over_rt" : float [unitless]
            Ideal-gas contribution to h/(R*T) via tau*d(alpha0)/d(tau)
            (Table-1 style, PDF Sec. on enthalpy), NOT alpha0 itself.
        "a01", "a02" : float [unitless]
            Pure-fluid ideal-gas alpha0 for fluid 1 / fluid 2, each
            evaluated at that fluid's own (tau_i, delta_i) via the IDAES
            external function for that fluid (analytic, not FD).
        "a01_tau_i", "a02_tau_i" : float [unitless]
            d(alpha0_i)/d(tau_i) for fluid 1 / fluid 2 (analytic, IDAES).
        "ar1_tau_i", "ar2_tau_i" : float [unitless]
            d(alphar_i)/d(tau_i) for fluid 1 / fluid 2 (analytic, IDAES).
        "ar_tau_mix", "ar_del_mix" : float [unitless]
            Mixture d(alphar)/d(tau) and d(alphar)/d(delta), assembled by
            chain-rule mapping each pure fluid's own tau_i/delta_i
            derivative back onto the mixture's tau/delta (via c_i, k_i;
            PDF Eqs. 12-13) and adding the departure-function tau/delta
            derivatives dar_tau/dar_del.

    Assumptions
    -----------
    - x1 is strictly between 0 and 1 (x1*ln(x1) diverges to -inf at the
      pure-component limits; callers should special-case x1 in {0, 1} if
      ever needed, though this module is written for the binary mixture
      case and does not currently guard against that).
    - d1, d2 use the same reference-state convention for alpha0 as encoded
      in their respective IDAES JSON files; no reference-state shift is
      applied here beyond what's already baked into alpha0_idaes_with_derivs.

    Failure modes
    -------------
    - x1 <= 0 or x1 >= 1 -> np.log(x1) or np.log(x2) is -inf/NaN.
    - T or rho_mol outside the range the IDAES ancillary/EOS correlations
      were fit over -> alpha0/alphar evaluators may return inaccurate or
      NaN values (no explicit bounds-checking is done here).

    References
    ----------
    - Helmholtz_R515B_mixutre_math.pdf (user's reference doc), Eqs. 4,
      14-16.
    - Bell (2023), J. Phys. Chem. Ref. Data 52(1), 013101, for the
      departure-function form evaluated via bell2023_departure_alphar.

    Notes
    -----
    - c1, c2, k1, k2 (temperature-/density-side chain-rule scale factors,
      PDF Eqs. 12-13) are recomputed locally here from _Tred/_vred rather
      than returned by `_mixture_reduced_state`, since that helper already
      applies them internally to build tau1/tau2/delta1/delta2 and only
      returns the resulting reduced states, not the scale factors
      themselves. Recomputing them here is a (deliberate) minor duplication
      needed because ar1_tau_i/ar2_tau_i/ar1_del_i/ar2_del_i (returned by
      the IDAES evaluators as d/d(tau_i), d/d(delta_i)) still need to be
      rescaled onto mixture tau/delta via the same c_i, k_i.
    """
    x2 = 1.0 - x1
    (
        tau,     ## mixture reduced inverse temperature (PDF Eq. 4), tau = Tred/T
        delta,   ## mixture reduced density (PDF Eq. 4), delta = rho_mol*vred
        tau1,    ## fluid-1 own reduced inverse temperature (= Tc1/T)
        tau2,    ## fluid-2 own reduced inverse temperature (= Tc2/T)
        delta1,  ## fluid-1 own reduced density
        delta2,  ## fluid-2 own reduced density
        _Tred,   ## mixture REDUCING temperature (PDF Eq. 9) -- a mixing-rule
                 ## construct, NOT the mixture's true critical temperature
                 ## except at the pure-component limits x1 in {0, 1}.
        _vred,   ## mixture REDUCING volume (PDF Eq. 10) -- same caveat as _Tred.
    ) = _mixture_reduced_state(d1, d2, x1, T, rho_mol)

    a01, a01_tau_i = alpha0_idaes_with_derivs(d1["eos"], tau1, delta1)      ## fluid-1 alpha0 and d(alpha0)/d(tau1), at fluid-1's own reduced state
    a02, a02_tau_i = alpha0_idaes_with_derivs(d2["eos"], tau2, delta2)      ## fluid-2 alpha0 and d(alpha0)/d(tau2)
    ar1, ar1_tau_i, ar1_del_i = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)  ## fluid-1 alphar, d/d(tau1), d/d(delta1)
    ar2, ar2_tau_i, ar2_del_i = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)  ## fluid-2 alphar, d/d(tau2), d/d(delta2)
    dar, dar_tau, dar_del = bell2023_departure_alphar(x1, x2, tau, delta)   ## Bell (2023) departure alphar and its d/d(tau), d/d(delta) (PDF Eq. 14/16)

    # Ideal mixture entropy-of-mixing contribution in alpha^0.
    a0_mix = x1 * a01 + x2 * a02 + x1 * np.log(x1) + x2 * np.log(x2)  ## mole-fraction-weighted pure alpha0 sum + ideal entropy-of-mixing term (PDF Eq. 16); this is alpha0, not an entropy
    # Table-1 consistent ideal tau contribution.
    h0_over_rt = 1.0 + x1 * tau1 * a01_tau_i + x2 * tau2 * a02_tau_i  ## ideal-gas h/(R*T) contribution: 1 + x1*tau1*d(a01)/d(tau1) + x2*tau2*d(a02)/d(tau2)

    # Chain-rule mapped residual derivatives wrt mixture tau/delta.
    MW1 = mw_from_json(d1)                                    ## fluid-1 molar mass [kg/mol]
    MW2 = mw_from_json(d2)                                    ## fluid-2 molar mass [kg/mol]
    rhoc1_mol = float(d1["basic"]["rhoc"]) / MW1               ## fluid-1 critical density [mol/m^3]
    rhoc2_mol = float(d2["basic"]["rhoc"]) / MW2               ## fluid-2 critical density [mol/m^3]
    rho_red_mol = 1.0 / _vred                                  ## mixture reducing molar density (reciprocal of reducing volume)
    c1 = float(d1["basic"]["Tc"]) / _Tred                      ## temperature-side chain-rule scale factor, c1 = Tc1/Tred (PDF Eq. 12-13)
    c2 = float(d2["basic"]["Tc"]) / _Tred                      ## c2 = Tc2/Tred
    k1 = rho_red_mol / rhoc1_mol                                ## density-side chain-rule scale factor, k1 = rho_red/rhoc1
    k2 = rho_red_mol / rhoc2_mol                                ## k2 = rho_red/rhoc2
    ar_tau_mix = x1 * (ar1_tau_i * c1) + x2 * (ar2_tau_i * c2) + dar_tau  ## mixture d(alphar)/d(tau): rescale each pure fluid's d/d(tau_i) by c_i, weight by mole fraction, add departure d/d(tau)
    ar_del_mix = x1 * (ar1_del_i * k1) + x2 * (ar2_del_i * k2) + dar_del  ## mixture d(alphar)/d(delta): rescale each pure fluid's d/d(delta_i) by k_i, weight by mole fraction, add departure d/d(delta)
    ar_mix = x1 * ar1 + x2 * ar2 + dar                          ## mixture alphar itself (PDF Eq. 14): mole-fraction-weighted pure alphar sum + departure term

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
    Analytic derivatives of the Bell (2023) binary reducing functions
    Tred(x1), vred(x1) with respect to x1, holding x2 = 1 - x1 (i.e. the
    total derivative along the binary composition line, NOT a partial
    derivative at fixed x2). These are the x1-derivatives of the same
    Tred/vred that `bell2023_Tred_vred` computes as values; this function
    differentiates PDF Eqs. 9-10 in closed form rather than by finite
    difference.

    Purpose
    -------
    Used wherever the mixture model needs d(Tred)/d(x1) or d(vred)/d(x1)
    directly (e.g. composition-derivative/fugacity-type formulas), as an
    analytic alternative to differencing `bell2023_Tred_vred` numerically.

    Inputs
    ------
    x1 : float [mol/mol]
        Mole fraction of component 1. x2 = 1 - x1.
    Tc1, Tc2 : float [K]
        Pure-component critical temperatures (real critical properties).
    vc1, vc2 : float [m^3/mol]
        Pure-component critical molar volumes (real critical properties).
    params : Bell2023PairParams [unitless]
        Fitted binary interaction parameters (beta_T, beta_v, gamma_T,
        gamma_v) for this fluid pair (PDF Eqs. 5-8; paper Table 2).

    Outputs
    -------
    dTred_dx1 : float [K]
        d(Tred)/d(x1) along x2 = 1 - x1.
    dvred_dx1 : float [m^3/mol]
        d(vred)/d(x1) along x2 = 1 - x1.

    Assumptions
    -----------
    - Same asymmetric-quadratic reducing-function form as
      `bell2023_Tred_vred` (PDF Eqs. 9-10, after Lemmon & Jacobsen 2004);
      this function must stay in lockstep with that one if the functional
      form ever changes.
    - beta_T, beta_v are assumed nonzero and DT, DV are assumed nonzero
      over the composition range of interest (no explicit guard here).

    Failure modes
    -------------
    - DT or DV -> 0 would blow up theta_T/theta_v and their derivatives;
      not expected for physically reasonable beta_T, beta_v values but not
      explicitly checked.

    References
    ----------
    - Helmholtz_R515B_mixutre_math.pdf (user's reference doc), Eqs. 9-10
      (values) differentiated here wrt x1.
    - Bell (2023), J. Phys. Chem. Ref. Data 52(1), 013101, Table 2 for the
      beta_T/beta_v/gamma_T/gamma_v parameter source (see also the
      CONFIRMED beta_v transcription-error note on
      BELL_2023_R1234ZE_R227EA).

    Notes
    -----
    - theta_T, theta_v, Tc12, vc12 here are recomputed locally rather than
      reused from `bell2023_Tred_vred`, since that function returns only
      the final Tred/vred values, not these intermediate quantities.
    """
    x2 = 1.0 - x1
    beta_T, beta_v, gamma_T, gamma_v = params.beta_T, params.beta_v, params.gamma_T, params.gamma_v

    DT = (beta_T ** 2) * x1 + x2   ## denominator of theta_T = 1/(beta_T^2*x1 + x2) (PDF Eq. 9 asymmetry weight)
    DV = (beta_v ** 2) * x1 + x2   ## denominator of theta_v = 1/(beta_v^2*x1 + x2) (PDF Eq. 10 asymmetry weight)
    theta_T = 1.0 / DT             ## composition-dependent asymmetry weight on the cross (12) term in Tred
    theta_v = 1.0 / DV             ## composition-dependent asymmetry weight on the cross (12) term in vred
    dtheta_T = -((beta_T ** 2) - 1.0) / (DT * DT)   ## d(theta_T)/d(x1), from d/dx1[1/DT] with dDT/dx1 = beta_T^2 - 1
    dtheta_v = -((beta_v ** 2) - 1.0) / (DV * DV)   ## d(theta_v)/d(x1)

    Tc12 = beta_T * gamma_T * np.sqrt(Tc1 * Tc2)    ## fitted cross (1-2) critical temperature (PDF Eq. 5/7; paper Table 2)
    vc12 = beta_v * gamma_v * ((vc1 ** (1.0 / 3.0) + vc2 ** (1.0 / 3.0)) ** 3) / 8.0  ## fitted cross (1-2) critical molar volume (PDF Eq. 6/8)

    ## d/dx1[x1*x2*theta] = (x2-x1)*theta + x1*x2*d(theta)/dx1, applied to
    ## the x1*x2*theta_T and x1*x2*theta_v cross-term factors in Tred/vred.
    dxx_theta_T = (x2 - x1) * theta_T + x1 * x2 * dtheta_T
    dxx_theta_v = (x2 - x1) * theta_v + x1 * x2 * dtheta_v

    ## d(Tred)/dx1 from differentiating Tred = x1^2*Tc1 + x2^2*Tc2 +
    ## 2*x1*x2*theta_T*Tc12 (PDF Eq. 9) term by term wrt x1 (x2 = 1-x1).
    dTred_dx1 = 2.0 * x1 * Tc1 - 2.0 * x2 * Tc2 + 2.0 * Tc12 * dxx_theta_T
    ## d(vred)/dx1, same pattern applied to vred = x1^2*vc1 + x2^2*vc2 +
    ## 2*x1*x2*theta_v*vc12 (PDF Eq. 10).
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
    - Central Finite-difference approximation; not exact analytic derivatives.
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
    d1 = load_idaes_helmholtz_json(fluid1) ## Loading fluid-1 properties
    d2 = load_idaes_helmholtz_json(fluid2) ## Loading fluid-2 properties
    x1 = x1_from_w1(d1, d2, w1) ## Mole fraction of fluid-1
    x2 = 1.0 - x1 ## Mole fraction of fluid-2

    MW1 = mw_from_json(d1)
    MW2 = mw_from_json(d2)
    MWmix = x1 * MW1 + x2 * MW2 ## Mixture molecular weight
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
        # This is the residual Helmholtz energy
        _tau, _delta, _tau1, _tau2, _d1, _d2, Tred, vred = _mixture_reduced_state(d1, d2, x1, T, rho_mol)
        _ = _tau, _delta, _tau1, _tau2, _d1, _d2
        rho_from_delta = delta_v / vred
        # Scale T to honor tau variation: tau = Tred / T_eval -> T_eval = Tred/tau_v.
        T_eval = Tred / tau_v
        return _mixture_alpha_eval(d1, d2, x1, T_eval, rho_from_delta)["ar_mix"]

    ar_tautau, ar_deldel, ar_taudel = _fd_2d(ar_func, tau, delta, rel=fd_rel)

    # Table-1 real-gas relations.
    Z = 1.0 + delta * ar_del
    print("Z")
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

    f1_Pa = x1 * rho_mol * R_u * T * np.exp(d_na_dn1) # Fugacity of fluid-1
    f2_Pa = x2 * rho_mol * R_u * T * np.exp(d_na_dn2) # Fugacity of fluid-2
    phi1 = f1_Pa / max(1e-300, x1 * p_Pa) ## coefficient of fugacity for fluid-1
    phi2 = f2_Pa / max(1e-300, x2 * p_Pa) ## coefficient of fugacity for fluid-2

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
