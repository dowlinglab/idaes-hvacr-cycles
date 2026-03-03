#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: Compute true binary VLE envelopes (bubble/dew) for a
fixed-overall-composition refrigerant blend using Helmholtz EOS mixture model.
Dependencies: numpy, scipy, matplotlib, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-02
# Assumptions:
# - Binary VLE solved from P equality and component chemical-potential equality.
# - Chemical potentials are evaluated from fugacity using analytic composition
#   derivatives of n*alpha^r.
# - Bubble curve fixes liquid composition x=z; dew curve fixes vapor composition y=z.
# - This module is isolated from compute_pressure_enthalpy workflow.
# TODO: Add optional finite-difference mu_i cross-check for debugging.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from linear_model_codex import (
    BELL_2023_R1234ZE_R227EA,
    R_u,
    alphar_idaes_with_derivs,
    bell2023_departure_alphar,
    bell2023_departure_base,
    bell2023_Tred_vred,
    load_idaes_helmholtz_json,
    mixture_alpha0_alphar_derivs,
    mw_from_json,
)


EPS_X = 1e-12
RHO_MIN_MOLM3 = 1e-9
RHO_MAX_MOLM3 = 2.0e4
LSQ_DIFF_STEP = 1e-7
LSQ_MAX_NFEV = 800
LSQ_METHOD = "trf"
LSQ_X_SCALE = 1.0


@dataclass
class MixState:
    """
    Purpose
    -------
    Container for one phase state returned by mixture EOS evaluations.

    Inputs
    ------
    p_pa : float [Pa]
    h_jmol : float [J/mol]
    g_jmol : float [J/mol]
    rho_mol : float [mol/m^3]
    rho_mass : float [kg/m^3]
    x1 : float [mol/mol]

    Outputs
    -------
    MixState [dataclass]
      Immutable-style record of phase thermodynamic state.

    Assumptions
    -----------
    - Values are computed at thermodynamically meaningful states.

    Failure modes
    -------------
    - No internal validation; invalid values can be stored if provided.

    References
    ----------
    - Lemmon and Tillner-Roth Helmholtz mixture framework.

    Notes on numerical stability
    ----------------------------
    - Not applicable; this class only stores scalar outputs.
    """

    p_pa: float
    h_jmol: float
    g_jmol: float
    rho_mol: float
    rho_mass: float
    x1: float


def w1_to_x1(w1: float, mw1: float, mw2: float) -> float:
    """
    Purpose
    -------
    Convert component-1 mass fraction to mole fraction.

    Inputs
    ------
    w1 : float [kg/kg]
    mw1 : float [kg/mol]
    mw2 : float [kg/mol]

    Outputs
    -------
    x1 : float [mol/mol]

    Assumptions
    -----------
    - 0 <= w1 <= 1 and mw1,mw2 > 0.

    Failure modes
    -------------
    - Division by zero if mw1/mw2 invalid.
    - Returns non-physical values if w1 is outside [0,1].

    References
    ----------
    - Standard mass-to-mole fraction conversion relation.

    Notes on numerical stability
    ----------------------------
    - Stable for positive molecular weights and bounded mass fractions.
    """
    w2 = 1.0 - w1
    n1 = w1 / mw1
    n2 = w2 / mw2
    return float(n1 / (n1 + n2))


def _clip_x(x1: float) -> float:
    """
    Purpose
    -------
    Clamp composition into an open interval to avoid log singularities.

    Inputs
    ------
    x1 : float [mol/mol]

    Outputs
    -------
    x1_clipped : float [mol/mol]

    Assumptions
    -----------
    - Composition is intended to be near [0,1].

    Failure modes
    -------------
    - No explicit error; silently clips out-of-range inputs.

    References
    ----------
    - Numerical guard for log(x) terms in mixture Helmholtz ideal mixing.

    Notes on numerical stability
    ----------------------------
    - Prevents overflow/NaN from log(0) in entropy-of-mixing terms.
    """
    return float(min(max(x1, EPS_X), 1.0 - EPS_X))


def _mix_alpha_and_derivs(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> Tuple[float, float, float]:
    """
    Evaluate mixture alpha and reduced derivatives at fixed composition.

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]

    Outputs
    -------
    alpha_mix : float [unitless]
    alpha_tau_mix : float [unitless]
    ar_del_mix : float [unitless]

    Assumptions
    -----------
    - d1 and d2 match expected IDAES Helmholtz JSON schema.
    - Reducing/departure parameters are for R1234ze(E)/R227ea.

    Failure modes
    -------------
    - Raises KeyError/ValueError for malformed JSON data.
    - Numeric overflow/underflow possible at extreme reduced states.

    References
    ----------
    - Bell (2023) reducing and departure formulation.
    - Lemmon and Tillner-Roth Helmholtz mixture identities.

    Notes on numerical stability
    ----------------------------
    - Uses composition clipping and delegated EOS derivative guards.
    """
    x1 = _clip_x(x1)
    x2 = 1.0 - x1

    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k
    delta = rho_mol * vred
    rho_red_mol = 1.0 / vred

    a0, a0_tau, ar, ar_tau, ar_del = mixture_alpha0_alphar_derivs(
        d1=d1,
        d2=d2,
        x1=x1,
        x2=x2,
        tau=tau,
        delta=delta,
        Tred=tred,
        rho_red_mol=rho_red_mol,
        pair_key="r1234ze|r227ea",
    )

    # Add ideal mixing entropy contribution in Helmholtz form.
    a0_mix = a0 + x1 * np.log(x1) + x2 * np.log(x2)
    return float(a0_mix + ar), float(a0_tau + ar_tau), float(ar_del)


def mix_state(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> MixState:
    """
    Compute pressure/enthalpy/Gibbs state at fixed composition.

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]

    Outputs
    -------
    MixState
      p [Pa], h [J/mol], g [J/mol], rho_mol [mol/m^3], rho_mass [kg/m^3], x1 [mol/mol]

    Assumptions
    -----------
    - Composition is fixed for this phase-state evaluation.
    - EOS identity forms for p,h,g from Helmholtz are valid in the evaluated region.

    Failure modes
    -------------
    - Propagates exceptions from reducing functions or alpha evaluations.
    - May return non-physical results if called in unstable/two-phase states.

    References
    ----------
    - Helmholtz identities: Z = 1 + delta*alphar_delta and h/(RT), g/(RT) forms.

    Notes on numerical stability
    ----------------------------
    - Accuracy depends on reduced derivative quality from delegated functions.
    """
    x1 = _clip_x(x1)
    x2 = 1.0 - x1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = x1 * mw1 + x2 * mw2

    alpha, alpha_tau, ar_del = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1)
    # Recover tau by invert identity alpha_tau contribution requires reduced tau.
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k
    delta = rho_mol * vred

    z = 1.0 + delta * ar_del
    p_pa = rho_mol * R_u * t_k * z
    h_jmol = R_u * t_k * (1.0 + tau * alpha_tau + delta * ar_del)
    g_jmol = R_u * t_k * (1.0 + alpha + delta * ar_del)
    return MixState(
        p_pa=float(p_pa),
        h_jmol=float(h_jmol),
        g_jmol=float(g_jmol),
        rho_mol=float(rho_mol),
        rho_mass=float(rho_mol * mw_mix),
        x1=float(x1),
    )


def total_helmholtz_a(d1: Dict, d2: Dict, t_k: float, v_m3: float, n1_mol: float, n2_mol: float) -> float:
    """
    Evaluate total Helmholtz energy A(T,V,n1,n2).

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    v_m3 : float [m^3]
    n1_mol : float [mol]
    n2_mol : float [mol]

    Outputs
    -------
    A : float [J]

    Assumptions
    -----------
    - Homogeneous phase representation at the supplied state.

    Failure modes
    -------------
    - Division by zero if v_m3 <= 0.
    - Propagates EOS evaluation exceptions.

    References
    ----------
    - A = n*R_u*T*alpha_mix.

    Notes on numerical stability
    ----------------------------
    - Small composition floors avoid divide-by-zero in x1 = n1/(n1+n2).
    """
    n1 = max(n1_mol, EPS_X)
    n2 = max(n2_mol, EPS_X)
    n = n1 + n2
    x1 = _clip_x(n1 / n)
    rho_mol = n / v_m3
    alpha, _, _ = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1)
    return float(n * R_u * t_k * alpha)


def _bell2023_reducing_derivs_binary_local(
    x1: float, tc1: float, tc2: float, vc1: float, vc2: float
) -> Tuple[float, float]:
    """
    Purpose
    -------
    Compute analytic composition derivatives of binary reducing functions.

    Inputs
    ------
    x1 : float [mol/mol]
    tc1, tc2 : float [K]
    vc1, vc2 : float [m^3/mol]

    Outputs
    -------
    dTred_dx1 : float [K]
    dvred_dx1 : float [m^3/mol]

    Assumptions
    -----------
    - Binary system with x2 = 1 - x1.
    - Bell parameter set corresponds to R1234ze(E)/R227ea.

    Failure modes
    -------------
    - Numeric issues if inputs are non-physical.

    References
    ----------
    - Bell (2023) corresponding-states reducing function derivatives.

    Notes on numerical stability
    ----------------------------
    - Closed-form derivatives avoid finite-difference noise in composition terms.
    """
    x2 = 1.0 - x1
    p = BELL_2023_R1234ZE_R227EA

    dt = (p.beta_T ** 2) * x1 + x2
    dv = (p.beta_v ** 2) * x1 + x2
    theta_t = 1.0 / dt
    theta_v = 1.0 / dv
    dtheta_t = -((p.beta_T ** 2) - 1.0) / (dt * dt)
    dtheta_v = -((p.beta_v ** 2) - 1.0) / (dv * dv)

    tc12 = p.beta_T * p.gamma_T * np.sqrt(tc1 * tc2)
    vc12 = p.beta_v * p.gamma_v * ((vc1 ** (1.0 / 3.0) + vc2 ** (1.0 / 3.0)) ** 3) / 8.0
    dxx_theta_t = (x2 - x1) * theta_t + x1 * x2 * dtheta_t
    dxx_theta_v = (x2 - x1) * theta_v + x1 * x2 * dtheta_v
    dtred_dx1 = 2.0 * x1 * tc1 - 2.0 * x2 * tc2 + 2.0 * tc12 * dxx_theta_t
    dvred_dx1 = 2.0 * x1 * vc1 - 2.0 * x2 * vc2 + 2.0 * vc12 * dxx_theta_v
    return float(dtred_dx1), float(dvred_dx1)


def chemical_potentials_analytic(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> Tuple[float, float]:
    """
    Compute component chemical potentials using analytic fugacity expressions.

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]

    Outputs
    -------
    mu1, mu2 : float [J/mol]

    Assumptions
    -----------
    - Fugacity relation f_i = x_i*rho*R*T*exp(d(n*alpha^r)/dn_i) is valid.
    - Composition derivatives follow fixed T,rho formulation.

    Failure modes
    -------------
    - Propagates EOS/reducing exceptions.
    - May overflow if exponent arguments are extreme.

    References
    ----------
    - Lemmon and Tillner-Roth mixture fugacity/chemical potential framework.
    - Bell (2023) reducing and departure derivatives.

    Notes on numerical stability
    ----------------------------
    - Uses lower bounds on rho and fugacity log arguments.
    """
    x1 = _clip_x(x1)
    x2 = 1.0 - x1
    rho_mol = max(float(rho_mol), RHO_MIN_MOLM3)

    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k
    delta = rho_mol * vred
    rho_red_mol = 1.0 / vred

    _, _, ar_mix, _, ar_del_mix = mixture_alpha0_alphar_derivs(
        d1=d1,
        d2=d2,
        x1=x1,
        x2=x2,
        tau=tau,
        delta=delta,
        Tred=tred,
        rho_red_mol=rho_red_mol,
        pair_key="r1234ze|r227ea",
    )

    c1 = tc1 / tred
    c2 = tc2 / tred
    k1 = rho_red_mol / rhoc1_mol
    k2 = rho_red_mol / rhoc2_mol
    tau1 = c1 * tau
    tau2 = c2 * tau
    delta1 = k1 * delta
    delta2 = k2 * delta
    ar1, _, _ = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)
    ar2, _, _ = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)

    dtred_dx1, dvred_dx1 = _bell2023_reducing_derivs_binary_local(x1, tc1, tc2, vc1, vc2)
    dtau_dx1 = dtred_dx1 / t_k
    ddelta_dx1 = rho_mol * dvred_dx1
    dep_base = bell2023_departure_base(tau, delta)
    _, dep_tau, dep_del = bell2023_departure_alphar(x1, x2, tau, delta)
    ddep_dx1 = (x2 - x1) * dep_base + dep_tau * dtau_dx1 + dep_del * ddelta_dx1
    dar_dx1 = (ar1 - ar2) + ddep_dx1
    dar_drho = ar_del_mix * vred

    d_na_dn1 = ar_mix + rho_mol * dar_drho + x2 * dar_dx1
    d_na_dn2 = ar_mix + rho_mol * dar_drho - x1 * dar_dx1
    f1_pa = x1 * rho_mol * R_u * t_k * np.exp(d_na_dn1)
    f2_pa = x2 * rho_mol * R_u * t_k * np.exp(d_na_dn2)
    mu1 = R_u * t_k * np.log(max(1e-300, f1_pa))
    mu2 = R_u * t_k * np.log(max(1e-300, f2_pa))
    return float(mu1), float(mu2)


def _sigmoid(z: float) -> float:
    """
    Purpose
    -------
    Map unconstrained scalar to (0,1) interval for phase compositions.

    Inputs
    ------
    z : float [unitless]

    Outputs
    -------
    y : float [unitless]

    Assumptions
    -----------
    - Input can span wide real values.

    Failure modes
    -------------
    - None expected for finite z.

    References
    ----------
    - Logistic transform.

    Notes on numerical stability
    ----------------------------
    - Uses branch form to reduce overflow for large |z|.
    """
    if z >= 0:
        ez = np.exp(-z)
        return float(1.0 / (1.0 + ez))
    ez = np.exp(z)
    return float(ez / (1.0 + ez))


def _safe_residual_vector(vals: np.ndarray) -> np.ndarray:
    """
    Purpose
    -------
    Convert non-finite residuals into large finite penalties for solver safety.

    Inputs
    ------
    vals : np.ndarray [unitless]

    Outputs
    -------
    res : np.ndarray [unitless]

    Assumptions
    -----------
    - Residual vector size is 3 for this VLE system.

    Failure modes
    -------------
    - None; always returns a finite vector.

    References
    ----------
    - Standard robust least-squares penalty strategy.

    Notes on numerical stability
    ----------------------------
    - Prevents NaN/Inf from breaking nonlinear iterations.
    """
    if np.all(np.isfinite(vals)):
        return vals
    return np.array([1.0e6, 1.0e6, 1.0e6], dtype=float)


def _rho_param_to_states(u0: float, u1: float) -> Tuple[float, float]:
    """
    Purpose
    -------
    Map unconstrained variables to physically ordered vapor/liquid densities.

    Inputs
    ------
    u0, u1 : float [unitless]

    Outputs
    -------
    rho_l : float [mol/m^3]
    rho_v : float [mol/m^3]

    Assumptions
    -----------
    - u0 controls rho_v via exp transform.
    - u1 controls rho_l-rho_v via exp transform.

    Failure modes
    -------------
    - None; values are clipped into bounded intervals.

    References
    ----------
    - Bounded reparameterization for nonlinear VLE solves.

    Notes on numerical stability
    ----------------------------
    - Enforces positivity and strict ordering rho_l > rho_v.
    """
    rho_v = float(np.exp(u0))
    drho = float(np.exp(u1))
    rho_v = min(max(rho_v, RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.95)
    drho = min(max(drho, RHO_MIN_MOLM3), RHO_MAX_MOLM3)
    rho_l = min(max(rho_v + drho, rho_v * (1.0 + 1.0e-8)), RHO_MAX_MOLM3)
    rho_v = min(rho_v, rho_l * (1.0 - 1.0e-8))
    return rho_l, rho_v


def solve_bubble_at_t(
    d1: Dict,
    d2: Dict,
    t_k: float,
    z1: float,
    rho_l0: float,
    rho_v0: float,
    y10: float,
) -> Dict:
    """
    Solve bubble state at fixed liquid composition x=z.

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    z1 : float [mol/mol]
        Fixed liquid composition for bubble solve.
    rho_l0, rho_v0 : float [mol/m^3]
        Initial liquid/vapor density guesses.
    y10 : float [mol/mol]
        Initial vapor composition guess.

    Outputs
    -------
    row : dict [mixed units]
      Contains status, pressure, densities, compositions, enthalpies, residuals.

    Assumptions
    -----------
    - Unknown vector is (rho_l, rho_v, y1).
    - Equilibrium equations are P_l=P_v and mu_i^l=mu_i^v for i=1,2.

    Failure modes
    -------------
    - May return DIVERGED if strict residual gates are not met.
    - Can hit penalized residual mode for non-finite trial states.

    References
    ----------
    - VLE equilibrium conditions from Helmholtz-fugacity formulation.

    Notes on numerical stability
    ----------------------------
    - Uses bounded least-squares and ordered density mapping.
    """
    z1 = _clip_x(z1)

    def res(u):
        rho_l, rho_v = _rho_param_to_states(float(u[0]), float(u[1]))
        y1 = _clip_x(_sigmoid(float(u[2])))
        try:
            st_l = mix_state(d1, d2, t_k, rho_l, z1)
            st_v = mix_state(d1, d2, t_k, rho_v, y1)
            mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, z1)
            mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, y1)
            r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
            r2 = (mu1_l - mu1_v) / (R_u * t_k)
            r3 = (mu2_l - mu2_v) / (R_u * t_k)
            return _safe_residual_vector(np.array([r1, r2, r3], dtype=float))
        except Exception:
            return np.array([1.0e6, 1.0e6, 1.0e6], dtype=float)

    rho_v_seed = min(max(float(rho_v0), RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.9)
    rho_l_seed = min(max(float(rho_l0), rho_v_seed * (1.0 + 1.0e-6)), RHO_MAX_MOLM3)
    dr_seed = max(rho_l_seed - rho_v_seed, 1.0e-6)
    y10 = _clip_x(y10)
    u0 = np.array([np.log(rho_v_seed), np.log(dr_seed), np.log(y10 / (1.0 - y10))], dtype=float)
    lb = np.array([np.log(RHO_MIN_MOLM3), np.log(1.0e-9), -30.0], dtype=float)
    ub = np.array([np.log(RHO_MAX_MOLM3), np.log(RHO_MAX_MOLM3), 30.0], dtype=float)
    sol = least_squares(
        res,
        u0,
        bounds=(lb, ub),
        method=LSQ_METHOD,
        ftol=1.0e-12,
        xtol=1.0e-12,
        gtol=1.0e-12,
        max_nfev=LSQ_MAX_NFEV,
        x_scale=LSQ_X_SCALE,
        diff_step=LSQ_DIFF_STEP,
    )

    rho_l, rho_v = _rho_param_to_states(float(sol.x[0]), float(sol.x[1]))
    y1 = _clip_x(_sigmoid(float(sol.x[2])))
    st_l = mix_state(d1, d2, t_k, rho_l, z1)
    st_v = mix_state(d1, d2, t_k, rho_v, y1)
    mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, z1)
    mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, y1)
    r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
    r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)
    ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * (1.0 + 1e-8))
    return {
        "status": "CONVERGED" if ok else "DIVERGED",
        "T_K": float(t_k),
        "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
        "rho_l_molm3": float(rho_l),
        "rho_v_molm3": float(rho_v),
        "x1_liq": float(z1),
        "y1_vap": float(y1),
        "h_l_Jmol": float(st_l.h_jmol),
        "h_v_Jmol": float(st_v.h_jmol),
        "r_P": float(r_p),
        "r_mu": float(r_mu),
        "iterations": int(sol.nfev),
        "notes": str(sol.message),
    }


def solve_dew_at_t(
    d1: Dict,
    d2: Dict,
    t_k: float,
    z1: float,
    rho_l0: float,
    rho_v0: float,
    x10: float,
) -> Dict:
    """
    Solve dew state at fixed vapor composition y=z.

    Inputs
    ------
    d1,d2 : dict [unitless]
    t_k : float [K]
    z1 : float [mol/mol]
        Fixed vapor composition for dew solve.
    rho_l0, rho_v0 : float [mol/m^3]
        Initial liquid/vapor density guesses.
    x10 : float [mol/mol]
        Initial liquid composition guess.

    Outputs
    -------
    row : dict [mixed units]
      Contains status, pressure, densities, compositions, enthalpies, residuals.

    Assumptions
    -----------
    - Unknown vector is (rho_l, rho_v, x1).
    - Equilibrium equations are P_l=P_v and mu_i^l=mu_i^v for i=1,2.

    Failure modes
    -------------
    - May return DIVERGED if strict residual gates are not met.
    - Can hit penalized residual mode for non-finite trial states.

    References
    ----------
    - VLE equilibrium conditions from Helmholtz-fugacity formulation.

    Notes on numerical stability
    ----------------------------
    - Uses bounded least-squares and ordered density mapping.
    """
    z1 = _clip_x(z1)

    def res(u):
        rho_l, rho_v = _rho_param_to_states(float(u[0]), float(u[1]))
        x1 = _clip_x(_sigmoid(float(u[2])))
        try:
            st_l = mix_state(d1, d2, t_k, rho_l, x1)
            st_v = mix_state(d1, d2, t_k, rho_v, z1)
            mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, x1)
            mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, z1)
            r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
            r2 = (mu1_l - mu1_v) / (R_u * t_k)
            r3 = (mu2_l - mu2_v) / (R_u * t_k)
            return _safe_residual_vector(np.array([r1, r2, r3], dtype=float))
        except Exception:
            return np.array([1.0e6, 1.0e6, 1.0e6], dtype=float)

    rho_v_seed = min(max(float(rho_v0), RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.9)
    rho_l_seed = min(max(float(rho_l0), rho_v_seed * (1.0 + 1.0e-6)), RHO_MAX_MOLM3)
    dr_seed = max(rho_l_seed - rho_v_seed, 1.0e-6)
    x10 = _clip_x(x10)
    u0 = np.array([np.log(rho_v_seed), np.log(dr_seed), np.log(x10 / (1.0 - x10))], dtype=float)
    lb = np.array([np.log(RHO_MIN_MOLM3), np.log(1.0e-9), -30.0], dtype=float)
    ub = np.array([np.log(RHO_MAX_MOLM3), np.log(RHO_MAX_MOLM3), 30.0], dtype=float)
    sol = least_squares(
        res,
        u0,
        bounds=(lb, ub),
        method=LSQ_METHOD,
        ftol=1.0e-12,
        xtol=1.0e-12,
        gtol=1.0e-12,
        max_nfev=LSQ_MAX_NFEV,
        x_scale=LSQ_X_SCALE,
        diff_step=LSQ_DIFF_STEP,
    )

    rho_l, rho_v = _rho_param_to_states(float(sol.x[0]), float(sol.x[1]))
    x1 = _clip_x(_sigmoid(float(sol.x[2])))
    st_l = mix_state(d1, d2, t_k, rho_l, x1)
    st_v = mix_state(d1, d2, t_k, rho_v, z1)
    mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, x1)
    mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, z1)
    r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
    r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)
    ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * (1.0 + 1e-8))
    return {
        "status": "CONVERGED" if ok else "DIVERGED",
        "T_K": float(t_k),
        "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
        "rho_l_molm3": float(rho_l),
        "rho_v_molm3": float(rho_v),
        "x1_liq": float(x1),
        "y1_vap": float(z1),
        "h_l_Jmol": float(st_l.h_jmol),
        "h_v_Jmol": float(st_v.h_jmol),
        "r_P": float(r_p),
        "r_mu": float(r_mu),
        "iterations": int(sol.nfev),
        "notes": str(sol.message),
    }


def run_true_vle_envelope(
    fluid1: str,
    fluid2: str,
    w1: float,
    t_vals: np.ndarray,
) -> Tuple[List[Dict], List[Dict], float]:
    """
    Purpose
    -------
    Run bubble and dew VLE solves across a temperature grid with continuation.

    Inputs
    ------
    fluid1, fluid2 : str [unitless]
    w1 : float [kg/kg]
        Component-1 mass fraction.
    t_vals : np.ndarray [K]

    Outputs
    -------
    bubble_rows : list[dict] [mixed units]
    dew_rows : list[dict] [mixed units]
    z1 : float [mol/mol]
        Overall mole fraction converted from input mass fraction.

    Assumptions
    -----------
    - Composition is fixed overall for the envelope run.
    - Continuation updates guesses from previously converged states.

    Failure modes
    -------------
    - Individual temperature points may fail and be marked DIVERGED.

    References
    ----------
    - Predictor-corrector continuation strategy for nonlinear phase-equilibrium traces.

    Notes on numerical stability
    ----------------------------
    - Continuation greatly improves robustness versus independent per-point solves.
    """
    d1 = load_idaes_helmholtz_json(fluid1)
    d2 = load_idaes_helmholtz_json(fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = w1_to_x1(w1, mw1, mw2)

    rhoc1 = float(d1["basic"]["rhoc"]) / mw1
    rhoc2 = float(d2["basic"]["rhoc"]) / mw2
    rho_l_guess = 0.8 * (z1 * rhoc1 + (1.0 - z1) * rhoc2)
    rho_v_guess = 0.01 * rho_l_guess
    y1_guess = z1
    x1_guess = z1

    bubble_rows: List[Dict] = []
    dew_rows: List[Dict] = []
    for t_k in np.asarray(t_vals, dtype=float):
        b = solve_bubble_at_t(d1, d2, float(t_k), z1, rho_l_guess, rho_v_guess, y1_guess)
        bubble_rows.append(b)
        if b["status"] == "CONVERGED":
            rho_l_guess = b["rho_l_molm3"]
            rho_v_guess = b["rho_v_molm3"]
            y1_guess = b["y1_vap"]

        d = solve_dew_at_t(d1, d2, float(t_k), z1, rho_l_guess, rho_v_guess, x1_guess)
        dew_rows.append(d)
        if d["status"] == "CONVERGED":
            rho_l_guess = d["rho_l_molm3"]
            rho_v_guess = d["rho_v_molm3"]
            x1_guess = d["x1_liq"]
    return bubble_rows, dew_rows, z1


def save_csv(rows: List[Dict], path: str | Path) -> None:
    """
    Purpose
    -------
    Persist branch solve rows to CSV.

    Inputs
    ------
    rows : list[dict] [mixed units]
    path : str | Path [filesystem path]

    Outputs
    -------
    None [unitless]

    Assumptions
    -----------
    - All rows share the same schema/keys.

    Failure modes
    -------------
    - File I/O exceptions propagate.

    References
    ----------
    - Standard CSV serialization.

    Notes on numerical stability
    ----------------------------
    - Not applicable.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def plot_envelope(
    bubble_rows: List[Dict],
    dew_rows: List[Dict],
    mw_mix: float,
    out_fig: str | Path,
) -> None:
    """
    Purpose
    -------
    Plot p-h envelope using converged bubble/dew branch points.

    Inputs
    ------
    bubble_rows, dew_rows : list[dict] [mixed units]
    mw_mix : float [kg/mol]
    out_fig : str | Path [filesystem path]

    Outputs
    -------
    None [unitless]
      Writes figure to disk.

    Assumptions
    -----------
    - Enthalpy conversion uses h[J/mol] -> h[kJ/kg] via mw_mix.

    Failure modes
    -------------
    - File I/O/plotting backend exceptions propagate.

    References
    ----------
    - p-h diagram convention for refrigeration cycle review.

    Notes on numerical stability
    ----------------------------
    - Filters to CONVERGED points to avoid plotting numerically invalid states.
    """
    b = [r for r in bubble_rows if r["status"] == "CONVERGED"]
    d = [r for r in dew_rows if r["status"] == "CONVERGED"]

    hb = np.array([(r["h_l_Jmol"] / mw_mix) * 1e-3 for r in b], dtype=float)
    pb = np.array([r["P_Pa"] * 1e-5 for r in b], dtype=float)
    hd = np.array([(r["h_v_Jmol"] / mw_mix) * 1e-3 for r in d], dtype=float)
    pd = np.array([r["P_Pa"] * 1e-5 for r in d], dtype=float)

    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if len(b):
        ax.plot(hb, pb, lw=2.0, label="Bubble line (liq)")
    if len(d):
        ax.plot(hd, pd, lw=2.0, label="Dew line (vap)")
    if len(b) and len(d):
        m = min(len(b), len(d))
        ax.fill_betweenx(pb[:m], hb[:m], hd[:m], alpha=0.15, label="Two-phase region")
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("R515B true VLE envelope (mu-equality solve)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")


def _cli() -> None:
    """
    Purpose
    -------
    Command-line entrypoint for true-VLE envelope generation.

    Inputs
    ------
    CLI args:
    --fluid1, --fluid2 [unitless]
    --w1 [kg/kg]
    --Tmin, --Tmax [K]
    --n [count]
    output file paths [filesystem paths]

    Outputs
    -------
    None [unitless]
      Writes branch CSVs, plot PNG, and metadata JSON.

    Assumptions
    -----------
    - User-provided range and composition are physically meaningful.

    Failure modes
    -------------
    - Raises if solver, plotting, or file operations fail.

    References
    ----------
    - Uses run_true_vle_envelope workflow in this module.

    Notes on numerical stability
    ----------------------------
    - Stability characteristics follow the underlying continuation and nonlinear solve.
    """
    p = argparse.ArgumentParser(description="Compute binary true VLE envelope from chemical-potential equality.")
    p.add_argument("--fluid1", default="r1234ze")
    p.add_argument("--fluid2", default="r227ea")
    p.add_argument("--w1", required=True, type=float, help="mass fraction fluid1 [kg/kg]")
    p.add_argument("--Tmin", required=True, type=float, help="minimum temperature [K]")
    p.add_argument("--Tmax", required=True, type=float, help="maximum temperature [K]")
    p.add_argument("--n", default=80, type=int, help="temperature points")
    p.add_argument("--bubble-csv", default="verification/r515b_true_vle_bubble.csv")
    p.add_argument("--dew-csv", default="verification/r515b_true_vle_dew.csv")
    p.add_argument("--fig", default="verification/r515b_true_vle_envelope.png")
    p.add_argument("--metadata", default="verification/r515b_true_vle_metadata.json")
    args = p.parse_args()

    t_vals = np.linspace(args.Tmin, args.Tmax, args.n)
    bubble_rows, dew_rows, z1 = run_true_vle_envelope(args.fluid1, args.fluid2, args.w1, t_vals)
    save_csv(bubble_rows, args.bubble_csv)
    save_csv(dew_rows, args.dew_csv)

    d1 = load_idaes_helmholtz_json(args.fluid1)
    d2 = load_idaes_helmholtz_json(args.fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = z1 * mw1 + (1.0 - z1) * mw2
    plot_envelope(bubble_rows, dew_rows, mw_mix, args.fig)

    nb = int(sum(r["status"] == "CONVERGED" for r in bubble_rows))
    nd = int(sum(r["status"] == "CONVERGED" for r in dew_rows))
    meta = {
        "fluid1": args.fluid1,
        "fluid2": args.fluid2,
        "w1_kgkg": float(args.w1),
        "z1_molmol": float(z1),
        "Tmin_K": float(args.Tmin),
        "Tmax_K": float(args.Tmax),
        "n_points": int(args.n),
        "bubble_converged": nb,
        "bubble_failed": int(args.n - nb),
        "dew_converged": nd,
        "dew_failed": int(args.n - nd),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "P equality + mu1/mu2 equality, analytic fugacity-based chemical potentials",
    }
    mpath = Path(args.metadata)
    mpath.parent.mkdir(parents=True, exist_ok=True)
    with mpath.open("w") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved bubble CSV: {args.bubble_csv}")
    print(f"Saved dew CSV: {args.dew_csv}")
    print(f"Saved figure: {args.fig}")
    print(f"Saved metadata: {args.metadata}")


if __name__ == "__main__":
    _cli()
