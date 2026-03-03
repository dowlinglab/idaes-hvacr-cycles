#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pressure-validated binary Helmholtz model (from user-provided implementation),
with chart-basis enthalpy support.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from idaes.models.properties.general_helmholtz import get_parameter_path

try:
    from scipy.optimize import brentq as _brentq
except Exception:
    _brentq = None

MW_GMOL_TO_KGMOL = 1e-3
R_KJkgK_TO_JkgK = 1000.0
KPA_TO_PA = 1000.0
BTU_LBM_PER_KJ_KG = 0.429922614


@dataclass(frozen=True)
class InteractionParams:
    betaT12: float
    betarho12: float
    gammaT12: float
    gammanu12: float
    F12: float


@dataclass(frozen=True)
class StatePoint:
    T: float
    rho: float
    w1: float
    w2: float
    alpha: float
    tau: float
    delta: float
    p: float
    h: float


@dataclass(frozen=True)
class ChartReference:
    T_ref_K: float = 273.15
    h_ref_kJkg: float = 200.0
    rho_ref_kgm3: float = 1258.4
    p_ref_kPa: Optional[float] = None


def load_json(refrigerant_name: str, parameter_path: Optional[str | Path] = None):
    if parameter_path is None:
        parameter_path = get_parameter_path()
    parameter_path = Path(parameter_path)
    with open(parameter_path / refrigerant_name, "r") as f:
        return json.load(f)


def _mw_mix_from_mass_fractions(comp1_data, comp2_data, w1, w2):
    MW1 = comp1_data["basic"]["MW"] * MW_GMOL_TO_KGMOL
    MW2 = comp2_data["basic"]["MW"] * MW_GMOL_TO_KGMOL
    MWmix = 1.0 / (w1 / MW1 + w2 / MW2)
    return MW1, MW2, MWmix


def compute_rho_mol(rho_mass, MWmix_kg_per_mol):
    return rho_mass / MWmix_kg_per_mol


def pull_ideal_params_from_json(data):
    eos = data["eos"]
    return {
        "last_term_ideal": eos["last_term_ideal"],
        "n0": eos["n0"],
        "g0": eos["g0"],
        "phi_ideal_type": eos["phi_ideal_type"],
    }


def pull_residual_params_from_json(data):
    eos = data["eos"]
    return {
        "last_term_residual": eos["last_term_residual"],
        "c": eos["c"],
        "d": eos["d"],
        "t": eos["t"],
        "n": eos["n"],
        "a": eos["a"],
        "b": eos["b"],
        "g": eos["g"],
        "e": eos["e"],
        "phi_residual_type": eos["phi_residual_type"],
    }


def compute_alpha_ideal_pure(data, tau, delta):
    alpha0_pure = np.log(delta) + data["n0"]["1"] + data["n0"]["2"] * tau + data["n0"]["3"] * np.log(tau)
    if data["phi_ideal_type"] == 1:
        h = int(data["last_term_ideal"])
        for index in range(4, h + 1):
            ni = float(data["n0"][str(index)])
            gi = float(data["g0"][str(index)])
            alpha0_pure += ni * np.log(1.0 - np.exp(-gi * tau))
    elif data["phi_ideal_type"] == 2:
        h1 = int(data["last_term_ideal"][0])
        h2 = int(data["last_term_ideal"][1])
        for index in range(4, h1 + 1):
            ni = float(data["n0"][str(index)])
            gi = float(data["g0"][str(index)])
            alpha0_pure += ni * (tau**gi)
        for index in range(h1 + 1, h2 + 1):
            ni = float(data["n0"][str(index)])
            gi = float(data["g0"][str(index)])
            alpha0_pure += ni * np.log(1.0 - np.exp(-gi * tau))
    elif data["phi_ideal_type"] == 3:
        h = int(data["last_term_ideal"])
        for index in range(4, h + 1):
            ni = float(data["n0"][str(index)])
            gi = float(data["g0"][str(index)])
            alpha0_pure += ni * (tau**gi)
    else:
        raise ValueError("Unexpected ideal type")
    return alpha0_pure


def compute_alpha_res_pure(data, tau, delta):
    alphares_pure = 0.0
    phi = data["phi_residual_type"]
    hlist = data["last_term_residual"]
    t1 = float(data["t"]["1"])
    h1 = int(hlist[0])
    h2 = int(hlist[1])
    if phi == 1:
        for index in range(1, h1 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            ti = float(data["t"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**ti)
        for index in range(h1 + 1, h2 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            ti = float(data["t"][str(index)])
            ci = float(data["c"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**ti) * np.exp(-(delta**ci))
    elif phi == 2:
        h3 = int(hlist[2])
        for index in range(1, h1 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**t1)
        for index in range(h1 + 1, h2 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            ci = float(data["c"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**t1) * np.exp(-(delta**ci))
        for index in range(h2 + 1, h3 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            ai = float(data["a"][str(index)])
            bi = float(data["b"][str(index)])
            ei = float(data["e"][str(index)])
            gi = float(data["g"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**t1) * np.exp(-ai * ((delta - ei) ** 2) - bi * ((tau - gi) ** 2))
    elif phi == 3:
        h3 = int(hlist[2])
        for index in range(1, h1 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**t1)
        for index in range(h1 + 1, h2 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            ci = float(data["c"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**t1) * np.exp(-(delta**ci))
        for index in range(h2 + 1, h3 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            ci = float(data["c"][str(index)])
            bi = float(data["b"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**t1) * np.exp(-(delta**ci)) * np.exp(-(tau**bi))
    elif phi == 4:
        h0 = int(hlist[0])
        for index in range(1, h0 + 1):
            ni = float(data["n"][str(index)])
            di = float(data["d"][str(index)])
            alphares_pure += ni * (delta**di) * (tau**t1)
        m = len(hlist) - 1
        for j in range(1, m + 1):
            hjm1 = int(hlist[j - 1])
            hj = int(hlist[j])
            summation = 0.0
            for i in range(hjm1 + 1, hj + 1):
                ni = float(data["n"][str(i)])
                di = float(data["d"][str(i)])
                summation += ni * (delta**di) * (tau**t1)
            alphares_pure += np.exp(-(delta**j)) * summation
    else:
        raise ValueError("Unexpected residual type")
    return alphares_pure


def _pure_reduced_vars(data, T, rho_mol):
    b = data["basic"]
    Tc = b["Tc"]
    rhoc_mass = b["rhoc"]
    MW_kg_per_mol = b["MW"] * MW_GMOL_TO_KGMOL
    rhoc_mol = rhoc_mass / MW_kg_per_mol
    return Tc / T, rho_mol / rhoc_mol


def reducing_functions(data1, data2, w1, w2, betaT12, betarho12, gammaT12, gammanu12):
    b1 = data1["basic"]
    b2 = data2["basic"]
    Tc1, Tc2 = b1["Tc"], b2["Tc"]
    rhoc1, rhoc2 = b1["rhoc"], b2["rhoc"]
    MW1, MW2 = b1["MW"] * MW_GMOL_TO_KGMOL, b2["MW"] * MW_GMOL_TO_KGMOL
    n1 = w1 / MW1
    n2 = w2 / MW2
    x1 = n1 / (n1 + n2)
    x2 = n2 / (n1 + n2)

    rhoc1_mol = rhoc1 / MW1
    rhoc2_mol = rhoc2 / MW2
    nu1 = 1.0 / rhoc1_mol
    nu2 = 1.0 / rhoc2_mol

    T12 = betaT12 * gammaT12 * np.sqrt(Tc1 * Tc2)
    Tred = (x1**2) * Tc1 + (x2**2) * Tc2 + 2 * x1 * x2 * ((x1 + x2) / ((betaT12**2) * x1 + x2)) * T12

    nu12 = betarho12 * gammanu12 * ((nu1 ** (1 / 3) + nu2 ** (1 / 3)) ** 3)
    nured = (x1**2) * nu1 + (x2**2) * nu2 + 2 * x1 * x2 * ((x1 + x2) / ((betarho12**2) * x1 + x2)) * nu12
    return x1, x2, Tred, nured


def compute_alpha_mix(data1, data2, T, rho_mol, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12):
    ideal_data1 = pull_ideal_params_from_json(data1)
    ideal_data2 = pull_ideal_params_from_json(data2)
    residual_data1 = pull_residual_params_from_json(data1)
    residual_data2 = pull_residual_params_from_json(data2)

    x1, x2, Tred, nured = reducing_functions(data1, data2, w1, w2, betaT12, betarho12, gammaT12, gammanu12)
    tau = Tred / T
    delta = rho_mol * nured

    tau1, delta1 = _pure_reduced_vars(data1, T, rho_mol)
    tau2, delta2 = _pure_reduced_vars(data2, T, rho_mol)
    alpha0_pure_1 = compute_alpha_ideal_pure(ideal_data1, tau1, delta1)
    alpha0_pure_2 = compute_alpha_ideal_pure(ideal_data2, tau2, delta2)
    alpha0_mix = x1 * alpha0_pure_1 + x2 * alpha0_pure_2 + x1 * np.log(x1) + x2 * np.log(x2)

    alphares_pure_1 = compute_alpha_res_pure(residual_data1, tau1, delta1)
    alphares_pure_2 = compute_alpha_res_pure(residual_data2, tau2, delta2)
    alphares_mix = x1 * alphares_pure_1 + x2 * alphares_pure_2

    alphar12 = (
        -0.057178 * (tau**1.290298) * delta * np.exp(-delta)
        + 0.031318 * (tau**0.038796) * (delta**2) * np.exp(-delta)
        - 0.027496 * (tau**2.640532) * (delta**3) * np.exp(-delta)
    )
    alpha_dep = x1 * x2 * F12 * alphar12
    return alpha0_mix + alphares_mix + alpha_dep, tau, delta


def alpha_at_tau_delta(data1, data2, tau, delta, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12):
    x1, x2, Tred, nured = reducing_functions(data1, data2, w1, w2, betaT12, betarho12, gammaT12, gammanu12)
    T = Tred / tau
    rho_mol = delta / nured
    alpha, _, _ = compute_alpha_mix(data1, data2, T, rho_mol, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12)
    return alpha


def alpha_and_partials_tau_delta(data1, data2, T, rho_mol, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12, dtau=1e-6, ddel=1e-6):
    _, _, Tred, nured = reducing_functions(data1, data2, w1, w2, betaT12, betarho12, gammaT12, gammanu12)
    tau = Tred / T
    delta = rho_mol * nured

    a0 = alpha_at_tau_delta(data1, data2, tau, delta, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12)
    ap = alpha_at_tau_delta(data1, data2, tau + dtau, delta, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12)
    am = alpha_at_tau_delta(data1, data2, tau - dtau, delta, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12)
    a_tau = (ap - am) / (2 * dtau)

    dp = alpha_at_tau_delta(data1, data2, tau, delta + ddel, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12)
    dm = alpha_at_tau_delta(data1, data2, tau, delta - ddel, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12)
    a_del = (dp - dm) / (2 * ddel)
    return a0, a_tau, a_del, tau, delta


def p_h_from_T_rho(data1, data2, T, rho_mass, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12, Ru=8.314462618):
    MW1, MW2, MWmix = _mw_mix_from_mass_fractions(data1, data2, w1, w2)
    rho_mol = compute_rho_mol(rho_mass, MWmix)
    _a, a_tau, a_del, tau, delta = alpha_and_partials_tau_delta(
        data1, data2, T, rho_mol, w1, w2, betaT12, F12, betarho12, gammaT12, gammanu12
    )
    a_del_res = a_del - 1.0 / delta
    p_Pa = rho_mol * Ru * T * (1.0 + delta * a_del_res)
    h_molar = Ru * T * (1.0 + tau * a_tau + delta * a_del_res)
    return p_Pa / KPA_TO_PA, h_molar


def compute_props(
    comp1_json: str,
    comp2_json: str,
    T: float,
    rho: float,
    w1: float,
    w2: float,
    interaction: InteractionParams,
    parameter_path: Optional[str | Path] = None,
) -> StatePoint:
    comp1_data = load_json(comp1_json, parameter_path)
    comp2_data = load_json(comp2_json, parameter_path)
    MW1, MW2, MWmix = _mw_mix_from_mass_fractions(comp1_data, comp2_data, w1, w2)
    rho_mol = compute_rho_mol(rho, MWmix)
    alpha, tau, delta = compute_alpha_mix(
        comp1_data, comp2_data, T, rho_mol, w1, w2,
        interaction.betaT12, interaction.F12, interaction.betarho12, interaction.gammaT12, interaction.gammanu12
    )
    p, h_molar = p_h_from_T_rho(
        comp1_data, comp2_data, T, rho, w1, w2,
        interaction.betaT12, interaction.F12, interaction.betarho12, interaction.gammaT12, interaction.gammanu12
    )
    h_kJkg = (h_molar / MWmix) / R_KJkgK_TO_JkgK
    return StatePoint(T, rho, w1, w2, alpha, tau, delta, p, h_kJkg)


def default_interaction() -> InteractionParams:
    return InteractionParams(betaT12=1.001247, betarho12=0.99290, gammaT12=0.989180, gammanu12=1.001581, F12=1.0)


def solve_rho_mass_for_P(T_K: float, P_target_kPa: float, comp1_json="r1234ze.json", comp2_json="r227ea.json", w1=0.911, w2=0.089, interaction: Optional[InteractionParams] = None, phase_hint="vapor", rho_min=1e-3, rho_max=2e3, ngrid=600):
    if interaction is None:
        interaction = default_interaction()
    rhos = np.logspace(np.log10(rho_min), np.log10(rho_max), ngrid)
    vals = []
    for rr in rhos:
        st = compute_props(comp1_json, comp2_json, T_K, rr, w1, w2, interaction)
        vals.append(st.p - P_target_kPa)
    vals = np.array(vals)

    brackets = []
    for i in range(len(rhos) - 1):
        f1, f2 = vals[i], vals[i + 1]
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
            f = lambda rr: compute_props(comp1_json, comp2_json, T_K, rr, w1, w2, interaction).p - P_target_kPa
            root = _brentq(f, a, b, maxiter=250)
        else:
            fa = compute_props(comp1_json, comp2_json, T_K, a, w1, w2, interaction).p - P_target_kPa
            fb = compute_props(comp1_json, comp2_json, T_K, b, w1, w2, interaction).p - P_target_kPa
            for _ in range(250):
                m = 0.5 * (a + b)
                fm = compute_props(comp1_json, comp2_json, T_K, m, w1, w2, interaction).p - P_target_kPa
                if fa * fm <= 0.0:
                    b, fb = m, fm
                else:
                    a, fa = m, fm
            root = 0.5 * (a + b)
        roots.append(float(root))

    roots = sorted(set(roots))
    return roots[-1] if phase_hint == "liquid" else roots[0]


def chart_offset(
    ref: ChartReference,
    comp1_json="r1234ze.json",
    comp2_json="r227ea.json",
    w1=0.911,
    w2=0.089,
    interaction: Optional[InteractionParams] = None,
) -> Tuple[float, float]:
    if interaction is None:
        interaction = default_interaction()
    rho_ref = ref.rho_ref_kgm3
    if ref.p_ref_kPa is not None:
        rho_ref_solved = solve_rho_mass_for_P(
            T_K=ref.T_ref_K,
            P_target_kPa=ref.p_ref_kPa,
            comp1_json=comp1_json,
            comp2_json=comp2_json,
            w1=w1,
            w2=w2,
            interaction=interaction,
            phase_hint="liquid",
        )
        if rho_ref_solved is None:
            raise RuntimeError("Could not solve rho_sat_liq from (T_ref, P_ref).")
        rho_ref = rho_ref_solved
    st_ref = compute_props(comp1_json, comp2_json, ref.T_ref_K, rho_ref, w1, w2, interaction)
    return ref.h_ref_kJkg - st_ref.h, rho_ref


def compute_props_with_chart_h(
    T: float,
    rho: float,
    ref: ChartReference,
    comp1_json="r1234ze.json",
    comp2_json="r227ea.json",
    w1=0.911,
    w2=0.089,
    interaction: Optional[InteractionParams] = None,
):
    if interaction is None:
        interaction = default_interaction()
    st = compute_props(comp1_json, comp2_json, T, rho, w1, w2, interaction)
    hoff, rho_ref_used = chart_offset(ref, comp1_json, comp2_json, w1, w2, interaction)
    h_chart = st.h + hoff
    return st, h_chart, hoff, rho_ref_used

