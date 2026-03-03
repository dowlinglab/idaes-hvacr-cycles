#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: Solve pure-fluid saturation states with Helmholtz EOS and
produce P-h dome data using IDAES JSON parameters and derivative helpers.
Dependencies: numpy, matplotlib, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-02
# Assumptions:
# - Pure-fluid saturation only (no mixture flash in this module).
# - Molar basis is used internally: rho [mol/m^3], p [Pa], h/g [J/mol].
# - Newton solve enforces P_l=P_v and g_l=g_v with damping and continuation.
# - Above-critical temperatures are not solved.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

from linear_model_codex import (
    R_u,
    alpha0_idaes_with_derivs,
    load_idaes_helmholtz_json,
    mw_from_json,
)


@dataclass
class PureState:
    rho_mol: float
    p_pa: float
    h_jmol: float
    g_jmol: float
    dpdrho: float
    dgdrho: float


@dataclass
class SaturationResult:
    T_K: float
    status: str
    rho_l_molm3: float
    rho_v_molm3: float
    p_Pa: float
    h_l_Jmol: float
    h_v_Jmol: float
    r_P: float
    r_mu: float
    iterations: int
    fallback_used: bool
    notes: str
    jacobian_cond: float
    step_norm: float


def _sat_density_guess_from_aux(data: Dict, T: float) -> Tuple[float, float]:
    """
    Estimate liquid/vapor saturation densities from IDAES auxiliary fits.

    Returns
    -------
    rho_l_guess, rho_v_guess : float [mol/m^3]
    """
    Tc = float(data["basic"]["Tc"])
    MW = mw_from_json(data)
    rhoc_mol = float(data["basic"]["rhoc"]) / MW
    theta = max(1e-12, 1.0 - T / Tc)

    aux = data.get("aux", {})
    dl = aux.get("delta_l_sat_approx")
    dv = aux.get("delta_v_sat_approx")

    if not dl or not dv:
        return 2.0 * rhoc_mol, 0.02 * rhoc_mol

    def _series(spec: Dict) -> float:
        n = spec.get("n", {})
        t = spec.get("t", {})
        idx = sorted(n.keys(), key=lambda k: int(k))
        return sum(float(n[k]) * (theta ** float(t[k])) for k in idx)

    sl = _series(dl)
    sv = _series(dv)

    # IDAES ancillary convention: type 1 -> delta = 1 + series, type 2 -> ln(delta) = series
    delta_l = 1.0 + sl if int(dl.get("type", 1)) == 1 else float(np.exp(sl))
    delta_v = float(np.exp(sv)) if int(dv.get("type", 2)) == 2 else max(1e-12, 1.0 + sv)

    rho_l = max(1e-9, delta_l * rhoc_mol)
    rho_v = max(1e-12, delta_v * rhoc_mol)
    if rho_v >= rho_l:
        rho_v = 0.1 * rho_l
    return float(rho_l), float(rho_v)


def _exp_d1_d2(delta: float, cval: float) -> Tuple[float, float, float]:
    """Return exp(-(delta**c)), first and second derivatives wrt delta."""
    extra = np.exp(-(delta**cval))
    d1 = -cval * (delta ** (cval - 1.0))
    d2 = -cval * (cval - 1.0) * (delta ** (cval - 2.0))
    extra_del = extra * d1
    extra_deldel = extra * (d2 + d1 * d1)
    return float(extra), float(extra_del), float(extra_deldel)


def alphar_idaes_with_second_derivs(eos: Dict, tau: float, delta: float) -> Tuple[float, float, float, float]:
    """
    Evaluate residual alpha and derivatives wrt tau, delta, and delta-delta.

    Inputs
    ------
    eos : dict [unitless]
    tau : float [unitless]
    delta : float [unitless]

    Outputs
    -------
    alphar, alphar_tau, alphar_del, alphar_deldel : float [unitless]

    Notes
    -----
    Uses analytic term-by-term differentiation for supported IDAES residual
    families (phi_residual_type 1..4).
    """
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

    t1 = float(t["1"])

    alphar = 0.0
    alphar_tau = 0.0
    alphar_del = 0.0
    alphar_deldel = 0.0

    def add_term(poly_coeff: float, di: float, ti: float, extra: float, extra_tau: float, extra_del: float, extra_deldel: float):
        nonlocal alphar, alphar_tau, alphar_del, alphar_deldel
        pref = poly_coeff * (delta ** di) * (tau ** ti)
        base = pref * extra
        alphar += base
        alphar_tau += pref * extra_tau + base * (ti / tau)
        alphar_del += pref * (extra_del + extra * (di / delta))
        alphar_deldel += pref * (
            extra_deldel
            + 2.0 * (di / delta) * extra_del
            + extra * (di * (di - 1.0) / (delta * delta))
        )

    if phi == 1:
        h1 = int(hlist[0]); h2 = int(hlist[1])
        for i in range(1, h1 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            add_term(ni, di, ti, 1.0, 0.0, 0.0, 0.0)
        for i in range(h1 + 1, h2 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ti = float(t[str(i)])
            ci = float(c[str(i)])
            extra, extra_del, extra_deldel = _exp_d1_d2(delta, ci)
            add_term(ni, di, ti, extra, 0.0, extra_del, extra_deldel)

    elif phi == 2:
        h1 = int(hlist[0]); h2 = int(hlist[1]); h3 = int(hlist[2])
        for i in range(1, h1 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)])
            add_term(ni, di, t1, 1.0, 0.0, 0.0, 0.0)
        for i in range(h1 + 1, h2 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ci = float(c[str(i)])
            extra, extra_del, extra_deldel = _exp_d1_d2(delta, ci)
            add_term(ni, di, t1, extra, 0.0, extra_del, extra_deldel)
        for i in range(h2 + 1, h3 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)])
            ai = float(a[str(i)]); bi = float(b[str(i)])
            ei = float(e[str(i)]); gi = float(g[str(i)])
            extra = np.exp(-ai * ((delta - ei) ** 2) - bi * ((tau - gi) ** 2))
            d1 = -(2.0 * ai * (delta - ei))
            d2 = -2.0 * ai
            extra_del = extra * d1
            extra_deldel = extra * (d2 + d1 * d1)
            extra_tau = extra * (-(2.0 * bi * (tau - gi)))
            add_term(ni, di, t1, float(extra), float(extra_tau), float(extra_del), float(extra_deldel))

    elif phi == 3:
        h1 = int(hlist[0]); h2 = int(hlist[1]); h3 = int(hlist[2])
        for i in range(1, h1 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)])
            add_term(ni, di, t1, 1.0, 0.0, 0.0, 0.0)
        for i in range(h1 + 1, h2 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)]); ci = float(c[str(i)])
            extra, extra_del, extra_deldel = _exp_d1_d2(delta, ci)
            add_term(ni, di, t1, extra, 0.0, extra_del, extra_deldel)
        for i in range(h2 + 1, h3 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)])
            ci = float(c[str(i)]); bi = float(b[str(i)])
            expt = np.exp(-(tau ** bi))
            expd, expd_del, expd_deldel = _exp_d1_d2(delta, ci)
            extra = expd * expt
            extra_del = expd_del * expt
            extra_deldel = expd_deldel * expt
            extra_tau = extra * (-(bi * (tau ** (bi - 1.0))))
            add_term(ni, di, t1, float(extra), float(extra_tau), float(extra_del), float(extra_deldel))

    elif phi == 4:
        h0 = int(hlist[0])
        for i in range(1, h0 + 1):
            ni = float(n[str(i)]); di = float(d[str(i)])
            add_term(ni, di, t1, 1.0, 0.0, 0.0, 0.0)

        m = len(hlist) - 1
        for j in range(1, m + 1):
            hjm1 = int(hlist[j - 1])
            hj = int(hlist[j])
            extra, extra_del, extra_deldel = _exp_d1_d2(delta, float(j))
            for i in range(hjm1 + 1, hj + 1):
                ni = float(n[str(i)]); di = float(d[str(i)])
                add_term(ni, di, t1, extra, 0.0, extra_del, extra_deldel)

    else:
        raise ValueError(f"Unexpected phi_residual_type={phi}")

    return float(alphar), float(alphar_tau), float(alphar_del), float(alphar_deldel)


def pure_state_properties(data: Dict, T: float, rho_mol: float) -> PureState:
    """
    Compute pure-fluid P, h, g and rho-derivatives for saturation solve.

    Inputs
    ------
    data : dict [unitless]
    T : float [K]
    rho_mol : float [mol/m^3]

    Outputs
    -------
    PureState
        p [Pa], h [J/mol], g [J/mol], dp/drho [Pa/(mol/m^3)],
        dg/drho [J/mol per (mol/m^3)].
    """
    Tc = float(data["basic"]["Tc"])
    MW = mw_from_json(data)
    rhoc_mol = float(data["basic"]["rhoc"]) / MW

    tau = Tc / T
    delta = rho_mol / rhoc_mol
    if delta <= 0.0:
        raise ValueError("delta must be positive")

    a0, a0_tau = alpha0_idaes_with_derivs(data["eos"], tau, delta)
    ar, ar_tau, ar_del, ar_deldel = alphar_idaes_with_second_derivs(data["eos"], tau, delta)

    alpha = a0 + ar
    Z = 1.0 + delta * ar_del

    p_pa = rho_mol * R_u * T * Z
    h_jmol = R_u * T * (1.0 + tau * (a0_tau + ar_tau) + delta * ar_del)
    g_jmol = R_u * T * (1.0 + alpha + delta * ar_del)

    dpdrho = R_u * T * (1.0 + 2.0 * delta * ar_del + (delta * delta) * ar_deldel)
    dgdrho = (R_u * T / rhoc_mol) * ((1.0 / delta) + 2.0 * ar_del + delta * ar_deldel)

    return PureState(
        rho_mol=float(rho_mol),
        p_pa=float(p_pa),
        h_jmol=float(h_jmol),
        g_jmol=float(g_jmol),
        dpdrho=float(dpdrho),
        dgdrho=float(dgdrho),
    )


def _residual_and_jac(data: Dict, T: float, rho_l: float, rho_v: float):
    st_l = pure_state_properties(data, T, rho_l)
    st_v = pure_state_properties(data, T, rho_v)

    r1 = st_l.p_pa - st_v.p_pa
    r2 = st_l.g_jmol - st_v.g_jmol
    r = np.array([r1, r2], dtype=float)

    J = np.array(
        [[st_l.dpdrho, -st_v.dpdrho], [st_l.dgdrho, -st_v.dgdrho]],
        dtype=float,
    )
    return r, J, st_l, st_v


def _solve_with_scipy(data: Dict, T: float, rho_l0: float, rho_v0: float, tol: float):
    """Fallback nonlinear solve in log-density space using scipy (never auto-accepted)."""
    MW = mw_from_json(data)
    rhoc = float(data["basic"]["rhoc"]) / MW
    min_gap = max(1e-8 * rhoc, 1e-6)

    def fun(y):
        rv = float(np.exp(y[0]))
        rl = rv + min_gap + float(np.exp(y[1]))
        r, _, _, _ = _residual_and_jac(data, T, rl, rv)
        return r

    def jac(y):
        rv = float(np.exp(y[0]))
        eg = float(np.exp(y[1]))
        rl = rv + min_gap + eg
        _, J, _, _ = _residual_and_jac(data, T, rl, rv)
        d_rl = J[:, 0]
        d_rv = J[:, 1]
        Jy = np.zeros_like(J)
        Jy[:, 0] = d_rl * rv + d_rv * rv
        Jy[:, 1] = d_rl * eg
        return Jy

    rv0 = max(min(rho_v0, 0.95 * rho_l0), 1e-12)
    gap0 = max(rho_l0 - rv0 - min_gap, 1e-12)
    y0 = np.array([np.log(rv0), np.log(gap0)], dtype=float)
    sol = root(fun, y0, jac=jac, method="hybr", tol=tol)
    if (not sol.success) or (not np.all(np.isfinite(sol.x))):
        raise RuntimeError(f"scipy root failed: {sol.message}")

    rho_v = float(np.exp(sol.x[0]))
    rho_l = rho_v + min_gap + float(np.exp(sol.x[1]))
    return rho_l, rho_v


def _acceptance_metrics(st_l: PureState, st_v: PureState, T: float) -> Tuple[float, float, bool]:
    """
    Compute strict acceptance metrics.

    r_P = |P_l - P_v| / max(1.0, 0.5*(P_l + P_v)) <= 1e-6
    r_mu = |g_l - g_v| / (R_u*T) <= 1e-6
    rho_l > rho_v*(1+1e-8)
    """
    denom_p = max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
    r_p = abs(st_l.p_pa - st_v.p_pa) / denom_p
    r_mu = abs(st_l.g_jmol - st_v.g_jmol) / (R_u * T)
    rho_ok = st_l.rho_mol > st_v.rho_mol * (1.0 + 1e-8)
    return float(r_p), float(r_mu), bool(rho_ok)


def _newton_attempt(
    data: Dict,
    T: float,
    rho_l0: float,
    rho_v0: float,
    maxiter: int,
    base_lambda: float,
) -> Tuple[float, float, int, bool, str, float, float, List[str]]:
    """
    Run one damped Newton attempt.

    Returns:
    rho_l, rho_v, iterations, accepted_step, notes, jacobian_cond, step_norm, trace
    """
    rho_l = max(float(rho_l0), 1e-9)
    rho_v = max(min(float(rho_v0), 0.95 * rho_l), 1e-12)
    trace: List[str] = []
    jac_cond = np.nan
    step_norm = np.nan

    for it in range(1, maxiter + 1):
        r, J, st_l, st_v = _residual_and_jac(data, T, rho_l, rho_v)
        rnorm = float(np.linalg.norm(r))
        trace.append(f"iter={it} rho_l={rho_l:.8e} rho_v={rho_v:.8e} r1={r[0]:.8e} r2={r[1]:.8e} norm={rnorm:.8e}")

        try:
            jac_cond = float(np.linalg.cond(J))
        except np.linalg.LinAlgError:
            jac_cond = np.inf

        if (not np.isfinite(jac_cond)) or jac_cond > 1e12:
            eps_l = max(1e-9, 1e-6 * rho_l)
            eps_v = max(1e-12, 1e-6 * rho_v)
            r_lp, _, _, _ = _residual_and_jac(data, T, rho_l + eps_l, rho_v)
            r_lm, _, _, _ = _residual_and_jac(data, T, max(1e-12, rho_l - eps_l), rho_v)
            r_vp, _, _, _ = _residual_and_jac(data, T, rho_l, rho_v + eps_v)
            r_vm, _, _, _ = _residual_and_jac(data, T, rho_l, max(1e-12, rho_v - eps_v))
            J = np.column_stack(((r_lp - r_lm) / (2.0 * eps_l), (r_vp - r_vm) / (2.0 * eps_v)))
            trace.append("  finite-difference Jacobian used")

        try:
            dx = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            return rho_l, rho_v, it, False, "Jacobian singular", jac_cond, float("nan"), trace

        dx[0] = np.clip(dx[0], -0.5 * rho_l, 0.5 * rho_l)
        dx[1] = np.clip(dx[1], -0.5 * rho_v, 0.5 * rho_v)
        step_norm = float(np.linalg.norm(dx))

        lam = float(base_lambda)
        accepted = False
        for _ in range(24):
            rho_l_try = rho_l + lam * dx[0]
            rho_v_try = rho_v + lam * dx[1]
            if rho_l_try <= 0.0 or rho_v_try <= 0.0 or rho_v_try >= rho_l_try:
                lam *= 0.5
                continue
            r_try, _, _, _ = _residual_and_jac(data, T, rho_l_try, rho_v_try)
            if np.linalg.norm(r_try) < (1.0 - 1e-4 * lam) * rnorm:
                rho_l = float(rho_l_try)
                rho_v = float(rho_v_try)
                accepted = True
                break
            lam *= 0.5

        if not accepted:
            return rho_l, rho_v, it, False, "line-search rejected step", jac_cond, step_norm, trace

    return rho_l, rho_v, maxiter, True, "maxiter reached", jac_cond, step_norm, trace


def _solve_saturation_with_status(
    data: Dict,
    T: float,
    maxiter: int,
    rho_l_guess: Optional[float],
    rho_v_guess: Optional[float],
    log_handle,
) -> SaturationResult:
    """
    Solve one T point with strict acceptance gates and explicit status.
    """
    Tc = float(data["basic"]["Tc"])
    if T >= Tc:
        return SaturationResult(
            T_K=float(T),
            status="DIVERGED",
            rho_l_molm3=float("nan"),
            rho_v_molm3=float("nan"),
            p_Pa=float("nan"),
            h_l_Jmol=float("nan"),
            h_v_Jmol=float("nan"),
            r_P=float("inf"),
            r_mu=float("inf"),
            iterations=0,
            fallback_used=False,
            notes=f"T >= Tc ({Tc:.6g} K)",
            jacobian_cond=float("nan"),
            step_norm=float("nan"),
        )

    rho_l_aux, rho_v_aux = _sat_density_guess_from_aux(data, T)
    rho_l0 = float(rho_l_guess if rho_l_guess is not None else rho_l_aux)
    rho_v0 = float(rho_v_guess if rho_v_guess is not None else rho_v_aux)
    rho_l0 = max(rho_l0, 1e-9)
    rho_v0 = max(min(rho_v0, 0.95 * rho_l0), 1e-12)

    attempts = [1.0, 0.5, 0.25]
    best = None
    used_fallback_seed = (rho_l_guess is None or rho_v_guess is None)
    failed_reason = "DIVERGED"
    all_traces: List[str] = []

    for idx, lam0 in enumerate(attempts, start=1):
        rho_l, rho_v, iters, _, note, jac_cond, step_norm, trace = _newton_attempt(
            data=data,
            T=T,
            rho_l0=rho_l0,
            rho_v0=rho_v0,
            maxiter=maxiter,
            base_lambda=lam0,
        )
        st_l = pure_state_properties(data, T, rho_l)
        st_v = pure_state_properties(data, T, rho_v)
        r_p, r_mu, rho_ok = _acceptance_metrics(st_l, st_v, T)

        all_traces.extend([f"attempt={idx} base_lambda={lam0:.3f}"] + trace + [
            f"attempt_summary r_P={r_p:.8e} r_mu={r_mu:.8e} rho_ok={rho_ok} jac_cond={jac_cond:.8e} step_norm={step_norm:.8e} note={note}",
            "",
        ])

        current = SaturationResult(
            T_K=float(T),
            status="CONVERGED" if (r_p <= 1e-6 and r_mu <= 1e-6 and rho_ok) else "FAILED_CHECKS",
            rho_l_molm3=float(rho_l),
            rho_v_molm3=float(rho_v),
            p_Pa=float(0.5 * (st_l.p_pa + st_v.p_pa)),
            h_l_Jmol=float(st_l.h_jmol),
            h_v_Jmol=float(st_v.h_jmol),
            r_P=float(r_p),
            r_mu=float(r_mu),
            iterations=int(iters),
            fallback_used=bool(used_fallback_seed),
            notes=f"attempt={idx}; {note}",
            jacobian_cond=float(jac_cond),
            step_norm=float(step_norm),
        )

        if best is None or (current.r_P + current.r_mu) < (best.r_P + best.r_mu):
            best = current
        if current.status == "CONVERGED":
            break

        failed_reason = "FAILED_CHECKS"

    # Optional scipy retry for recovery only; still must pass hard gates.
    if best is None or best.status != "CONVERGED":
        try:
            rl_s, rv_s = _solve_with_scipy(data, T, best.rho_l_molm3 if best else rho_l0, best.rho_v_molm3 if best else rho_v0, tol=1e-12)
            st_l = pure_state_properties(data, T, rl_s)
            st_v = pure_state_properties(data, T, rv_s)
            r_p, r_mu, rho_ok = _acceptance_metrics(st_l, st_v, T)
            all_traces.append(f"scipy_retry r_P={r_p:.8e} r_mu={r_mu:.8e} rho_ok={rho_ok}")
            if r_p <= 1e-6 and r_mu <= 1e-6 and rho_ok:
                best = SaturationResult(
                    T_K=float(T),
                    status="CONVERGED",
                    rho_l_molm3=float(rl_s),
                    rho_v_molm3=float(rv_s),
                    p_Pa=float(0.5 * (st_l.p_pa + st_v.p_pa)),
                    h_l_Jmol=float(st_l.h_jmol),
                    h_v_Jmol=float(st_v.h_jmol),
                    r_P=float(r_p),
                    r_mu=float(r_mu),
                    iterations=int(best.iterations + 1 if best else 1),
                    fallback_used=True,
                    notes="scipy_retry_converged",
                    jacobian_cond=float(best.jacobian_cond if best else np.nan),
                    step_norm=float(best.step_norm if best else np.nan),
                )
            else:
                failed_reason = "FAILED_CHECKS"
        except Exception as err:
            all_traces.append(f"scipy_retry_failed: {err}")
            failed_reason = "DIVERGED"

    if best is None:
        best = SaturationResult(
            T_K=float(T),
            status="DIVERGED",
            rho_l_molm3=float(rho_l_aux),
            rho_v_molm3=float(rho_v_aux),
            p_Pa=float("nan"),
            h_l_Jmol=float("nan"),
            h_v_Jmol=float("nan"),
            r_P=float("inf"),
            r_mu=float("inf"),
            iterations=0,
            fallback_used=True,
            notes="no candidate produced",
            jacobian_cond=float("nan"),
            step_norm=float("nan"),
        )
    elif best.status != "CONVERGED":
        # Keep diagnostic point from aux estimates for plotting failures.
        st_l_aux = pure_state_properties(data, T, rho_l_aux)
        st_v_aux = pure_state_properties(data, T, rho_v_aux)
        r_p_aux, r_mu_aux, _ = _acceptance_metrics(st_l_aux, st_v_aux, T)
        best = SaturationResult(
            T_K=float(T),
            status=failed_reason,
            rho_l_molm3=float(rho_l_aux),
            rho_v_molm3=float(rho_v_aux),
            p_Pa=float(0.5 * (st_l_aux.p_pa + st_v_aux.p_pa)),
            h_l_Jmol=float(st_l_aux.h_jmol),
            h_v_Jmol=float(st_v_aux.h_jmol),
            r_P=float(r_p_aux),
            r_mu=float(r_mu_aux),
            iterations=int(best.iterations),
            fallback_used=True,
            notes=f"{best.notes}; fallback_seed_used",
            jacobian_cond=float(best.jacobian_cond),
            step_norm=float(best.step_norm),
        )

    if best.status != "CONVERGED":
        log_handle.write(f"T={T:.6f} status={best.status}\n")
        for ln in all_traces:
            log_handle.write(f"{ln}\n")
        log_handle.write(
            f"final r_P={best.r_P:.8e} r_mu={best.r_mu:.8e} jac_cond={best.jacobian_cond:.8e} step_norm={best.step_norm:.8e} notes={best.notes}\n\n"
        )

    return best


def saturation_point_at_T(
    fluid: str | Dict,
    T: float,
    tol: float = 1e-10,
    maxiter: int = 50,
    rho_l_guess: Optional[float] = None,
    rho_v_guess: Optional[float] = None,
) -> Tuple[float, float, float, float, float]:
    """
    Solve pure-fluid saturation point at fixed T by equating P and g.

    Inputs
    ------
    fluid : str | dict [unitless]
        Fluid name or already loaded JSON dict.
    T : float [K]
    tol : float [unitless]
        Relative residual tolerance on max(|r_i| / scale_i).
    maxiter : int [unitless]
        Maximum Newton iterations.
    rho_l_guess : float | None [mol/m^3]
    rho_v_guess : float | None [mol/m^3]

    Outputs
    -------
    rho_l, rho_v : float [mol/m^3]
    p_sat : float [Pa]
    h_l, h_v : float [J/mol]

    Failure modes
    -------------
    - RuntimeError if strict acceptance gates are not satisfied.
    """
    data = load_idaes_helmholtz_json(fluid) if isinstance(fluid, str) else fluid
    with Path("/tmp/helmholtz_sat_null.log").open("w") as null_log:
        result = _solve_saturation_with_status(
            data=data,
            T=float(T),
            maxiter=maxiter,
            rho_l_guess=rho_l_guess,
            rho_v_guess=rho_v_guess,
            log_handle=null_log,
        )
    if result.status != "CONVERGED":
        raise RuntimeError(
            f"saturation solve failed at T={T:.6f} K: status={result.status}, "
            f"r_P={result.r_P:.3e}, r_mu={result.r_mu:.3e}, notes={result.notes}"
        )
    return result.rho_l_molm3, result.rho_v_molm3, result.p_Pa, result.h_l_Jmol, result.h_v_Jmol


def compute_saturation_dome(
    fluid: str,
    T_vals: np.ndarray,
    tol: float = 1e-10,
    maxiter: int = 50,
) -> Dict[str, np.ndarray]:
    """
    Compute saturation dome arrays over increasing temperature grid.

    Outputs use mixed reporting units for convenience:
    - p_kPa [kPa]
    - h_l_kJkg, h_v_kJkg [kJ/kg]
    - rho_l_molm3, rho_v_molm3 [mol/m^3]
    """
    _ = tol  # kept for API compatibility; strict acceptance gates are fixed at 1e-6.
    data = load_idaes_helmholtz_json(fluid)
    MW = mw_from_json(data)
    Tc = float(data["basic"]["Tc"])

    rows: List[SaturationResult] = []
    rho_l_guess = None
    rho_v_guess = None
    critical_reached = False

    log_path = Path("verification/log_saturation_solver.txt")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as logf:
        logf.write(f"Saturation solver diagnostics ({datetime.now(timezone.utc).isoformat()})\n")
        logf.write(f"fluid={fluid}\n\n")
        for T in np.asarray(T_vals, dtype=float):
            if T >= Tc:
                break
            result = _solve_saturation_with_status(
                data=data,
                T=float(T),
                maxiter=maxiter,
                rho_l_guess=rho_l_guess,
                rho_v_guess=rho_v_guess,
                log_handle=logf,
            )
            rows.append(result)

            if result.status == "CONVERGED":
                rho_l_guess = result.rho_l_molm3
                rho_v_guess = result.rho_v_molm3
                mean_rho = 0.5 * (result.rho_l_molm3 + result.rho_v_molm3)
                if abs(result.rho_l_molm3 - result.rho_v_molm3) / max(mean_rho, 1.0) < 1e-4:
                    critical_reached = True
                    break
            else:
                rho_l_guess, rho_v_guess = _sat_density_guess_from_aux(data, float(T))

    # Full run rows (converged + failed) for diagnostics.
    all_T = np.array([r.T_K for r in rows], dtype=float)
    all_status = np.array([r.status for r in rows], dtype=object)
    all_rho_l = np.array([r.rho_l_molm3 for r in rows], dtype=float)
    all_rho_v = np.array([r.rho_v_molm3 for r in rows], dtype=float)
    all_p = np.array([r.p_Pa for r in rows], dtype=float)
    all_h_l = np.array([r.h_l_Jmol for r in rows], dtype=float)
    all_h_v = np.array([r.h_v_Jmol for r in rows], dtype=float)
    all_r_p = np.array([r.r_P for r in rows], dtype=float)
    all_r_mu = np.array([r.r_mu for r in rows], dtype=float)
    all_iter = np.array([r.iterations for r in rows], dtype=int)
    all_fallback = np.array([r.fallback_used for r in rows], dtype=bool)
    all_notes = np.array([r.notes for r in rows], dtype=object)
    all_jac_cond = np.array([r.jacobian_cond for r in rows], dtype=float)
    all_step_norm = np.array([r.step_norm for r in rows], dtype=float)

    mask = all_status == "CONVERGED"
    return {
        "fluid": np.array([fluid], dtype=object),
        "T_K": all_T[mask],
        "p_kPa": all_p[mask] * 1e-3,
        "h_l_kJkg": (all_h_l[mask] / MW) * 1e-3,
        "h_v_kJkg": (all_h_v[mask] / MW) * 1e-3,
        "rho_l_molm3": all_rho_l[mask],
        "rho_v_molm3": all_rho_v[mask],
        "critical_reached": np.array([critical_reached], dtype=bool),
        "Tc_K": np.array([Tc], dtype=float),
        "all_T_K": all_T,
        "all_status": all_status,
        "all_rho_l_molm3": all_rho_l,
        "all_rho_v_molm3": all_rho_v,
        "all_p_Pa": all_p,
        "all_h_l_Jmol": all_h_l,
        "all_h_v_Jmol": all_h_v,
        "all_r_P": all_r_p,
        "all_r_mu": all_r_mu,
        "all_iterations": all_iter,
        "all_fallback_used": all_fallback,
        "all_notes": all_notes,
        "all_jacobian_cond": all_jac_cond,
        "all_step_norm": all_step_norm,
        "log_path": np.array([str(log_path)], dtype=object),
    }


def save_dome_csv(dome: Dict[str, np.ndarray], path: str | Path) -> None:
    """Save per-temperature saturation solver status CSV."""
    path = Path(path)
    fieldnames = [
        "T_K",
        "status",
        "rho_l_molm3",
        "rho_v_molm3",
        "p_Pa",
        "h_l_Jmol",
        "h_v_Jmol",
        "r_P",
        "r_mu",
        "iterations",
        "fallback_used",
        "notes",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        n = len(dome["all_T_K"])
        for i in range(n):
            w.writerow(
                {
                    "T_K": float(dome["all_T_K"][i]),
                    "status": str(dome["all_status"][i]),
                    "rho_l_molm3": float(dome["all_rho_l_molm3"][i]),
                    "rho_v_molm3": float(dome["all_rho_v_molm3"][i]),
                    "p_Pa": float(dome["all_p_Pa"][i]),
                    "h_l_Jmol": float(dome["all_h_l_Jmol"][i]),
                    "h_v_Jmol": float(dome["all_h_v_Jmol"][i]),
                    "r_P": float(dome["all_r_P"][i]),
                    "r_mu": float(dome["all_r_mu"][i]),
                    "iterations": int(dome["all_iterations"][i]),
                    "fallback_used": bool(dome["all_fallback_used"][i]),
                    "notes": str(dome["all_notes"][i]),
                }
            )


def plot_ph_dome_from_saturation_data(
    h_l_kJkg: np.ndarray,
    h_v_kJkg: np.ndarray,
    p_kPa: np.ndarray,
    outpath: str | Path,
) -> None:
    """Plot clean P-h saturation dome using CONVERGED points only."""
    outpath = Path(outpath)
    p_bar = np.asarray(p_kPa, dtype=float) / 100.0
    h_l = np.asarray(h_l_kJkg, dtype=float)
    h_v = np.asarray(h_v_kJkg, dtype=float)

    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    ax.plot(h_l, p_bar, lw=2.0, label="Saturated liquid")
    ax.plot(h_v, p_bar, lw=2.0, label="Saturated vapor")
    ax.fill_betweenx(p_bar, h_l, h_v, alpha=0.18, label="Two-phase region")
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("Pure-fluid P-h Saturation Dome")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    if len(p_bar) > 0:
        ic = np.argmax(p_bar)
        ax.plot([h_l[ic], h_v[ic]], [p_bar[ic], p_bar[ic]], "o", ms=4)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")


def plot_ph_dome_with_failures(dome: Dict[str, np.ndarray], outpath: str | Path) -> None:
    """Plot clean dome and overlay FAILED/DIVERGED points as red X markers."""
    outpath = Path(outpath)
    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)

    if len(dome["T_K"]) > 0:
        p_bar = dome["p_kPa"] / 100.0
        ax.plot(dome["h_l_kJkg"], p_bar, lw=2.0, label="Saturated liquid (CONVERGED)")
        ax.plot(dome["h_v_kJkg"], p_bar, lw=2.0, label="Saturated vapor (CONVERGED)")
        ax.fill_betweenx(p_bar, dome["h_l_kJkg"], dome["h_v_kJkg"], alpha=0.15, label="Two-phase region (CONVERGED)")

    failed_mask = dome["all_status"] != "CONVERGED"
    if np.any(failed_mask):
        MW = mw_from_json(load_idaes_helmholtz_json(str(dome["fluid"][0])))
        p_fail_bar = (dome["all_p_Pa"][failed_mask] * 1e-5)
        h_l_fail = (dome["all_h_l_Jmol"][failed_mask] / MW) * 1e-3
        h_v_fail = (dome["all_h_v_Jmol"][failed_mask] / MW) * 1e-3
        ax.plot(h_l_fail, p_fail_bar, "rx", ms=5, label="Failed liquid estimate")
        ax.plot(h_v_fail, p_fail_bar, "rx", ms=5, label="Failed vapor estimate")

    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("Pure-fluid P-h Saturation Dome (with failures)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")


def _cli() -> None:
    """CLI for dome computation and plotting."""
    parser = argparse.ArgumentParser(description="Compute pure-fluid saturation dome from Helmholtz EOS.")
    parser.add_argument("--fluid", required=True, help="fluid stem, e.g. r1234ze")
    parser.add_argument("--Tmin", required=True, type=float, help="minimum temperature [K]")
    parser.add_argument("--Tmax", required=True, type=float, help="maximum temperature [K], must be below Tc")
    parser.add_argument("--n", default=200, type=int, help="number of temperature points")
    parser.add_argument("--tol", default=1e-10, type=float, help="legacy arg; strict acceptance gates are fixed at 1e-6")
    parser.add_argument("--maxiter", default=50, type=int, help="max Newton iterations per T")
    parser.add_argument("--out", default="verification/saturation_dome_run.csv", help="output status CSV path")
    parser.add_argument("--csv", default=None, help="alias for --out")
    parser.add_argument("--fig", default="verification/ph_dome_clean.png", help="clean converged dome plot")
    parser.add_argument("--fig-fail", default="verification/ph_dome_with_failures.png", help="dome plot with failed points")
    parser.add_argument("--metadata", default="verification/saturation_run_metadata.json", help="metadata JSON output")
    args = parser.parse_args()

    out_csv = args.csv if args.csv else args.out
    T_vals = np.linspace(args.Tmin, args.Tmax, args.n)
    dome = compute_saturation_dome(args.fluid, T_vals, tol=args.tol, maxiter=args.maxiter)

    save_dome_csv(dome, out_csv)
    plot_ph_dome_from_saturation_data(dome["h_l_kJkg"], dome["h_v_kJkg"], dome["p_kPa"], args.fig)
    plot_ph_dome_with_failures(dome, args.fig_fail)

    n_attempted = int(len(dome["all_T_K"]))
    n_converged = int(np.sum(dome["all_status"] == "CONVERGED"))
    n_failed = int(n_attempted - n_converged)

    meta = {
        "fluid": args.fluid,
        "Tmin_K": float(args.Tmin),
        "Tmax_K": float(args.Tmax),
        "n_points_requested": int(args.n),
        "n_points_attempted": n_attempted,
        "n_points_converged": n_converged,
        "n_points_failed": n_failed,
        "tol": float(args.tol),
        "acceptance_r_P": 1e-6,
        "acceptance_r_mu": 1e-6,
        "maxiter": int(args.maxiter),
        "critical_reached": bool(dome["critical_reached"][0]),
        "Tc_K": float(dome["Tc_K"][0]),
        "log_file": str(dome["log_path"][0]),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    meta_path = Path(args.metadata)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    with meta_path.open("w") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved CSV: {out_csv}")
    print(f"Saved figure: {args.fig}")
    print(f"Saved failure figure: {args.fig_fail}")
    print(f"Saved metadata: {args.metadata}")
    print(f"Saved solver log: {dome['log_path'][0]}")


if __name__ == "__main__":
    _cli()
