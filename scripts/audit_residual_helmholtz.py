#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-03
Purpose of file: Audit residual Helmholtz implementation and derivative consistency
for pure-fluid, departure-term, and mixture-residual pathways.
Dependencies: numpy, csv, pathlib, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-03
# Assumptions:
# - Audit targets the currently hardwired R1234ze(E)/R227ea Bell mixture model.
# - Finite-difference checks are local consistency checks, not external validation.
# - Composition is fixed for tau/delta derivative checks unless explicitly noted.
# TODO: Add external-reference validation audit once REFPROP/NIST baseline table is finalized.
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from linear_model_codex import (  # noqa: E402
    BELL_2023_R1234ZE_R227EA,
    _nares_from_n,
    _mixture_alpha_eval,
    _mixture_reduced_state,
    _bell2023_reducing_derivs_binary,
    alphar_idaes_with_derivs,
    bell2023_departure_alphar,
    bell2023_departure_base,
    bell2023_Tred_vred,
    load_idaes_helmholtz_json,
    mixture_alpha0_alphar_derivs,
    mw_from_json,
)


EPS = 1e-12


@dataclass
class AuditRow:
    """
    Purpose
    -------
    Record one derivative audit check result.

    Inputs
    ------
    check : str [unitless]
    id : str [unitless]
    metric : str [unitless]
    analytic : float [various]
    finite_diff : float [various]
    abs_err : float [same as metric]
    rel_err : float [unitless]

    Outputs
    -------
    AuditRow [dataclass]

    Assumptions
    -----------
    - analytic and finite_diff correspond to the same mathematical quantity.

    Failure modes
    -------------
    - None; storage-only container.

    References
    ----------
    - Internal audit schema for residual Helmholtz checks.

    Notes on numerical stability
    ----------------------------
    - Not applicable.
    """

    check: str
    id: str
    metric: str
    analytic: float
    finite_diff: float
    abs_err: float
    rel_err: float


def _fd1(func: Callable[[float], float], x: float, h_rel: float = 1e-6) -> float:
    """
    Purpose
    -------
    Compute first derivative by central finite difference.

    Inputs
    ------
    func : callable [unitless]
    x : float [function argument units]
    h_rel : float [unitless]

    Outputs
    -------
    dfunc_dx : float [func units per x units]

    Assumptions
    -----------
    - func is smooth in a neighborhood around x.

    Failure modes
    -------------
    - Numerical cancellation if h is too small.

    References
    ----------
    - Standard central finite-difference derivative.

    Notes on numerical stability
    ----------------------------
    - Step uses max(1e-8, h_rel*|x|) floor.
    """
    h = max(1e-8, h_rel * abs(x))
    xp = x + h
    xm = max(EPS, x - h)
    return (func(xp) - func(xm)) / (xp - xm)


def _rel_err(a: float, b: float) -> float:
    """
    Purpose
    -------
    Compute robust relative error.

    Inputs
    ------
    a : float [units of compared quantity]
    b : float [units of compared quantity]

    Outputs
    -------
    rel : float [unitless]

    Assumptions
    -----------
    - At least one of a or b is finite.

    Failure modes
    -------------
    - Returns NaN if either value is NaN.

    References
    ----------
    - Relative error definition with floor denominator.

    Notes on numerical stability
    ----------------------------
    - Denominator floor avoids blow-up near zero reference values.
    """
    denom = max(abs(a), abs(b), 1e-14)
    return abs(a - b) / denom


def _pure_residual_audit() -> List[AuditRow]:
    """
    Purpose
    -------
    Audit pure-fluid alphar first derivatives against finite differences.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    rows : list[AuditRow]

    Assumptions
    -----------
    - States sampled are away from singular/critical boundaries.

    Failure modes
    -------------
    - Propagates exceptions from EOS evaluation.

    References
    ----------
    - alphar_idaes_with_derivs implementation in linear_model_codex.py.

    Notes on numerical stability
    ----------------------------
    - Uses moderate tau/delta sample points to limit FD noise.
    """
    rows: List[AuditRow] = []
    fluids = ["r1234ze", "r227ea"]
    taus = [0.75, 0.95, 1.15]
    deltas = [0.05, 0.3, 0.8, 1.4]

    for fluid in fluids:
        d = load_idaes_helmholtz_json(fluid)
        eos = d["eos"]
        for tau in taus:
            for delta in deltas:
                ar, ar_tau, ar_del = alphar_idaes_with_derivs(eos, tau, delta)
                fd_tau = _fd1(lambda tt: alphar_idaes_with_derivs(eos, tt, delta)[0], tau)
                fd_del = _fd1(lambda dd: alphar_idaes_with_derivs(eos, tau, dd)[0], delta)
                rows.append(
                    AuditRow(
                        check="pure_alphar",
                        id=f"{fluid}|tau={tau}|delta={delta}",
                        metric="d_ar_d_tau",
                        analytic=float(ar_tau),
                        finite_diff=float(fd_tau),
                        abs_err=abs(float(ar_tau) - float(fd_tau)),
                        rel_err=_rel_err(float(ar_tau), float(fd_tau)),
                    )
                )
                rows.append(
                    AuditRow(
                        check="pure_alphar",
                        id=f"{fluid}|tau={tau}|delta={delta}",
                        metric="d_ar_d_delta",
                        analytic=float(ar_del),
                        finite_diff=float(fd_del),
                        abs_err=abs(float(ar_del) - float(fd_del)),
                        rel_err=_rel_err(float(ar_del), float(fd_del)),
                    )
                )
                _ = ar
    return rows


def _departure_residual_audit() -> List[AuditRow]:
    """
    Purpose
    -------
    Audit Bell departure term derivatives against finite differences.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    rows : list[AuditRow]

    Assumptions
    -----------
    - Uses R1234ze(E)/R227ea departure coefficients currently hardwired.

    Failure modes
    -------------
    - Propagates exceptions from departure evaluation.

    References
    ----------
    - bell2023_departure_alphar implementation in linear_model_codex.py.

    Notes on numerical stability
    ----------------------------
    - Composition values avoid x=0 or x=1 endpoint singular behavior.
    """
    rows: List[AuditRow] = []
    x1_vals = [0.2, 0.5, 0.8]
    taus = [0.75, 0.95, 1.15]
    deltas = [0.05, 0.3, 0.8, 1.4]

    for x1 in x1_vals:
        x2 = 1.0 - x1
        for tau in taus:
            for delta in deltas:
                ar, ar_tau, ar_del = bell2023_departure_alphar(x1, x2, tau, delta)
                fd_tau = _fd1(lambda tt: bell2023_departure_alphar(x1, x2, tt, delta)[0], tau)
                fd_del = _fd1(lambda dd: bell2023_departure_alphar(x1, x2, tau, dd)[0], delta)
                rows.append(
                    AuditRow(
                        check="departure_alphar",
                        id=f"x1={x1}|tau={tau}|delta={delta}",
                        metric="d_dep_d_tau",
                        analytic=float(ar_tau),
                        finite_diff=float(fd_tau),
                        abs_err=abs(float(ar_tau) - float(fd_tau)),
                        rel_err=_rel_err(float(ar_tau), float(fd_tau)),
                    )
                )
                rows.append(
                    AuditRow(
                        check="departure_alphar",
                        id=f"x1={x1}|tau={tau}|delta={delta}",
                        metric="d_dep_d_delta",
                        analytic=float(ar_del),
                        finite_diff=float(fd_del),
                        abs_err=abs(float(ar_del) - float(fd_del)),
                        rel_err=_rel_err(float(ar_del), float(fd_del)),
                    )
                )
                _ = ar
    return rows


def _mixture_residual_tau_delta_audit() -> List[AuditRow]:
    """
    Purpose
    -------
    Audit mixture residual ar_tau/ar_del mapping vs finite differences.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    rows : list[AuditRow]

    Assumptions
    -----------
    - Composition fixed during tau/delta finite-difference checks.

    Failure modes
    -------------
    - Propagates exceptions from mixture derivative evaluation.

    References
    ----------
    - mixture_alpha0_alphar_derivs in linear_model_codex.py.

    Notes on numerical stability
    ----------------------------
    - Uses fixed reducing scales at constant composition.
    """
    rows: List[AuditRow] = []
    d1 = load_idaes_helmholtz_json("r1234ze")
    d2 = load_idaes_helmholtz_json("r227ea")
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    x1_vals = [0.2, 0.5, 0.8]
    taus = [0.75, 0.95, 1.15]
    deltas = [0.05, 0.3, 0.8, 1.4]

    for x1 in x1_vals:
        x2 = 1.0 - x1
        tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
        rho_red = 1.0 / vred

        def ar_func(tt: float, dd: float) -> float:
            return mixture_alpha0_alphar_derivs(
                d1=d1,
                d2=d2,
                x1=x1,
                x2=x2,
                tau=tt,
                delta=dd,
                Tred=tred,
                rho_red_mol=rho_red,
                pair_key="r1234ze|r227ea",
            )[2]

        for tau in taus:
            for delta in deltas:
                _, _, ar_mix, ar_tau, ar_del = mixture_alpha0_alphar_derivs(
                    d1=d1,
                    d2=d2,
                    x1=x1,
                    x2=x2,
                    tau=tau,
                    delta=delta,
                    Tred=tred,
                    rho_red_mol=rho_red,
                    pair_key="r1234ze|r227ea",
                )
                fd_tau = _fd1(lambda tt: ar_func(tt, delta), tau)
                fd_del = _fd1(lambda dd: ar_func(tau, dd), delta)
                rows.append(
                    AuditRow(
                        check="mixture_ar_tau_delta",
                        id=f"x1={x1}|tau={tau}|delta={delta}",
                        metric="d_ar_mix_d_tau",
                        analytic=float(ar_tau),
                        finite_diff=float(fd_tau),
                        abs_err=abs(float(ar_tau) - float(fd_tau)),
                        rel_err=_rel_err(float(ar_tau), float(fd_tau)),
                    )
                )
                rows.append(
                    AuditRow(
                        check="mixture_ar_tau_delta",
                        id=f"x1={x1}|tau={tau}|delta={delta}",
                        metric="d_ar_mix_d_delta",
                        analytic=float(ar_del),
                        finite_diff=float(fd_del),
                        abs_err=abs(float(ar_del) - float(fd_del)),
                        rel_err=_rel_err(float(ar_del), float(fd_del)),
                    )
                )
                _ = ar_mix
    return rows


def _composition_derivative_audit() -> List[AuditRow]:
    """
    Purpose
    -------
    Audit analytic d(n*ar)/dn_i against finite differences at fixed T,V,n_j.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    rows : list[AuditRow]

    Assumptions
    -----------
    - Binary system and analytic relations from current compute_table1_properties.

    Failure modes
    -------------
    - Propagates exceptions from EOS and derivative evaluations.

    References
    ----------
    - linear_model_codex.compute_table1_properties composition derivative block.

    Notes on numerical stability
    ----------------------------
    - Uses small but bounded finite-difference steps on mole numbers.
    """
    rows: List[AuditRow] = []
    d1 = load_idaes_helmholtz_json("r1234ze")
    d2 = load_idaes_helmholtz_json("r227ea")
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    test_states = [
        (280.0, 150.0, 0.30),
        (320.0, 2200.0, 0.50),
        (360.0, 3500.0, 0.80),
    ]

    for T, rho_mol, x1 in test_states:
        x2 = 1.0 - x1
        vals = _mixture_alpha_eval(d1, d2, x1, T, rho_mol)
        tau = vals["tau"]
        delta = vals["delta"]
        ar = vals["ar_mix"]
        ar_del = vals["ar_del_mix"]

        tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
        dTred_dx1, dvred_dx1 = _bell2023_reducing_derivs_binary(x1, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
        dtau_dx1 = dTred_dx1 / T
        ddelta_dx1 = rho_mol * dvred_dx1

        _, _, tau1, tau2, delta1, delta2, _, _ = _mixture_reduced_state(d1, d2, x1, T, rho_mol)
        ar1, _, _ = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)
        ar2, _, _ = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)
        dep_base = bell2023_departure_base(tau, delta)
        _, dep_tau, dep_del = bell2023_departure_alphar(x1, x2, tau, delta)
        ddep_dx1 = (x2 - x1) * dep_base + dep_tau * dtau_dx1 + dep_del * ddelta_dx1
        dar_dx1 = (ar1 - ar2) + ddep_dx1
        dar_drho = ar_del * vred

        d_na_dn1 = ar + rho_mol * dar_drho + x2 * dar_dx1
        d_na_dn2 = ar + rho_mol * dar_drho - x1 * dar_dx1

        # Finite-difference reference at constant T, V, n_j.
        n_tot = 1.0
        n1 = x1 * n_tot
        n2 = x2 * n_tot
        V = n_tot / rho_mol
        dn1 = max(1e-8, 1e-6 * n1)
        dn2 = max(1e-8, 1e-6 * n2)

        fd_dn1 = (_nares_from_n(d1, d2, T, V, n1 + dn1, n2) - _nares_from_n(d1, d2, T, V, max(EPS, n1 - dn1), n2)) / (
            (n1 + dn1) - max(EPS, n1 - dn1)
        )
        fd_dn2 = (_nares_from_n(d1, d2, T, V, n1, n2 + dn2) - _nares_from_n(d1, d2, T, V, n1, max(EPS, n2 - dn2))) / (
            (n2 + dn2) - max(EPS, n2 - dn2)
        )

        rows.append(
            AuditRow(
                check="composition_d_nar_dn",
                id=f"T={T}|rho={rho_mol}|x1={x1}",
                metric="d_nar_dn1",
                analytic=float(d_na_dn1),
                finite_diff=float(fd_dn1),
                abs_err=abs(float(d_na_dn1) - float(fd_dn1)),
                rel_err=_rel_err(float(d_na_dn1), float(fd_dn1)),
            )
        )
        rows.append(
            AuditRow(
                check="composition_d_nar_dn",
                id=f"T={T}|rho={rho_mol}|x1={x1}",
                metric="d_nar_dn2",
                analytic=float(d_na_dn2),
                finite_diff=float(fd_dn2),
                abs_err=abs(float(d_na_dn2) - float(fd_dn2)),
                rel_err=_rel_err(float(d_na_dn2), float(fd_dn2)),
            )
        )

    return rows


def _summarize(rows: List[AuditRow]) -> Dict[str, Dict[str, float]]:
    """
    Purpose
    -------
    Aggregate error statistics per audit check.

    Inputs
    ------
    rows : list[AuditRow]

    Outputs
    -------
    summary : dict[str, dict[str, float]]

    Assumptions
    -----------
    - rows contain finite error values.

    Failure modes
    -------------
    - Empty group returns NaN statistics.

    References
    ----------
    - Internal audit reporting format.

    Notes on numerical stability
    ----------------------------
    - Uses nan-safe NumPy aggregations.
    """
    summary: Dict[str, Dict[str, float]] = {}
    checks = sorted({r.check for r in rows})
    for chk in checks:
        rsub = [r for r in rows if r.check == chk]
        rel = np.array([r.rel_err for r in rsub], dtype=float)
        abs_e = np.array([r.abs_err for r in rsub], dtype=float)
        summary[chk] = {
            "n": float(len(rsub)),
            "mean_rel_err": float(np.nanmean(rel)),
            "max_rel_err": float(np.nanmax(rel)),
            "mean_abs_err": float(np.nanmean(abs_e)),
            "max_abs_err": float(np.nanmax(abs_e)),
        }
    return summary


def main() -> None:
    """
    Purpose
    -------
    Execute residual Helmholtz audit and save detailed and summary CSV outputs.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    Writes:
      verification/residual_helmholtz_audit_details.csv
      verification/residual_helmholtz_audit_summary.csv

    Assumptions
    -----------
    - Required fluid parameter files are available in IDAES installation.

    Failure modes
    -------------
    - Raises if EOS calls fail unexpectedly.

    References
    ----------
    - linear_model_codex residual Helmholtz and mixture derivative implementations.

    Notes on numerical stability
    ----------------------------
    - This is an internal consistency audit against finite differences.
    """
    rows: List[AuditRow] = []
    rows.extend(_pure_residual_audit())
    rows.extend(_departure_residual_audit())
    rows.extend(_mixture_residual_tau_delta_audit())
    rows.extend(_composition_derivative_audit())

    out_dir = Path("verification")
    out_dir.mkdir(parents=True, exist_ok=True)
    details_path = out_dir / "residual_helmholtz_audit_details.csv"
    summary_path = out_dir / "residual_helmholtz_audit_summary.csv"

    with details_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["check", "id", "metric", "analytic", "finite_diff", "abs_err", "rel_err"],
        )
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "check": r.check,
                    "id": r.id,
                    "metric": r.metric,
                    "analytic": f"{r.analytic:.17e}",
                    "finite_diff": f"{r.finite_diff:.17e}",
                    "abs_err": f"{r.abs_err:.17e}",
                    "rel_err": f"{r.rel_err:.17e}",
                }
            )

    summary = _summarize(rows)
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["check", "n", "mean_rel_err", "max_rel_err", "mean_abs_err", "max_abs_err"],
        )
        w.writeheader()
        for chk, stats in summary.items():
            w.writerow({"check": chk, **stats})

    print(f"Saved details: {details_path}")
    print(f"Saved summary: {summary_path}")
    for chk, stats in summary.items():
        print(
            f"{chk}: n={int(stats['n'])}, mean_rel={stats['mean_rel_err']:.3e}, "
            f"max_rel={stats['max_rel_err']:.3e}, max_abs={stats['max_abs_err']:.3e}"
        )


if __name__ == "__main__":
    main()
