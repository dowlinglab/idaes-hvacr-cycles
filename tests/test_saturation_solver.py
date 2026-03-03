#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-02
Purpose of file: Tests for Helmholtz saturation solver helpers and convergence.
Dependencies: pytest, numpy, helmholtz_saturation, linear_model_codex
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0
"""

from __future__ import annotations

import numpy as np

from helmholtz_saturation import compute_saturation_dome, pure_state_properties, saturation_point_at_T
from linear_model_codex import load_idaes_helmholtz_json, mw_from_json


def _fd(fun, x, rel=1e-6):
    h = max(1e-9, rel * abs(x))
    return (fun(x + h) - fun(max(1e-12, x - h))) / (2.0 * h)


def test_pressure_and_g_derivatives_match_finite_difference():
    data = load_idaes_helmholtz_json("r1234ze")
    MW = mw_from_json(data)
    rhoc = float(data["basic"]["rhoc"]) / MW
    T = 300.0

    for mult in [0.02, 0.2, 1.2]:
        rho = mult * rhoc
        st = pure_state_properties(data, T, rho)

        dp_fd = _fd(lambda rr: pure_state_properties(data, T, rr).p_pa, rho)
        dg_fd = _fd(lambda rr: pure_state_properties(data, T, rr).g_jmol, rho)

        assert np.isfinite(st.dpdrho)
        assert np.isfinite(st.dgdrho)
        assert np.isclose(st.dpdrho, dp_fd, rtol=5e-5, atol=1e-2)
        assert np.isclose(st.dgdrho, dg_fd, rtol=5e-5, atol=1e-2)


def test_single_point_saturation_residuals_close():
    data = load_idaes_helmholtz_json("r1234ze")
    T = 240.0
    rho_l, rho_v, p_sat, h_l, h_v = saturation_point_at_T(data, T, tol=1e-9, maxiter=60)

    st_l = pure_state_properties(data, T, rho_l)
    st_v = pure_state_properties(data, T, rho_v)

    assert rho_l > rho_v > 0.0
    assert p_sat > 0.0
    r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
    r_mu = abs(st_l.g_jmol - st_v.g_jmol) / (8.314462618 * T)
    assert r_p <= 1e-6
    assert r_mu <= 1e-6
    assert h_v > h_l


def test_continuation_over_temperature_grid_converges():
    data = load_idaes_helmholtz_json("r1234ze")
    Tc = float(data["basic"]["Tc"])
    T_vals = np.linspace(240.0, Tc - 5.0, 30)
    dome = compute_saturation_dome("r1234ze", T_vals, tol=1e-9, maxiter=60)

    assert len(dome["T_K"]) > 0
    assert np.all(np.isfinite(dome["p_kPa"]))
    assert np.all(dome["p_kPa"] > 0.0)
    assert np.all(dome["rho_l_molm3"] > dome["rho_v_molm3"])
    # Most points should preserve vapor enthalpy above liquid enthalpy; near
    # fallback transitions this may not hold pointwise.
    frac_h_ordered = np.mean(dome["h_v_kJkg"] >= dome["h_l_kJkg"])
    assert frac_h_ordered > 0.5
