# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa Narasimhan
Creation date: 2026-03-03
Purpose of file: Regression guard to detect accidental edits in frozen true VLE reference implementation.
Dependencies: pytest, mixture_vle_true_reference, linear_model_codex
Context reference: PROJECT_CONTEXT.md
"""

from __future__ import annotations

import math

import mixture_vle_true_reference as ref
from linear_model_codex import load_idaes_helmholtz_json, mw_from_json


def _overall_x1_from_w1(w1: float, mw1: float, mw2: float) -> float:
    """
    Purpose
    -------
    Convert mass fraction to mole fraction for test setup consistency.

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
    - Positive molecular weights and 0 <= w1 <= 1.

    Failure modes
    -------------
    - Division by zero if an invalid molecular weight is supplied.

    References
    ----------
    - Standard mass-to-mole conversion relation.

    Notes on numerical stability
    ----------------------------
    - Stable for bounded, positive inputs.
    """
    n1 = w1 / mw1
    n2 = (1.0 - w1) / mw2
    return n1 / (n1 + n2)


def test_frozen_reference_regression_points() -> None:
    """
    Purpose
    -------
    Validate frozen reference outputs at three fixed points under recorded
    solver settings to catch accidental implementation drift.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    None [unitless]

    Assumptions
    -----------
    - JSON parameter files and frozen solver code are unchanged.

    Failure modes
    -------------
    - AssertionError if any regression value drifts beyond tolerance.

    References
    ----------
    - diagnostics generated during frozen-module creation (2026-03-03).

    Notes on numerical stability
    ----------------------------
    - Uses tight tolerances with deterministic seeds and fixed solver profile.
    """
    d1 = load_idaes_helmholtz_json("r1234ze")
    d2 = load_idaes_helmholtz_json("r227ea")
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = _overall_x1_from_w1(0.911, mw1, mw2)

    # Bubble profile used in preserved diagnostics
    ref.LSQ_METHOD = "trf"
    ref.LSQ_DIFF_STEP = 1e-6
    ref.LSQ_MAX_NFEV = 800
    ref.LSQ_X_SCALE = 1.0
    ref.LSQ_FTOL = 1e-14
    ref.LSQ_XTOL = 1e-14
    ref.LSQ_GTOL = 1e-14
    ref.LSQ_LOSS = "cauchy"
    ref.LSQ_F_SCALE = 10.0

    b300 = ref.solve_bubble_at_t(d1, d2, 300.0, z1, 5000.0, 50.0, z1)
    b340 = ref.solve_bubble_at_t(
        d1,
        d2,
        340.0,
        z1,
        b300["rho_l_molm3"],
        b300["rho_v_molm3"],
        b300["y1_vap"],
    )

    # Dew profile used in preserved diagnostics
    ref.LSQ_METHOD = "trf"
    ref.LSQ_DIFF_STEP = 1e-8
    ref.LSQ_MAX_NFEV = 800
    ref.LSQ_X_SCALE = 2.0
    ref.LSQ_FTOL = 1e-10
    ref.LSQ_XTOL = 1e-10
    ref.LSQ_GTOL = 1e-10
    ref.LSQ_LOSS = "huber"
    ref.LSQ_F_SCALE = 3.0

    d300 = ref.solve_dew_at_t(d1, d2, 300.0, z1, 5000.0, 50.0, z1)

    assert b300["status"] == "CONVERGED"
    assert b340["status"] == "CONVERGED"
    assert d300["status"] == "CONVERGED"

    assert math.isclose(b300["P_Pa"], 683941.3593653208, rel_tol=1e-10, abs_tol=1e-6)
    assert math.isclose(b340["P_Pa"], 1582906.5135530196, rel_tol=1e-10, abs_tol=1e-6)
    assert math.isclose(d300["P_Pa"], 681251.2806882628, rel_tol=1e-10, abs_tol=1e-6)

    assert math.isclose(b300["y1_vap"], 0.9504784105956765, rel_tol=1e-11, abs_tol=1e-12)
    assert math.isclose(b340["y1_vap"], 0.9441794516234248, rel_tol=1e-11, abs_tol=1e-12)
    assert math.isclose(d300["x1_liq"], 0.9239396916178307, rel_tol=1e-11, abs_tol=1e-12)
