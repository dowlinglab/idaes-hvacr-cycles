"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5), Claude AI
QA/Testing Responsibility: Shilpa Narasimhan
Creation date: 2026-03-02
Edit date: 2026-08-13
Purpose of file: Compute true binary VLE envelopes (bubble/dew) for a
fixed-overall-composition refrigerant blend using Helmholtz EOS mixture model.
Dependencies: numpy, scipy, matplotlib, linear_model_codex

This file will serve as the validation for Vanilla code. Meaning that
if the dew point and bubble points are validated against Honeywell data,
the code is assumed to be validated.

Version: v0.2.0

This file computes the saturation dome ONLY (bubble line = saturated
liquid, dew line = saturated vapor), from the self-consistent 3-equation
VLE solve (pressure equality + component chemical-potential equality;
see the KEY EQUATIONS REFERENCE block below, tag M6). It does not compute
or represent interior two-phase (quality/lever-rule) states -- that is a
distinct, out-of-scope calculation belonging to a flash/cycle model, not
to dome validation.

# BREADCRUMB:
# Date: 2026-03-02
# Assumptions:
# - Binary VLE solved from P equality and component chemical-potential equality.
# - Chemical potentials are evaluated from fugacity using analytic composition
#   derivatives of n*alpha^r.
# - Bubble curve fixes liquid composition x=z; dew curve fixes vapor composition y=z.
# - This module is isolated from compute_pressure_enthalpy workflow.
# TODO: Add optional finite-difference mu_i cross-check for debugging.
# TODO: Keep thermodynamic identities cited in docstrings when equations are added.
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
from scipy.optimize import brentq, least_squares, root

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


# =============================================================================
# KEY EQUATIONS REFERENCE (read this once; every function below tags its lines
# back to one of these, e.g. "# (M6)")
#
# Notation: subscript 1 = R-1234ze(E), subscript 2 = R-227ea. x1,x2 = liquid
# mole fractions (x1+x2=1); y1,y2 = vapor mole fractions. T [K] = temperature.
# rho [mol/m^3] = molar density. R_u [J/(mol.K)] = universal gas constant.
#
# (M1) Reduced variables (inputs to every alpha evaluation):
#        tau   = T_red(x1) / T          (reduced inverse temperature)
#        delta = rho * v_red(x1)        (reduced density; v_red = 1/rho_red)
#      T_red(x1), v_red(x1) are the Bell (2023) binary reducing functions
#      (paper Eq. 3-5): quadratic mixing rule in x1,x2 plus a beta/gamma
#      "cross" term that lets the mixture's reducing point deviate from a
#      simple mole-fraction average. Implemented in bell2023_Tred_vred()
#      (imported from linear_model_codex, already validated separately).
#
# (M2) Non-dimensional Helmholtz energy of the mixture (Bell 2023 Eq. 1):
#        alpha_mix(tau,delta,x1) = alpha0_mix + alphar_mix
#      alpha0_mix = ideal-gas part, INCLUDING the entropy-of-mixing term
#        alpha0_mix = alpha0_CS + x1*ln(x1) + x2*ln(x2)          (M2a)
#      NOTE: x1 -> 0 or x1 -> 1 makes ln(x1) or ln(x2) blow up. Every call
#      site clips x1 into (EPS_X, 1-EPS_X) via _clip_x() BEFORE it ever
#      reaches (M2a), so ln(0) never actually happens -- EPS_X itself is
#      just a module-level constant defined right below this block (search
#      "EPS_X = "); it is not yet in scope at THIS point in the file, only
#      at runtime once the module has finished loading top to bottom.
#      alphar_mix = residual part = corresponding-states (CS) contribution
#        from the two pure fluids PLUS a fitted departure function
#        (Bell 2023 Eq. 6-7, Table 7 coefficients):
#        alphar_mix = alphar_CS(tau,delta,x1) + alphar_dep(tau,delta,x1,x2)
#      Both a0/ar and their tau- and delta-derivatives come out of
#      mixture_alpha0_alphar_derivs() (imported, already validated).
#
# (M3) Compressibility factor and pressure (standard Helmholtz-EOS identity):
#        Z = 1 + delta * (d alphar_mix / d delta)_tau,x            [alphar_del]
#        P = rho * R_u * T * Z
#      NOTE: on the deep-liquid branch, Z is a near-total cancellation
#      between "1" and a large negative alphar_del*delta term (Z ~ 0.2 at
#      R-515B's real liquid state) -- this is why pressure is extremely
#      sensitive to any error in alphar_del there. Never substitute a real
#      / externally-measured rho into this formula directly; always let the
#      3-equation solve below determine rho self-consistently first.
#
# (M4) Molar enthalpy and Gibbs energy identities (from alpha and its
#      tau-derivative alpha_tau = d(alpha_mix)/d(tau) at fixed delta,x):
#        h/(R_u T) = 1 + tau*alpha_tau + delta*alphar_del
#        g/(R_u T) = 1 + alpha_mix     + delta*alphar_del
#
# (M5) Fugacity / chemical potential of each component (Lemmon &
#      Tillner-Roth mixture fugacity form). d(n*alphar)/dn_i is the
#      composition (mole-number) derivative of the *residual* Helmholtz
#      energy at fixed T,V, obtained by the chain rule through
#      tau(x1), delta(x1), and the direct x1-dependence of alphar itself:
#        d(n alphar)/dn1 = alphar_mix + rho*(d alphar/d rho) + x2*(d alphar/d x1)
#        d(n alphar)/dn2 = alphar_mix + rho*(d alphar/d rho) - x1*(d alphar/d x1)
#      then:
#        f_i = x_i * rho * R_u * T * exp[ d(n alphar)/dn_i ]
#        mu_i = R_u * T * ln(f_i)
#
# (M6) Bubble/dew 3x3 residual system solved by least_squares, at FIXED T
#      and fixed overall/feed composition z1. Unknowns: rho_l, rho_v, and
#      the free-phase composition (y1 for bubble, x1 for dew).
#        r1 = [P(T,rho_l,x_liq) - P(T,rho_v,x_vap)] / avg(P)   mechanical eq.
#        r2 = [mu1_liq - mu1_vap] / (R_u T)                    component-1 eq.
#        r3 = [mu2_liq - mu2_vap] / (R_u T)                    component-2 eq.
#      The P reported for the dome at this T is P from the CONVERGED
#      solution -- never a P computed by plugging in an assumed rho.
# =============================================================================

# =============================================================================
# FULL REFERENCES (cited by tag, e.g. "[B23]", in every function's docstring
# References section below -- this is the one place the complete bibliographic
# entry lives, so no docstring below needs to repeat it in full).
#
# [B23]   Bell, I. H. (2023). "Mixture Model for Refrigerant Pairs
#         R-32/1234yf, R-32/1234ze(E), R-1234ze(E)/227ea, R-1234yf/152a, and
#         R-125/1234yf." Journal of Physical and Chemical Reference Data,
#         52(1), 013101. https://doi.org/10.1063/5.0135368
#         (Source of the binary reducing functions (M1, Eq. 3-5), the
#         departure function (M2, Eq. 6-7), and the Table 7 / Table 2
#         coefficients used for the R-1234ze(E)/R-227ea pair -- the surrogate
#         binary this file uses to represent pseudo-pure R-515B.)
#
# [LT99]  Lemmon, E. W., & Tillner-Roth, R. (1999). "A Helmholtz energy
#         equation of state for calculating the thermodynamic properties of
#         fluid mixtures." Fluid Phase Equilibria, 165(1), 1-21.
#         https://doi.org/10.1016/S0378-3812(99)00262-9
#         (General corresponding-states-plus-departure-function Helmholtz
#         mixture framework this file's EOS structure follows: the P/h/g
#         identities in M3-M4, and the fugacity/chemical-potential form used
#         to build the M5/M6 residuals for the VLE solve.)
#
# [BCL99] Branch, M. A., Coleman, T. F., & Li, Y. (1999). "A Subspace,
#         Interior, and Conjugate Gradient Method for Large-Scale
#         Bound-Constrained Minimization Problems." SIAM Journal on
#         Scientific Computing, 21(1), 1-23.
#         https://doi.org/10.1137/S1064827595289108
#         (The Trust Region Reflective ["trf"] algorithm used by
#         scipy.optimize.least_squares below -- see LSQ_METHOD.)
# =============================================================================

EPS_X = 1e-12
RHO_MIN_MOLM3 = 1e-9  # hard floor on molar density [mol/m^3]; keeps log/exp args finite
RHO_MAX_MOLM3 = 2.0e4  # hard ceiling on molar density [mol/m^3]; well above any physical liquid rho for this pair
# Scaled-sigmoid bounds for density parameterization.
# 0.001..2.0 g/cc ~= 1..2000 kg/m^3 mapped to molar-density bounds using
# conservative component MW limits for R1234ze/R227ea.
# WHY: least_squares works in an UNCONSTRAINED variable space (u0,u1,u2 in
# solve_bubble_at_t/solve_dew_at_t below). Rather than adding explicit
# inequality constraints (rho_v>0, rho_l>rho_v, 0<y1<1) to the optimizer,
# we map each unconstrained real number through a sigmoid into the physically
# valid range -- so *every* trial point the solver tries is automatically
# physical, with no separate feasibility check needed. RHO_MAP_* bounds the
# vapor density; DRHO_MAP_* bounds the *increment* (rho_l - rho_v), which is
# what guarantees rho_l > rho_v everywhere (see _rho_param_to_states below).
RHO_MAP_MIN_MOLM3 = 5.0
RHO_MAP_MAX_MOLM3 = 2.0e4
DRHO_MAP_MIN_MOLM3 = 1.0
DRHO_MAP_MAX_MOLM3 = 2.0e4
# scipy.optimize.least_squares tuning (kept fixed so every dome point uses
# an identical, already-validated solver configuration -- see the
# 2026-03-03 "bubble-tuned" breadcrumbs for why these specific values were
# chosen: trf handles the explicit bounds cleanly; diff_step sets the
# finite-difference step for the numerically-differentiated Jacobian since
# no analytic Jacobian of the residuals is provided).
LSQ_DIFF_STEP = 1e-6  # finite-difference step size used to build the Jacobian of res(u)
LSQ_MAX_NFEV = 800  # max residual-function evaluations before giving up
LSQ_METHOD = "trf"  # Trust Region Reflective -- supports the box bounds (lb,ub) below
LSQ_X_SCALE = 1.0  # relative scaling between the 3 unknowns (kept uniform; they're already O(1) after the sigmoid map)
LSQ_FTOL = 1.0e-14  # convergence tolerance on residual (cost) change
LSQ_XTOL = 1.0e-14  # convergence tolerance on solution-vector change
LSQ_GTOL = 1.0e-14  # convergence tolerance on gradient norm
LSQ_LOSS = "huber"  # switched from "cauchy" 2026-08-13: cauchy's aggressive down-weighting let the bubble solver settle early on a nearby-but-wrong basin near T~278K, producing a non-monotonic h_l turnover (S-kink). huber matches the dew branch's already-working profile.
LSQ_F_SCALE = 3.0  # matched to dew's LSQ_F_SCALE_DEW; scale parameter for the loss function


LSQ_METHOD_DEW = "trf"
LSQ_DIFF_STEP_DEW = 1e-8
LSQ_MAX_NFEV_DEW = 800
LSQ_X_SCALE_DEW = 2.0
LSQ_FTOL_DEW = 1.0e-10
LSQ_XTOL_DEW = 1.0e-10
LSQ_GTOL_DEW = 1.0e-10
LSQ_LOSS_DEW = "huber"
LSQ_F_SCALE_DEW = 3.0

# -----------------------------------------------------------------------
# Tapered density-separation acceptance gate (added 2026-08-13).
#
# WHY: the acceptance check requires rho_l > rho_v * RATIO to reject the
# degenerate "fake" root where the solver finds two nearly-identical
# densities instead of true phase separation (see the 2026-08-13
# breadcrumb entries on the bubble-branch S-kink). A fixed RATIO=3.0
# works everywhere EXCEPT near the critical point, where rho_l -> rho_v
# is the correct physical behavior (that is the literal definition of a
# critical point) -- a fixed ratio of 3.0 there would reject the TRUE
# solution, not just fake ones, silently truncating the dome before
# bubble and dew can meet.
#
# FIX: taper the required ratio linearly from RHO_SEP_RATIO_FAR (far from
# critical) down to RHO_SEP_RATIO_NEAR_TC (at/above the mixture's
# reducing temperature Tred_mix, used here as a model-consistent proxy
# for the mixture critical temperature -- already confirmed elsewhere in
# this project to match the Honeywell datasheet's 382.04 K almost
# exactly at this composition). The floor stays strictly > 1.0 so the
# gate still rejects a literal rho_l==rho_v degenerate root even at Tc.
# -----------------------------------------------------------------------
RHO_SEP_RATIO_FAR = 3.0  # required rho_l/rho_v ratio at/below (Tred_mix - RHO_SEP_TAPER_START_K)
RHO_SEP_RATIO_NEAR_TC = 1.05  # relaxed floor allowed right at/above Tred_mix; still > 1.0 to reject exact degeneracy
RHO_SEP_TAPER_START_K = 30.0  # begin relaxing the ratio this many K below Tred_mix


def _rho_separation_min_ratio(t_k: float, tred_mix: float) -> float:
    """
    Purpose
    -------
    Compute the temperature-dependent minimum acceptable rho_l/rho_v ratio
    for the bubble/dew acceptance gate, tapering toward 1.0 near the
    mixture's critical region so true near-critical solutions are not
    rejected alongside fake degenerate ones.

    Inputs
    ------
    t_k : float [K]
      Current solve temperature.
    tred_mix : float [K]
      Mixture reducing temperature at the fixed feed composition (from
      bell2023_Tred_vred), used as a proxy for the mixture critical
      temperature.

    Outputs
    -------
    min_ratio : float [unitless]
      Required rho_l/rho_v ratio for this T; RHO_SEP_RATIO_FAR far from
      Tred_mix, linearly relaxing to RHO_SEP_RATIO_NEAR_TC within
      RHO_SEP_TAPER_START_K of it.

    Assumptions
    -----------
    - Linear taper is a simple, monotonic, easily-tuned approximation --
      not a literature critical-scaling-law fit (real near-critical
      density separation follows (Tc-T)^beta with beta~0.325 for real
      fluids, but a linear ratio taper is sufficient here since this gate
      is a numerical-acceptance heuristic, not a physical model).

    Failure modes
    -------------
    - None (pure arithmetic, always returns a finite value).
    """
    dt = tred_mix - t_k  # distance below the mixture reducing/critical temperature [K]; can go negative if t_k > tred_mix
    if dt >= RHO_SEP_TAPER_START_K:
        return RHO_SEP_RATIO_FAR
    if dt <= 0.0:
        return RHO_SEP_RATIO_NEAR_TC
    frac = dt / RHO_SEP_TAPER_START_K  # 1.0 far from Tred_mix -> 0.0 right at Tred_mix
    return RHO_SEP_RATIO_NEAR_TC + frac * (RHO_SEP_RATIO_FAR - RHO_SEP_RATIO_NEAR_TC)

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
    - [LT99] (see file-level FULL REFERENCES block) -- this container holds
      the P/h/g outputs of that general Helmholtz-mixture EOS framework.

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
    - Not an external literature citation -- this is a direct algebraic
      identity (n_i = w_i/MW_i, then x1 = n1/(n1+n2)), not drawn from any of
      the file-level FULL REFERENCES.

    Notes on numerical stability
    ----------------------------
    - Stable for positive molecular weights and bounded mass fractions.
    """
    # Standard mass-fraction -> mole-fraction conversion: moles per unit mass
    # of each component (n_i = w_i/MW_i for 1 kg total), then normalize.
    w2 = 1.0 - w1  # mass fraction of fluid 2
    n1 = w1 / mw1  # "moles" of fluid 1 per kg of mixture (w1 [kg/kg] / mw1 [kg/mol] = [mol/kg])
    n2 = w2 / mw2  # "moles" of fluid 2 per kg of mixture
    return float(n1 / (n1 + n2))  # x1 = n1 / (n1+n2)


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
    - Not an external literature citation -- this is an internal numerical
      guard for the x1*ln(x1)+x2*ln(x2) entropy-of-mixing term (M2a, part of
      the [LT99]/[B23] ideal-mixing formulation), not a citation itself.

    Notes on numerical stability
    ----------------------------
    - Prevents overflow/NaN from log(0) in entropy-of-mixing terms.
    """
    return float(min(max(x1, EPS_X), 1.0 - EPS_X))  # clamp to (EPS_X, 1-EPS_X) so ln(x1) and ln(x2)=ln(1-x1) in (M2a) never hit ln(0)


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
    - [B23] -- binary reducing functions (M1, Eq. 3-5) and departure
      function (M2, Eq. 6-7).
    - [LT99] -- overall corresponding-states + departure-function Helmholtz
      mixture structure (alpha_mix = alpha0_mix + alphar_mix, M2).

    Notes on numerical stability
    ----------------------------
    - Uses composition clipping and delegated EOS derivative guards.
    """
    x1 = _clip_x(x1)  # liquid-or-vapor mole fraction of R-1234ze(E), clamped off {0,1} to keep ln(x1) finite
    x2 = 1.0 - x1  # mole fraction of R-227ea

    mw1 = mw_from_json(d1)  # molecular weight of fluid 1 [kg/mol], read from its IDAES-Helmholtz-JSON
    mw2 = mw_from_json(d2)  # molecular weight of fluid 2 [kg/mol]
    tc1 = float(d1["basic"]["Tc"])  # pure-fluid critical temperature of fluid 1 [K]
    tc2 = float(d2["basic"]["Tc"])  # pure-fluid critical temperature of fluid 2 [K]
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1  # pure-fluid critical molar density of fluid 1 [mol/m^3] (rhoc is stored as mass density, convert via MW)
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2  # pure-fluid critical molar density of fluid 2 [mol/m^3]
    vc1 = 1.0 / rhoc1_mol  # pure-fluid critical molar volume of fluid 1 [m^3/mol]
    vc2 = 1.0 / rhoc2_mol  # pure-fluid critical molar volume of fluid 2 [m^3/mol]

    # (M1): Bell (2023) Eq. 3-5 binary reducing functions -- how the mixture's
    # own "critical-like" reducing temperature/volume depend on composition.
    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k  # (M1) reduced inverse temperature = T_red(x1) / T [unitless]
    delta = rho_mol * vred  # (M1) reduced density = rho * v_red(x1) [unitless]
    rho_red_mol = 1.0 / vred  # reducing molar density [mol/m^3], = 1/v_red; used elsewhere to map pure-fluid tau_i/delta_i

    # (M2): non-dimensional Helmholtz energy of the mixture and its tau/delta
    # derivatives, split into ideal-gas (a0) and residual (ar) parts. ar
    # already includes both the corresponding-states term AND the fitted
    # Bell (2023) departure function (Eq. 6-7) -- this single call is the
    # one place both pieces of physics enter the whole file.
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
    # a0    = alpha0_CS(tau,delta,x1)          ideal-gas alpha, corresponding-states part only
    # a0_tau= d(alpha0_CS)/d(tau)  at fixed delta,x
    # ar    = alphar_mix(tau,delta,x1)         residual alpha = alphar_CS + alphar_dep, (M2)
    # ar_tau= d(alphar_mix)/d(tau) at fixed delta,x
    # ar_del= d(alphar_mix)/d(delta) at fixed tau,x  -- this is the quantity that enters Z in (M3)

    # (M2a): add the ideal-gas entropy-of-mixing term explicitly. This term is
    # NOT part of mixture_alpha0_alphar_derivs()'s output above -- it has to be
    # added here because it depends only on composition (x1*ln(x1)+x2*ln(x2)),
    # not on tau or delta, and is delta-independent (so it never affects
    # pressure via (M3), only h/g/mu through alpha_mix itself).
    a0_mix = a0 + x1 * np.log(x1) + x2 * np.log(x2)
    # Returns: total alpha_mix, total d(alpha_mix)/d(tau), and ar_del (needed
    # standalone for the Z/pressure identity in mix_state's (M3) below).
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
    - [LT99] -- Z = 1 + delta*alphar_delta (M3) and the h/(RT), g/(RT)
      identities (M4) are standard multiparameter-Helmholtz-EOS forms from
      this framework.

    Notes on numerical stability
    ----------------------------
    - Accuracy depends on reduced derivative quality from delegated functions.
    """
    x1 = _clip_x(x1)  # phase composition this state is evaluated at (liquid x1 for bubble's liquid side, y1 for its vapor side, etc.)
    x2 = 1.0 - x1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = x1 * mw1 + x2 * mw2  # mole-fraction-weighted mixture molecular weight [kg/mol], for converting rho_mol -> rho_mass below

    # alpha = total non-dim Helmholtz energy alpha_mix (M2); alpha_tau = its
    # tau-derivative (needed for h in (M4)); ar_del = d(alphar_mix)/d(delta)
    # (needed for both Z in (M3) and h,g in (M4)).
    alpha, alpha_tau, ar_del = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1)
    # Recompute tau/delta independently here (rather than threading them
    # through from _mix_alpha_and_derivs) purely so this function's P/h/g
    # identities below are self-contained and easy to read against (M3)/(M4).
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)  # (M1)
    tau = tred / t_k  # (M1)
    delta = rho_mol * vred  # (M1)

    # (M3): compressibility factor. On the deep-liquid branch this is a
    # near-total cancellation between 1.0 and (delta*ar_del) ~ -0.8 --
    # Z ends up small (~0.2 at R-515B's real liquid state), which is why P
    # is so sensitive to any error in ar_del there. This is fine here because
    # rho_mol at this point in the code is always something the solver
    # produced itself (self-consistent), never a real/externally-measured rho.
    z = 1.0 + delta * ar_del
    p_pa = rho_mol * R_u * t_k * z  # (M3): P = rho*R*T*Z
    h_jmol = R_u * t_k * (1.0 + tau * alpha_tau + delta * ar_del)  # (M4): h/(R T) = 1 + tau*alpha_tau + delta*ar_del
    g_jmol = R_u * t_k * (1.0 + alpha + delta * ar_del)  # (M4): g/(R T) = 1 + alpha_mix + delta*ar_del
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
    - [LT99] -- A = n*R_u*T*alpha_mix is the direct extensive-Helmholtz-energy
      identity for the corresponding-states-plus-departure alpha_mix (M2)
      this framework defines; not a separate citation of its own.

    Notes on numerical stability
    ----------------------------
    - Small composition floors avoid divide-by-zero in x1 = n1/(n1+n2).
    """
    n1 = max(n1_mol, EPS_X)  # moles of fluid 1, floored away from 0 so x1=n1/(n1+n2) never divides by zero
    n2 = max(n2_mol, EPS_X)  # moles of fluid 2
    n = n1 + n2  # total moles [mol]
    x1 = _clip_x(n1 / n)  # overall mole fraction of fluid 1
    rho_mol = n / v_m3  # total molar density = n/V [mol/m^3]
    alpha, _, _ = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1)  # alpha_mix, (M2)
    return float(n * R_u * t_k * alpha)  # A = n*R_u*T*alpha_mix -- total (extensive) Helmholtz free energy [J]


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
    - [B23] -- analytic x1-derivatives of the binary reducing functions
      T_red(x1), v_red(x1) (M1, Eq. 3-5).

    Notes on numerical stability
    ----------------------------
    - Closed-form derivatives avoid finite-difference noise in composition terms.
    """
    # Analytic d/dx1 of Bell (2023) Eq. 3-5:
    #   T_red(x1) = x1^2*Tc1 + x2^2*Tc2 + 2*x1*x2*tc12*theta_T(x1)
    #   v_red(x1) = x1^2*vc1 + x2^2*vc2 + 2*x1*x2*vc12*theta_v(x1)
    # where theta_T(x1) = (x1+x2)/(beta_T^2*x1+x2) = 1/(beta_T^2*x1+x2)
    # (the (x1+x2) numerator collapses to 1 since x1+x2=1 always), and
    # tc12 = beta_T*gamma_T*sqrt(Tc1*Tc2), vc12 = beta_v*gamma_v*(vc1^(1/3)+vc2^(1/3))^3/8
    # are the "cross" reducing parameters from Table 2 of the paper.
    # Differentiating these closed forms in x1 (rather than using finite
    # differences) avoids introducing FD noise into the chemical-potential
    # calculation that consumes this derivative.
    x2 = 1.0 - x1
    p = BELL_2023_R1234ZE_R227EA  # Table 2 interaction parameters: beta_T, gamma_T, beta_v, gamma_v

    dt = (p.beta_T ** 2) * x1 + x2  # denominator of theta_T = beta_T^2*x1 + x2
    dv = (p.beta_v ** 2) * x1 + x2  # denominator of theta_v = beta_v^2*x1 + x2
    theta_t = 1.0 / dt  # theta_T(x1), the temperature "combining" weight from Eq. 3
    theta_v = 1.0 / dv  # theta_v(x1), the volume "combining" weight from Eq. 5
    # d(theta_T)/dx1: quotient rule on 1/dt, with d(dt)/dx1 = beta_T^2 - 1
    # (since dx2/dx1 = -1, so d(dt)/dx1 = beta_T^2 - 1).
    dtheta_t = -((p.beta_T ** 2) - 1.0) / (dt * dt)
    dtheta_v = -((p.beta_v ** 2) - 1.0) / (dv * dv)

    tc12 = p.beta_T * p.gamma_T * np.sqrt(tc1 * tc2)  # cross reducing-temperature parameter [K]
    vc12 = p.beta_v * p.gamma_v * ((vc1 ** (1.0 / 3.0) + vc2 ** (1.0 / 3.0)) ** 3) / 8.0  # cross reducing-volume parameter [m^3/mol]; the /8 divisor is part of Eq. 5 itself
    # d/dx1 of [x1*x2*theta(x1)] via the product rule: d(x1*x2)/dx1 = x2 - x1
    # (since x2 = 1-x1), plus x1*x2 * d(theta)/dx1.
    dxx_theta_t = (x2 - x1) * theta_t + x1 * x2 * dtheta_t
    dxx_theta_v = (x2 - x1) * theta_v + x1 * x2 * dtheta_v
    # d(T_red)/dx1 = d(x1^2)/dx1 * Tc1 + d(x2^2)/dx1 * Tc2 + 2*tc12*d(x1*x2*theta_T)/dx1
    #              = 2*x1*Tc1 - 2*x2*Tc2 + 2*tc12*dxx_theta_t
    dtred_dx1 = 2.0 * x1 * tc1 - 2.0 * x2 * tc2 + 2.0 * tc12 * dxx_theta_t
    dvred_dx1 = 2.0 * x1 * vc1 - 2.0 * x2 * vc2 + 2.0 * vc12 * dxx_theta_v
    return float(dtred_dx1), float(dvred_dx1)  # [K], [m^3/mol] -- consumed by chemical_potentials_analytic's (M5) composition-derivative chain rule


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
    - [LT99] -- mixture fugacity/chemical-potential form f_i, mu_i (M5),
      built from the composition (mole-number) derivative of n*alphar.
    - [B23] -- the alphar_mix, alphar_del, and composition-derivative terms
      (dar_dx1) this M5 formula is evaluated on.

    Notes on numerical stability
    ----------------------------
    - Uses lower bounds on rho and fugacity log arguments.
    """
    x1 = _clip_x(x1)  # composition of THIS phase (liquid x1 or vapor y1, depending on caller)
    x2 = 1.0 - x1
    rho_mol = max(float(rho_mol), RHO_MIN_MOLM3)  # guard against rho=0, which would make ln(f_i) diverge below

    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)  # (M1)
    tau = tred / t_k  # (M1)
    delta = rho_mol * vred  # (M1)
    rho_red_mol = 1.0 / vred  # reducing molar density [mol/m^3]; used to convert mixture (tau,delta) into each pure fluid's own (tau_i,delta_i) below

    # ar_mix = alphar_mix(tau,delta,x1) (M2); ar_del_mix = d(alphar_mix)/d(delta)
    # at fixed tau,x -- needed for the rho*(d alphar/d rho) term of (M5) via
    # the chain rule d(alphar)/d(rho) = d(alphar)/d(delta) * d(delta)/d(rho) = ar_del_mix * vred.
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

    # --- Composition derivative d(alphar_mix)/d(x1), needed for the x2*(d
    # alphar/d x1) term of (M5). alphar_mix = alphar_CS + alphar_dep, so this
    # splits into two pieces computed separately below: (ar1-ar2) from the
    # corresponding-states part, and ddep_dx1 from the departure function.

    # Map the MIXTURE's (tau,delta) into each PURE fluid's own reduced
    # variables via the standard corresponding-states chain rule:
    #   tau_i   = Tc_i/T = (Tc_i/T_red)*(T_red/T) = c_i * tau
    #   delta_i = rho/rhoc_i = (rho_red/rhoc_i)*(rho*v_red) = k_i * delta
    c1 = tc1 / tred  # Tc1/T_red
    c2 = tc2 / tred  # Tc2/T_red
    k1 = rho_red_mol / rhoc1_mol  # rho_red/rhoc1
    k2 = rho_red_mol / rhoc2_mol  # rho_red/rhoc2
    tau1 = c1 * tau  # pure fluid 1's own reduced inverse temperature = Tc1/T
    tau2 = c2 * tau  # pure fluid 2's own reduced inverse temperature = Tc2/T
    delta1 = k1 * delta  # pure fluid 1's own reduced density = rho/rhoc1
    delta2 = k2 * delta  # pure fluid 2's own reduced density = rho/rhoc2
    ar1, _, _ = alphar_idaes_with_derivs(d1["eos"], tau1, delta1)  # alphar of PURE fluid 1 evaluated at (tau1,delta1)
    ar2, _, _ = alphar_idaes_with_derivs(d2["eos"], tau2, delta2)  # alphar of PURE fluid 2 evaluated at (tau2,delta2)

    # d(T_red)/dx1, d(v_red)/dx1 -- analytic derivatives of the Eq. 3/5
    # reducing functions (see _bell2023_reducing_derivs_binary_local above).
    dtred_dx1, dvred_dx1 = _bell2023_reducing_derivs_binary_local(x1, tc1, tc2, vc1, vc2)
    dtau_dx1 = dtred_dx1 / t_k  # d(tau)/d(x1) at fixed T, via tau = T_red(x1)/T
    ddelta_dx1 = rho_mol * dvred_dx1  # d(delta)/d(x1) at fixed rho, via delta = rho*v_red(x1)
    dep_base = bell2023_departure_base(tau, delta)  # the departure sum evaluated WITHOUT its x1*x2 prefactor (Eq. 6), needed since d(x1*x2)/dx1 = x2-x1
    _, dep_tau, dep_del = bell2023_departure_alphar(x1, x2, tau, delta)  # dep_tau, dep_del = d(alphar_dep)/d(tau), d(alphar_dep)/d(delta) at fixed x1
    # d(alphar_dep)/dx1, by the product+chain rule on
    # alphar_dep = x1*x2*dep_base(tau(x1),delta(x1)):
    #   d/dx1 = d(x1*x2)/dx1 * dep_base + x1*x2*[dep_tau*dtau_dx1 + dep_del*ddelta_dx1]
    #         = (x2-x1)*dep_base       + dep_tau*dtau_dx1 + dep_del*ddelta_dx1
    # (the x1*x2 factor is already folded into dep_tau/dep_del as returned by
    # bell2023_departure_alphar, so it is NOT reapplied to the second term here)
    ddep_dx1 = (x2 - x1) * dep_base + dep_tau * dtau_dx1 + dep_del * ddelta_dx1
    # Total d(alphar_mix)/dx1 at fixed tau,delta: pure-fluid CS part (ar1-ar2,
    # from alphar_CS = x1*ar1(tau1,delta1) + x2*ar2(tau2,delta2), differentiated
    # holding tau1,delta1,tau2,delta2 fixed -- i.e. only the explicit x1,x2
    # prefactors vary) plus the departure-function part.
    dar_dx1 = (ar1 - ar2) + ddep_dx1
    dar_drho = ar_del_mix * vred  # (M5): d(alphar_mix)/d(rho) = d(alphar_mix)/d(delta) * d(delta)/d(rho), chain rule since delta=rho*vred

    # (M5): composition (mole-number) derivatives of n*alphar at fixed T,V.
    # d_na_dn1 = alphar_mix + rho*(d alphar/d rho) + x2*(d alphar/d x1)
    # d_na_dn2 = alphar_mix + rho*(d alphar/d rho) - x1*(d alphar/d x1)
    d_na_dn1 = ar_mix + rho_mol * dar_drho + x2 * dar_dx1
    d_na_dn2 = ar_mix + rho_mol * dar_drho - x1 * dar_dx1
    # (M5): fugacity f_i = x_i*rho*R*T*exp[d(n alphar)/dn_i]  [Pa]
    f1_pa = x1 * rho_mol * R_u * t_k * np.exp(d_na_dn1)
    f2_pa = x2 * rho_mol * R_u * t_k * np.exp(d_na_dn2)
    # (M5): chemical potential mu_i = R*T*ln(f_i)  [J/mol]. The max(1e-300,.)
    # floor guards ln() against a numerically-zero fugacity at extreme states.
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
    - Not an external literature citation -- the logistic function is a
      standard mathematical form (sigma(z) = 1/(1+e^-z)), used here purely
      as the bounded-reparameterization trick paired with the [BCL99]
      solver below (see _scaled_sigmoid/_rho_param_to_states).

    Notes on numerical stability
    ----------------------------
    - Uses branch form to reduce overflow for large |z|.
    """
    # Standard logistic sigmoid, sigmoid(z) = 1/(1+exp(-z)), maps any real z
    # into the open interval (0,1). The branch on sign(z) is purely for
    # numerical stability: exp(-z) for z>>0, or exp(z) for z<<0, keeps the
    # exponent argument from overflowing either way.
    if z >= 0:
        ez = np.exp(-z)
        return float(1.0 / (1.0 + ez))
    ez = np.exp(z)
    return float(ez / (1.0 + ez))


def _logit_from_unit_interval(y: float) -> float:
    """
    Purpose
    -------
    Map y in (0,1) to unconstrained real via logit.

    Inputs
    ------
    y : float [unitless]

    Outputs
    -------
    z : float [unitless]

    Assumptions
    -----------
    - y represents a bounded normalized variable.

    Failure modes
    -------------
    - None; y is clipped into open interval to avoid log singularities.

    References
    ----------
    - Not an external literature citation -- exact algebraic inverse of the
      logistic sigmoid used in _sigmoid above (logit(y) = ln(y/(1-y))).

    Notes on numerical stability
    ----------------------------
    - Clips away from exactly 0/1 to avoid inf values.
    """
    # Inverse of _sigmoid: logit(y) = ln(y/(1-y)). Used only when SEEDING the
    # solver -- given a physical initial guess in (0,1), find the
    # unconstrained u that _sigmoid maps back to it.
    yc = min(max(float(y), 1.0e-12), 1.0 - 1.0e-12)  # clip away from exactly 0 or 1 to avoid log(0)/div-by-0
    return float(np.log(yc / (1.0 - yc)))


def _scaled_sigmoid(z: float, lo: float, hi: float) -> float:
    """
    Purpose
    -------
    Map unconstrained z to bounded interval [lo, hi] using sigmoid scaling.

    Inputs
    ------
    z : float [unitless]
    lo : float [physical units]
    hi : float [physical units]

    Outputs
    -------
    x : float [physical units]

    Assumptions
    -----------
    - hi > lo.

    Failure modes
    -------------
    - Raises ValueError if interval is invalid.

    References
    ----------
    - [BCL99] -- this scaled-sigmoid is the specific reparameterization used
      to keep every trial point inside the box bounds that paper's
      Trust-Region-Reflective ("trf") method expects (see LSQ_METHOD).

    Notes on numerical stability
    ----------------------------
    - Stable for large |z| due sigmoid implementation in _sigmoid.
    """
    # x = lo + (hi-lo)*sigmoid(z): rescales the (0,1) output of _sigmoid onto
    # [lo,hi]. This is how the solver's unconstrained u0,u1 (in
    # _rho_param_to_states below) are turned into physically bounded
    # densities without the optimizer ever needing an explicit constraint.
    if not (hi > lo):
        raise ValueError("scaled-sigmoid requires hi > lo")
    return float(lo + (hi - lo) * _sigmoid(z))


def _inverse_scaled_sigmoid(x: float, lo: float, hi: float) -> float:
    """
    Purpose
    -------
    Map bounded physical x in [lo, hi] to unconstrained z for solver seeding.

    Inputs
    ------
    x : float [physical units]
    lo : float [physical units]
    hi : float [physical units]

    Outputs
    -------
    z : float [unitless]

    Assumptions
    -----------
    - hi > lo.

    Failure modes
    -------------
    - Raises ValueError if interval is invalid.

    References
    ----------
    - Not an external literature citation -- exact algebraic inverse of
      _scaled_sigmoid above, used by [BCL99]'s "trf" solver's initial-guess
      encoding (see LSQ_METHOD).

    Notes on numerical stability
    ----------------------------
    - Clips normalized coordinate to avoid infinities at interval bounds.
    """
    # Inverse of _scaled_sigmoid: normalize x into (0,1) via (x-lo)/(hi-lo),
    # then apply the logit. Used only for turning a physical seed guess
    # (e.g. rho_v_seed in mol/m^3) into the unconstrained u fed to least_squares.
    if not (hi > lo):
        raise ValueError("inverse scaled-sigmoid requires hi > lo")
    x_clip = min(max(float(x), lo), hi)
    y = (x_clip - lo) / (hi - lo)
    return _logit_from_unit_interval(y)


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
    - Not an external literature citation -- internal safeguard that
      replaces non-finite residual entries (NaN/Inf from a bad trial state)
      with a large finite penalty, so [BCL99]'s "trf" solver always sees a
      well-posed residual vector and can reject that trial point normally.

    Notes on numerical stability
    ----------------------------
    - Prevents NaN/Inf from breaking nonlinear iterations.
    """
    # If any of (r1,r2,r3) came out NaN/Inf (e.g. an overflowed exp() in
    # chemical_potentials_analytic at an extreme trial state), replace the
    # WHOLE vector with a large finite penalty so least_squares treats this
    # as "very wrong" and steers away from it, instead of crashing on NaN.
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
    - u0 controls rho_v via scaled-sigmoid transform.
    - u1 controls rho_l-rho_v via scaled-sigmoid transform.

    Failure modes
    -------------
    - None; values are clipped into bounded intervals.

    References
    ----------
    - Not an external literature citation -- combines _scaled_sigmoid (rho_v
      bound) and the rho_l=rho_v+increment trick above to guarantee
      rho_l > rho_v > 0 for every (u0,u1) the [BCL99] solver proposes.

    Notes on numerical stability
    ----------------------------
    - Enforces positivity and strict ordering rho_l > rho_v.
    """
    # u0 controls rho_v directly (bounded to [RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3]).
    rho_v = _scaled_sigmoid(float(u0), RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3)
    # u1 controls the GAP (rho_l - rho_v), not rho_l itself, and that gap is
    # bounded to be strictly positive (DRHO_MAP_MIN_MOLM3 > 0). This is what
    # guarantees rho_l > rho_v for every possible (u0,u1) the solver tries --
    # there is no way to express "liquid denser than vapor" being violated.
    drho = _scaled_sigmoid(float(u1), DRHO_MAP_MIN_MOLM3, DRHO_MAP_MAX_MOLM3)
    rho_l = min(max(rho_v + drho, rho_v * (1.0 + 1.0e-8)), RHO_MAX_MOLM3)  # rho_l = rho_v + drho, with a tiny floor above rho_v and a hard ceiling
    rho_v = min(rho_v, rho_l * (1.0 - 1.0e-8))  # re-clip rho_v slightly below rho_l in case the RHO_MAX_MOLM3 ceiling above changed rho_l
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
    - [LT99] -- the P/mu equilibrium conditions (M6) this residual system
      solves for the bubble point.
    - [BCL99] -- the "trf" least_squares algorithm this solves them with.

    Notes on numerical stability
    ----------------------------
    - Uses bounded least-squares and ordered density mapping.
    """
    z1 = _clip_x(z1)  # FIXED liquid-phase (feed) composition for this bubble solve -- z1 is never adjusted by the solver
    z2 = 1.0 - z1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    # Mixture reducing temperature at this fixed feed composition -- used ONLY
    # to taper the rho-separation acceptance gate near critical (see
    # _rho_separation_min_ratio above); not used anywhere else in this solve.
    tred_mix, _vred_mix = bell2023_Tred_vred(z1, z2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)

    def res(u):
        """(M6) residual vector [r1,r2,r3] for the unknowns u=(u0,u1,u2),
        where u0,u1 encode (rho_l,rho_v) and u2 encodes y1 (vapor comp)."""
        rho_l, rho_v = _rho_param_to_states(float(u[0]), float(u[1]))  # decode this trial's liquid/vapor density [mol/m^3]
        y1 = _clip_x(_sigmoid(float(u[2])))  # decode this trial's vapor mole fraction
        try:
            st_l = mix_state(d1, d2, t_k, rho_l, z1)  # liquid-phase state: P_l, h_l, ... at (T, rho_l, z1)
            st_v = mix_state(d1, d2, t_k, rho_v, y1)  # vapor-phase state: P_v, h_v, ... at (T, rho_v, y1)
            mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, z1)  # liquid-phase mu1, mu2 (M5)
            mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, y1)  # vapor-phase mu1, mu2 (M5)
            # (M6) r1: mechanical equilibrium, P_l == P_v, normalized by the
            # average pressure so the residual is a relative (not absolute) error.
            r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
            # (M6) r2, r3: phase equilibrium for each component, mu_i^liq ==
            # mu_i^vap, normalized by R*T so all three residuals are dimensionless
            # and on comparable scale for the optimizer.
            r2 = (mu1_l - mu1_v) / (R_u * t_k)
            r3 = (mu2_l - mu2_v) / (R_u * t_k)
            return _safe_residual_vector(np.array([r1, r2, r3], dtype=float))
        except Exception:
            # Any exception (e.g. overflow deep in the EOS at a wild trial
            # state) is treated the same as a non-finite residual: a large
            # penalty steers the solver away without crashing the whole run.
            return np.array([1.0e6, 1.0e6, 1.0e6], dtype=float)

    def _attempt(rho_l_seed_in: float, rho_v_seed_in: float, y1_seed_in: float) -> Dict:
        """One full least_squares solve of (M6) from one specific starting
        guess (rho_l_seed_in, rho_v_seed_in, y1_seed_in). Called multiple
        times with different seeds by the retry loop below."""
        rho_v_seed = min(max(float(rho_v_seed_in), RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.9)  # clamp seed into the valid rho_v range
        rho_l_seed = min(max(float(rho_l_seed_in), rho_v_seed * (1.0 + 1.0e-6)), RHO_MAX_MOLM3)  # clamp seed liquid density, keeping it above rho_v_seed
        dr_seed = max(rho_l_seed - rho_v_seed, 1.0e-6)  # seed density gap, must stay strictly positive
        y1_seed = _clip_x(y1_seed_in)
        # Convert the physical seed (rho_l_seed, rho_v_seed, y1_seed) into the
        # solver's unconstrained u-space via the inverse sigmoid maps, so
        # least_squares starts its search already at the physically-motivated guess.
        u0 = np.array(
            [
                _inverse_scaled_sigmoid(rho_v_seed, RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3),
                _inverse_scaled_sigmoid(dr_seed, DRHO_MAP_MIN_MOLM3, DRHO_MAP_MAX_MOLM3),
                np.log(y1_seed / (1.0 - y1_seed)),  # inline logit(y1_seed); same formula as _logit_from_unit_interval
            ],
            dtype=float,
        )
        lb = np.array([-30.0, -30.0, -30.0], dtype=float)  # generous bounds on the UNCONSTRAINED u-space itself (the physical bounds come from the sigmoid maps, not these)
        ub = np.array([30.0, 30.0, 30.0], dtype=float)
        # Solve (M6): drive res(u) -> [0,0,0] starting from u0.
        sol = least_squares(
            res,
            u0,
            bounds=(lb, ub),
            method=LSQ_METHOD,
            ftol=LSQ_FTOL,
            xtol=LSQ_XTOL,
            gtol=LSQ_GTOL,
            max_nfev=LSQ_MAX_NFEV,
            x_scale=LSQ_X_SCALE,
            diff_step=LSQ_DIFF_STEP,
            loss=LSQ_LOSS,
            f_scale=LSQ_F_SCALE,
        )
        rho_l, rho_v = _rho_param_to_states(float(sol.x[0]), float(sol.x[1]))  # decode the CONVERGED (or best-effort) solution's densities
        y1 = _clip_x(_sigmoid(float(sol.x[2])))  # decode the converged vapor composition
        # Re-evaluate everything at the final (rho_l, rho_v, y1) to get the
        # reported P/h and to independently verify how well (M6) is actually
        # satisfied (r_p, r_mu below), rather than trusting sol.success alone.
        st_l = mix_state(d1, d2, t_k, rho_l, z1)
        st_v = mix_state(d1, d2, t_k, rho_v, y1)
        mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, z1)
        mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, y1)
        r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))  # final mechanical-equilibrium residual (M6, r1), absolute value
        r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)  # worst of the two chemical-potential residuals (M6, r2/r3)
        # Strict convergence gate: scipy reports success AND both physical
        # residuals are below 1e-6 AND the liquid density is genuinely above
        # the vapor density (guards against a degenerate rho_l==rho_v root).
        # The required ratio itself is TAPERED near the mixture critical
        # temperature (see _rho_separation_min_ratio) so true near-critical
        # solutions -- where rho_l and rho_v are physically close -- are not
        # rejected alongside fake degenerate ones.
        min_ratio = _rho_separation_min_ratio(t_k, tred_mix)
        ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * min_ratio)
        return {
            "status": "CONVERGED" if ok else "DIVERGED",
            "T_K": float(t_k),
            "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),  # (M6): dome pressure at this T = the solved P_l~=P_v, NEVER computed by feeding in an assumed rho
            "rho_l_molm3": float(rho_l),  # solved liquid molar density [mol/m^3] -- this is the quantity that inherits the ~2% Bell-2023 extrapolation bias at R-515B's real composition
            "rho_v_molm3": float(rho_v),  # solved vapor molar density [mol/m^3]
            "x1_liq": float(z1),  # liquid composition (fixed input for the bubble solve)
            "y1_vap": float(y1),  # solved equilibrium vapor composition
            "h_l_Jmol": float(st_l.h_jmol),  # saturated-liquid molar enthalpy [J/mol], via (M4) at (T,rho_l,z1)
            "h_v_Jmol": float(st_v.h_jmol),  # saturated-vapor molar enthalpy [J/mol], via (M4) at (T,rho_v,y1)
            "r_P": float(r_p),  # final mechanical-equilibrium residual, should be <=1e-6 for a trustworthy point
            "r_mu": float(r_mu),  # final chemical-potential residual, should be <=1e-6 for a trustworthy point
            "iterations": int(sol.nfev),
            "notes": str(sol.message),
            "rho_sep_min_ratio": float(min_ratio),  # required ratio actually used at this T (3.0 far from critical, tapering to 1.05 near Tred_mix) -- lets downstream plotting flag near-critical points distinctly
            "near_critical": bool(min_ratio < RHO_SEP_RATIO_FAR),  # True once this point is inside the taper window; use to overlay/mark near-critical dome points as lower-confidence rather than plotting them identically to the rest
        }

    y10 = _clip_x(y10)
    # Four different starting guesses, perturbing both densities and the
    # vapor-composition seed, tried in order until one converges. This exists
    # because a single fixed seed can converge to a spurious/nonphysical root
    # at some (T,composition) combinations even while scipy reports success
    # (see the 2026-03-03 breadcrumbs on solver seed-robustness) -- retrying
    # from nearby seeds is cheap insurance against that failure mode. NOTE:
    # when called from run_true_vle_envelope's temperature sweep below, the
    # FIRST seed here is already the previous temperature's converged
    # solution (continuation), so these four attempts are a second layer of
    # robustness on top of that, not a replacement for it.
    attempts = [
        (rho_l0, rho_v0, y10),
        (1.1 * rho_l0, 0.7 * rho_v0, y10),
        (0.9 * rho_l0, 0.5 * rho_v0, _clip_x(y10 + 0.01)),
        (1.2 * rho_l0, 0.5 * rho_v0, _clip_x(y10 - 0.01)),
    ]
    best = None
    best_score = np.inf
    for idx, (rl, rv, yy) in enumerate(attempts):
        row = _attempt(rl, rv, yy)
        if row["status"] == "CONVERGED":
            row["notes"] = f"{row['notes']} | retry={idx}"
            return row
        # None converged yet: keep the least-bad attempt (smallest combined
        # residual) in case every seed ultimately fails, so the caller still
        # gets a "best effort" DIVERGED row with diagnostic residuals rather
        # than nothing.
        score = float(row["r_P"] + row["r_mu"])
        if score < best_score:
            best = row
            best_score = score
    assert best is not None
    best["notes"] = f"{best['notes']} | retry=best_failed"
    return best


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
    - [LT99] -- the P/mu equilibrium conditions (M6) this residual system
      solves for the dew point (mirrors solve_bubble_at_t with x1 free).
    - [BCL99] -- the "trf" least_squares algorithm this solves them with.

    Notes on numerical stability
    ----------------------------
    - Uses bounded least-squares and ordered density mapping.
    """
    z1 = _clip_x(z1)  # FIXED vapor-phase (feed) composition for this dew solve -- mirror image of solve_bubble_at_t's fixed z1 (there it's the liquid comp)
    z2 = 1.0 - z1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    # Mixture reducing temperature at this fixed feed composition -- used ONLY
    # to taper the rho-separation acceptance gate near critical (mirrors
    # solve_bubble_at_t; see _rho_separation_min_ratio above).
    tred_mix, _vred_mix = bell2023_Tred_vred(z1, z2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)

    def res(u):
        """Same (M6) residual system as solve_bubble_at_t.res(), but with the
        roles of liquid/vapor composition swapped: here the VAPOR composition
        z1 is fixed and the LIQUID composition x1 is the third unknown
        (decoded from u[2]), since a dew point means "first drop of liquid
        forms from a saturated vapor of known (feed) composition z1"."""
        rho_l, rho_v = _rho_param_to_states(float(u[0]), float(u[1]))
        x1 = _clip_x(_sigmoid(float(u[2])))  # solved liquid composition (unknown here, vs. y1 being unknown in the bubble solve)
        try:
            st_l = mix_state(d1, d2, t_k, rho_l, x1)  # liquid state at the SOLVED x1
            st_v = mix_state(d1, d2, t_k, rho_v, z1)  # vapor state at the FIXED feed composition z1
            mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, x1)
            mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, z1)
            r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))  # (M6) r1: mechanical equilibrium
            r2 = (mu1_l - mu1_v) / (R_u * t_k)  # (M6) r2: component-1 phase equilibrium
            r3 = (mu2_l - mu2_v) / (R_u * t_k)  # (M6) r3: component-2 phase equilibrium
            return _safe_residual_vector(np.array([r1, r2, r3], dtype=float))
        except Exception:
            return np.array([1.0e6, 1.0e6, 1.0e6], dtype=float)

    def _attempt(rho_l_seed_in: float, rho_v_seed_in: float, x1_seed_in: float) -> Dict:
        """One full least_squares solve of (M6) from one starting guess --
        identical structure to solve_bubble_at_t's _attempt(), just solving
        for x1 instead of y1."""
        rho_v_seed = min(max(float(rho_v_seed_in), RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.9)
        rho_l_seed = min(max(float(rho_l_seed_in), rho_v_seed * (1.0 + 1.0e-6)), RHO_MAX_MOLM3)
        dr_seed = max(rho_l_seed - rho_v_seed, 1.0e-6)
        x1_seed = _clip_x(x1_seed_in)
        u0 = np.array(
            [
                _inverse_scaled_sigmoid(rho_v_seed, RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3),
                _inverse_scaled_sigmoid(dr_seed, DRHO_MAP_MIN_MOLM3, DRHO_MAP_MAX_MOLM3),
                np.log(x1_seed / (1.0 - x1_seed)),
            ],
            dtype=float,
        )
        lb = np.array([-30.0, -30.0, -30.0], dtype=float)
        ub = np.array([30.0, 30.0, 30.0], dtype=float)
        sol = least_squares(
            res,
            u0,
            bounds=(lb, ub),
            method=LSQ_METHOD_DEW,
            ftol=LSQ_FTOL_DEW,
            xtol=LSQ_XTOL_DEW,
            gtol=LSQ_GTOL_DEW,
            max_nfev=LSQ_MAX_NFEV_DEW,
            x_scale=LSQ_X_SCALE_DEW,
            diff_step=LSQ_DIFF_STEP_DEW,
            loss=LSQ_LOSS_DEW,
            f_scale=LSQ_F_SCALE_DEW,
        )
        rho_l, rho_v = _rho_param_to_states(float(sol.x[0]), float(sol.x[1]))
        x1 = _clip_x(_sigmoid(float(sol.x[2])))
        st_l = mix_state(d1, d2, t_k, rho_l, x1)
        st_v = mix_state(d1, d2, t_k, rho_v, z1)
        mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, x1)
        mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, z1)
        r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
        r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)
        # Same tapered gate as solve_bubble_at_t -- see _rho_separation_min_ratio.
        min_ratio = _rho_separation_min_ratio(t_k, tred_mix)
        ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * min_ratio)
        return {
            "status": "CONVERGED" if ok else "DIVERGED",
            "T_K": float(t_k),
            "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),  # dome pressure at this T from the dew branch; compare to solve_bubble_at_t's P at the same T as an azeotrope consistency check
            "rho_l_molm3": float(rho_l),  # solved (incipient) liquid density [mol/m^3]
            "rho_v_molm3": float(rho_v),  # solved saturated-vapor density [mol/m^3]
            "x1_liq": float(x1),  # solved equilibrium liquid composition (the unknown for a dew solve)
            "y1_vap": float(z1),  # vapor composition (fixed input for the dew solve)
            "h_l_Jmol": float(st_l.h_jmol),  # incipient-liquid molar enthalpy [J/mol]
            "h_v_Jmol": float(st_v.h_jmol),  # saturated-vapor molar enthalpy [J/mol]
            "r_P": float(r_p),
            "r_mu": float(r_mu),
            "iterations": int(sol.nfev),
            "notes": str(sol.message),
            "rho_sep_min_ratio": float(min_ratio),  # required ratio actually used at this T -- see solve_bubble_at_t for details
            "near_critical": bool(min_ratio < RHO_SEP_RATIO_FAR),  # True inside the taper window; use to overlay/mark near-critical dome points distinctly
        }

    x10 = _clip_x(x10)
    # Same four-seed retry strategy as solve_bubble_at_t, perturbing
    # densities and the (here: liquid) composition seed.
    attempts = [
        (rho_l0, rho_v0, x10),
        (1.1 * rho_l0, 0.7 * rho_v0, x10),
        (0.9 * rho_l0, 0.5 * rho_v0, _clip_x(x10 + 0.01)),
        (1.2 * rho_l0, 0.5 * rho_v0, _clip_x(x10 - 0.01)),
    ]
    best = None
    best_score = np.inf
    for idx, (rl, rv, xx) in enumerate(attempts):
        row = _attempt(rl, rv, xx)
        if row["status"] == "CONVERGED":
            row["notes"] = f"{row['notes']} | retry={idx}"
            return row
        score = float(row["r_P"] + row["r_mu"])
        if score < best_score:
            best = row
            best_score = score
    assert best is not None
    best["notes"] = f"{best['notes']} | retry=best_failed"
    return best


def _pressure_rho_derivatives_fd(
    d1: Dict,
    d2: Dict,
    t_k: float,
    rho_mol: float,
    x1: float,
    rel_step: float = 1.0e-4,
) -> Tuple[float, float]:
    """
    Purpose
    -------
    Estimate (dP/drho)_T and (d2P/drho2)_T at fixed (T, x1) via central
    finite differences on mix_state's own pressure evaluation, for use in
    the mixture critical-point solve (see solve_mixture_critical_point).

    Inputs
    ------
    d1, d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]
    rel_step : float [unitless]
      Relative finite-difference step size (h = rel_step * rho_mol).

    Outputs
    -------
    dp_drho : float [Pa/(mol/m^3)]
    d2p_drho2 : float [Pa/(mol/m^3)^2]

    Assumptions
    -----------
    - P(rho) is smooth enough at fixed (T,x1) for a 3-point central
      finite-difference stencil to be accurate near the critical region.
    - No analytic third-derivative-of-alphar machinery exists in this
      file (only first derivatives are used, via mix_state), so finite
      differences on P itself -- not analytic differentiation of
      alphar's own derivatives -- is the simplest correct path to
      (d^2 P/d rho^2)_T here.

    Failure modes
    -------------
    - Propagates exceptions from mix_state (e.g. EOS overflow) if
      rho+/-h lands in an unstable/invalid region.

    References
    ----------
    - Standard central finite-difference formulas; not a specific
      external citation.

    Notes on numerical stability
    ----------------------------
    - rel_step=1e-4 balances truncation error (favors larger h) against
      floating-point cancellation error (favors smaller h) for P values
      that are O(1e6-1e7) Pa in this near-critical region.
    """
    h = max(rel_step * rho_mol, 1.0e-6)  # absolute FD step [mol/m^3]; floored so it never collapses to zero at tiny rho
    p_plus = mix_state(d1, d2, t_k, rho_mol + h, x1).p_pa
    p_minus = mix_state(d1, d2, t_k, rho_mol - h, x1).p_pa
    p_mid = mix_state(d1, d2, t_k, rho_mol, x1).p_pa
    dp_drho = (p_plus - p_minus) / (2.0 * h)  # standard central first derivative
    d2p_drho2 = (p_plus - 2.0 * p_mid + p_minus) / (h * h)  # standard central second derivative
    return float(dp_drho), float(d2p_drho2)


def _dp_drho_fd(
    d1: Dict,
    d2: Dict,
    t_k: float,
    rho_mol: float,
    x1: float,
    rel_step: float = 1.0e-4,
) -> float:
    """
    Purpose
    -------
    Estimate (dP/drho)_T ONLY (first derivative) via a 2-point central
    finite difference on mix_state's pressure. Lighter-weight than
    _pressure_rho_derivatives_fd's 3-point stencil; used repeatedly by
    the inner bracketed spinodal-density root-find in
    solve_mixture_critical_point (called many times per outer step, so
    the cheaper 2-evaluation form matters).

    Inputs
    ------
    d1, d2 : dict [unitless]
    t_k : float [K]
    rho_mol : float [mol/m^3]
    x1 : float [mol/mol]
    rel_step : float [unitless]

    Outputs
    -------
    dp_drho : float [Pa/(mol/m^3)]

    Assumptions
    -----------
    - Same smoothness assumption as _pressure_rho_derivatives_fd.

    Failure modes
    -------------
    - Propagates mix_state exceptions (caller wraps in try/except).

    References
    ----------
    - Standard central finite-difference formula.

    Notes on numerical stability
    ----------------------------
    - Not applicable beyond the general FD step-size discussion in
      _pressure_rho_derivatives_fd.
    """
    h = max(rel_step * rho_mol, 1.0e-6)
    p_plus = mix_state(d1, d2, t_k, rho_mol + h, x1).p_pa
    p_minus = mix_state(d1, d2, t_k, rho_mol - h, x1).p_pa
    return float((p_plus - p_minus) / (2.0 * h))


def _find_sign_change_bracket(
    func,
    x_lo: float,
    x_hi: float,
    n_scan: int = 30,
) -> Optional[Tuple[float, float, float, float]]:
    """
    Purpose
    -------
    Scan [x_lo, x_hi] at n_scan evenly-spaced points and return the first
    bracket (a, b, f(a), f(b)) where func changes sign -- a prerequisite
    for scipy.optimize.brentq, which requires a verified sign change and
    does not search for one itself.

    Inputs
    ------
    func : callable float -> float
    x_lo, x_hi : float
    n_scan : int [count]

    Outputs
    -------
    (a, b, fa, fb) : tuple[float,float,float,float], or None if no sign
      change was found (or func raised/returned non-finite everywhere).

    Assumptions
    -----------
    - A single relevant sign change exists within the scanned window and
      is coarse enough to be caught by n_scan points; the caller chooses
      the window/resolution based on how far from the true root the
      window's endpoints are expected to be.

    Failure modes
    -------------
    - Returns None (not an exception) if func raises or is non-finite at
      every scan point, or if no sign change is found -- callers must
      check for None.

    References
    ----------
    - Standard root-bracketing technique underlying Brent's method
      (scipy.optimize.brentq); matches the "nested and bounded
      iterations of Brent's method" approach validated for exactly this
      mixture-critical-point problem by Hoteit et al., as cited in Bell
      & Jager (2017), Fluid Phase Equilibria 433, 159-173 (the paper
      this whole nested-bracketed-solve design is based on).

    Notes on numerical stability
    ----------------------------
    - Non-finite/exception points are treated as gaps (reset the
      previous-point tracker) rather than false sign changes.
    """
    xs = np.linspace(x_lo, x_hi, n_scan)
    f_prev: Optional[float] = None
    x_prev: Optional[float] = None
    for x in xs:
        try:
            f = float(func(float(x)))
        except Exception:
            f_prev, x_prev = None, None
            continue
        if not np.isfinite(f):
            f_prev, x_prev = None, None
            continue
        if f_prev is not None and f != 0.0 and np.sign(f) != np.sign(f_prev):
            return (float(x_prev), float(x), float(f_prev), float(f))
        f_prev, x_prev = f, float(x)
    return None


# Search-window and scan-resolution constants for
# solve_mixture_critical_point's nested search. Windows are deliberately
# NARROW relative to Bell & Jager (2017)'s from-scratch algorithm (which
# searches from a cold start with no prior information) because we
# already have an excellent initial guess (Tred_mix, 1/vred_mix from the
# Bell 2023 mixing rule, or better: real converged bubble/dew densities
# from the sweep -- see run_true_vle_envelope).
CRIT_RHO_SCAN_LO_FACTOR = 0.7  # inner rho scan window default (only used when no explicit rho_scan_lo/hi override is passed): [factor*rho_guess, ...]
CRIT_RHO_SCAN_HI_FACTOR = 1.6  # inner rho scan window default: [..., factor*rho_guess]
CRIT_RHO_SCAN_N = 150  # inner scan resolution (evaluations per T trial) -- raised from an initial 40 while chasing the (now-abandoned) d2P/drho2-based design; kept high here too since resolving a narrow spinodal-density gap still benefits from it
CRIT_T_SCAN_LO_K = 6.0  # outer T scan window default (only used when no explicit t_scan_lo/hi override is passed): [t_guess-this, ...]
CRIT_T_SCAN_HI_K = 3.0  # outer T scan window default: [..., t_guess+this]
# CRIT_BISECT_MAX_ITERS/CRIT_BISECT_T_TOL_K/CRIT_WIDTH_FRAC_GATE: the
# bisection-on-dip-existence design's tuning (see solve_mixture_critical_point
# docstring for why this replaced a d2P/drho2-based outer root-find).
CRIT_BISECT_MAX_ITERS = 50  # hard cap on bisection steps (T-tolerance below is reached in ~log2((t_hi-t_lo)/tol) steps well under this)
CRIT_BISECT_T_TOL_K = 1.0e-4  # bisection stops once the T bracket is narrower than this
CRIT_WIDTH_FRAC_GATE = 0.02  # converged requires the final spinodal-density gap to be below this fraction of rho_scale (a physically meaningful "the two branches have essentially merged" check)


def _find_all_sign_changes(
    func,
    x_lo: float,
    x_hi: float,
    n_scan: int = 30,
) -> List[Tuple[float, float, float, float]]:
    """
    Purpose
    -------
    Like _find_sign_change_bracket, but returns EVERY bracket found while
    scanning [x_lo, x_hi], not just the first. Used to find BOTH spinodal
    branches (vapor-side and liquid-side) from a single scan of
    (dP/drho)_T at fixed T, since the classic sub-critical isotherm shape
    has dP/drho go positive -> negative (vapor-side spinodal) -> positive
    again (liquid-side spinodal) as density increases.

    Inputs
    ------
    func : callable float -> float
    x_lo, x_hi : float
    n_scan : int [count]

    Outputs
    -------
    brackets : list[tuple[float,float,float,float]]
      Each entry is (a, b, f(a), f(b)) for one detected sign change, in
      scan order (increasing x). Empty list if none found.

    Assumptions
    -----------
    - Same as _find_sign_change_bracket, generalized to multiple crossings.

    Failure modes
    -------------
    - Returns an empty list (not an exception) if func raises or is
      non-finite at every scan point, or no sign change is found.

    References
    ----------
    - Same basis as _find_sign_change_bracket (Brent's-method bracketing;
      Bell & Jager 2017 / Hoteit et al.'s nested-bounded-Brent approach).

    Notes on numerical stability
    ----------------------------
    - Non-finite/exception points are treated as gaps (reset the
      previous-point tracker) rather than false sign changes.
    """
    xs = np.linspace(x_lo, x_hi, n_scan)
    brackets: List[Tuple[float, float, float, float]] = []
    f_prev: Optional[float] = None
    x_prev: Optional[float] = None
    for x in xs:
        try:
            f = float(func(float(x)))
        except Exception:
            f_prev, x_prev = None, None
            continue
        if not np.isfinite(f):
            f_prev, x_prev = None, None
            continue
        if f_prev is not None and f != 0.0 and np.sign(f) != np.sign(f_prev):
            brackets.append((float(x_prev), float(x), float(f_prev), float(f)))
        f_prev, x_prev = f, float(x)
    return brackets


def _spinodal_pair_at_t(
    d1: Dict,
    d2: Dict,
    t_k: float,
    x1: float,
    rho_scan_lo: float,
    rho_scan_hi: float,
    rel_step: float = 1.0e-4,
    n_scan: int = CRIT_RHO_SCAN_N,
) -> Optional[Tuple[float, float]]:
    """
    Purpose
    -------
    At FIXED T, find BOTH spinodal densities (vapor-side rho_v_sp,
    liquid-side rho_l_sp, with rho_v_sp < rho_l_sp) bounding the
    mechanically unstable ((dP/drho)_T<0) density region within
    [rho_scan_lo, rho_scan_hi], via _find_all_sign_changes + brentq
    refinement on the first and last detected brackets.

    Inputs
    ------
    d1, d2 : dict [unitless]
    t_k : float [K]
    x1 : float [mol/mol]
    rho_scan_lo, rho_scan_hi : float [mol/m^3]
    rel_step : float [unitless]
    n_scan : int [count]

    Outputs
    -------
    (rho_v_sp, rho_l_sp) : tuple[float,float] [mol/m^3], or None if fewer
      than 2 sign changes were found -- i.e. no "dip" exists in this
      window at this T, which happens when T is at/above the mixture's
      critical temperature (confirmed empirically for this
      mixture/composition via diagnose_spinodal.py: a clear two-root dip
      is present at T=381.5K and absent by T=382.0K).

    Assumptions
    -----------
    - Exactly one dip (two spinodal roots) is expected in the scanned
      window below Tc; if more than 2 sign changes are found (unlikely
      for a well-behaved EOS in a window this narrow), only the FIRST
      and LAST are used, which still correctly bounds the full unstable
      region.

    Failure modes
    -------------
    - Returns None rather than raising if fewer than 2 brackets are found.

    References
    ----------
    - Bell & Jager (2017), Fluid Phase Equilibria 433, 159-173 -- Section
      3.2's L1=0 spinodal contour, fixed-composition-only reduction.

    Notes on numerical stability
    ----------------------------
    - This function's outputs feed a BISECTION on whether it returns
      None or a real pair (see solve_mixture_critical_point) rather than
      feeding a second-derivative-based root-find. A prior design used
      (d2P/drho2)_T at a single spinodal point as an outer root-finding
      target, but diagnose_spinodal.py showed that quantity flips sign
      essentially at random near critical (values like -5.9e-3, -7.3e-3,
      +4.0e-3, -8.1e-4, +1.9e-3 across consecutive 0.5-1K steps) --
      finite-difference noise dominating a signal that is genuinely
      shrinking toward zero. The spinodal-density GAP used here only
      needs first derivatives, which were independently confirmed
      reliable (bracketed roots agreeing with brentq to ~1e-13).
    """
    def f(rho: float) -> float:
        return _dp_drho_fd(d1, d2, t_k, rho, x1, rel_step=rel_step)

    brackets = _find_all_sign_changes(f, rho_scan_lo, rho_scan_hi, n_scan=n_scan)
    if len(brackets) < 2:
        return None
    a0, b0, _fa0, _fb0 = brackets[0]
    a1, b1, _fa1, _fb1 = brackets[-1]
    rho_v_sp = float(brentq(f, a0, b0, xtol=1.0e-6, rtol=1.0e-12, maxiter=100))
    rho_l_sp = float(brentq(f, a1, b1, xtol=1.0e-6, rtol=1.0e-12, maxiter=100))
    return rho_v_sp, rho_l_sp


def solve_mixture_critical_point(
    d1: Dict,
    d2: Dict,
    z1: float,
    t_guess: float,
    rho_guess: float,
    rho_scan_lo: Optional[float] = None,
    rho_scan_hi: Optional[float] = None,
    t_scan_lo: Optional[float] = None,
    t_scan_hi: Optional[float] = None,
) -> Dict:
    """
    Purpose
    -------
    Solve for the mixture's true critical point (Tc_mix, rhoc_mix) at
    fixed feed composition z1 -- the actual point where bubble and dew
    genuinely coincide. Distinct from Tred_mix/vred_mix (the Bell 2023
    mixing-rule reducing point used elsewhere in this file for the
    near-critical acceptance-gate taper), which is only an approximation
    to this true EOS critical point, close enough to seed the search.

    THIRD design iteration for this function (see prior two attempts
    documented in PROJECT_CONTEXT.md, both of which failed on real runs):
    (1) a flat 2D least_squares solve on [(dP/drho)_T, (d2P/drho2)_T]
    stalled near critical (near-singular Jacobian, expected at a genuine
    critical point); (2) a nested bracketed-Brent redesign using
    (d2P/drho2)_T as the OUTER root-finding target found no sign-change
    bracket even with a correctly-placed window and a 4-5x resolution
    increase. Diagnosed via a standalone script (diagnose_spinodal.py)
    that printed raw (dP/drho)_T values across the search grid: the
    density gap between the two spinodal branches DOES shrink cleanly
    toward zero as T approaches the true Tc (a real, well-behaved
    physical signal), but (d2P/drho2)_T evaluated there is essentially
    finite-difference NOISE at that point -- its sign flips
    unpredictably across consecutive small T steps instead of shrinking
    monotonically, because the true curvature genuinely IS shrinking
    toward zero (a weak signal) while the FD noise floor does not, so
    noise dominates exactly where precision is needed most.

    This THIRD design avoids the noisy second derivative entirely:
    Stage 1 (_spinodal_pair_at_t) finds BOTH spinodal densities at a
    given T using only first derivatives (already independently
    confirmed reliable to ~1e-13). Stage 2 is a BISECTION on T: at each
    trial T, check only whether a two-root "dip" exists at all
    (_spinodal_pair_at_t returns a pair) or not (returns None) -- this is
    a robust, discrete, always-well-defined test (no noisy quantity to
    root-find), and the T where the dip disappears converges to Tc as
    the bisection bracket narrows.

    Inputs
    ------
    d1, d2 : dict [unitless]
    z1 : float [mol/mol]
      Fixed overall (azeotrope) composition.
    t_guess, rho_guess : float [K], [mol/m^3]
      Initial guess -- Tred_mix and 1/vred_mix from bell2023_Tred_vred,
      or (preferably) values derived from real converged bubble/dew
      densities (see run_true_vle_envelope). Used to build default
      search windows (CRIT_*_SCAN_* factors/offsets) UNLESS explicit
      window bounds are provided below.
    rho_scan_lo, rho_scan_hi : float | None [mol/m^3]
      Optional explicit override for the inner (density) search window.
      If omitted, defaults to [CRIT_RHO_SCAN_LO_FACTOR,
      CRIT_RHO_SCAN_HI_FACTOR]*rho_guess.
    t_scan_lo, t_scan_hi : float | None [K]
      Optional explicit override for the outer (temperature) bisection
      window (absolute K values, not offsets). If omitted, defaults to
      [t_guess-CRIT_T_SCAN_LO_K, t_guess+CRIT_T_SCAN_HI_K]. REQUIRES a
      dip to exist at t_scan_lo and NOT exist at t_scan_hi (checked
      explicitly before bisecting; returns converged=False with a
      specific diagnostic note if either precondition fails).

    Outputs
    -------
    result : dict [mixed units]
      T_K, rho_molm3, P_Pa, h_Jmol, x1, converged (bool), dp_drho_resid,
      d2p_drho2_resid (both now DIAGNOSTIC ONLY, not the convergence
      gate -- see above), spinodal_width_molm3, iterations, notes. This
      single physical point is used as the shared closing vertex for
      BOTH the bubble and dew branches in plot_envelope.

    Assumptions
    -----------
    - Critical point solved at FIXED composition z1 (this file's
      pseudo-pure/azeotrope treatment) -- the standard fixed-N
      (dP/dV)_T=(d2P/dV2)_T=0 criticality condition, not the full
      multicomponent Legendre-transform/matrix-determinant machinery of
      Bell & Jager (2017) (needed only when composition is a free
      variable; it is not here).
    - Bisection precondition: a genuine two-root dip exists at t_scan_lo
      and has vanished by t_scan_hi. If the true Tc is outside this
      window, both preconditions are checked explicitly and the function
      fails loudly (specific note) rather than bisecting on a false
      assumption.

    Failure modes
    -------------
    - Returns converged=False if either bisection-window precondition
      fails, or if bisection completes but the final spinodal-density
      gap (as a fraction of rho_scale) still exceeds CRIT_WIDTH_FRAC_GATE.

    References
    ----------
    - Bell, I.H., Jager, A. (2017). "Calculation of critical points from
      Helmholtz-energy-explicit mixture models." Fluid Phase Equilibria,
      433, 159-173. https://doi.org/10.1016/j.fluid.2016.10.030 --
      source of the general nested-search strategy (Section 3.2) and the
      observation that Newton-Raphson requires "quite good estimates"
      near critical (Section 3.1); also citing, within that paper, Hoteit
      et al.'s validated "nested and bounded iterations of Brent's
      method" approach for critical-point finding, which this
      bisection-on-dip-existence design is a further, noise-robust
      adaptation of.
    - [LT99] -- P(T,rho,x) Helmholtz-EOS identity (M3) this solve's
      finite-difference derivatives are taken of.

    Notes on numerical stability
    ----------------------------
    - dp_drho_resid/d2p_drho2_resid are still computed and reported for
      diagnostic continuity with the prior design, but are NOT used to
      determine `converged` -- see the Purpose section for why the
      second-derivative quantity is unreliable near critical.
    - Residuals/width are nondimensionalized by FIXED p_scale/rho_scale
      constants computed once from the seed, not recomputed per-iteration.
    """
    if rho_scan_lo is None:
        rho_scan_lo = CRIT_RHO_SCAN_LO_FACTOR * rho_guess
    if rho_scan_hi is None:
        rho_scan_hi = CRIT_RHO_SCAN_HI_FACTOR * rho_guess
    if t_scan_lo is None:
        t_scan_lo = t_guess - CRIT_T_SCAN_LO_K
    if t_scan_hi is None:
        t_scan_hi = t_guess + CRIT_T_SCAN_HI_K

    p_scale = max(abs(mix_state(d1, d2, t_guess, rho_guess, z1).p_pa), 1.0)
    rho_scale = max(rho_guess, 1.0)

    def _fallback(notes: str, iters: int) -> Dict:
        return {
            "T_K": float(t_guess),
            "rho_molm3": float(rho_guess),
            "P_Pa": float(p_scale),
            "h_Jmol": float("nan"),
            "x1": float(z1),
            "converged": False,
            "dp_drho_resid": float("nan"),
            "d2p_drho2_resid": float("nan"),
            "spinodal_width_molm3": float("nan"),
            "iterations": iters,
            "notes": notes,
        }

    pair_lo = _spinodal_pair_at_t(d1, d2, t_scan_lo, z1, rho_scan_lo, rho_scan_hi)
    if pair_lo is None:
        return _fallback(
            f"no spinodal dip found at t_scan_lo={t_scan_lo:.3f} K -- lower T bound is already at/above the mixture critical temperature, or rho window [{rho_scan_lo:.1f},{rho_scan_hi:.1f}] mol/m^3 doesn't bracket it there",
            0,
        )
    pair_hi = _spinodal_pair_at_t(d1, d2, t_scan_hi, z1, rho_scan_lo, rho_scan_hi)
    if pair_hi is not None:
        return _fallback(
            f"spinodal dip STILL present at t_scan_hi={t_scan_hi:.3f} K (rho_v_sp={pair_hi[0]:.2f}, rho_l_sp={pair_hi[1]:.2f}) -- upper T bound is not high enough to be past the true critical temperature; widen t_scan_hi",
            0,
        )

    # Bisect on dip existence: t_lo always has a confirmed dip, t_hi never does.
    t_lo, t_hi = float(t_scan_lo), float(t_scan_hi)
    last_pair = pair_lo
    n_iter = 0
    for n_iter in range(1, CRIT_BISECT_MAX_ITERS + 1):
        if (t_hi - t_lo) < CRIT_BISECT_T_TOL_K:
            break
        t_mid = 0.5 * (t_lo + t_hi)
        pair_mid = _spinodal_pair_at_t(d1, d2, t_mid, z1, rho_scan_lo, rho_scan_hi)
        if pair_mid is not None:
            t_lo = t_mid
            last_pair = pair_mid
        else:
            t_hi = t_mid

    t_c = t_lo  # last T with a confirmed dip; bracket (t_hi-t_lo) is now tiny
    rho_v_sp, rho_l_sp = last_pair
    rho_c = 0.5 * (rho_v_sp + rho_l_sp)
    width = rho_l_sp - rho_v_sp

    dpdrho_final = _dp_drho_fd(d1, d2, t_c, rho_c, z1)
    _dp_check, d2p_final = _pressure_rho_derivatives_fd(d1, d2, t_c, rho_c, z1)
    r1 = abs(float(dpdrho_final) * rho_scale / p_scale)
    r2 = abs(float(d2p_final) * rho_scale * rho_scale / p_scale)
    width_frac = width / rho_scale
    ok = bool((t_hi - t_lo) < CRIT_BISECT_T_TOL_K * 2.0 and width_frac < CRIT_WIDTH_FRAC_GATE)
    st = mix_state(d1, d2, t_c, rho_c, z1)
    return {
        "T_K": t_c,
        "rho_molm3": float(rho_c),
        "P_Pa": float(st.p_pa),
        "h_Jmol": float(st.h_jmol),
        "x1": float(z1),
        "converged": ok,
        "dp_drho_resid": float(r1),
        "d2p_drho2_resid": float(r2),
        "spinodal_width_molm3": float(width),
        "iterations": int(n_iter),
        "notes": f"bisection on spinodal-dip existence (width-based, avoids noisy d2P/drho2 -- see diagnose_spinodal.py); final T bracket={t_hi - t_lo:.2e} K, spinodal width={width:.3f} mol/m^3 ({width_frac * 100.0:.3f}% of rho_scale)",
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
    crit_point : dict [mixed units]
        The mixture's true critical point at fixed composition z1, from
        solve_mixture_critical_point -- the single physical state where
        bubble and dew genuinely coincide, used by plot_envelope to
        close the dome rather than leave a gap at the top.

    Assumptions
    -----------
    - Composition is fixed overall for the envelope run.
    - Continuation updates guesses from previously converged states.

    Failure modes
    -------------
    - Individual temperature points may fail and be marked DIVERGED.
    - crit_point["converged"] may be False if the 2-equation
      criticality solve itself fails; callers must check this before
      trusting/plotting it.

    References
    ----------
    - Not a citation of a single external paper -- continuation/homotopy
      (seeding each new temperature from the previous point's converged
      solution) is a standard, generic numerical-continuation technique.
      See [BCL99] for the per-point solver each seed feeds into.

    Notes on numerical stability
    ----------------------------
    - Continuation greatly improves robustness versus independent per-point solves.
    """
    d1 = load_idaes_helmholtz_json(fluid1)  # pure-fluid IDAES-Helmholtz-JSON dict for fluid 1 (r1234ze)
    d2 = load_idaes_helmholtz_json(fluid2)  # pure-fluid IDAES-Helmholtz-JSON dict for fluid 2 (r227ea)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = w1_to_x1(w1, mw1, mw2)  # convert the FIXED overall mass fraction (e.g. w1=0.911 for R-515B) into the mole fraction used everywhere else in this file

    # Compute Tred_mix/vred_mix (Bell 2023 mixing-rule reducing point) now
    # -- used as the FALLBACK critical-point guess/window if the sweep
    # below doesn't produce usable near-critical bubble/dew data (see the
    # crit_point computation AFTER the loop). Deliberately not calling
    # solve_mixture_critical_point yet: a first attempt calling it here,
    # up front, with a window derived purely from rho_guess=1/vred_mix,
    # found a spurious root sitting right at the window's lower edge
    # (2959.94 mol/m^3 against a 2960.0 bound) -- BELOW the real
    # converged dew branch's own vapor density at T=380K (2864.2
    # mol/m^3). The window was clipping the very branch it needed to
    # resolve. Waiting until after the sweep lets us build a
    # window/guess from ACTUAL solved near-critical densities instead.
    tc1_g = float(d1["basic"]["Tc"])
    tc2_g = float(d2["basic"]["Tc"])
    rhoc1_mol_g = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol_g = float(d2["basic"]["rhoc"]) / mw2
    vc1_g = 1.0 / rhoc1_mol_g
    vc2_g = 1.0 / rhoc2_mol_g
    tred_mix_g, vred_mix_g = bell2023_Tred_vred(z1, 1.0 - z1, tc1_g, tc2_g, vc1_g, vc2_g, BELL_2023_R1234ZE_R227EA)

    # Rough initial density guesses for the very first temperature point:
    # 80% of a mole-fraction-weighted average of the two pure critical
    # densities for the liquid guess, and 1% of that for the vapor guess.
    # These are only ever used ONCE -- every subsequent temperature reuses
    # the previous point's converged solution as its seed (continuation,
    # see the loop below), which is far more robust than reseeding from this
    # rough guess at every point independently.
    rhoc1 = float(d1["basic"]["rhoc"]) / mw1
    rhoc2 = float(d2["basic"]["rhoc"]) / mw2
    rho_l_seed = 0.8 * (z1 * rhoc1 + (1.0 - z1) * rhoc2)
    rho_v_seed = 0.01 * rho_l_seed
    bubble_rho_l_guess = rho_l_seed
    bubble_rho_v_guess = rho_v_seed
    dew_rho_l_guess = rho_l_seed
    dew_rho_v_guess = rho_v_seed
    y1_guess = z1  # first vapor-composition guess for the bubble solve: assume near-ideal (y1~z1), reasonable for a near-azeotropic blend like R-515B
    x1_guess = z1  # first liquid-composition guess for the dew solve, same reasoning

    bubble_rows: List[Dict] = []
    dew_rows: List[Dict] = []
    for t_k in np.asarray(t_vals, dtype=float):
        # Bubble solve at this T, seeded from the PREVIOUS temperature's
        # converged bubble state (continuation).
        b = solve_bubble_at_t(d1, d2, float(t_k), z1, bubble_rho_l_guess, bubble_rho_v_guess, y1_guess)
        bubble_rows.append(b)
        if b["status"] == "CONVERGED":
            # Advance the seed for the NEXT temperature only if this one
            # actually converged -- otherwise keep reusing the last known-good
            # state rather than propagating a bad/DIVERGED guess forward.
            bubble_rho_l_guess = b["rho_l_molm3"]
            bubble_rho_v_guess = b["rho_v_molm3"]
            y1_guess = b["y1_vap"]

        # Dew solve at this T, seeded independently from its own continuation
        # chain (kept separate from the bubble chain since bubble and dew are
        # different equilibrium loci, even though R-515B's near-azeotropic
        # behavior means they end up numerically close).
        d = solve_dew_at_t(d1, d2, float(t_k), z1, dew_rho_l_guess, dew_rho_v_guess, x1_guess)
        dew_rows.append(d)
        if d["status"] == "CONVERGED":
            dew_rho_l_guess = d["rho_l_molm3"]
            dew_rho_v_guess = d["rho_v_molm3"]
            x1_guess = d["x1_liq"]

    # Now solve the mixture's true critical point, using the sweep's OWN
    # near-critical converged densities (if available) to build a much
    # better-grounded search window than the Bell-mixing-rule guess alone
    # (see the comment above the Tred_mix/vred_mix computation for why).
    # Spinodal points (where the critical-point solve searches) lie
    # BETWEEN the equilibrium bubble/dew densities at a given T, so the
    # last converged bubble rho_l (highest-T liquid) and dew rho_v
    # (highest-T vapor) bracket where the true spinodal branches -- and
    # thus the critical point -- should be found, with a modest margin
    # added on each side.
    last_bubble_conv = next((r for r in reversed(bubble_rows) if r["status"] == "CONVERGED"), None)
    last_dew_conv = next((r for r in reversed(dew_rows) if r["status"] == "CONVERGED"), None)
    if last_bubble_conv is not None and last_dew_conv is not None:
        rho_l_last = float(last_bubble_conv["rho_l_molm3"])
        rho_v_last = float(last_dew_conv["rho_v_molm3"])
        t_last = max(float(last_bubble_conv["T_K"]), float(last_dew_conv["T_K"]))
        crit_rho_guess = 0.5 * (rho_l_last + rho_v_last)  # midpoint of the closest real converged liquid/vapor densities -- a much better-grounded guess than 1/vred_mix alone
        crit_rho_scan_lo = 0.85 * rho_v_last  # margin BELOW the real vapor density, so the vapor-side spinodal (which sits above rho_v_last) is never clipped
        crit_rho_scan_hi = 1.15 * rho_l_last  # margin ABOVE the real liquid density, so the liquid-side spinodal (which sits below rho_l_last) is never clipped
        crit_t_scan_lo = t_last - 2.0  # small margin below the highest T actually reached, so the search doesn't start narrower than where real data already exists
        crit_t_scan_hi = tred_mix_g + CRIT_T_SCAN_HI_K  # keep the upper bound tied to the mixing-rule estimate, since the sweep itself never goes supercritical
        crit_point = solve_mixture_critical_point(
            d1, d2, z1,
            t_guess=tred_mix_g,
            rho_guess=crit_rho_guess,
            rho_scan_lo=crit_rho_scan_lo,
            rho_scan_hi=crit_rho_scan_hi,
            t_scan_lo=crit_t_scan_lo,
            t_scan_hi=crit_t_scan_hi,
        )
    else:
        # No usable near-critical convergence from the sweep (e.g. Tmax
        # was set far below critical) -- fall back to the pure
        # mixing-rule-based guess/default windows.
        crit_point = solve_mixture_critical_point(d1, d2, z1, t_guess=tred_mix_g, rho_guess=1.0 / vred_mix_g)

    return bubble_rows, dew_rows, z1, crit_point


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
    - Not an external literature citation -- uses Python's built-in `csv`
      module to serialize the dicts run_true_vle_envelope returns.

    Notes on numerical stability
    ----------------------------
    - Not applicable.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())  # column order = insertion order of the first row's dict (status, T_K, P_Pa, rho_l_molm3, ... -- see solve_bubble_at_t/solve_dew_at_t's return dicts)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)

QUALITY_LINE_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]  # matches the Honeywell TDS p-h chart's printed x=0.1...0.9 lines


def compute_quality_lines(
        bubble_rows: List[Dict],
        dew_rows: List[Dict],
        qualities: List[float] = QUALITY_LINE_VALUES,
)-> Dict[float,List[Dict]]:
    """
        Lever-rule interior two-phase states: at each T where BOTH bubble and
        dew converged, h(q) = h_l + q*(h_v - h_l), q in (0,1), at the shared
        P = avg(P_bubble, P_dew). Pure post-processing of already-solved
        rows -- never touches the solve itself, no new EOS calls.
    """
    lines: Dict[float, List[Dict]] = {q:[] for q in qualities}
    for b,d in zip(bubble_rows,dew_rows):
        if b["status"] != "CONVERGED" or d["status"] != "CONVERGED":
            continue
        if abs(b["T_K"] - d["T_K"])>1.0e-6:
            continue
        p_pa = 0.5*(b["P_Pa"] + d["P_Pa"])
        h_l = b["h_l_Jmol"]
        h_v = d["h_v_Jmol"]
        for q in qualities:
            lines[q].append(
                {
                    "T_K": float(b["T_K"]),
                    "P_Pa": float(p_pa),
                    "quality": float(q),
                    "h_Jmol": float(h_l + q*(h_v-h_l)),
                }
            )
    return lines

def plot_envelope(
    bubble_rows: List[Dict],
    dew_rows: List[Dict],
    mw_mix: float,
    out_fig: str | Path,
    crit_point: Optional[Dict] = None,
    quality_lines: Optional[Dict[float, List[Dict]]] = None,
) -> None:
    """
    Purpose
    -------
    Plot mixture VLE envelope in p-h space.

    Inputs
    ------
    bubble_rows, dew_rows : list[dict] [mixed units]
    mw_mix : float [kg/mol]
    out_fig : str | Path [filesystem path]
    crit_point : dict | None [mixed units]
        Output of solve_mixture_critical_point (via
        run_true_vle_envelope). If provided and crit_point["converged"]
        is True, its (h,P) coordinate is appended as the shared final
        vertex of BOTH the bubble and dew lines, so the two lines
        terminate at the exact same point instead of leaving a gap --
        this is a genuine solved closure (the true EOS critical point),
        not a cosmetic extrapolation. Rendered as one continuous solid
        black line per branch (no dashed near-critical overlay, no
        marker at the closure point) per explicit user preference
        (2026-08-13): "just need a smooth connection."

    Outputs
    -------
    None [unitless]
      Writes figure to disk.

    Assumptions
    -----------
    - Enthalpy conversion uses h[J/mol] -> h[kJ/kg] via mw_mix.
    - Bubble and dew are the only two loci plotted; this is a saturation-dome
      plot, not a flash/quality diagram.

    Failure modes
    -------------
    - File I/O/plotting backend exceptions propagate.
    - If crit_point is None or not converged, the dome is plotted exactly
      as before (near-critical segments end wherever the tapered gate
      stopped converging) -- no silent fabrication of a closure point.

    References
    ----------
    - [B23] and [LT99] -- the bubble/dew states plotted here come from that
      Helmholtz mixture EOS framework (see module KEY EQUATIONS REFERENCE
      block, tags M1-M6); this function is presentation-only.

    Notes on numerical stability
    ----------------------------
    - Filters to CONVERGED points to avoid plotting numerically invalid states.
    """
    # Only plot points that actually satisfied (M6) to the strict tolerance --
    # a DIVERGED row's densities/enthalpies are not meaningful and would
    # visually corrupt the dome if plotted.
    b = [r for r in bubble_rows if r["status"] == "CONVERGED"]
    d = [r for r in dew_rows if r["status"] == "CONVERGED"]

    # 2026-08-13: NOT split into far/near-critical subsets here (unlike the
    # design this file was copied from) -- per explicit user preference,
    # both branches are drawn as ONE continuous solid line each, no dashed
    # overlay distinguishing the tapered-gate closing segment. The taper
    # itself (see _rho_separation_min_ratio) still governs which points
    # are accepted as CONVERGED in the first place; only the VISUAL
    # demarcation is removed.

    # h_l_Jmol/h_v_Jmol are MOLAR enthalpies [J/mol]; divide by mw_mix
    # [kg/mol] to get MASS-specific enthalpy [J/kg], then *1e-3 for [kJ/kg]
    # (the conventional refrigerant p-H chart unit). P_Pa*1e-5 converts
    # [Pa] -> [bar] for the y-axis.
    hb = np.array([(r["h_l_Jmol"] / mw_mix) * 1e-3 for r in b], dtype=float)
    pb = np.array([r["P_Pa"] * 1e-5 for r in b], dtype=float)
    hd = np.array([(r["h_v_Jmol"] / mw_mix) * 1e-3 for r in d], dtype=float)
    pd = np.array([r["P_Pa"] * 1e-5 for r in d], dtype=float)

    # Append the solved critical point as the shared final vertex of BOTH
    # lines -- this is what actually closes the dome into a smooth,
    # unbroken curve (genuine solved coincidence, not a cosmetic patch).
    # Only done if the critical-point solve itself converged; otherwise
    # the dome is left exactly as the tapered gate produced it, with no
    # fabricated closure.
    have_crit = bool(crit_point is not None and crit_point.get("converged", False))
    if have_crit:
        hc = (crit_point["h_Jmol"] / mw_mix) * 1e-3  # [kJ/kg]
        pc = crit_point["P_Pa"] * 1e-5  # [bar]
        hb = np.append(hb, hc)
        pb = np.append(pb, pc)
        hd = np.append(hd, hc)
        pd = np.append(pd, pc)

    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if len(hb):
        ax.plot(hb, pb, lw=2.0, color="black", label="Bubble line (liq)")
    if len(hd):
        ax.plot(hd, pd, lw=2.0, color="black", label="Dew line (vap)")
    if quality_lines:
        # At Tc, h_l=h_v (that's the literal definition of the critical
        # point), so h(x)=h_l+x*(h_v-h_l) collapses to the SAME value for
        # every quality x there -- every quality line should terminate at
        # exactly the same (hc,pc) coordinate as the bubble/dew boundary's
        # own closure above, not wherever the T-sweep's last converged
        # point happened to land (see 2026-08-13 breadcrumb: "shouldn't
        # they meet at the critical point?", confirmed against Honeywell's
        # own chart, where the quality lines visibly converge at the apex).
        for q in sorted(quality_lines.keys()):
            rows_q = quality_lines[q]
            if not rows_q:
                continue
            hq = np.array([(r["h_Jmol"] / mw_mix) * 1e-3 for r in rows_q], dtype=float)
            pq = np.array([r["P_Pa"] * 1e-5 for r in rows_q], dtype=float)
            if have_crit:
                hq = np.append(hq, hc)
                pq = np.append(pq, pc)
            ax.plot(hq, pq, lw=0.75, color="black")
            ax.annotate(f"x={q:.1f}", (hq[0], pq[0]), fontsize=7, color="black", ha="right", va="center")
    # NOTE: no interior "quality line" / lever-rule shading is drawn here --
    # this plot is the saturation dome only (bubble + dew loci). Interior
    # two-phase state representation is out of scope for this file.
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [kJ/kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("R515B true VLE envelope (mu-equality solve)")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")


def plot_envelope_honeywell_units(
    bubble_rows: List[Dict],
    dew_rows: List[Dict],
    mw_mix: float,
    out_fig: str | Path,
    crit_point: Optional[Dict] = None,
    quality_lines: Optional[Dict[float, List[Dict]]] = None,
) -> None:
    """
    Purpose
    -------
    Validation-only twin of plot_envelope: same SI solve results, same
    single-black-line/no-marker presentation, but converted to the
    Honeywell TDS p-h chart's IP units (Btu/lbm, psia) at the very end via
    to_honeywell_units, for direct visual comparison against page 2 of the
    Honeywell Solstice N15 TDS.

    Inputs
    ------
    bubble_rows, dew_rows : list[dict] [mixed units]
    mw_mix : float [kg/mol]
    out_fig : str | Path [filesystem path]
    crit_point : dict | None [mixed units]

    Outputs
    -------
    None [unitless]
      Writes figure to disk.

    Assumptions
    -----------
    - This function exists ONLY for Honeywell-chart validation. The
      DVCT-facing pipeline uses plot_envelope (SI, kJ/kg / bar) and is
      completely unaffected by this function's existence -- per explicit
      user instruction (2026-08-13), unit conversion must never be baked
      into the core SI solve or into plot_envelope itself.
    - Same CONVERGED-only filtering and critical-point closure logic as
      plot_envelope; only the final unit conversion (via
      to_honeywell_units) and axis labels differ.
    - Does NOT reconcile the Honeywell chart's reference-state footnote
      (see to_honeywell_units docstring) -- shape/width comparisons are
      valid, absolute enthalpy position may still be offset.

    Failure modes
    -------------
    - File I/O/plotting backend exceptions propagate.
    - If crit_point is None or not converged, the dome is plotted exactly
      as the tapered gate produced it, with no fabricated closure.

    References
    ----------
    - Honeywell Solstice N15 (R-515B) Technical Data Sheet, p.2.
    - [B23] and [LT99] -- see plot_envelope; identical underlying data.

    Notes on numerical stability
    ----------------------------
    - Filters to CONVERGED points to avoid plotting numerically invalid states.
    """
    b = [r for r in bubble_rows if r["status"] == "CONVERGED"]
    d = [r for r in dew_rows if r["status"] == "CONVERGED"]

    # Same SI computation as plot_envelope (kJ/kg, Pa) -- conversion to
    # Honeywell's IP units happens ONLY in the final to_honeywell_units
    # call below, never upstream of this point.
    hb_si = [(r["h_l_Jmol"] / mw_mix) * 1e-3 for r in b]
    pb_si = [r["P_Pa"] for r in b]
    hd_si = [(r["h_v_Jmol"] / mw_mix) * 1e-3 for r in d]
    pd_si = [r["P_Pa"] for r in d]

    have_crit = bool(crit_point is not None and crit_point.get("converged", False))
    if have_crit:
        hc_si = (crit_point["h_Jmol"] / mw_mix) * 1e-3
        pc_si = crit_point["P_Pa"]
        hb_si = hb_si + [hc_si]
        pb_si = pb_si + [pc_si]
        hd_si = hd_si + [hc_si]
        pd_si = pd_si + [pc_si]

    def _convert_all(h_si_list: List[float], p_si_list: List[float]) -> Tuple[np.ndarray, np.ndarray]:
        """Apply to_honeywell_units point-by-point, at the very end -- the
        only place SI->IP conversion happens in this whole function."""
        h_ip: List[float] = []
        p_ip: List[float] = []
        for h_si, p_si in zip(h_si_list, p_si_list):
            h_conv, p_conv = to_honeywell_units(h_si, p_si)
            h_ip.append(h_conv)
            p_ip.append(p_conv)
        return np.array(h_ip, dtype=float), np.array(p_ip, dtype=float)

    hb, pb = _convert_all(hb_si, pb_si)
    hd, pd = _convert_all(hd_si, pd_si)

    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=160)
    if len(hb):
        ax.plot(hb, pb, lw=2.0, color="black", label="Bubble line (liq)")
    if len(hd):
        ax.plot(hd, pd, lw=2.0, color="black", label="Dew line (vap)")
    if quality_lines:
        # Same critical-point closure as plot_envelope -- all quality
        # lines terminate at the exact same (hc_si,pc_si) coordinate in SI,
        # BEFORE the final to_honeywell_units conversion (conversion stays
        # final-step-only, per standing instruction).
        for q in sorted(quality_lines.keys()):
            rows_q = quality_lines[q]
            if not rows_q:
                continue
            hq_si = [(r["h_Jmol"] / mw_mix) * 1e-3 for r in rows_q]
            pq_si = [r["P_Pa"] for r in rows_q]
            if have_crit:
                hq_si = hq_si + [hc_si]
                pq_si = pq_si + [pc_si]
            hq, pq = _convert_all(hq_si, pq_si)
            ax.plot(hq, pq, lw=0.75, color="black")
            ax.annotate(f"x={q:.1f}", (hq[0], pq[0]), fontsize=7, color="black", ha="right", va="center")
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy [Btu/lbm]")
    ax.set_ylabel("Pressure [psia]")
    ax.set_title("R515B true VLE envelope -- Honeywell TDS units (validation only)")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")


# =============================================================================
# HONEYWELL SOLSTICE N15 (R-515B) TECHNICAL DATA SHEET REVALIDATION
# Added 2026-08-13. Ground-truth reference table transcribed directly (via
# pdftotext -layout, not manual retyping) from the "PRESSURE AND ENTHALPY"
# table on page 3 of the Honeywell Solstice N15 (R-515B) Technical Data
# Sheet (filename 8077001-ras-tds-solstice-n15-ltr-en): Temperature (F) vs
# Pressure (psig), 0-170F in 2F steps, 86 points total. The datasheet lists
# this blend as "Zero glide" (see PHYSICAL PROPERTIES / BENEFITS on page 1),
# i.e. Honeywell itself treats bubble and dew pressure at a given T as
# coincident -- this is the basis for the x1=y1 simplification used in
# run_honeywell_pt_comparison below: the model's bubble and dew pressures
# are averaged into a single P_model per T and compared against this one
# chart pressure, rather than requiring a separate reference for each branch
# (the datasheet only publishes one P-T curve for this blend, consistent
# with its own "Zero glide" claim).
# =============================================================================
HONEYWELL_TDS_PT_TABLE_F_PSIG: List[Tuple[float, float]] = [
    (0.0, 0.7), (2.0, 1.5), (4.0, 2.3), (6.0, 3.1), (8.0, 3.9),
    (10.0, 4.8), (12.0, 5.7), (14.0, 6.6), (16.0, 7.6), (18.0, 8.6),
    (20.0, 9.6), (22.0, 10.7), (24.0, 11.8), (26.0, 13.0), (28.0, 14.2),
    (30.0, 15.4), (32.0, 16.6), (34.0, 17.9), (36.0, 19.3), (38.0, 20.7),
    (40.0, 22.1), (42.0, 23.6), (44.0, 25.1), (46.0, 26.7), (48.0, 28.3),
    (50.0, 29.9), (52.0, 31.6), (54.0, 33.4), (56.0, 35.2),
    (58.0, 37.1), (60.0, 39.0), (62.0, 40.9), (64.0, 42.9), (66.0, 45.0),
    (68.0, 47.2), (70.0, 49.3), (72.0, 51.6), (74.0, 53.9), (76.0, 56.2),
    (78.0, 58.7), (80.0, 61.2), (82.0, 63.7), (84.0, 66.3), (86.0, 69.0),
    (88.0, 71.7), (90.0, 74.6), (92.0, 77.4), (94.0, 80.4), (96.0, 83.4),
    (98.0, 86.5), (100.0, 89.7), (102.0, 92.9), (104.0, 96.2), (106.0, 99.6),
    (108.0, 103.1), (110.0, 106.6), (112.0, 110.3), (114.0, 114.0),
    (116.0, 117.8), (118.0, 121.6), (120.0, 125.6), (122.0, 129.6), (124.0, 133.8),
    (126.0, 138.0), (128.0, 142.3), (130.0, 146.7), (132.0, 151.2), (134.0, 155.7),
    (136.0, 160.4), (138.0, 165.2), (140.0, 170.1), (142.0, 175.0), (144.0, 180.1),
    (146.0, 185.3), (148.0, 190.5), (150.0, 195.9), (152.0, 201.4), (154.0, 207.0),
    (156.0, 212.7), (158.0, 218.5), (160.0, 224.4), (162.0, 230.4), (164.0, 236.5),
    (166.0, 242.8), (168.0, 249.2), (170.0, 255.7),
]
assert len(HONEYWELL_TDS_PT_TABLE_F_PSIG) == 86, "Honeywell TDS table should have 86 points (0-170F, 2F steps)"

PSI_TO_PA = 6894.757293168361  # exact SI conversion factor, 1 psi = 6894.757293168361 Pa
PSIG_TO_PSIA_OFFSET = 14.696  # standard atmospheric pressure [psia] used by the datasheet's gauge-pressure convention
BTU_LBM_TO_KJ_KG = 2.326  # exact: 1 Btu_IT/lbm = 2.326 kJ/kg (used to match the Honeywell TDS p-h chart's IP-unit axes, page 2)


def to_honeywell_units(h_kj_kg: float, p_pa: float) -> Tuple[float, float]:
    """
    Purpose
    -------
    Convert a single already-computed SI (h [kJ/kg], P [Pa]) state to the
    IP units used by the Honeywell Solstice N15 TDS p-h chart (page 2):
    enthalpy in Btu/lbm, pressure in psia (absolute, no gauge offset --
    unlike the page-3 P-T table's psig convention, the p-h chart's y-axis
    is explicitly labeled "Pressure (psia)").

    Inputs
    ------
    h_kj_kg : float [kJ/kg]
    p_pa : float [Pa]

    Outputs
    -------
    h_btu_lbm : float [Btu/lbm]
    p_psia : float [psia]

    Assumptions
    -----------
    - h_kj_kg, p_pa are FINAL, already-computed values from the SI solve
      (mix_state / run_true_vle_envelope / solve_mixture_critical_point).
      This function performs unit conversion ONLY -- it must never be
      called anywhere inside the solve itself, and plot_envelope (the
      DVCT-facing SI plot) must never call it. Per explicit user
      instruction (2026-08-13): "The conversion should be applied on
      final values only. We will use the code as is for DVCT, but
      validation should be against Honeywell."
    - Does NOT reconcile the Honeywell chart's reference-state footnote
      ("h = 200 kJ/kg, s = 1.00 kJ/kg-K; sat. liq. at 0C") against
      whatever zero-point is baked into the IDAES Helmholtz JSON this
      file's EOS pulls from -- units only, not reference-state alignment.
      See the 2026-08-13 breadcrumb entry for why this matters (absolute
      enthalpy values may still be offset even after this conversion).

    Failure modes
    -------------
    - None; pure arithmetic.

    References
    ----------
    - Honeywell Solstice N15 (R-515B) Technical Data Sheet, p.2
      ("PRESSURE AND ENTHALPY" p-h chart, axis units).
    """
    h_btu_lbm = h_kj_kg / BTU_LBM_TO_KJ_KG
    p_psia = p_pa / PSI_TO_PA
    return float(h_btu_lbm), float(p_psia)


def _f_to_k(t_f: float) -> float:
    """Convert Fahrenheit to Kelvin: K = (F-32)*5/9 + 273.15."""
    return float((t_f - 32.0) * 5.0 / 9.0 + 273.15)


def _psig_to_pa(p_psig: float) -> float:
    """Convert gauge pressure [psig] to absolute pressure [Pa] via +14.696 psia offset, then the exact psi->Pa factor."""
    return float((p_psig + PSIG_TO_PSIA_OFFSET) * PSI_TO_PA)


def run_honeywell_pt_comparison(
    fluid1: str,
    fluid2: str,
    w1: float,
    out_csv: str | Path = "verification/r515b_honeywell_pt_revalidation.csv",
    out_summary_json: str | Path = "verification/r515b_honeywell_pt_revalidation_summary.json",
) -> Dict:
    """
    Purpose
    -------
    Run the true-VLE solve at exactly the 86 temperatures in the Honeywell
    Solstice N15 TDS P-T table, compare the model's saturation pressure at
    each point against the datasheet's chart pressure, and report MAPE/bias/
    max-error accuracy metrics -- the real, current-file Honeywell
    revalidation (distinct from the older 1.84% MAPE result, which was from
    a different file, mixture_dome_validation.py).

    Inputs
    ------
    fluid1, fluid2 : str [unitless]
    w1 : float [kg/kg]
        Component-1 (R-1234ze(E)) mass fraction; 0.911 for R-515B.
    out_csv : str | Path [filesystem path]
    out_summary_json : str | Path [filesystem path]

    Outputs
    -------
    summary : dict [mixed units]
      MAPE_pct, bias_pct, max_abs_err_pct, convergence counts, and metadata.
      Also written to out_summary_json; the full per-point comparison is
      written to out_csv.

    Assumptions
    -----------
    - x1=y1 simplification: bubble and dew pressure at a given T are
      averaged into one P_model and compared against the datasheet's single
      published P-T curve, consistent with Honeywell's own "Zero glide"
      characterization of this blend (see module-level comment above this
      function). If only one of bubble/dew converges at a given T, that
      branch's pressure is used alone rather than discarding the point.
    - The 86-point table (HONEYWELL_TDS_PT_TABLE_F_PSIG) is already ascending
      in T, so passing it straight into run_true_vle_envelope's continuation
      sweep is valid (each point seeds from the previous, cooler point).

    Failure modes
    -------------
    - Individual T points may DIVERGE (see bubble_status/dew_status columns
      in the output CSV); MAPE/bias are computed over converged points only,
      and n_converged_points/n_total_points in the summary reports how many
      that was.

    References
    ----------
    - Honeywell Solstice N15 (R-515B) Technical Data Sheet, p.3 (ground-truth
      reference table).
    - [B23], [LT99] -- the model this table is validating (see module KEY
      EQUATIONS REFERENCE block).

    Notes on numerical stability
    ----------------------------
    - Inherits run_true_vle_envelope's continuation-seeded robustness; no
      additional solver logic is introduced here.
    """
    t_f_vals = np.array([p[0] for p in HONEYWELL_TDS_PT_TABLE_F_PSIG], dtype=float)
    p_psig_vals = np.array([p[1] for p in HONEYWELL_TDS_PT_TABLE_F_PSIG], dtype=float)
    t_k_vals = np.array([_f_to_k(t) for t in t_f_vals], dtype=float)
    p_ref_pa_vals = np.array([_psig_to_pa(p) for p in p_psig_vals], dtype=float)

    bubble_rows, dew_rows, z1, crit_point = run_true_vle_envelope(fluid1, fluid2, w1, t_k_vals)

    rows_out: List[Dict] = []
    abs_err_pct: List[float] = []
    signed_err_pct: List[float] = []
    for i in range(len(t_f_vals)):
        b = bubble_rows[i]
        d = dew_rows[i]
        b_ok = b["status"] == "CONVERGED"
        d_ok = d["status"] == "CONVERGED"
        p_bubble = float(b["P_Pa"]) if b_ok else float("nan")
        p_dew = float(d["P_Pa"]) if d_ok else float("nan")
        if b_ok and d_ok:
            p_avg = 0.5 * (p_bubble + p_dew)
        elif b_ok:
            p_avg = p_bubble
        elif d_ok:
            p_avg = p_dew
        else:
            p_avg = float("nan")

        p_ref = float(p_ref_pa_vals[i])
        err_pct = float((p_avg - p_ref) / p_ref * 100.0) if np.isfinite(p_avg) else float("nan")
        row = {
            "T_F": float(t_f_vals[i]),
            "T_K": float(t_k_vals[i]),
            "P_chart_psig": float(p_psig_vals[i]),
            "P_chart_Pa": p_ref,
            "P_model_bubble_Pa": p_bubble if np.isfinite(p_bubble) else None,
            "P_model_dew_Pa": p_dew if np.isfinite(p_dew) else None,
            "P_model_avg_Pa": p_avg if np.isfinite(p_avg) else None,
            "err_pct": err_pct if np.isfinite(err_pct) else None,
            "bubble_status": b["status"],
            "dew_status": d["status"],
            "bubble_notes": b["notes"],
            "dew_notes": d["notes"],
        }
        rows_out.append(row)
        if np.isfinite(err_pct):
            abs_err_pct.append(abs(err_pct))
            signed_err_pct.append(err_pct)

    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows_out[0].keys())
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)

    n_total = len(rows_out)
    n_converged = len(abs_err_pct)
    mape = float(np.mean(abs_err_pct)) if abs_err_pct else float("nan")
    bias = float(np.mean(signed_err_pct)) if signed_err_pct else float("nan")
    max_abs_err = float(np.max(abs_err_pct)) if abs_err_pct else float("nan")
    converged_rows = [r for r in rows_out if r["err_pct"] is not None]
    worst_t_f = None
    if abs_err_pct:
        worst_idx = int(np.argmax(abs_err_pct))
        worst_t_f = converged_rows[worst_idx]["T_F"]

    summary = {
        "fluid1": fluid1,
        "fluid2": fluid2,
        "w1_kgkg": float(w1),
        "z1_molmol": float(z1),
        "source": "Honeywell Solstice N15 (R-515B) Technical Data Sheet, p.3, Pressure/Enthalpy table (0-170F, 2F steps, psig)",
        "n_total_points": n_total,
        "n_converged_points": n_converged,
        "n_diverged_points": n_total - n_converged,
        "MAPE_pct": mape,
        "bias_pct": bias,
        "max_abs_err_pct": max_abs_err,
        "max_abs_err_at_T_F": worst_t_f,
        "x1_equals_y1_simplification": True,
        "critical_point": crit_point,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    out_summary_json = Path(out_summary_json)
    out_summary_json.parent.mkdir(parents=True, exist_ok=True)
    with out_summary_json.open("w") as f:
        json.dump(summary, f, indent=2)

    print(f"Honeywell P-T comparison: {n_converged}/{n_total} points converged.")
    if abs_err_pct:
        print(f"MAPE={mape:.4f}%  bias={bias:.4f}%  max_abs_err={max_abs_err:.4f}% (at T={worst_t_f}F)")
    print(f"Saved comparison CSV: {out_csv}")
    print(f"Saved summary JSON: {out_summary_json}")
    return summary


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
    - Not an external literature citation -- internal cross-reference: this
      is the argparse entrypoint that wires run_true_vle_envelope, save_csv,
      and plot_envelope together.

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
    p.add_argument("--fig-honeywell-units", default="verification/r515b_true_vle_envelope_honeywell_units.png", help="Validation-only twin of --fig, converted to the Honeywell TDS p-h chart's IP units (Btu/lbm, psia) at the final step only; DVCT-facing --fig stays SI.")
    p.add_argument("--metadata", default="verification/r515b_true_vle_metadata.json")
    p.add_argument(
        "--honeywell-compare",
        action="store_true",
        help="Run the 86-point Honeywell TDS P-T revalidation instead of a linspace dome sweep; ignores --Tmin/--Tmax/--n and the dome CSV/fig/metadata args.",
    )
    p.add_argument("--honeywell-csv", default="verification/r515b_honeywell_pt_revalidation.csv")
    p.add_argument("--honeywell-summary", default="verification/r515b_honeywell_pt_revalidation_summary.json")
    args = p.parse_args()

    if args.honeywell_compare:
        run_honeywell_pt_comparison(
            args.fluid1,
            args.fluid2,
            args.w1,
            out_csv=args.honeywell_csv,
            out_summary_json=args.honeywell_summary,
        )
        return

    t_vals = np.linspace(args.Tmin, args.Tmax, args.n)  # evenly-spaced temperature grid [K] the dome will be evaluated at
    # Run the full continuation sweep -- this is the single call that
    # produces the entire dome (bubble_rows and dew_rows, one entry per T)
    # plus the solved mixture critical point used to close the top of the dome.
    bubble_rows, dew_rows, z1, crit_point = run_true_vle_envelope(args.fluid1, args.fluid2, args.w1, t_vals)
    save_csv(bubble_rows, args.bubble_csv)
    save_csv(dew_rows, args.dew_csv)

    d1 = load_idaes_helmholtz_json(args.fluid1)
    d2 = load_idaes_helmholtz_json(args.fluid2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = z1 * mw1 + (1.0 - z1) * mw2  # overall mixture molecular weight [kg/mol], for the J/mol -> kJ/kg conversion in plot_envelope
    quality_lines = compute_quality_lines(bubble_rows, dew_rows)  # lever-rule interior states, pure post-processing of already-solved rows -- no new EOS calls
    plot_envelope(bubble_rows, dew_rows, mw_mix, args.fig, crit_point=crit_point, quality_lines=quality_lines)  # SI, DVCT-facing -- unaffected by the Honeywell-units conversion below
    plot_envelope_honeywell_units(bubble_rows, dew_rows, mw_mix, args.fig_honeywell_units, crit_point=crit_point, quality_lines=quality_lines)  # validation-only, IP units, converted at the final step only

    nb = int(sum(r["status"] == "CONVERGED" for r in bubble_rows))  # count of bubble points that met the strict (M6) convergence gate
    nd = int(sum(r["status"] == "CONVERGED" for r in dew_rows))  # count of dew points that met the strict (M6) convergence gate
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
        "critical_point": crit_point,  # from solve_mixture_critical_point: T_K, rho_molm3, P_Pa, h_Jmol, x1, converged, residuals, notes -- the solved (dP/drho)_T=(d2P/drho2)_T=0 state used to close the dome in the figure
    }
    mpath = Path(args.metadata)
    mpath.parent.mkdir(parents=True, exist_ok=True)
    with mpath.open("w") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved bubble CSV: {args.bubble_csv}")
    print(f"Saved dew CSV: {args.dew_csv}")
    print(f"Saved figure (SI, DVCT-facing): {args.fig}")
    print(f"Saved figure (Honeywell IP units, validation-only): {args.fig_honeywell_units}")
    print(f"Saved metadata: {args.metadata}")


if __name__ == "__main__":
    _cli()
