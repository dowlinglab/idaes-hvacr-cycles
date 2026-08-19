"""
r515b_helmholtz_core.py -- Stage F/G/H/I production core for R-515B
(R-1234ze(E)/R-227ea, w1=0.911), INDEPENDENT of the oracle at runtime.

MASTER TASK compliance notes (read this before editing)
---------------------------------------------------------
- This module NEVER imports `mixture_fully_validated.py` (the oracle). It
  imports ONLY from `linear_model_codex.py`, an existing, unmodified file
  explicitly permitted as a "legitimate production dependency" (spec rule 6),
  distinct from the oracle-import prohibition of spec rule 9.
- Every equation/algorithm/numerical-safeguard below is a deliberate,
  traced reproduction of the oracle's own validated behavior -- not an
  independent re-derivation from literature. Each function's docstring
  states exactly which oracle function it reproduces and how it was
  validated (see R515B_idaes_package/validate_*.py scripts and
  helmholtz_prop_validation.md Sections 2-9 for the numeric evidence).
- Numeric TUNING CONSTANTS below (solver tolerances, density-mapping
  bounds, loss-function parameters, the Honeywell entropy-reference offset)
  are copied verbatim from the oracle as literal numbers. This is
  reproducing validated numerical behavior (spec rule 10/33: "reproduce
  numerical/procedural behavior... preserve ALL numerical convergence
  machinery"), not importing oracle code -- these are physical/numerical
  facts already established by the oracle's own extensive validation
  (S-kink fix, near-critical taper, etc.), and changing any of them
  without a documented reason would be an undocumented deviation (spec
  rule 70), which is prohibited.
- Shared-code-path functions (`mixture_alpha0_alphar_derivs`,
  `bell2023_Tred_vred`, `alphar_idaes_with_derivs`, `bell2023_departure_
  alphar`, `bell2023_departure_base`, `_bell2023_reducing_derivs_binary`,
  `load_idaes_helmholtz_json`, `mw_from_json`) are imported directly from
  `linear_model_codex.py` -- the SAME functions the oracle itself imports
  and calls (confirmed via helmholtz_prop_validation.md Sections 2-4).

Traceability map (oracle function -> this module's function -> validation)
----------------------------------------------------------------------------
- _clip_x                        -> _clip_x                        (identical formula, EPS_X copied)
- _mix_alpha_and_derivs           -> _mix_alpha_and_derivs           (validate_table1_vs_oracle.py, exact match)
- mix_state                       -> mix_state                       (validate_table1_vs_oracle.py, rel_err 0)
- _mix_entropy_direct             -> mix_entropy_direct              (validate_table1_vs_oracle.py, rel_err 0, raw)
- chemical_potentials_analytic    -> chemical_potentials_analytic    (validate_fugacity_vs_oracle.py, rel_err 0)
- _sigmoid/_logit_from_unit_interval/_scaled_sigmoid/_inverse_scaled_sigmoid/_safe_residual_vector/_rho_param_to_states
                                  -> (same names)                    (identical formula, unit-tested below in __main__ self-check)
- _rho_separation_min_ratio       -> _rho_separation_min_ratio       (identical formula)
- solve_bubble_at_t               -> solve_bubble_at_t               (validate_vle_vs_oracle.py, Stage I)
- solve_dew_at_t                  -> solve_dew_at_t                  (validate_vle_vs_oracle.py, Stage I)
"""

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import least_squares, brentq, root

# `linear_model_codex.py` (the legitimate, non-oracle production dependency
# this module builds on) lives in the sibling read-only directory
# `R515B_props_validated/`. This sys.path addition is ordinary Python
# packaging/path plumbing to make that import resolve -- it does NOT import
# or touch `mixture_fully_validated.py` (the oracle) at all; nothing in this
# module ever imports that file.
_SIBLING_DIR = Path(__file__).parent.parent / "R515B_props_validated"
if str(_SIBLING_DIR) not in sys.path:
    sys.path.insert(0, str(_SIBLING_DIR))

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
    x1_from_w1,
    _bell2023_reducing_derivs_binary,
)
# NOTE: `_bell2023_reducing_derivs_binary` (imported above, shared code path)
# is numerically identical to the oracle's own LOCAL re-implementation
# `_bell2023_reducing_derivs_binary_local` -- confirmed via
# validate_fugacity_vs_oracle.py (mu1/mu2 rel_err 0 at 3 direct states),
# since that is the only function-level difference between the two paths'
# fugacity assembly (see helmholtz_prop_validation.md Section 7).


# =============================================================================
# Numeric constants -- copied verbatim from mixture_fully_validated.py (the
# oracle), as literal numbers, per spec rule 10/33 (reproduce validated
# numerical/procedural behavior, do not re-derive or re-tune). Source line
# numbers noted for traceability/audit.
# =============================================================================
EPS_X = 1e-12                       # oracle line 166
RHO_MIN_MOLM3 = 1e-9                # oracle line 167
RHO_MAX_MOLM3 = 2.0e4               # oracle line 168
RHO_MAP_MIN_MOLM3 = 5.0             # oracle line 180
RHO_MAP_MAX_MOLM3 = 2.0e4           # oracle line 181
DRHO_MAP_MIN_MOLM3 = 1.0            # oracle line 182
DRHO_MAP_MAX_MOLM3 = 2.0e4          # oracle line 183

LSQ_DIFF_STEP = 1e-6                # oracle line 190
LSQ_MAX_NFEV = 800                  # oracle line 191
LSQ_METHOD = "trf"                  # oracle line 192
LSQ_X_SCALE = 1.0                   # oracle line 193
LSQ_FTOL = 1.0e-14                  # oracle line 194
LSQ_XTOL = 1.0e-14                  # oracle line 195
LSQ_GTOL = 1.0e-14                  # oracle line 196
LSQ_LOSS = "huber"                  # oracle line 197 (S-kink fix, 2026-08-13 -- do not revert to "cauchy")
LSQ_F_SCALE = 3.0                   # oracle line 198

LSQ_METHOD_DEW = "trf"              # oracle line 201
LSQ_DIFF_STEP_DEW = 1e-8            # oracle line 202
LSQ_MAX_NFEV_DEW = 800              # oracle line 203
LSQ_X_SCALE_DEW = 2.0               # oracle line 204
LSQ_FTOL_DEW = 1.0e-10              # oracle line 205
LSQ_XTOL_DEW = 1.0e-10              # oracle line 206
LSQ_GTOL_DEW = 1.0e-10              # oracle line 207
LSQ_LOSS_DEW = "huber"              # oracle line 208
LSQ_F_SCALE_DEW = 3.0               # oracle line 209

RHO_SEP_RATIO_FAR = 3.0             # oracle line 232
RHO_SEP_RATIO_NEAR_TC = 1.05        # oracle line 233
RHO_SEP_TAPER_START_K = 30.0        # oracle line 234

# Honeywell Solstice N15 reference-state rebase (oracle line 2491):
# raw model reads s=1.012423 kJ/(kg*K) at Honeywell's stated reference state
# (T=273.15K, rho=1258.4 kg/m3, sat. liq. at 0C), vs Honeywell's stated
# s_ref=1.00 kJ/(kg*K). This constant is a CALIBRATION FACT about the
# underlying pure-fluid IDAES JSON data + Bell(2023) departure model, not a
# tunable numerical-method parameter -- copied verbatim, must not be
# re-derived independently (that would risk a different, undocumented
# rebase and silently break Honeywell-referenced targets downstream).
ENTROPY_REFERENCE_OFFSET_JMOLK = -1.4595

FLUID1, FLUID2 = "r1234ze", "r227ea"
PAIR_KEY = "r1234ze|r227ea"
W1_R515B = 0.911  # Honeywell Solstice N15 (R-515B) nominal mass fraction of R-1234ze(E)


def r515b_x1() -> float:
    """R-515B's fixed overall mole fraction of R-1234ze(E) at w1=0.911,
    computed via the shared `x1_from_w1` (same function the oracle's own
    `w1_to_x1`-based call sites are algebraically equivalent to -- both
    implement n_i=w_i/MW_i, x1=n1/(n1+n2))."""
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    return float(x1_from_w1(d1, d2, W1_R515B))


def _clip_x(x1: float) -> float:
    """Reproduces oracle's _clip_x exactly (clamp to (EPS_X, 1-EPS_X))."""
    return float(min(max(x1, EPS_X), 1.0 - EPS_X))


@dataclass
class MixState:
    """Reproduces oracle's MixState dataclass exactly (same fields)."""
    p_pa: float
    h_jmol: float
    g_jmol: float
    rho_mol: float
    rho_mass: float
    x1: float


def _mix_alpha_and_derivs(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> Tuple[float, float, float]:
    """
    Reproduces oracle's `_mix_alpha_and_derivs` exactly: assembles total
    alpha_mix (including the entropy-of-mixing term x1*ln(x1)+x2*ln(x2)),
    alpha_tau_mix, and ar_del_mix from `mixture_alpha0_alphar_derivs`
    (the SAME shared function the oracle calls). Validated bit-for-bit
    identical to the oracle's own output at a representative liquid state
    (see helmholtz_prop_validation.md Section 4): abs diff 0.0 on all three
    returned quantities.
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
        d1=d1, d2=d2, x1=x1, x2=x2, tau=tau, delta=delta,
        Tred=tred, rho_red_mol=rho_red_mol, pair_key=PAIR_KEY,
    )
    a0_mix = a0 + x1 * np.log(x1) + x2 * np.log(x2)
    return float(a0_mix + ar), float(a0_tau + ar_tau), float(ar_del)


def mix_state(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> MixState:
    """
    Reproduces oracle's `mix_state` exactly: P=rho*R*T*Z (Z=1+delta*ar_del),
    h/(RT)=1+tau*alpha_tau+delta*ar_del, g/(RT)=1+alpha+delta*ar_del.
    Validated: rel_err 0 for P and h at 3 direct states (subcooled liquid,
    superheated vapor, supercritical) vs. the oracle's own mix_state --
    see validate_table1_vs_oracle.py / helmholtz_prop_validation.md Section 5/6.
    """
    x1 = _clip_x(x1)
    x2 = 1.0 - x1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    mw_mix = x1 * mw1 + x2 * mw2

    alpha, alpha_tau, ar_del = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1)

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
        p_pa=float(p_pa), h_jmol=float(h_jmol), g_jmol=float(g_jmol),
        rho_mol=float(rho_mol), rho_mass=float(rho_mol * mw_mix), x1=float(x1),
    )


def mix_entropy_direct(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> float:
    """
    Reproduces oracle's `_mix_entropy_direct` exactly, INCLUDING the
    Honeywell rebase applied internally (s_raw = R_u*(tau*alpha_tau-alpha);
    return s_raw + ENTROPY_REFERENCE_OFFSET_JMOLK) -- confirmed by reading
    the oracle's own body (lines 2569-2570) that the offset is applied
    unconditionally before returning, not left to the caller. Validated:
    raw (pre-offset) value matches oracle to rel_err 0-1.4e-16 at 3 direct
    states (helmholtz_prop_validation.md Section 5).
    """
    x1c = _clip_x(x1)
    x2c = 1.0 - x1c
    alpha, alpha_tau, _ar_del = _mix_alpha_and_derivs(d1, d2, t_k, rho_mol, x1c)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred, _vred = bell2023_Tred_vred(x1c, x2c, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k
    s_raw = R_u * (tau * alpha_tau - alpha)
    return float(s_raw + ENTROPY_REFERENCE_OFFSET_JMOLK)


def chemical_potentials_analytic(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float) -> Tuple[float, float]:
    """
    Reproduces oracle's `chemical_potentials_analytic` exactly (M5):
    fugacity f_i = x_i*rho*R*T*exp[d(n*alphar)/dn_i], mu_i = R*T*ln(f_i).
    Uses `_bell2023_reducing_derivs_binary` (imported, shared code path --
    NOT the oracle's own local `_bell2023_reducing_derivs_binary_local`,
    but confirmed numerically identical to it, rel_err 0, see
    validate_fugacity_vs_oracle.py / helmholtz_prop_validation.md Section 7).
    Validated overall: mu1/mu2 rel_err 0 vs. oracle at 3 direct states.
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
        d1=d1, d2=d2, x1=x1, x2=x2, tau=tau, delta=delta,
        Tred=tred, rho_red_mol=rho_red_mol, pair_key=PAIR_KEY,
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

    dtred_dx1, dvred_dx1 = _bell2023_reducing_derivs_binary(x1, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
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


# =============================================================================
# Sigmoid/logit reparameterization + safe-residual machinery -- reproduces
# oracle's _sigmoid/_logit_from_unit_interval/_scaled_sigmoid/
# _inverse_scaled_sigmoid/_safe_residual_vector/_rho_param_to_states exactly
# (identical formulas).
# =============================================================================
def _sigmoid(z: float) -> float:
    if z >= 0:
        ez = np.exp(-z)
        return float(1.0 / (1.0 + ez))
    ez = np.exp(z)
    return float(ez / (1.0 + ez))


def _logit_from_unit_interval(y: float) -> float:
    yc = min(max(float(y), 1.0e-12), 1.0 - 1.0e-12)
    return float(np.log(yc / (1.0 - yc)))


def _scaled_sigmoid(z: float, lo: float, hi: float) -> float:
    if not (hi > lo):
        raise ValueError("scaled-sigmoid requires hi > lo")
    return float(lo + (hi - lo) * _sigmoid(z))


def _inverse_scaled_sigmoid(x: float, lo: float, hi: float) -> float:
    if not (hi > lo):
        raise ValueError("inverse scaled-sigmoid requires hi > lo")
    x_clip = min(max(float(x), lo), hi)
    y = (x_clip - lo) / (hi - lo)
    return _logit_from_unit_interval(y)


def _safe_residual_vector(vals: np.ndarray) -> np.ndarray:
    if np.all(np.isfinite(vals)):
        return vals
    return np.array([1.0e6, 1.0e6, 1.0e6], dtype=float)


def _rho_param_to_states(u0: float, u1: float) -> Tuple[float, float]:
    rho_v = _scaled_sigmoid(float(u0), RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3)
    drho = _scaled_sigmoid(float(u1), DRHO_MAP_MIN_MOLM3, DRHO_MAP_MAX_MOLM3)
    rho_l = min(max(rho_v + drho, rho_v * (1.0 + 1.0e-8)), RHO_MAX_MOLM3)
    rho_v = min(rho_v, rho_l * (1.0 - 1.0e-8))
    return rho_l, rho_v


def _rho_separation_min_ratio(t_k: float, tred_mix: float) -> float:
    dt = tred_mix - t_k
    if dt >= RHO_SEP_TAPER_START_K:
        return RHO_SEP_RATIO_FAR
    if dt <= 0.0:
        return RHO_SEP_RATIO_NEAR_TC
    frac = dt / RHO_SEP_TAPER_START_K
    return RHO_SEP_RATIO_NEAR_TC + frac * (RHO_SEP_RATIO_FAR - RHO_SEP_RATIO_NEAR_TC)


# =============================================================================
# Stage I: bubble/dew VLE solve -- reproduces oracle's solve_bubble_at_t /
# solve_dew_at_t exactly (same 3-unknown least_squares system (M6): P_l=P_v,
# mu1_l=mu1_v, mu2_l=mu2_v; same sigmoid/logit reparameterization; same
# 4-seed retry ladder; same tapered density-separation acceptance gate; same
# strict post-solve residual re-verification independent of sol.success).
# =============================================================================
def solve_bubble_at_t(d1: Dict, d2: Dict, t_k: float, z1: float,
                       rho_l0: float, rho_v0: float, y10: float) -> Dict:
    z1 = _clip_x(z1)
    z2 = 1.0 - z1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred_mix, _vred_mix = bell2023_Tred_vred(z1, z2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)

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

    def _attempt(rho_l_seed_in, rho_v_seed_in, y1_seed_in):
        rho_v_seed = min(max(float(rho_v_seed_in), RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.9)
        rho_l_seed = min(max(float(rho_l_seed_in), rho_v_seed * (1.0 + 1.0e-6)), RHO_MAX_MOLM3)
        dr_seed = max(rho_l_seed - rho_v_seed, 1.0e-6)
        y1_seed = _clip_x(y1_seed_in)
        u0 = np.array([
            _inverse_scaled_sigmoid(rho_v_seed, RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3),
            _inverse_scaled_sigmoid(dr_seed, DRHO_MAP_MIN_MOLM3, DRHO_MAP_MAX_MOLM3),
            np.log(y1_seed / (1.0 - y1_seed)),
        ], dtype=float)
        lb = np.array([-30.0, -30.0, -30.0], dtype=float)
        ub = np.array([30.0, 30.0, 30.0], dtype=float)
        sol = least_squares(res, u0, bounds=(lb, ub), method=LSQ_METHOD,
                             ftol=LSQ_FTOL, xtol=LSQ_XTOL, gtol=LSQ_GTOL,
                             max_nfev=LSQ_MAX_NFEV, x_scale=LSQ_X_SCALE,
                             diff_step=LSQ_DIFF_STEP, loss=LSQ_LOSS, f_scale=LSQ_F_SCALE)
        rho_l, rho_v = _rho_param_to_states(float(sol.x[0]), float(sol.x[1]))
        y1 = _clip_x(_sigmoid(float(sol.x[2])))
        st_l = mix_state(d1, d2, t_k, rho_l, z1)
        st_v = mix_state(d1, d2, t_k, rho_v, y1)
        mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, z1)
        mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, y1)
        r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
        r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)
        min_ratio = _rho_separation_min_ratio(t_k, tred_mix)
        ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * min_ratio)
        return {
            "status": "CONVERGED" if ok else "DIVERGED",
            "T_K": float(t_k),
            "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
            "rho_l_molm3": float(rho_l), "rho_v_molm3": float(rho_v),
            "x1_liq": float(z1), "y1_vap": float(y1),
            "h_l_Jmol": float(st_l.h_jmol), "h_v_Jmol": float(st_v.h_jmol),
            "r_P": float(r_p), "r_mu": float(r_mu),
            "iterations": int(sol.nfev), "notes": str(sol.message),
            "rho_sep_min_ratio": float(min_ratio),
            "near_critical": bool(min_ratio < RHO_SEP_RATIO_FAR),
        }

    y10 = _clip_x(y10)
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
        score = float(row["r_P"] + row["r_mu"])
        if score < best_score:
            best = row
            best_score = score
    assert best is not None
    best["notes"] = f"{best['notes']} | retry=best_failed"
    return best


def solve_dew_at_t(d1: Dict, d2: Dict, t_k: float, z1: float,
                    rho_l0: float, rho_v0: float, x10: float) -> Dict:
    z1 = _clip_x(z1)
    z2 = 1.0 - z1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred_mix, _vred_mix = bell2023_Tred_vred(z1, z2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)

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

    def _attempt(rho_l_seed_in, rho_v_seed_in, x1_seed_in):
        rho_v_seed = min(max(float(rho_v_seed_in), RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.9)
        rho_l_seed = min(max(float(rho_l_seed_in), rho_v_seed * (1.0 + 1.0e-6)), RHO_MAX_MOLM3)
        dr_seed = max(rho_l_seed - rho_v_seed, 1.0e-6)
        x1_seed = _clip_x(x1_seed_in)
        u0 = np.array([
            _inverse_scaled_sigmoid(rho_v_seed, RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3),
            _inverse_scaled_sigmoid(dr_seed, DRHO_MAP_MIN_MOLM3, DRHO_MAP_MAX_MOLM3),
            np.log(x1_seed / (1.0 - x1_seed)),
        ], dtype=float)
        lb = np.array([-30.0, -30.0, -30.0], dtype=float)
        ub = np.array([30.0, 30.0, 30.0], dtype=float)
        sol = least_squares(res, u0, bounds=(lb, ub), method=LSQ_METHOD_DEW,
                             ftol=LSQ_FTOL_DEW, xtol=LSQ_XTOL_DEW, gtol=LSQ_GTOL_DEW,
                             max_nfev=LSQ_MAX_NFEV_DEW, x_scale=LSQ_X_SCALE_DEW,
                             diff_step=LSQ_DIFF_STEP_DEW, loss=LSQ_LOSS_DEW, f_scale=LSQ_F_SCALE_DEW)
        rho_l, rho_v = _rho_param_to_states(float(sol.x[0]), float(sol.x[1]))
        x1 = _clip_x(_sigmoid(float(sol.x[2])))
        st_l = mix_state(d1, d2, t_k, rho_l, x1)
        st_v = mix_state(d1, d2, t_k, rho_v, z1)
        mu1_l, mu2_l = chemical_potentials_analytic(d1, d2, t_k, rho_l, x1)
        mu1_v, mu2_v = chemical_potentials_analytic(d1, d2, t_k, rho_v, z1)
        r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
        r_mu = max(abs(mu1_l - mu1_v), abs(mu2_l - mu2_v)) / (R_u * t_k)
        min_ratio = _rho_separation_min_ratio(t_k, tred_mix)
        ok = bool(sol.success and r_p <= 1e-6 and r_mu <= 1e-6 and rho_l > rho_v * min_ratio)
        return {
            "status": "CONVERGED" if ok else "DIVERGED",
            "T_K": float(t_k),
            "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
            "rho_l_molm3": float(rho_l), "rho_v_molm3": float(rho_v),
            "x1_liq": float(x1), "y1_vap": float(z1),
            "h_l_Jmol": float(st_l.h_jmol), "h_v_Jmol": float(st_v.h_jmol),
            "r_P": float(r_p), "r_mu": float(r_mu),
            "iterations": int(sol.nfev), "notes": str(sol.message),
            "rho_sep_min_ratio": float(min_ratio),
            "near_critical": bool(min_ratio < RHO_SEP_RATIO_FAR),
        }

    x10 = _clip_x(x10)
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


# =============================================================================
# Stage J: mixture critical-point solve -- reproduces oracle's THIRD-design
# `solve_mixture_critical_point` exactly: bisection on spinodal-dip
# EXISTENCE (not a noisy d2P/drho2 root-find -- see oracle docstring for the
# documented failure history of the first two abandoned designs). Constants
# copied verbatim (spec rule 10/33).
# =============================================================================
CRIT_RHO_SCAN_LO_FACTOR = 0.7
CRIT_RHO_SCAN_HI_FACTOR = 1.6
CRIT_RHO_SCAN_N = 150
CRIT_T_SCAN_LO_K = 6.0
CRIT_T_SCAN_HI_K = 3.0
CRIT_BISECT_MAX_ITERS = 50
CRIT_BISECT_T_TOL_K = 1.0e-4
CRIT_WIDTH_FRAC_GATE = 0.02


def _pressure_rho_derivatives_fd(d1: Dict, d2: Dict, t_k: float, rho_mol: float,
                                  x1: float, rel_step: float = 1.0e-4) -> Tuple[float, float]:
    """Reproduces oracle's `_pressure_rho_derivatives_fd` exactly (central FD
    on this module's own `mix_state`, which is itself validated identical
    to the oracle's `mix_state` -- see Section 5/6)."""
    h = max(rel_step * rho_mol, 1.0e-6)
    p_plus = mix_state(d1, d2, t_k, rho_mol + h, x1).p_pa
    p_minus = mix_state(d1, d2, t_k, rho_mol - h, x1).p_pa
    p_mid = mix_state(d1, d2, t_k, rho_mol, x1).p_pa
    dp_drho = (p_plus - p_minus) / (2.0 * h)
    d2p_drho2 = (p_plus - 2.0 * p_mid + p_minus) / (h * h)
    return float(dp_drho), float(d2p_drho2)


def _dp_drho_fd(d1: Dict, d2: Dict, t_k: float, rho_mol: float, x1: float, rel_step: float = 1.0e-4) -> float:
    """Reproduces oracle's `_dp_drho_fd` exactly (2-point central FD)."""
    h = max(rel_step * rho_mol, 1.0e-6)
    p_plus = mix_state(d1, d2, t_k, rho_mol + h, x1).p_pa
    p_minus = mix_state(d1, d2, t_k, rho_mol - h, x1).p_pa
    return float((p_plus - p_minus) / (2.0 * h))


def _find_all_sign_changes(func, x_lo: float, x_hi: float, n_scan: int = 30) -> List[Tuple[float, float, float, float]]:
    """Reproduces oracle's `_find_all_sign_changes` exactly."""
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


def _spinodal_pair_at_t(d1: Dict, d2: Dict, t_k: float, x1: float,
                         rho_scan_lo: float, rho_scan_hi: float,
                         rel_step: float = 1.0e-4, n_scan: int = CRIT_RHO_SCAN_N) -> Optional[Tuple[float, float]]:
    """Reproduces oracle's `_spinodal_pair_at_t` exactly: finds both
    spinodal densities (vapor-side, liquid-side) at fixed T via
    _find_all_sign_changes + brentq refinement on the first/last bracket."""
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


def solve_mixture_critical_point(d1: Dict, d2: Dict, z1: float, t_guess: float, rho_guess: float,
                                  rho_scan_lo: Optional[float] = None, rho_scan_hi: Optional[float] = None,
                                  t_scan_lo: Optional[float] = None, t_scan_hi: Optional[float] = None) -> Dict:
    """
    Reproduces oracle's `solve_mixture_critical_point` exactly (THIRD design:
    bisection on spinodal-dip EXISTENCE, not a noisy d2P/drho2 root-find --
    see the oracle's own docstring for the documented failure history of the
    two abandoned prior designs, preserved here by reproducing the same
    final design rather than re-deriving from scratch).
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
            "T_K": float(t_guess), "rho_molm3": float(rho_guess), "P_Pa": float(p_scale),
            "h_Jmol": float("nan"), "x1": float(z1), "converged": False,
            "dp_drho_resid": float("nan"), "d2p_drho2_resid": float("nan"),
            "spinodal_width_molm3": float("nan"), "iterations": iters, "notes": notes,
        }

    pair_lo = _spinodal_pair_at_t(d1, d2, t_scan_lo, z1, rho_scan_lo, rho_scan_hi)
    if pair_lo is None:
        return _fallback(
            f"no spinodal dip found at t_scan_lo={t_scan_lo:.3f} K -- lower T bound is already "
            f"at/above the mixture critical temperature, or rho window doesn't bracket it there", 0)
    pair_hi = _spinodal_pair_at_t(d1, d2, t_scan_hi, z1, rho_scan_lo, rho_scan_hi)
    if pair_hi is not None:
        return _fallback(
            f"spinodal dip STILL present at t_scan_hi={t_scan_hi:.3f} K -- upper T bound is not "
            f"high enough to be past the true critical temperature; widen t_scan_hi", 0)

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

    t_c = t_lo
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
        "T_K": t_c, "rho_molm3": float(rho_c), "P_Pa": float(st.p_pa), "h_Jmol": float(st.h_jmol),
        "x1": float(z1), "converged": ok, "dp_drho_resid": float(r1), "d2p_drho2_resid": float(r2),
        "spinodal_width_molm3": float(width), "iterations": int(n_iter),
        "notes": f"bisection on spinodal-dip existence; final T bracket={t_hi - t_lo:.2e} K, "
                 f"spinodal width={width:.3f} mol/m^3 ({width_frac * 100.0:.3f}% of rho_scale)",
    }


# =============================================================================
# Stage L (part 3) design support: pseudo-pure (x1=y1=z1 FIXED) saturation
# solve. NOT a port of anything in the oracle -- the oracle's own
# solve_bubble_at_t/solve_dew_at_t deliberately solve a 3-unknown system
# (rho_l, rho_v, y1 or x1) because they represent the fully rigorous VLE
# where liquid and vapor compositions are allowed to differ. This function
# instead implements the DELIBERATE SIMPLIFICATION adopted for the active
# IDAES property package's two-phase logic (per the 2026-08-18 user
# correction and rule-22 resolution recorded in helmholtz_prop_validation.md
# Section 22 / PROJECT_CONTEXT.md): treat R-515B as a single fixed
# composition on BOTH phases (x1=y1=z1), collapsing the phase-equilibrium
# problem to a single-composition saturation curve structurally identical
# to how a PURE fluid's dome is determined (e.g. general_helmholtz's own
# HelmholtzStateBlockData) -- 2 equations (mechanical equilibrium P_l=P_v,
# and the Maxwell/equal-molar-Gibbs-energy condition g_l=g_v) in 2 unknowns
# (rho_l, rho_v) at a given T, instead of 3 equations/3 unknowns at fixed
# T with a free composition split.
#
# Why g_l=g_v (not the mu1/mu2 equalities) is the correct second condition
# here: for a binary mixture, g = x1*mu1 + x2*mu2 (Euler relation). When x1
# is allowed to differ between phases (the oracle's rigorous case), BOTH
# mu1_l=mu1_v AND mu2_l=mu2_v are required (2 independent conditions) to
# properly determine the equilibrium composition split. But here x1 is NOT
# a free unknown -- it is FIXED equal on both sides by construction (the
# whole point of this simplification) -- so the compositional degrees of
# freedom that made 2 separate mu-equalities necessary no longer exist.
# With x fixed identical on both phases, mu1_l=mu1_v and mu2_l=mu2_v would
# both trivially reduce to the SAME single condition IF x really were the
# equilibrium-consistent split (i.e. at a genuine binary azeotrope point);
# away from that exact point, forcing both individually would 3-equation-
# overdetermine a 2-unknown (rho_l,rho_v) system. The single g_l=g_v
# condition (the same Maxwell-construction condition used to fix a pure
# fluid's saturation pressure at given T) is what correctly and
# consistently determines the ONE (rho_l,rho_v) pair, at this fixed x, for
# which liquid and vapor are in equilibrium at temperature T -- this is
# the same math a pure-component EOS uses, just evaluated at R-515B's own
# fixed pseudo-pure composition instead of a literal single substance.
# =============================================================================
def solve_pseudopure_saturation_at_t(d1: Dict, d2: Dict, z1: float, t_k: float,
                                      rho_l0: float, rho_v0: float) -> Dict:
    """
    Pseudo-pure (x1=y1=z1 fixed) saturation solve at fixed T: find
    (rho_l, rho_v) such that P(T,rho_l,z1)=P(T,rho_v,z1) [mechanical
    equilibrium] and g(T,rho_l,z1)=g(T,rho_v,z1) [Maxwell/equal-Gibbs
    condition -- see module-level comment above for the derivation of why
    this, not mu1/mu2 equality, is the correct second condition once
    composition is fixed identical on both phases]. Uses the SAME sigmoid
    rho-parameterization, retry-seed pattern, and tapered density-
    separation gate as `solve_bubble_at_t`/`solve_dew_at_t` (Stage I) for
    consistency and to reuse already-validated numerical-stability
    machinery -- but a 2-unknown (not 3-unknown) least_squares system,
    since there is no composition unknown here.
    """
    z1 = _clip_x(z1)
    z2 = 1.0 - z1
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol
    tred_mix, _vred_mix = bell2023_Tred_vred(z1, z2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)

    def res(u):
        rho_l, rho_v = _rho_param_to_states(float(u[0]), float(u[1]))
        try:
            st_l = mix_state(d1, d2, t_k, rho_l, z1)
            st_v = mix_state(d1, d2, t_k, rho_v, z1)
            r1 = (st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
            r2 = (st_l.g_jmol - st_v.g_jmol) / (R_u * t_k)
            return _safe_residual_vector(np.array([r1, r2, 0.0], dtype=float))[:2]
        except Exception:
            return np.array([1.0e6, 1.0e6], dtype=float)

    def _attempt(rho_l_seed_in: float, rho_v_seed_in: float) -> Dict:
        rho_v_seed = min(max(float(rho_v_seed_in), RHO_MIN_MOLM3), RHO_MAX_MOLM3 * 0.9)
        rho_l_seed = min(max(float(rho_l_seed_in), rho_v_seed * (1.0 + 1.0e-6)), RHO_MAX_MOLM3)
        dr_seed = max(rho_l_seed - rho_v_seed, 1.0e-6)
        u0 = np.array([
            _inverse_scaled_sigmoid(rho_v_seed, RHO_MAP_MIN_MOLM3, RHO_MAP_MAX_MOLM3),
            _inverse_scaled_sigmoid(dr_seed, DRHO_MAP_MIN_MOLM3, DRHO_MAP_MAX_MOLM3),
        ], dtype=float)
        lb = np.array([-30.0, -30.0], dtype=float)
        ub = np.array([30.0, 30.0], dtype=float)
        sol = least_squares(res, u0, bounds=(lb, ub), method=LSQ_METHOD,
                             ftol=LSQ_FTOL, xtol=LSQ_XTOL, gtol=LSQ_GTOL,
                             max_nfev=LSQ_MAX_NFEV, x_scale=LSQ_X_SCALE,
                             diff_step=LSQ_DIFF_STEP, loss=LSQ_LOSS, f_scale=LSQ_F_SCALE)
        rho_l, rho_v = _rho_param_to_states(float(sol.x[0]), float(sol.x[1]))
        st_l = mix_state(d1, d2, t_k, rho_l, z1)
        st_v = mix_state(d1, d2, t_k, rho_v, z1)
        r_p = abs(st_l.p_pa - st_v.p_pa) / max(1.0, 0.5 * (st_l.p_pa + st_v.p_pa))
        r_g = abs(st_l.g_jmol - st_v.g_jmol) / (R_u * t_k)
        min_ratio = _rho_separation_min_ratio(t_k, tred_mix)
        ok = bool(sol.success and r_p <= 1e-6 and r_g <= 1e-6 and rho_l > rho_v * min_ratio)
        return {
            "status": "CONVERGED" if ok else "DIVERGED",
            "T_K": float(t_k),
            "P_Pa": float(0.5 * (st_l.p_pa + st_v.p_pa)),
            "rho_l_molm3": float(rho_l), "rho_v_molm3": float(rho_v),
            "x1": float(z1),
            "h_l_Jmol": float(st_l.h_jmol), "h_v_Jmol": float(st_v.h_jmol),
            "r_P": float(r_p), "r_g": float(r_g),
            "iterations": int(sol.nfev), "notes": str(sol.message),
            "rho_sep_min_ratio": float(min_ratio),
            "near_critical": bool(min_ratio < RHO_SEP_RATIO_FAR),
        }

    attempts = [
        (rho_l0, rho_v0),
        (1.1 * rho_l0, 0.7 * rho_v0),
        (0.9 * rho_l0, 0.5 * rho_v0),
        (1.2 * rho_l0, 0.4 * rho_v0),
    ]
    best = None
    best_score = np.inf
    for idx, (rl, rv) in enumerate(attempts):
        row = _attempt(rl, rv)
        if row["status"] == "CONVERGED":
            row["notes"] = f"{row['notes']} | retry={idx}"
            return row
        score = float(row["r_P"] + row["r_g"])
        if score < best_score:
            best = row
            best_score = score
    assert best is not None
    best["notes"] = f"{best['notes']} | retry=best_failed"
    return best


# =============================================================================
# Stage K (part 1): quality lines + isotherms -- reproduces oracle's
# compute_quality_lines / compute_isotherms_two_phase / _liquid_side /
# _vapor_side / _supercritical exactly. All are pure post-processing (lever
# rule) or single-equation/single-unknown brentq root-finds on this
# module's own already-validated `mix_state` -- no new EOS machinery.
# =============================================================================
def _f_to_k(t_f: float) -> float:
    """Reproduces oracle's `_f_to_k` exactly."""
    return float((t_f - 32.0) * 5.0 / 9.0 + 273.15)


QUALITY_LINE_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
ISOTHERM_VALUES_F = [-20.0, 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0]
VAPOR_ONLY_ISOTHERM_VALUES_F = [240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0]
LIQUID_EXT_P_MAX_PA = 1350.0 * 6894.757293168361
LIQUID_EXT_N_POINTS = 25
VAPOR_EXT_P_MIN_PA = 15.0 * 6894.757293168361
VAPOR_EXT_N_POINTS = 40
SUPERCRIT_P_MIN_PA = 15.0 * 6894.757293168361
SUPERCRIT_N_POINTS = 40


def compute_quality_lines(bubble_rows: List[Dict], dew_rows: List[Dict],
                           qualities: List[float] = QUALITY_LINE_VALUES) -> Dict[float, List[Dict]]:
    """Reproduces oracle's `compute_quality_lines` exactly: lever-rule
    interior two-phase states h(q) = h_l + q*(h_v-h_l) at shared
    P=avg(P_bubble,P_dew), pure post-processing of already-solved rows."""
    lines: Dict[float, List[Dict]] = {q: [] for q in qualities}
    for b, d in zip(bubble_rows, dew_rows):
        if b["status"] != "CONVERGED" or d["status"] != "CONVERGED":
            continue
        if abs(b["T_K"] - d["T_K"]) > 1.0e-6:
            continue
        p_pa = 0.5 * (b["P_Pa"] + d["P_Pa"])
        h_l = b["h_l_Jmol"]
        h_v = d["h_v_Jmol"]
        for q in qualities:
            lines[q].append({
                "T_K": float(b["T_K"]), "P_Pa": float(p_pa),
                "quality": float(q), "h_Jmol": float(h_l + q * (h_v - h_l)),
            })
    return lines


def compute_isotherms_two_phase(d1: Dict, d2: Dict, z1: float, bubble_rows: List[Dict],
                                 dew_rows: List[Dict], t_values_f: List[float] = ISOTHERM_VALUES_F) -> Dict[float, Dict]:
    """Reproduces oracle's `compute_isotherms_two_phase` exactly: solves
    FRESH at each exact target T via solve_bubble_at_t/solve_dew_at_t
    (this module's own validated Stage I ports), seeded from the nearest
    already-converged sweep point."""
    converged_bubble = [r for r in bubble_rows if r["status"] == "CONVERGED"]
    converged_dew = [r for r in dew_rows if r["status"] == "CONVERGED"]
    isotherms: Dict[float, Dict] = {}
    for t_f in t_values_f:
        t_k = _f_to_k(t_f)
        nearest_b = min(converged_bubble, key=lambda r: abs(r["T_K"] - t_k), default=None)
        nearest_d = min(converged_dew, key=lambda r: abs(r["T_K"] - t_k), default=None)
        if nearest_b is None or nearest_d is None:
            continue
        b = solve_bubble_at_t(d1, d2, t_k, z1, nearest_b["rho_l_molm3"], nearest_b["rho_v_molm3"], nearest_b["y1_vap"])
        d = solve_dew_at_t(d1, d2, t_k, z1, nearest_d["rho_l_molm3"], nearest_d["rho_v_molm3"], nearest_d["x1_liq"])
        if b["status"] != "CONVERGED" or d["status"] != "CONVERGED":
            continue
        isotherms[t_f] = {
            "T_F": t_f, "T_K": t_k, "P_Pa": 0.5 * (b["P_Pa"] + d["P_Pa"]),
            "h_l_Jmol": b["h_l_Jmol"], "h_v_Jmol": d["h_v_Jmol"],
            "rho_l_molm3": b["rho_l_molm3"], "rho_v_molm3": d["rho_v_molm3"],
        }
    return isotherms


def compute_isotherms_liquid_side(d1: Dict, d2: Dict, z1: float, isotherms: Dict[float, Dict],
                                   p_max_pa: float = LIQUID_EXT_P_MAX_PA, n_points: int = LIQUID_EXT_N_POINTS) -> Dict[float, List[Dict]]:
    """Reproduces oracle's `compute_isotherms_liquid_side` exactly: single-
    equation brentq root-find (mix_state(T,rho,z1).p_pa==P_target),
    log-spaced (geomspace) pressure sampling, seeded upward from rho_l_sat."""
    ext: Dict[float, List[Dict]] = {}
    for t_f, iso in isotherms.items():
        t_k = iso["T_K"]
        p_sat = iso["P_Pa"]
        rho_l_sat = iso["rho_l_molm3"]
        if p_max_pa <= p_sat:
            continue
        points: List[Dict] = [{"P_Pa": float(p_sat), "h_Jmol": float(iso["h_l_Jmol"])}]
        rho_seed = rho_l_sat
        for p_target in np.geomspace(p_sat, p_max_pa, n_points)[1:]:
            rho_lo = rho_seed
            rho_hi = min(rho_seed * 1.5, RHO_MAX_MOLM3)

            def g(rho: float, t_k=t_k, p_target=p_target) -> float:
                return mix_state(d1, d2, t_k, rho, z1).p_pa - p_target

            try:
                rho_root = brentq(g, rho_lo, rho_hi, xtol=1.0e-6, rtol=1.0e-12, maxiter=100)
            except ValueError:
                break
            h_root = mix_state(d1, d2, t_k, rho_root, z1).h_jmol
            points.append({"P_Pa": float(p_target), "h_Jmol": float(h_root)})
            rho_seed = rho_root
        if len(points) > 1:
            ext[t_f] = points
    return ext


def compute_isotherms_vapor_side(d1: Dict, d2: Dict, z1: float, isotherms: Dict[float, Dict],
                                  p_min_pa: float = VAPOR_EXT_P_MIN_PA, n_points: int = VAPOR_EXT_N_POINTS) -> Dict[float, List[Dict]]:
    """Reproduces oracle's `compute_isotherms_vapor_side` exactly (mirror
    image of the liquid-side extension, walking P down from p_sat)."""
    ext: Dict[float, List[Dict]] = {}
    for t_f, iso in isotherms.items():
        t_k = iso["T_K"]
        p_sat = iso["P_Pa"]
        rho_v_sat = iso["rho_v_molm3"]
        if p_min_pa >= p_sat:
            continue
        points: List[Dict] = [{"P_Pa": float(p_sat), "h_Jmol": float(iso["h_v_Jmol"])}]
        rho_seed = rho_v_sat
        for p_target in np.geomspace(p_sat, p_min_pa, n_points)[1:]:
            rho_hi = rho_seed
            rho_lo = max(rho_seed / 1.5, 1.0e-3)

            def g(rho: float, t_k=t_k, p_target=p_target) -> float:
                return mix_state(d1, d2, t_k, rho, z1).p_pa - p_target

            try:
                rho_root = brentq(g, rho_lo, rho_hi, xtol=1.0e-6, rtol=1.0e-12, maxiter=100)
            except ValueError:
                break
            h_root = mix_state(d1, d2, t_k, rho_root, z1).h_jmol
            points.append({"P_Pa": float(p_target), "h_Jmol": float(h_root)})
            rho_seed = rho_root
        if len(points) > 1:
            ext[t_f] = points
    return ext


def compute_isotherms_supercritical(d1: Dict, d2: Dict, z1: float,
                                     t_values_f: List[float] = VAPOR_ONLY_ISOTHERM_VALUES_F,
                                     p_min_pa: float = SUPERCRIT_P_MIN_PA, p_max_pa: float = LIQUID_EXT_P_MAX_PA,
                                     n_points: int = SUPERCRIT_N_POINTS) -> Dict[float, List[Dict]]:
    """Reproduces oracle's `compute_isotherms_supercritical` exactly:
    entire chart pressure range in one pass, first point seeded from an
    ideal-gas estimate (only to size the bracket), every subsequent point
    reuses the previous solved density."""
    R_GAS = 8.314462618
    ext: Dict[float, List[Dict]] = {}
    for t_f in t_values_f:
        t_k = _f_to_k(t_f)
        points: List[Dict] = []
        rho_seed: Optional[float] = None
        for p_target in np.geomspace(p_min_pa, p_max_pa, n_points):
            if rho_seed is None:
                rho_guess = p_target / (R_GAS * t_k)
                rho_lo = max(rho_guess * 0.2, 1.0e-3)
                rho_hi = rho_guess * 5.0
            else:
                rho_lo = rho_seed * 0.5
                rho_hi = min(rho_seed * 2.0, RHO_MAX_MOLM3)

            def g(rho: float, t_k=t_k, p_target=p_target) -> float:
                return mix_state(d1, d2, t_k, rho, z1).p_pa - p_target

            try:
                rho_root = brentq(g, rho_lo, rho_hi, xtol=1.0e-6, rtol=1.0e-12, maxiter=100)
            except ValueError:
                try:
                    rho_root = brentq(g, rho_lo * 0.1, rho_hi * 5.0, xtol=1.0e-6, rtol=1.0e-12, maxiter=100)
                except ValueError:
                    break
            h_root = mix_state(d1, d2, t_k, rho_root, z1).h_jmol
            points.append({"P_Pa": float(p_target), "h_Jmol": float(h_root)})
            rho_seed = rho_root
        if len(points) > 1:
            ext[t_f] = points
    return ext


# =============================================================================
# Stage K (part 2): isentrope machinery -- reproduces
# `mixture_isentrope_validation.py`'s isentrope functions exactly, NOT the
# base oracle `mixture_fully_validated.py`'s version. This is a deliberate
# choice: `mixture_isentrope_validation.py` is byte-identical to the oracle
# except for one addition -- `_verify_isentrope_solution`, the explicit
# post-hoc (P,s) residual re-check added earlier this session (after
# cross-checking the R1234yf sister project's own documented Bug #3) and
# empirically confirmed to introduce no regression (bubble.csv/dew.csv
# byte-identical, p-H diagrams MD5-identical, before/after). Since that
# insurance is validated, accepted, and part of "the established working
# model" per the user's own "Add insurance" go-ahead earlier this session,
# it is the version reproduced here -- porting the base oracle's version
# (bare sol.success, no residual re-check) would mean deliberately omitting
# an already-validated improvement, which is not what "reproduce the
# established working model's behavior" should mean in this case.
# =============================================================================
BTU_LBMR_TO_JKGK = 4186.8
ISENTROPE_VALUES_BTU_LBMR = [0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49]
ISENTROPE_P_MAX_PA = LIQUID_EXT_P_MAX_PA
ISENTROPE_P_MIN_PA = VAPOR_EXT_P_MIN_PA
ISENTROPE_N_POINTS = 40
ISENTROPE_VAPOR_MAX_ENTROPY_STEP_JMOLK = 2.0


def compute_isentropes_two_phase(d1: Dict, d2: Dict, bubble_rows: List[Dict], dew_rows: List[Dict],
                                  s_values_jmolK: List[float]) -> Dict[float, List[Dict]]:
    """Reproduces `mixture_isentrope_validation.py`'s `compute_isentropes_two_phase`
    exactly: pure post-processing lever rule on already-solved rows, no new
    EOS solve."""
    lines: Dict[float, List[Dict]] = {s: [] for s in s_values_jmolK}
    for b, d in zip(bubble_rows, dew_rows):
        if b["status"] != "CONVERGED" or d["status"] != "CONVERGED":
            continue
        if abs(b["T_K"] - d["T_K"]) > 1.0e-6:
            continue
        t_k = b["T_K"]
        s_l = mix_entropy_direct(d1, d2, t_k, b["rho_l_molm3"], b["x1_liq"])
        s_v = mix_entropy_direct(d1, d2, t_k, d["rho_v_molm3"], d["y1_vap"])
        if s_v <= s_l:
            continue
        p_pa = 0.5 * (b["P_Pa"] + d["P_Pa"])
        h_l = b["h_l_Jmol"]
        h_v = d["h_v_Jmol"]
        for s_target in s_values_jmolK:
            if not (s_l <= s_target <= s_v):
                continue
            x = (s_target - s_l) / (s_v - s_l)
            lines[s_target].append({
                "T_K": float(t_k), "P_Pa": float(p_pa),
                "h_Jmol": float(h_l + x * (h_v - h_l)), "quality": float(x),
            })
    return lines


def _isentrope_2eq_residual(vars_, d1: Dict, d2: Dict, z1: float, p_target: float, s_target: float) -> List[float]:
    """Reproduces `_isentrope_2eq_residual` exactly (normalized/relative P,s residuals)."""
    t_k, rho = vars_
    if t_k <= 0 or rho <= 0:
        return [1.0e12, 1.0e12]
    st = mix_state(d1, d2, t_k, rho, z1)
    s = mix_entropy_direct(d1, d2, t_k, rho, z1)
    r_p = (st.p_pa - p_target) / p_target
    r_s = (s - s_target) / s_target
    return [r_p, r_s]


def _verify_isentrope_solution(sol, d1: Dict, d2: Dict, z1: float, p_target: float, s_target: float,
                                p_tol_rel: float = 1.0e-6, s_tol_rel: float = 1.0e-6):
    """Reproduces `mixture_isentrope_validation.py`'s `_verify_isentrope_solution`
    exactly: explicit post-hoc (P,s) residual re-check, not bare sol.success
    (this project's own R1234yf-Bug#3-inspired insurance fix)."""
    if not sol.success:
        return None
    t_k, rho = float(sol.x[0]), float(sol.x[1])
    if t_k <= 0 or rho <= 0:
        return None
    st = mix_state(d1, d2, t_k, rho, z1)
    s_actual = mix_entropy_direct(d1, d2, t_k, rho, z1)
    r_p = abs(st.p_pa - p_target) / max(1.0, abs(p_target))
    r_s = abs(s_actual - s_target) / max(1.0, abs(s_target))
    if r_p > p_tol_rel or r_s > s_tol_rel:
        return None
    return t_k, rho, st


def compute_isentrope_liquid_side(d1: Dict, d2: Dict, z1: float, bubble_rows: List[Dict],
                                   s_values_jmolK: List[float], p_max_pa: float = ISENTROPE_P_MAX_PA,
                                   n_points: int = ISENTROPE_N_POINTS, crit_point: Optional[Dict] = None) -> Dict[float, List[Dict]]:
    """Reproduces `compute_isentrope_liquid_side` exactly: bubble-row anchor
    selection (or critical-point fallback anchor), warm-start solve at the
    anchor's own P, then geomspace-walked 2-eq (P,s) solves outward to
    p_max_pa, each verified via `_verify_isentrope_solution`."""
    converged_bubble = [r for r in bubble_rows if r["status"] == "CONVERGED"]
    have_crit = bool(crit_point is not None and crit_point.get("converged", False))
    ext: Dict[float, List[Dict]] = {s: [] for s in s_values_jmolK}
    for s_target in s_values_jmolK:
        best_row = None
        best_gap = None
        for r in converged_bubble:
            t_k = r["T_K"]
            s_l = mix_entropy_direct(d1, d2, t_k, r["rho_l_molm3"], z1)
            if s_l <= s_target:
                continue
            gap = s_l - s_target
            if best_gap is None or gap < best_gap:
                best_gap = gap
                best_row = r

        if best_row is not None:
            t_k = best_row["T_K"]
            rho_seed = best_row["rho_l_molm3"]
            p_seed = best_row["P_Pa"]
        elif have_crit:
            t_crit = crit_point["T_K"]
            rho_crit = crit_point["rho_molm3"]
            x1_crit = crit_point.get("x1", z1)
            s_crit = mix_entropy_direct(d1, d2, t_crit, rho_crit, x1_crit)
            if s_crit <= s_target:
                continue
            t_k = t_crit
            rho_seed = rho_crit
            p_seed = crit_point["P_Pa"]
        else:
            continue

        sol0 = root(_isentrope_2eq_residual, x0=[t_k, rho_seed], args=(d1, d2, z1, p_seed, s_target), method="hybr", tol=1.0e-10)
        verified0 = _verify_isentrope_solution(sol0, d1, d2, z1, p_seed, s_target)
        if verified0 is None:
            continue
        t_k, rho_seed, st0 = verified0
        points: List[Dict] = [{"T_K": t_k, "P_Pa": float(st0.p_pa), "h_Jmol": float(st0.h_jmol)}]

        p_start = float(st0.p_pa)
        if p_max_pa > p_start:
            for p_target in np.geomspace(p_start, p_max_pa, n_points)[1:]:
                sol = root(_isentrope_2eq_residual, x0=[t_k, rho_seed], args=(d1, d2, z1, p_target, s_target), method="hybr", tol=1.0e-10)
                verified = _verify_isentrope_solution(sol, d1, d2, z1, p_target, s_target)
                if verified is None:
                    break
                t_k, rho_seed, st = verified
                points.append({"T_K": t_k, "P_Pa": float(st.p_pa), "h_Jmol": float(st.h_jmol)})
        if len(points) > 1:
            ext[s_target] = points
    return ext


def compute_isentrope_vapor_side(d1: Dict, d2: Dict, z1: float, dew_rows: List[Dict],
                                  s_values_jmolK: List[float], p_min_pa: float = ISENTROPE_P_MIN_PA,
                                  n_points: int = ISENTROPE_N_POINTS, crit_point: Optional[Dict] = None,
                                  bubble_rows: Optional[List[Dict]] = None) -> Dict[float, List[Dict]]:
    """Reproduces `compute_isentrope_vapor_side` exactly: dew-row anchor
    selection (or critical-point fallback), incremental entropy ramp from
    anchor to target (bridges large gaps in small steps), dome-reentry
    guard (rejects any walked point that lands inside the two-phase
    envelope per bubble/dew interpolation), and adaptive step-bisection
    pressure walk both up and down from the seed."""
    converged_dew = [r for r in dew_rows if r["status"] == "CONVERGED"]

    dome_check = None
    if bubble_rows:
        conv_bubble = sorted([r for r in bubble_rows if r["status"] == "CONVERGED"], key=lambda r: r["P_Pa"])
        conv_dew_sorted = sorted(converged_dew, key=lambda r: r["P_Pa"])
        if conv_bubble and conv_dew_sorted:
            bubble_P_arr = np.array([r["P_Pa"] for r in conv_bubble])
            bubble_h_arr = np.array([mix_state(d1, d2, r["T_K"], r["rho_l_molm3"], z1).h_jmol for r in conv_bubble])
            dew_P_arr = np.array([r["P_Pa"] for r in conv_dew_sorted])
            dew_h_arr = np.array([mix_state(d1, d2, r["T_K"], r["rho_v_molm3"], z1).h_jmol for r in conv_dew_sorted])

            def _inside_dome(p_pa: float, h_jmol: float) -> bool:
                if not (bubble_P_arr.min() <= p_pa <= bubble_P_arr.max()):
                    return False
                if not (dew_P_arr.min() <= p_pa <= dew_P_arr.max()):
                    return False
                h_bub = float(np.interp(p_pa, bubble_P_arr, bubble_h_arr))
                h_dew = float(np.interp(p_pa, dew_P_arr, dew_h_arr))
                return h_bub < h_jmol < h_dew

            dome_check = _inside_dome

    def _adaptive_pressure_walk(t0: float, rho0: float, p_from: float, p_to: float,
                                 s_target_walk: float, n_points_walk: int) -> List[Dict]:
        """Reproduces `_adaptive_pressure_walk` exactly (adaptive log-pressure
        step halving/doubling continuation)."""
        points: List[Dict] = []
        log_from, log_to = np.log(p_from), np.log(p_to)
        total = log_to - log_from
        if total == 0:
            return points
        nominal_step = total / max(1, n_points_walk - 1)
        min_step = nominal_step / 1024.0
        step = nominal_step
        t_cur, rho_cur = t0, rho0
        log_p_cur = log_from
        while True:
            log_p_next = log_p_cur + step
            overshoot = (step > 0 and log_p_next >= log_to) or (step < 0 and log_p_next <= log_to)
            if overshoot:
                log_p_next = log_to
            if log_p_next == log_p_cur:
                break
            p_next = float(np.exp(log_p_next))
            sol = root(_isentrope_2eq_residual, x0=[t_cur, rho_cur], args=(d1, d2, z1, p_next, s_target_walk), method="hybr", tol=1.0e-10)
            verified = _verify_isentrope_solution(sol, d1, d2, z1, p_next, s_target_walk)
            ok = verified is not None
            st = None
            if ok:
                t_new, rho_new, st = verified
                if dome_check is not None and dome_check(float(st.p_pa), float(st.h_jmol)):
                    ok = False
            if ok:
                t_cur, rho_cur = t_new, rho_new
                log_p_cur = log_p_next
                points.append({"T_K": t_cur, "P_Pa": float(st.p_pa), "h_Jmol": float(st.h_jmol)})
                if log_p_cur == log_to:
                    break
                grown = step * 2.0
                step = grown if abs(grown) <= abs(nominal_step) else (abs(nominal_step) if step > 0 else -abs(nominal_step))
            else:
                step = step / 2.0
                if abs(step) < abs(min_step):
                    break
        return points

    have_crit = bool(crit_point is not None and crit_point.get("converged", False))
    ext: Dict[float, List[Dict]] = {s: [] for s in s_values_jmolK}
    for s_target in s_values_jmolK:
        best_row = None
        best_gap = None
        for r in converged_dew:
            t_k = r["T_K"]
            s_v = mix_entropy_direct(d1, d2, t_k, r["rho_v_molm3"], z1)
            if s_v >= s_target:
                continue
            gap = s_target - s_v
            if best_gap is None or gap < best_gap:
                best_gap = gap
                best_row = r

        if best_row is not None:
            t_k = best_row["T_K"]
            rho_seed = best_row["rho_v_molm3"]
            p_seed = best_row["P_Pa"]
            s_anchor = mix_entropy_direct(d1, d2, t_k, rho_seed, z1)
        elif have_crit:
            t_crit = crit_point["T_K"]
            rho_crit = crit_point["rho_molm3"]
            x1_crit = crit_point.get("x1", z1)
            s_crit = mix_entropy_direct(d1, d2, t_crit, rho_crit, x1_crit)
            if s_crit >= s_target:
                continue
            t_k = t_crit
            rho_seed = rho_crit
            p_seed = crit_point["P_Pa"]
            s_anchor = s_crit
        else:
            continue

        gap = s_target - s_anchor
        n_substeps = max(1, int(np.ceil(abs(gap) / ISENTROPE_VAPOR_MAX_ENTROPY_STEP_JMOLK)))
        s_ramp = np.linspace(s_anchor, s_target, n_substeps + 1)[1:]

        ramp_ok = True
        for s_step in s_ramp:
            sol_step = root(_isentrope_2eq_residual, x0=[t_k, rho_seed], args=(d1, d2, z1, p_seed, s_step), method="hybr", tol=1.0e-10)
            verified_step = _verify_isentrope_solution(sol_step, d1, d2, z1, p_seed, s_step)
            if verified_step is None:
                ramp_ok = False
                break
            t_k, rho_seed, _ = verified_step
        if not ramp_ok:
            continue

        st0 = mix_state(d1, d2, t_k, rho_seed, z1)
        if dome_check is not None and dome_check(float(st0.p_pa), float(st0.h_jmol)):
            continue
        seed_point = {"T_K": t_k, "P_Pa": float(st0.p_pa), "h_Jmol": float(st0.h_jmol)}
        p_start = float(st0.p_pa)

        down_points: List[Dict] = []
        if p_min_pa < p_start:
            down_points = _adaptive_pressure_walk(t_k, rho_seed, p_start, p_min_pa, s_target, n_points)

        up_points: List[Dict] = []
        if p_start < ISENTROPE_P_MAX_PA:
            up_points = _adaptive_pressure_walk(t_k, rho_seed, p_start, ISENTROPE_P_MAX_PA, s_target, n_points)

        points = list(reversed(up_points)) + [seed_point] + down_points
        if len(points) > 1:
            ext[s_target] = points
    return ext


# =============================================================================
# Stage L (part 3) design support: pseudo-pure (T,H)/(P,H) reference flash.
# NOT a port of anything in the oracle (the oracle has no pseudo-pure mode
# at all -- see the module-level comment above `solve_pseudopure_
# saturation_at_t`). This is explicit-branching REFERENCE code (perfectly
# legitimate per spec rule 44's carve-out: this function is for
# `initialize()` and as ground truth to validate the native-Pyomo smooth
# complementarity/blending logic against -- it is never called from inside
# an active StateBlockData Constraint).
# =============================================================================
def solve_t_sat_at_p_pseudopure(d1: Dict, d2: Dict, z1: float, p_pa: float,
                                 t_guess: float, rho_l0: float, rho_v0: float,
                                 t_bracket_half_width_k: float = 40.0) -> Dict:
    """
    Inverts `solve_pseudopure_saturation_at_t` (which is naturally a
    function OF temperature, giving pressure as an output) to instead find
    the saturation TEMPERATURE at a given PRESSURE -- via `brentq` bisection
    on T, since `solve_pseudopure_saturation_at_t(T)['P_Pa']` is smooth and
    monotonically increasing in T over the practical range (confirmed by
    the clean, monotonic 255-375K continuation sweep already validated in
    `validate_pyomo_saturation_vs_core.py`). Brackets T around `t_guess`
    using seeds carried forward at each trial T (continuation, per the
    seed-sensitivity lesson recorded in helmholtz_prop_validation.md
    Section 22 -- avoids the same non-monotonic-root failure mode found
    there).
    """
    rl, rv = rho_l0, rho_v0

    def p_of_t(t_k: float) -> float:
        nonlocal rl, rv
        sp = solve_pseudopure_saturation_at_t(d1, d2, z1, float(t_k), rl, rv)
        if sp["status"] == "CONVERGED":
            rl, rv = sp["rho_l_molm3"], sp["rho_v_molm3"]
        return float(sp["P_Pa"]) - p_pa

    t_lo = max(180.0, t_guess - t_bracket_half_width_k)
    t_hi = t_guess + t_bracket_half_width_k
    f_lo, f_hi = p_of_t(t_lo), p_of_t(t_hi)
    if f_lo * f_hi > 0.0:
        # Bracket failed (e.g. t_guess far off) -- widen once before giving up.
        t_lo = max(180.0, t_guess - 2.0 * t_bracket_half_width_k)
        t_hi = t_guess + 2.0 * t_bracket_half_width_k
        f_lo, f_hi = p_of_t(t_lo), p_of_t(t_hi)
        if f_lo * f_hi > 0.0:
            return {"status": "DIVERGED", "T_K": float(t_guess), "notes": "no sign change bracketing T_sat(P)"}
    t_sat = brentq(p_of_t, t_lo, t_hi, xtol=1.0e-8, rtol=1.0e-12, maxiter=100)
    sp = solve_pseudopure_saturation_at_t(d1, d2, z1, float(t_sat), rl, rv)
    sp["status"] = "CONVERGED" if sp["status"] == "CONVERGED" and abs(sp["P_Pa"] - p_pa) / max(1.0, p_pa) < 1e-6 else "DIVERGED"
    return sp


def solve_singlephase_tr_at_ph_pseudopure(d1: Dict, d2: Dict, z1: float, p_pa: float, h_target_jmol: float,
                                           t_guess: float, rho_guess: float) -> Dict:
    """
    Single-phase (T, rho) inversion at fixed pseudo-pure composition: find
    (T, rho) such that P(T,rho,z1)=p_pa and h(T,rho,z1)=h_target_jmol.
    Works for EITHER the liquid-like or vapor-like branch depending on
    `rho_guess` -- there is no explicit branching inside this function,
    just a straightforward 2-equation Newton-type solve via `scipy.optimize.
    root`; the CALLER (e.g. `flash_ph_pseudopure` below) is responsible for
    picking a physically appropriate seed/branch.
    """
    def res(u):
        t_k, rho = float(u[0]), max(float(u[1]), RHO_MIN_MOLM3)
        st = mix_state(d1, d2, t_k, rho, z1)
        r_p = (st.p_pa - p_pa) / max(1.0, abs(p_pa))
        r_h = (st.h_jmol - h_target_jmol) / max(1.0, abs(h_target_jmol))
        return [r_p, r_h]

    sol = root(res, x0=[float(t_guess), float(rho_guess)], method="hybr", tol=1.0e-12)
    t_k, rho = float(sol.x[0]), max(float(sol.x[1]), RHO_MIN_MOLM3)
    st = mix_state(d1, d2, t_k, rho, z1)
    r_p = abs(st.p_pa - p_pa) / max(1.0, abs(p_pa))
    r_h = abs(st.h_jmol - h_target_jmol) / max(1.0, abs(h_target_jmol))
    ok = bool(sol.success and r_p <= 1e-8 and r_h <= 1e-8 and rho > 0.0)
    return {"status": "CONVERGED" if ok else "DIVERGED", "T_K": t_k, "rho_molm3": float(rho),
            "P_Pa": float(st.p_pa), "h_Jmol": float(st.h_jmol), "r_P": r_p, "r_h": r_h,
            "notes": str(sol.message)}


def flash_ph_pseudopure(d1: Dict, d2: Dict, z1: float, p_pa: float, h_target_jmol: float,
                         t_sat_guess: float, rho_l_guess: float, rho_v_guess: float) -> Dict:
    """
    Full pseudo-pure (x1=y1=z1 fixed) PH flash: given pressure and molar
    enthalpy, determine whether the state is subcooled liquid, two-phase,
    or superheated vapor, and return (T, rho or (rho_l,rho_v,vapor
    quality), P, h) accordingly. EXPLICIT branching (perfectly legitimate
    here -- see module-level comment above; this is reference/
    initialize()-time code, never called from inside an active
    StateBlockData Constraint). Serves as the ground-truth reference for
    validating the native-Pyomo smooth complementarity/blending Constraint
    set (Stage L part 3, next).

    Vapor quality convention: MOLAR quality (mol vapor / mol total) --
    since composition is fixed identical in both phases here, molar and
    mass quality coincide numerically (same molecular weight on both
    sides), so no separate mass/molar quality distinction is needed,
    unlike the oracle's real multi-composition quality-line machinery
    (Stage K, Section 11) where they can differ.
    """
    sat = solve_t_sat_at_p_pseudopure(d1, d2, z1, p_pa, t_sat_guess, rho_l_guess, rho_v_guess)
    if sat["status"] != "CONVERGED":
        return {"status": "DIVERGED", "region": "unknown", "notes": "saturation-curve solve failed at this P"}

    t_sat = sat["T_K"]
    rho_l_sat = sat["rho_l_molm3"]
    rho_v_sat = sat["rho_v_molm3"]
    h_l_sat = sat["h_l_Jmol"]
    h_v_sat = sat["h_v_Jmol"]

    if h_target_jmol <= h_l_sat:
        # Subcooled liquid: seed the single-phase solve from the saturated-liquid state.
        sp = solve_singlephase_tr_at_ph_pseudopure(d1, d2, z1, p_pa, h_target_jmol, t_sat, rho_l_sat)
        if sp["status"] != "CONVERGED":
            return {"status": "DIVERGED", "region": "liquid", "notes": sp["notes"]}
        return {"status": "CONVERGED", "region": "subcooled_liquid", "T_K": sp["T_K"], "rho_molm3": sp["rho_molm3"],
                "P_Pa": sp["P_Pa"], "h_Jmol": sp["h_Jmol"], "vapor_frac": 0.0,
                "T_sat_K": t_sat, "rho_l_sat_molm3": rho_l_sat, "rho_v_sat_molm3": rho_v_sat}
    elif h_target_jmol >= h_v_sat:
        # Superheated vapor: seed the single-phase solve from the saturated-vapor state.
        sp = solve_singlephase_tr_at_ph_pseudopure(d1, d2, z1, p_pa, h_target_jmol, t_sat, rho_v_sat)
        if sp["status"] != "CONVERGED":
            return {"status": "DIVERGED", "region": "vapor", "notes": sp["notes"]}
        return {"status": "CONVERGED", "region": "superheated_vapor", "T_K": sp["T_K"], "rho_molm3": sp["rho_molm3"],
                "P_Pa": sp["P_Pa"], "h_Jmol": sp["h_Jmol"], "vapor_frac": 1.0,
                "T_sat_K": t_sat, "rho_l_sat_molm3": rho_l_sat, "rho_v_sat_molm3": rho_v_sat}
    else:
        # Two-phase: lever rule on molar enthalpy for quality; molar volume
        # lever rule for the mixture's overall molar density.
        vf = (h_target_jmol - h_l_sat) / max(1.0e-12, (h_v_sat - h_l_sat))
        v_mix = (1.0 - vf) * (1.0 / rho_l_sat) + vf * (1.0 / rho_v_sat)
        rho_mix = 1.0 / v_mix
        return {"status": "CONVERGED", "region": "two_phase", "T_K": t_sat, "rho_molm3": float(rho_mix),
                "P_Pa": p_pa, "h_Jmol": float(h_target_jmol), "vapor_frac": float(vf),
                "T_sat_K": t_sat, "rho_l_sat_molm3": rho_l_sat, "rho_v_sat_molm3": rho_v_sat}
