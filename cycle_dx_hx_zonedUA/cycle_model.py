# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Cycle solver with standardized HX call sites.

Context
-------
Variant C (`cycle_dx_hx_zonedUA`): compressor/expansion models are unchanged,
while HX duties come from zoned finite UA. Superheat/subcooling are computed outputs.
"""

from typing import Dict, Tuple

import CoolProp.CoolProp as CP

from .config import CycleConfig
from .hx_models import solve_condenser, solve_evaporator
from .types import AirStream, CycleResult, RefrigerantState


def _compute_plf(plr: float, cd: float) -> float:
    """Return part-load fraction using the legacy PLR/CD correlation."""
    plf = 1.0 - cd * (1.0 - plr)
    return max(0.0, min(1.0, plf))


def _compute_evap_outputs(p: float, h: float, fluid: str) -> Dict[str, float]:
    """Compute actual evaporator outlet SH and quality from ``(P, h)``."""
    h_f = CP.PropsSI("H", "P", p, "Q", 0, fluid)
    h_g = CP.PropsSI("H", "P", p, "Q", 1, fluid)
    t_sat = CP.PropsSI("T", "P", p, "Q", 1, fluid)

    if h_f < h < h_g:
        x = (h - h_f) / max(h_g - h_f, 1.0e-9)
        sh = 0.0
    elif h >= h_g:
        try:
            t_out = CP.PropsSI("T", "P", p, "H", h, fluid)
            sh = max(0.0, t_out - t_sat)
        except ValueError:
            sh = float("nan")
        x = float("nan")
    else:
        x = 0.0
        sh = 0.0

    return {"SH_actual": sh, "x_evap_out": x}


def _compute_cond_outputs(p: float, h: float, fluid: str) -> Dict[str, float]:
    """Compute actual condenser outlet SC and quality from ``(P, h)``."""
    h_f = CP.PropsSI("H", "P", p, "Q", 0, fluid)
    h_g = CP.PropsSI("H", "P", p, "Q", 1, fluid)
    t_sat = CP.PropsSI("T", "P", p, "Q", 0, fluid)

    if h_f < h < h_g:
        x = (h - h_f) / max(h_g - h_f, 1.0e-9)
        sc = 0.0
    elif h <= h_f:
        try:
            t_out = CP.PropsSI("T", "P", p, "H", h, fluid)
            sc = max(0.0, t_sat - t_out)
        except ValueError:
            sc = float("nan")
        x = float("nan")
    else:
        x = 1.0
        sc = 0.0

    return {"SC_actual": sc, "x_cond_out": x}


def _compressor_step(h1: float, p_evap: float, p_cond: float, eta_isen: float, fluid: str) -> Tuple[float, float, float]:
    """Return compressor outlet enthalpy, inlet entropy, and effective inlet h.

    If the evaporator outlet is still in two-phase region, use saturated vapor
    at evaporator pressure as a surrogate compressor inlet state for stability.
    """
    h_g = CP.PropsSI("H", "P", p_evap, "Q", 1, fluid)
    h1_eff = max(h1, h_g + 1.0e-3)
    s1 = CP.PropsSI("S", "P", p_evap, "H", h1_eff, fluid)
    try:
        h2s = CP.PropsSI("H", "P", p_cond, "S", s1, fluid)
    except ValueError:
        h2s = CP.PropsSI("H", "P", p_cond, "Q", 1, fluid)
    h2 = h1_eff + (h2s - h1_eff) / eta_isen
    return h2, s1, h1_eff


def solve_cycle_point(
    fluid: str,
    t_evap_sat_c: float,
    t_cond_sat_c: float,
    cfg: CycleConfig,
    t_air_evap_in_c: float,
    t_air_cond_in_c: float,
) -> CycleResult:
    """Solve one cycle point using fixed-point closure on compressor inlet enthalpy."""
    t_evap_sat_k = t_evap_sat_c + 273.15
    t_cond_sat_k = t_cond_sat_c + 273.15

    p_evap = CP.PropsSI("P", "T", t_evap_sat_k, "Q", 1, fluid)
    p_cond = CP.PropsSI("P", "T", t_cond_sat_k, "Q", 0, fluid)

    h1 = CP.PropsSI("H", "P", p_evap, "Q", 1, fluid)
    state1_out = RefrigerantState(p=p_evap, h=h1, m_dot=cfg.m_dot_ref)
    state3 = RefrigerantState(p=p_cond, h=CP.PropsSI("H", "P", p_cond, "Q", 0, fluid), m_dot=cfg.m_dot_ref)
    h4_value = state3.h
    q_evap = 0.0
    q_cond = 0.0
    aux_evap: Dict[str, float] = {}
    aux_cond: Dict[str, float] = {}

    air_cond = AirStream(t_in=t_air_cond_in_c + 273.15, m_dot=cfg.m_dot_air_cond, cp=cfg.cp_air_cond)
    air_evap = AirStream(t_in=t_air_evap_in_c + 273.15, m_dot=cfg.m_dot_air_evap, cp=cfg.cp_air_evap)

    for _ in range(25):
        h2, s1, h1_comp = _compressor_step(h1, p_evap, p_cond, cfg.eta_isentropic, fluid)
        state2 = RefrigerantState(p=p_cond, h=h2, m_dot=cfg.m_dot_ref)

        state3, q_cond, aux_cond = solve_condenser(state2, air_cond, cfg, fluid)
        state4 = RefrigerantState(p=p_evap, h=state3.h, m_dot=cfg.m_dot_ref)
        h4_value = state4.h
        state1_out, q_evap, aux_evap = solve_evaporator(state4, air_evap, cfg, fluid)

        h1_new = state1_out.h
        if abs(h1_new - h1) < 1.0e-3:
            h1 = h1_new
            break
        h1 = 0.5 * h1 + 0.5 * h1_new

    h2, s1, h1_comp = _compressor_step(h1, p_evap, p_cond, cfg.eta_isentropic, fluid)
    q_evap_full = q_evap
    q_cond_full = q_cond
    w_comp_full = cfg.m_dot_ref * (h2 - h1_comp)
    cop_full = q_evap_full / w_comp_full if w_comp_full > 0 else 0.0
    plf = _compute_plf(cfg.plr, cfg.cd)
    cop_part = plf * cop_full
    q_evap_part = cfg.plr * q_evap_full
    q_cond_part = cfg.plr * q_cond_full
    w_comp_part = q_evap_part / cop_part if cop_part > 0 else 0.0

    evap_out = _compute_evap_outputs(p_evap, state1_out.h, fluid)
    cond_out = _compute_cond_outputs(p_cond, state3.h, fluid)

    diagnostics: Dict[str, float] = {
        "h1_in": h1,
        "h1_comp_effective": h1_comp,
        "h1_out": state1_out.h,
        "h2": h2,
        "h3": state3.h,
        "h4": h4_value,
        "s1": s1,
        "p_ratio": p_cond / p_evap,
        "PLR": cfg.plr,
        "CD": cfg.cd,
        "PLF": plf,
        "COP_full": cop_full,
        "COP_part": cop_part,
        "m_dot_effective": cfg.m_dot_ref * cfg.plr,
        **evap_out,
        **cond_out,
    }
    diagnostics.update({f"evap_{k}": float(v) for k, v in aux_evap.items()})
    diagnostics.update({f"cond_{k}": float(v) for k, v in aux_cond.items()})

    return CycleResult(
        fluid=fluid,
        t_evap_sat_c=float(t_evap_sat_c),
        t_cond_sat_c=float(t_cond_sat_c),
        p_evap_pa=float(p_evap),
        p_cond_pa=float(p_cond),
        m_dot_ref=float(cfg.m_dot_ref * cfg.plr),
        q_evap_w=float(q_evap_part),
        q_cond_w=float(q_cond_part),
        w_comp_w=float(w_comp_part),
        cop=float(cop_part),
        diagnostics=diagnostics,
    )
