# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Cycle solver with standardized HX call sites.

Context
-------
This module is the solver scaffold used by all variants. In this directory,
its behavior differs from the project baseline only by delegating HX work to
`hx_models.py`, which here implements ideal HX behavior.
"""

from typing import Dict, Tuple

import CoolProp.CoolProp as CP

from .config import CycleConfig
from .hx_models import solve_condenser, solve_evaporator
from .types import AirStream, CycleResult, RefrigerantState


def _state_from_p_t(p: float, t: float, fluid: str) -> Tuple[float, float]:
    """Return ``(h, s)`` at pressure ``p`` and temperature ``t``."""
    h = CP.PropsSI("H", "P", p, "T", t, fluid)
    s = CP.PropsSI("S", "P", p, "T", t, fluid)
    return h, s


def _compute_plf(plr: float, cd: float) -> float:
    """Return part-load fraction using the legacy PLR/CD correlation."""
    plf = 1.0 - cd * (1.0 - plr)
    return max(0.0, min(1.0, plf))


def solve_cycle_point(
    fluid: str,
    t_evap_sat_c: float,
    t_cond_sat_c: float,
    cfg: CycleConfig,
    t_air_evap_in_c: float,
    t_air_cond_in_c: float,
) -> CycleResult:
    """Solve one vapor-compression cycle point.

    The compressor and expansion models are intentionally identical across all
    variants; only evaporator/condenser behavior is swapped through
    ``solve_evaporator`` and ``solve_condenser``.
    """

    t_evap_sat_k = t_evap_sat_c + 273.15
    t_cond_sat_k = t_cond_sat_c + 273.15

    p_evap = CP.PropsSI("P", "T", t_evap_sat_k, "Q", 1, fluid)
    p_cond = CP.PropsSI("P", "T", t_cond_sat_k, "Q", 0, fluid)

    t1 = t_evap_sat_k + cfg.superheat_target_K
    h1, s1 = _state_from_p_t(p_evap, t1, fluid)

    h2s = CP.PropsSI("H", "P", p_cond, "S", s1, fluid)
    h2 = h1 + (h2s - h1) / cfg.eta_isentropic
    state2 = RefrigerantState(p=p_cond, h=h2, m_dot=cfg.m_dot_ref)

    air_cond = AirStream(
        t_in=t_air_cond_in_c + 273.15,
        m_dot=cfg.m_dot_air_cond,
        cp=cfg.cp_air_cond,
    )
    state3, q_cond, aux_cond = solve_condenser(state2, air_cond, cfg, fluid)

    state4 = RefrigerantState(p=p_evap, h=state3.h, m_dot=cfg.m_dot_ref)

    air_evap = AirStream(
        t_in=t_air_evap_in_c + 273.15,
        m_dot=cfg.m_dot_air_evap,
        cp=cfg.cp_air_evap,
    )
    state1_out, q_evap, aux_evap = solve_evaporator(state4, air_evap, cfg, fluid)

    q_evap_full = q_evap
    q_cond_full = q_cond
    w_comp_full = cfg.m_dot_ref * (h2 - h1)
    cop_full = q_evap_full / w_comp_full if w_comp_full > 0 else 0.0
    plf = _compute_plf(cfg.plr, cfg.cd)
    cop_part = plf * cop_full
    q_evap_part = cfg.plr * q_evap_full
    q_cond_part = cfg.plr * q_cond_full
    w_comp_part = q_evap_part / cop_part if cop_part > 0 else 0.0

    diagnostics: Dict[str, float] = {
        "h1_in": h1,
        "h1_out": state1_out.h,
        "h2": h2,
        "h3": state3.h,
        "h4": state4.h,
        "p_ratio": p_cond / p_evap,
        "PLR": cfg.plr,
        "CD": cfg.cd,
        "PLF": plf,
        "COP_full": cop_full,
        "COP_part": cop_part,
        "m_dot_effective": cfg.m_dot_ref * cfg.plr,
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
