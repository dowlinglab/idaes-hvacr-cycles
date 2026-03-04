# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Idealized HX models with common solver interface.

Context
-------
This directory is the HX-ideal variant. The evaporator and condenser are
modeled as ideal components meeting superheat/subcool targets without UA
limitation. This module exists to keep the same API used by finite-UA variants.
"""

from typing import Dict, Tuple

import CoolProp.CoolProp as CP

from .config import CycleConfig
from .types import AirStream, RefrigerantState


def solve_evaporator(
    state_in: RefrigerantState,
    air_in: AirStream,
    cfg: CycleConfig,
    refrigerant: str,
) -> Tuple[RefrigerantState, float, Dict[str, float]]:
    """Solve ideal evaporator duty and outlet state.

    Parameters follow the cross-variant standardized HX API.
    """
    t_sat = CP.PropsSI("T", "P", state_in.p, "Q", 1, refrigerant)
    t_out = t_sat + cfg.superheat_target_K
    h_out = CP.PropsSI("H", "P", state_in.p, "T", t_out, refrigerant)

    q = state_in.m_dot * (h_out - state_in.h)
    t_air_out = air_in.t_in - q / (air_in.m_dot * air_in.cp)

    state_out = RefrigerantState(p=state_in.p, h=h_out, m_dot=state_in.m_dot)
    aux = {
        "UA_used": float("inf"),
        "t_sat": t_sat,
        "t_ref_out": t_out,
        "t_air_out": t_air_out,
        "q_total": q,
        "limited": 0.0,
    }
    return state_out, q, aux


def solve_condenser(
    state_in: RefrigerantState,
    air_in: AirStream,
    cfg: CycleConfig,
    refrigerant: str,
) -> Tuple[RefrigerantState, float, Dict[str, float]]:
    """Solve ideal condenser duty and outlet state.

    Parameters follow the cross-variant standardized HX API.
    """
    t_sat = CP.PropsSI("T", "P", state_in.p, "Q", 0, refrigerant)
    t_out = t_sat - cfg.subcool_target_K
    h_out = CP.PropsSI("H", "P", state_in.p, "T", t_out, refrigerant)

    q = state_in.m_dot * (state_in.h - h_out)
    t_air_out = air_in.t_in + q / (air_in.m_dot * air_in.cp)

    state_out = RefrigerantState(p=state_in.p, h=h_out, m_dot=state_in.m_dot)
    aux = {
        "UA_used": float("inf"),
        "t_sat": t_sat,
        "t_ref_out": t_out,
        "t_air_out": t_air_out,
        "q_total": q,
        "limited": 0.0,
    }
    return state_out, q, aux
