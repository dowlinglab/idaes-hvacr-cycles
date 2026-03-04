# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Finite-UA lumped HX models with standardized cycle interface.

Context
-------
Variant B (`cycle_dx_hx_lumpedUA`): evaporator and condenser each use a
single-zone phase-change-dominant epsilon-NTU formulation where UA directly
determines heat transfer; no SH/SC outlet targets are enforced inside HX.

This implementation adds physical bounds so UA-driven duties remain
thermodynamically consistent.
"""

from typing import Dict, Tuple

import CoolProp.CoolProp as CP

from .config import CycleConfig
from .types import AirStream, RefrigerantState


def _eps_cr0(ntu: float) -> float:
    """Counterflow effectiveness for ``Cr -> 0`` phase-change approximation."""
    if ntu <= 0.0:
        return 0.0
    return 1.0 - pow(2.718281828459045, -ntu)


def _h_at_p_t(p: float, t: float, fluid: str) -> float:
    """Safe pressure-temperature flash helper."""
    return CP.PropsSI("H", "P", p, "T", t, fluid)


def solve_evaporator(
    state_in: RefrigerantState,
    air_in: AirStream,
    cfg: CycleConfig,
    refrigerant: str,
) -> Tuple[RefrigerantState, float, Dict[str, float]]:
    """Compute evaporator duty from UA and driving temperature only."""
    c_air = max(air_in.m_dot * air_in.cp, 1.0e-9)
    t_sat = CP.PropsSI("T", "P", state_in.p, "Q", 1, refrigerant)

    ntu = cfg.UA_evap_total / c_air
    eps = _eps_cr0(ntu)

    q_raw = max(eps * c_air * (air_in.t_in - t_sat), 0.0)
    h_raw = state_in.h + q_raw / state_in.m_dot

    # Physical cap: refrigerant outlet cannot exceed inlet-air temperature.
    t_ref_out_max = max(t_sat, air_in.t_in - 0.5)
    h_max = _h_at_p_t(state_in.p, t_ref_out_max, refrigerant)

    h_out = min(max(h_raw, state_in.h), h_max)
    q = state_in.m_dot * (h_out - state_in.h)
    t_air_out = air_in.t_in - q / c_air

    state_out = RefrigerantState(p=state_in.p, h=h_out, m_dot=state_in.m_dot)
    aux = {
        "UA_used": cfg.UA_evap_total,
        "NTU": ntu,
        "epsilon": eps,
        "t_sat": t_sat,
        "t_air_out": t_air_out,
        "q_total": q,
        "q_raw": q_raw,
        "limited": 1.0 if abs(q - q_raw) > 1.0e-9 else 0.0,
    }
    return state_out, q, aux


def solve_condenser(
    state_in: RefrigerantState,
    air_in: AirStream,
    cfg: CycleConfig,
    refrigerant: str,
) -> Tuple[RefrigerantState, float, Dict[str, float]]:
    """Compute condenser duty from UA and driving temperature only."""
    c_air = max(air_in.m_dot * air_in.cp, 1.0e-9)
    t_sat = CP.PropsSI("T", "P", state_in.p, "Q", 0, refrigerant)

    ntu = cfg.UA_cond_total / c_air
    eps = _eps_cr0(ntu)

    q_raw = max(eps * c_air * (t_sat - air_in.t_in), 0.0)
    h_raw = state_in.h - q_raw / state_in.m_dot

    # Physical cap: refrigerant outlet cannot cool below entering air temperature.
    t_ref_out_min = min(t_sat - 1.0e-3, air_in.t_in + 0.5)
    h_min = _h_at_p_t(state_in.p, t_ref_out_min, refrigerant)

    h_out = max(min(h_raw, state_in.h), h_min)
    q = state_in.m_dot * (state_in.h - h_out)
    t_air_out = air_in.t_in + q / c_air

    state_out = RefrigerantState(p=state_in.p, h=h_out, m_dot=state_in.m_dot)
    aux = {
        "UA_used": cfg.UA_cond_total,
        "NTU": ntu,
        "epsilon": eps,
        "t_sat": t_sat,
        "t_air_out": t_air_out,
        "q_total": q,
        "q_raw": q_raw,
        "limited": 1.0 if abs(q - q_raw) > 1.0e-9 else 0.0,
    }
    return state_out, q, aux
