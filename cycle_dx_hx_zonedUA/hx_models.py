# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""DX zoned finite-UA HX models with standardized cycle interface.

Context
-------
Variant C (`cycle_dx_hx_zonedUA`): sequential evaporator 2-zone and condenser
3-zone models where each zone duty is UA-driven only (no target-driven Q_req),
with phase-boundary and air-approach bounds for thermodynamic consistency.
"""

from typing import Dict, Tuple

import CoolProp.CoolProp as CP

from .config import CycleConfig
from .types import AirStream, RefrigerantState


def _eps_counterflow(ntu: float, c_r: float) -> float:
    """Counterflow epsilon-NTU effectiveness for ``0 <= c_r < 1``."""
    if ntu <= 0.0:
        return 0.0
    c_r = max(min(c_r, 0.999999), 0.0)
    if abs(1.0 - c_r) < 1.0e-8:
        return ntu / (1.0 + ntu)
    exp_term = pow(2.718281828459045, -ntu * (1.0 - c_r))
    return (1.0 - exp_term) / (1.0 - c_r * exp_term)


def _q_single_phase_ua(ua: float, c_air: float, c_ref: float, dt_hot_cold: float) -> Tuple[float, float, float]:
    """Return UA-driven heat transfer and epsilon-NTU details for one zone."""
    c_min = max(min(c_air, c_ref), 1.0e-9)
    c_max = max(c_air, c_ref)
    c_r = c_min / max(c_max, 1.0e-9)
    ntu = ua / c_min
    eps = _eps_counterflow(ntu, c_r)
    q = eps * c_min * max(dt_hot_cold, 0.0)
    return q, ntu, eps


def solve_evaporator(
    state_in: RefrigerantState,
    air_in: AirStream,
    cfg: CycleConfig,
    refrigerant: str,
) -> Tuple[RefrigerantState, float, Dict[str, float]]:
    """Solve 2-zone evaporator sequentially with UA-driven zone duties."""
    c_air = max(air_in.m_dot * air_in.cp, 1.0e-9)
    t_sat = CP.PropsSI("T", "P", state_in.p, "Q", 1, refrigerant)
    h_g = CP.PropsSI("H", "P", state_in.p, "Q", 1, refrigerant)

    ntu_tp = cfg.UA_evap_tp / c_air
    eps_tp = 1.0 - pow(2.718281828459045, -max(ntu_tp, 0.0))
    q_tp_raw = eps_tp * c_air * max(air_in.t_in - t_sat, 0.0)
    q_tp_cap = max(state_in.m_dot * (h_g - state_in.h), 0.0)
    q_tp = min(q_tp_raw, q_tp_cap)

    h_after_tp = state_in.h + q_tp / state_in.m_dot
    t_air_after_tp = air_in.t_in - q_tp / c_air

    if h_after_tp >= h_g - 1.0e-6:
        t_ref_in_sh = CP.PropsSI("T", "P", state_in.p, "H", max(h_after_tp, h_g), refrigerant)
        cp_ref_sh = CP.PropsSI("C", "P", state_in.p, "T", max(t_ref_in_sh, t_sat + 1.0e-3), refrigerant)
        c_ref_sh = max(state_in.m_dot * cp_ref_sh, 1.0e-9)
        q_sh_raw, ntu_sh, eps_sh = _q_single_phase_ua(
            cfg.UA_evap_sh, c_air, c_ref_sh, t_air_after_tp - t_ref_in_sh
        )
        t_ref_out_max = max(t_ref_in_sh, t_air_after_tp - 0.5)
        h_max = CP.PropsSI("H", "P", state_in.p, "T", t_ref_out_max, refrigerant)
        q_sh_cap = max(state_in.m_dot * (h_max - h_after_tp), 0.0)
        q_sh = min(q_sh_raw, q_sh_cap)
    else:
        q_sh = 0.0
        ntu_sh = 0.0
        eps_sh = 0.0

    h_out = h_after_tp + q_sh / state_in.m_dot
    t_air_out = t_air_after_tp - q_sh / c_air
    q_total = q_tp + q_sh

    state_out = RefrigerantState(p=state_in.p, h=h_out, m_dot=state_in.m_dot)
    aux = {
        "UA_used": cfg.UA_evap_tp + cfg.UA_evap_sh,
        "q_total": q_total,
        "zone_q_tp": q_tp,
        "zone_q_sh": q_sh,
        "zone_q_tp_raw": q_tp_raw,
        "zone_q_sh_raw": q_sh if ntu_sh == 0 else q_sh_raw,
        "NTU_tp": ntu_tp,
        "NTU_sh": ntu_sh,
        "eps_tp": eps_tp,
        "eps_sh": eps_sh,
        "t_air_out": t_air_out,
        "t_sat": t_sat,
    }
    return state_out, q_total, aux


def solve_condenser(
    state_in: RefrigerantState,
    air_in: AirStream,
    cfg: CycleConfig,
    refrigerant: str,
) -> Tuple[RefrigerantState, float, Dict[str, float]]:
    """Solve 3-zone condenser sequentially with UA-driven zone duties."""
    c_air = max(air_in.m_dot * air_in.cp, 1.0e-9)
    t_sat = CP.PropsSI("T", "P", state_in.p, "Q", 0, refrigerant)
    h_v = CP.PropsSI("H", "P", state_in.p, "Q", 1, refrigerant)
    h_l = CP.PropsSI("H", "P", state_in.p, "Q", 0, refrigerant)

    t_ref_in_ds = CP.PropsSI("T", "P", state_in.p, "H", state_in.h, refrigerant)
    cp_ref_ds = CP.PropsSI("C", "P", state_in.p, "T", max(t_ref_in_ds, t_sat + 1.0e-3), refrigerant)
    c_ref_ds = max(state_in.m_dot * cp_ref_ds, 1.0e-9)
    q_ds_raw, ntu_ds, eps_ds = _q_single_phase_ua(cfg.UA_cond_ds, c_air, c_ref_ds, t_ref_in_ds - air_in.t_in)
    q_ds_cap = max(state_in.m_dot * (state_in.h - h_v), 0.0)
    q_ds = min(q_ds_raw, q_ds_cap)

    h_after_ds = state_in.h - q_ds / state_in.m_dot
    t_air_after_ds = air_in.t_in + q_ds / c_air

    ntu_tp = cfg.UA_cond_tp / c_air
    eps_tp = 1.0 - pow(2.718281828459045, -max(ntu_tp, 0.0))
    q_tp_raw = eps_tp * c_air * max(t_sat - t_air_after_ds, 0.0)
    q_tp_cap = max(state_in.m_dot * (h_after_ds - h_l), 0.0) if h_after_ds <= h_v + 1.0e-6 else 0.0
    q_tp = min(q_tp_raw, q_tp_cap)

    h_after_tp = h_after_ds - q_tp / state_in.m_dot
    t_air_after_tp = t_air_after_ds + q_tp / c_air

    if h_after_tp <= h_l + 1.0e-6:
        t_ref_in_sc = CP.PropsSI("T", "P", state_in.p, "H", min(h_after_tp, h_l), refrigerant)
        t_cp_sc = min(max(t_ref_in_sc, 200.0), t_sat - 1.0e-3)
        cp_ref_sc = CP.PropsSI("C", "P", state_in.p, "T", t_cp_sc, refrigerant)
        c_ref_sc = max(state_in.m_dot * cp_ref_sc, 1.0e-9)
        q_sc_raw, ntu_sc, eps_sc = _q_single_phase_ua(
            cfg.UA_cond_sc, c_air, c_ref_sc, t_ref_in_sc - t_air_after_tp
        )
        t_ref_out_min = min(t_ref_in_sc, t_air_after_tp + 0.5)
        h_min = CP.PropsSI("H", "P", state_in.p, "T", t_ref_out_min, refrigerant)
        q_sc_cap = max(state_in.m_dot * (h_after_tp - h_min), 0.0)
        q_sc = min(q_sc_raw, q_sc_cap)
    else:
        q_sc = 0.0
        ntu_sc = 0.0
        eps_sc = 0.0

    h_out = h_after_tp - q_sc / state_in.m_dot
    t_air_out = t_air_after_tp + q_sc / c_air
    q_total = q_ds + q_tp + q_sc

    state_out = RefrigerantState(p=state_in.p, h=h_out, m_dot=state_in.m_dot)
    aux = {
        "UA_used": cfg.UA_cond_ds + cfg.UA_cond_tp + cfg.UA_cond_sc,
        "q_total": q_total,
        "zone_q_ds": q_ds,
        "zone_q_tp": q_tp,
        "zone_q_sc": q_sc,
        "zone_q_ds_raw": q_ds_raw,
        "zone_q_tp_raw": q_tp_raw,
        "zone_q_sc_raw": q_sc if ntu_sc == 0 else q_sc_raw,
        "NTU_ds": ntu_ds,
        "NTU_tp": ntu_tp,
        "NTU_sc": ntu_sc,
        "eps_ds": eps_ds,
        "eps_tp": eps_tp,
        "eps_sc": eps_sc,
        "t_air_out": t_air_out,
        "t_sat": t_sat,
    }
    return state_out, q_total, aux
