"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/14/2026

What this does
    Reads the solved model and returns its heat balance and state points:
    mass flow, every pipe's T, P, h, and every unit's heat or work.

Why
    The COP alone hides how it was reached. Collaborators (Zhuorui) need the
    duties and states behind each COP, in one row per ambient.

Basis
    Cascade: CO2 (cold) loop flow fixed at 1 kg/s; hot-loop flow solved.
    Single stage: refrigerant flow fixed at 1 kg/s.
    So duties/works are in kW per 1 kg/s of the fixed loop. COP and all
    ratios (W/Q_box, Q_air/Q_box) do not depend on this basis.

Columns (prefix hot_ / cold_ for the cascade, none for single stage)
    m_kg_s                         loop mass flow
    <unit>_in/out_T_C, _P_kPa, _h_kJkg    for evap, comp, cond, valve
    Q_evap_kW, W_comp_kW, Q_cond_kW       unit duties (Q_cond negative = heat out)
Summary (cascade_balance / single_balance)
    Q_box_kW, W_total_kW, Q_air_kW, (Q_cascade_kW), energy_residual_kW

Use
    from heat_balance import cascade_balance, single_balance
    row.update(cascade_balance(cycle))
    cascade_balance_keys() gives the column names (for runners that
    declare their CSV columns up front, e.g. run_cascade_sweep_r515b.py).
"""

from pyomo.environ import value

K = 273.15
UNITS = (("evap", "evaporator"), ("comp", "compressor"), ("cond", "condenser"), ("valve", "expansion_valve"))


def _state(cv, side):
    st = cv.properties_in[0] if side == "in" else cv.properties_out[0]
    return dict(T_C=value(st.temperature) - K, P_kPa=value(st.pressure) / 1e3,
                h_kJkg=value(st.enth_mass) / 1e3)


def loop_report(fs, prefix=""):
    """Mass flow, states and duties for one loop (prefix "hot_", "cold_" or "")."""
    out = {}
    evap = getattr(fs, f"{prefix}evaporator")
    out[f"{prefix}m_kg_s"] = value(evap.control_volume.properties_in[0].flow_mass)
    for short, name in UNITS:
        unit = getattr(fs, f"{prefix}{name}")
        for side in ("in", "out"):
            for k, v in _state(unit.control_volume, side).items():
                out[f"{prefix}{short}_{side}_{k}"] = v
    out[f"{prefix}Q_evap_kW"] = value(evap.heat_duty[0]) / 1e3
    out[f"{prefix}W_comp_kW"] = value(getattr(fs, f"{prefix}compressor").work_mechanical[0]) / 1e3
    out[f"{prefix}Q_cond_kW"] = value(getattr(fs, f"{prefix}condenser").heat_duty[0]) / 1e3
    return out


def cascade_balance(cycle):
    fs = cycle.model.fs
    out = {}
    out.update(loop_report(fs, "cold_"))
    out.update(loop_report(fs, "hot_"))
    q_box, q_cas = out["cold_Q_evap_kW"], out["hot_Q_evap_kW"]
    w = out["cold_W_comp_kW"] + out["hot_W_comp_kW"]
    q_air = -out["hot_Q_cond_kW"]
    out.update(Q_box_kW=q_box, Q_cascade_kW=q_cas, W_total_kW=w, Q_air_kW=q_air,
               energy_residual_kW=q_box + w - q_air)
    return out


def single_balance(cycle):
    fs = cycle.model.fs
    out = loop_report(fs, "")
    q_box, w, q_air = out["Q_evap_kW"], out["W_comp_kW"], -out["Q_cond_kW"]
    out.update(Q_box_kW=q_box, W_total_kW=w, Q_air_kW=q_air, energy_residual_kW=q_box + w - q_air)
    return out


def loop_keys(prefix=""):
    """Column names loop_report() produces, in order."""
    keys = [f"{prefix}m_kg_s"]
    for short, _ in UNITS:
        for side in ("in", "out"):
            keys += [f"{prefix}{short}_{side}_{k}" for k in ("T_C", "P_kPa", "h_kJkg")]
    return keys + [f"{prefix}Q_evap_kW", f"{prefix}W_comp_kW", f"{prefix}Q_cond_kW"]


def cascade_balance_keys():
    """Column names cascade_balance() produces, in order."""
    return (loop_keys("cold_") + loop_keys("hot_")
            + ["Q_box_kW", "Q_cascade_kW", "W_total_kW", "Q_air_kW", "energy_residual_kW"])
