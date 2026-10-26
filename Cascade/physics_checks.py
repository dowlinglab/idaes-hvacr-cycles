"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/11/2026

What this does
    Physics checks on one solved cycle. Returns a dict of numbers plus a
    list of failed checks. Empty list = physically sane.

Why
    converged=True only means IPOPT is happy. It does not mean the cycle
    is physical. Every COP we report goes through these first.

Checks (cascade)
    energy      |Q_box + W_hot + W_cold - Q_air| / Q_air        < 1e-4
    carnot      COP < T_box_sat / (T_cond_sat - T_box_sat)
    comp_dh     compressor h_out - h_in > 0 (both loops)
    comp_vap    compressor outlet is vapor (x = 1, T > T_sat)
    valve_h     valve h_out = h_in (isenthalpic)
    valve_2ph   valve outlet two-phase (0 < x < 1)
    sh_sc       superheat / subcool inside their windows
    ratio       pressure ratio <= cap
    cross_hot   cascade HX hot end: CO2 in > hot-loop vapor out
    cross_cold  cascade HX cold end: CO2 out - T_CAS,E >= cold_end_dT
    cross_air   air condenser: refrigerant out >= ambient the model sees
    cross_box   CO2 evaporator out < -20 C box air
    co2_subcrit CO2 condensing below its critical point
    p_bound     no port pressure sits on its lower/upper bound (a bound
                that binds means the "optimum" is set by a safety bound,
                not physics; found 10/11: the engine's derived hot-loop
                pressure floor pinned T_CAS,E at about -26 C)
Info only (no pass/fail)
    t_dis_hot / t_dis_cold   discharge temperatures, C
    vacuum                   hot-loop evaporating pressure < 1 atm

Checks (single stage): same idea, one loop: energy, carnot, comp_dh,
comp_vap, valve_h, valve_2ph, sh_sc, ratio, cross_air, cross_box.

Use
    from physics_checks import cascade_checks, single_checks
    info, failed = cascade_checks(cycle, ambient_model_C)
"""

from pyomo.environ import value

K = 273.15
BOX_AIR_C = -20.0                   # cold-storage box air
CO2_P_CRIT_PA = 7.3773e6
SH_SC_TOL_K = 0.05                  # slack on the SH/SC windows


def _st(unit, side):
    """State block at a unit's inlet ('in') or outlet ('out')."""
    cv = unit.control_volume
    return cv.properties_in[0] if side == "in" else cv.properties_out[0]


def _x(state):
    """Vapor fraction, or None if this package doesn't have it."""
    try:
        return value(state.vapor_frac)
    except Exception:
        return None


def _pressure_on_bound(fs, units, rel_tol=1e-6):
    """Names of port pressures sitting on a bound (relative tol)."""
    hits = []
    for u in units:
        unit = getattr(fs, u)
        for side in ("inlet", "outlet"):
            v = getattr(unit, side).pressure[0]
            p = value(v)
            for bnd in (v.lb, v.ub):
                if bnd is not None and abs(p - bnd) <= rel_tol * max(abs(bnd), 1.0):
                    hits.append(f"{u}.{side}")
    return hits


def _loop(fs, prefix, info, failed, sh_min, sh_max, sc_min, sc_max, ratio_cap, tag):
    """Checks shared by any loop: compressor, valve, SH/SC, ratio."""
    comp = getattr(fs, f"{prefix}compressor")
    valve = getattr(fs, f"{prefix}expansion_valve")
    evap = getattr(fs, f"{prefix}evaporator")
    cond = getattr(fs, f"{prefix}condenser")

    ci, co = _st(comp, "in"), _st(comp, "out")
    dh = (value(co.enth_mass) - value(ci.enth_mass)) / 1e3
    info[f"{tag}comp_dh_kJkg"] = dh
    if dh <= 0:
        failed.append(f"{tag}comp_dh")

    x_out = _x(co)
    t_out, t_sat_out = value(co.temperature), value(co.temperature_sat)
    info[f"{tag}comp_x_out"] = x_out
    info[f"{tag}t_dis_C"] = t_out - K
    if (x_out is not None and x_out < 0.999) or t_out < t_sat_out - 0.01:
        failed.append(f"{tag}comp_vap")

    vi, vo = _st(valve, "in"), _st(valve, "out")
    dh_v = abs(value(vo.enth_mass) - value(vi.enth_mass)) / 1e3
    info[f"{tag}valve_dh_kJkg"] = dh_v
    if dh_v > 1e-3:
        failed.append(f"{tag}valve_h")
    xv = _x(vo)
    info[f"{tag}valve_x_out"] = xv
    if xv is not None and not (0.0 < xv < 1.0):
        failed.append(f"{tag}valve_2ph")

    eo = _st(evap, "out")
    sh = value(eo.temperature) - value(eo.temperature_sat)
    co_ = _st(cond, "out")
    sc = value(co_.temperature_sat) - value(co_.temperature)
    info[f"{tag}sh_K"], info[f"{tag}sc_K"] = sh, sc
    if not (sh_min - SH_SC_TOL_K <= sh <= sh_max + SH_SC_TOL_K):
        failed.append(f"{tag}sh")
    if not (sc_min - SH_SC_TOL_K <= sc <= sc_max + SH_SC_TOL_K):
        failed.append(f"{tag}sc")

    ratio = value(co.pressure) / value(ci.pressure)
    info[f"{tag}ratio"] = ratio
    if ratio > ratio_cap * 1.001:
        failed.append(f"{tag}ratio")


def cascade_checks(cycle, ambient_model_C, sh=(3, 8), sc_hot=(3, 8), sc_cold=(0, 8),
                   ratio_cap=20.0, cold_end_dT=1.0):
    """All cascade checks. ambient_model_C = ambient the model sees
    (after head pressure control). Returns (info dict, failed list)."""
    fs = cycle.model.fs
    info, failed = {}, []

    q_box = value(fs.cold_evaporator.heat_duty[0])
    q_air = -value(fs.hot_condenser.heat_duty[0])
    w_hot = value(fs.hot_compressor.work_mechanical[0])
    w_cold = value(fs.cold_compressor.work_mechanical[0])
    info["energy_rel"] = abs(q_box + w_hot + w_cold - q_air) / abs(q_air)
    if info["energy_rel"] > 1e-4:
        failed.append("energy")

    cop = q_box / (w_hot + w_cold)
    t_box = value(_st(fs.cold_evaporator, "out").temperature_sat)
    t_cond = value(_st(fs.hot_condenser, "out").temperature_sat)
    info["cop_check"] = cop
    info["cop_carnot"] = t_box / (t_cond - t_box)
    if cop >= info["cop_carnot"]:
        failed.append("carnot")

    _loop(fs, "hot_", info, failed, *sh, *sc_hot, ratio_cap, "hot_")
    _loop(fs, "cold_", info, failed, *sh, *sc_cold, ratio_cap, "cold_")

    # Cascade HX: CO2 condenses (hot side), hot loop evaporates (cold side)
    co2_in = value(_st(fs.cold_condenser, "in").temperature)
    co2_out = value(_st(fs.cold_condenser, "out").temperature)
    he_in = value(_st(fs.hot_evaporator, "in").temperature)
    he_out = value(_st(fs.hot_evaporator, "out").temperature)
    info["cross_hot_K"] = co2_in - he_out
    info["cross_cold_K"] = co2_out - he_in
    info["T_cas_e_C"] = value(_st(fs.hot_evaporator, "out").temperature_sat) - K
    if info["cross_hot_K"] <= 0:
        failed.append("cross_hot")
    if info["cross_cold_K"] < cold_end_dT - 0.01:
        failed.append("cross_cold")

    t_air_out = value(_st(fs.hot_condenser, "out").temperature)
    info["cross_air_K"] = t_air_out - (ambient_model_C + K)
    if info["cross_air_K"] < 0:
        failed.append("cross_air")

    t_box_out = value(_st(fs.cold_evaporator, "out").temperature)
    info["cross_box_K"] = (BOX_AIR_C + K) - t_box_out
    if info["cross_box_K"] <= 0:
        failed.append("cross_box")

    p_co2 = value(_st(fs.cold_condenser, "in").pressure)
    info["co2_P_MPa"] = p_co2 / 1e6
    if p_co2 >= CO2_P_CRIT_PA:
        failed.append("co2_subcrit")

    units = [f"{r}_{u}" for r in ("hot", "cold")
             for u in ("evaporator", "compressor", "condenser", "expansion_valve")]
    info["p_on_bound"] = " ".join(_pressure_on_bound(fs, units))
    if info["p_on_bound"]:
        failed.append("p_bound")

    p_low_hot = value(_st(fs.hot_evaporator, "out").pressure)
    info["hot_P_evap_kPa"] = p_low_hot / 1e3
    info["vacuum"] = p_low_hot < 101325.0
    return info, failed


def single_checks(cycle, ambient_model_C, sh=(3, 8), sc=(3, 8), ratio_cap=20.0):
    """All single-stage checks. Returns (info dict, failed list)."""
    fs = cycle.model.fs
    info, failed = {}, []

    q_box = value(fs.evaporator.heat_duty[0])
    q_air = -value(fs.condenser.heat_duty[0])
    w = value(fs.compressor.work_mechanical[0])
    info["energy_rel"] = abs(q_box + w - q_air) / abs(q_air)
    if info["energy_rel"] > 1e-4:
        failed.append("energy")

    cop = q_box / w
    t_box = value(_st(fs.evaporator, "out").temperature_sat)
    t_cond = value(_st(fs.condenser, "out").temperature_sat)
    info["cop_check"] = cop
    info["cop_carnot"] = t_box / (t_cond - t_box)
    if cop >= info["cop_carnot"]:
        failed.append("carnot")

    _loop(fs, "", info, failed, *sh, *sc, ratio_cap, "")

    t_air_out = value(_st(fs.condenser, "out").temperature)
    info["cross_air_K"] = t_air_out - (ambient_model_C + K)
    if info["cross_air_K"] < 0:
        failed.append("cross_air")
    t_box_out = value(_st(fs.evaporator, "out").temperature)
    info["cross_box_K"] = (BOX_AIR_C + K) - t_box_out
    if info["cross_box_K"] <= 0:
        failed.append("cross_box")

    info["p_on_bound"] = " ".join(_pressure_on_bound(
        fs, ["evaporator", "compressor", "condenser", "expansion_valve"]))
    if info["p_on_bound"]:
        failed.append("p_bound")

    p_low = value(_st(fs.evaporator, "out").pressure)
    info["P_evap_kPa"] = p_low / 1e3
    info["vacuum"] = p_low < 101325.0
    return info, failed
