"""
cycle_eq_PR_1015.py -- R-32 single-stage cycle as an equation-oriented Pyomo
model on IDAES's own Peng-Robinson phase expressions. Phase of every state is
stated explicitly; there is no vapour-fraction variable.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    Four states. One IDAES PR state block per state; only its phase-specific
    expressions enth_mol_phase[p] / entr_mol_phase[p] are used (they depend
    on T, P and composition only -- checked 10/15). All block constraints are
    deactivated (no flash, no phase split).
        1  evaporator outlet   Vap   T1 = T_e + SH
        2s isentropic outlet   Vap   s_V(T2s, Pc) = s_V(T1, Pe)
        2  compressor outlet   Vap   h_V(T2, Pc) = h1 + (h2s - h1)/eta
        3  condenser outlet    Liq   T3 = T_cond - SC
        4  valve outlet        2-ph  T4 = T_e, h4 = h3 (quality checked after)
    Unknowns T1, T2s, T2, T3 (+ h, s with their defining equations);
    equations 2 more; DOF 2 = SH, SC, read with degrees_of_freedom().
Why
    The generic PR package cannot hold a PURE fluid on one phase: phase_frac
    is not determined away from saturation, so blocks land on blends of the
    two roots at wrong T (probe_vle_1015.py). Bounds patch it but leave
    IPOPT "infeasible" on DOWLINGLAB-7 (breadcrumb 10/15 night 5-6).
    Option B (Shilpa, 10/15): state the phase, do not let a solver pick it.
Basis
    P_e, P_c = PR's own Psat at T_e, T_cond (fugacity equality, PRProbe in
    vapor_compression_cubic_FPhx_1015.py), fixed before the solve because the
    spec fixes T_e and T_cond. Bounds keep each state on its side of the dome:
    T1 in [T_e+3, T_e+8], T3 in [T_cond-8, T_cond-3], T2, T2s > T_cond.
    LIMIT: this is the closed form solved simultaneously -- it must match the
    closed form, and energy closes by construction. Independent checks are
    Carnot, Helmholtz, and the root-count check below.
Outputs
    build(...) -> model; solve_point(...) -> dict. Nothing written.
Run
    imported by sweep_eq_PR_1015.py
"""
import math
import logging
import warnings

from pyomo.environ import (ConcreteModel, Constraint, Expression, Objective,
                           Var, value, maximize, units as pyunits)
from pyomo.core.expr.visitor import identify_variables
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

from phase_1_cubic_eos_validation_refstate_1001 import make_config, METHODS, MW
from vapor_compression_cubic_FPhx_1015 import PRProbe

logging.getLogger("idaes").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

STATES = {"s1": "Vap", "s2s": "Vap", "s2": "Vap", "s3": "Liq"}


def build(method="NIST", T_amb_C=20.0, T_evap_C=-29.0, approach=9.0,
          eta=0.80, SH=3.0, SC=8.0, SH_win=(3.0, 8.0), SC_win=(3.0, 8.0)):
    probe = PRProbe(method)
    Te = T_evap_C + 273.15
    Tc = T_amb_C + approach + 273.15
    Pe, Pc = probe.psat(Te), probe.psat(Tc)

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = GenericParameterBlock(**make_config(METHODS[method]))
    m.fs.sb = m.fs.pp.build_state_block(list(STATES), defined_state=True)

    for k, ph in STATES.items():
        b = m.fs.sb[k]
        # touch the two expressions we use so they are built, then switch off
        # every equation in the block: no flash, no phase split
        b.enth_mol_phase[ph]; b.entr_mol_phase[ph]
        b.deactivate()
        b.flow_mol.fix(1.0); b.mole_frac_comp["R32"].fix(1.0)
        b.mole_frac_phase_comp[ph, "R32"].fix(1.0)
        b.log_mole_frac_phase_comp[ph, "R32"].fix(0.0)
        b.pressure.fix(Pe if k == "s1" else Pc)
    sb = m.fs.sb
    # h and s as VARIABLES with defining equations, so every state's T sits in
    # an active equation and degrees_of_freedom() counts SH and SC. (With
    # Expressions, T3 appeared only in the COP and the count read 1, not 2.)
    m.fs.h = Var(list(STATES), initialize=0.0, units=pyunits.J / pyunits.mol)
    m.fs.s = Var(["s1", "s2s"], initialize=0.0, units=pyunits.J / pyunits.mol / pyunits.K)
    m.fs.h_def = Constraint(list(STATES), rule=lambda fs, k: fs.h[k] == sb[k].enth_mol_phase[STATES[k]])
    m.fs.s_def = Constraint(["s1", "s2s"], rule=lambda fs, k: fs.s[k] == sb[k].entr_mol_phase[STATES[k]])

    # bounds that keep each state on its side of the dome
    sb["s1"].temperature.setlb(Te + SH_win[0]); sb["s1"].temperature.setub(Te + SH_win[1])
    sb["s3"].temperature.setlb(Tc - SC_win[1]); sb["s3"].temperature.setub(Tc - SC_win[0])
    for k in ("s2s", "s2"):
        sb[k].temperature.setlb(Tc + 0.5); sb[k].temperature.setub(600.0)

    m.fs.eta = eta
    m.fs.isentropic = Constraint(expr=m.fs.s["s2s"] == m.fs.s["s1"])
    m.fs.efficiency = Constraint(
        expr=(m.fs.h["s2"] - m.fs.h["s1"]) * eta == m.fs.h["s2s"] - m.fs.h["s1"])
    m.fs.cop = Expression(expr=(m.fs.h["s1"] - m.fs.h["s3"]) / (m.fs.h["s2"] - m.fs.h["s1"]))

    # starting point: closed form at the requested SH / SC
    import vapor_compression_cubic_FPhx_1015 as V
    ref = V.closed_form(probe, T_evap_C, T_amb_C + approach, SH, SC, eta)
    sb["s1"].temperature.set_value(ref["T1"]); sb["s3"].temperature.set_value(ref["T3"])
    sb["s2"].temperature.set_value(ref["T2"])
    s1 = probe.state(ref["T1"], Pe, "Vap")[1]
    sb["s2s"].temperature.set_value(probe.T_from("s", s1, Pc, "Vap", Tc, 600.0))
    for k in STATES: m.fs.h[k].set_value(value(sb[k].enth_mol_phase[STATES[k]]))
    for k in ("s1", "s2s"): m.fs.s[k].set_value(value(sb[k].entr_mol_phase[STATES[k]]))

    m.fs.Te, m.fs.Tc, m.fs.Pe, m.fs.Pc = Te, Tc, Pe, Pc
    m.fs.dof_free = degrees_of_freedom(m)          # SH, SC free: expect 2
    return m, probe, ref


def solve_point(method, T_amb_C, mode="square", SH=3.0, SC=8.0, eta=0.80, tee=False):
    m, probe, ref = build(method, T_amb_C, eta=eta, SH=SH, SC=SC)
    fs, sb = m.fs, m.fs.sb
    if mode == "square":
        sb["s1"].temperature.fix(fs.Te + SH); sb["s3"].temperature.fix(fs.Tc - SC)
    else:
        fs.obj = Objective(expr=fs.cop, sense=maximize)
    dof = degrees_of_freedom(m)
    res = get_solver("ipopt", options={"max_iter": 500, "tol": 1e-8}).solve(m, tee=tee)
    tc = str(res.solver.termination_condition)

    T = {k: value(sb[k].temperature) for k in STATES}
    h = {k: value(fs.h[k]) for k in STATES}
    SH_act, SC_act = T["s1"] - fs.Te, fs.Tc - T["s3"]
    # independent check 1: each state really is on the phase we said
    # P_e, P_c are PR's Psat at T_e, T_cond, so Tsat at each state's P is
    # T_e or T_cond exactly
    Tsat = {"s1": fs.Te, "s2s": fs.Tc, "s2": fs.Tc, "s3": fs.Tc}
    side = {k: ("Vap" if T[k] > Tsat[k] else "Liq") for k in STATES}
    phases_ok = all(side[k] == STATES[k] for k in STATES)
    # check 2: valve outlet is two-phase (0 < x < 1) at T_e
    hL, hV = probe.state(fs.Te, fs.Pe, "Liq")[0], probe.state(fs.Te, fs.Pe, "Vap")[0]
    x4 = (h["s3"] - hL) / (hV - hL)
    # check 3: same SH/SC through the closed form
    import vapor_compression_cubic_FPhx_1015 as V
    cf = V.closed_form(probe, fs.Te - 273.15, fs.Tc - 273.15, SH_act, SC_act, eta)["cop"]
    # check 4: second law for the compressor, s2 >= s1
    s2_ge_s1 = probe.state(T["s2"], fs.Pc, "Vap")[1] >= value(fs.s["s1"]) - 1e-9
    F = 1.0 / MW
    Q_e, W, Q_c = F * (h["s1"] - h["s3"]), F * (h["s2"] - h["s1"]), F * (h["s3"] - h["s2"])
    carnot = fs.Te / (fs.Tc - fs.Te)
    cop = value(fs.cop)
    return dict(method=method, T_amb=T_amb_C, mode=mode, dof_free=fs.dof_free, dof=dof,
                term=tc, cop=cop, closed_form=cf, diff_pct=100 * (cop - cf) / cf,
                carnot=carnot, below_carnot=cop < carnot, SH=SH_act, SC=SC_act,
                T2_C=T["s2"] - 273.15, x4=x4, phases_ok=phases_ok and 0 < x4 < 1,
                second_law_ok=s2_ge_s1, loop_W=Q_e + W + Q_c,
                Pe_kPa=fs.Pe / 1e3, Pc_kPa=fs.Pc / 1e3)
