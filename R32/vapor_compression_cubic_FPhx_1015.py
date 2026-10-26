"""
vapor_compression_cubic_FPhx_1015.py -- R-32 single-stage cycle on the generic
Peng-Robinson package with the FPhx state definition (flow, x, P, h).

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    Builds evaporator -> compressor -> condenser -> valve on a PR
    GenericParameterBlock whose state variables are (flow_mol, mole_frac_comp,
    pressure, enth_mol). Seeds every state block deterministically, solves,
    and returns the model for cycle_guards.check_cycle().

Why
    FTPx ports carry T and P, which do not fix the state of a pure fluid at
    saturation; the loop did not conserve energy (10/08). FPhx ports carry
    enth_mol, so expand_arcs writes the enthalpy equality itself and (P, h)
    fixes quality everywhere. See BREADCRUMB_09-28.md, 10/15 evening + late.

    This is a NEW lean class, not a patch of _1001. The _1001 initialization
    is ~600 lines of FTPx-specific workarounds (phase_frac fixes, temperature
    fences, sum_mole_frac deactivations) that do not apply under FPhx.

Basis
    Seeds: pr_phase_state() reads enth_mol_phase straight off a PR block for
    the phase chosen by the sign of ln(phi_L) - ln(phi_V). No solve, so no
    wrong-root (10/15 late: the SmoothVLE flash returned a 61 % vapour
    "subcooled liquid" and reported optimal).
    Saturation pressures are PR's OWN (fugacity equality), not CoolProp's.

Outputs
    build_cycle(...) -> (model, info dict). Nothing written to disk.

Status 10/15
    DOF verified: 0 square, 2 with SH/SC free. The cycle does NOT solve:
    the PR block can land on a blend of liquid and vapour roots at a wrong T.
    PRProbe.psat/tsat bracket on nonzero ln(phi_L)-ln(phi_V) first (a plain
    bisection drifted to Pc for SPGP). See BREADCRUMB_09-28.md, 10/15 night.

Run
    imported by run_FPhx_cycle_1015.py and closedform_4methods_1015.py
"""
import math
import logging

from pyomo.environ import (ConcreteModel, Constraint, Objective, Var, value,
                           TransformationFactory, maximize, units as pyunits)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Heater, Compressor, PressureChanger
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.state_definitions import FPhx
from idaes.models.properties.modular_properties.eos.ceos import Cubic

from phase_1_cubic_eos_validation_refstate_1001 import make_config, METHODS, MW

logging.getLogger("idaes").setLevel(logging.ERROR)


def make_config_FPhx(p):
    cfg = make_config(p)
    cfg["state_definition"] = FPhx
    cfg["state_bounds"] = {
        "flow_mol": (0.0, 19.2, 1000.0, pyunits.mol / pyunits.s),
        "enth_mol": (-6.0e4, -5.0e3, 6.0e4, pyunits.J / pyunits.mol),
        "temperature": (180.0, 300.0, 500.0, pyunits.K),
        "pressure": (1.0e4, 1.0e5, 1.0e7, pyunits.Pa),
    }
    return cfg


# ---------------------------------------------------------------- PR probe --
class PRProbe:
    """Evaluates PR expressions on a scratch FTPx block. No solver."""

    def __init__(self, method):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = GenericParameterBlock(**make_config(METHODS[method]))
        m.fs.sb = m.fs.pp.build_state_block([0], defined_state=True)
        self.b = m.fs.sb[0]
        self.m = m
        self.Tc = METHODS[method]["Tc"]; self.Pc = METHODS[method]["Pc"]

    def _set(self, T, P):
        b = self.b
        b.temperature.set_value(T); b.pressure.set_value(P)
        b.flow_mol.set_value(1.0); b.mole_frac_comp["R32"].set_value(1.0)
        for ph in ("Liq", "Vap"):
            b.mole_frac_phase_comp[ph, "R32"].set_value(1.0)

    def dlnphi(self, T, P):
        self._set(T, P)
        return (math.log(value(Cubic.fug_coeff_phase_comp(self.b, "Liq", "R32")))
                - math.log(value(Cubic.fug_coeff_phase_comp(self.b, "Vap", "R32"))))

    def tsat(self, P, lo=200.0, hi=None):
        # same bracketing logic as psat: dlnphi < 0 below Tsat (liquid)
        import numpy as np
        hi = hi or self.Tc - 0.01
        grid = np.linspace(lo, hi, 120)
        vals = [self.dlnphi(T, P) for T in grid]
        a = b = None
        for i in range(len(grid) - 1):
            if vals[i] < 0 and vals[i + 1] > 0: a, b = grid[i], grid[i + 1]; break
        if a is None: raise ValueError(f"no PR Tsat bracket at P={P:.0f} Pa")
        for _ in range(60):
            mid = 0.5 * (a + b); f = self.dlnphi(mid, P)
            if f < 0: a = mid
            elif f > 0: b = mid
            else: return mid  # exact equality = converged
        return 0.5 * (a + b)

    def psat(self, T, lo=1.0e4, hi=None):
        # dlnphi < 0 -> liquid stable -> P above Psat. Outside the 3-root
        # window both phases take the same root and dlnphi == 0 exactly, so
        # bracket on a scan of NONZERO values first (10/15: a plain bisection
        # drifted to Pc for SPGP, whose Tc is only 52 C).
        import numpy as np
        hi = hi or self.Pc * 0.999
        grid = np.geomspace(lo, hi, 120)
        vals = [self.dlnphi(T, P) for P in grid]
        a = b = None
        for i in range(len(grid) - 1):
            if vals[i] > 0 and vals[i + 1] < 0: a, b = grid[i], grid[i + 1]; break
        if a is None: raise ValueError(f"no PR Psat bracket at T={T:.2f} K")
        for _ in range(60):
            mid = math.sqrt(a * b); f = self.dlnphi(T, mid)
            if f < 0: b = mid
            elif f > 0: a = mid
            else: return mid  # exact equality = converged
        return math.sqrt(a * b)

    def state(self, T, P, phase=None):
        """(h, s, phase) on the stable root unless phase is forced."""
        if phase is None:
            phase = "Vap" if T > self.tsat(P) else "Liq"
        self._set(T, P)
        return (value(self.b.enth_mol_phase[phase]),
                value(self.b.entr_mol_phase[phase]), phase)

    def T_from(self, prop, target, P, phase, lo, hi):
        idx = 0 if prop == "h" else 1
        f = lambda T: self.state(T, P, phase)[idx] - target
        flo = f(lo)
        for _ in range(80):
            mid = 0.5 * (lo + hi); fm = f(mid)
            if (fm > 0) == (flo > 0): lo, flo = mid, fm
            else: hi = mid
        return 0.5 * (lo + hi)


def closed_form(probe, T_evap_C, T_cond_C, SH, SC, eta):
    """Four state points from the same PR expressions. Reference for the solve."""
    Pe = probe.psat(T_evap_C + 273.15); Pc = probe.psat(T_cond_C + 273.15)
    T1 = T_evap_C + 273.15 + SH
    h1, s1, _ = probe.state(T1, Pe, "Vap")
    T2s = probe.T_from("s", s1, Pc, "Vap", T_cond_C + 273.15, 600.0)
    h2s = probe.state(T2s, Pc, "Vap")[0]
    h2 = h1 + (h2s - h1) / eta
    T2 = probe.T_from("h", h2, Pc, "Vap", T_cond_C + 273.15, 700.0)
    T3 = T_cond_C + 273.15 - SC
    h3 = probe.state(T3, Pc, "Liq")[0]
    return dict(Pe=Pe, Pc=Pc, T1=T1, h1=h1, T2=T2, h2=h2, T3=T3, h3=h3,
                h4=h3, cop=(h1 - h3) / (h2 - h1))


# ---------------------------------------------------------------- flowsheet --
def build_cycle(method="NIST", T_amb_C=20.0, T_evap_C=-29.0, approach=9.0,
                SH=3.0, SC=8.0, eta=0.80, flow_kg_s=1.0):
    probe = PRProbe(method)
    T_cond_C = T_amb_C + approach
    ref = closed_form(probe, T_evap_C, T_cond_C, SH, SC, eta)

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**make_config_FPhx(METHODS[method]))
    fs = m.fs
    fs.evaporator = Heater(property_package=fs.properties)
    fs.compressor = Compressor(property_package=fs.properties,
                               thermodynamic_assumption=ThermodynamicAssumption.isentropic)
    fs.condenser = Heater(property_package=fs.properties)
    fs.expansion_valve = PressureChanger(property_package=fs.properties, compressor=False,
                                         thermodynamic_assumption=ThermodynamicAssumption.adiabatic)
    fs.a1 = Arc(source=fs.evaporator.outlet, destination=fs.compressor.inlet)
    fs.a2 = Arc(source=fs.compressor.outlet, destination=fs.condenser.inlet)
    fs.a3 = Arc(source=fs.condenser.outlet, destination=fs.expansion_valve.inlet)
    fs.a4 = Arc(source=fs.expansion_valve.outlet, destination=fs.evaporator.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    sb = {  # the eight state blocks
        "evap_in": fs.evaporator.control_volume.properties_in[0],
        "evap_out": fs.evaporator.control_volume.properties_out[0],
        "comp_in": fs.compressor.control_volume.properties_in[0],
        "comp_out": fs.compressor.control_volume.properties_out[0],
        "cond_in": fs.condenser.control_volume.properties_in[0],
        "cond_out": fs.condenser.control_volume.properties_out[0],
        "valve_in": fs.expansion_valve.control_volume.properties_in[0],
        "valve_out": fs.expansion_valve.control_volume.properties_out[0],
    }
    info = dict(method=method, T_amb_C=T_amb_C, T_cond_C=T_cond_C, ref=ref,
                dof_unspecified=degrees_of_freedom(m))

    # ---- specifications (closed form uses the same ones) ----
    Ts_e = probe.tsat(ref["Pe"]); Ts_c = probe.tsat(ref["Pc"])
    fs.Tsat_evap = Ts_e; fs.Tsat_cond = Ts_c   # plain floats, PR's own Tsat
    sb["evap_in"].flow_mol.fix(flow_kg_s / MW)
    sb["evap_in"].mole_frac_comp["R32"].fix(1.0)
    sb["evap_in"].pressure.fix(ref["Pe"])
    fs.compressor.outlet.pressure[0].fix(ref["Pc"])
    fs.compressor.efficiency_isentropic.fix(eta)
    # 2 design variables: outlet temperatures (superheat, subcool)
    fs.T_evap_out = Var(initialize=Ts_e + SH, units=pyunits.K)
    fs.T_cond_out = Var(initialize=Ts_c - SC, units=pyunits.K)
    fs.T_evap_out_def = Constraint(expr=sb["evap_out"].temperature == fs.T_evap_out)
    fs.T_cond_out_def = Constraint(expr=sb["cond_out"].temperature == fs.T_cond_out)
    fs.T_evap_out.fix(Ts_e + SH); fs.T_cond_out.fix(Ts_c - SC)
    info["dof_before_loop_fixes"] = degrees_of_freedom(m)
    # Closed loop: flow and composition are each conserved by every unit, so
    # one arc's flow and composition equalities are implied by the other
    # three. Remove exactly those two; cycle_guards check 4 verifies flow.
    fs.a4_expanded.flow_mol_equality.deactivate()
    fs.a4_expanded.mole_frac_comp_equality.deactivate()
    info["dof_square"] = degrees_of_freedom(m)
    fs.T_evap_out.unfix(); fs.T_cond_out.unfix()
    info["dof_design"] = degrees_of_freedom(m)
    fs.T_evap_out.fix(); fs.T_cond_out.fix()
    return m, sb, info, probe


def seed(m, sb, info, probe):
    """Deterministic seeds from the closed form; phase chosen by fugacity."""
    r = info["ref"]
    pts = {"evap_out": (r["T1"], r["Pe"]), "comp_in": (r["T1"], r["Pe"]),
           "comp_out": (r["T2"], r["Pc"]), "cond_in": (r["T2"], r["Pc"]),
           "cond_out": (r["T3"], r["Pc"]), "valve_in": (r["T3"], r["Pc"])}
    F = 1.0 / MW
    for k, (T, P) in pts.items():
        h, _, ph = probe.state(T, P)
        _seed_block(sb[k], F, P, T, h, 1.0 if ph == "Vap" else 0.0)
    # two-phase valve outlet / evap inlet: lever rule at PR saturation
    Ts = probe.tsat(r["Pe"])
    hL = probe.state(Ts, r["Pe"], "Liq")[0]; hV = probe.state(Ts, r["Pe"], "Vap")[0]
    q = (r["h4"] - hL) / (hV - hL)
    for k in ("valve_out", "evap_in"):
        _seed_block(sb[k], F, r["Pe"], Ts, r["h4"], q)
    info["seed_quality_valve_out"] = q


def _seed_block(b, F, P, T, h, vf):
    b.flow_mol.set_value(F); b.mole_frac_comp["R32"].set_value(1.0)
    b.pressure.set_value(P); b.enth_mol.set_value(h); b.temperature.set_value(T)
    b.phase_frac["Vap"].set_value(vf); b.phase_frac["Liq"].set_value(1 - vf)
    for ph in ("Liq", "Vap"):
        b.flow_mol_phase[ph].set_value(F * (vf if ph == "Vap" else 1 - vf))
        b.mole_frac_phase_comp[ph, "R32"].set_value(1.0)
    for name in ("_teq", "_t1"):
        c = getattr(b, name, None)
        if c is not None:
            for i in c: c[i].set_value(T)
    for name in ("temperature_bubble", "temperature_dew"):
        c = getattr(b, name, None)
        if c is not None:
            for i in c: c[i].set_value(T)


def phase_check(sb, probe, tol_K=0.05):
    """State-block root check: phase_frac must match the side of PR's own
    saturation line. Returns list of (block, T, Tsat, vap_frac, verdict)."""
    out = []
    for k, b in sb.items():
        T, P, vf = value(b.temperature), value(b.pressure), value(b.phase_frac["Vap"])
        Ts = probe.tsat(P)
        if T > Ts + tol_K: ok = vf > 0.99
        elif T < Ts - tol_K: ok = vf < 0.01
        else: ok = True  # on the line: two-phase allowed
        out.append((k, T - 273.15, Ts - 273.15, vf, "ok" if ok else "WRONG ROOT"))
    return out


def solve(m, tee=False):
    solver = get_solver("ipopt", options={"max_iter": 300, "tol": 1e-8})
    res = solver.solve(m, tee=tee)
    return str(res.solver.termination_condition)


def seed_blockwise(m, sb, info, probe, tee=False):
    """Seed v2. For each state block (incl. the isentropic ghost) fix the state
    vars (F, x, P, h) at the closed-form values, initialize the block's own
    flash in isolation, then release. Under FPhx the flash at fixed (P, h) is
    well posed -- h decides quality -- unlike (T, P) on the dome."""
    from pyomo.environ import Var
    r = info["ref"]; F = 1.0 / MW
    seed(m, sb, info, probe)  # T, phase_frac guesses
    blocks = dict(sb)
    iso = m.fs.compressor.properties_isentropic[0]
    T2s = probe.T_from("s", probe.state(r["T1"], r["Pe"], "Vap")[1], r["Pc"], "Vap",
                       r["Pc"] and info["T_cond_C"] + 273.15, 600.0)
    h2s = probe.state(T2s, r["Pc"], "Vap")[0]
    _seed_block(iso, F, r["Pc"], T2s, h2s, 1.0)
    blocks["comp_isentropic"] = iso
    report = {}
    for k, b in blocks.items():
        b.parent_component().initialize(state_args={"flow_mol": F, "mole_frac_comp": {"R32": 1.0},
                                 "pressure": value(b.pressure),
                                 "enth_mol": value(b.enth_mol)},
                     hold_state=False)
        report[k] = (value(b.temperature) - 273.15, value(b.phase_frac["Vap"]))
    info["seed_report"] = report
    return report


# phase of each single-phase block, in flow order. Two-phase blocks
# (valve_out, evap_in) are left free: there T is pinned to Tsat anyway.
SINGLE_PHASE = {"evap_out": "Vap", "comp_in": "Vap", "comp_out": "Vap",
                "cond_in": "Vap", "cond_out": "Liq", "valve_in": "Liq"}


def bound_phases(m, sb, eps=1e-5):
    """Bound phase_frac to the physical side on single-phase blocks.

    Why: for a PURE component nothing in SmoothVLE forces phase_frac to 0 or
    1 away from saturation, so a block can sit at a wrong T as a blend of the
    liquid and vapour roots (probe_vle_1015.py). A bound does not change DOF.
    probe_bound_1015.py: the T error scales with eps, 10x per decade --
    1e-3 -> 0.32 K, 1e-4 -> 0.034 K, 1e-5 -> 0.004 K; 1e-6 is sometimes
    infeasible. 1e-5 is the default."""
    blocks = dict(sb)
    blocks["comp_isentropic"] = m.fs.compressor.properties_isentropic[0]
    phases = dict(SINGLE_PHASE, comp_isentropic="Vap")
    for k, ph in phases.items():
        if ph == "Vap": blocks[k].phase_frac["Vap"].setlb(1 - eps)
        else: blocks[k].phase_frac["Vap"].setub(eps)


def seed_bounded(m, sb, info, probe, eps=1e-5):
    """Seed v3. Closed-form guesses, phase bounds on, then solve each state
    block ALONE at fixed (F, x, P, h) -- the same set-up probe_bound_1015.py
    verified to within 0.004 K. IDAES's own block initialize() is not used:
    it ignores the bounds and fails on the two-phase blocks."""
    from pyomo.environ import Block
    r = info["ref"]; F = 1.0 / MW
    seed(m, sb, info, probe)
    iso = m.fs.compressor.properties_isentropic[0]
    s1 = probe.state(r["T1"], r["Pe"], "Vap")[1]
    T2s = probe.T_from("s", s1, r["Pc"], "Vap", info["T_cond_C"] + 273.15, 600.0)
    _seed_block(iso, F, r["Pc"], T2s, probe.state(T2s, r["Pc"], "Vap")[0], 1.0)
    bound_phases(m, sb, eps)
    blocks = dict(sb); blocks["comp_isentropic"] = iso
    solver = get_solver("ipopt", options={"max_iter": 500})
    report = {}
    for k, b in blocks.items():
        sv = [b.flow_mol, b.mole_frac_comp["R32"], b.pressure, b.enth_mol]
        was = [v.fixed for v in sv]
        for v in sv: v.fix()
        tmp = Block(concrete=True)   # solve just this block
        res = solver.solve(b)
        for v, w in zip(sv, was):
            if not w: v.unfix()
        report[k] = (str(res.solver.termination_condition),
                     value(b.temperature) - 273.15, value(b.phase_frac["Vap"]))
    info["seed_report"] = report
    return report


def bound_two_phase_T(m, sb, probe, dT=0.01):
    """Two-phase blocks (valve_out, evap_in) must sit AT PR's boiling
    temperature. Same pure-component degeneracy as the single-phase blocks:
    without this the solve put them at -26.3 C, 24 % vapour, P = Psat(-29 C)
    (10/15, NIST 20 C) -- right enthalpy, impossible state. The pressure is
    fixed, so Tsat is known; a bound does not change DOF."""
    for k in ("valve_out", "evap_in"):
        b = sb[k]
        Ts = probe.tsat(value(b.pressure))
        b.temperature.setlb(Ts - dT); b.temperature.setub(Ts + dT)
