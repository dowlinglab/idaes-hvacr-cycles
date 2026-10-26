"""
probe_bound_1015.py -- does a PR state block STAY at the true state, and does
a phase_frac bound force it there?

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    One PR state block (FPhx, P and h fixed) at four known states. Three
    starts: seeded at truth; seeded at truth + phase_frac bound; default
    start + bound. SmoothVLE and CubicComplementarityVLE.
Why
    probe_vle_1015.py showed blocks landing on impossible states. This says
    whether that is a bad START (fixable by seeding) or a NON-UNIQUE model
    (needs an extra condition such as the bound).
Basis
    Truth from PRProbe (fugacity sign, no solver).
Outputs
    table on screen.
Run
    python probe_bound_1015.py 2>&1 | tee probe_bound_1015.log
"""
import logging, warnings
logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
from pyomo.environ import ConcreteModel, value, units as pyunits
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.state_definitions import FPhx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.smooth_VLE_2 import CubicComplementarityVLE
import vapor_compression_cubic_FPhx_1015 as V
from phase_1_cubic_eos_validation_refstate_1001 import make_config, METHODS

probe = V.PRProbe("NIST")
cases = [("evap_out SH3", 247.15, 283.3e3), ("cond_out 21C", 294.15, 1878e3),
         ("cond_out 26C", 299.15, 1878e3), ("cond_out 14C", 287.15, 1878e3)]
solver = get_solver("ipopt", options={"max_iter": 500})
print(f"{'case':14s}{'VLE':9s}{'start':15s}{'term':16s}{'T C':>9s}{'true':>8s}{'vap':>8s}  verdict")
for name, T, P in cases:
    h_true, _, ph = probe.state(T, P); vf_true = 1.0 if ph == "Vap" else 0.0
    for vle_name, vle in [("Smooth", SmoothVLE), ("CubicCC", CubicComplementarityVLE)]:
        for mode in ("seed", "seed+bound", "default+bound"):
            cfg = make_config(METHODS["NIST"]); cfg["phase_equilibrium_state"] = {("Vap", "Liq"): vle}
            cfg["state_definition"] = FPhx
            cfg["state_bounds"] = {"flow_mol": (0, 1, 1000, pyunits.mol/pyunits.s),
                "enth_mol": (-6e4, -5e3, 6e4, pyunits.J/pyunits.mol),
                "temperature": (180, 300, 500, pyunits.K), "pressure": (1e4, 1e5, 1e7, pyunits.Pa)}
            m = ConcreteModel(); m.fs = FlowsheetBlock(dynamic=False)
            m.fs.pp = GenericParameterBlock(**cfg)
            m.fs.sb = m.fs.pp.build_state_block([0], defined_state=True); b = m.fs.sb[0]
            b.flow_mol.fix(1); b.mole_frac_comp["R32"].fix(1); b.pressure.fix(P); b.enth_mol.fix(h_true)
            try: m.fs.sb.initialize()
            except Exception: pass
            if mode.startswith("seed"):
                b.temperature.set_value(T)
                b.phase_frac["Vap"].set_value(vf_true); b.phase_frac["Liq"].set_value(1 - vf_true)
                b.flow_mol_phase["Vap"].set_value(vf_true); b.flow_mol_phase["Liq"].set_value(1 - vf_true)
            if "bound" in mode:
                if vf_true > 0.5: b.phase_frac["Vap"].setlb(0.999)
                else: b.phase_frac["Vap"].setub(0.001)
            try: tc = str(solver.solve(m).solver.termination_condition)
            except Exception as e: tc = "EXC"
            Tm = value(b.temperature); vf = value(b.phase_frac["Vap"])
            print(f"{name:14s}{vle_name:9s}{mode:15s}{tc:16s}{Tm-273.15:9.2f}{T-273.15:8.2f}{vf:8.4f}  {'ok' if abs(Tm-T) < 0.05 else 'WRONG'}")
