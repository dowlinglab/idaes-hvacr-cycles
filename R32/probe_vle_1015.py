"""
probe_vle_1015.py -- does an isolated PR state block land on the right phase?

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    Solves ONE state block at a known state (T,P or P,h) for four setups:
    SmoothVLE / CubicComplementarityVLE  x  FTPx / FPhx. Compares T, vapour
    fraction and h with the true state from the fugacity-sign closed form.
Why
    To find out whether wrong cycle states come from the cycle or from the
    property block itself.
Basis
    Truth: PRProbe in vapor_compression_cubic_FPhx_1015.py (no solver).
Outputs
    table printed to screen. NOTE: "WRONG" on comp_out is a false alarm --
    above Tc there is one root, so vapour fraction has no meaning there.
Run
    cd ~/Desktop/DVCT_Project/property_packages_DVRT_code/idaes-hvacr-cycles/R32
    conda activate myidaesenv
    python probe_vle_1015.py 2>&1 | tee probe_vle_1015.log
"""
import logging; logging.disable(logging.WARNING)
import warnings; warnings.filterwarnings("ignore")
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.state_definitions import FTPx, FPhx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.smooth_VLE_2 import CubicComplementarityVLE
from pyomo.environ import units as pyunits
import vapor_compression_cubic_FPhx_1015 as V
from phase_1_cubic_eos_validation_refstate_1001 import make_config, METHODS
probe = V.PRProbe("NIST")
cases = [("evap_out SH3", 247.15, 283.3e3), ("cond_out 21C", 294.15, 1878e3),
         ("cond_out 26C", 299.15, 1878e3), ("cond_out 14C", 287.15, 1878e3),
         ("comp_out", 387.4, 1892.6e3)]
solver = get_solver("ipopt", options={"max_iter": 500})
print(f"{'case':14s}{'VLE':10s}{'state':6s}{'term':12s}{'T C':>9s}{'vap':>8s}{'h':>10s}{'h true':>10s}  verdict")
for name, T, P in cases:
    h_true, _, ph = probe.state(T, P); vf_true = 1.0 if ph == "Vap" else 0.0
    for vle_name, vle in [("Smooth", SmoothVLE), ("CubicCC", CubicComplementarityVLE)]:
        for sd_name, sd in [("FTPx", FTPx), ("FPhx", FPhx)]:
            cfg = make_config(METHODS["NIST"]); cfg["phase_equilibrium_state"] = {("Vap","Liq"): vle}
            cfg["state_definition"] = sd
            if sd is FPhx:
                cfg["state_bounds"] = {"flow_mol": (0, 1, 1000, pyunits.mol/pyunits.s),
                    "enth_mol": (-6e4, -5e3, 6e4, pyunits.J/pyunits.mol),
                    "temperature": (180, 300, 500, pyunits.K), "pressure": (1e4, 1e5, 1e7, pyunits.Pa)}
            m = ConcreteModel(); m.fs = FlowsheetBlock(dynamic=False)
            m.fs.pp = GenericParameterBlock(**cfg)
            m.fs.sb = m.fs.pp.build_state_block([0], defined_state=True); b = m.fs.sb[0]
            b.flow_mol.fix(1); b.mole_frac_comp["R32"].fix(1); b.pressure.fix(P)
            if sd is FTPx: b.temperature.fix(T)
            else: b.enth_mol.fix(h_true)
            try:
                m.fs.sb.initialize()
                tc = str(solver.solve(m).solver.termination_condition)
            except Exception as e:
                tc = "EXC " + type(e).__name__
            Tm = value(b.temperature); vf = value(b.phase_frac["Vap"]); hm = value(b.enth_mol)
            ok = abs(vf - vf_true) < 0.01 and abs(Tm - T) < 0.05 and abs(hm - h_true) < 5
            print(f"{name:14s}{vle_name:10s}{sd_name:6s}{tc:12s}{Tm-273.15:9.2f}{vf:8.4f}{hm:10.1f}{h_true:10.1f}  {'ok' if ok else 'WRONG'}")
