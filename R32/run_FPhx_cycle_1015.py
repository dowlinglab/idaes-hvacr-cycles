"""
run_FPhx_cycle_1015.py -- build, seed, solve the FPhx cycle at one ambient.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    Builds the FPhx cycle (NIST), prints DOF, seeds every block, solves,
    bounds phases (bound_phases, bound_two_phase_T),
    prints a per-block phase check and cycle_guards.check_cycle.
Why
    First solve of the FPhx cycle.
Basis
    vapor_compression_cubic_FPhx_1015.py
Outputs
    printed. As of 10/15 night: NIST 20 C solves, COP 2.6149 = closed form,
    all 8 blocks physical, cycle_guards pass.
Run
    python run_FPhx_cycle_1015.py 20 NIST 2>&1 | tee run_FPhx_cycle_1015.log
"""
import sys
from pyomo.environ import value
from idaes.core.util.model_statistics import degrees_of_freedom
import vapor_compression_cubic_FPhx_1015 as V
from cycle_guards import check_cycle
Tamb=float(sys.argv[1]) if len(sys.argv)>1 else 20.0
method=sys.argv[2] if len(sys.argv)>2 else "NIST"
m, sb, info, probe = V.build_cycle(method, T_amb_C=Tamb)
r=info["ref"]
rep=V.seed_bounded(m, sb, info, probe, eps=1e-5)
V.bound_two_phase_T(m, sb, probe)
print("seed (block-wise flash at fixed P,h):")
for k,(tc,T,vf) in rep.items(): print(f"  {k:16s} {tc:12s} T {T:8.2f} C  vap {vf:.6f}")
print("DOF at solve:", degrees_of_freedom(m))
tc = V.solve(m)
fs=m.fs
Q=value(fs.evaporator.heat_duty[0]); W=value(fs.compressor.work_mechanical[0]); K=value(fs.condenser.heat_duty[0])
cop=Q/W
print("termination:", tc, " COP %.4f  closed form %.4f  diff %.4f %%"%(cop, r["cop"], 100*(cop-r["cop"])/r["cop"]))
print(f"{'block':10s}{'T C':>9s}{'PR Tsat':>9s}{'vap':>9s}{'h':>11s}  check")
for k,T,Ts,vf,v in V.phase_check(sb, probe):
    print(f"{k:10s}{T:9.2f}{Ts:9.2f}{vf:9.4f}{value(sb[k].enth_mol):11.1f}  {v}")
rpt = check_cycle(m, cop, T_evap_C=-29.0, T_cond_C=info["T_cond_C"], fluid=None, strict=False)
print(rpt.get("summary", rpt))
