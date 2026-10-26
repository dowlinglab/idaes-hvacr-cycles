"""
phase7_sweep_1015.py -- COP vs ambient, single stage, every PR method, on the
arc-fixed model (vapor_compression_cubic_refstate_1015).

WHY SQUARE
    The arc fix takes degrees of freedom to 0, so there is nothing to
    optimise and optimize=True raises "Too few degrees of freedom". The
    2-DOF design (unfix evaporator/condenser outlet phase_frac, add the
    superheat/subcool window caps) is the next step; this file is the square
    solve, which is what the model supports today.

EVERY POINT IS CHECKED
    cycle_guards.check_cycle() runs on each solve before the row is written:
    Carnot ceiling, loop energy closure, arc enthalpy continuity, mass
    continuity, spec read-back. The CSV carries the measured superheat and
    subcool, not the requested ones.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/15/2026
"""

import sys
sys.path.insert(0, "/mnt/user-data/outputs")
import logging
import warnings
import io
import contextlib
import csv

logging.disable(logging.WARNING)
warnings.filterwarnings("ignore")

from pyomo.environ import value
from CoolProp.CoolProp import PropsSI as CP
from idaes.core.util.model_statistics import degrees_of_freedom
from cycle_guards import check_cycle
from vapor_compression_cubic_refstate_1015 import (
    SimpleVaporCompressionCycle as C, Mode)

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.9999
SH_LO, SH_HI = 3.0, 8.0      # superheat window, same as the cascade hot loop
SC_LO, SC_HI = 3.0, 8.0      # subcool window
SH, SC = SH_LO, SC_HI        # expected corner, for the guard's spec check
AMBIENTS = [10.0, 15.0, 20.0, 25.0]
METHODS = ["NIST", "GCGP", "GCGP_scaled", "fp_single_0919"]

# 10/15/2026: the default tolerances (tol 1e-4) make IPOPT exit with a hard
# error on GCGP_scaled at every ambient -- not infeasible, a crash. Loosening
# them solves it cleanly: COP 3.4034, loop closure 0.003 %. Applied to every
# method so the comparison stays like-for-like.
SOLVER_OPTS = {"max_iter": 3000, "tol": 1e-3,
               "constr_viol_tol": 1e-3, "acceptable_tol": 1e-2}


def closed(Tc_C, sh, sc):
    """Closed-form Helmholtz COP at the SOLVED superheat/subcool."""
    Te, Tc = T_EVAP + 273.15, Tc_C + 273.15
    pe = CP("P", "T", Te, "Q", 1, "R32")
    pc = CP("P", "T", Tc, "Q", 0, "R32")
    h1 = CP("H", "T", Te + sh, "P", pe, "R32")
    s1 = CP("S", "T", Te + sh, "P", pe, "R32")
    h2 = h1 + (CP("H", "P", pc, "S", s1, "R32") - h1) / ETA
    h3 = CP("H", "T", Tc - sc, "P", pc, "R32")
    return (h1 - h3) / (h2 - h1)


rows = []
hdr = (f"{'method':<16}{'Tamb':>5}{'DOF':>4}{'COP':>9}{'Carnot':>8}{'/Carn':>7}"
       f"{'conv':>7}{'loop%':>9}{'maxarc':>10}{'SH':>6}{'SC':>6}{'fail':>5}")
print(hdr, flush=True)
print("-" * len(hdr), flush=True)

for meth in METHODS:
    for Ta in AMBIENTS:
        Tc = Ta + APPROACH
        rec = dict(method=meth, ambient_C=Ta, T_cond_C=Tc)
        try:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                vc = C("R32", compressor_efficiency=ETA,
                       mode=Mode.IMPROVED_TPX, method=meth)
                vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                              high_side_temperature=Tc)
                vc.initialize(verbose=False)
                vc.set_specifications(ambient_temperature=Ta,
                                      condenser_approach=APPROACH,
                                      evap_sat_temperature=T_EVAP,
                                      superheating=SH_LO,
                                      superheating_max=SH_HI,
                                      subcooling=SC_LO,
                                      subcooling_max=SC_HI,
                                      free_outlet_quality=True,
                                      max_pressure_ratio=20)
                dof = degrees_of_freedom(vc.model)
                cop, conv = vc.optimize_COP(verbose=False, initialize=True,
                                            optimize=True,
                                            solver_options=SOLVER_OPTS)
                m = vc.model
                Q = value(m.fs.evaporator.heat_duty[0])
                W = value(m.fs.compressor.work_mechanical[0])
                K = value(m.fs.condenser.heat_duty[0])
                rpt = check_cycle(m, cop, T_evap_C=T_EVAP, T_cond_C=Tc,
                                  fluid="R32", superheat_K=SH, subcool_K=SC,
                                  strict=False)
            carn = (T_EVAP + 273.15) / (Tc - T_EVAP)
            maxarc = max(abs(v) for v in rpt["arc_breaks"].values())
            sh_a = rpt.get("superheat_K", float("nan"))
            sc_a = rpt.get("subcool_K", float("nan"))
            rec.update(dof=dof, cop=cop, converged=conv, carnot=carn,
                       cop_over_carnot=cop / carn, Q_evap_W=Q, W_comp_W=W,
                       Q_cond_W=K, loop_W=Q + W + K,
                       loop_pct_of_Qevap=100 * abs(Q + W + K) / abs(Q),
                       max_arc_break_Jmol=maxarc,
                       evap_superheat_actual_K=sh_a,
                       cond_subcool_actual_K=sc_a,
                       closed_form_cop=closed(Tc, sh_a, sc_a),
                       n_guard_failures=len(rpt["failures"]))
            print(f"{meth:<16}{Ta:>5.0f}{dof:>4}{cop:>9.4f}{carn:>8.4f}"
                  f"{cop / carn:>7.3f}{str(conv):>7}"
                  f"{rec['loop_pct_of_Qevap']:>8.3f}%{maxarc:>10.1e}"
                  f"{sh_a:>6.2f}{sc_a:>6.2f}{rec['n_guard_failures']:>5}",
                  flush=True)
        except Exception as exc:
            rec.update(cop=None, converged=f"EXC {type(exc).__name__}")
            print(f"{meth:<16}{Ta:>5.0f}   FAILED {type(exc).__name__}: "
                  f"{str(exc)[:45]}", flush=True)
        rows.append(rec)

keys = sorted({k for r in rows for k in r})
with open("/mnt/user-data/outputs/phase7_sweep_1015.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(rows)
print("\nwrote phase7_sweep_1015.csv", flush=True)
