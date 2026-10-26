"""
sweep_FPhx_1015.py -- R-32 COP vs ambient on the IDAES PR (FPhx) cycle, all
three PR methods, with every state checked.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    For NIST, GCGP (GCGP_scaled), SPGP (fp_single_0919) at T_amb 10/15/20/25:
    build FPhx cycle, seed (seed_bounded), bound phases, solve.
    mode "square":   SH 3 K, SC 8 K fixed.
    mode "optimize": SH, SC free in 3..8 K, maximise COP (2 DOF).
    Then: phase check on all 8 blocks, cycle_guards.check_cycle, and the
    closed-form COP from the same PR equations.
Why
    First IDAES PR sweep where every state is physical (10/15 night).
Basis
    vapor_compression_cubic_FPhx_1015.py; spec T_evap -29 C, T_cond T_amb+9,
    eta 0.80, 1 kg/s.
Outputs
    sweep_FPhx_<mode>_1015.csv
Run
    python sweep_FPhx_1015.py square   2>&1 | tee sweep_FPhx_square_1015.log
    python sweep_FPhx_1015.py optimize 2>&1 | tee sweep_FPhx_optimize_1015.log
"""
import sys, csv, logging, warnings
logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
from pyomo.environ import value, Objective, Constraint, maximize
from idaes.core.util.model_statistics import degrees_of_freedom
import vapor_compression_cubic_FPhx_1015 as V
from cycle_guards import check_cycle

MODE = sys.argv[1] if len(sys.argv) > 1 else "square"
METHODS = [("NIST", "NIST"), ("GCGP", "GCGP_scaled"), ("SPGP", "fp_single_0919")]
AMB = [10.0, 15.0, 20.0, 25.0]
rows = []
print(f"{'method':6s}{'Tamb':>5s}{'DOF':>4s}{'term':>11s}{'COP':>8s}{'closed':>8s}{'diff%':>9s}"
      f"{'SH':>6s}{'SC':>6s}{'loopW':>8s}{'phases':>8s}{'guards':>8s}")
for label, key in METHODS:
    for a in AMB:
        row = dict(method=label, T_amb=a)
        try:
            m, sb, info, probe = V.build_cycle(key, T_amb_C=a)
            V.seed_bounded(m, sb, info, probe, eps=1e-5)
            V.bound_two_phase_T(m, sb, probe)
            fs = m.fs
            if MODE == "optimize":
                fs.T_evap_out.unfix(); fs.T_cond_out.unfix()
                fs.T_evap_out.setlb(fs.Tsat_evap + 3); fs.T_evap_out.setub(fs.Tsat_evap + 8)
                fs.T_cond_out.setlb(fs.Tsat_cond - 8); fs.T_cond_out.setub(fs.Tsat_cond - 3)
                fs.obj = Objective(expr=fs.evaporator.heat_duty[0] / fs.compressor.work_mechanical[0],
                                   sense=maximize)
            dof = degrees_of_freedom(m)
            tc = V.solve(m)
            Q = value(fs.evaporator.heat_duty[0]); W = value(fs.compressor.work_mechanical[0])
            cop = Q / W
            SH = value(fs.T_evap_out) - fs.Tsat_evap; SC = fs.Tsat_cond - value(fs.T_cond_out)
            ref = V.closed_form(probe, -29.0, info["T_cond_C"], SH, SC, 0.80)["cop"]
            pc = V.phase_check(sb, probe)
            ph_ok = all(v == "ok" for *_, v in pc)
            rpt = check_cycle(m, cop, T_evap_C=-29.0, T_cond_C=info["T_cond_C"], fluid=None, strict=False)
            loop = Q + W + value(fs.condenser.heat_duty[0])
            row.update(dof=dof, term=tc, cop=cop, closed_form=ref, diff_pct=100*(cop-ref)/ref,
                       SH=SH, SC=SC, loop_W=loop, phases_ok=ph_ok, guards_ok=rpt["ok"],
                       carnot=(244.15)/(info["T_cond_C"]+29.0),
                       T_dis_C=value(sb["comp_out"].temperature)-273.15)
            if not ph_ok:
                row["bad_blocks"] = ";".join(k for k, *_, v in pc if v != "ok")
            print(f"{label:6s}{a:5.0f}{dof:4d}{tc:>11s}{cop:8.4f}{ref:8.4f}{row['diff_pct']:9.4f}"
                  f"{SH:6.2f}{SC:6.2f}{loop:8.2f}{str(ph_ok):>8s}{str(rpt['ok']):>8s}")
        except Exception as e:
            row.update(term="EXC " + type(e).__name__ + ": " + str(e)[:60])
            print(f"{label:6s}{a:5.0f}  FAILED {row['term']}")
        rows.append(row)
keys = ["method","T_amb","dof","term","cop","closed_form","diff_pct","carnot","SH","SC",
        "loop_W","T_dis_C","phases_ok","guards_ok","bad_blocks"]
with open(f"sweep_FPhx_{MODE}_1015.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
