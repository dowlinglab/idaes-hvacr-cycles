"""
sweep_eq_PR_1015.py -- R-32 COP vs ambient, equation-oriented PR cycle
(cycle_eq_PR_1015.py), three PR methods, square and optimised.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    NIST, GCGP (GCGP_scaled), SPGP (fp_single_0919) x T_amb 10/15/20/25.
      square:   SH 3 K, SC 8 K fixed (DOF 0)
      optimize: SH, SC free in 3..8 K, maximise COP (DOF 2), run twice:
                started at SH 8 / SC 3 and at SH 3 / SC 8. Two starts because
                for SPGP COP vs SH has a MINIMUM near 6 K (closed form, SC 8,
                T_amb 10: 3.205768 @3, 3.205427 @6, 3.205689 @8), so both
                ends are local maxima.
    Checks per point: IPOPT result, COP vs closed form, phase of every
    state against PR's own saturation line, valve outlet 0 < x < 1, below
    Carnot, compressor s2 >= s1.
Why
    Option B (10/15): no vapour-fraction variable, so no fake states.
Basis
    T_evap -29 C, T_cond = T_amb + 9, eta 0.80.
Outputs
    sweep_eq_PR_1015.csv
Run
    python sweep_eq_PR_1015.py 2>&1 | tee sweep_eq_PR_1015.log
"""
import csv, logging, warnings
logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
import cycle_eq_PR_1015 as C

METHODS = [("NIST", "NIST"), ("GCGP", "GCGP_scaled"), ("SPGP", "fp_single_0919")]
AMB = [10.0, 15.0, 20.0, 25.0]
rows = []
print(f"{'method':6s}{'start':>6s}{'mode':>9s}{'Tamb':>5s}{'DOF':>4s}{'IPOPT':>11s}{'COP':>9s}{'closed':>9s}"
      f"{'diff%':>9s}{'SH':>6s}{'SC':>6s}{'x4':>7s}{'T2 C':>7s}  checks")
for mode, start in (("square", (3.0, 8.0)), ("optimize", (8.0, 3.0)), ("optimize", (3.0, 8.0))):
    for label, key in METHODS:
        for a in AMB:
            try:
                r = C.solve_point(key, a, mode, SH=start[0], SC=start[1])
                r["method"] = label; r["start"] = f"{start[0]:.0f}/{start[1]:.0f}"
                ok = r["phases_ok"] and r["below_carnot"] and r["second_law_ok"]
                print(f"{label:6s}{r['start']:>6s}{mode:>9s}{a:5.0f}{r['dof']:4d}{r['term']:>11s}{r['cop']:9.4f}"
                      f"{r['closed_form']:9.4f}{r['diff_pct']:9.4f}{r['SH']:6.2f}{r['SC']:6.2f}"
                      f"{r['x4']:7.3f}{r['T2_C']:7.1f}  {'ok' if ok else 'FAIL'}")
            except Exception as e:
                r = dict(method=label, T_amb=a, mode=mode, term="EXC " + str(e)[:60])
                print(f"{label:6s}{mode:>9s}{a:5.0f}  FAILED {r['term']}")
            rows.append(r)
keys = ["method", "mode", "start", "T_amb", "dof", "term", "cop", "closed_form", "diff_pct", "carnot",
        "SH", "SC", "x4", "T2_C", "Pe_kPa", "Pc_kPa", "phases_ok", "below_carnot", "second_law_ok"]
with open("sweep_eq_PR_1015.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
