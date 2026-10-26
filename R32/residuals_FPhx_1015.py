"""
residuals_FPhx_1015.py -- why does IPOPT say "infeasible" when the COP matches?

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    One cycle (method, T_amb), set up exactly as sweep_FPhx_1015.py square.
    Prints: IPOPT result, the 10 largest equation residuals, every variable
    sitting on a bound. Then repeats with the two-phase temperature bound
    loosened (0.01 K, 0.1 K, off) to see whether that bound is the cause.
Why
    5 of 12 square-sweep points end "infeasible" while matching the closed
    form and passing every check (10/15 night).
Basis
    vapor_compression_cubic_FPhx_1015.py
Outputs
    printed.
Run
    python residuals_FPhx_1015.py NIST 10 2>&1 | tee residuals_NIST_10.log
"""
import sys, logging, warnings
logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
from pyomo.environ import value, Constraint, Var
import vapor_compression_cubic_FPhx_1015 as V

KEYS = {"NIST": "NIST", "GCGP": "GCGP_scaled", "SPGP": "fp_single_0919"}
label = sys.argv[1] if len(sys.argv) > 1 else "NIST"
T_amb = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0


def residual(c):
    body = value(c.body)
    r = 0.0
    if c.lower is not None: r = max(r, value(c.lower) - body)
    if c.upper is not None: r = max(r, body - value(c.upper))
    return r


def run(dT):
    m, sb, info, probe = V.build_cycle(KEYS[label], T_amb_C=T_amb)
    V.seed_bounded(m, sb, info, probe, eps=1e-5)
    if dT is not None:
        V.bound_two_phase_T(m, sb, probe, dT=dT)
    tc = V.solve(m)
    fs = m.fs
    cop = value(fs.evaporator.heat_duty[0]) / value(fs.compressor.work_mechanical[0])
    ref = info["ref"]["cop"]
    T2ph = [value(sb[k].temperature) - 273.15 for k in ("valve_out", "evap_in")]
    return m, sb, tc, cop, ref, T2ph, probe


print(f"=== {label}  T_amb {T_amb:.0f} C ===\n")
m, sb, tc, cop, ref, T2ph, probe = run(0.01)
print(f"IPOPT: {tc}   COP {cop:.5f}   closed form {ref:.5f}\n")

res = []
for c in m.component_data_objects(Constraint, active=True):
    try: res.append((residual(c), c.name))
    except Exception: res.append((float("nan"), c.name))
res.sort(key=lambda x: -x[0] if x[0] == x[0] else 1e300)
print("10 largest equation residuals (raw units of each equation):")
for r, n in res[:10]:
    print(f"  {r:11.3e}  {n}")

print("\nvariables on a bound (within 1e-8 relative):")
n_at = 0
for v in m.component_data_objects(Var):
    if v.fixed or v.value is None: continue
    for side, b in (("lb", v.lb), ("ub", v.ub)):
        if b is not None and abs(v.value - b) <= 1e-8 * max(1.0, abs(b)):
            print(f"  {side} {b:12.6g}  {v.name}"); n_at += 1
if n_at == 0: print("  none")

print("\ntwo-phase temperature bound loosened:")
print(f"  {'bound':>8s}{'IPOPT':>13s}{'COP':>10s}{'closed':>10s}{'T valve_out':>13s}{'T evap_in':>11s}  phases")
for dT in (0.01, 0.1, None):
    m2, sb2, tc2, cop2, ref2, T2, pr2 = run(dT)
    ok = all(v == "ok" for *_, v in V.phase_check(sb2, pr2))
    lab = "off" if dT is None else f"{dT} K"
    print(f"  {lab:>8s}{tc2:>13s}{cop2:10.5f}{ref2:10.5f}{T2[0]:13.3f}{T2[1]:11.3f}  {ok}")
