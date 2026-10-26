"""
scaling_check_FPhx_1015.py -- is "infeasible" only a scaling/tolerance artifact?

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    One cycle (method, T_amb), set up as sweep_FPhx_1015.py square, solved
    three ways:
      raw      tol 1e-8, no user scaling        (what the sweep did)
      tol1e-6  IDAES default tolerance
      scaled   energy balances x 1e-5, enthalpies x 1e-4, flows x 1e-1,
               IPOPT nlp_scaling_method = user-scaling, tol 1e-8
    Prints IPOPT result, COP vs closed form, largest RELATIVE residual
    (residual / largest term in that equation), phase check.
Why
    On the Mac, NIST 10 C ends "infeasible" with raw residuals of 3.5e-3 W on
    ~2.9e5 W balances (~1e-8 relative). See breadcrumb 10/15 night (5).
    This tests whether IPOPT is judging big equations in absolute units.
    It does NOT remove the two-phase degeneracy; that is option B.
Basis
    vapor_compression_cubic_FPhx_1015.py
Outputs
    printed.
Run
    python scaling_check_FPhx_1015.py NIST 10 2>&1 | tee scaling_NIST_10.log
"""
import sys, logging, warnings
logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
from pyomo.environ import value, Constraint, Suffix
from pyomo.core.expr.visitor import identify_variables
from idaes.core.solvers import get_solver
import vapor_compression_cubic_FPhx_1015 as V

KEYS = {"NIST": "NIST", "GCGP": "GCGP_scaled", "SPGP": "fp_single_0919"}
label = sys.argv[1] if len(sys.argv) > 1 else "NIST"
T_amb = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0


def setup():
    m, sb, info, probe = V.build_cycle(KEYS[label], T_amb_C=T_amb)
    V.seed_bounded(m, sb, info, probe, eps=1e-5)
    V.bound_two_phase_T(m, sb, probe, dT=0.01)
    return m, sb, info, probe


def add_scaling(m):
    """Suffix-based user scaling: only the large-magnitude quantities."""
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)
    fs = m.fs
    for u in (fs.evaporator, fs.compressor, fs.condenser, fs.expansion_valve):
        for c in u.control_volume.enthalpy_balances.values():
            m.scaling_factor[c] = 1e-5
        for t in (u.control_volume.properties_in, u.control_volume.properties_out):
            b = t[0]
            m.scaling_factor[b.enth_mol] = 1e-4
            m.scaling_factor[b.flow_mol] = 1e-1
            m.scaling_factor[b.pressure] = 1e-5
    for k in ("heat_duty",):
        for u in (fs.evaporator, fs.condenser):
            m.scaling_factor[getattr(u, k)[0]] = 1e-5
    m.scaling_factor[fs.compressor.work_mechanical[0]] = 1e-5


def rel_resid(c):
    body = value(c.body)
    r = 0.0
    if c.lower is not None: r = max(r, value(c.lower) - body)
    if c.upper is not None: r = max(r, body - value(c.upper))
    size = max([1.0] + [abs(value(v)) for v in identify_variables(c.body) if v.value is not None])
    return r / size, r


print(f"=== {label}  T_amb {T_amb:.0f} C ===")
print(f"{'setup':>9s}{'IPOPT':>13s}{'COP':>10s}{'closed':>10s}{'max rel resid':>15s}{'(raw)':>11s}  phases  worst equation")
for name in ("raw", "tol1e-6", "scaled"):
    m, sb, info, probe = setup()
    opts = {"max_iter": 500, "tol": 1e-8}
    if name == "tol1e-6": opts["tol"] = 1e-6
    if name == "scaled":
        add_scaling(m); opts["nlp_scaling_method"] = "user-scaling"
    solver = get_solver("ipopt", options=opts)
    try: tc = str(solver.solve(m).solver.termination_condition)
    except Exception as e: tc = "EXC"
    fs = m.fs
    cop = value(fs.evaporator.heat_duty[0]) / value(fs.compressor.work_mechanical[0])
    worst = max(((rel_resid(c), c.name) for c in m.component_data_objects(Constraint, active=True)),
                key=lambda x: x[0][0])
    ok = all(v == "ok" for *_, v in V.phase_check(sb, probe))
    short = ".".join(worst[1].split(".")[1:])[-55:]
    print(f"{name:>9s}{tc:>13s}{cop:10.5f}{info['ref']['cop']:10.5f}{worst[0][0]:15.2e}{worst[0][1]:11.2e}  {str(ok):6s}  {short}")
