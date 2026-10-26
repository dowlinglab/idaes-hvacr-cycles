"""
nist_fix_1008.py -- why NIST's seeded-from-cold re-solve gave 2.2837 at 15 C,
and the fix.

In coldseed_1008.py the re-solve seeds only T and P on the four outlets
(8 numbers). Everything else (flows, phase fractions, enthalpies) comes from
a fresh vc.initialize(), so the starting point is not the converged state.
For NIST at 15 C that mixed start converged to a different solution
(2.2837 vs 3.3974).

This script runs NIST at each ambient three ways:
    cold        cold start, unchanged from coldseed_1008.py
    tp_seed     T,P-only seed from the converged cold state (reproduces 2.2837)
    full_seed   EVERY free variable copied from the converged cold model,
                applied after set_specifications so specs are never touched

and prints T, P, vapour fraction and h at all four state points for each,
so a phase flip or other change in the 15 C TP-seeded solution is visible.

If the cold solution is a genuine converged point, full_seed must return the
cold COP. If it does not, the problem has more than one solution near there.

Writes nist_fix_1008.csv and cop_vs_ambient_nist_fix_1008.png.
Same cycle, solver call and conditions as coldseed_1008.py.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""
import logging, warnings, io, contextlib, csv
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyomo.environ import value, Var
from vapor_compression_cubic_refstate_1001 import SimpleVaporCompressionCycle as C, Mode

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.9999
AMB = [10, 15, 20, 25]
METHOD = "NIST"
UNITS = ("evap_out", "comp_out", "cond_out", "valve_out")
HELM = [3.948947, 3.526194, 3.190784, 2.914391]


def _units(model):
    return dict(zip(UNITS, (model.fs.evaporator, model.fs.compressor,
                            model.fs.condenser, model.fs.expansion_valve)))


def _v(x):
    try:
        return value(x)
    except Exception:
        return float("nan")


def states(model):
    out = {}
    for s, u in _units(model).items():
        st = u.control_volume.properties_out[0]
        out[s] = {"T": _v(st.temperature), "P": _v(st.pressure),
                  "x": _v(st.phase_frac["Vap"]), "h": _v(st.enth_mol)}
    return out


def tp_seed(model, src):
    """The coldseed_1008 seed: T and P only, applied BEFORE set_specifications."""
    for s, u in _units(model).items():
        st = u.control_volume.properties_out[0]
        st.temperature.set_value(src[s]["T"]); st.pressure.set_value(src[s]["P"])


def full_seed(model, src_model):
    """Copy every free variable's value from the converged model.
    Applied AFTER set_specifications, and fixed variables are skipped, so the
    seed can only change the starting guess, never a specification."""
    n = 0
    for v in model.component_data_objects(Var, descend_into=True):
        if v.fixed:
            continue
        src = src_model.find_component(v.name)
        if src is None or src.value is None:
            continue
        v.set_value(src.value, skip_validation=True)
        n += 1
    return n


def run(Tamb, mode=None, src=None):
    """mode: None (cold), 'tp' (src = states dict), 'full' (src = converged model)."""
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc = C('R32', compressor_efficiency=ETA, mode=Mode.IMPROVED_TPX, method=METHOD)
            vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                          high_side_temperature=Tamb+APPROACH)
            vc.initialize(verbose=False)
            if mode == "tp":
                tp_seed(vc.model, src)
            vc.set_specifications(ambient_temperature=Tamb, condenser_approach=APPROACH,
                                  evap_sat_temperature=T_EVAP, superheating=0,
                                  subcooling=0, max_pressure_ratio=10)
            if mode == "full":
                full_seed(vc.model, src)
            cop, conv = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
        return cop, bool(conv), vc
    except Exception as e:
        print(f"    exception: {type(e).__name__}: {e}")
        return None, False, None


def show(tag, cop, conv, st):
    c = "None" if cop is None else f"{cop:.6f}"
    print(f"  {tag:<10} COP {c}  conv={conv}")
    if st is None:
        return
    for s in UNITS:
        d = st[s]
        print(f"      {s:<10} T {d['T']-273.15:8.3f} C   P {d['P']/1e5:8.4f} bar   "
              f"x_vap {d['x']:8.5f}   h {d['h']:11.2f} J/mol")


rows = []
res = {"cold": [], "tp_seed": [], "full_seed": []}
for Tamb in AMB:
    print(f"\n=== NIST  T_amb = {Tamb} C ===")
    c0, v0, vc0 = run(Tamb)
    s0 = states(vc0.model) if vc0 is not None else None
    show("cold", c0, v0, s0)
    out = {"cold": (c0, v0)}
    if v0:
        c1, v1, vc1 = run(Tamb, "tp", s0)
        show("tp_seed", c1, v1, states(vc1.model) if vc1 is not None else None)
        c2, v2, vc2 = run(Tamb, "full", vc0.model)
        show("full_seed", c2, v2, states(vc2.model) if vc2 is not None else None)
        out["tp_seed"], out["full_seed"] = (c1, v1), (c2, v2)
    else:
        out["tp_seed"] = out["full_seed"] = (None, False)
    for k, (c, v) in out.items():
        res[k].append(c if (v and c is not None) else np.nan)
        rows.append(dict(method=METHOD, ambient_C=Tamb, run=k, cop=c, converged=v))

with open("nist_fix_1008.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

cv = lambda y: np.round(np.diff(np.diff(np.asarray(y, float))), 4)
print("\nSummary")
print(f"{'run':<10}" + "".join(f"{a:>10}" for a in AMB) + "   curvature")
for k, y in res.items():
    print(f"{k:<10}" + "".join(f"{x:>10.4f}" if np.isfinite(x) else f"{'--':>10}" for x in y)
          + f"   {cv(y)}")

amb = np.array(AMB, float)
fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.plot(amb, HELM, "-", color="k", marker="o", ms=5, lw=2.4, label="Helmholtz EoS / CoolProp")
for k, ls, mfc, lab in (("cold", ":", "none", "PR EoS / NIST, cold start"),
                        ("tp_seed", "--", "none", "PR EoS / NIST, T,P seed from cold"),
                        ("full_seed", "-", "tab:blue", "PR EoS / NIST, full-state seed from cold")):
    y = np.asarray(res[k], float); ok = np.isfinite(y)
    ax.plot(amb[ok], y[ok], ls, color="tab:blue", marker="s", ms=4.5, lw=1.3,
            mfc=mfc, mew=1.1, label=lab)
ax.set_xlabel("Ambient temperature [deg C]"); ax.set_ylabel("COP"); ax.set_xticks(amb)
ax.set_xlim(amb.min()-1, amb.max()+1); ax.grid(True, ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
fig.tight_layout(); fig.savefig("cop_vs_ambient_nist_fix_1008.png", dpi=200)
print("\nwrote nist_fix_1008.csv, cop_vs_ambient_nist_fix_1008.png")
