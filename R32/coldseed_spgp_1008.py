"""
coldseed_spgp_1008.py -- SPGP (fp_single_0919): solve each ambient from a cold
start, then re-solve warm-started from that converged cold-start state.

Same cycle, solve call and warmstart_from_state as algseed_test.py; only the
seed differs (own converged cold state instead of the algebraic points).

Writes coldseed_spgp_1008.csv and cop_vs_ambient_coldseed_1008.png.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""
import logging, warnings, io, contextlib, csv
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyomo.environ import value
from vapor_compression_cubic_refstate_1001 import SimpleVaporCompressionCycle as C, Mode

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.9999
AMB = [10, 15, 20, 25]
METHOD = "fp_single_0919"
UNITS = ("evap_out", "comp_out", "cond_out", "valve_out")

def _units(model):
    return dict(zip(UNITS, (model.fs.evaporator, model.fs.compressor,
                            model.fs.condenser, model.fs.expansion_valve)))

def warmstart_from_state(model, pts):
    for s, u in _units(model).items():
        st = u.control_volume.properties_out[0]
        st.temperature.set_value(pts[s]["T"]); st.pressure.set_value(pts[s]["P"])

def grab(model):
    return {s: {"T": value(u.control_volume.properties_out[0].temperature),
                "P": value(u.control_volume.properties_out[0].pressure)}
            for s, u in _units(model).items()}

def run(Tamb, seed):
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc = C('R32', compressor_efficiency=ETA, mode=Mode.IMPROVED_TPX, method=METHOD)
            vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                          high_side_temperature=Tamb+APPROACH)
            vc.initialize(verbose=False)
            if seed is not None:
                warmstart_from_state(vc.model, seed)
            vc.set_specifications(ambient_temperature=Tamb, condenser_approach=APPROACH,
                                  evap_sat_temperature=T_EVAP, superheating=0,
                                  subcooling=0, max_pressure_ratio=10)
            cop, conv = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
            pts = grab(vc.model) if conv else None
        return cop, bool(conv), pts
    except Exception as e:
        return None, False, None

POSTER = [3.739092, 3.549206, 3.184143, 3.041566]   # phase6_0924 chase result
HELM   = [3.948947, 3.526194, 3.190784, 2.914391]
ALGSEED= [3.4878, 3.2201, 2.9714, 2.8036]          # algseed_test.py, 10/08

rows = []
for Tamb in AMB:
    c0, v0, p0 = run(Tamb, None)
    c1, v1 = (None, False)
    if v0:
        c1, v1, _ = run(Tamb, p0)
    rows.append(dict(ambient_C=Tamb, cold_cop=c0, cold_converged=v0,
                     coldseed_cop=c1, coldseed_converged=v1))
    f = lambda c: "None" if c is None else f"{c:.6f}"
    print(f"Tamb={Tamb:>3}  cold {f(c0)} conv={v0}   seeded-from-cold {f(c1)} conv={v1}", flush=True)

with open("coldseed_spgp_1008.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def arr(key, convkey):
    return np.array([r[key] if (r[convkey] and r[key] is not None) else np.nan for r in rows], float)
cold, seeded = arr("cold_cop", "cold_converged"), arr("coldseed_cop", "coldseed_converged")
cv = lambda y: np.round(np.diff(np.diff(y)), 4)
print("\ncurvature  poster", cv(POSTER), " cold", cv(cold), " seeded-from-cold", cv(seeded),
      " algseed", cv(ALGSEED))

amb = np.array(AMB, float)
fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.plot(amb, HELM, "-", color="k", marker="o", ms=5, lw=2.4, label="Helmholtz EoS / CoolProp")
ax.plot(amb, POSTER, "--", color="tab:red", marker="D", ms=4.5, lw=1.3, mfc="none", mew=1.1,
        label="PR EoS / SPGP, phase6 chase (poster)")
ok = np.isfinite(cold)
ax.plot(amb[ok], cold[ok], ":", color="tab:orange", marker="s", ms=4.5, lw=1.3, mfc="none",
        mew=1.1, label="PR EoS / SPGP, cold start")
ok = np.isfinite(seeded)
ax.plot(amb[ok], seeded[ok], "-", color="tab:red", marker="D", ms=4.5, lw=1.6,
        label="PR EoS / SPGP, seeded from converged cold start")
ax.set_xlabel("Ambient temperature [deg C]"); ax.set_ylabel("COP"); ax.set_xticks(amb)
ax.set_xlim(amb.min()-1, amb.max()+1); ax.grid(True, ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
fig.tight_layout(); fig.savefig("cop_vs_ambient_coldseed_1008.png", dpi=200)
print("wrote coldseed_spgp_1008.csv, cop_vs_ambient_coldseed_1008.png")
