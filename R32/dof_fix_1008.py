"""
dof_fix_1008.py -- close the two free temperatures and check the answer no
longer depends on the starting point.

ROOT CAUSE (from nist_fix_1008.py output, 10/08)
    With superheating=0 and subcooling=0, set_specifications() in
    vapor_compression_cubic_refstate_1001.py:
        - fixes evap_out phase_frac[Vap] = 1 and DEACTIVATES superheating_constraint
        - fixes cond_out phase_frac[Vap] = 0 and DEACTIVATES subcooling_constraint
    The pressures are pinned (temperature_bubble == T_evap and == T_amb+approach),
    but nothing pins evap_out.temperature or cond_out.temperature. They are only
    bounded to +/-10 K of saturation. With optimize=False the objective is off,
    so IPOPT solves a feasibility problem with two free variables and returns
    whichever feasible point is nearest the starting guess. Evidence, same P in
    every run:
        cond_out T at 20 C ambient: 29.31 / 34.80 / 38.50 C, all x_vap = 0
        evap_out T at 15 C ambient: -21.9 to -25.5 C, all x_vap = 1
    COP follows those free temperatures, hence the seed dependence.

THE FIX (applied here, model file untouched)
    After set_specifications, add two equalities:
        evap_out.T == temperature_bubble + SH      (SH = 0: saturated vapour)
        cond_out.T == temperature_bubble - SC      (SC = 0: saturated liquid)
    and check degrees_of_freedom == 0 before solving.

    If exact saturation (SH = SC = 0) does not converge (phase_frac fixed at a
    dome boundary is ill-conditioned, see the 07-23 comments in the model),
    set OFFSET_K to 0.1 to step just off the boundary.

TEST
    Each method x ambient solved cold, then T,P-seeded and full-state-seeded
    from that cold solution. If the fix works, all three agree.

Writes dof_fix_1008.csv and cop_vs_ambient_dof_fix_1008.png.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""
import logging, warnings, io, contextlib, csv, sys
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyomo.environ import value, Var, Constraint
from idaes.core.util.model_statistics import degrees_of_freedom
from vapor_compression_cubic_refstate_1001 import SimpleVaporCompressionCycle as C, Mode

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.9999
AMB = [10, 15, 20, 25]
OFFSET_K = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0   # SH = SC = OFFSET_K
METHODS = {   # label: (config name, colour, marker, poster values)
    "NIST": ("NIST",           "tab:blue",  "s", [3.567653, 3.397422, 3.131608, 2.800516]),
    "GCGP": ("GCGP",           "tab:green", "^", [3.509379, 3.258122, 3.054898, 2.890954]),
    "SPGP": ("fp_single_0919", "tab:red",   "D", [3.739092, 3.549206, 3.184143, 3.041566]),
}
HELM = [3.948947, 3.526194, 3.190784, 2.914391]
UNITS = ("evap_out", "comp_out", "cond_out", "valve_out")


def _units(model):
    return dict(zip(UNITS, (model.fs.evaporator, model.fs.compressor,
                            model.fs.condenser, model.fs.expansion_valve)))


def pin_saturation(model):
    ev = model.fs.evaporator.control_volume.properties_out[0]
    cd = model.fs.condenser.control_volume.properties_out[0]
    model.fs.evap_out_T_pin = Constraint(
        expr=ev.temperature == ev.temperature_bubble["Vap", "Liq"] + OFFSET_K)
    model.fs.cond_out_T_pin = Constraint(
        expr=cd.temperature == cd.temperature_bubble["Vap", "Liq"] - OFFSET_K)


def tp_seed(model, src):
    for s, u in _units(model).items():
        st = u.control_volume.properties_out[0]
        st.temperature.set_value(src[s]["T"]); st.pressure.set_value(src[s]["P"])


def full_seed(model, src_model):
    for v in model.component_data_objects(Var, descend_into=True):
        if v.fixed:
            continue
        src = src_model.find_component(v.name)
        if src is not None and src.value is not None:
            v.set_value(src.value, skip_validation=True)


def states(model):
    return {s: {"T": value(u.control_volume.properties_out[0].temperature),
                "P": value(u.control_volume.properties_out[0].pressure)}
            for s, u in _units(model).items()}


def run(cfg, Tamb, mode=None, src=None):
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc = C('R32', compressor_efficiency=ETA, mode=Mode.IMPROVED_TPX, method=cfg)
            vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                          high_side_temperature=Tamb+APPROACH)
            vc.initialize(verbose=False)
            if mode == "tp":
                tp_seed(vc.model, src)
            vc.set_specifications(ambient_temperature=Tamb, condenser_approach=APPROACH,
                                  evap_sat_temperature=T_EVAP, superheating=0,
                                  subcooling=0, max_pressure_ratio=10)
            pin_saturation(vc.model)
            dof = degrees_of_freedom(vc.model)
            if mode == "full":
                full_seed(vc.model, src)
            cop, conv = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
        return cop, bool(conv), vc, dof
    except Exception as e:
        print(f"    exception: {type(e).__name__}: {e}")
        return None, False, None, None


f = lambda c: "  None  " if c is None else f"{c:.6f}"
cv = lambda y: np.round(np.diff(np.diff(np.asarray(y, float))), 4)
print(f"saturation pin offset = {OFFSET_K} K\n")

rows, res = [], {}
for lab, (cfg, *_) in METHODS.items():
    res[lab] = {"cold": [], "tp_seed": [], "full_seed": []}
    for Tamb in AMB:
        c0, v0, vc0, dof = run(cfg, Tamb)
        out = {"cold": (c0, v0)}
        if v0:
            c1, v1, *_ = run(cfg, Tamb, "tp", states(vc0.model))
            c2, v2, *_ = run(cfg, Tamb, "full", vc0.model)
            out["tp_seed"], out["full_seed"] = (c1, v1), (c2, v2)
            ev = vc0.model.fs.evaporator.control_volume.properties_out[0]
            cd = vc0.model.fs.condenser.control_volume.properties_out[0]
            chk = (f"evap_out T-Tsat {value(ev.temperature)-value(ev.temperature_bubble['Vap','Liq']):+.4f} K  "
                   f"cond_out T-Tsat {value(cd.temperature)-value(cd.temperature_bubble['Vap','Liq']):+.4f} K")
        else:
            out["tp_seed"] = out["full_seed"] = (None, False)
            chk = ""
        print(f"{lab:<5} Tamb={Tamb:>3}  DOF={dof}  cold {f(c0)} {v0!s:<5}  "
              f"tp {f(out['tp_seed'][0])} {out['tp_seed'][1]!s:<5}  "
              f"full {f(out['full_seed'][0])} {out['full_seed'][1]!s:<5}  {chk}", flush=True)
        for k, (c, v) in out.items():
            res[lab][k].append(c if (v and c is not None) else np.nan)
            rows.append(dict(method=lab, config=cfg, ambient_C=Tamb, offset_K=OFFSET_K,
                             dof=dof, run=k, cop=c, converged=v))

with open("dof_fix_1008.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

print("\nSummary (cold start, pinned)        spread across cold/tp/full      curvature")
print(f"  Helmholtz poster {np.round(HELM, 4)}                                   {cv(HELM)}")
for lab, (_, _, _, poster) in METHODS.items():
    y = np.array(res[lab]["cold"])
    allr = np.vstack([res[lab][k] for k in ("cold", "tp_seed", "full_seed")])
    spread = np.nanmax(allr, 0) - np.nanmin(allr, 0)
    print(f"  {lab:<5} {np.round(y, 4)}   {np.round(spread, 6)}   {cv(y)}")
    print(f"        poster {np.round(poster, 4)}")

amb = np.array(AMB, float)
fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.plot(amb, HELM, "-", color="k", marker="o", ms=5, lw=2.4, label="Helmholtz EoS / CoolProp")
for lab, (_, col, mk, _) in METHODS.items():
    y = np.array(res[lab]["cold"]); ok = np.isfinite(y)
    ax.plot(amb[ok], y[ok], "--", color=col, marker=mk, ms=4.5, lw=1.3, mfc="none", mew=1.1,
            label=f"PR EoS / {lab}")
ax.set_xlabel("Ambient temperature [deg C]"); ax.set_ylabel("COP"); ax.set_xticks(amb)
ax.set_xlim(amb.min()-1, amb.max()+1); ax.grid(True, ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
fig.tight_layout(); fig.savefig("cop_vs_ambient_dof_fix_1008.png", dpi=200)
print("\nwrote dof_fix_1008.csv, cop_vs_ambient_dof_fix_1008.png")
