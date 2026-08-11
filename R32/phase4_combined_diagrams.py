"""
phase4_combined_diagrams.py -- ONE p-h diagram, ONE p-T diagram, and ONE T-s
diagram, each combining all 4 ambients (10/15/20/25 C) and all 3 converged
methods (Helmholtz, NIST, GCGP -- SPGP excluded, does not converge at any
ambient, see BREADCRUMB_07-20.md 2026-08-10 update) on a single set of axes.

Runs the REAL cycle model (vapor_compression_cubic.py / vapor_compression_plr.py)
directly -- this must be run in myidaesenv (the environment with the compiled
cubic-root extension), NOT in a generic sandbox. Produces 3 PNG files in the
same directory as this script.

Reference-state note (see chat discussion 2026-08-10): Helmholtz and the real
CoolProp dome share the identical reference state (both built from
Tillner-Roth & Yokozeki 1997 for R32 -- IIR convention, h=200 kJ/kg,
s=1.00 kJ/kg-K at 0 C sat. liquid), so in principle Helmholtz needs NO shift
against the dome. NIST and GCGP share a *different* common reference with
EACH OTHER (both pinned to ideal-gas h=0 at T_ref=298.15K in the IDAES
generic-property config), but that reference is not the IIR one, so they
need a large, expected shift against the dome. Per Shilpa's explicit
instruction, this script does NOT change the anchoring approach -- it uses
the same one constant-offset-per-method scheme as before (computed once at
T_amb=20 by matching that method's own evap_out property to the real fluid's
saturated-vapor property at the shared -29 C setpoint), applied identically
to h and s. The offsets are printed so any drift between methods is visible,
not hidden.

Usage: python3 phase4_combined_diagrams.py   (run in myidaesenv)
"""
import os
from pyomo.environ import value, Var
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from vapor_compression_cubic import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
from vapor_compression_plr import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]
M_R32 = 0.052024  # kg/mol
H_MAX = 700e3
T_EVAP_SET_K = 244.15  # -29 C, shared evaporator setpoint
STREAM_POINTS = ["evap_out", "comp_out", "cond_out", "valve_out"]
OUTDIR = os.path.dirname(os.path.abspath(__file__))


def relax_enth_bounds(vc, hmax=H_MAX):
    for v in vc.model.component_data_objects(Var, active=True, descend_into=True):
        if v.parent_component().local_name == "enth_mass":
            if v.ub is not None and v.ub < hmax:
                v.setub(hmax)


def _cubic_point(state):
    return {
        "T": value(state.temperature),
        "P": value(state.pressure),
        "h": value(state.enth_mol) / M_R32,
        "s": value(state.entr_mol) / M_R32,
        "x": value(state.phase_frac["Vap"]),
    }


def _helm_point(state):
    return {
        "T": value(state.temperature),
        "P": value(state.pressure),
        "h": value(state.enth_mass),
        "s": value(state.entr_mass),
        "x": value(state.phase_frac["Vap"]),
    }


def run_cubic(method, Tamb):
    Tcond_sat = Tamb + 9
    vc = CubicCycle(FLUID, compressor_efficiency=0.9999, mode=CubicMode.IMPROVED_TPX, method=method)
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)
    vc.set_specifications(ambient_temperature=Tamb, condenser_approach=9,
                           evap_sat_temperature=-29, superheating=0, subcooling=0,
                           max_pressure_ratio=10)
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    m = vc.model
    points = {
        "evap_out": _cubic_point(m.fs.evaporator.control_volume.properties_out[0]),
        "comp_out": _cubic_point(m.fs.compressor.control_volume.properties_out[0]),
        "cond_out": _cubic_point(m.fs.condenser.control_volume.properties_out[0]),
        "valve_out": _cubic_point(m.fs.expansion_valve.control_volume.properties_out[0]),
    }
    return {"cop": cop, "converged": converged, "points": points}


def run_helm(Tamb):
    Tcond_sat = Tamb + 9
    vc = HelmCycle(FLUID, compressor_efficiency=0.9999, mode=HelmMode.IMPROVED_TPX)
    relax_enth_bounds(vc)
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)
    vc.set_specifications(ambient_temperature=Tamb, condenser_approach=9,
                           evap_sat_temperature=-29, superheating=0, subcooling=0,
                           max_pressure_ratio=10)
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    m = vc.model
    points = {
        "evap_out": _helm_point(m.fs.evaporator.control_volume.properties_out[0]),
        "comp_out": _helm_point(m.fs.compressor.control_volume.properties_out[0]),
        "cond_out": _helm_point(m.fs.condenser.control_volume.properties_out[0]),
        "valve_out": _helm_point(m.fs.expansion_valve.control_volume.properties_out[0]),
    }
    return {"cop": cop, "converged": converged, "points": points}


METHODS = ["Helmholtz", "NIST", "GCGP"]
results = {m: {} for m in METHODS}

for Tamb in AMBIENTS:
    print(f"\n{'='*60}\n  T_amb = {Tamb} C\n{'='*60}")
    for method in ("NIST", "GCGP"):
        try:
            results[method][Tamb] = run_cubic(method, Tamb)
            print(f"  {method}: converged={results[method][Tamb]['converged']}, COP={results[method][Tamb]['cop']:.4f}")
        except Exception as e:
            results[method][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
            print(f"  {method}: FAILED -- {type(e).__name__}: {e}")
    try:
        results["Helmholtz"][Tamb] = run_helm(Tamb)
        print(f"  Helmholtz: converged={results['Helmholtz'][Tamb]['converged']}, COP={results['Helmholtz'][Tamb]['cop']:.4f}")
    except Exception as e:
        results["Helmholtz"][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
        print(f"  Helmholtz: FAILED -- {type(e).__name__}: {e}")

# --- real dome from CoolProp ---
Tc = CP.PropsSI("Tcrit", FLUID)
Ttrip = CP.PropsSI("Ttriple", FLUID)
T_dome = np.linspace(Ttrip + 0.5, Tc - 0.05, 400)
h_liq = np.array([CP.PropsSI("H", "T", T, "Q", 0, FLUID) for T in T_dome]) / 1e3
h_vap = np.array([CP.PropsSI("H", "T", T, "Q", 1, FLUID) for T in T_dome]) / 1e3
s_liq = np.array([CP.PropsSI("S", "T", T, "Q", 0, FLUID) for T in T_dome]) / 1e3
s_vap = np.array([CP.PropsSI("S", "T", T, "Q", 1, FLUID) for T in T_dome]) / 1e3
P_liq = np.array([CP.PropsSI("P", "T", T, "Q", 0, FLUID) for T in T_dome]) / 1e5
P_vap = np.array([CP.PropsSI("P", "T", T, "Q", 1, FLUID) for T in T_dome]) / 1e5
T_dome_C = T_dome - 273.15

h_vap_real_at_evap = CP.PropsSI("H", "T", T_EVAP_SET_K, "Q", 1, FLUID) / 1e3
s_vap_real_at_evap = CP.PropsSI("S", "T", T_EVAP_SET_K, "Q", 1, FLUID) / 1e3

# --- per-method constant offsets, anchored at T_amb=20 evap_out (same scheme as before) ---
h_offsets, s_offsets = {}, {}
for m in METHODS:
    r20 = results[m].get(20, {})
    if r20.get("converged"):
        h_offsets[m] = h_vap_real_at_evap - r20["points"]["evap_out"]["h"] / 1e3
        s_offsets[m] = s_vap_real_at_evap - r20["points"]["evap_out"]["s"] / 1e3
    else:
        h_offsets[m] = 0.0
        s_offsets[m] = 0.0
print("\nEnthalpy offsets (kJ/kg):", h_offsets)
print("Entropy offsets (kJ/kg-K):", s_offsets)

COLORS = {"Helmholtz": "black", "NIST": "tab:blue", "GCGP": "tab:red"}
LINESTYLES = {10: "dotted", 15: "dashed", 20: "solid", 25: (0, (3, 1, 1, 1))}


def _loop(method, Tamb, key, offset=0.0, scale=1e-3):
    r = results[method].get(Tamb, {})
    if not r.get("converged"):
        return None, None
    pts = r["points"]
    vals = [pts[sp][key] * scale + offset for sp in STREAM_POINTS]
    vals.append(vals[0])
    return vals


# ============================== Figure 1: p-h ==============================
fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(h_liq, P_liq, color="tab:blue", lw=1.5, label="R32 sat. liquid (CoolProp)")
ax.plot(h_vap, P_vap, color="tab:orange", lw=1.5, label="R32 sat. vapor (CoolProp)")
for method in METHODS:
    for Tamb in AMBIENTS:
        hs = _loop(method, Tamb, "h", offset=h_offsets[method])
        r = results[method].get(Tamb, {})
        if hs is None or not r.get("converged"):
            continue
        Ps = [r["points"][sp]["P"] / 1e5 for sp in STREAM_POINTS]
        Ps.append(Ps[0])
        lbl = f"{method}, T_amb={Tamb}C" if Tamb == 20 else None
        ax.plot(hs, Ps, "o-", color=COLORS[method], ls=LINESTYLES[Tamb],
                 lw=1.8 if Tamb == 20 else 1.2, ms=4, alpha=1.0 if Tamb == 20 else 0.6, label=lbl)
ax.set_yscale("log")
ax.set_xlabel("Enthalpy, h (kJ/kg) -- shifted per method to common real-fluid basis")
ax.set_ylabel("Pressure, P (bar, log scale)")
ax.set_title("R32 p-h diagram: real dome + all methods, all ambients")
ax.legend(loc="lower right", fontsize=8)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "combined_ph_diagram.png"), dpi=150)
plt.close(fig)

# ============================== Figure 2: p-T ==============================
fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(T_dome_C, P_liq, color="tab:blue", lw=1.5, label="R32 sat. liquid (CoolProp)")
ax.plot(T_dome_C, P_vap, color="tab:orange", lw=1.5, label="R32 sat. vapor (CoolProp)")
for method in METHODS:
    for Tamb in AMBIENTS:
        r = results[method].get(Tamb, {})
        if not r.get("converged"):
            continue
        Ts = [r["points"][sp]["T"] - 273.15 for sp in STREAM_POINTS]
        Ts.append(Ts[0])
        Ps = [r["points"][sp]["P"] / 1e5 for sp in STREAM_POINTS]
        Ps.append(Ps[0])
        lbl = f"{method}, T_amb={Tamb}C" if Tamb == 20 else None
        ax.plot(Ts, Ps, "o-", color=COLORS[method], ls=LINESTYLES[Tamb],
                 lw=1.8 if Tamb == 20 else 1.2, ms=4, alpha=1.0 if Tamb == 20 else 0.6, label=lbl)
ax.set_yscale("log")
ax.set_xlabel("Temperature, T (C) -- NOT shifted (native reported values, incl. any known artifacts)")
ax.set_ylabel("Pressure, P (bar, log scale)")
ax.set_title("R32 p-T diagram: real dome + all methods, all ambients")
ax.legend(loc="upper left", fontsize=8)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "combined_pT_diagram.png"), dpi=150)
plt.close(fig)

# ============================== Figure 3: T-s ==============================
fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(s_liq, T_dome_C, color="tab:blue", lw=1.5, label="R32 sat. liquid (CoolProp)")
ax.plot(s_vap, T_dome_C, color="tab:orange", lw=1.5, label="R32 sat. vapor (CoolProp)")
for method in METHODS:
    for Tamb in AMBIENTS:
        ss = _loop(method, Tamb, "s", offset=s_offsets[method])
        r = results[method].get(Tamb, {})
        if ss is None or not r.get("converged"):
            continue
        Ts = [r["points"][sp]["T"] - 273.15 for sp in STREAM_POINTS]
        Ts.append(Ts[0])
        lbl = f"{method}, T_amb={Tamb}C" if Tamb == 20 else None
        ax.plot(ss, Ts, "o-", color=COLORS[method], ls=LINESTYLES[Tamb],
                 lw=1.8 if Tamb == 20 else 1.2, ms=4, alpha=1.0 if Tamb == 20 else 0.6, label=lbl)
ax.set_xlabel("Entropy, s (kJ/kg-K) -- shifted per method to common real-fluid basis")
ax.set_ylabel("Temperature, T (C) -- NOT shifted (native reported values, incl. any known artifacts)")
ax.set_title("R32 T-s diagram: real dome + all methods, all ambients")
ax.legend(loc="upper left", fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "combined_Ts_diagram.png"), dpi=150)
plt.close(fig)

print("\nSaved:")
print(" ", os.path.join(OUTDIR, "combined_ph_diagram.png"))
print(" ", os.path.join(OUTDIR, "combined_pT_diagram.png"))
print(" ", os.path.join(OUTDIR, "combined_Ts_diagram.png"))
