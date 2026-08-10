"""
phase4_full_pfd_comparison.py -- full PFD-style state table (T, P, h, s,
vapor quality) at all four cycle stream points, for all four methods
(Helmholtz baseline, NIST, GCGP, SPGP), at every ambient in the sweep.

Extends phase4_property_state_comparison.py (which only reported T/h at
evap_out/comp_isen/comp_out/cond_out for Helmholtz/NIST/GCGP) by adding:
  - the expansion valve outlet (the 4th stream point -- previously missing)
  - pressure and entropy at every point, not just T and h
  - vapor fraction (quality) at every point
  - SPGP, which is expected to fail/show large residuals -- included
    explicitly rather than skipped, since the point of this run is to show
    ALL cases considered, including the ones that don't converge cleanly

All cubic-PR properties are molar (enth_mol in J/mol, entr_mol in J/mol/K);
converted to mass basis (M_R32 = 0.052024 kg/mol) to compare directly
against Helmholtz's native mass-basis enth_mass/entr_mass.

Run this on the machine where idaes's compiled cubic-root extension is
installed (myidaesenv) -- it will NOT run in a sandbox without that
extension (RuntimeError: "Cubic root external functions are not
available").

Usage: python3 phase4_full_pfd_comparison.py
"""
from pyomo.environ import value, Var
from vapor_compression_cubic import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
from vapor_compression_plr import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]
M_R32 = 0.052024  # kg/mol
H_MAX = 700e3      # J/kg -- same enth_mass relaxation phase_3a_helmholtz_cop.py uses

STREAM_POINTS = ["evap_out", "comp_out", "cond_out", "valve_out"]


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
    return {
        "cop": cop, "converged": converged,
        "ratioP": value(m.fs.compressor.ratioP[0]),
        "comp_isen_T": value(m.fs.compressor.properties_isentropic[0].temperature),
        "comp_isen_h": value(m.fs.compressor.properties_isentropic[0].enth_mol) / M_R32,
        "points": points,
    }


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
    return {
        "cop": cop, "converged": converged,
        "ratioP": value(m.fs.compressor.outlet.pressure[0]) / value(m.fs.compressor.inlet.pressure[0]),
        "comp_isen_T": value(m.fs.compressor.properties_isentropic[0].temperature),
        "comp_isen_h": value(m.fs.compressor.properties_isentropic[0].enth_mass),
        "points": points,
    }


METHODS = ["Helmholtz", "NIST", "GCGP", "SPGP"]
results = {m: {} for m in METHODS}

for Tamb in AMBIENTS:
    print(f"\n{'='*70}\n  Running T_amb = {Tamb} C\n{'='*70}")
    for method in ("NIST", "GCGP", "SPGP"):
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

# --- Summary: COP, pressure ratio, isentropic reference block ---
print(f"\n{'='*100}\n  Summary: COP, pressure ratio, isentropic reference (T,h)\n{'='*100}")
print(f"{'T_amb':>7}{'method':>12}{'COP':>9}{'ratioP':>9}{'T_isen(K)':>11}{'h_isen(kJ/kg)':>15}")
for Tamb in AMBIENTS:
    for method in METHODS:
        r = results[method].get(Tamb, {})
        if not r.get("converged"):
            err = r.get("error", "not converged")
            print(f"{Tamb:>7}{method:>12}  FAILED: {err}")
            continue
        print(f"{Tamb:>7}{method:>12}{r['cop']:>9.4f}{r['ratioP']:>9.4f}"
              f"{r['comp_isen_T']:>11.3f}{r['comp_isen_h']/1e3:>15.2f}")

# --- Full PFD table: one block per (T_amb, method), all 4 stream points ---
print(f"\n{'='*100}\n  Full PFD state table: T(K), P(bar), h(kJ/kg), s(kJ/kg-K), quality\n{'='*100}")
for Tamb in AMBIENTS:
    for method in METHODS:
        r = results[method].get(Tamb, {})
        print(f"\n--- T_amb={Tamb}C, {method} ---")
        if not r.get("converged"):
            print(f"  NOT CONVERGED: {r.get('error', 'unknown')}")
            continue
        print(f"  {'stream':>10}{'T(K)':>10}{'P(bar)':>10}{'h(kJ/kg)':>11}{'s(kJ/kg-K)':>12}{'quality':>9}")
        for sp in STREAM_POINTS:
            pt = r["points"][sp]
            print(f"  {sp:>10}{pt['T']:>10.3f}{pt['P']/1e5:>10.4f}{pt['h']/1e3:>11.3f}{pt['s']/1e3:>12.4f}{pt['x']:>9.4f}")
