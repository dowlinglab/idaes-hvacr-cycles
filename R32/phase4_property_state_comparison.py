"""
phase4_property_state_comparison.py -- direct answer to "is something
happening to the PREDICTED PROPERTIES at T_amb=20 specifically, or is
this purely a compressor-solve numerics bug?" (task #37 follow-up).

Everything diagnosed so far (compressor_fix_regression_check.py,
phase4_gcgp_sweep.py --diagnose) has looked at the SYMPTOM: the
compressor's real outlet landing on a wrong branch relative to its own
isentropic reference. This script instead compares the actual predicted
STATE properties (T, P, h at each of the 4 cycle streams, plus pressure
ratio) side by side for GCGP, NIST, and the Phase 3a Helmholtz baseline,
at every ambient in the sweep -- to see whether GCGP's own property
predictions diverge from NIST/Helmholtz specifically at T_amb=20 (e.g. a
notably different pressure ratio, or a jump/kink other ambients don't
show), or whether the properties track smoothly across all three methods
and the compressor bug really is isolated numerics.

All cubic-PR enthalpies are molar (J/mol); converted to J/kg (M_R32 =
0.052024 kg/mol) so they're directly comparable to Helmholtz's native
mass-basis enth_mass.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-27
"""
from pyomo.environ import value
from vapor_compression_cubic import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
from vapor_compression_plr import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode
from pyomo.environ import Var

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]
M_R32 = 0.052024  # kg/mol
H_MAX = 700e3      # J/kg -- same enth_mass relaxation phase_3a_helmholtz_cop.py uses


def relax_enth_bounds(vc, hmax=H_MAX):
    for v in vc.model.component_data_objects(Var, active=True, descend_into=True):
        if v.parent_component().local_name == "enth_mass":
            if v.ub is not None and v.ub < hmax:
                v.setub(hmax)


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
    e_out = m.fs.evaporator.control_volume.properties_out[0]
    comp_out = m.fs.compressor.control_volume.properties_out[0]
    comp_isen = m.fs.compressor.properties_isentropic[0]
    c_out = m.fs.condenser.control_volume.properties_out[0]
    return {
        "cop": cop, "converged": converged,
        "ratioP": value(m.fs.compressor.ratioP[0]),
        "P_low": value(m.fs.compressor.inlet.pressure[0]),
        "P_high": value(m.fs.compressor.outlet.pressure[0]),
        "evap_out_T": value(e_out.temperature), "evap_out_h": value(e_out.enth_mol) / M_R32,
        "comp_isen_T": value(comp_isen.temperature), "comp_isen_h": value(comp_isen.enth_mol) / M_R32,
        "comp_out_T": value(comp_out.temperature), "comp_out_h": value(comp_out.enth_mol) / M_R32,
        "cond_out_T": value(c_out.temperature), "cond_out_h": value(c_out.enth_mol) / M_R32,
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
    return {
        "cop": cop, "converged": converged,
        "ratioP": value(m.fs.compressor.outlet.pressure[0]) / value(m.fs.compressor.inlet.pressure[0]),
        "P_low": value(m.fs.compressor.inlet.pressure[0]),
        "P_high": value(m.fs.compressor.outlet.pressure[0]),
        "evap_out_T": value(m.fs.evaporator.outlet.temperature[0]), "evap_out_h": value(m.fs.evaporator.outlet.enth_mass[0]),
        "comp_isen_T": value(m.fs.compressor.properties_isentropic[0].temperature), "comp_isen_h": value(m.fs.compressor.properties_isentropic[0].enth_mass),
        "comp_out_T": value(m.fs.compressor.outlet.temperature[0]), "comp_out_h": value(m.fs.compressor.outlet.enth_mass[0]),
        "cond_out_T": value(m.fs.condenser.outlet.temperature[0]), "cond_out_h": value(m.fs.condenser.outlet.enth_mass[0]),
    }


results = {"NIST": {}, "GCGP": {}, "Helmholtz": {}}
for Tamb in AMBIENTS:
    print(f"\n{'='*70}\n  T_amb = {Tamb} C\n{'='*70}")
    for method in ("NIST", "GCGP"):
        try:
            results[method][Tamb] = run_cubic(method, Tamb)
        except Exception as e:
            results[method][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
    try:
        results["Helmholtz"][Tamb] = run_helm(Tamb)
    except Exception as e:
        results["Helmholtz"][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}

# --- pressure ratio + P_low/P_high table ---
print(f"\n{'='*100}\n  Pressure ratio and suction/discharge pressure by method\n{'='*100}")
print(f"{'T_amb':>7}{'method':>12}{'P_low(bar)':>12}{'P_high(bar)':>13}{'ratioP':>9}")
for Tamb in AMBIENTS:
    for method in ("Helmholtz", "NIST", "GCGP"):
        r = results[method].get(Tamb, {})
        if not r.get("converged"):
            print(f"{Tamb:>7}{method:>12}{'FAILED':>12}")
            continue
        print(f"{Tamb:>7}{method:>12}{r['P_low']/1e5:>12.4f}{r['P_high']/1e5:>13.4f}{r['ratioP']:>9.4f}")

# --- state property table: T (K) at each stream ---
print(f"\n{'='*100}\n  Temperature (K) at each stream, by method\n{'='*100}")
print(f"{'T_amb':>7}{'method':>12}{'evap_out':>11}{'comp_isen':>11}{'comp_out':>11}{'cond_out':>11}")
for Tamb in AMBIENTS:
    for method in ("Helmholtz", "NIST", "GCGP"):
        r = results[method].get(Tamb, {})
        if not r.get("converged"):
            continue
        print(f"{Tamb:>7}{method:>12}{r['evap_out_T']:>11.3f}{r['comp_isen_T']:>11.3f}"
              f"{r['comp_out_T']:>11.3f}{r['cond_out_T']:>11.3f}")

# --- state property table: h (J/kg) at each stream ---
print(f"\n{'='*100}\n  Enthalpy (kJ/kg) at each stream, by method\n{'='*100}")
print(f"{'T_amb':>7}{'method':>12}{'evap_out':>11}{'comp_isen':>11}{'comp_out':>11}{'cond_out':>11}")
for Tamb in AMBIENTS:
    for method in ("Helmholtz", "NIST", "GCGP"):
        r = results[method].get(Tamb, {})
        if not r.get("converged"):
            continue
        print(f"{Tamb:>7}{method:>12}{r['evap_out_h']/1e3:>11.2f}{r['comp_isen_h']/1e3:>11.2f}"
              f"{r['comp_out_h']/1e3:>11.2f}{r['cond_out_h']/1e3:>11.2f}")
