"""
debug_nist_in_sweep.py -- Test #2 of the history-dependency check (see
BREADCRUMB_07-20.md, and debug_nist_isolated.py for Test #1).

Runs the FULL method/ambient sweep in the EXACT SAME ORDER as
phase4_full_pfd_comparison.py (all methods x all ambients), so NIST@T10 is
solved with exactly the same preceding call history it would have had in
that original run. Prints ONLY the NIST@T_amb=15 result at the end.

Compare this script's NIST@T15 output directly against three separate runs
of `python3 debug_nist_isolated.py` (which solves NIST@T10 completely
fresh, nothing run before it):
  - If they all match -> NIST@T15 is deterministic regardless of context;
    history/state leakage is ruled out as the explanation for the earlier
    324.571K vs 310.175K discrepancy.
  - If this in-sweep result differs from the isolated result -> some
    solver/model state is leaking between iterations (a fixable ordering
    bug), which would need to be ruled out before crediting "uncertainty
    propagation" through the EOS parameters themselves.

    Author: Claude AI
Testing and Validation: Shilpa Narasimhan

Date created: 08/11/2026
"""
from pyomo.environ import value, Var
from vapor_compression_cubic_refstate import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode
from vapor_compression_plr import SimpleVaporCompressionCycle as HelmCycle, Mode as HelmMode

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]
M_R32 = 0.052024
H_MAX = 700e3
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


METHODS = ["NIST", "GCGP", "SPGP"]
results = {m: {} for m in METHODS}
results["Helmholtz"] = {}

for Tamb in AMBIENTS:
    print(f"[running] T_amb={Tamb}C ...")
    for method in ("NIST", "GCGP", "SPGP"):
        try:
            results[method][Tamb] = run_cubic(method, Tamb)
        except Exception as e:
            results[method][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}
    try:
        results["Helmholtz"][Tamb] = run_helm(Tamb)
    except Exception as e:
        results["Helmholtz"][Tamb] = {"converged": False, "error": f"{type(e).__name__}: {e}"}

r = results["NIST"][15]
print(f"\nIN-SWEEP RESULT -- NIST @ T_amb=15C -- "
      f"converged={r.get('converged')}, COP={r.get('cop', float('nan')):.4f}")
print(f"  {'stream':>10}{'T(K)':>10}{'P(bar)':>10}{'h(kJ/kg)':>11}{'s(kJ/kg-K)':>12}{'quality':>9}")
for sp in STREAM_POINTS:
    pt = r["points"][sp]
    print(f"  {sp:>10}{pt['T']:>10.3f}{pt['P']/1e5:>10.4f}{pt['h']/1e3:>11.3f}{pt['s']/1e3:>12.4f}{pt['x']:>9.4f}")
