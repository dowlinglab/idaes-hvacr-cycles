"""
debug_nist_isolated.py -- Test #1 of the history-dependency check (see
BREADCRUMB_07-20.md). Builds and solves ONLY NIST @ T_amb=15C, nothing
else -- a completely fresh process, no other method/ambient solved before
it. Run this as `python3 debug_nist_isolated.py` THREE SEPARATE TIMES
(three separate invocations, not a loop inside one process) and compare
the printed comp_out/cond_out T,h,s across the three runs.

If all three isolated runs agree with each other, NIST@T15 is
deterministic in isolation -- the next question is whether it ALSO
matches what debug_nist_in_sweep.py gets for the same point when solved
as part of the full method/ambient sweep (same order as
phase4_full_pfd_comparison.py). If isolated-run and in-sweep results
differ, that's history/state leakage between iterations -- a fixable
bug, not "uncertainty propagation" through the EOS itself.

Author: Claude AI
Testing and Validation: Shilpa Narasimhan

Date created: 08/11/2026
"""
from pyomo.environ import value
from vapor_compression_cubic_refstate import SimpleVaporCompressionCycle as CubicCycle, Mode as CubicMode

FLUID = "R32"
M_R32 = 0.052024
STREAM_POINTS = ["evap_out", "comp_out", "cond_out", "valve_out"]


def _cubic_point(state):
    return {
        "T": value(state.temperature),
        "P": value(state.pressure),
        "h": value(state.enth_mol) / M_R32,
        "s": value(state.entr_mol) / M_R32,
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


if __name__ == "__main__":
    r = run_cubic("NIST", 15)
    print(f"\nISOLATED RUN -- NIST @ T_amb=15C -- converged={r['converged']}, COP={r['cop']:.4f}")
    print(f"  {'stream':>10}{'T(K)':>10}{'P(bar)':>10}{'h(kJ/kg)':>11}{'s(kJ/kg-K)':>12}{'quality':>9}")
    for sp in STREAM_POINTS:
        pt = r["points"][sp]
        print(f"  {sp:>10}{pt['T']:>10.3f}{pt['P']/1e5:>10.4f}{pt['h']/1e3:>11.3f}{pt['s']/1e3:>12.4f}{pt['x']:>9.4f}")
