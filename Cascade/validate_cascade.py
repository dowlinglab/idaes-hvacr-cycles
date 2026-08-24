"""
validate_cascade.py

Validation suite for the R134a/CO2 cascade cycle in
`vapor_compression_cascade.py`. Run this after ANY change to that file.

Three stages, cheapest first:

  Stage 1 -- structural self-test. Same fluid (R134a) on BOTH loops. This
             isolates the cascade plumbing (energy balance, approach
             constraint, free mass flow) from any one fluid's physics.
             If this fails, the coupling is broken, not the fluid.

  Stage 2 -- invariant checks on the real R134a/CO2 pair at one point.
             Verifies the things that must hold in ANY valid solution:
               * cascade approach == cascade_approach_dT exactly
               * cascade energy balance closes to ~0
               * CO2 stays subcritical
               * compressor pressure ratios inside their cap
               * NOTHING pinned at a derived pressure bound -- if a
                 variable sits on a bound, the "solution" is being set by
                 our constants rather than by physics. This check is the
                 direct regression test for the bug that motivated this
                 file (see BREADCRUMB.md): hardcoded R134a-shaped bounds
                 (50-2000 / 100-5000 kPa, condenser outlet 30-50 C) made
                 the CO2 loop infeasible by construction, with constant
                 residuals of 18.7158 K and 0.6487 MPa that never moved
                 with operating conditions.

  Stage 3 -- full ambient sweep, at conditions matched to the existing
             single-stage benchmark so the two are directly comparable.

Usage:  python3 validate_cascade.py
"""

import logging
import numpy as np
from pyomo.environ import value

from vapor_compression_cascade import CascadeCycle, DEFAULT_FLUIDS, C_to_K

logging.disable(logging.CRITICAL)

# Matched to R1234yf/cop_vs_ambient_plr_benchmark_combined.csv so the cascade
# numbers can be compared against the single-stage R134a/R1234yf/R515B study.
BENCHMARK_PLR = 0.75
BENCHMARK_CD = 0.13
BENCHMARK_EVAP_SAT_C = -29.0
BENCHMARK_CONDENSER_APPROACH_C = 9.0
BENCHMARK_AMBIENTS_C = np.arange(10, 46, 5)

CASCADE_APPROACH_DT = 3.0

# ---------------------------------------------------------------------------
# Two DIFFERENT fluid configs, and it matters which one a comparison uses.
#
# DEFAULT_FLUIDS (in vapor_compression_cascade.py) is the cascade's own DESIGN:
# compressor efficiency 0.75, superheat/subcool capped at 3 K (R134a) and 5 K
# (CO2) -- the ranges chosen deliberately for this project.
#
# BENCHMARK_FLUIDS below instead mirrors the single-stage benchmark's settings
# exactly (run_plr_benchmark_r1234yf_vs_r134a.py): COMPRESSOR_EFFICIENCY = 0.80,
# SUPERHEAT_MAX_C = SUBCOOL_MAX_C = 8.0. Use THIS for any cascade-vs-single-stage
# comparison.
#
# Why this matters (see BREADCRUMB.md): comparing the as-designed cascade
# (eta 0.75, caps 3/5 K) against the single-stage benchmark (eta 0.80, cap 8 K)
# makes the cascade look ~7% worse at ambient 35 C. Matching the settings flips
# it to ~5% BETTER. The entire apparent deficit was the settings mismatch, not
# cascade physics.
BENCHMARK_FLUIDS = {
    "hot": {
        "gh_component": "r134a", "coolprop_name": "R134a", "custom_json": None,
        "efficiency": 0.80, "superheat_max": 8.0, "subcool_max": 8.0, "Tmax_C": None,
    },
    "cold": {
        "gh_component": "co2", "coolprop_name": "CO2", "custom_json": None,
        "efficiency": 0.80, "superheat_max": 8.0, "subcool_max": 8.0, "Tmax_C": None,
    },
}


def _build(fluids, ambient_C, cold_evap_C, plr=BENCHMARK_PLR, cd=BENCHMARK_CD,
           cascade_approach_dT=CASCADE_APPROACH_DT):
    cycle = CascadeCycle(fluids=fluids, PLR=plr, CD=cd)
    cycle.specify_initial_conditions(hot_ambient_C=ambient_C, cold_evap_C=cold_evap_C)
    cycle.initialize(verbose=False)
    cycle.set_specifications(
        hot={"ambient_temperature": ambient_C,
             "condenser_approach": BENCHMARK_CONDENSER_APPROACH_C},
        cold={"evap_sat_temperature": cold_evap_C},
        cascade_approach_dT=cascade_approach_dT,
    )
    return cycle


def _pinned_at_pressure_bound(cycle, rtol=1e-6):
    """Return a list of (name, value_kPa, which_bound) for any pressure
    variable sitting on its bound. Empty list = bounds are permissive."""
    fs = cycle.model.fs
    pinned = []
    for role in cycle.ROLES:
        for uname in ("evaporator", "compressor", "condenser", "expansion_valve"):
            unit = getattr(fs, f"{role}_{uname}")
            for pname, port in (("in", unit.inlet), ("out", unit.outlet)):
                P = port.pressure[0]
                v = value(P)
                if P.lb is not None and abs(v - P.lb) / max(abs(P.lb), 1.0) < rtol:
                    pinned.append((f"{role}_{uname}.{pname}", v / 1000.0, "lb"))
                elif P.ub is not None and abs(v - P.ub) / max(abs(P.ub), 1.0) < rtol:
                    pinned.append((f"{role}_{uname}.{pname}", v / 1000.0, "ub"))
    return pinned


def stage1_structural_self_test():
    print("=" * 72)
    print("STAGE 1 -- structural self-test (R134a on BOTH loops)")
    print("=" * 72)
    same_fluid = {
        role: {"gh_component": "r134a", "coolprop_name": "R134a", "custom_json": None,
               "efficiency": 0.75, "superheat_max": 3.0, "subcool_max": 3.0,
               "Tmax_C": None}
        for role in ("hot", "cold")
    }
    # Mild conditions: this stage tests plumbing, not envelope limits.
    cycle = _build(same_fluid, ambient_C=25.0, cold_evap_C=-10.0)
    cop, converged = cycle.optimize_COP(verbose=False)
    print(f"  converged = {converged}   COP = {cop if cop is None else round(cop, 4)}")
    if not converged:
        print("  FAIL: cascade coupling itself is broken (fluid-independent).")
        return False
    print("  PASS: two-loop coupling works.\n")
    return True


def stage2_invariants(ambient_C=35.0, cold_evap_C=BENCHMARK_EVAP_SAT_C):
    print("=" * 72)
    print(f"STAGE 2 -- invariants, real R134a/CO2 @ ambient={ambient_C}C, "
          f"evap={cold_evap_C}C")
    print("=" * 72)
    cycle = _build(DEFAULT_FLUIDS, ambient_C, cold_evap_C)

    print("  derived bounds (option c -- from each fluid's own sat curve):")
    for role, b in cycle._derived_bounds.items():
        lo = tuple(round(x, 1) for x in b["low_side_pressure_kPa"])
        hi = tuple(round(x, 1) for x in b["high_side_pressure_kPa"])
        print(f"    {role:4s} ({b['fluid']:5s}) evap faces {b['evaporator_faces']:7s} "
              f"cond faces {b['condenser_faces']:7s}")
        print(f"           low  P {lo} kPa   high P {hi} kPa")

    cop, converged = cycle.optimize_COP(verbose=False)
    print(f"\n  converged = {converged}   COP = {cop if cop is None else round(cop, 4)}")
    if not converged:
        print("  FAIL: did not converge.")
        return False

    fs = cycle.model.fs
    ok = True

    he = value(fs.hot_evaporator.control_volume.properties_out[0].temperature_sat) - C_to_K
    cc = value(fs.cold_condenser.control_volume.properties_out[0].temperature_sat) - C_to_K
    approach = cc - he
    good = abs(approach - CASCADE_APPROACH_DT) < 1e-4
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] cascade approach = {approach:.6f} K "
          f"(target {CASCADE_APPROACH_DT})")

    qh = value(fs.hot_evaporator.heat_duty[0])
    qc = value(fs.cold_condenser.heat_duty[0])
    good = abs(qh + qc) < 1e-3 * max(abs(qh), 1.0)
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] energy balance: {qh:.1f} + ({qc:.1f}) "
          f"= {qh + qc:.3e} W")

    good = cc < 30.978
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] CO2 condensing {cc:.2f} C < Tcrit 30.98 C "
          f"(subcritical)")

    for role in cycle.ROLES:
        rp = value(getattr(fs, f"{role}_compressor").ratioP[0])
        good = 1.1 - 1e-6 <= rp <= 8.0 + 1e-6
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] {role} pressure ratio = {rp:.3f} "
              f"(cap 8)")

    pinned = _pinned_at_pressure_bound(cycle)
    ok &= not pinned
    if pinned:
        print("  [FAIL] variables pinned at a derived pressure bound -- the "
              "solution is set by our constants, not physics:")
        for name, v, which in pinned:
            print(f"           {name} = {v:.1f} kPa at {which}")
    else:
        print("  [PASS] nothing pinned at a derived pressure bound")

    print()
    return ok


def stage3_ambient_sweep(cold_evap_C=BENCHMARK_EVAP_SAT_C, fluids=None,
                          label="as-designed"):
    fluids = fluids if fluids is not None else DEFAULT_FLUIDS
    print("=" * 72)
    print(f"STAGE 3 -- ambient sweep [{label}] @ evap={cold_evap_C}C, "
          f"PLR={BENCHMARK_PLR}, CD={BENCHMARK_CD}, "
          f"eta={fluids['hot']['efficiency']}, SH/SC cap={fluids['hot']['subcool_max']}K")
    print("=" * 72)
    T_cold_K = cold_evap_C + C_to_K
    rows = []
    for ambient_C in BENCHMARK_AMBIENTS_C:
        cop_carnot = T_cold_K / ((ambient_C + C_to_K) - T_cold_K)
        try:
            cycle = _build(fluids, float(ambient_C), cold_evap_C)
            cop, converged = cycle.optimize_COP(verbose=False)
        except Exception as exc:
            print(f"  ambient={ambient_C:3.0f} C  EXCEPTION {type(exc).__name__}")
            rows.append({"ambient_C": float(ambient_C), "converged": False})
            continue
        if not converged:
            print(f"  ambient={ambient_C:3.0f} C  FAILED")
            rows.append({"ambient_C": float(ambient_C), "converged": False})
            continue
        fs = cycle.model.fs
        he = value(fs.hot_evaporator.control_volume.properties_out[0].temperature_sat) - C_to_K
        rows.append({
            "ambient_C": float(ambient_C),
            "converged": True,
            "cop_full": cycle._last_cop_full,
            "cop_part": cycle._last_cop_part,
            "cop_carnot": cop_carnot,
            "second_law_efficiency": cycle._last_cop_full / cop_carnot,
            "intermediate_T_C": he,
            "hot_ratioP": value(fs.hot_compressor.ratioP[0]),
            "cold_ratioP": value(fs.cold_compressor.ratioP[0]),
            "hot_mass_flow": value(fs.hot_evaporator.inlet.flow_mass[0]),
        })

    print()
    print(f"  {'amb':>4} {'COP_full':>9} {'COP_part':>9} {'Carnot':>7} {'eta_II':>7} "
          f"{'T_int':>7} {'rP_hot':>7} {'rP_cold':>8} {'mdot_hot':>9}")
    for r in rows:
        if not r["converged"]:
            print(f"  {r['ambient_C']:4.0f}  FAILED")
            continue
        print(f"  {r['ambient_C']:4.0f} {r['cop_full']:9.4f} {r['cop_part']:9.4f} "
              f"{r['cop_carnot']:7.3f} {r['second_law_efficiency']:7.3f} "
              f"{r['intermediate_T_C']:7.2f} {r['hot_ratioP']:7.3f} "
              f"{r['cold_ratioP']:8.3f} {r['hot_mass_flow']:9.3f}")

    n_ok = sum(1 for r in rows if r["converged"])
    print(f"\n  {n_ok}/{len(rows)} converged\n")
    return rows


if __name__ == "__main__":
    s1 = stage1_structural_self_test()
    s2 = stage2_invariants()
    rows = stage3_ambient_sweep()
    n_ok = sum(1 for r in rows if r["converged"])
    print("=" * 72)
    print(f"SUMMARY: stage1={'PASS' if s1 else 'FAIL'}  "
          f"stage2={'PASS' if s2 else 'FAIL'}  "
          f"stage3={n_ok}/{len(rows)} converged")
    print("=" * 72)
