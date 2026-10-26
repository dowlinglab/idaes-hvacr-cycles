"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/11/2026

What this does
    Second chance for R515B/CO2 cascade points that failed in
    run_cascade_sweep_r515b.py (even after its warm-start retry).
    Fixes the failed rows in sweep_cascade_R515B_<case>.csv in place.

Why
    No-limit, cold ambients: the engine derives pressure bounds from
    R515B's saturation curve, and that derivation fails below about -25 C
    ("pressure-bound derivation did not converge"). Warm starts can also
    be blocked by bounds derived at the neighbour's ambient.

Tiers (stops at the first that works)
    A  bounds    fresh build, explicit wide pressure bounds (no derivation)
    B  cont      build + solve at the nearest good ambient (wide bounds),
                 then walk the condenser ambient to the failed point in
                 1 K steps, re-solving warm each step (continuation)
    Pin-and-verify (main engine's _pin_and_verify, borrowed as-is) runs
    whenever a solve ends not-ok: pin T_CAS,E, Newton on the pin's
    multiplier (= dCOP/dT_CAS,E) until |g| < 1e-5 COP/K.
Every rescued point goes through physics_checks.py like the rest.

Wide bounds (kPa)   low side 30..600, high side 50..3000. These are safety
                    bounds only; the saturation constraints set the pressures.

Run
    python retry_r515b_points.py nolimit     (or 21C)
    python retry_r515b_points.py nolimit -30 -25 -20 -15
        also redo these ambients (10/11: their T_CAS,E sat on the engine's
        derived pressure floor, ~-26 C, so their COP was bound-limited)
"""

import csv
import sys

CASE = sys.argv[1] if len(sys.argv) > 1 else "21C"
FORCE = [float(a) for a in sys.argv[2:]]   # ambients to redo even if they "passed"
sys.argv = [sys.argv[0], CASE]          # the runner reads the case from argv
import run_cascade_sweep_r515b as r
from physics_checks import cascade_checks
from pyomo.environ import value
from vapor_compression_cascade import CascadeCycle as MainCascade

K = 273.15
WIDE_LOW_KPA = (30.0, 600.0)
WIDE_HIGH_KPA = (50.0, 3000.0)
CSV_PATH = f"sweep_cascade_R515B_{CASE}.csv"


def build(amb_model):
    """Fresh R515B cascade at amb_model, wide explicit pressure bounds."""
    cycle, _ = r._build_cycle("R515B")
    cycle.specify_initial_conditions(hot_ambient_C=amb_model, cold_evap_C=r.EVAP_SAT_TARGET_C,
                                     cond_approach_C=r.COND_APPROACH_C)
    cycle.initialize(verbose=False)
    hot = dict(ambient_temperature=amb_model, condenser_approach=r.COND_APPROACH_C,
               superheating=r.SUPERHEAT_C, superheating_max=r.SUPERHEAT_MAX_C,
               subcooling=r.SUBCOOL_C, subcooling_max=r.SUBCOOL_MAX_C,
               max_pressure_ratio=r.MAX_PRESSURE_RATIO,
               low_side_pressure=WIDE_LOW_KPA, high_side_pressure=WIDE_HIGH_KPA)
    cold = dict(evap_sat_temperature=r.EVAP_SAT_TARGET_C, superheating=r.SUPERHEAT_C,
                superheating_max=r.SUPERHEAT_MAX_C, subcooling=r.SUBCOOL_COLD_C,
                subcooling_max=r.SUBCOOL_MAX_C, max_pressure_ratio=r.MAX_PRESSURE_RATIO)
    cycle.set_specifications(hot=hot, cold=cold, cascade_approach_dT=r.CASCADE_APPROACH_DT_C,
                             cascade_cold_end_dT=r.CASCADE_COLD_END_DT_C,
                             plr=r.PLR_VALUE, cd=r.CD_VALUE)
    return cycle


def solve(cycle, initialize):
    """optimize_COP, then pin-and-verify if it did not certify."""
    cop, ok = cycle.optimize_COP(verbose=False, initialize=initialize, optimize=True)
    if not ok:
        ok, T, g = MainCascade._pin_and_verify(cycle)
        if ok:
            cop = value(cycle.model.fs.cop)
            print(f"   pin-and-verify certified: T_CAS,E {T - K:.2f} C, dCOP/dT {g:.1e}", flush=True)
    return cop, ok


def tier_a(amb):
    cycle = build(r.effective_ambient(amb))
    cop, ok = solve(cycle, True)
    return cycle, cop, ok


def tier_b(amb, nb):
    a0, a1 = r.effective_ambient(nb), r.effective_ambient(amb)
    cycle = build(a0)
    cop, ok = solve(cycle, True)
    if not ok:
        return cycle, cop, False
    step = 1.0 if a1 > a0 else -1.0
    t = a0
    while abs(t - a1) > 1e-9:
        t = t + step if abs(a1 - t) > 1.0 else a1
        cycle.model.fs.hot_condenser.ambient_T.set_value(t + K)
        cop, ok = solve(cycle, False)
        print(f"   cont {nb:.0f} -> {amb:.0f} C: at {t:6.1f} C ok {ok} COP {cop:.4f}", flush=True)
        if not ok:
            return cycle, cop, False
    return cycle, cop, True


def main():
    with open(CSV_PATH) as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        rows = list(reader)
    def ok_row(x):
        return (x["converged"] in ("1", 1) and not x["checks_failed"]
                and float(x["ambient_C"]) not in FORCE)
    good = [float(x["ambient_C"]) for x in rows if ok_row(x)]
    todo = [x for x in rows if not ok_row(x)]
    while todo and good:
        # failed point closest to a good one first, so the next can start from it
        x = min(todo, key=lambda y: min(abs(float(y["ambient_C"]) - g) for g in good))
        todo.remove(x)
        amb = float(x["ambient_C"])
        result = None
        attempts = [("bounds", lambda: tier_a(amb))]
        for nb in sorted(good, key=lambda a: (abs(a - amb), -a))[:2]:
            attempts.append((f"cont from {nb:.0f} C", lambda nb=nb: tier_b(amb, nb)))
        for how, fn in attempts:
            try:
                cycle, cop, ok = fn()
            except Exception as exc:
                print(f"RETRY {amb:6.1f} C {how}: ERROR {type(exc).__name__}: {str(exc)[:70]}", flush=True)
                continue
            if not ok:
                print(f"RETRY {amb:6.1f} C {how}: not converged", flush=True)
                continue
            info, failed = cascade_checks(cycle, r.effective_ambient(amb),
                                          sc_cold=(r.SUBCOOL_COLD_C, r.SUBCOOL_MAX_C),
                                          cold_end_dT=r.CASCADE_COLD_END_DT_C)
            print(f"RETRY {amb:6.1f} C {how}: COP {cop:.4f} checks {failed or 'ok'}", flush=True)
            if not failed:
                result = (cycle, cop, info, how)
                break
        if result is None:
            continue
        cycle, cop, info, how = result
        ccarnot = r._carnot_cop(amb)
        x.update(cop_full=cop, cop_part=cycle.get_part_load_cop(), cop_carnot=ccarnot,
                 second_law_efficiency=cop / ccarnot, converged=1, warm_start_from_C=how,
                 checks_failed="", T_cas_e_C=info["T_cas_e_C"], hot_P_evap_kPa=info["hot_P_evap_kPa"],
                 vacuum=info["vacuum"], energy_rel=info["energy_rel"])
        x.update(r._report_unit_ops(cycle))
        good.append(amb)
    with open(CSV_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    n = sum(1 for x in rows if x["converged"] in ("1", 1) and not x["checks_failed"])
    print(f"SUMMARY R515B {CASE}: {n}/{len(rows)} points count -> {CSV_PATH}", flush=True)


if __name__ == "__main__":
    main()
