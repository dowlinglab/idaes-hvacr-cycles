"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/11/2026

What this does
    COP vs ambient (-30..50 C) for the R134a/CO2 and R1234yf/CO2 cascades,
    main engine (vapor_compression_cascade.py). One case per run:
        21C      head pressure control, condensing never below 21 C
        nolimit  condensing = ambient + 9 C at every ambient

Setpoints
    Same as diagnose_45C.py (imported as d): CO2 evap -29 C, ambient + 9 C
    condenser, SH/SC 3..8 K (CO2 SC 0..8 K), 3 K cascade gap, 1 K cold end,
    PLR 0.75, CD 0.13, ratio cap 20.

Head pressure control
    The model sees effective_ambient = max(ambient, 21 - 9) = max(ambient, 12 C).
    Every ambient <= 12 C is the same cycle. CSV keeps the real ambient.

Cascade window
    cascade_T_window_C = (-26, 20) C instead of the engine default (-15, 20).
    The default's derived pressure floor pinned T_CAS,E at ~-26 C at cold
    no-limit ambients. physics_checks.py now flags any pressure on a bound.

Per point
    1. Fresh build, seeded and specified at the ambient the model sees.
    2. optimize_COP. It already tries, in order: default solve + retries,
       scaled fallback, pin-and-verify, feasibility acceptance.
    3. Physics checks (physics_checks.py). A point counts only if it is
       certified (not feasibility-accepted) AND passes every check.
Then a retry pass
    Each point that did not count is re-solved warm, from the nearest point
    that did (same idea as the R515B 15 C fix).

Outputs (this folder)
    sweep_cascade_<fluid>_<case>.csv   (COP, checks, heat balance, every stream)

Run
    python run_cascade_sweep_hpc.py 21C        (or nolimit)
    python run_cascade_sweep_hpc.py 21C R134a  (one fluid only)
"""

import csv
import sys
import time

import numpy as np
from pyomo.environ import value

import diagnose_45C as d
from physics_checks import cascade_checks
from heat_balance import cascade_balance

K = 273.15
AMBIENT_SWEEP_C = np.arange(-30, 55, 5)
FLUIDS = ["R134a", "R1234yf"]

CASE = sys.argv[1] if len(sys.argv) > 1 else "21C"
assert CASE in ("21C", "nolimit"), "case must be 21C or nolimit"
# 21 C = Copeland's lowest condenser temperature
# (https://copeland.custhelp.com/app/answers/detail/a_id/3278/)
T_COND_MIN_C = 21.0 if CASE == "21C" else None
if len(sys.argv) > 2:
    FLUIDS = sys.argv[2:]


# Cascade window for the engine's derived pressure bounds. Default (-15, 20) C
# puts the hot-loop evaporating floor at 0.7*P_sat(-18 C) = T_CAS,E ~ -26 C,
# which BINDS at cold no-limit ambients (10/11). -26 C = CO2 evap -29 + 3 K gap,
# the lowest the cascade can physically go, so the bound never binds.
CASCADE_T_WINDOW_C = (-26.0, 20.0)


def effective_ambient(ambient_c):
    """Ambient the model sees. Head pressure control holds condensing at
    or above T_COND_MIN_C; nolimit returns the real ambient."""
    if T_COND_MIN_C is None:
        return ambient_c
    return max(ambient_c, T_COND_MIN_C - d.COND_APPROACH_C)


def build(fluid, amb_model):
    """Fresh cascade, seeded + specified at amb_model. Not solved."""
    cycle, _ = d._build_cycle(fluid)
    cycle.specify_initial_conditions(hot_ambient_C=amb_model, cold_evap_C=d.EVAP_SAT_TARGET_C)
    cycle.initialize(verbose=False)
    cold = dict(evap_sat_temperature=d.EVAP_SAT_TARGET_C, superheating=d.SUPERHEAT_C,
                superheating_max=d.SUPERHEAT_MAX_C, subcooling=d.SUBCOOL_COLD_C,
                subcooling_max=d.SUBCOOL_MAX_C, max_pressure_ratio=d.MAX_PRESSURE_RATIO)
    hot = dict(ambient_temperature=amb_model, condenser_approach=d.COND_APPROACH_C,
               superheating=d.SUPERHEAT_C, superheating_max=d.SUPERHEAT_MAX_C,
               subcooling=d.SUBCOOL_C, subcooling_max=d.SUBCOOL_MAX_C,
               max_pressure_ratio=d.MAX_PRESSURE_RATIO)
    cycle.set_specifications(hot=hot, cold=cold, cascade_approach_dT=d.CASCADE_APPROACH_DT_C,
                             cascade_cold_end_dT=d.CASCADE_COLD_END_DT_C,
                             cascade_T_window_C=CASCADE_T_WINDOW_C,
                             plr=d.PLR_VALUE, cd=d.CD_VALUE)
    return cycle


def record(cycle, ambient, amb_model, cop, ok, how):
    """One result row: COP, which tier certified it, physics checks."""
    row = dict(ambient_C=ambient, ambient_model_C=amb_model, cop_full=cop, solver_ok=ok, how=how,
               scaled=getattr(cycle, "used_scaling_fallback", False),
               pin_verify=getattr(cycle, "used_pin_verify", False),
               feas_accepted=getattr(cycle, "accepted_on_feasibility", False))
    if ok:
        row["cop_part"] = cycle.get_part_load_cop()
        info, failed = cascade_checks(cycle, amb_model, cold_end_dT=d.CASCADE_COLD_END_DT_C)
        row.update(info)
        row["checks_failed"] = " ".join(failed)
        row.update(cascade_balance(cycle))          # heat balance + every stream (T, P, h, flow)
    else:
        row["checks_failed"] = "not converged"
    row["counts"] = bool(ok and not row["feas_accepted"] and not row["checks_failed"])
    return row


def show(fluid, row):
    tag = f" ({row['how']})" if row["how"] != "cold" else ""
    if row["counts"]:
        print(f"ROW {fluid:8s} {CASE:7s} {row['ambient_C']:6.1f} C{tag} | COP {row['cop_full']:.4f} | "
              f"T_CAS,E {row['T_cas_e_C']:6.2f} | scaled {row['scaled']} pin {row['pin_verify']} | checks ok",
              flush=True)
    else:
        print(f"ROW {fluid:8s} {CASE:7s} {row['ambient_C']:6.1f} C{tag} | DOES NOT COUNT | "
              f"ok {row['solver_ok']} feas {row['feas_accepted']} | {row['checks_failed']}", flush=True)


def run_fluid(fluid):
    rows, cycles = {}, {}
    for a in AMBIENT_SWEEP_C:
        a = float(a)
        amb_model = effective_ambient(a)
        t0 = time.time()
        try:
            cycle = build(fluid, amb_model)
            cop, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=True)
            row = record(cycle, a, amb_model, cop, ok, "cold")
            if row["counts"]:
                cycles[a] = cycle
        except Exception as exc:
            row = dict(ambient_C=a, ambient_model_C=amb_model, solver_ok=False, counts=False, how="cold",
                       checks_failed=f"ERROR {type(exc).__name__}: {str(exc)[:60]}",
                       feas_accepted=False)
        row["seconds"] = round(time.time() - t0)
        rows[a] = row
        show(fluid, row)

    # Retry pass: warm start from the nearest point that counts
    for a in sorted(rows):
        if rows[a]["counts"] or not cycles:
            continue
        nb = min(cycles, key=lambda b: (abs(b - a), -b))
        amb_model = effective_ambient(a)
        t0 = time.time()
        try:
            cycle = build(fluid, effective_ambient(nb))
            cop, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=True)
            cycle.model.fs.hot_condenser.ambient_T.set_value(amb_model + K)
            cop, ok = cycle.optimize_COP(verbose=False, initialize=False, optimize=True)
            row = record(cycle, a, amb_model, cop, ok, f"warm from {nb:.0f} C")
        except Exception as exc:
            row = dict(rows[a], how=f"warm from {nb:.0f} C failed: {type(exc).__name__}")
        row["seconds"] = round(time.time() - t0)
        if row["counts"]:
            rows[a] = row
        show(fluid, row)

    out = f"sweep_cascade_{fluid}_{CASE}.csv"
    keys = []
    for r in rows.values():
        keys += [k for k in r if k not in keys]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for a in sorted(rows):
            w.writerow(rows[a])
    n = sum(r["counts"] for r in rows.values())
    print(f"SUMMARY {fluid} {CASE}: {n}/{len(rows)} points count -> {out}", flush=True)


if __name__ == "__main__":
    for fl in FLUIDS:
        run_fluid(fl)
