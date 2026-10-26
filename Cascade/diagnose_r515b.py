"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/09/2026

What this does
    Runs the R515B/CO2 cascade one ambient point at a time, stage by stage,
    to find where (and at which ambient) it fails.

Why
    The 09/15 sweep never converged for R515B. Checking every point here,
    with a stage label and physics checks, is faster and clearer than
    running the full sweep blind.

Retry
    After all points, each failed point is re-solved starting from the
    nearest converged point's solution (warm start). Lines say "(warm from X C)".

Stages (per point)
    seed    starting guesses at this ambient
    init    initialize every unit
    specs   set the specs (setpoints + cascade coupling). DOF should be 5.
    solve   warm-up (no objective), then maximize COP. ok should be True.

Physics checks (printed per point, after a solve)
    comp_dh  R515B compressor h_out - h_in, kJ/kg   (must be > 0)
    x_out    R515B compressor outlet vapor fraction (must be 1)
    co2_sh   CO2 evaporator superheat, K             (3 to 8)
    cold_end CO2 out - T_CAS,E, K                    (must be >= 1)

Switches (bottom of file)
    RUN_ALL_AMBIENTS  R515B at every ambient in the sweep (-30..50 C)
    RUN_ONE_AMBIENT   R515B at AMBIENT only, full IPOPT output
    RUN_R134A_PH      R134a/CO2 forced into P-h (the 10/10 isolation test)

Settings
    All setpoints come from run_cascade_sweep_r515b.py (imported as r),
    so this always matches the sweep -- including head pressure control
    (the model sees r.effective_ambient(ambient); lines show the real ambient).

Run
    python diagnose_r515b.py 2>&1 | grep "^DIAG\\|^SUMMARY"
"""

from pyomo.environ import value
from idaes.core.util.model_statistics import degrees_of_freedom

import run_cascade_sweep_r515b as r
from vapor_compression_cascade_r515b import CascadeCycle, Mode

AMBIENT = 30.0      # used by RUN_ONE_AMBIENT
K = 273.15


def setup(cycle, ambient, stage):
    """seed -> init -> specs for one ambient. `stage` is a 1-item list that
    always holds the stage being worked on, so a crash can be labeled."""
    stage[0] = "seed"
    cycle.specify_initial_conditions(hot_ambient_C= ambient, cold_evap_C= r.EVAP_SAT_TARGET_C,
                                     cond_approach_C= r.COND_APPROACH_C)
    stage[0] = "init"
    cycle.initialize(verbose = False)

    stage[0] = "specs"
    hot = dict(ambient_temperature = ambient, condenser_approach = r.COND_APPROACH_C,
               superheating = r.SUPERHEAT_C, superheating_max = r.SUPERHEAT_MAX_C,
               subcooling = r.SUBCOOL_C, subcooling_max = r.SUBCOOL_MAX_C,
               max_pressure_ratio = r.MAX_PRESSURE_RATIO)
    cold = dict(evap_sat_temperature = r.EVAP_SAT_TARGET_C,
                superheating = r.SUPERHEAT_C, superheating_max = r.SUPERHEAT_MAX_C,
                subcooling = r.SUBCOOL_COLD_C, subcooling_max = r.SUBCOOL_MAX_C,
                max_pressure_ratio = r.MAX_PRESSURE_RATIO)
    cycle.set_specifications(hot = hot, cold = cold,
                             cascade_approach_dT = r.CASCADE_APPROACH_DT_C,
                             cascade_cold_end_dT = r.CASCADE_COLD_END_DT_C,
                             plr = r.PLR_VALUE, cd = r.CD_VALUE)
    return degrees_of_freedom(cycle.model)


def checks(cycle):
    """The four physics checks listed in the module docstring."""
    fs = cycle.model.fs
    comp = fs.hot_compressor.control_volume
    ce = fs.cold_evaporator.control_volume.properties_out[0]
    cc = fs.cold_condenser.control_volume.properties_out[0]
    t_cas_e = value(fs.hot_evaporator.control_volume.properties_in[0].temperature)
    return dict(
        comp_dh = (value(comp.properties_out[0].enth_mass) - value(comp.properties_in[0].enth_mass)) / 1e3,
        x_out = value(comp.properties_out[0].vapor_frac),
        co2_sh = value(ce.temperature - ce.temperature_sat),
        cold_end = value(cc.temperature) - t_cas_e,
        T_cas_e_C = t_cas_e - K,
    )


def run_point(ambient, verbose = False):
    """All four stages at one ambient. Never raises: returns one result row."""
    stage = ["build"]
    row = dict(ambient = ambient, stage = "build", dof = None, cop = None, ok = False)
    try:
        cycle, _ = r._build_cycle("R515B")
        row["dof"] = setup(cycle, r.effective_ambient(ambient), stage)   # head pressure control
        stage[0] = "solve"
        cop, ok = cycle.optimize_COP(verbose = verbose, initialize = True, optimize = True)
        row.update(cop = cop, ok = ok, stage = "done")
        row.update(checks(cycle))
    except Exception as exc:
        row.update(stage = stage[0], error = f"{type(exc).__name__}: {str(exc)[:80]}")
    return row


def run_point_warm(ambient, neighbour):
    """Retry for a failed point: solve the converged neighbour ambient first,
    then move the condenser ambient to `ambient` and re-solve from there
    (no warm-up pre-solve). Never raises: returns one result row."""
    stage = ["build"]
    row = dict(ambient = ambient, stage = "build", dof = None, cop = None, ok = False, warm_from = neighbour)
    try:
        cycle, _ = r._build_cycle("R515B")
        row["dof"] = setup(cycle, r.effective_ambient(neighbour), stage)
        stage[0] = "solve neighbour"
        cop, ok = cycle.optimize_COP(verbose = False, initialize = True, optimize = True)
        if not ok:
            raise RuntimeError(f"neighbour {neighbour} C did not converge")
        stage[0] = "warm solve"
        cycle.model.fs.hot_condenser.ambient_T.set_value(r.effective_ambient(ambient) + K)
        cop, ok = cycle.optimize_COP(verbose = False, initialize = False, optimize = True)
        row.update(cop = cop, ok = ok, stage = "done")
        row.update(checks(cycle))
    except Exception as exc:
        row.update(stage = stage[0], error = f"{type(exc).__name__}: {str(exc)[:80]}")
    return row


def show(row):
    """One line per point, grep-able by DIAG."""
    if "error" in row:
        print(f"DIAG {row['ambient']:6.1f} C | FAILED at {row['stage']:5s} | {row['error']}", flush = True)
        return
    tag = f" (warm from {row['warm_from']:.0f} C)" if row.get("warm_from") is not None else ""
    print(f"DIAG {row['ambient']:6.1f} C{tag} | DOF {row['dof']} | ok {str(row['ok']):5s} | COP {row['cop']:8.4f} | "
          f"T_CAS,E {row['T_cas_e_C']:6.2f} C | comp_dh {row['comp_dh']:7.2f} | x_out {row['x_out']:.3f} | "
          f"co2_sh {row['co2_sh']:5.2f} | cold_end {row['cold_end']:5.2f}", flush = True)


# --- What to run (switch on what you need) --------------------------------
RUN_ALL_AMBIENTS = True    # R515B at every sweep ambient, one line each
RUN_ONE_AMBIENT = False    # R515B at AMBIENT, full IPOPT output
RUN_R134A_PH = False       # R134a/CO2 forced into P-h, no R515B

if RUN_ALL_AMBIENTS:
    rows = []
    for a in r.AMBIENT_SWEEP_C:
        row = run_point(float(a))
        show(row)
        rows.append(row)
    # Retry pass: failed points, warm-started from the nearest converged one
    good = [x["ambient"] for x in rows if x["ok"]]
    for i, x in enumerate(rows):
        if x["ok"] or not good:
            continue
        nb = min(good, key = lambda a: (abs(a - x["ambient"]), -a))   # nearest; warmer wins a tie
        rows[i] = run_point_warm(x["ambient"], nb)
        show(rows[i])
    n_ok = sum(1 for x in rows if x["ok"])
    print(f"SUMMARY {n_ok}/{len(rows)} points converged", flush = True)

if RUN_ONE_AMBIENT:
    show(run_point(AMBIENT, verbose = True))

if RUN_R134A_PH:
    base, _ = r._build_cycle("R134a")          # only to reuse its fluid settings
    cycle = CascadeCycle(fluids = base.fluids, mode = Mode.PH)
    stage = ["build"]
    print("DIAG R134a P-h DOF =", setup(cycle, AMBIENT, stage))
    print("DIAG R134a P-h", cycle.optimize_COP(verbose = True, initialize = True, optimize = True))
