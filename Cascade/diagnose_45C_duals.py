"""
diagnose_45C_duals.py

Author:  Claude (Anthropic), for Shilpa Narasimhan

Solve the 45 C cascade once per fluid (same settings as diagnose_45C.py)
and print, for each superheat/subcool constraint in both loops:
    slack -- distance from the limit (~0 means the limit is active)
    dual  -- IPOPT's multiplier (how hard the limit is pushing back)
An active limit with a dual ~0 (touching, not pushing) would point to a
degenerate solution. Inactive limits should show duals ~1e-8.

Run:
    cd ~/idaes-hvacr-cycles/Cascade
    python diagnose_45C_duals.py R134a R1234yf
"""
import os
import sys
import time

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
os.chdir(THIS_DIR)

import diagnose_45C as d
from pyomo.environ import Suffix, value
from idaes.core.solvers import get_solver

CONSTRAINTS = [
    "evaporator.superheating_constraint",
    "evaporator.superheating_upper_constraint",
    "compressor.vapor_constraint",
    "condenser.subcooling_constraint",
    "condenser.subcooling_lower_constraint",
]


def run(fluid, ambient=45.0):
    # Build, initialize and specify exactly as diagnose_45C.py does.
    cycle, _ = d._build_cycle(fluid)
    cycle.specify_initial_conditions(hot_ambient_C=ambient, cold_evap_C=d.EVAP_SAT_TARGET_C)
    cycle.initialize(verbose=False)
    cold = dict(evap_sat_temperature=d.EVAP_SAT_TARGET_C, superheating=d.SUPERHEAT_C,
                superheating_max=d.SUPERHEAT_MAX_C, subcooling=d.SUBCOOL_C,
                subcooling_max=d.SUBCOOL_MAX_C, max_pressure_ratio=d.MAX_PRESSURE_RATIO)
    hot = dict(condenser_approach=d.COND_APPROACH_C, superheating=d.SUPERHEAT_C,
               superheating_max=d.SUPERHEAT_MAX_C, subcooling=d.SUBCOOL_C,
               subcooling_max=d.SUBCOOL_MAX_C, max_pressure_ratio=d.MAX_PRESSURE_RATIO,
               ambient_temperature=ambient)
    cycle.set_specifications(hot=hot, cold=cold, cascade_approach_dT=d.CASCADE_APPROACH_DT_C,
                             plr=d.PLR_VALUE, cd=d.CD_VALUE)

    m = cycle.model
    fs = m.fs
    m.dual = Suffix(direction=Suffix.IMPORT)   # ask Pyomo to bring back IPOPT's multipliers

    solver = get_solver()
    solver.options = {"max_iter": 1000, "tol": 1e-6, "linear_solver": "ma57"}

    # Warm-start pre-solve (objective off), as optimize_COP does. Advisory only.
    fs.compute_cop.deactivate()
    fs.obj.deactivate()
    try:
        solver.solve(m)
    except Exception as e:  # noqa: BLE001
        print("pre-solve:", e)
    work = value(fs.hot_compressor.work_mechanical[0]) + value(fs.cold_compressor.work_mechanical[0])
    fs.cop.set_value(value(fs.cold_evaporator.heat_duty[0]) / work)

    # The real (objective) solve -- one attempt, no retries.
    fs.compute_cop.activate()
    fs.obj.activate()
    t0 = time.time()
    res = solver.solve(m)
    tc = res.solver.termination_condition

    print(f"\n=== {fluid}/CO2 at {ambient} C: {tc}, COP = {value(fs.cop):.4f}, {time.time() - t0:.0f} s ===")
    print(f"{'constraint':<46}{'slack':>10}{'dual':>12}")
    for loop in ("hot", "cold"):
        for name in CONSTRAINTS:
            unit, cname = name.split(".")
            con = getattr(getattr(fs, f"{loop}_{unit}"), cname)
            body = value(con.body)
            slack = body - value(con.lower) if con.has_lb() else value(con.upper) - body
            print(f"{loop + '_' + name:<46}{slack:10.4f}{m.dual.get(con, float('nan')):12.3e}")
    T_split = value(fs.hot_evaporator.control_volume.properties_out[0].temperature_sat) - 273.15
    print(f"hot evaporator T_sat (split) = {T_split:.2f} C")


if __name__ == "__main__":
    for fluid in (sys.argv[1:] or ["R134a", "R1234yf"]):
        run(fluid)
