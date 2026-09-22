"""
Author: Shilpa Narasimhan

Support: Claude AI

Date Created: 09/22/2026

PATH-TRACING diagnostic for the R1234yf @ 45C cascade problem.

Context: diagnose_45C.py showed the free (searching) solve at R1234yf @ 45C
never confirms optimal -- it hits the 1000-iteration cap while superheat/
subcool bounce around. diagnose_45C_pinned.py showed that if you REMOVE the
search (pin superheat=subcool=8K directly), the model solves that exact
point cleanly in ~17 iterations. So the destination (8K/8K) is fine; the
problem is somewhere in how the free search GETS there.

This script makes that "getting there" visible as a plot. IPOPT doesn't
expose a live per-iteration callback through the plain Pyomo solver
interface used here, so instead this reconstructs the path by RE-SOLVING
the exact same problem from scratch, over and over, each time capping
max_iter at a different value (1, 3, 5, ... up to 1000) and recording
where superheat/subcool/COP actually are at that cutoff. Because every
run starts from an identical fresh build + identical warm start (IPOPT is
deterministic given the same starting point and options), the max_iter=50
run's ending point IS where the real, uncapped run actually was after its
50th iteration too -- so stringing these together across increasing
max_iter values reconstructs the real trajectory.

Only the MAIN (COP-maximizing) solve's max_iter is varied -- the warm-start
presolve always gets a full, uncapped 1000 iterations first, exactly like
the real code, so every checkpoint starts from the same place.

Output: diag45_path_r1234yf.csv (one row per checkpoint) and
diag45_path_r1234yf.png (SH/SC vs. iteration, with COP as a second panel).
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from idaes.core.solvers import get_solver
from pyomo.environ import value

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT_DIR = os.path.dirname(_THIS_DIR)
if _REPO_ROOT_DIR not in sys.path:
    sys.path.insert(0, _REPO_ROOT_DIR)

from vapor_compression_cascade import CascadeCycle

# ---- Same setpoints as the real sweep / diagnose_45C.py (REAL free range,
# NOT pinned -- floor=3, cap=8, the actual problematic search window) ----
AMBIENT_C = 45.0
EVAP_SAT_TARGET_C = -29.0
COND_APPROACH_C = 9.0
SUPERHEAT_C = 3.0
SUBCOOL_C = 3.0
SUPERHEAT_MAX_C = 8.0
SUBCOOL_MAX_C = 8.0
COMPRESSOR_EFFICIENCY = 0.80
MAX_PRESSURE_RATIO = 20.0
CASCADE_APPROACH_DT_C = 3.0
PLR_VALUE = 0.75
CD_VALUE = 0.13

# Dense early (the path moves fast early on), coarser later, up through the
# real 1000-iteration cap.
CHECKPOINTS = sorted(set(
    list(range(1, 20, 2))
    + list(range(20, 100, 10))
    + list(range(100, 1001, 50))
))


def _build_r1234yf_cycle():
    cold_fluid = {
        "gh_component": "co2", "coolprop_name": "CO2", "custom_json": None,
        "efficiency": COMPRESSOR_EFFICIENCY,
        "superheat_max": SUPERHEAT_MAX_C, "subcool_max": SUBCOOL_MAX_C,
        "Tmax_C": None,
    }
    hot_fluid = {
        "gh_component": "r1234yf", "coolprop_name": "R1234yf",
        "custom_json": os.path.join(_REPO_ROOT_DIR, "R1234yf", "r1234yf_general_helmholtz.json"),
        "efficiency": COMPRESSOR_EFFICIENCY,
        "superheat_max": SUPERHEAT_MAX_C, "subcool_max": SUBCOOL_MAX_C,
        "Tmax_C": None,
    }
    return CascadeCycle(fluids={"hot": hot_fluid, "cold": cold_fluid})


def _fresh_cycle_at_start():
    """Build + initialize + specify a brand-new cycle, identical every call
    -- this is the common starting point every checkpoint solves from."""
    cycle = _build_r1234yf_cycle()
    cycle.specify_initial_conditions(hot_ambient_C=AMBIENT_C, cold_evap_C=EVAP_SAT_TARGET_C)
    cycle.initialize(verbose=False)

    cold_kwargs = dict(
        evap_sat_temperature=EVAP_SAT_TARGET_C, superheating=SUPERHEAT_C,
        superheating_max=SUPERHEAT_MAX_C, subcooling=SUBCOOL_C,
        subcooling_max=SUBCOOL_MAX_C, max_pressure_ratio=MAX_PRESSURE_RATIO,
    )
    hot_kwargs = dict(
        condenser_approach=COND_APPROACH_C, superheating=SUPERHEAT_C,
        superheating_max=SUPERHEAT_MAX_C, subcooling=SUBCOOL_C,
        subcooling_max=SUBCOOL_MAX_C, max_pressure_ratio=MAX_PRESSURE_RATIO,
        ambient_temperature=AMBIENT_C,
    )
    cycle.set_specifications(hot=hot_kwargs, cold=cold_kwargs,
                              cascade_approach_dT=CASCADE_APPROACH_DT_C,
                              plr=PLR_VALUE, cd=CD_VALUE)
    return cycle


def _sh_sc_cop(fs):
    """Cold-loop evaporator superheat, hot-loop condenser subcool, COP --
    same three fields the console reports use, read straight off the model."""
    sh = (value(fs.cold_evaporator.control_volume.properties_out[0].temperature)
          - value(fs.cold_evaporator.control_volume.properties_out[0].temperature_sat))
    sc = (value(fs.hot_condenser.control_volume.properties_out[0].temperature_sat)
          - value(fs.hot_condenser.control_volume.properties_out[0].temperature))
    try:
        cop = value(fs.cop)
    except Exception:
        cop = float("nan")
    return sh, sc, cop


def _solve_capped(cycle, max_iter_cap):
    """Warm-start (uncapped, 1000 iters, identical every call) then the
    real COP-maximizing solve, capped at max_iter_cap iterations. Returns
    the termination condition, ok flag, and SH/SC/COP at that cutoff."""
    solver = get_solver()
    fs = cycle.model.fs

    solver.options = {"max_iter": 1000, "tol": 1e-6, "linear_solver": "ma57"}
    fs.compute_cop.deactivate()
    fs.obj.deactivate()
    CascadeCycle._safe_solve(solver, cycle.model, tee=False)
    try:
        work_total = (value(fs.hot_compressor.work_mechanical[0])
                      + value(fs.cold_compressor.work_mechanical[0]))
        if work_total:
            fs.cop.set_value(value(fs.cold_evaporator.heat_duty[0]) / work_total)
    except Exception:  # noqa: BLE001 -- warm start only
        pass

    fs.compute_cop.activate()
    fs.obj.activate()
    solver.options = {"max_iter": max_iter_cap, "tol": 1e-6, "linear_solver": "ma57"}
    tc, ok = CascadeCycle._safe_solve(solver, cycle.model, tee=False)

    sh, sc, cop = _sh_sc_cop(fs)
    return tc, ok, cop, sh, sc


def main():
    rows = []
    for k in CHECKPOINTS:
        cycle = _fresh_cycle_at_start()
        tc, ok, cop, sh, sc = _solve_capped(cycle, k)
        rows.append((k, tc, ok, cop, sh, sc))
        print(f"max_iter={k:>4d} | tc={tc:<28s} | ok={int(ok)} | "
              f"COP={cop:8.4f} | SH={sh:6.3f}K | SC={sc:6.3f}K")

    with open("diag45_path_r1234yf.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["max_iter_cap", "termination_condition", "ok", "cop_full", "sh_actual_K", "sc_actual_K"])
        for k, tc, ok, cop, sh, sc in rows:
            w.writerow([k, tc, int(ok), cop, sh, sc])
    print("Saved: diag45_path_r1234yf.csv")

    iters = [r[0] for r in rows]
    sh_vals = [r[4] for r in rows]
    sc_vals = [r[5] for r in rows]
    cop_vals = [r[3] for r in rows]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 7.5), dpi=160, sharex=True)

    ax1.plot(iters, sh_vals, marker="o", markersize=3, linewidth=1.5, label="SH_actual (cold_evap)")
    ax1.plot(iters, sc_vals, marker="s", markersize=3, linewidth=1.5, label="SC_actual (hot_cond)")
    ax1.axhline(8.0, linestyle="--", color="black", linewidth=1.2, alpha=0.7, label="8K (confirmed correct corner)")
    ax1.axhline(3.0, linestyle=":", color="gray", linewidth=1.0, alpha=0.6, label="3K (floor)")
    ax1.set_ylabel("Superheat / Subcool (K)")
    ax1.set_title("R1234yf @ 45C: free-search optimizer path (reconstructed by capped reruns)")
    ax1.legend(loc="best", fontsize=8)
    ax1.grid(True, alpha=0.25)

    ax2.plot(iters, cop_vals, marker="d", markersize=3, linewidth=1.5, color="darkred", label="COP_full")
    ax2.set_xlabel("Main-solve iteration cap (max_iter)")
    ax2.set_ylabel("COP_full")
    ax2.legend(loc="best", fontsize=8)
    ax2.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig("diag45_path_r1234yf.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved: diag45_path_r1234yf.png")


if __name__ == "__main__":
    main()
