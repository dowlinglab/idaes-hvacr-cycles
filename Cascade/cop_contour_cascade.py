"""
cop_contour_cascade.py

Author:  Shilpa Narasimhan
Support: Claude (Anthropic)

Why this file exists
--------------------
IPOPT never certifies the free R-1234yf/CO2 cascade at 45 C. The cause is
the cascade-HX level T_CAS,E (R-1234yf evaporating temperature): COP is
nearly flat along it, so IPOPT zigzags. R-134a/CO2 zigzags too and only
escapes by luck. Pinning T_CAS,E makes every solve certify in ~11 IPOPT
iterations, so we sweep it ourselves and pick the best -- the same method
as Yilmaz et al. (2014, Purdue IRACC paper 1541).

What is in here
---------------
build            set up one cascade (same specs as diagnose_45C.py)
cop_at           COP with T_CAS,E pinned at one value
sweep_45C_both   COP vs T_CAS,E at 45 C for R-1234yf and R-134a + plot
free_solve_log   free (unpinned) 45 C solve, IPOPT log saved to a file
read_ipopt_log   pull iteration, COP, inf_pr, inf_du out of a log
plot_ipopt_logs  plot those for both fluids (the "zigzag" figure)
best_cop         certified best COP at one ambient: coarse + fine sweep

Run
---
Pick what to run with the switches at the bottom (RUN_...), then:
    cd ~/idaes-hvacr-cycles/Cascade
    python cop_contour_cascade.py 2>&1 | grep -v WARNING
"""

import os
import sys
from pyomo.environ import value
import re
import numpy as np
import matplotlib.pyplot as plt
from idaes.core.solvers import get_solver
import CoolProp.CoolProp as CP
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
import diagnose_45C as d

CO2_T_CRIT_C = CP.PropsSI("Tcrit", "CO2") - 273.15   # ~31 C; caps T_CAS,E in best_cop

def build(fluid, ambient = 45.0):
    """
    Set up one cascade, ready to solve. Nothing is solved here.

    Same four calls and settings as diagnose_45C.py, so results are
    directly comparable to the real sweeps:
      1. _build_cycle               hot loop = fluid, cold loop = CO2
      2. specify_initial_conditions CoolProp starting guesses
      3. initialize                 each unit solved alone, one by one
      4. set_specifications         -29 C CO2 evaporator, ambient + 9 C
                                    condenser, SH/SC limits, 3 K cascade gap
    Careful: _build_cycle treats ANY name it doesn't know as R-134a,
    silently. Spell fluid exactly "R1234yf" or "R134a".
    """
    """
    Set up one cascade, ready to solve. Nothing is solved here.

    Same four calls and same settings as diagnose_45C.py, so results here
    are directly comparable to the real sweeps:
      1. _build_cycle            -> hot loop = fluid, cold loop = CO2
      2. specify_initial_conditions -> CoolProp starting guesses
      3. initialize              -> each unit solved on its own, one by one
      4. set_specifications      -> -29 C CO2 evaporator, ambient + 9 C
                                    condenser, SH/SC limits, 3 K cascade gap
    Careful: _build_cycle treats ANY name it doesn't know as R-134a,
    silently. Spell fluid exactly "R1234yf" or "R134a".
    """
    cycle, _ = d._build_cycle(fluid)
    cycle.specify_initial_conditions(hot_ambient_C= ambient, cold_evap_C= d.EVAP_SAT_TARGET_C)
    cycle.initialize(verbose = False)
    cold = dict(evap_sat_temperature = d.EVAP_SAT_TARGET_C,
                superheating = d.SUPERHEAT_C, superheating_max = d.SUPERHEAT_MAX_C,
                subcooling = d.SUBCOOL_COLD_C, subcooling_max = d.SUBCOOL_MAX_C,
                max_pressure_ratio = d.MAX_PRESSURE_RATIO)
    hot = dict(ambient_temperature = ambient, condenser_approach = d.COND_APPROACH_C,
               superheating = d.SUPERHEAT_C, superheating_max = d.SUPERHEAT_MAX_C,
               subcooling = d.SUBCOOL_C, subcooling_max = d.SUBCOOL_MAX_C,
               max_pressure_ratio = d.MAX_PRESSURE_RATIO)

    cycle.set_specifications(hot = hot, cold = cold,
                             cascade_approach_dT= d.CASCADE_APPROACH_DT_C,
                             cascade_cold_end_dT = d.CASCADE_COLD_END_DT_C,
                             plr = d.PLR_VALUE, cd = d.CD_VALUE)

    return cycle

from idaes.core.util.model_statistics import degrees_of_freedom
# cycle = build("R1234yf", ambient = 45.0)
# fs = cycle.model.fs
# print(type(fs.hot_evaporator))
# evap_out = fs.hot_evaporator.control_volume.properties_out[0]
# print("hot evap outlet T(K):", value(evap_out.temperature))
# print("hot evap T_sat(K):", value(evap_out.temperature_sat))
# print("degrees of freedom:", degrees_of_freedom(cycle.model))

# evap = fs.hot_evaporator
# evap.T_sat_set.set_value(5.0 + 273.15)
# evap.evap_sat_constraint.activate()
# print("degrees of freedom after pinning:", degrees_of_freedom(cycle.model))

# cop, ok = cycle.optimize_COP(verbose = True, initialize= True, optimize = True)
# print("T_CAS,E pinned at 5 C -> COP:", cop, "ok:", ok)

def cop_at(fluid, ambient, T_cas_e_C):
    """
    COP with the cascade-HX level pinned at T_CAS,E = T_cas_e_C.

    Pinning = switching on hot_evaporator.evap_sat_constraint
    (T_sat == T_sat_set). The model builds it but leaves it off for the
    hot loop. On -> degrees of freedom 5 -> 4; IPOPT only picks SH/SC
    and certifies in ~11 iterations.
    Returns (COP, ok). ok = True also covers optimize_COP's
    feasibility-accepted case; pinned solves normally certify outright.
    """
    """
    COP with the cascade-HX level pinned: T_CAS,E = T_cas_e_C.

    Pinning = switching on hot_evaporator.evap_sat_constraint
    (T_sat == T_sat_set). The model builds it but leaves it off for the
    hot loop. On -> one fewer degree of freedom (5 -> 4); IPOPT then only
    picks SH/SC, and certifies in ~11 iterations.
    Returns (COP, ok). ok = True also covers optimize_COP's
    feasibility-accepted case; pinned solves normally certify outright.
    """
    cycle = build(fluid, ambient = ambient)
    evap = cycle.model.fs.hot_evaporator
    evap.T_sat_set.set_value(T_cas_e_C + 273.15)
    evap.evap_sat_constraint.activate()
    cop, ok = cycle.optimize_COP(verbose = False, initialize = True, optimize = True)
    return cop, ok



def sweep_45C_both(fluids=("R1234yf", "R134a"), T_list=range(-10, 16)):
    """
    COP vs pinned T_CAS,E at 45 C for each fluid, then one plot.

    Plot shows each fluid's COP as % below its own peak, so the two
    curves can be compared by shape even though their COP levels differ.
    Saves: cop_vs_TcasE_45C_both.png
    """
    results = {}
    for fluid in fluids:
        results[fluid] = []
        for T in T_list:
            cop, ok = cop_at(fluid, 45.0, T)
            results[fluid].append((T, cop, ok))
            print(f"{fluid}  T_CAS,E = {T:4d} C   COP = {cop:.5f}   ok = {ok}")

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for fluid in fluids:
        T_vals = [r[0] for r in results[fluid]]
        cop_vals = [r[1] for r in results[fluid]]
        peak = max(cop_vals)
        rel = [100 * (c / peak - 1) for c in cop_vals]      # % below own peak
        ax.plot(T_vals, rel, "o-", label=f"{fluid}/CO$_2$")
    ax.set_xlabel("T_CAS,E (deg C)")
    ax.set_ylabel("COP relative to own peak (%)")
    ax.set_title("45 deg C ambient")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig("cop_vs_TcasE_45C_both.png", dpi=200)
    print("saved cop_vs_TcasE_45C_both.png")
    return results



def free_solve_log(fluid,logfile, ambient = 45.0):
    """
    One free solve (T_CAS,E NOT pinned) with IPOPT's iteration table saved.

    Same two stages as optimize_COP -- warm-up with no objective, then
    maximize COP -- but no retries, so the log shows one clean run.
    "output_file" is overwritten by each solve, so the file ends up holding
    the real (objective) solve, not the warm-up.
    """
    cycle = build(fluid, ambient = ambient)
    fs = cycle.model.fs
    solver = get_solver()
    solver.options = {"max_iter": 1000, "tol": 1e-6, "linear_solver": "ma57",
                      "output_file": logfile, "print_level": 5 }

    fs.compute_cop.deactivate()
    fs.obj.deactivate()
    try:
        solver.solve(cycle.model)
    except Exception as e:
        print("warm-up:", e)
    w = value(fs.hot_compressor.work_mechanical[0]) + value(fs.cold_compressor.work_mechanical[0])
    fs.cop.set_value(value(fs.cold_evaporator.heat_duty[0]) / w)

    fs.compute_cop.activate()
    fs.obj.activate()
    res = solver.solve(cycle.model)
    from idaes.core.util.model_diagnostics import DiagnosticsToolbox
    dt = DiagnosticsToolbox(cycle.model)
    dt.report_numerical_issues()
    print(fluid,res.solver.termination_condition, value(fs.cop))


def read_ipopt_log(logfile):
    """
    Pull the iteration table out of an IPOPT output_file.

    Returns four lists: iteration, COP, inf_pr, inf_du.
    - Iteration lines start with the iteration number (restoration-phase
      lines end it with "r", e.g. "12r"); every other line is skipped.
    - IPOPT minimizes, so it logs the objective as -COP; the sign is
      flipped back here.
    """
    pattern = re.compile(r"^\s*(\d+)r?\s+(\S+)\s+(\S+)\s+(\S+)\s")
    it, cop, inf_pr, inf_du = [], [], [], []
    with open(logfile) as f:
        for line in f:
            m = pattern.match(line)
            if m is None:
                continue
            try:
                obj, pr, du = (float(m.group(k)) for k in (2, 3, 4))
            except ValueError:          # header or text line that happened to match
                continue
            it.append(int(m.group(1)))
            cop.append(-obj)
            inf_pr.append(pr)
            inf_du.append(du)
    return it, cop, inf_pr, inf_du


def plot_ipopt_logs(fluids=("R134a", "R1234yf"), outfile="ipopt_paths_45C.png"):
    """
    Two panels from the free-solve logs: COP (top) and inf_du (bottom) per
    IPOPT iteration.

    What to look for: both fluids zigzag in inf_du (0.01-1) instead of
    dropping to tol = 1e-6. R-134a escapes with one big step near
    iteration ~347; R-1234yf never does.
    Needs: ipopt_45C_<fluid>.log from free_solve_log.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    for fluid in fluids:
        it, cop, inf_pr, inf_du = read_ipopt_log(f"ipopt_45C_{fluid}.log")
        ax1.plot(it, cop, lw=1.2, label=f"{fluid}/CO$_2$")
        ax2.semilogy(it, inf_du, lw=0.6)
    ax1.set_ylabel("COP")
    ax1.set_ylim(1.67, 1.745)           # zoom on the crawl; early iterations are far below
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.4)
    ax2.axhline(1e-6, color="k", ls="--", lw=1)
    ax2.text(ax2.get_xlim()[1], 1.5e-6, "tol = 1e-6 ", ha="right")
    ax2.set_ylabel("dual infeasibility (inf_du)")
    ax2.set_xlabel("IPOPT iteration")
    ax2.grid(True, alpha=0.4)
    ax1.set_title("Free 45 deg C cascade solves: IPOPT progress")
    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    print("saved", outfile)


def best_cop(fluid, ambient, coarse_step=2.0, fine_step=0.25):
    """
    Certified best COP at one ambient: sweep T_CAS,E, keep the best point.

    Why: letting IPOPT choose T_CAS,E zigzags (see plot_ipopt_logs).
    Every pinned solve certifies, so we do the 1-D search ourselves.

    Steps
    1. Range: from 4 K above the CO2 evaporator (-29 C) up to the lower of
       (a) 8 K below the hot condenser (ambient + 9 C), or
       (b) 5 K below where CO2 would condense at its critical point
           (31 C - 3 K gap - 5 K = ~23 C). Above that the CO2 side cannot
           condense, and solves fail (ValueError).
       The peak moves with ambient (about -14 C at 10 C, about 2 C at
       45 C), so a fixed range would miss it.
    2. Coarse sweep every coarse_step K.
    3. Fine sweep, every fine_step K, within +/- coarse_step of the
       coarse best.
    4. Points with ok = False are dropped, never trusted.

    Returns (T_CAS,E in C, COP), or None if nothing solved.
    Cost: ~35-45 pinned solves, ~1 s each; failed points are slower
    (optimize_COP retries them).
    """
    T_lo = d.EVAP_SAT_TARGET_C + 4.0
    T_hi = min(ambient + d.COND_APPROACH_C - 8.0,
               CO2_T_CRIT_C - d.CASCADE_APPROACH_DT_C - 5.0)

    def sweep(T_values):
        points = []
        for T in T_values:
            cop, ok = cop_at(fluid, ambient, float(T))
            print(f"  {fluid} {ambient:5.1f} C ambient  T_CAS,E = {T:6.2f} C  "
                  f"COP = {cop:.5f}  ok = {ok}")
            if ok:
                points.append((float(T), cop))
        return points

    coarse = sweep(np.arange(T_lo, T_hi + 1e-9, coarse_step))
    if not coarse:
        return None
    T_best, _ = max(coarse, key=lambda p: p[1])
    if T_best in (coarse[0][0], coarse[-1][0]):
        print(f"  WARNING: best point is at the edge of the range ({T_best} C) -- "
              "the true peak may lie outside it")

    fine = sweep(np.arange(T_best - coarse_step, T_best + coarse_step + 1e-9, fine_step))
    return max(coarse + fine, key=lambda p: p[1])


def audit(fluid, ambient=45.0):
    """Scaling audit: solve with all 5 DOF free, then ask IDAES what is
    numerically wrong with the model at that answer.

    Why: the free R1234yf solve at 45 C never certifies, but every pinned
    solve does. If the model is badly scaled or has near-duplicate
    constraints, IPOPT's step math gets shaky near the flat COP peak.
    This checks for that.

    Reads the printout:
        Jacobian condition number   big (>1e10) = ill-conditioned
        WARNINGS                    serious, look at these first
        Cautions                    worth a look, often harmless

    Returns the toolbox, so you can call its display_* methods after.
    """
    cycle = build(fluid, ambient=ambient)
    cop, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=True)
    print(f"AUDIT {fluid} at {ambient} C: COP = {cop:.5f}, ok = {ok}")
    dt = DiagnosticsToolbox(cycle.model)
    dt.report_numerical_issues()
    return dt


# --- What to run (switch on what you need) ------------------------------
RUN_SWEEP_45C = False      # COP vs T_CAS,E at 45 C, both fluids (~2 min)
RUN_FREE_LOGS = False      # free 45 C solves with IPOPT logs (~10 s)
RUN_PLOT_LOGS = False      # plot the logs (needs RUN_FREE_LOGS once)
RUN_BEST_COP = False       # certified best COP at 45 C, both fluids (~1.5 min)
RUN_AUDIT = True           # scaling audit of the free 45 C R1234yf solve (~15 s)

if RUN_SWEEP_45C:
    sweep_45C_both()

if RUN_FREE_LOGS:
    for fluid in ["R134a", "R1234yf"]:
        free_solve_log(fluid, f"ipopt_45C_{fluid}.log")

if RUN_PLOT_LOGS:
    plot_ipopt_logs()

if RUN_BEST_COP:
    for fluid in ["R1234yf", "R134a"]:
        result = best_cop(fluid, 45.0)
        print(f"BEST {fluid} at 45 C: T_CAS,E = {result[0]:.2f} C, COP = {result[1]:.5f}")

if RUN_AUDIT:
    dt = audit("R1234yf")
    dt.display_near_parallel_constraints()       # which equations are duplicates (LICQ check)
    dt.display_variables_at_or_outside_bounds()  # which variables sit on a bound
    print("\n=== EXTREME JACOBIAN COLUMNS ===")
    dt.display_variables_with_extreme_jacobians()
    print("\n=== EXTREME JACOBIAN ROWS ===")
    dt.display_constraints_with_extreme_jacobians()
    print("\n=== EXTREME JACOBIAN ENTRIES ===")
    dt.display_extreme_jacobian_entries()
