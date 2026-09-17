"""
diagnose_one_sensitivity_point_0917.py -- run ONE Tc/Pc sensitivity case in
isolation and SHOW what the cycle solve is actually doing.

Why this file exists (2026-09-17): the Tc/Pc sensitivity sweep
(phase6_final_sensitivity_Tc_*.py) reports 251 of 900 case-ambient solves as
non-converged, and the failures are not randomly placed -- every case at
Tc >= +3% fails at every ambient (180 of the 251), while small deviations
(|dev| <= 1%) come back scattered between 75% and 98% with no physical
pattern. The sweep only records a pass/fail boolean, so there is no way to
see WHY from its output.

This script runs a single case end to end, stage by stage, and reports:

  1. WHERE it stops -- construction, specify_initial_conditions,
     vc.initialize(), set_specifications, or the solve itself. These are
     different failures with different causes, and the sweep cannot tell
     them apart.

  2. WHETHER the primal solution is actually fine. A "failed" case may
     still have satisfied every constraint -- see vapor_compression_cubic_
     refstate_0903.py lines 1623-1632, which documents that the coupled
     solve at SH=SC=0 can reach a point with ZERO large constraint
     residuals and still be reported infeasible, because the
     log_mole_frac_tbub/tdew variables sit EXACTLY at 0.0 against a
     (None, 0] bound (correct for a pure fluid: log(1.0) = 0) which trips
     Ipopt's dual-feasibility check. This script prints the largest
     residual and lists the variables sitting on their bounds, so that
     situation is visible rather than inferred.

  3. WHERE THE FOUR CYCLE POINTS LANDED, plotted. The state values are read
     out whether or not the solve converged -- on a failure they show
     where the solver got stuck, which is the whole point.

The P-T plot is the useful one for the Tc-perturbation question, because
pressure and temperature need no reference state, so nothing can be
misaligned by an enthalpy datum. It draws:

  - the perturbed fluid's true saturation curve (Ambrose-Walton, from this
    case's own Tc/Pc/omega),
  - the Antoine curve that IDAES actually uses to SEED the bubble/dew
    calculation -- a single fit to real R32, shared unchanged by every
    perturbed case (phase_1_cubic_eos_validation_refstate_0903.py lines
    113-115), and
  - the four converged/stuck cycle points.

If the hypothesis is right, at Tc >= +3% those two curves separate sharply,
and the solver is being seeded at a pressure that belongs to a different
fluid than the one it is being asked to solve.

Nothing here writes to the sweep's own outputs, and nothing is modified in
p1cev.METHODS beyond adding one uniquely-named entry, same convention as
register_case() in the sweep.

Usage:
    python3 diagnose_one_sensitivity_point_0917.py

Edit the CONFIG block below to change which case is examined.

Author: Shilpa Narasimhan. Support: Claude AI.
"""

import os
import traceback

import numpy as np
import matplotlib.pyplot as plt
from pyomo.environ import value, Var, Constraint

import phase_1_cubic_eos_validation_refstate_0903 as p1cev
from vapor_compression_cubic_refstate_0903 import (
    SimpleVaporCompressionCycle as CubicCycle,
    Mode as CubicMode,
)

# =============================================================================
# CONFIG -- the one case to examine
# =============================================================================

TC_DEV_PCT = 3.0    # +3% is in the block that fails at EVERY ambient
PC_DEV_PCT = 0.0
TAMB = 20           # deg C

# Also run the unperturbed baseline (0%, 0%) for side-by-side contrast.
# Set False for a faster run if you only care about the failing case.
RUN_BASELINE_TOO = True

VERBOSE_SOLVER = True   # True -> print Ipopt's own iteration log

# NIST baseline, same literals the sweep pins itself to (see
# phase6_final_sensitivity_Tc_0904.py NIST_BASE). F/G are zeroed here to
# match the _0903 property file, which zeroes F/G for every method.
NIST_BASE = {
    "Pc": 57.82e5, "Tc": 351.3, "omega": 0.2769,
    "A": -6.098682, "B": 179.2200, "C": -122.3682, "D": 32.30207, "E": 0.491361,
    "F": 0.0, "G": 0.0,
}

FLUID = "R32"
M_R32 = 0.052024        # kg/mol
EVAP_SAT_C = -29.0      # evaporator saturation setpoint
COND_APPROACH = 9.0     # condenser sat = T_amb + 9

HERE = os.path.dirname(os.path.abspath(__file__))


# =============================================================================
# Saturation curves -- pure numpy, no IDAES, no reference state needed
# =============================================================================

def psat_ambrose_walton_bar(T_K, Tc, Pc_bar, omega):
    """Ambrose-Walton corresponding-states saturation pressure (Poling
    sec. 7-4) -- the perturbed fluid's TRUE vapor-pressure curve, i.e. what
    this case's own Tc/Pc/omega imply. Same correlation compare_cp_methods_*
    uses. Returns NaN at or above Tc."""
    T_K = np.asarray(T_K, dtype=float)
    Tr = T_K / Tc
    out = np.full(Tr.shape, np.nan)
    ok = Tr < 1.0
    tr = Tr[ok]
    tau = 1.0 - tr
    f0 = (-5.97616*tau + 1.29874*tau**1.5 - 0.60394*tau**2.5 - 1.06841*tau**5)/tr
    f1 = (-5.03365*tau + 1.11505*tau**1.5 - 5.41217*tau**2.5 - 7.46628*tau**5)/tr
    f2 = (-0.64771*tau + 2.41539*tau**1.5 - 4.26979*tau**2.5 + 3.25259*tau**5)/tr
    out[ok] = Pc_bar*np.exp(f0 + omega*f1 + omega**2*f2)
    return out


def psat_antoine_bar(T_K):
    """The Antoine curve IDAES actually uses as the bubble/dew INITIAL GUESS
    (pressure_sat_comp: NIST). Coefficients come from the property file, are
    a single fit to REAL R32, and are identical for every perturbed case --
    they are never re-fit when Tc/Pc move. That is the point of plotting it."""
    return 10.0**(p1cev.Antoine_A - p1cev.Antoine_B/(np.asarray(T_K, float) + p1cev.Antoine_C))


# =============================================================================
# Model inspection helpers
# =============================================================================

def max_constraint_residual(model):
    """Largest |body - target| over all active constraints. Same diagnostic
    the sweep uses. A SMALL value on a 'failed' case means the primal
    solution was fine and the failure was a certification failure."""
    worst, worst_name = 0.0, None
    for con in model.component_data_objects(Constraint, active=True, descend_into=True):
        try:
            body = value(con.body, exception=False)
            if body is None:
                continue
            lo = value(con.lower, exception=False) if con.lower is not None else None
            up = value(con.upper, exception=False) if con.upper is not None else None
            target = lo if (lo is not None and up is not None and lo == up) else 0.0
            resid = abs(body - target)
        except Exception:
            continue
        if resid > worst:
            worst, worst_name = resid, con.name
    return worst, worst_name


def vars_on_bounds(model, tol=1e-7, limit=25):
    """Unfixed variables sitting ON a bound at the current point. These are
    what make this problem hard for an interior-point solver: Ipopt wants to
    approach bounds from strictly inside, and here the SOLUTION is the
    boundary (log_mole_frac = log(1.0) = 0 against an upper bound of 0)."""
    hits = []
    for v in model.component_data_objects(Var, active=True, descend_into=True):
        if v.fixed or v.value is None:
            continue
        if v.lb is not None and abs(v.value - v.lb) <= tol:
            hits.append((v.name, v.value, "lb", v.lb))
        elif v.ub is not None and abs(v.value - v.ub) <= tol:
            hits.append((v.name, v.value, "ub", v.ub))
        if len(hits) >= limit:
            break
    return hits


def read_points(m):
    """The four cycle state points, read whether or not the solve converged.
    On a failure these show where the solver got stuck. Mirrors
    _cubic_point() in the sweep (h/s divided by M_R32 -> mass basis)."""
    units = {
        "evap_out":  m.fs.evaporator,
        "comp_out":  m.fs.compressor,
        "cond_out":  m.fs.condenser,
        "valve_out": m.fs.expansion_valve,
    }
    pts = {}
    for label, unit in units.items():
        st = unit.control_volume.properties_out[0]
        def g(expr):
            try:
                v = value(expr, exception=False)
                return float(v) if v is not None else np.nan
            except Exception:
                return np.nan
        pts[label] = {
            "T_K":   g(st.temperature),
            "P_Pa":  g(st.pressure),
            "h":     g(st.enth_mol)/M_R32 if not np.isnan(g(st.enth_mol)) else np.nan,
            "s":     g(st.entr_mol)/M_R32 if not np.isnan(g(st.entr_mol)) else np.nan,
            "x":     g(st.phase_frac["Vap"]),
        }
    return pts


# =============================================================================
# Run one case, stage by stage
# =============================================================================

def run_one(tc_dev, pc_dev, Tamb, verbose=VERBOSE_SOLVER):
    """Returns a dict describing everything that happened. Never raises --
    a failure at any stage is captured and reported, because WHICH stage
    fails is the information this script exists to produce."""
    Tc = NIST_BASE["Tc"] * (1 + tc_dev/100.0)
    Pc = NIST_BASE["Pc"] * (1 + pc_dev/100.0)
    name = f"DIAG_Tc{tc_dev:+.2f}pct_Pc{pc_dev:+.2f}pct"
    p1cev.METHODS[name] = {**NIST_BASE, "Tc": Tc, "Pc": Pc}

    out = {
        "name": name, "tc_dev": tc_dev, "pc_dev": pc_dev, "Tamb": Tamb,
        "Tc": Tc, "Pc": Pc, "omega": NIST_BASE["omega"],
        "stage_failed": None, "error": None,
        "cop": None, "converged": False,
        "points": None, "max_resid": None, "worst_con": None, "on_bounds": [],
    }

    Tcond_sat = Tamb + COND_APPROACH
    banner = f"  CASE  Tc {tc_dev:+.2f}%  Pc {pc_dev:+.2f}%   T_amb = {Tamb} C"
    print("\n" + "=" * 78)
    print(banner)
    print("=" * 78)
    print(f"  perturbed Tc = {Tc:.3f} K ({Tc-273.15:.2f} C),  Pc = {Pc/1e5:.4f} bar,  omega = {NIST_BASE['omega']}")

    # What the two saturation curves say at the two setpoints. If these
    # disagree badly, the solver is being seeded for the wrong fluid.
    for label, T_C in [("evaporator", EVAP_SAT_C), ("condenser", Tcond_sat)]:
        T_K = T_C + 273.15
        true_p = float(psat_ambrose_walton_bar(np.array([T_K]), Tc, Pc/1e5, NIST_BASE["omega"])[0])
        seed_p = float(psat_antoine_bar(np.array([T_K]))[0])
        ratio = seed_p/true_p if true_p and not np.isnan(true_p) else np.nan
        print(f"  {label:<11} setpoint {T_C:6.1f} C:  true Psat = {true_p:8.4f} bar | "
              f"Antoine seed = {seed_p:8.4f} bar | seed/true = {ratio:5.2f}x")

    vc = None
    try:
        print("\n  [1/5] constructing flowsheet + property package ...")
        vc = CubicCycle(FLUID, compressor_efficiency=0.9999,
                        mode=CubicMode.IMPROVED_TPX, method=name)
        print("        ok")

        print("  [2/5] specify_initial_conditions ...")
        vc.specify_initial_conditions(low_side_temperature=EVAP_SAT_C,
                                      high_side_temperature=Tcond_sat)
        print("        ok")

        print("  [3/5] vc.initialize() ...")
        vc.initialize(verbose=False)
        print("        ok")

        print("  [4/5] set_specifications ...")
        vc.set_specifications(ambient_temperature=Tamb, condenser_approach=COND_APPROACH,
                              evap_sat_temperature=EVAP_SAT_C, superheating=0, subcooling=0,
                              max_pressure_ratio=10)
        print("        ok")

        print("  [5/5] optimize_COP(initialize=True, optimize=False) -- square solve ...")
        cop, converged = vc.optimize_COP(verbose=verbose, initialize=True, optimize=False)
        out["cop"], out["converged"] = cop, converged
        print(f"        returned COP = {cop}, converged = {converged}")

    except Exception as e:
        stage = ("construction" if vc is None else "initialize/specify/solve")
        out["stage_failed"] = stage
        out["error"] = f"{type(e).__name__}: {e}"
        print(f"\n  !! STOPPED during {stage}")
        print(f"     {out['error']}")
        traceback.print_exc()

    # Post-mortem -- runs whether or not anything above succeeded.
    if vc is not None and getattr(vc, "model", None) is not None:
        m = vc.model
        try:
            out["max_resid"], out["worst_con"] = max_constraint_residual(m)
        except Exception:
            pass
        try:
            out["on_bounds"] = vars_on_bounds(m)
        except Exception:
            pass
        try:
            out["points"] = read_points(m)
        except Exception:
            pass

    print("\n  --- post-mortem -------------------------------------------------")
    if out["max_resid"] is not None:
        print(f"  largest constraint residual : {out['max_resid']:.3e}")
        print(f"  worst constraint            : {out['worst_con']}")
        if out["max_resid"] < 1e-4 and not out["converged"]:
            print("  NOTE: residual is tiny but the case is flagged NOT converged --")
            print("        the primal solution is essentially satisfied. This is a")
            print("        certification failure, not a failure to find the point.")
    if out["on_bounds"]:
        print(f"  variables sitting ON a bound (first {len(out['on_bounds'])}):")
        for nm, val, side, b in out["on_bounds"]:
            print(f"      {nm:<62} = {val: .6e}  ({side} = {b})")
    if out["points"]:
        print("  cycle points:")
        print(f"      {'point':<11}{'T [C]':>10}{'P [bar]':>11}{'h [kJ/kg]':>12}"
              f"{'s [kJ/kg/K]':>13}{'vap frac':>10}")
        for label in ["evap_out", "comp_out", "cond_out", "valve_out"]:
            p = out["points"][label]
            print(f"      {label:<11}{p['T_K']-273.15:10.2f}{p['P_Pa']/1e5:11.4f}"
                  f"{p['h']/1e3:12.2f}{p['s']/1e3:13.4f}{p['x']:10.4f}")
    return out


# =============================================================================
# Plots
# =============================================================================

def plot_case(res, fname_prefix):
    """P-T diagram (reference-state-free, so nothing can be misaligned) plus
    the cycle polygon in p-h."""
    Tc, Pc_bar, omega = res["Tc"], res["Pc"]/1e5, res["omega"]
    Tamb = res["Tamb"]

    T_C = np.linspace(-60.0, min(Tc - 273.15 - 0.5, 130.0), 400)
    T_K = T_C + 273.15
    p_true = psat_ambrose_walton_bar(T_K, Tc, Pc_bar, omega)
    p_seed = psat_antoine_bar(T_K)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(T_C, p_true, "-",  color="tab:blue", lw=2,
            label=f"true saturation curve (this case: Tc={Tc-273.15:.1f} C)")
    ax.plot(T_C, p_seed, "--", color="tab:orange", lw=2,
            label="Antoine seed curve (fit to REAL R32, never perturbed)")
    ax.axvline(Tc - 273.15, color="tab:blue", ls=":", alpha=0.6)
    ax.annotate("Tc", (Tc - 273.15, ax.get_ylim()[1]), color="tab:blue",
                ha="center", va="top", fontsize=9)

    # the two target setpoints the cycle must hit
    for T_set, lab in [(EVAP_SAT_C, "evaporator setpoint"),
                       (Tamb + COND_APPROACH, "condenser setpoint")]:
        ax.axvline(T_set, color="grey", ls=":", alpha=0.8)
        ax.text(T_set, 0.04, f" {lab} ({T_set:.0f} C)", rotation=90, fontsize=8,
                color="grey", ha="right", va="bottom", transform=
                ax.get_xaxis_transform())

    if res["points"]:
        order = ["evap_out", "comp_out", "cond_out", "valve_out"]
        xs = [res["points"][k]["T_K"] - 273.15 for k in order]
        ys = [res["points"][k]["P_Pa"]/1e5 for k in order]
        ax.plot(xs + [xs[0]], ys + [ys[0]], "-o", color="tab:red", ms=7,
                lw=1.5, label="cycle points (where the solver ended up)")
        # NOTE: evap_out and valve_out both sit at the evaporator pressure and
        # saturation temperature, so they land on TOP of each other in P-T --
        # that is real, not a plotting bug. Labels are pushed in different
        # directions so both stay readable; the p-h plot separates them
        # properly, since they differ in enthalpy, not in T or P.
        offsets = {"evap_out": (6, 10), "comp_out": (6, 8),
                   "cond_out": (6, 8), "valve_out": (6, -14)}
        for k, x, y in zip(order, xs, ys):
            ax.annotate(k, (x, y), textcoords="offset points",
                        xytext=offsets[k], fontsize=8, color="tab:red")

    ax.set_yscale("log")
    ax.set_xlabel("T [deg C]")
    ax.set_ylabel("P [bar]")
    status = "CONVERGED" if res["converged"] else "NOT CONVERGED"
    ax.set_title(f"P-T: Tc {res['tc_dev']:+.2f}%, Pc {res['pc_dev']:+.2f}%, "
                 f"T_amb={Tamb} C  [{status}]")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    plt.tight_layout()
    p1 = os.path.join(HERE, f"{fname_prefix}_PT.png")
    plt.savefig(p1, dpi=150)
    print(f"  wrote {os.path.basename(p1)}")

    # p-h cycle polygon (IDAES's own enthalpy basis, so no dome overlaid --
    # the shape alone shows whether the cycle is sane)
    if res["points"] and not any(np.isnan(res["points"][k]["h"]) for k in res["points"]):
        order = ["evap_out", "comp_out", "cond_out", "valve_out"]
        hs = [res["points"][k]["h"]/1e3 for k in order]
        ps = [res["points"][k]["P_Pa"]/1e5 for k in order]
        fig, ax = plt.subplots(figsize=(9, 6.5))
        ax.plot(hs + [hs[0]], ps + [ps[0]], "-o", color="tab:red", ms=7, lw=1.5)
        for k, x, y in zip(order, hs, ps):
            ax.annotate(f"{k}\nx={res['points'][k]['x']:.3f}", (x, y),
                        textcoords="offset points", xytext=(8, 6), fontsize=8)
        ax.set_yscale("log")
        ax.set_xlabel("h [kJ/kg]  (IDAES reference state)")
        ax.set_ylabel("P [bar]")
        ax.set_title(f"p-h cycle points: Tc {res['tc_dev']:+.2f}%, "
                     f"T_amb={Tamb} C  [{status}]")
        ax.grid(True, which="both", ls=":", alpha=0.6)
        plt.tight_layout()
        p2 = os.path.join(HERE, f"{fname_prefix}_ph.png")
        plt.savefig(p2, dpi=150)
        print(f"  wrote {os.path.basename(p2)}")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    results = []

    res = run_one(TC_DEV_PCT, PC_DEV_PCT, TAMB)
    plot_case(res, f"diag_Tc{TC_DEV_PCT:+.2f}_Pc{PC_DEV_PCT:+.2f}_Tamb{TAMB}")
    results.append(res)

    if RUN_BASELINE_TOO and (TC_DEV_PCT, PC_DEV_PCT) != (0.0, 0.0):
        base = run_one(0.0, 0.0, TAMB)
        plot_case(base, f"diag_Tc+0.00_Pc+0.00_Tamb{TAMB}")
        results.append(base)

    print("\n" + "=" * 78)
    print("  SUMMARY")
    print("=" * 78)
    print(f"  {'case':<26}{'converged':>11}{'COP':>10}{'max resid':>13}  stage failed")
    for r in results:
        tag = f"Tc{r['tc_dev']:+.2f}% Pc{r['pc_dev']:+.2f}% @{r['Tamb']}C"
        cop = f"{r['cop']:.4f}" if isinstance(r["cop"], float) else "--"
        mr = f"{r['max_resid']:.2e}" if r["max_resid"] is not None else "--"
        print(f"  {tag:<26}{str(r['converged']):>11}{cop:>10}{mr:>13}  {r['stage_failed'] or '-'}")
    print()
