"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 09/14/2026

What this does
    Sweeps ambient temperature for the cascade with R515B in the hot loop
    and CO2 in the cold loop. R134a/CO2 runs first as the baseline. Both use
    the same setpoints and the same solve steps, so only the hot-loop fluid
    differs.

Why it imports a different engine
    R515B is a blend. It is not in CoolProp or IDAES general_helmholtz.
    It uses this repo's own package (R515BParameterBlock), which only has
    P-h state variables. So CascadeCycle comes from
    vapor_compression_cascade_r515b.py, not vapor_compression_cascade.py.

Setpoints
    PLR 0.75, CD 0.13       part-load correction (CD > 0 keeps it on)
    CO2 evaporating         -29 C
    Hot condenser           ambient + 9 C
    Cascade gap             3 C   (T_CAS,C - T_CAS,E)
    Superheat / subcool     3 C floor, 8 C cap (CO2 subcool floor 0 C)
    Cascade cold end        CO2 out >= T_CAS,E + 1 C (no temperature cross)
    Head pressure control   condensing never below 21 C (T_COND_MIN_C), so
                            every ambient <= 12 C runs as a 12 C point
    Compressor efficiency   0.80, both loops

Reading the results
    - COP is maximized, so superheat/subcool can differ between fluids.
      That is expected.
    - Cold-loop flow is fixed at 1 kg/s. Hot-loop flow is solved, so
      hot-loop duty and work are NOT per kg.
    - Every ambient point is built fresh. No point depends on another.
    - Carnot COP uses -29 C and ambient + 9 C, so it is the same for both
      fluids. It is compared to full-load COP, not part-load.

Outputs (written to this folder; <case> = 21C or nolimit)
    sweep_cascade_R515B_<case>.csv
    cop_vs_ambient_plr_benchmark_combined_r515b_<case>.csv
    cop_vs_ambient_plr_benchmark_comparison_r515b_<case>.png / .pdf

Run
    python run_cascade_sweep_r515b.py            (21 C limit)
    python run_cascade_sweep_r515b.py nolimit    (no condensing limit)
"""

import csv
import sys

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from pyomo.environ import value

# This file lives in Cascade/, alongside vapor_compression_cascade_r515b.py
# -- that import is a plain local one, no extra sys.path setup needed for
# it. R515B's own custom property package lives under repo-root's
# R515B_idaes_package/ instead, but vapor_compression_cascade_r515b.py
# adds that to sys.path itself (relative to its own location), so this
# runner doesn't need to duplicate that setup.
from vapor_compression_cascade_r515b import CascadeCycle
from physics_checks import cascade_checks
from heat_balance import cascade_balance, cascade_balance_keys


# ---- Setpoints (collaborator-sourced) ----
PLR_VALUE = 0.75
CD_VALUE = 0.13
EVAP_SAT_TARGET_C = -29.0
COND_APPROACH_C = 9.0
SUPERHEAT_C = 3.0
SUBCOOL_C = 3.0
SUPERHEAT_MAX_C = 8.0
SUBCOOL_MAX_C = 8.0
COMPRESSOR_EFFICIENCY = 0.80
MAX_PRESSURE_RATIO = 20.0
# Internal hot-evap/cold-cond approach inside the cascade HX itself --
# distinct from COND_APPROACH_C (the hot loop's ambient-side approach).
CASCADE_APPROACH_DT_C = 3.0
# CO2 side of the cascade HX: subcool floor 0 K and 1 K cold-end approach,
# so CO2 never leaves colder than T_CAS,E (10/10).
SUBCOOL_COLD_C = 0.0
CASCADE_COLD_END_DT_C = 1.0

# Ascending order, matching validate_cascade.py's BENCHMARK_AMBIENTS_C
# convention. Order no longer affects convergence -- each point is solved
# from its own fresh, independently-seeded CascadeCycle.
AMBIENT_SWEEP_C = np.arange(-30, 55, 5)   # -30 to 50 C (collaborator range)

# R134a baseline OFF: certified R134a/CO2 numbers already exist from the main
# engine (sweep_amb_R134a.csv, 10/10). This engine lacks the main engine's
# fallbacks, so R134a stalls here at low ambient. Turn on only to compare engines.
RUN_R134A_BASELINE = False

# Head pressure control: real systems hold a minimum condensing temperature.
# 21 C = Copeland's lowest condenser temperature
# (https://copeland.custhelp.com/app/answers/detail/a_id/3278/; also the
# 70 F floating-head minimum in Copeland's low-condensing webinar).
# Zhuorui (10/10): use a typical minimum condensing temperature.
T_COND_MIN_C = 21.0

# Case from the command line: "21C" (default, head pressure control) or
# "nolimit" (condensing = ambient + 9 C everywhere). Tags the output files.
CASE = sys.argv[1] if len(sys.argv) > 1 else "21C"
assert CASE in ("21C", "nolimit"), "case must be 21C or nolimit"
if CASE == "nolimit":
    T_COND_MIN_C = None


def effective_ambient(ambient_c):
    """Ambient the model sees, so condensing (ambient + approach) never
    drops below T_COND_MIN_C. Below 12 C every point looks like 12 C.
    No limit (T_COND_MIN_C = None): the real ambient.
    Results and plots still use the real ambient."""
    if T_COND_MIN_C is None:
        return ambient_c
    return max(ambient_c, T_COND_MIN_C - COND_APPROACH_C)

def _build_cycle(fluid_name):
    """Build one empty cascade for the chosen hot-loop fluid.

    fluid_name: "R515B" or "R134a". The cold loop is always CO2.
    "property_block": "r515b" tells the engine to use R515B's own package
    (and P-h mode). Anything else falls through to R134a.
    Returns (cycle, fluid_name).
    """
    cold_fluid = {
        "gh_component": "co2", "coolprop_name": "CO2", "custom_json": None,
        "efficiency": COMPRESSOR_EFFICIENCY,
        "superheat_max": SUPERHEAT_MAX_C, "subcool_max": SUBCOOL_MAX_C,
        "Tmax_C": None,
    }

    if fluid_name == "R515B":
        # Not a CoolProp/general_helmholtz fluid -- "property_block":
        # "r515b" tells CascadeCycle to build R515BParameterBlock() for
        # this role instead, forcing the whole cascade into PH mode.
        # gh_component/coolprop_name are unused placeholders here (both
        # only get read on the non-R515B branch of the engine's property
        # block construction and CoolProp-bound-derivation code paths).
        hot_fluid = {
            "property_block": "r515b",
            "gh_component": "r515b", "coolprop_name": "R515B", "custom_json": None,
            "efficiency": COMPRESSOR_EFFICIENCY,
            "superheat_max": SUPERHEAT_MAX_C, "subcool_max": SUBCOOL_MAX_C,
            "Tmax_C": None,
        }
    else:  # "R134a" baseline
        hot_fluid = {
            "gh_component": "r134a", "coolprop_name": "R134a", "custom_json": None,
            "efficiency": COMPRESSOR_EFFICIENCY,
            "superheat_max": SUPERHEAT_MAX_C, "subcool_max": SUBCOOL_MAX_C,
            "Tmax_C": None,
        }

    cycle = CascadeCycle(fluids={"hot": hot_fluid, "cold": cold_fluid})
    return cycle, fluid_name


C_to_K = 273.15


def _report_loop_ops(fs, role):
    """Read one loop's solved state into a flat dict.

    role: "hot" or "cold". Gives inlet/outlet T and P for all four units,
    duty and work, and the actual superheat/subcool.
    Why: so every COP in the CSV comes with the operating point behind it.
    """
    evaporator = getattr(fs, f"{role}_evaporator")
    compressor = getattr(fs, f"{role}_compressor")
    condenser = getattr(fs, f"{role}_condenser")
    expansion_valve = getattr(fs, f"{role}_expansion_valve")

    # NOTE: in PH mode, ports only expose the actual state variables
    # (pressure, enth_mass, flow_mass) -- NOT temperature. Temperature is a
    # derived quantity, only available via control_volume.properties_in/out.
    def T_in(unit):
        return value(unit.control_volume.properties_in[0].temperature) - C_to_K

    def T_out(unit):
        return value(unit.control_volume.properties_out[0].temperature) - C_to_K

    def P(port):
        return value(port.pressure[0]) / 1000.0  # kPa

    evap_T_out = value(evaporator.control_volume.properties_out[0].temperature)
    evap_T_sat = value(evaporator.control_volume.properties_out[0].temperature_sat)
    cond_T_out = value(condenser.control_volume.properties_out[0].temperature)
    cond_T_sat = value(condenser.control_volume.properties_out[0].temperature_sat)

    return {
        "evap_in_T_C": T_in(evaporator), "evap_in_P_kPa": P(evaporator.inlet),
        "evap_out_T_C": T_out(evaporator), "evap_out_P_kPa": P(evaporator.outlet),
        "evap_T_sat_C": evap_T_sat - C_to_K,
        "evap_superheat_actual_K": evap_T_out - evap_T_sat,
        "evap_heat_duty_per_kg": value(evaporator.heat_duty[0]),
        "comp_in_T_C": T_in(compressor), "comp_in_P_kPa": P(compressor.inlet),
        "comp_out_T_C": T_out(compressor), "comp_out_P_kPa": P(compressor.outlet),
        "comp_ratioP": value(compressor.ratioP[0]),
        "comp_work_per_kg": value(compressor.work_mechanical[0]),
        "cond_in_T_C": T_in(condenser), "cond_in_P_kPa": P(condenser.inlet),
        "cond_out_T_C": T_out(condenser), "cond_out_P_kPa": P(condenser.outlet),
        "cond_T_sat_C": cond_T_sat - C_to_K,
        "cond_subcool_actual_K": cond_T_sat - cond_T_out,
        "cond_heat_duty_per_kg": value(condenser.heat_duty[0]),
        "valve_in_T_C": T_in(expansion_valve), "valve_in_P_kPa": P(expansion_valve.inlet),
        "valve_out_T_C": T_out(expansion_valve), "valve_out_P_kPa": P(expansion_valve.outlet),
    }


def _report_unit_ops(cycle):
    """Both loops' states in one dict, keys prefixed hot_ / cold_.

    Why: one flat dict = one CSV row.
    """
    fs = cycle.model.fs
    out = {}
    for role in ("hot", "cold"):
        for k, v in _report_loop_ops(fs, role).items():
            out[f"{role}_{k}"] = v
    out.update(cascade_balance(cycle))      # mass flows, h of every stream, heat balance (10/14)
    return out


def _carnot_cop(ambient_c):
    """Carnot COP = T_cold / (T_hot - T_cold), in K.

    T_cold = -29 C (CO2 evaporating). T_hot = ambient + 9 C.
    Same for both fluids. It is the reversible limit, not "efficiency = 1".
    """
    t_cold_k = EVAP_SAT_TARGET_C + C_to_K
    t_hot_k = effective_ambient(ambient_c) + COND_APPROACH_C + C_to_K   # condensing T (head pressure control)
    return t_cold_k / (t_hot_k - t_cold_k)


_UNIT_OP_FIELDS_BASE = [
    "evap_in_T_C", "evap_in_P_kPa", "evap_out_T_C", "evap_out_P_kPa", "evap_T_sat_C",
    "evap_superheat_actual_K", "evap_heat_duty_per_kg",
    "comp_in_T_C", "comp_in_P_kPa", "comp_out_T_C", "comp_out_P_kPa", "comp_ratioP", "comp_work_per_kg",
    "cond_in_T_C", "cond_in_P_kPa", "cond_out_T_C", "cond_out_P_kPa", "cond_T_sat_C",
    "cond_subcool_actual_K", "cond_heat_duty_per_kg",
    "valve_in_T_C", "valve_in_P_kPa", "valve_out_T_C", "valve_out_P_kPa",
]
# Doubled for both loops -- matches _report_unit_ops's hot_/cold_ prefixed keys.
_UNIT_OP_FIELDS = [f"{role}_{f}" for role in ("hot", "cold") for f in _UNIT_OP_FIELDS_BASE]
# + heat balance columns (mass flows, enthalpies, duties in kW); skip names already present
_UNIT_OP_FIELDS += [k for k in cascade_balance_keys() if k not in _UNIT_OP_FIELDS]


def _run_fluid(fluid_name, out_stem):
    """Run the full ambient sweep for one hot-loop fluid.

    Steps per ambient point:
        1. Build a fresh cascade and seed it at this ambient.
        2. Initialize, then set specs (hot loop gets this ambient).
        3. optimize_COP. A crash counts as not converged.
        4. Record COPs and the operating point (NaN if not converged).
    Then a retry pass: each failed point is re-solved starting from the
    nearest converged point's solution (warm start). Why: at 15 C R515B got
    stuck from a cold start (condenser outlet on the saturation line) but
    converged in seconds from the 20 C solution (10/10). Points that
    converged the first time are left exactly as they are.
    Every converged point also goes through physics_checks.py; the CSV
    column checks_failed lists any check it fails (empty = all pass).
    Writes <out_stem>.csv and prints one line per point.
    Returns (ambients, full COPs, part COPs, converged flags, solved dict).
    """
    # Cold loop's spec dict: fixed load-side target (evap_sat_temperature),
    # no ambient_temperature -- that belongs to the hot loop only.
    cold_kwargs = dict(
        evap_sat_temperature=EVAP_SAT_TARGET_C, superheating=SUPERHEAT_C,
        superheating_max=SUPERHEAT_MAX_C, subcooling=SUBCOOL_COLD_C,
        subcooling_max=SUBCOOL_MAX_C, max_pressure_ratio=MAX_PRESSURE_RATIO,
    )
    # Hot loop's spec dict base -- ambient_temperature gets added per-point
    # inside the sweep loop below, since that's what's actually swept.
    hot_kwargs_base = dict(
        condenser_approach=COND_APPROACH_C, superheating=SUPERHEAT_C,
        superheating_max=SUPERHEAT_MAX_C, subcooling=SUBCOOL_C,
        subcooling_max=SUBCOOL_MAX_C, max_pressure_ratio=MAX_PRESSURE_RATIO,
    )

    solved = {}
    good_cycles = {}     # ambient -> converged cycle, for the warm-start retry
    warm_from = {}       # ambient -> neighbour ambient it was warm-started from
    checks = {}          # ambient -> (info, failed) from physics_checks.py
    same_cycle = {}      # ambient the model sees -> first ambient solved there (head pressure control)
    for ambient_c in AMBIENT_SWEEP_C:
        # Fresh CascadeCycle for every point, seeded at THAT point's own
        # ambient temperature -- no dependency on any other point's
        # solution (matches validate_cascade.py's _build pattern).
        amb_eff = effective_ambient(float(ambient_c))      # head pressure control
        # Head pressure control: every ambient <= 12 C is the SAME cycle (model sees 12 C).
        # Reuse the first solved one instead of solving it again (saves ~8 slow solves).
        if amb_eff in same_cycle:
            a0 = same_cycle[amb_eff]
            ok0, c0, cp0, _, _, uo0 = solved[a0]
            cc = _carnot_cop(ambient_c)
            solved[ambient_c] = (ok0, c0, cp0, cc, c0 / cc, uo0)
            checks[ambient_c] = checks[a0]
            print(f"{fluid_name:>10s} | Ambient {ambient_c:>5.1f} C | converged=1 | same cycle as {a0:.0f} C "
                  f"(model sees {amb_eff:.0f} C) | COP_full={c0:.4f}")
            continue
        cycle, fluid_used = _build_cycle(fluid_name)
        hot_kwargs = dict(hot_kwargs_base)
        hot_kwargs["ambient_temperature"] = amb_eff

        # Whole point in one try: a crash in seeding, init, specs or solve
        # marks THIS point as not converged instead of stopping the sweep.
        try:
            cycle.specify_initial_conditions(hot_ambient_C=amb_eff, cold_evap_C=EVAP_SAT_TARGET_C,
                                             cond_approach_C=COND_APPROACH_C)
            cycle.initialize(verbose=False)
            cycle.set_specifications(hot=hot_kwargs, cold=cold_kwargs,
                                      cascade_approach_dT=CASCADE_APPROACH_DT_C,
                                      cascade_cold_end_dT=CASCADE_COLD_END_DT_C,
                                      plr=PLR_VALUE, cd=CD_VALUE)
            cop_full, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=True)
        except Exception as exc:
            print(f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | ERROR: {type(exc).__name__}: {exc}")
            cop_full, ok = float("nan"), False

        cop_part = cycle.get_part_load_cop() if ok else float("nan")
        cop_carnot = _carnot_cop(ambient_c)
        second_law_eff = (cop_full / cop_carnot) if ok else float("nan")
        unit_ops = _report_unit_ops(cycle) if ok else {k: float("nan") for k in _UNIT_OP_FIELDS}
        solved[ambient_c] = (ok, cop_full if ok else float("nan"), cop_part, cop_carnot, second_law_eff, unit_ops)
        if ok:
            good_cycles[ambient_c] = cycle
            same_cycle.setdefault(amb_eff, ambient_c)
            checks[ambient_c] = cascade_checks(cycle, amb_eff, sc_cold=(SUBCOOL_COLD_C, SUBCOOL_MAX_C),
                                               cold_end_dT=CASCADE_COLD_END_DT_C)

        if ok:
            print(
                f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | converged=1 | "
                f"COP_full={cop_full:.4f} | COP_part={cop_part:.4f} | COP_carnot={cop_carnot:.4f} | "
                f"2nd_law_eff={second_law_eff*100:.1f}% | "
                f"SH_actual(cold_evap)={unit_ops['cold_evap_superheat_actual_K']:.2f}K | "
                f"SC_actual(hot_cond)={unit_ops['hot_cond_subcool_actual_K']:.2f}K | "
                f"hot_ratioP={unit_ops['hot_comp_ratioP']:.2f} | cold_ratioP={unit_ops['cold_comp_ratioP']:.2f}"
            )
        else:
            print(f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | converged=0 | "
                  f"COP_carnot={cop_carnot:.4f} (thermodynamic bound is defined regardless of convergence)")

    # --- Retry pass: warm start failed points from the nearest converged one ---
    for ambient_c in AMBIENT_SWEEP_C:
        if solved[ambient_c][0] or not good_cycles:
            continue
        nb = min(good_cycles, key=lambda a: (abs(a - ambient_c), -a))   # nearest; warmer wins a tie
        cycle = good_cycles[nb]      # note: this neighbour's model is moved to the new ambient
        try:
            cycle.model.fs.hot_condenser.ambient_T.set_value(effective_ambient(float(ambient_c)) + C_to_K)
            cop_full, ok = cycle.optimize_COP(verbose=False, initialize=False, optimize=True)
        except Exception as exc:
            print(f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | warm-start ERROR: {type(exc).__name__}: {exc}")
            cop_full, ok = float("nan"), False
        if not ok:
            print(f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | converged=0 (warm start from {nb:.1f} C also failed)")
            continue
        cop_part = cycle.get_part_load_cop()
        cop_carnot = _carnot_cop(ambient_c)
        unit_ops = _report_unit_ops(cycle)
        solved[ambient_c] = (True, cop_full, cop_part, cop_carnot, cop_full / cop_carnot, unit_ops)
        warm_from[ambient_c] = nb
        amb_w = effective_ambient(float(ambient_c))
        checks[ambient_c] = cascade_checks(cycle, amb_w, sc_cold=(SUBCOOL_COLD_C, SUBCOOL_MAX_C),
                                           cold_end_dT=CASCADE_COLD_END_DT_C)
        print(f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | converged=1 (warm start from {nb:.1f} C) | "
              f"COP_full={cop_full:.4f} | COP_part={cop_part:.4f}")

    ambient_sorted = sorted(solved.keys())
    cop_full_vals = [solved[a][1] for a in ambient_sorted]
    cop_part_vals = [solved[a][2] for a in ambient_sorted]
    converged_vals = [solved[a][0] for a in ambient_sorted]

    out_csv = f"{out_stem}.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ambient_C", "cop_full", "cop_part", "cop_carnot", "second_law_efficiency",
                     "converged", "warm_start_from_C", "fluid_name_used",
                     "checks_failed", "T_cas_e_C", "hot_P_evap_kPa", "vacuum", "energy_rel"] + _UNIT_OP_FIELDS)
        for ta in ambient_sorted:
            ok, cfull, cpart, ccarnot, seff, unit_ops = solved[ta]
            info, failed = checks.get(ta, ({}, ["not converged"]))
            w.writerow([ta, cfull, cpart, ccarnot, seff, int(ok), warm_from.get(ta, ""), fluid_used,
                        " ".join(failed), info.get("T_cas_e_C", ""), info.get("hot_P_evap_kPa", ""),
                        info.get("vacuum", ""), info.get("energy_rel", "")]
                       + [unit_ops[k] for k in _UNIT_OP_FIELDS])

    print(f"Saved: {out_csv}")
    print(f"Converged points: {sum(converged_vals)}/{len(ambient_sorted)}")
    bad = {a: f for a, (i, f) in checks.items() if f}
    print(f"Physics checks failed at: {bad if bad else 'none'}")
    return np.array(ambient_sorted), cop_full_vals, cop_part_vals, converged_vals, solved


def _print_full_report(fluid_solved, tol_K=0.5):
    """Print COP and operating conditions for every ambient, every fluid.

    fluid_solved: {fluid label: results dict from _run_fluid}.
    The first fluid is the baseline. Superheat/subcool differences vs the
    baseline are printed as facts, not errors (each fluid optimizes its own).
    """
    labels = list(fluid_solved.keys())
    baseline_label = labels[0]
    all_ambients = sorted(set().union(*[set(d.keys()) for d in fluid_solved.values()]))

    print("\n" + "=" * 100)
    print(f"FULL OPERATING-CONDITION REPORT -- every point, {len(labels)} fluid(s): {', '.join(labels)}")
    print("=" * 100)

    empty_row = (False, float("nan"), float("nan"), float("nan"), float("nan"), {})

    for amb in all_ambients:
        rows = {}
        for label in labels:
            ok, cfull, cpart, ccarnot, seff, uo = fluid_solved[label].get(
                amb, (False, float("nan"), float("nan"), _carnot_cop(amb), float("nan"), {}))
            rows[label] = (ok, cfull, cpart, ccarnot, seff, uo)

        any_carnot = next((r[3] for r in rows.values() if not np.isnan(r[3])), _carnot_cop(amb))
        print(f"\n--- Ambient = {amb:.1f} C --- COP_carnot = {any_carnot:.4f} (same for all fluids) ---")

        for label in labels:
            ok, cfull, cpart, ccarnot, ceff, uo = rows[label]
            if not ok:
                print(f"  {label:>8s} | NOT CONVERGED")
                continue
            print(
                f"  {label:>8s} | COP_full={cfull:.4f} COP_part={cpart:.4f} "
                f"2nd_law_eff={ceff*100:.1f}% | "
                f"cold_evap[in {uo['cold_evap_in_T_C']:.2f}C/{uo['cold_evap_in_P_kPa']:.1f}kPa -> "
                f"out {uo['cold_evap_out_T_C']:.2f}C/{uo['cold_evap_out_P_kPa']:.1f}kPa, "
                f"Tsat={uo['cold_evap_T_sat_C']:.2f}C, SH={uo['cold_evap_superheat_actual_K']:.2f}K] | "
                f"hot_cond[in {uo['hot_cond_in_T_C']:.2f}C/{uo['hot_cond_in_P_kPa']:.1f}kPa -> "
                f"out {uo['hot_cond_out_T_C']:.2f}C/{uo['hot_cond_out_P_kPa']:.1f}kPa, "
                f"Tsat={uo['hot_cond_T_sat_C']:.2f}C, SC={uo['hot_cond_subcool_actual_K']:.2f}K] | "
                f"hot_comp[ratioP={uo['hot_comp_ratioP']:.2f}] "
                f"cold_comp[ratioP={uo['cold_comp_ratioP']:.2f}]"
            )

        ok_base, _, _, _, _, uo_base = rows[baseline_label]
        if ok_base and len(labels) > 1:
            for label in labels[1:]:
                ok_l, _, _, _, _, uo_l = rows[label]
                if not ok_l:
                    print(f"  -> {label} vs {baseline_label}: comparison not possible ({label} did not converge)")
                    continue
                dsh = uo_l["cold_evap_superheat_actual_K"] - uo_base["cold_evap_superheat_actual_K"]
                dsc = uo_l["hot_cond_subcool_actual_K"] - uo_base["hot_cond_subcool_actual_K"]
                print(f"  -> {label} vs {baseline_label}: each fluid's own COP-maximizing point: "
                      f"dSH={dsh:+.2f}K, dSC={dsc:+.2f}K (both allowed to independently land anywhere in [3,8]K)")
        elif len(labels) > 1:
            print(f"  -> comparison vs {baseline_label} not possible ({baseline_label} did not converge)")

    print("\n" + "=" * 100)


def _write_combined_csv(fluid_solved, out_path="cop_vs_ambient_plr_benchmark_combined.csv"):
    """One CSV row per ambient, all fluids side by side.

    Columns per fluid: converged flag, COPs, second-law efficiency, and the
    full operating point. Last columns: superheat/subcool difference vs the
    baseline fluid (informational).
    """
    labels = list(fluid_solved.keys())
    baseline_label = labels[0]
    all_ambients = sorted(set().union(*[set(d.keys()) for d in fluid_solved.values()]))
    empty_uo = {k: float("nan") for k in _UNIT_OP_FIELDS}
    empty_row = (False, float("nan"), float("nan"), float("nan"), float("nan"), empty_uo)

    header = ["ambient_C", "cop_carnot"]
    for label in labels:
        prefix = f"{label.lower()}_"
        header += [f"{prefix}converged", f"{prefix}cop_full", f"{prefix}cop_part", f"{prefix}second_law_efficiency"]
        header += [f"{prefix}{f}" for f in _UNIT_OP_FIELDS]
    for label in labels[1:]:
        header += [f"superheat_delta_K_{label.lower()}_vs_{baseline_label.lower()}",
                   f"subcool_delta_K_{label.lower()}_vs_{baseline_label.lower()}"]

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for amb in all_ambients:
            rows = {}
            for label in labels:
                ok, cfull, cpart, ccarnot, seff, uo = fluid_solved[label].get(
                    amb, (False, float("nan"), float("nan"), _carnot_cop(amb), float("nan"), empty_uo))
                rows[label] = (ok, cfull, cpart, ccarnot, seff, uo)

            any_carnot = next((r[3] for r in rows.values() if not np.isnan(r[3])), _carnot_cop(amb))
            row = [amb, any_carnot]
            for label in labels:
                ok, cfull, cpart, _, seff, uo = rows[label]
                row += [int(ok), cfull, cpart, seff] + [uo[k] for k in _UNIT_OP_FIELDS]

            ok_base, _, _, _, _, uo_base = rows[baseline_label]
            for label in labels[1:]:
                ok_l, _, _, _, _, uo_l = rows[label]
                if ok_base and ok_l:
                    dsh = uo_l["cold_evap_superheat_actual_K"] - uo_base["cold_evap_superheat_actual_K"]
                    dsc = uo_l["hot_cond_subcool_actual_K"] - uo_base["hot_cond_subcool_actual_K"]
                else:
                    dsh, dsc = float("nan"), float("nan")
                row += [dsh, dsc]

            w.writerow(row)

    print(f"Saved: {out_path}")


def _plot_comparison(results):
    """Plot part-load COP vs ambient for each fluid, plus Carnot.

    Only converged points are drawn. Full-load COP is left off the plot
    (it is in the CSV and the report).
    """
    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=160)
    for fluid_used, (ambient_c, cop_full_vals, cop_part_vals, converged_vals, _solved) in results.items():
        ok_mask = np.array(converged_vals, dtype=bool) & np.isfinite(cop_part_vals)
        if np.any(ok_mask):
            # Part-load COP only -- full-load COP is still computed and
            # reported (console report + combined CSV, and used for
            # second_law_efficiency = COP_full/COP_carnot), just no longer
            # plotted here: the full-load dashed lines made the chart hard
            # to read (per user feedback) without adding plot-only value,
            # since the same numbers are already in the report/CSV.
            ax.plot(
                ambient_c[ok_mask],
                np.array(cop_part_vals)[ok_mask],
                marker="s",
                linewidth=2.2,
                label=f"{fluid_used} (part-load)",
            )

    # Carnot COP is fluid-independent (same T_evap_sat/T_cond_sat targets
    # for both fluids), so plot it once using either result's ambient axis.
    any_ambient_c = next(iter(results.values()))[0]
    carnot_vals = [_carnot_cop(a) for a in any_ambient_c]
    ax.plot(any_ambient_c, carnot_vals, linestyle=":", linewidth=2.0, color="black", label="Carnot (reversible limit)")

    ax.set_title("PLR benchmark (fresh build per point): R134a vs R515B vs Carnot")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"cop_vs_ambient_plr_benchmark_comparison_r515b_{CASE}.png", bbox_inches="tight")
    fig.savefig(f"cop_vs_ambient_plr_benchmark_comparison_r515b_{CASE}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved: cop_vs_ambient_plr_benchmark_comparison_r515b.png")
    print("Saved: cop_vs_ambient_plr_benchmark_comparison_r515b.pdf")


def main():
    """R134a baseline first, then R515B. Then report, combined CSV, plot."""
    results = {}
    if RUN_R134A_BASELINE:
        results["R134a"] = _run_fluid("R134a", "cop_vs_ambient_plr_benchmark_r134a_for_r515b")
    results["R515B"] = _run_fluid("R515B", f"sweep_cascade_R515B_{CASE}")

    fluid_solved = {label: r[4] for label, r in results.items()}
    _print_full_report(fluid_solved)
    _write_combined_csv(fluid_solved, out_path=f"cop_vs_ambient_plr_benchmark_combined_r515b_{CASE}.csv")
    _plot_comparison(results)


if __name__ == "__main__":
    main()


