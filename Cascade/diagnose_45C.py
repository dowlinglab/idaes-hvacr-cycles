"""
Author: Shilpa Narasimhan

Support: Claude AI

Date Created: 09/22/2026

SINGLE-POINT DIAGNOSTIC COPY of run_cascade_sweep_r1234yf.py, for digging
into the 45C convergence issue without waiting on/rerunning the other 7
ambient points every time. The ONLY functional change from the real sweep
script is AMBIENT_SWEEP_C below (now just [45.0]); everything else --
CascadeCycle setup, solver options (imported live from
vapor_compression_cascade.py), setpoints -- is identical to the real
sweep, so whatever happens here is exactly what the real sweep would do
at 45C. Output file names are prefixed "diag45_" so this never
overwrites the real sweep's CSVs/plots.

Because AMBIENT_SWEEP_C now has only one point, that point is always
"the last point" -- so the existing is_last_point verbose flag below
means every solve here prints the full IPOPT iteration log automatically,
with no code change needed.

Ambient-temperature sweep runner for the two-loop CASCADE cycle defined
in Cascade/vapor_compression_cascade.py (CascadeCycle). Compares two
variants of the cascade -- hot=R134a/cold=CO2 vs hot=R1234yf/cold=CO2 --
across the same ambient sweep, using the same setpoints and the same
solve procedure for both, so the comparison is apples-to-apples by
construction. The cold loop (CO2) is identical in both variants; only
the hot loop's fluid changes.

Setpoints:
    PLR = 0.75
    CD = 0.13            <- kept nonzero on purpose: CD=0 assumes the part
                            -load correction is a no-op (PLF=1), i.e. a
                            perfect cycle. 0.13 keeps the part-load
                            degradation modeled instead of assumed away.
    cold loop evap_sat_temperature = -29 C   (fixed load-side target)
    hot loop condenser_approach    = 9 C     (ambient-side approach)
    cascade_approach_dT            = 3 C     (internal hot-evap/cold-cond
                                              approach inside the cascade
                                              HX itself)
    superheat/subcool = 3 C floor, 8 C cap, both loops
    compressor_efficiency = 0.80, both loops

Carnot comparison:
    COP_carnot = T_cold_K / (T_hot_K - T_cold_K), using T_cold = the cold
    loop's evap_sat_temperature (-29C, fixed) and T_hot = ambient +
    condenser_approach (the hot loop's condenser saturation temp). This
    is the reversible-cycle limit -- NOT the same as compressor_efficiency
    =1 in this model (which would still have superheat/subcool/pressure-
    drop irreversibilities); it depends only on the two true reservoir
    temperatures, so it is IDENTICAL for both fluid variants at a given
    ambient, and unaffected by the cascade's internal coupling step.
    Compared against COP_full (not COP_part): Carnot is a thermodynamic-
    cycle bound, while PLR/CD part-load derating is a separate, empirical
    equipment/controls correction unrelated to reversibility -- mixing
    the two would conflate different kinds of "loss." second_law_efficiency
    = COP_full / COP_carnot is reported per variant per point.

Fresh build per ambient point: a brand-new CascadeCycle is constructed
    and initialized independently for every point in AMBIENT_SWEEP_C, each
    one seeded via specify_initial_conditions(hot_ambient_C=<that point's
    own temperature>, ...) -- no point depends on any other point's
    solution. This matches validate_cascade.py's _build/stage3_ambient_sweep
    pattern. (Earlier versions of this runner used a single persistent
    CascadeCycle, warm-started sequentially point-to-point; that approach
    could leave the coldest points several warm-start hops removed from a
    starting guess tailored to them, occasionally failing to converge. The
    fresh-build approach gives every point its own tailored guess instead.)

What the reported COP represents:
    `optimize_COP(optimize=True)` MAXIMIZES COP subject to the active
    constraints, including the superheat/subcool floor+cap -- it does NOT
    solve at one fixed operating point. Each variant independently picks
    whichever point inside [floor, cap] is best for it at that ambient, so
    the two variants' actual superheat/subcool can legitimately differ
    point to point; that is expected optimizer behavior, not a bug. Every
    unit operation's actual inlet/outlet temperature and pressure, plus
    heat duty / work / actual superheat / actual subcool, for BOTH loops
    (hot_ and cold_ prefixed keys), is recorded for every converged point
    -- console report and CSVs -- so the underlying operating point is
    always visible alongside the COP.

    Mass flow: only ONE loop's mass flow is fixed to 1 kg/s
    (`CascadeCycle`'s `flow_fixed_role`, default "cold"); the other
    loop's mass flow is a free variable solved by the cascade energy
    balance. So the cold loop's heat_duty/work_mechanical values are
    literally per 1 kg of cold-loop refrigerant, but the hot loop's are
    per whatever mass flow the solver finds for THAT loop -- not
    necessarily 1 kg. Fine for COP (an intensive ratio either way), but
    the hot-loop duty/work numbers should not be read as normalized to
    the same basis as the cold loop's.
"""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from pyomo.environ import value

# This file lives in Cascade/, alongside vapor_compression_cascade.py --
# that import is a plain local one. R1234yf's actual parameter JSON lives
# one level up at repo root (R1234yf/r1234yf_general_helmholtz.json,
# referenced directly in _build_cycle's custom_json path), which is why
# the repo root needs to be on sys.path too.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT_DIR = os.path.dirname(_THIS_DIR)
if _REPO_ROOT_DIR not in sys.path:
    sys.path.insert(0, _REPO_ROOT_DIR)

from vapor_compression_cascade import CascadeCycle


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

# DIAGNOSTIC: just the one problem point. (Real sweep script uses
# np.arange(10, 46, 5) -- 10 through 45 in 5C steps.)
AMBIENT_SWEEP_C = np.array([45.0])

def _build_cycle(fluid_name):
    """The only fluid-specific branch in this whole runner. `fluid_name`
    selects the HOT loop's fluid ("R134a" or "R1234yf"); the cold loop is
    always CO2. Returns a full two-loop CascadeCycle, not a single-stage
    cycle -- everything downstream (_run_fluid, _report_unit_ops) has to
    match CascadeCycle's real API, not the old single-stage one."""
    cold_fluid = {
        "gh_component": "co2", "coolprop_name": "CO2", "custom_json": None,
        "efficiency": COMPRESSOR_EFFICIENCY,
        "superheat_max": SUPERHEAT_MAX_C, "subcool_max": SUBCOOL_MAX_C,
        "Tmax_C": None,
    }

    if fluid_name == "R1234yf":
        hot_fluid = {
            "gh_component": "r1234yf", "coolprop_name": "R1234yf",
            "custom_json": os.path.join(_REPO_ROOT_DIR, "R1234yf", "r1234yf_general_helmholtz.json"),
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
    """Snapshot ONE loop's (role="hot" or "cold") actual operating
    condition after a solve -- inlet/outlet T & P, duty/work, and actual
    superheat/subcool. Same fields _report_unit_ops always returned, just
    scoped to one loop's four units instead of a single unnamed loop
    (the cascade has two loops, named fs.hot_evaporator/fs.cold_evaporator
    etc., not fs.evaporator)."""
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
    """Both loops' unit-op snapshot, merged into one flat dict with
    hot_/cold_ prefixes (e.g. "cold_evap_superheat_actual_K",
    "hot_cond_subcool_actual_K") -- suitable for a single CSV row."""
    fs = cycle.model.fs
    out = {}
    for role in ("hot", "cold"):
        for k, v in _report_loop_ops(fs, role).items():
            out[f"{role}_{k}"] = v
    return out


def _carnot_cop(ambient_c):
    """Reversible-cycle COP limit: T_cold / (T_hot - T_cold), both in K.

    T_cold = evap_sat_temperature (fixed, -29C). T_hot = condenser
    saturation temp = ambient + condenser_approach. Both are
    equality-constrained targets, independent of which fluid is running --
    so this number is IDENTICAL for both fluids at a given ambient. This
    is NOT compressor_efficiency=1 in our model (that would still have
    superheat/subcool/pressure-drop irreversibilities); it's the pure
    thermodynamic ceiling, unrelated to any of our modeling assumptions.
    """
    t_cold_k = EVAP_SAT_TARGET_C + C_to_K
    t_hot_k = ambient_c + COND_APPROACH_C + C_to_K
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


def _run_fluid(fluid_name, out_stem):
    # Cold loop's spec dict: fixed load-side target (evap_sat_temperature),
    # no ambient_temperature -- that belongs to the hot loop only.
    cold_kwargs = dict(
        evap_sat_temperature=EVAP_SAT_TARGET_C, superheating=SUPERHEAT_C,
        superheating_max=SUPERHEAT_MAX_C, subcooling=SUBCOOL_C,
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
    for ambient_c in AMBIENT_SWEEP_C:
        # Fresh CascadeCycle for every point, seeded at THAT point's own
        # ambient temperature -- no dependency on any other point's
        # solution (matches validate_cascade.py's _build pattern).
        cycle, fluid_used = _build_cycle(fluid_name)
        cycle.specify_initial_conditions(hot_ambient_C=float(ambient_c), cold_evap_C=EVAP_SAT_TARGET_C)
        cycle.initialize(verbose=False)

        hot_kwargs = dict(hot_kwargs_base)
        hot_kwargs["ambient_temperature"] = float(ambient_c)
        cycle.set_specifications(hot=hot_kwargs, cold=cold_kwargs,
                                  cascade_approach_dT=CASCADE_APPROACH_DT_C,
                                  plr=PLR_VALUE, cd=CD_VALUE)

        # Temporary diagnostic: turn on solver-verbose output for the last
        # (hottest) point in the sweep, and print any swallowed exception's
        # actual message -- added to find out why R1234yf failed at 45C
        # while R134a converged fine there (see BREADCRUMB.md). Safe to
        # remove once that's diagnosed.
        is_last_point = (ambient_c == AMBIENT_SWEEP_C[-1])
        try:
            cop_full, ok = cycle.optimize_COP(verbose=is_last_point, initialize=True, optimize=True)
            if is_last_point and not ok:
                print(f"  -> {fluid_used} at {ambient_c}C: optimize_COP returned ok=False (no exception raised)")
        except Exception as exc:
            print(f"  -> {fluid_used} at {ambient_c}C raised exception: {exc!r}")
            cop_full, ok = float("nan"), False

        cop_part = cycle.get_part_load_cop() if ok else float("nan")
        cop_carnot = _carnot_cop(ambient_c)
        second_law_eff = (cop_full / cop_carnot) if ok else float("nan")
        unit_ops = _report_unit_ops(cycle) if ok else {k: float("nan") for k in _UNIT_OP_FIELDS}
        solved[ambient_c] = (ok, cop_full if ok else float("nan"), cop_part, cop_carnot, second_law_eff, unit_ops)

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

    ambient_sorted = sorted(solved.keys())
    cop_full_vals = [solved[a][1] for a in ambient_sorted]
    cop_part_vals = [solved[a][2] for a in ambient_sorted]
    converged_vals = [solved[a][0] for a in ambient_sorted]

    out_csv = f"{out_stem}.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ambient_C", "cop_full", "cop_part", "cop_carnot", "second_law_efficiency",
                     "converged", "fluid_name_used"] + _UNIT_OP_FIELDS)
        for ta in ambient_sorted:
            ok, cfull, cpart, ccarnot, seff, unit_ops = solved[ta]
            w.writerow([ta, cfull, cpart, ccarnot, seff, int(ok), fluid_used]
                       + [unit_ops[k] for k in _UNIT_OP_FIELDS])

    print(f"Saved: {out_csv}")
    print(f"Converged points: {sum(converged_vals)}/{len(ambient_sorted)}")
    return np.array(ambient_sorted), cop_full_vals, cop_part_vals, converged_vals, solved


def _print_full_report(fluid_solved, tol_K=0.5):
    """Print COP and every unit operation's actual operating condition for
    EVERY ambient point, every fluid in `fluid_solved` (a {label: solved}
    dict, N fluids -- not hardcoded to two).

    NOTE: each fluid independently solves "maximize COP subject to the
    [floor, cap] superheat/subcool window" -- so different fluids landing on
    DIFFERENT actual superheat/subcool at the same ambient is an EXPECTED
    outcome of the optimization, not an error. The SH/SC delta (vs. the
    first fluid listed, used only as a reporting baseline) is reported as a
    fact, not flagged as a mismatch/defect.
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
    """One row per ambient point, every fluid in `fluid_solved` side by side
    -- COP + full operating conditions for each fluid's own independently-
    optimized solve, plus the SH/SC delta of each non-baseline fluid vs. the
    first-listed (baseline) fluid, as a plain fact (each fluid is free to
    land anywhere in [3,8]K; differing is expected, not an error -- so this
    is informational, not a pass/fail flag).
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

    ax.set_title("DIAGNOSTIC: 45C only -- R134a vs R1234yf vs Carnot")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig("diag45_comparison.png", bbox_inches="tight")
    fig.savefig("diag45_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved: diag45_comparison.png")
    print("Saved: diag45_comparison.pdf")


def main():
    results = {}
    results["R134a"] = _run_fluid("R134a", "diag45_cop_vs_ambient_plr_benchmark_r134a")
    results["R1234yf"] = _run_fluid("R1234yf", "diag45_cop_vs_ambient_plr_benchmark_r1234yf")

    fluid_solved = {label: r[4] for label, r in results.items()}
    _print_full_report(fluid_solved)
    _write_combined_csv(fluid_solved, out_path="diag45_combined.csv")
    _plot_comparison(results)


if __name__ == "__main__":
    main()


