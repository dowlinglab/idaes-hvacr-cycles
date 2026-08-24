"""
Author: Shilpa Narasimhan

Support: Claude AI

Date Created: 08/24/2026

WORK IN PROGRESS -- intended end state: an ambient-temperature sweep
runner for the two-loop CASCADE cycle in
Cascade/vapor_compression_cascade.py (CascadeCycle, hot=R134a/cold=CO2 by
default), NOT a single-stage fluid benchmark.

This file started life as a copy of R1234yf/'s single-stage R134a vs
R1234yf vs R515B PLR benchmark runner (see BREADCRUMB.md's 2026-08-24
entry for the full history). Two cleanup passes so far:
  1. R515B has been removed from this file entirely (import, sys.path
     setup, the `_build_cycle` branch, the sweep call in `main()`, and
     the plot title/legend). R515B is a fixed-composition blend with its
     own hand-built property package -- it doesn't apply to a
     general_helmholtz cascade run, and PLR-benchmarking multiple
     standalone fluids side-by-side isn't this file's job anymore.
  2. Docstring corrected to describe what this file actually is (this
     block).

STILL NOT DONE (see BREADCRUMB.md for the full section-by-section plan):
`_build_cycle`, `_run_fluid`, `_report_unit_ops`, and the
`set_specifications`/`specify_initial_conditions` calls in the sweep loop
still reflect the OLD single-stage API (`SimpleVaporCompressionCyclePLRFixed`,
flat kwargs, one fluid at a time) -- none of that matches `CascadeCycle`'s
real API yet (`hot=`/`cold=` spec dicts, `hot_ambient_C=`/`cold_evap_C=`,
etc.). As of this note, running this file will still fail on the R134a
branch (`SimpleVaporCompressionCyclePLRFixed`/`Mode` are referenced but
never imported here). The `CascadeCycle` import below is in place for
that conversion but not yet wired into the sweep logic.

Setpoints (collaborator-sourced, unchanged from prior benchmarks):
    PLR = 0.75
    CD = 0.13            <- kept nonzero on purpose: CD=0 assumes the part
                            -load correction is a no-op (PLF=1), i.e. a
                            perfect cycle. 0.13 keeps the part-load
                            degradation modeled instead of assumed away.
    evap_sat_temperature = -29 C
    condenser_approach   = 9 C
    superheat/subcool     = 3 C floor, 8 C cap (see BREADCRUMB.md for the
                            floor+cap fix history)
    compressor_efficiency = 0.80

Carnot comparison:
    COP_carnot = T_cold_K / (T_hot_K - T_cold_K), using T_cold=evap_sat
    (-29C, fixed) and T_hot=ambient+condenser_approach. This is the
    reversible-cycle limit -- NOT the same as compressor_efficiency=1 in
    this model (which would still have superheat/subcool/pressure-drop
    irreversibilities); it depends only on the two saturation
    temperatures, so it is IDENTICAL for both fluids at a given ambient.
    Compared against COP_full (not COP_part): Carnot is a thermodynamic-
    cycle bound, while PLR/CD part-load derating is a separate, empirical
    equipment/controls correction unrelated to reversibility -- mixing
    the two would conflate different kinds of "loss." second_law_efficiency
    = COP_full / COP_carnot is reported per fluid per point.

Sequential warm-start (the fix for the earlier 10C/15C non-convergence):
    Diagnosed directly -- at 10C/15C ambient, the ONLY constraint with a
    large residual was the subcool floor (actual subcool stuck at exactly
    0K while everything else -- superheat, pressure ratio, COP -- looked
    physically normal). `initialize()` always cold-starts from a fixed
    high_side_temperature=30C guess regardless of ambient; that's close
    enough to converge near ambient=20-25C but IPOPT (a local solver) got
    stuck at the extremes. This is a warm-start problem, not a physical
    infeasibility: matching the initial guess to each ambient's own target
    state converges cleanly, and so does marching through the sweep in
    order and reusing each converged solution as the next point's warm
    start (what this runner does). No bound, cap, or CD change ever fixed
    this -- all were tested and ruled out first.

IMPORTANT -- what the reported COP actually represents:
    `optimize_COP(optimize=True)` MAXIMIZES COP subject to the active
    constraints, including the superheat/subcool floor+cap. It does NOT
    solve at one fixed operating point per fluid -- it picks whichever
    point INSIDE the [floor, cap] window is best for THAT fluid at THAT
    ambient. This means the two fluids' actual superheat/subcool CAN and
    SHOULD legitimately differ point to point -- each fluid is finding its
    own best operating condition independently, not being forced onto a
    shared one. That is expected optimizer behavior, not a bug (confirmed:
    at ambient=10C, R134a's COP-maximizing point is superheat=3K/subcool=8K
    while R1234yf's is superheat=8K/subcool=8K -- both are each fluid's
    genuine optimum, not an error).

    Because of this, a bare "COP=X" per point does not tell you what
    condition that COP was achieved at, or whether the two fluids were
    compared at the same condition. So every unit operation's actual
    inlet/outlet temperature and pressure, plus heat duty / work / actual
    superheat / actual subcool, is recorded and reported for every
    converged point, for both fluids -- both in the console report and in
    the CSVs -- so the underlying operating point is always visible
    alongside the COP, not hidden behind it. Mass flow is fixed to 1 kg/s
    (a normalized basis -- see `specify_initial_conditions`/`initialize`
    in vapor_compression_plr_fixed.py), so heat_duty/work_mechanical here
    are PER UNIT MASS FLOW, not absolute capacities -- fine for COP (an
    intensive ratio) but not to be read as real duty magnitudes.
"""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from pyomo.environ import value

# This file now lives IN R1234yf/, so the R1234yf cycle module is a plain
# local import. The R134a-side module lives one level up at repo root, so
# that's the one that needs the sys.path insert (reversed from when this
# file lived at repo root).
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT_DIR = os.path.dirname(_THIS_DIR)
if _REPO_ROOT_DIR not in sys.path:
    sys.path.insert(0, _REPO_ROOT_DIR)

from vapor_compression_cascade import CascadeCycle
from vapor_compression_plr_r1234yf import (
    Mode as ModeR1234yf,
    SimpleVaporCompressionCyclePLRR1234yf,
)


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

# Sequential sweep order: cold-start at 20C (known to converge reliably
# from the default initial guess), then march outward in both directions,
# warm-starting each point from the previous point's converged solution.
SWEEP_ORDER_C = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 15.0, 10.0]


def _build_cycle(fluid_name):
    """The only fluid-specific branch in this whole runner."""
    if fluid_name == "R1234yf":
        cycle = SimpleVaporCompressionCyclePLRR1234yf(
            compressor_efficiency=COMPRESSOR_EFFICIENCY,
            PLR=PLR_VALUE,
            CD=CD_VALUE,
            mode=ModeR1234yf.PH,
        )
        return cycle, "R1234yf"

    cycle = SimpleVaporCompressionCyclePLRFixed(
        fluid_name=fluid_name,
        compressor_efficiency=COMPRESSOR_EFFICIENCY,
        PLR=PLR_VALUE,
        CD=CD_VALUE,
        mode=Mode.PH,
    )
    return cycle, fluid_name


C_to_K = 273.15


def _report_unit_ops(cycle):
    """Snapshot every unit operation's actual operating condition after a
    solve -- inlet/outlet T & P, duty/work, and actual superheat/subcool.
    Returns a flat dict suitable for a CSV row. Mass flow is fixed to
    1 kg/s (normalized basis), so heat_duty/work_mechanical are per unit
    mass flow, not absolute capacities.
    """
    m = cycle.model.fs

    # NOTE: in PH mode, ports only expose the actual state variables
    # (pressure, enth_mass, flow_mass) -- NOT temperature. Temperature is a
    # derived quantity, only available via control_volume.properties_in/out.
    def T_in(unit):
        return value(unit.control_volume.properties_in[0].temperature) - C_to_K

    def T_out(unit):
        return value(unit.control_volume.properties_out[0].temperature) - C_to_K

    def P(port):
        return value(port.pressure[0]) / 1000.0  # kPa

    evap_T_out = value(m.evaporator.control_volume.properties_out[0].temperature)
    evap_T_sat = value(m.evaporator.control_volume.properties_out[0].temperature_sat)
    cond_T_out = value(m.condenser.control_volume.properties_out[0].temperature)
    cond_T_sat = value(m.condenser.control_volume.properties_out[0].temperature_sat)

    return {
        "evap_in_T_C": T_in(m.evaporator), "evap_in_P_kPa": P(m.evaporator.inlet),
        "evap_out_T_C": T_out(m.evaporator), "evap_out_P_kPa": P(m.evaporator.outlet),
        "evap_T_sat_C": evap_T_sat - C_to_K,
        "evap_superheat_actual_K": evap_T_out - evap_T_sat,
        "evap_heat_duty_per_kg": value(m.evaporator.heat_duty[0]),
        "comp_in_T_C": T_in(m.compressor), "comp_in_P_kPa": P(m.compressor.inlet),
        "comp_out_T_C": T_out(m.compressor), "comp_out_P_kPa": P(m.compressor.outlet),
        "comp_ratioP": value(m.compressor.ratioP[0]),
        "comp_work_per_kg": value(m.compressor.work_mechanical[0]),
        "cond_in_T_C": T_in(m.condenser), "cond_in_P_kPa": P(m.condenser.inlet),
        "cond_out_T_C": T_out(m.condenser), "cond_out_P_kPa": P(m.condenser.outlet),
        "cond_T_sat_C": cond_T_sat - C_to_K,
        "cond_subcool_actual_K": cond_T_sat - cond_T_out,
        "cond_heat_duty_per_kg": value(m.condenser.heat_duty[0]),
        "valve_in_T_C": T_in(m.expansion_valve), "valve_in_P_kPa": P(m.expansion_valve.inlet),
        "valve_out_T_C": T_out(m.expansion_valve), "valve_out_P_kPa": P(m.expansion_valve.outlet),
    }


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


_UNIT_OP_FIELDS = [
    "evap_in_T_C", "evap_in_P_kPa", "evap_out_T_C", "evap_out_P_kPa", "evap_T_sat_C",
    "evap_superheat_actual_K", "evap_heat_duty_per_kg",
    "comp_in_T_C", "comp_in_P_kPa", "comp_out_T_C", "comp_out_P_kPa", "comp_ratioP", "comp_work_per_kg",
    "cond_in_T_C", "cond_in_P_kPa", "cond_out_T_C", "cond_out_P_kPa", "cond_T_sat_C",
    "cond_subcool_actual_K", "cond_heat_duty_per_kg",
    "valve_in_T_C", "valve_in_P_kPa", "valve_out_T_C", "valve_out_P_kPa",
]


def _run_fluid(fluid_name, out_stem):
    cycle, fluid_used = _build_cycle(fluid_name)

    high_seed = 30.0
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=high_seed)
    cycle.initialize(verbose=False)

    base_kwargs = dict(
        low_side_pressure=(60.0, 200.0),
        high_side_pressure=(500.0, 4000.0),
        evaporator_temperature=(-55.0, 0.0),
        condenser_temperature=(15.0, 60.0),
        evap_sat_temperature=EVAP_SAT_TARGET_C,
        condenser_approach=COND_APPROACH_C,
        superheating=SUPERHEAT_C,
        subcooling=SUBCOOL_C,
        superheating_max=SUPERHEAT_MAX_C,
        subcooling_max=SUBCOOL_MAX_C,
        max_pressure_ratio=MAX_PRESSURE_RATIO,
        plr=PLR_VALUE,
        cd=CD_VALUE,
    )

    solved = {}
    first_point = True
    for ambient_c in SWEEP_ORDER_C:
        run_kwargs = dict(base_kwargs)
        run_kwargs["ambient_temperature"] = float(ambient_c)
        cycle.set_specifications(**run_kwargs)

        # Only the FIRST point in the sweep re-initializes from scratch.
        # Every later point warm-starts from the previous point's
        # converged solution (this is the fix for the 10C/15C failures).
        try:
            cop_full, ok = cycle.optimize_COP(verbose=False, initialize=first_point, optimize=True)
        except Exception:
            cop_full, ok = float("nan"), False
        first_point = False

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
                f"SH_actual={unit_ops['evap_superheat_actual_K']:.2f}K | "
                f"SC_actual={unit_ops['cond_subcool_actual_K']:.2f}K | "
                f"ratioP={unit_ops['comp_ratioP']:.2f}"
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
                f"evap[in {uo['evap_in_T_C']:.2f}C/{uo['evap_in_P_kPa']:.1f}kPa -> "
                f"out {uo['evap_out_T_C']:.2f}C/{uo['evap_out_P_kPa']:.1f}kPa, "
                f"Tsat={uo['evap_T_sat_C']:.2f}C, SH={uo['evap_superheat_actual_K']:.2f}K] | "
                f"comp[in {uo['comp_in_T_C']:.2f}C/{uo['comp_in_P_kPa']:.1f}kPa -> "
                f"out {uo['comp_out_T_C']:.2f}C/{uo['comp_out_P_kPa']:.1f}kPa, "
                f"ratioP={uo['comp_ratioP']:.2f}] | "
                f"cond[in {uo['cond_in_T_C']:.2f}C/{uo['cond_in_P_kPa']:.1f}kPa -> "
                f"out {uo['cond_out_T_C']:.2f}C/{uo['cond_out_P_kPa']:.1f}kPa, "
                f"Tsat={uo['cond_T_sat_C']:.2f}C, SC={uo['cond_subcool_actual_K']:.2f}K] | "
                f"valve[in {uo['valve_in_T_C']:.2f}C/{uo['valve_in_P_kPa']:.1f}kPa -> "
                f"out {uo['valve_out_T_C']:.2f}C/{uo['valve_out_P_kPa']:.1f}kPa]"
            )

        ok_base, _, _, _, _, uo_base = rows[baseline_label]
        if ok_base and len(labels) > 1:
            for label in labels[1:]:
                ok_l, _, _, _, _, uo_l = rows[label]
                if not ok_l:
                    print(f"  -> {label} vs {baseline_label}: comparison not possible ({label} did not converge)")
                    continue
                dsh = uo_l["evap_superheat_actual_K"] - uo_base["evap_superheat_actual_K"]
                dsc = uo_l["cond_subcool_actual_K"] - uo_base["cond_subcool_actual_K"]
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
                    dsh = uo_l["evap_superheat_actual_K"] - uo_base["evap_superheat_actual_K"]
                    dsc = uo_l["cond_subcool_actual_K"] - uo_base["cond_subcool_actual_K"]
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

    ax.set_title("PLR benchmark (sequential warm-start): R134a vs R1234yf vs Carnot")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig("cop_vs_ambient_plr_benchmark_comparison.png", bbox_inches="tight")
    fig.savefig("cop_vs_ambient_plr_benchmark_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved: cop_vs_ambient_plr_benchmark_comparison.png")
    print("Saved: cop_vs_ambient_plr_benchmark_comparison.pdf")


def main():
    results = {}
    results["R134a"] = _run_fluid("R134a", "cop_vs_ambient_plr_benchmark_r134a")
    results["R1234yf"] = _run_fluid("R1234yf", "cop_vs_ambient_plr_benchmark_r1234yf")

    fluid_solved = {label: r[4] for label, r in results.items()}
    _print_full_report(fluid_solved)
    _write_combined_csv(fluid_solved)
    _plot_comparison(results)


if __name__ == "__main__":
    main()


