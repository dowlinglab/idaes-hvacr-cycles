import csv
import itertools
import math
from dataclasses import dataclass

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr import Mode, SimpleVaporCompressionCyclePLR


@dataclass(frozen=True)
class CandidateBand:
    name: str
    evap_sat_candidates_c: tuple[float, ...]
    cond_approach_candidates_c: tuple[float, ...]
    low_side_pressure_kpa: tuple[float, float]
    high_side_pressure_kpa: tuple[float, float]
    subcooling_candidates_c: tuple[float, ...]
    superheating_candidates_c: tuple[float, ...]


def _solve_point(cycle, kwargs):
    try:
        cycle.set_specifications(**kwargs)
        _, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        cop_part_load = cycle.get_part_load_cop() if converged else float("nan")
        if not math.isfinite(cop_part_load):
            return False, float("nan")
        return bool(converged), float(cop_part_load)
    except Exception:
        try:
            retry_kwargs = dict(kwargs)
            retry_kwargs["debug_disable_arc_pressure_eq"] = True
            cycle.set_specifications(**retry_kwargs)
            _, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
            cop_part_load = cycle.get_part_load_cop() if converged else float("nan")
            if not math.isfinite(cop_part_load):
                return False, float("nan")
            return bool(converged), float(cop_part_load)
        except Exception:
            return False, float("nan")


def _evaluate_candidate_band(cycle, band, ambient_c, common_kwargs):
    attempts = 0
    feasible = 0
    best_cop = float("nan")

    for evap_sat_c, cond_app_c, subcool_c, superheat_c in itertools.product(
        band.evap_sat_candidates_c,
        band.cond_approach_candidates_c,
        band.subcooling_candidates_c,
        band.superheating_candidates_c,
    ):
        attempts += 1
        run_kwargs = dict(common_kwargs)
        run_kwargs.update(
            dict(
                low_side_pressure=band.low_side_pressure_kpa,
                high_side_pressure=band.high_side_pressure_kpa,
                ambient_temperature=float(ambient_c),
                evap_sat_temperature=float(evap_sat_c),
                condenser_approach=float(cond_app_c),
                subcooling=float(subcool_c),
                superheating=float(superheat_c),
            )
        )
        ok, cop = _solve_point(cycle, run_kwargs)
        if ok:
            feasible += 1
            if (not math.isfinite(best_cop)) or (cop > best_cop):
                best_cop = float(cop)

    return {
        "ambient_C": float(ambient_c),
        "case": band.name,
        "attempts": int(attempts),
        "feasible": int(feasible),
        "feasible_fraction": float(feasible / attempts if attempts else 0.0),
        "best_cop": float(best_cop) if math.isfinite(best_cop) else float("nan"),
    }


def _evaluate_repo_relaxed_case(cycle, ambient_c, common_kwargs):
    # Mirrors active outside-of-demo runner settings in run_plr_cold_storage_r134a_copy.py.
    attempts = 0
    feasible = 0
    best_cop = float("nan")
    for subcool_c, superheat_c in itertools.product((3.0, 5.0, 7.0), (3.0, 5.0, 7.0)):
        attempts += 1
        run_kwargs = dict(common_kwargs)
        run_kwargs.update(
            dict(
                low_side_pressure=(60.0, 120.0),
                high_side_pressure=(500.0, 1000.0),
                evaporator_temperature=(-55.0, -20.0),
                condenser_temperature=(15.0, 50.0),
                ambient_temperature=float(ambient_c),
                condenser_approach=20.0,
                subcooling=float(subcool_c),
                superheating=float(superheat_c),
            )
        )
        ok, cop = _solve_point(cycle, run_kwargs)
        if ok:
            feasible += 1
            if (not math.isfinite(best_cop)) or (cop > best_cop):
                best_cop = float(cop)

    return {
        "ambient_C": float(ambient_c),
        "case": "Repo relaxed (current run settings)",
        "attempts": int(attempts),
        "feasible": int(feasible),
        "feasible_fraction": float(feasible / attempts if attempts else 0.0),
        "best_cop": float(best_cop) if math.isfinite(best_cop) else float("nan"),
    }


def _plot_overlay(out_png, out_pdf, title, ambient, row_a, row_b):
    xa = np.array([r["ambient_C"] for r in row_a], dtype=float)
    xb = np.array([r["ambient_C"] for r in row_b], dtype=float)

    ya_feas = np.array([r["feasible_fraction"] for r in row_a], dtype=float)
    yb_feas = np.array([r["feasible_fraction"] for r in row_b], dtype=float)

    ya_cop = np.array([r["best_cop"] for r in row_a], dtype=float)
    yb_cop = np.array([r["best_cop"] for r in row_b], dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.2), dpi=170, sharex=True)

    ax0 = axes[0]
    ax0.plot(xa, ya_feas, marker="o", linewidth=2.0, label=row_a[0]["case"])
    ax0.plot(xb, yb_feas, marker="s", linewidth=2.0, label=row_b[0]["case"])
    ax0.set_ylabel("Feasible fraction")
    ax0.set_ylim(-0.02, 1.02)
    ax0.grid(True, alpha=0.25)
    ax0.legend(loc="best")
    ax0.set_title(title)

    ax1 = axes[1]
    oka = np.isfinite(ya_cop)
    okb = np.isfinite(yb_cop)
    if np.any(oka):
        ax1.plot(xa[oka], ya_cop[oka], marker="o", linewidth=2.0, label=f"{row_a[0]['case']}: best feasible COP")
    if np.any(~oka):
        ax1.plot(xa[~oka], np.zeros(np.sum(~oka)), "x", label=f"{row_a[0]['case']}: infeasible")

    if np.any(okb):
        ax1.plot(xb[okb], yb_cop[okb], marker="s", linewidth=2.0, label=f"{row_b[0]['case']}: best feasible COP")
    if np.any(~okb):
        ax1.plot(xb[~okb], np.zeros(np.sum(~okb)), "x", label=f"{row_b[0]['case']}: infeasible")

    ax1.set_xlabel("Ambient temperature (C)")
    ax1.set_ylabel("Best feasible COP")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    ambient_grid_c = np.arange(10.0, 46.0, 5.0)

    strict_band = CandidateBand(
        name="Strict band (collaborator)",
        evap_sat_candidates_c=(-30.0, -29.0, -28.0),
        cond_approach_candidates_c=(8.0, 9.0, 10.0),
        low_side_pressure_kpa=(60.0, 120.0),
        high_side_pressure_kpa=(500.0, 1000.0),
        subcooling_candidates_c=(3.0, 5.0, 7.0),
        superheating_candidates_c=(3.0, 5.0, 7.0),
    )
    relaxed_earlier = CandidateBand(
        name="Relaxed band (earlier demo)",
        evap_sat_candidates_c=(-35.0, -30.0, -25.0),
        cond_approach_candidates_c=(8.0, 14.0, 20.0),
        low_side_pressure_kpa=(60.0, 140.0),
        high_side_pressure_kpa=(450.0, 1200.0),
        subcooling_candidates_c=(3.0, 5.0, 7.0),
        superheating_candidates_c=(3.0, 5.0, 7.0),
    )

    common_kwargs = dict(
        evaporator_temperature=(-55.0, -20.0),
        condenser_temperature=(15.0, 60.0),
        max_pressure_ratio=20.0,
        plr=0.75,
        cd=0.13,
    )

    cycle = SimpleVaporCompressionCyclePLR(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        mode=Mode.PH,
    )
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    rows = []
    strict_rows = []
    relaxed1_rows = []
    repo_rows = []

    for amb in ambient_grid_c:
        r_strict = _evaluate_candidate_band(cycle, strict_band, float(amb), common_kwargs)
        r_relaxed1 = _evaluate_candidate_band(cycle, relaxed_earlier, float(amb), common_kwargs)
        r_repo = _evaluate_repo_relaxed_case(cycle, float(amb), common_kwargs)

        strict_rows.append(r_strict)
        relaxed1_rows.append(r_relaxed1)
        repo_rows.append(r_repo)
        rows.extend([r_strict, r_relaxed1, r_repo])

    out_csv = "plr_two_relaxed_overlays_r134a.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["ambient_C", "case", "attempts", "feasible", "feasible_fraction", "best_cop"],
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    _plot_overlay(
        out_png="overlay_1_strict_vs_earlier_relaxed_r134a.png",
        out_pdf="overlay_1_strict_vs_earlier_relaxed_r134a.pdf",
        title="Overlay 1: Strict vs Earlier Relaxed Band (R134a)",
        ambient=ambient_grid_c,
        row_a=strict_rows,
        row_b=relaxed1_rows,
    )
    _plot_overlay(
        out_png="overlay_2_strict_vs_repo_relaxed_r134a.png",
        out_pdf="overlay_2_strict_vs_repo_relaxed_r134a.pdf",
        title="Overlay 2: Strict vs Repo-Current Relaxed Setup (R134a)",
        ambient=ambient_grid_c,
        row_a=strict_rows,
        row_b=repo_rows,
    )

    print(f"Saved: {out_csv}")
    print("Saved: overlay_1_strict_vs_earlier_relaxed_r134a.png/.pdf")
    print("Saved: overlay_2_strict_vs_repo_relaxed_r134a.png/.pdf")
    print(
        "Repo-relaxed settings confirmed from run_plr_cold_storage_r134a_copy.py: "
        "evaporator_temperature=(-55,-20), condenser_temperature=(15,50), "
        "condenser_approach=20, low/high pressure=(60,120)/(500,1000), plr=0.75, cd=0.13, "
        "with subcooling and superheating swept over {3,5,7}"
    )


if __name__ == "__main__":
    main()
