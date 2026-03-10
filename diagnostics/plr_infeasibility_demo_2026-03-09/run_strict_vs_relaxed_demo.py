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
class SweepCase:
    name: str
    evap_sat_candidates_c: tuple[float, ...]
    cond_approach_candidates_c: tuple[float, ...]
    low_side_pressure_kpa: tuple[float, float]
    high_side_pressure_kpa: tuple[float, float]


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


def _evaluate_case(cycle, case, ambient_c, common_kwargs):
    attempts = 0
    feasible = 0
    best_cop = float("nan")
    best_evap = float("nan")
    best_cond_app = float("nan")

    for evap_sat_c, cond_app_c in itertools.product(case.evap_sat_candidates_c, case.cond_approach_candidates_c):
        attempts += 1
        run_kwargs = dict(common_kwargs)
        run_kwargs.update(
            dict(
                low_side_pressure=case.low_side_pressure_kpa,
                high_side_pressure=case.high_side_pressure_kpa,
                ambient_temperature=float(ambient_c),
                evap_sat_temperature=float(evap_sat_c),
                condenser_approach=float(cond_app_c),
            )
        )
        ok, cop = _solve_point(cycle, run_kwargs)
        if ok:
            feasible += 1
            if (not math.isfinite(best_cop)) or (cop > best_cop):
                best_cop = float(cop)
                best_evap = float(evap_sat_c)
                best_cond_app = float(cond_app_c)

    frac_feasible = feasible / attempts if attempts else 0.0
    return {
        "ambient_C": float(ambient_c),
        "case": case.name,
        "attempts": int(attempts),
        "feasible": int(feasible),
        "feasible_fraction": float(frac_feasible),
        "best_cop": float(best_cop) if math.isfinite(best_cop) else float("nan"),
        "best_evap_sat_C": best_evap,
        "best_cond_approach_C": best_cond_app,
    }


def main():
    tsp_c = -20.0
    ambient_grid_c = np.arange(10.0, 46.0, 5.0)

    strict_case = SweepCase(
        name="Strict band",
        evap_sat_candidates_c=(-30.0, -29.0, -28.0),
        cond_approach_candidates_c=(8.0, 9.0, 10.0),
        low_side_pressure_kpa=(60.0, 120.0),
        high_side_pressure_kpa=(500.0, 1000.0),
    )

    cycle = SimpleVaporCompressionCyclePLR(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        mode=Mode.PH,
    )
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    common_kwargs = dict(
        evaporator_temperature=(-55.0, -20.0),
        condenser_temperature=(15.0, 60.0),
        subcooling=1.0,
        superheating=3.0,
        max_pressure_ratio=20.0,
        plr=0.75,
        cd=0.13,
    )

    rows = []
    for ambient_c in ambient_grid_c:
        rows.append(_evaluate_case(cycle, strict_case, float(ambient_c), common_kwargs))

    out_csv = "plr_strict_band_feasibility_r134a.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "ambient_C",
                "case",
                "attempts",
                "feasible",
                "feasible_fraction",
                "best_cop",
                "best_evap_sat_C",
                "best_cond_approach_C",
            ],
        )
        w.writeheader()
        for row in rows:
            w.writerow(row)

    strict_rows = [r for r in rows if r["case"] == strict_case.name]

    x_strict = np.array([r["ambient_C"] for r in strict_rows], dtype=float)
    y_strict_feas = np.array([r["feasible_fraction"] for r in strict_rows], dtype=float)
    y_strict_cop = np.array([r["best_cop"] for r in strict_rows], dtype=float)

    out_png = "plr_strict_band_feasibility_r134a.png"
    out_pdf = "plr_strict_band_feasibility_r134a.pdf"

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.5), dpi=170, sharex=True)

    ax0 = axes[0]
    ax0.plot(x_strict, y_strict_feas, marker="o", linewidth=2.0, label="Strict band")
    ax0.set_ylabel("Feasible fraction")
    ax0.set_ylim(-0.02, 1.02)
    ax0.grid(True, alpha=0.25)
    ax0.legend(loc="best")
    ax0.set_title("R134a PLR Feasibility Demonstration (T_sp = -20 C)")

    ax1 = axes[1]
    strict_ok = np.isfinite(y_strict_cop)
    if np.any(strict_ok):
        ax1.plot(x_strict[strict_ok], y_strict_cop[strict_ok], marker="o", linewidth=2.0, label="Strict band: best feasible COP")
    if np.any(~strict_ok):
        ax1.plot(x_strict[~strict_ok], np.zeros(np.sum(~strict_ok)), "x", label="Strict infeasible")

    ax1.set_xlabel("Ambient temperature (C)")
    ax1.set_ylabel("Best feasible COP")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print(f"Ambient grid (C): {ambient_grid_c.tolist()}")
    print(
        "Strict feasible points by ambient:",
        [f"{int(r['feasible'])}/{int(r['attempts'])}" for r in strict_rows],
    )


if __name__ == "__main__":
    main()
