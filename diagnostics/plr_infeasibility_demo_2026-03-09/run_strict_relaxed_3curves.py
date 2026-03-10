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
class Band:
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
        return True, float(cop_part_load)
    except Exception:
        try:
            retry_kwargs = dict(kwargs)
            retry_kwargs["debug_disable_arc_pressure_eq"] = True
            cycle.set_specifications(**retry_kwargs)
            _, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
            cop_part_load = cycle.get_part_load_cop() if converged else float("nan")
            if not (converged and math.isfinite(cop_part_load)):
                return False, float("nan")
            return True, float(cop_part_load)
        except Exception:
            return False, float("nan")


def _evaluate_band_for_scsh(cycle, band, ambient_c, scsh_c, common_kwargs):
    attempts = 0
    feasible = 0
    best_cop = float("nan")

    for evap_sat_c, cond_app_c in itertools.product(
        band.evap_sat_candidates_c,
        band.cond_approach_candidates_c,
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
                subcooling=float(scsh_c),
                superheating=float(scsh_c),
            )
        )
        ok, cop = _solve_point(cycle, run_kwargs)
        if ok:
            feasible += 1
            if (not math.isfinite(best_cop)) or (cop > best_cop):
                best_cop = float(cop)

    return {
        "ambient_C": float(ambient_c),
        "band": band.name,
        "scsh_C": float(scsh_c),
        "attempts": int(attempts),
        "feasible": int(feasible),
        "feasible_fraction": float(feasible / attempts if attempts else 0.0),
        "best_cop": float(best_cop) if math.isfinite(best_cop) else float("nan"),
    }


def main():
    ambient_grid_c = np.arange(10.0, 46.0, 5.0)
    scsh_values_c = (3.0, 5.0, 7.0)

    strict = Band(
        name="Strict",
        evap_sat_candidates_c=(-30.0, -29.0, -28.0),
        cond_approach_candidates_c=(8.0, 9.0, 10.0),
        low_side_pressure_kpa=(60.0, 120.0),
        high_side_pressure_kpa=(500.0, 1000.0),
    )
    relaxed = Band(
        name="Relaxed",
        evap_sat_candidates_c=(-35.0, -30.0, -25.0),
        cond_approach_candidates_c=(8.0, 14.0, 20.0),
        low_side_pressure_kpa=(60.0, 140.0),
        high_side_pressure_kpa=(450.0, 1200.0),
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
    for amb in ambient_grid_c:
        for scsh in scsh_values_c:
            rows.append(_evaluate_band_for_scsh(cycle, strict, float(amb), float(scsh), common_kwargs))
            rows.append(_evaluate_band_for_scsh(cycle, relaxed, float(amb), float(scsh), common_kwargs))

    out_csv = "plr_strict_relaxed_scsh_3_5_7_r134a.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["ambient_C", "band", "scsh_C", "attempts", "feasible", "feasible_fraction", "best_cop"],
        )
        w.writeheader()
        for row in rows:
            w.writerow(row)

    out_png = "plr_strict_relaxed_scsh_3_5_7_r134a.png"
    out_pdf = "plr_strict_relaxed_scsh_3_5_7_r134a.pdf"

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.4), dpi=170, sharex=True)

    colors = {3.0: "#1f77b4", 5.0: "#ff7f0e", 7.0: "#2ca02c"}

    for idx, band_name in enumerate(("Strict", "Relaxed")):
        ax = axes[idx]
        band_rows = [r for r in rows if r["band"] == band_name]
        for scsh in scsh_values_c:
            cur = [r for r in band_rows if abs(r["scsh_C"] - scsh) < 1e-9]
            x = np.array([r["ambient_C"] for r in cur], dtype=float)
            y = np.array([r["best_cop"] for r in cur], dtype=float)
            ok = np.isfinite(y)
            label = f"subcool=superheat={int(scsh)} C"
            if np.any(ok):
                ax.plot(x[ok], y[ok], marker="o", linewidth=2.0, color=colors[scsh], label=label)
            if np.any(~ok):
                ax.plot(x[~ok], np.zeros(np.sum(~ok)), "x", color=colors[scsh])
        ax.set_ylabel("Best feasible COP")
        ax.grid(True, alpha=0.25)
        ax.set_title(f"{band_name} band")
        ax.legend(loc="best")

    axes[-1].set_xlabel("Ambient temperature (C)")
    fig.suptitle("R134a: Strict vs Relaxed bands with sc=sh in {3,5,7} C", y=0.995)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
