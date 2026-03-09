import csv
import math

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression import SimpleVaporCompressionCycle, Mode as ModeBase
from vapor_compression_plr import SimpleVaporCompressionCyclePLR, Mode as ModePLR


def _build_cycles_with_fallback():
    candidates = [
        ("r1234ze", "R1234ze(E)"),
        ("R1234ze(E)", None),
        ("R1234ZEE", None),
        ("R1234zeE", None),
        ("R1234zee", None),
    ]
    last_err = None
    for idaes_name, cp_name in candidates:
        try:
            base = SimpleVaporCompressionCycle(
                fluid_name=idaes_name,
                compressor_efficiency=0.75,
                mode=ModeBase.PH,
            )
            plr = SimpleVaporCompressionCyclePLR(
                fluid_name=idaes_name,
                compressor_efficiency=0.75,
                mode=ModePLR.PH,
            )
            if cp_name is not None:
                base.fluid_name = cp_name
                plr.fluid_name = cp_name
            return base, plr, idaes_name
        except Exception as exc:
            last_err = exc
    raise RuntimeError(f"Failed to build R1234ze cycles: {last_err}")


def _solve_with_retry(cycle, run_kwargs):
    try:
        cop_full, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        return cop_full, converged
    except Exception:
        try:
            retry_kwargs = dict(run_kwargs)
            retry_kwargs.update(dict(debug_disable_arc_pressure_eq=True))
            cycle.set_specifications(**retry_kwargs)
            cop_full, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
            return cop_full, converged
        except Exception:
            return float("nan"), False


def main():
    cold_storage_setpoint_c = -20.0
    ambient_temps = np.arange(10.0, 46.0, 5.0)
    condenser_approach_c = 20.0

    low_side_bounds = (40.0, 80.0)
    high_side_bounds = (350.0, 800.0)

    cycle_base, cycle_plr, fluid_name_used = _build_cycles_with_fallback()

    cycle_base.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle_base.initialize(verbose=False)
    cycle_plr.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle_plr.initialize(verbose=False)

    base_kwargs = dict(
        low_side_pressure=low_side_bounds,
        high_side_pressure=high_side_bounds,
        evaporator_temperature=(-55, -20),
        condenser_temperature=(15, 50),
        subcooling=1,
        superheating=3,
        max_pressure_ratio=20,
    )

    plr_kwargs = dict(base_kwargs)
    plr_kwargs.update(dict(plr=0.75, cd=0.13))

    cop_base_vals = []
    cop_plr_vals = []
    cop_carnot_vals = []
    converged_base = []
    converged_plr = []

    for ambient_c in ambient_temps:
        TL_K = cold_storage_setpoint_c + 273.15
        TH_K = ambient_c + condenser_approach_c + 273.15
        cop_carnot_vals.append(TL_K / (TH_K - TL_K))

        kwargs_base = dict(base_kwargs)
        kwargs_base.update(dict(
            ambient_temperature=float(ambient_c),
            condenser_approach=float(condenser_approach_c),
        ))
        cycle_base.set_specifications(**kwargs_base)
        cop_base, ok_base = _solve_with_retry(cycle_base, kwargs_base)

        kwargs_plr = dict(plr_kwargs)
        kwargs_plr.update(dict(
            ambient_temperature=float(ambient_c),
            condenser_approach=float(condenser_approach_c),
        ))
        cycle_plr.set_specifications(**kwargs_plr)
        _, ok_plr = _solve_with_retry(cycle_plr, kwargs_plr)
        cop_plr = cycle_plr.get_part_load_cop() if ok_plr else float("nan")

        cop_base_vals.append(float(cop_base) if (ok_base and math.isfinite(cop_base)) else float("nan"))
        cop_plr_vals.append(float(cop_plr) if (ok_plr and math.isfinite(cop_plr)) else float("nan"))
        converged_base.append(bool(ok_base))
        converged_plr.append(bool(ok_plr))

    out_csv = "cop_vs_ambient_plr_cold_storage_r1234zee_copy.csv"
    out_png = "cop_vs_ambient_plr_cold_storage_r1234zee_copy.png"
    out_pdf = "cop_vs_ambient_plr_cold_storage_r1234zee_copy.pdf"

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "ambient_C", "cop_vapor_compression", "cop_plr_adjusted", "cop_carnot_ref",
            "converged_vapor_compression", "converged_plr", "fluid_name_used",
        ])
        for Ta, c0, c1, cc, ok0, ok1 in zip(ambient_temps, cop_base_vals, cop_plr_vals, cop_carnot_vals, converged_base, converged_plr):
            w.writerow([float(Ta), c0, c1, cc, int(ok0), int(ok1), fluid_name_used])

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    ok0 = np.array(converged_base, dtype=bool) & np.isfinite(cop_base_vals)
    ok1 = np.array(converged_plr, dtype=bool) & np.isfinite(cop_plr_vals)

    if np.any(ok0):
        ax.plot(ambient_temps[ok0], np.array(cop_base_vals)[ok0], marker="o", linewidth=2.2, label="vapor_compression COP")
    if np.any(ok1):
        ax.plot(ambient_temps[ok1], np.array(cop_plr_vals)[ok1], marker="s", linewidth=2.2, label="vapor_compression + PLR COP")
    ax.plot(ambient_temps, cop_carnot_vals, linestyle="--", linewidth=2.0, color="black", label="Carnot COP")

    ax.set_title("R1234ze: COP Overlay (Same Run Grid)")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"fluid_name_used: {fluid_name_used}")
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print(f"Pressure bounds (kPa): low={low_side_bounds}, high={high_side_bounds}")
    print(f"Converged vapor_compression: {int(np.sum(ok0))}/{len(ambient_temps)}")
    print(f"Converged PLR: {int(np.sum(ok1))}/{len(ambient_temps)}")


if __name__ == "__main__":
    main()
