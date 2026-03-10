import csv
import math

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr_hx_0d import (
    Mode,
    SimpleVaporCompressionCyclePLRHX0D,
)


def _solve_with_retry(cycle, run_kwargs):
    try:
        cop_full, converged = cycle.optimize_COP(verbose=False, initialize=False, optimize=False)
        return cop_full, converged
    except Exception:
        try:
            retry_kwargs = dict(run_kwargs)
            retry_kwargs.update(dict(debug_disable_arc_pressure_eq=True))
            cycle.set_specifications(**retry_kwargs)
            cop_full, converged = cycle.optimize_COP(verbose=False, initialize=False, optimize=False)
            return cop_full, converged
        except Exception:
            return float("nan"), False


def main():
    cold_storage_setpoint_c = -20.0
    ambient_temps = np.arange(10.0, 46.0, 5.0)

    low_side_bounds = (60.0, 200.0)
    high_side_bounds = (500.0, 4000.0)

    cycle_lc = SimpleVaporCompressionCyclePLRHX0D(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        PLR=0.75,
        CD=0.13,
        mode=Mode.PH,
    )

    cycle_lc.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)

    lc_kwargs = dict(
        low_side_pressure=low_side_bounds,
        high_side_pressure=high_side_bounds,
        evaporator_temperature=(
            cold_storage_setpoint_c - 10.0,
            cold_storage_setpoint_c - 8.0,
        ),
        condenser_temperature=(18.0, 20.0),
        subcooling=3,
        superheating=3,
        max_pressure_ratio=20,
        plr=0.75,
        cd=0.13,
    )

    cop_full_vals = []
    cop_part_vals = []
    cop_carnot_vals = []
    converged_lc = []
    sh_vals = []
    sc_vals = []

    for ambient_c in ambient_temps:
        tl_k = cold_storage_setpoint_c + 273.15
        th_k = ambient_c + 9.0 + 273.15
        cop_carnot_vals.append(tl_k / (th_k - tl_k))

        kwargs_lc = dict(lc_kwargs)
        kwargs_lc.update(
            dict(
                ambient_temperature=float(ambient_c),
                condenser_temperature=(float(ambient_c + 8.0), float(ambient_c + 10.0)),
                cold_storage_setpoint=float(cold_storage_setpoint_c),
                evap_offset_bounds=(-10.0, -8.0),
            )
        )
        cycle_lc.set_specifications(**kwargs_lc)
        cop_full, ok_lc = _solve_with_retry(cycle_lc, kwargs_lc)
        cop_part = cycle_lc.get_part_load_cop() if ok_lc else float("nan")
        sh, sc = cycle_lc.get_actual_sh_sc() if ok_lc else (float("nan"), float("nan"))

        cop_full_vals.append(float(cop_full) if (ok_lc and math.isfinite(cop_full)) else float("nan"))
        cop_part_vals.append(float(cop_part) if (ok_lc and math.isfinite(cop_part)) else float("nan"))
        converged_lc.append(bool(ok_lc))
        sh_vals.append(float(sh))
        sc_vals.append(float(sc))

    out_csv = "cop_vs_ambient_plr_cold_storage_r134a_hx_0d.csv"
    out_png = "cop_vs_ambient_plr_cold_storage_r134a_hx_0d.png"
    out_pdf = "cop_vs_ambient_plr_cold_storage_r134a_hx_0d.pdf"

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "ambient_C",
                "cop_lc_full",
                "cop_lc_part",
                "cop_carnot_ref",
                "SH_actual_K",
                "SC_actual_K",
                "converged_lc",
            ]
        )
        for ta, cfull, cpart, cc, sh, sc, ok in zip(
            ambient_temps, cop_full_vals, cop_part_vals, cop_carnot_vals, sh_vals, sc_vals, converged_lc
        ):
            w.writerow([float(ta), cfull, cpart, cc, sh, sc, int(ok)])

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    ok = np.array(converged_lc, dtype=bool) & np.isfinite(cop_part_vals)
    if np.any(ok):
        ax.plot(
            ambient_temps[ok],
            np.array(cop_part_vals)[ok],
            marker="s",
            linewidth=2.2,
            label="PLR + HX (IDAES LC)",
        )
    ax.plot(
        ambient_temps,
        cop_carnot_vals,
        linestyle="--",
        linewidth=2.0,
        color="black",
        label="Carnot COP",
    )

    ax.set_title("R134a: PLR + HX (IDAES LC) vs Carnot")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print(f"Pressure bounds (kPa): low={low_side_bounds}, high={high_side_bounds}")
    print(f"Converged LC points: {int(np.sum(ok))}/{len(ambient_temps)}")


if __name__ == "__main__":
    main()
