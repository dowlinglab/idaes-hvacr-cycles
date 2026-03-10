import csv
import math

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr_epsntu_copy import (
    Mode,
    SimpleVaporCompressionCyclePLREpsNTU,
)


def _solve_with_retry(cycle, run_kwargs):
    try:
        cop_full, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        return cop_full, converged
    except Exception:
        return float("nan"), False


def main():
    cold_storage_setpoint_c = -20.0
    ambient_temps = np.arange(10.0, 46.0, 5.0)

    # Baseline: match non-IDAES lumped model defaults.
    cycle = SimpleVaporCompressionCyclePLREpsNTU(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        PLR=0.75,
        CD=0.13,
        mode=Mode.PH,
        UA_evap_total=1500.0,
        UA_cond_total=1800.0,
        UA_scale=1.0,
        m_dot_ref=0.02,
        m_dot_air_evap=1.2,
        m_dot_air_cond=1.5,
        cp_air_evap=1006.0,
        cp_air_cond=1006.0,
    )

    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    kwargs_base = dict(
        low_side_pressure=(60.0, 200.0),
        high_side_pressure=(500.0, 4000.0),
        evaporator_temperature=(
            cold_storage_setpoint_c - 10.0,
            cold_storage_setpoint_c - 8.0,
        ),
        condenser_temperature=(18.0, 20.0),
        plr=0.75,
        cd=0.13,
        UA_evap_total=1500.0,
        UA_cond_total=1800.0,
        UA_scale=1.0,
    )

    cop_full_vals = []
    cop_part_vals = []
    cop_carnot_vals = []
    converged_vals = []
    sh_vals = []
    sc_vals = []

    for ambient_c in ambient_temps:
        tl_k = cold_storage_setpoint_c + 273.15
        th_k = ambient_c + 9.0 + 273.15
        cop_carnot_vals.append(tl_k / (th_k - tl_k))

        run_kwargs = dict(kwargs_base)
        run_kwargs.update(
            dict(
                ambient_temperature=float(ambient_c),
                condenser_temperature=(float(ambient_c + 8.0), float(ambient_c + 10.0)),
                cold_storage_setpoint=float(cold_storage_setpoint_c),
            )
        )
        cycle.set_specifications(**run_kwargs)
        cop_full, ok = _solve_with_retry(cycle, run_kwargs)
        cop_part = cycle.get_part_load_cop() if ok else float("nan")
        diag = cycle._last_result.diagnostics if (ok and cycle._last_result is not None) else {}

        cop_full_vals.append(float(cop_full) if (ok and math.isfinite(cop_full)) else float("nan"))
        cop_part_vals.append(float(cop_part) if (ok and math.isfinite(cop_part)) else float("nan"))
        converged_vals.append(bool(ok))
        sh_vals.append(float(diag.get("SH_actual", float("nan"))))
        sc_vals.append(float(diag.get("SC_actual", float("nan"))))

    out_csv = "cop_vs_ambient_plr_cold_storage_r134a_epsntu_copy.csv"
    out_png = "cop_vs_ambient_plr_cold_storage_r134a_epsntu_copy.png"
    out_pdf = "cop_vs_ambient_plr_cold_storage_r134a_epsntu_copy.pdf"

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "ambient_C",
                "cop_epsntu_full",
                "cop_epsntu_part",
                "cop_carnot_ref",
                "SH_actual_K",
                "SC_actual_K",
                "converged",
            ]
        )
        for ta, cfull, cpart, cc, sh, sc, ok in zip(
            ambient_temps, cop_full_vals, cop_part_vals, cop_carnot_vals, sh_vals, sc_vals, converged_vals
        ):
            w.writerow([float(ta), cfull, cpart, cc, sh, sc, int(ok)])

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    ok_mask = np.array(converged_vals, dtype=bool) & np.isfinite(cop_part_vals)
    if np.any(ok_mask):
        ax.plot(
            ambient_temps[ok_mask],
            np.array(cop_part_vals)[ok_mask],
            marker="s",
            linewidth=2.2,
            label="PLR + HX (eps-NTU copy)",
        )
    ax.plot(
        ambient_temps,
        cop_carnot_vals,
        linestyle="--",
        linewidth=2.0,
        color="black",
        label="Carnot COP",
    )
    ax.set_title("R134a: Non-IDAES Epsilon-NTU Copy vs Carnot")
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
    print(f"Converged points: {int(np.sum(ok_mask))}/{len(ambient_temps)}")


if __name__ == "__main__":
    main()
