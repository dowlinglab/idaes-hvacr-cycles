import csv
import math

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import CoolProp.CoolProp as CP

from vapor_compression_plr import SimpleVaporCompressionCyclePLR, Mode


def main():
    cold_storage_setpoint_c = -20.0
    evap_sat_c = cold_storage_setpoint_c - 10.0  # -30 C
    condenser_approach_c = 10.0
    ambient_temps = np.arange(15.0, 46.0, 5.0)
    cp_fluid = "R134a"

    # Rounded fixed bounds requested by user (kPa).
    low_side_bounds = (60.0, 120.0)
    high_side_bounds = (500.0, 1000.0)

    cycle = SimpleVaporCompressionCyclePLR(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        mode=Mode.PH,
    )
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    base_kwargs = dict(
        low_side_pressure=low_side_bounds,
        high_side_pressure=high_side_bounds,
        evaporator_temperature=(-45, -15),
        condenser_temperature=(15, 50),
        subcooling=3,
        superheating=3,
        max_pressure_ratio=20,
        plr=0.75,
        cd=0.13,
    )

    cop_vals = []
    converged_vals = []
    cop_carnot_vals = []

    for ambient_c in ambient_temps:
        TL_K = evap_sat_c + 273.15
        TH_K = ambient_c + condenser_approach_c + 273.15
        cop_carnot_vals.append(TL_K / (TH_K - TL_K))

        cycle.set_specifications(
            ambient_temperature=float(ambient_c),
            condenser_approach=condenser_approach_c,
            evap_sat_temperature=evap_sat_c,
            **base_kwargs,
        )

        try:
            _, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
            cop_use = cycle.get_part_load_cop()
            if (not converged) or (not math.isfinite(cop_use)):
                cop_use = float("nan")
                converged = False
        except Exception:
            try:
                cycle.set_specifications(debug_disable_arc_pressure_eq=True)
                _, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
                cop_use = cycle.get_part_load_cop()
                if (not converged) or (not math.isfinite(cop_use)):
                    cop_use = float("nan")
                    converged = False
            except Exception:
                cop_use = float("nan")
                converged = False

        cop_vals.append(float(cop_use) if math.isfinite(cop_use) else float("nan"))
        converged_vals.append(bool(converged))

    out_csv = "cop_vs_ambient_plr_cold_storage_r134a_copy.csv"
    out_png = "cop_vs_ambient_plr_cold_storage_r134a_copy.png"
    out_pdf = "cop_vs_ambient_plr_cold_storage_r134a_copy.pdf"

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "ambient_C", "cop_part_load", "cop_carnot_ref", "converged",
            "cold_storage_setpoint_C", "evap_sat_setpoint_C", "cond_sat_target_C",
        ])
        for Ta, cop, cop_carnot, ok in zip(ambient_temps, cop_vals, cop_carnot_vals, converged_vals):
            w.writerow([
                float(Ta), cop, float(cop_carnot), int(ok),
                cold_storage_setpoint_c, evap_sat_c, float(Ta + condenser_approach_c),
            ])

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    ok_mask = np.array(converged_vals, dtype=bool) & np.isfinite(cop_vals)
    bad_mask = ~ok_mask

    if np.any(ok_mask):
        ax.plot(ambient_temps[ok_mask], np.array(cop_vals)[ok_mask], marker="o", linewidth=2.2, label="Model COP")
    if np.any(bad_mask):
        ax.plot(ambient_temps[bad_mask], np.array(cop_vals)[bad_mask], marker="x", linestyle="none", label="Not converged")

    ax.plot(ambient_temps, cop_carnot_vals, linestyle="--", linewidth=2.0, color="black", label="Carnot reference")
    ax.set_title("R134a: PLR COP vs Ambient (Cold Storage)")
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
    print(f"Converged points: {int(np.sum(ok_mask))}/{len(ambient_temps)}")


if __name__ == "__main__":
    main()
