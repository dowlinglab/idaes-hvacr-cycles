"""
R134a Ambient Sweep Runner for PLR + IDAES NTU HX Copy

Author: Shilpa Narasimhan
Codex Support: OpenAI Codex
QA/Testing: Shilpa Narasimhan

Description:
    Runs COP vs ambient (10-45 C) for the IDAES HeatExchanger1D PLR+HX copy.

Context Breadcrumb:
    Uses the frozen bounds policy:
    Tevap in [Tcold_sp-10, Tcold_sp-8], Tcond in [Tamb+8, Tamb+10].
"""

import csv
import math

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr_hx1d import (
    Mode,
    SimpleVaporCompressionCyclePLRHX1D,
)


def _solve_with_retry(cycle, run_kwargs):
    try:
        cop_full, converged = cycle.optimize_COP(verbose=False, initialize=False, optimize=False)
        return cop_full, converged
    except Exception:
        return float("nan"), False


def main():
    cold_storage_setpoint_c = -20.0
    ambient_temps = np.arange(10.0, 46.0, 5.0)

    # IDAES HX1D model copy with frozen PLR settings.
    cycle = SimpleVaporCompressionCyclePLRHX1D(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        PLR=0.75,
        CD=0.13,
        mode=Mode.PH,
        finite_elements=2,
    )

    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)

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
                evap_offset_bounds=(-10.0, -8.0),
            )
        )
        cycle.set_specifications(**run_kwargs)
        cop_full, ok = _solve_with_retry(cycle, run_kwargs)
        cop_part = cycle.get_part_load_cop() if ok else float("nan")
        sh, sc = cycle.get_actual_sh_sc() if ok else (float("nan"), float("nan"))

        cop_full_vals.append(float(cop_full) if (ok and math.isfinite(cop_full)) else float("nan"))
        cop_part_vals.append(float(cop_part) if (ok and math.isfinite(cop_part)) else float("nan"))
        converged_vals.append(bool(ok))
        sh_vals.append(float(sh))
        sc_vals.append(float(sc))

    out_csv = "cop_vs_ambient_plr_cold_storage_r134a_hx1d.csv"
    out_png = "cop_vs_ambient_plr_cold_storage_r134a_hx1d.png"
    out_pdf = "cop_vs_ambient_plr_cold_storage_r134a_hx1d.pdf"

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
            label="PLR + HX (IDAES HX1D)",
        )
    ax.plot(
        ambient_temps,
        cop_carnot_vals,
        linestyle="--",
        linewidth=2.0,
        color="black",
        label="Carnot COP",
    )
    ax.set_title("R134a: IDAES HX1D PLR + HX vs Carnot")
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
