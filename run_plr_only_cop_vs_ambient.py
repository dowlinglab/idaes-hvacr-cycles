"""
COP vs ambient runner for the plain IDAES PLR-only cycle model.

This runner intentionally uses `vapor_compression_plr.py` without HX-copy
blocks so the resulting curves reflect the current IDAES PLR-only path.
"""

import csv
import math

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr import Mode, SimpleVaporCompressionCyclePLR


AMBIENT_TEMPS_C = np.arange(10.0, 46.0, 5.0)
PLR_VALUE = 0.75
CD_VALUE = 0.13
EVAP_SAT_TARGET_C = -29.0
COND_APPROACH_C = 9.0
SUPERHEAT_C = 3.0
SUBCOOL_C = 3.0


def _build_cycle(fluid_name: str):
    candidates = [fluid_name]
    if fluid_name.lower() == "r1234ze(e)":
        candidates.extend(["R1234ZEE", "R1234zeE", "R1234zee", "r1234ze"])

    last_err = None
    for name in candidates:
        try:
            cycle = SimpleVaporCompressionCyclePLR(
                fluid_name=name,
                compressor_efficiency=0.75,
                PLR=PLR_VALUE,
                CD=CD_VALUE,
                mode=Mode.PH,
            )
            return cycle, name
        except Exception as exc:
            last_err = exc
    raise RuntimeError(f"Failed to build PLR-only cycle for {fluid_name}: {last_err}")


def _solve_with_retry(cycle, run_kwargs):
    try:
        cycle.set_specifications(**run_kwargs)
        cop_full, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        return cop_full, converged
    except Exception:
        try:
            retry_kwargs = dict(run_kwargs)
            retry_kwargs["debug_disable_arc_pressure_eq"] = True
            cycle.set_specifications(**retry_kwargs)
            cop_full, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
            return cop_full, converged
        except Exception:
            return float("nan"), False


def _run_fluid(fluid_name: str, out_stem: str):
    cycle, fluid_used = _build_cycle(fluid_name)
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    kwargs_base = dict(
        low_side_pressure=(60.0, 200.0),
        high_side_pressure=(500.0, 4000.0),
        evaporator_temperature=(-55.0, 0.0),
        condenser_temperature=(15.0, 60.0),
        evap_sat_temperature=EVAP_SAT_TARGET_C,
        condenser_approach=COND_APPROACH_C,
        superheating=SUPERHEAT_C,
        subcooling=SUBCOOL_C,
        max_pressure_ratio=20.0,
        plr=PLR_VALUE,
        cd=CD_VALUE,
    )

    cop_full_vals = []
    cop_part_vals = []
    cop_carnot_vals = []
    converged_vals = []

    for ambient_c in AMBIENT_TEMPS_C:
        tl_k = EVAP_SAT_TARGET_C + 273.15
        th_k = ambient_c + COND_APPROACH_C + 273.15
        cop_carnot_vals.append(tl_k / max(th_k - tl_k, 1.0e-9))

        run_kwargs = dict(kwargs_base)
        run_kwargs["ambient_temperature"] = float(ambient_c)

        cop_full, ok = _solve_with_retry(cycle, run_kwargs)
        cop_part = cycle.get_part_load_cop() if ok else float("nan")

        cop_full_vals.append(float(cop_full) if (ok and math.isfinite(cop_full)) else float("nan"))
        cop_part_vals.append(float(cop_part) if (ok and math.isfinite(cop_part)) else float("nan"))
        converged_vals.append(bool(ok))

        print(
            f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | "
            f"converged={int(ok)} | COP_part={cop_part_vals[-1]:.4f}"
        )

    out_csv = f"{out_stem}.csv"
    out_png = f"{out_stem}.png"
    out_pdf = f"{out_stem}.pdf"

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "ambient_C",
                "cop_plr_only_full",
                "cop_plr_only_part",
                "cop_carnot_ref",
                "converged",
                "fluid_name_used",
                "evap_sat_target_C",
                "cond_approach_C",
            ]
        )
        for ta, cfull, cpart, cc, ok in zip(
            AMBIENT_TEMPS_C, cop_full_vals, cop_part_vals, cop_carnot_vals, converged_vals
        ):
            w.writerow([float(ta), cfull, cpart, cc, int(ok), fluid_used, EVAP_SAT_TARGET_C, COND_APPROACH_C])

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    ok_mask = np.array(converged_vals, dtype=bool) & np.isfinite(cop_part_vals)
    if np.any(ok_mask):
        ax.plot(
            AMBIENT_TEMPS_C[ok_mask],
            np.array(cop_part_vals)[ok_mask],
            marker="s",
            linewidth=2.2,
            label="PLR-only COP",
        )
    ax.plot(
        AMBIENT_TEMPS_C,
        cop_carnot_vals,
        linestyle="--",
        linewidth=2.0,
        color="black",
        label="Carnot COP",
    )
    ax.set_title(f"{fluid_used}: IDAES PLR-only COP vs Ambient")
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
    print(f"Converged points: {int(np.sum(ok_mask))}/{len(AMBIENT_TEMPS_C)}")


def main():
    _run_fluid("R134a", "cop_vs_ambient_plr_only_r134a")
    try:
        _run_fluid("R1234ze(E)", "cop_vs_ambient_plr_only_r1234zee")
    except Exception as exc:
        print(f"Skipped R1234ze(E): {exc}")


if __name__ == "__main__":
    main()
