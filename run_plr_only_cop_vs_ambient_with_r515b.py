"""
COP vs ambient runner for the plain IDAES PLR-only cycle model --
extended to include R-515B alongside R134a, mirroring exactly how
`run_plr_only_cop_vs_ambient_with_r1234yf.py` extended the original
runner to include R1234yf.

This is a NEW file, not an edit of `run_plr_only_cop_vs_ambient.py` or
`run_plr_only_cop_vs_ambient_with_r1234yf.py` (both belong to the
existing PLR project and are left untouched). It reuses the EXACT same
setpoints, methodology, and output format as those runners -- same
AMBIENT_TEMPS_C sweep, same PLR/CD values, same evaporator saturation
target and condenser approach -- and adds R-515B as an additional
fluid, wired to the R515B_idaes_package/ project's own validated,
oracle-traced Helmholtz property package (via
`vapor_compression_plr_r515b.py` in that subdirectory).

R134a still goes through the original, unmodified `vapor_compression_
plr.py` exactly as the original runner does.

One deliberate, documented adaptation for the R-515B leg only: its
`specify_initial_conditions` generic warm-start seed uses
high_side_temperature=40 instead of the original runner's 30. This is
NOT a change to the benchmark itself (the actual swept/fixed setpoints
below -- ambient grid, PLR, CD, evap_sat target, condenser approach,
superheat/subcool -- are identical for both fluids); it only changes
the rough initial guess used to warm-start the flowsheet before
`set_specifications`/`optimize_COP` retarget it to the real conditions.
Confirmed empirically (see PROJECT_CONTEXT.md, 2026-08-19) that R-515B's
compressor unit-model initialization raises `InitializationError` at the
original runner's high_side_temperature=30 generic seed regardless of
low_side_temperature (reproduced identically against the UNMODIFIED
`vapor_compression_r515b_integration.py`, so this is not something the
PLR grafting introduced) but succeeds cleanly at 40 -- so 40 is used for
the R-515B leg's warm start only, exactly as Stage N's own validation
script (`validate_stage_n_integration.py`) already did before this file
existed.
"""

import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr import Mode, SimpleVaporCompressionCyclePLR

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_R515B_DIR = os.path.join(_THIS_DIR, "R515B_idaes_package")
if _R515B_DIR not in sys.path:
    sys.path.insert(0, _R515B_DIR)

from vapor_compression_plr_r515b import R515BVaporCompressionCyclePLR  # noqa: E402


# Setpoints -- identical to run_plr_only_cop_vs_ambient.py /
# run_plr_only_cop_vs_ambient_with_r1234yf.py (collaborator-sourced)
AMBIENT_TEMPS_C = np.arange(10.0, 46.0, 5.0)
PLR_VALUE = 0.75
CD_VALUE = 0.13
EVAP_SAT_TARGET_C = -29.0
COND_APPROACH_C = 9.0
SUPERHEAT_C = 3.0
SUBCOOL_C = 3.0


def _build_cycle(fluid_name: str):
    if fluid_name == "R515B":
        cycle = R515BVaporCompressionCyclePLR(
            compressor_efficiency=0.75,
            PLR=PLR_VALUE,
            CD=CD_VALUE,
        )
        return cycle, "R515B"

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
    if fluid_name == "R515B":
        # See module docstring: R-515B-only adaptation of the generic
        # warm-start seed (40 instead of the original runner's 30);
        # the actual benchmark setpoints below are unchanged.
        cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=40)
    else:
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

    # Solve in DESCENDING ambient order rather than the AMBIENT_TEMPS_C
    # array's own ascending order. Per the user's own diagnosis (this is a
    # warm-starting issue, not a physical infeasibility): the model/solver
    # state persists across `set_specifications`/`optimize_COP` calls on
    # the SAME `cycle` object between loop iterations (no rebuild), so
    # each point after the first is effectively warm-started from the
    # PREVIOUS point solved, not from the generic cold `specify_initial_
    # conditions` seed. In ascending order the two historically-hardest
    # points (10C, 15C) are solved FIRST, straight from the generic cold
    # seed, with nothing better to warm-start from. Solving descending
    # (45C -> 10C) instead means 10C/15C are only ever attempted after
    # their nearest neighbor (15C/20C) has already converged, so they
    # inherit a much closer-to-correct starting point. Results are
    # collected in a dict and re-assembled in the original ascending
    # order afterward, so the CSV/plot output shape is unchanged.
    results_by_ambient = {}
    for ambient_c in sorted(AMBIENT_TEMPS_C, reverse=True):
        run_kwargs = dict(kwargs_base)
        run_kwargs["ambient_temperature"] = float(ambient_c)

        cop_full, ok = _solve_with_retry(cycle, run_kwargs)
        cop_part = cycle.get_part_load_cop() if ok else float("nan")

        results_by_ambient[float(ambient_c)] = (
            float(cop_full) if (ok and math.isfinite(cop_full)) else float("nan"),
            float(cop_part) if (ok and math.isfinite(cop_part)) else float("nan"),
            bool(ok),
        )

        print(
            f"{fluid_used:>10s} | Ambient {ambient_c:>5.1f} C | "
            f"converged={int(ok)} | COP_part={results_by_ambient[float(ambient_c)][1]:.4f}",
            flush=True,
        )

    cop_full_vals = []
    cop_part_vals = []
    cop_carnot_vals = []
    converged_vals = []

    for ambient_c in AMBIENT_TEMPS_C:
        tl_k = EVAP_SAT_TARGET_C + 273.15
        th_k = ambient_c + COND_APPROACH_C + 273.15
        cop_carnot_vals.append(tl_k / max(th_k - tl_k, 1.0e-9))

        cfull, cpart, ok = results_by_ambient[float(ambient_c)]
        cop_full_vals.append(cfull)
        cop_part_vals.append(cpart)
        converged_vals.append(ok)

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

    return AMBIENT_TEMPS_C, cop_full_vals, cop_part_vals, converged_vals


def _plot_comparison(results):
    """Overlay all fluids' part-load COP curves on one figure."""
    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=160)
    for fluid_used, (ambient_c, cop_full_vals, cop_part_vals, converged_vals) in results.items():
        ok_mask = np.array(converged_vals, dtype=bool) & np.isfinite(cop_part_vals)
        if np.any(ok_mask):
            ax.plot(
                ambient_c[ok_mask],
                np.array(cop_part_vals)[ok_mask],
                marker="s",
                linewidth=2.2,
                label=fluid_used,
            )
    ax.set_title("PLR-only part-load COP vs Ambient -- R134a vs R-515B")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP (part-load)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig("cop_vs_ambient_plr_only_comparison_r515b.png", bbox_inches="tight")
    fig.savefig("cop_vs_ambient_plr_only_comparison_r515b.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved: cop_vs_ambient_plr_only_comparison_r515b.png")
    print("Saved: cop_vs_ambient_plr_only_comparison_r515b.pdf")


def main():
    results = {}
    results["R134a"] = _run_fluid("R134a", "cop_vs_ambient_plr_only_r134a_for_r515b_compare")
    try:
        results["R515B"] = _run_fluid("R515B", "cop_vs_ambient_plr_only_r515b")
    except Exception as exc:
        print(f"Skipped R515B: {exc}")

    if len(results) > 1:
        _plot_comparison(results)


if __name__ == "__main__":
    main()
