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
class StrictBand:
    evap_approach_points_c: tuple[tuple[float, float], ...]
    low_side_pressure_kpa: tuple[float, float]
    high_side_pressure_kpa: tuple[float, float]


def _solve(cycle, kwargs):
    try:
        cycle.set_specifications(**kwargs)
        _, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        cop = cycle.get_part_load_cop() if ok else float("nan")
        return bool(ok and math.isfinite(cop)), float(cop) if math.isfinite(cop) else float("nan")
    except Exception:
        try:
            kw = dict(kwargs)
            kw["debug_disable_arc_pressure_eq"] = True
            cycle.set_specifications(**kw)
            _, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
            cop = cycle.get_part_load_cop() if ok else float("nan")
            return bool(ok and math.isfinite(cop)), float(cop) if math.isfinite(cop) else float("nan")
        except Exception:
            return False, float("nan")


def eval_strict(cycle, strict, ambient_c, sc, sh, common):
    attempts = 0
    feasible = 0
    best_cop = float("nan")
    for evap_sat_c, cond_app_c in strict.evap_approach_points_c:
        attempts += 1
        kw = dict(common)
        kw.update(
            dict(
                low_side_pressure=strict.low_side_pressure_kpa,
                high_side_pressure=strict.high_side_pressure_kpa,
                ambient_temperature=float(ambient_c),
                evap_sat_temperature=float(evap_sat_c),
                condenser_approach=float(cond_app_c),
                subcooling=float(sc),
                superheating=float(sh),
            )
        )
        ok, cop = _solve(cycle, kw)
        if ok:
            feasible += 1
            if (not math.isfinite(best_cop)) or (cop > best_cop):
                best_cop = cop
    return attempts, feasible, best_cop


def eval_repo_relaxed(cycle, ambient_c, sc, sh, common, tsp_c):
    # Updated relaxed definition per user:
    # T_evap_sat in (Tsp-30, Tsp) with condenser_approach fixed at 20 C.
    evap_sat_candidates = np.arange(tsp_c - 30.0, tsp_c + 1.0, 1.0)
    attempts = 0
    feasible = 0
    best_cop = float("nan")
    for evap_sat_c in evap_sat_candidates:
        attempts += 1
        kw = dict(common)
        kw.update(
            dict(
                low_side_pressure=(60.0, 120.0),
                high_side_pressure=(500.0, 1000.0),
                # Relaxed sweep needs headroom above Tsp because outlet temp = T_sat + superheating.
                evaporator_temperature=(tsp_c - 30.0, tsp_c + 10.0),
                condenser_temperature=(15.0, 50.0),
                ambient_temperature=float(ambient_c),
                evap_sat_temperature=float(evap_sat_c),
                condenser_approach=20.0,
                subcooling=float(sc),
                superheating=float(sh),
            )
        )
        ok, cop = _solve(cycle, kw)
        if ok:
            feasible += 1
            if (not math.isfinite(best_cop)) or (cop > best_cop):
                best_cop = cop
    return attempts, feasible, best_cop


def main():
    tsp_c = -20.0
    ambient = np.arange(10.0, 46.0, 5.0)
    pairs = [
        (1.0, 1.0), (1.0, 3.0), (1.0, 5.0), (1.0, 7.0),
        (3.0, 3.0), (3.0, 5.0), (3.0, 7.0), (5.0, 5.0), (5.0, 7.0),
    ]

    strict = StrictBand(
        # Exact strict set requested by user for Tsp=-20 C:
        # (Tsp-10,8),(Tsp-9,8),(Tsp-8,8),
        # (Tsp-10,9),(Tsp-9,9),(Tsp-8,9),
        # (Tsp-10,10),(Tsp-9,10),(Tsp-8,10)
        evap_approach_points_c=(
            (-30.0, 8.0), (-29.0, 8.0), (-28.0, 8.0),
            (-30.0, 9.0), (-29.0, 9.0), (-28.0, 9.0),
            (-30.0, 10.0), (-29.0, 10.0), (-28.0, 10.0),
        ),
        low_side_pressure_kpa=(60.0, 120.0),
        high_side_pressure_kpa=(500.0, 1000.0),
    )

    common = dict(
        evaporator_temperature=(-55.0, -20.0),
        condenser_temperature=(15.0, 60.0),
        max_pressure_ratio=20.0,
        plr=0.75,
        cd=0.13,
    )

    cycle = SimpleVaporCompressionCyclePLR(fluid_name="R134a", compressor_efficiency=0.75, mode=Mode.PH)
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    rows = []
    for a in ambient:
        for sc, sh in pairs:
            at, fe, cp = eval_strict(cycle, strict, a, sc, sh, common)
            rows.append(dict(ambient_C=float(a), band="Strict", subcooling_C=sc, superheating_C=sh,
                             attempts=at, feasible=fe, feasible_fraction=fe/at if at else 0.0,
                             best_cop=cp if math.isfinite(cp) else float("nan")))

            at, fe, cp = eval_repo_relaxed(cycle, a, sc, sh, common, tsp_c=tsp_c)
            rows.append(dict(ambient_C=float(a), band="Relaxed(repo)", subcooling_C=sc, superheating_C=sh,
                             attempts=at, feasible=fe, feasible_fraction=fe/at if at else 0.0,
                             best_cop=cp if math.isfinite(cp) else float("nan")))

    out_csv = "plr_strict_vs_repo_relaxed_scsh_pairs_r134a.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ambient_C","band","subcooling_C","superheating_C","attempts","feasible","feasible_fraction","best_cop"])
        w.writeheader(); w.writerows(rows)

    # compact plot focused on (1,3) feasibility expectation + all pairs COP traces
    fig, axes = plt.subplots(2, 1, figsize=(10.2, 8.2), dpi=170, sharex=True)
    cmap = plt.get_cmap("tab10")
    colors = {p: cmap(i % 10) for i, p in enumerate(pairs)}

    for ax, band in zip(axes, ["Strict", "Relaxed(repo)"]):
        b = [r for r in rows if r["band"] == band]
        for p in pairs:
            cur = [r for r in b if abs(r["subcooling_C"]-p[0])<1e-9 and abs(r["superheating_C"]-p[1])<1e-9]
            x = np.array([r["ambient_C"] for r in cur])
            y = np.array([r["best_cop"] for r in cur], dtype=float)
            ok = np.isfinite(y)
            lab = f"({int(p[0])},{int(p[1])})"
            if np.any(ok):
                ax.plot(x[ok], y[ok], marker="o", linewidth=1.8, color=colors[p], label=lab)
            if np.any(~ok):
                ax.plot(x[~ok], np.zeros(np.sum(~ok)), "x", color=colors[p], alpha=0.75)
        ax.set_title(band)
        ax.set_ylabel("Best feasible COP")
        ax.grid(True, alpha=0.25)
        ax.legend(ncol=3, fontsize=8, title="(sc,sh)")

    axes[-1].set_xlabel("Ambient temperature (C)")
    fig.suptitle("R134a: Strict vs Repo-relaxed with requested (sc,sh) pairs", y=0.995)
    fig.tight_layout()
    out_png = "plr_strict_vs_repo_relaxed_scsh_pairs_r134a.png"
    out_pdf = "plr_strict_vs_repo_relaxed_scsh_pairs_r134a.pdf"
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")

    # Print specific check the user asked about
    check = [r for r in rows if r["band"]=="Relaxed(repo)" and abs(r["subcooling_C"]-1.0)<1e-9 and abs(r["superheating_C"]-3.0)<1e-9]
    print("Relaxed(repo) (1,3) feasible by ambient:", [int(r["feasible"]) for r in check])
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
