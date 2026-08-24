"""
cop_vs_ambient.py

COP vs. ambient temperature sweep for the R134a (hot) / CO2 (cold)
cascade cycle, using CascadeCycle from vapor_compression_cascade.py.
Same shape of study as the single-stage R134a/R1234yf/R515B benchmark
(R1234yf/cop_vs_ambient_plr_benchmark_combined.csv) -- ambient sweep,
per-point CSV with convergence + state details, Carnot reference line.

For each ambient temperature, a fresh CascadeCycle is built, initialized,
specified, and optimized -- rebuilding per point avoids carrying a stuck
solve/bad basis from one point into the next.

Fixed for the sweep (edit below to change):
    cold_evap_C               = -27 C  (cold loop delivery/evap temperature)
    hot_condenser_approach_C  = 9 C    (hot condenser sat = ambient + approach)
    cascade_approach_dT       = 3 K
    superheat/subcool          = fluid defaults (CO2 [0,5]K, R134a [0,3]K)
    max_pressure_ratio         = 8 (both loops)
    flow_fixed_role             = "cold" @ 1 kg/s (hot loop's flow floats)

Carnot COP here uses the SYSTEM's overall reservoirs (cold_evap_C on the
cold end, ambient_C on the hot end) -- the theoretical best case for the
cascade as a whole, not per-loop. Real cascade COP is necessarily below
this since it stacks two real compression stages plus the cascade HX's
own approach-temperature penalty on top of the ambient-side condenser
approach.
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from pyomo.environ import value

from vapor_compression_cascade import CascadeCycle, DEFAULT_FLUIDS, C_to_K


def _loop_diagnostics(cycle, role):
    fs = cycle.model.fs
    evaporator, compressor, condenser, expansion_valve = cycle._unit_operations[role]
    evap_out = evaporator.control_volume.properties_out[0]
    cond_out = condenser.control_volume.properties_out[0]

    return {
        f"{role}_mass_flow": value(evaporator.inlet.flow_mass[0]),
        f"{role}_evap_T_sat_C": value(evap_out.temperature_sat) - C_to_K,
        f"{role}_evap_superheat_actual_K": value(evap_out.temperature) - value(evap_out.temperature_sat),
        f"{role}_evap_duty": value(evaporator.heat_duty[0]),
        f"{role}_comp_ratioP": value(compressor.ratioP[0]),
        f"{role}_comp_work": value(compressor.work_mechanical[0]),
        f"{role}_cond_T_sat_C": value(cond_out.temperature_sat) - C_to_K,
        f"{role}_cond_subcool_actual_K": value(cond_out.temperature_sat) - value(cond_out.temperature),
        f"{role}_cond_duty": value(condenser.heat_duty[0]),
    }


def run_sweep(ambient_C_values, cold_evap_C=-27.0, hot_condenser_approach_C=9.0,
              cascade_approach_dT=3.0, verbose=False):
    results = []
    T_cold_K = cold_evap_C + C_to_K

    for ambient_C in ambient_C_values:
        print(f"--- ambient = {ambient_C} C ---")
        T_hot_K = ambient_C + C_to_K
        cop_carnot = T_cold_K / (T_hot_K - T_cold_K)

        row = {"ambient_C": ambient_C, "cop_carnot": cop_carnot,
               "converged": False, "cop_full": None, "cop_part": None,
               "second_law_efficiency": None,
               "cascade_approach_actual_K": None}

        cycle = CascadeCycle(fluids=DEFAULT_FLUIDS)
        try:
            cycle.specify_initial_conditions(hot_ambient_C=ambient_C, cold_evap_C=cold_evap_C)
            cycle.initialize(verbose=False)
            cycle.set_specifications(
                hot={"ambient_temperature": ambient_C, "condenser_approach": hot_condenser_approach_C},
                cold={"evap_sat_temperature": cold_evap_C},
                cascade_approach_dT=cascade_approach_dT,
            )
            cop, converged = cycle.optimize_COP(verbose=verbose)

            row["converged"] = bool(converged)
            if converged:
                row["cop_full"] = cycle._last_cop_full
                row["cop_part"] = cycle._last_cop_part
                row["second_law_efficiency"] = cycle._last_cop_full / cop_carnot
                for role in cycle.ROLES:
                    row.update(_loop_diagnostics(cycle, role))
                row["cascade_approach_actual_K"] = row["cold_cond_T_sat_C"] - row["hot_evap_T_sat_C"]
        except Exception as exc:
            print(f"  FAILED at ambient={ambient_C} C: {exc}")
            row["converged"] = False

        print(f"  converged={row['converged']}, cop_full={row['cop_full']}, cop_part={row['cop_part']}")
        results.append(row)

    return results


def write_csv(results, out_path="cop_vs_ambient_cascade.csv"):
    fieldnames = []
    for row in results:
        for k in row:
            if k not in fieldnames:
                fieldnames.append(k)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {out_path}")


def plot_results(results, save_path=None):
    converged = [r for r in results if r["converged"]]
    failed = [r for r in results if not r["converged"]]

    plt.figure()
    plt.plot([r["ambient_C"] for r in results], [r["cop_carnot"] for r in results],
              'k--', label="Carnot (system)")
    plt.plot([r["ambient_C"] for r in converged], [r["cop_full"] for r in converged],
              'o-', label="Cascade COP (full load)")
    plt.plot([r["ambient_C"] for r in converged], [r["cop_part"] for r in converged],
              's-', label="Cascade COP (PLR/CD part-load)")
    if failed:
        ylim = plt.gca().get_ylim()
        plt.plot([r["ambient_C"] for r in failed], [ylim[0]] * len(failed), 'rx', label="Did not converge")

    plt.xlabel("Ambient temperature (C)")
    plt.ylabel("COP")
    plt.title("R134a/CO2 cascade COP vs. ambient temperature")
    plt.grid(True)
    plt.legend()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved plot to {save_path}")
    plt.show()


if __name__ == "__main__":
    ambient_C_values = np.arange(10, 46, 5)  # 10 to 45 C, 5 C steps
    results = run_sweep(ambient_C_values)
    write_csv(results, out_path="cop_vs_ambient_cascade.csv")
    plot_results(results, save_path="cop_vs_ambient_cascade.png")
