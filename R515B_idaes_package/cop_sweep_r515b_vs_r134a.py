"""
COP sweep comparison: R-515B (our new validated integration) vs. R134a
(the original vapor_compression.py, using the real HelmholtzParameterBlock)
-- both run through the SAME simplified vapor-compression cycle model
(evaporator -> compressor -> condenser -> expansion valve), swept over
evaporating temperature at a fixed condensing temperature, so the
comparison isolates the refrigerant's own thermodynamic behavior rather
than any difference in cycle assumptions.

Condensing side held fixed via ambient_temperature=35 degC + condenser_
approach=5 degC (condensing sat T = 40 degC) for both fluids. Evaporating
side swept via evap_sat_temperature over a range comfortably inside R-515B's
already-validated saturation-curve grid (255-375K, i.e. -18.15 to 101.85
degC) so no new, unvalidated territory is explored for R-515B.

Not part of the master task's Stage A-O build -- this is a downstream
analysis using the already-validated model, run in response to a direct
user question ("how does the COP sweep compare to R134a's?").
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
PROJECT_DIR = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PROJECT_DIR))

import json  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")

from vapor_compression_r515b_integration import R515BVaporCompressionCycle  # noqa: E402
from vapor_compression import SimpleVaporCompressionCycle, Mode  # noqa: E402

T_EVAP_LIST_C = [-20.0, -15.0, -10.0, -5.0, 0.0, 5.0]
T_COND_AMBIENT_C = 35.0
T_COND_APPROACH_C = 5.0  # condensing sat T = 40 degC
SUBCOOL_C = 3.0
SUPERHEAT_C = 3.0
MAX_PRESSURE_RATIO = 8.0
COMPRESSOR_EFFICIENCY = 0.75


def run_r515b_point(t_evap_c):
    cyc = R515BVaporCompressionCycle(compressor_efficiency=COMPRESSOR_EFFICIENCY)
    cyc.specify_initial_conditions(low_side_temperature=t_evap_c, high_side_temperature=T_COND_AMBIENT_C + T_COND_APPROACH_C)
    cyc.initialize(verbose=False)
    cyc.set_specifications(
        evap_sat_temperature=t_evap_c,
        ambient_temperature=T_COND_AMBIENT_C, condenser_approach=T_COND_APPROACH_C,
        subcooling=SUBCOOL_C, superheating=SUPERHEAT_C, max_pressure_ratio=MAX_PRESSURE_RATIO,
    )
    cop, converged = cyc.optimize_COP(verbose=False, initialize=True, optimize=True)
    return float(cop), bool(converged)


def run_r134a_point(t_evap_c):
    cyc = SimpleVaporCompressionCycle("R134a", compressor_efficiency=COMPRESSOR_EFFICIENCY, mode=Mode.PH)
    cyc.specify_initial_conditions(low_side_temperature=t_evap_c, high_side_temperature=T_COND_AMBIENT_C + T_COND_APPROACH_C)
    cyc.initialize(verbose=False)
    cyc.set_specifications(
        evap_sat_temperature=t_evap_c,
        ambient_temperature=T_COND_AMBIENT_C, condenser_approach=T_COND_APPROACH_C,
        subcooling=SUBCOOL_C, superheating=SUPERHEAT_C, max_pressure_ratio=MAX_PRESSURE_RATIO,
    )
    cop, converged = cyc.optimize_COP(verbose=False, initialize=True, optimize=True)
    return float(cop), bool(converged)


def main():
    results = {"T_evap_C": T_EVAP_LIST_C, "R515B": [], "R134a": []}
    for t in T_EVAP_LIST_C:
        try:
            cop, ok = run_r515b_point(t)
        except Exception as e:
            cop, ok = None, False
            print(f"[R515B] T_evap={t}C EXCEPTION: {e}")
        results["R515B"].append({"cop": cop, "converged": ok})
        print(f"[R515B] T_evap={t:6.1f}C  COP={cop}  converged={ok}")

    for t in T_EVAP_LIST_C:
        try:
            cop, ok = run_r134a_point(t)
        except Exception as e:
            cop, ok = None, False
            print(f"[R134a] T_evap={t}C EXCEPTION: {e}")
        results["R134a"].append({"cop": cop, "converged": ok})
        print(f"[R134a] T_evap={t:6.1f}C  COP={cop}  converged={ok}")

    out_path = HERE / "cop_sweep_r515b_vs_r134a_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
