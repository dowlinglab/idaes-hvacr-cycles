import importlib
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

MODULE_NAME = "vapor_compression_copy"


def _disable_plots():
    try:
        plt.show = lambda *args, **kwargs: None
    except Exception:
        pass


def _write_kv(path, data):
    lines = []
    for k, v in data.items():
        lines.append(f"{k}: {v}")
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    _disable_plots()
    report_lines = []
    values = {}
    status = "PASS"

    try:
        mod = importlib.import_module(MODULE_NAME)
        SimpleVaporCompressionCycle = mod.SimpleVaporCompressionCycle
        Mode = mod.Mode
    except Exception as exc:
        report_lines.append(f"IMPORT ERROR: {exc}")
        status = "FAIL"
        Path("cop_test_report.txt").write_text("\n".join(report_lines) + "\n")
        _write_kv("cop_test_values.txt", values)
        return 1

    try:
        cycle = SimpleVaporCompressionCycle(
            "R134a", compressor_efficiency=0.75, mode=Mode.IMPROVED_TPX
        )
        cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
        cycle.initialize(verbose=False)

        def apply_specs(flag):
            cycle.set_specifications(
                ambient_temperature=35,
                condenser_approach=5,
                evap_sat_temperature=-10,
                superheating=3,
                subcooling=3,
                max_pressure_ratio=4,
                debug_disable_arc_pressure_eq=flag,
            )

        apply_specs(False)
        cop, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        if not converged:
            apply_specs(True)
            cop, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)

        m = cycle.model
        Q_evap = m.fs.evaporator.heat_duty[0].value
        Q_cond = m.fs.condenser.heat_duty[0].value
        W_comp = m.fs.compressor.work_mechanical[0].value

        values.update(
            {
                "termination": getattr(cycle, "last_termination_condition", None),
                "converged": converged,
                "Q_evap": Q_evap,
                "Q_cond": Q_cond,
                "W_comp": W_comp,
                "COP_var": m.fs.cop.value,
            }
        )

        if W_comp is None or abs(W_comp) < 1e-8 or Q_evap is None:
            cop_calc = 0.0
        else:
            cop_calc = Q_evap / W_comp
        values["COP_calc"] = cop_calc

        # Sign checks
        sign_errors = []
        if Q_evap is not None and Q_evap < 0:
            sign_errors.append("Evaporator heat duty is negative")
        if Q_cond is not None and Q_cond > 0:
            sign_errors.append("Condenser heat duty is not negative")
        if W_comp is not None and W_comp <= 0:
            sign_errors.append("Compressor work is non-positive")

        if sign_errors:
            status = "FAIL"
            report_lines.append("SIGN ERROR:")
            report_lines.extend([f"- {e}" for e in sign_errors])

        # Carnot comparison
        carnot = None
        try:
            carnot = cycle.compute_carnot()
        except Exception:
            carnot = None
        values["Carnot"] = carnot
        if carnot is not None and cop_calc is not None:
            if cop_calc > carnot * 1.01:
                status = "FAIL"
                report_lines.append(
                    f"COP exceeds Carnot: COP={cop_calc}, Carnot={carnot}"
                )

        report_lines.append(f"STATUS: {status}")

    except Exception as exc:
        status = "FAIL"
        report_lines.append(f"RUNTIME ERROR: {exc}")

    Path("cop_test_report.txt").write_text("\n".join(report_lines) + "\n")
    _write_kv("cop_test_values.txt", values)

    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
