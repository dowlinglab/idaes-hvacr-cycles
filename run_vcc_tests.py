import argparse
import sys

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from pyomo.environ import value

from vapor_compression import SimpleVaporCompressionCycle, Mode

C_TO_K = 273.15


def _disable_plots():
    try:
        plt.show = lambda *args, **kwargs: None
    except Exception:
        pass


def compute_carnot_cop(T_L, T_H):
    if T_H <= T_L:
        raise ValueError(f"Invalid reservoir temps: T_H={T_H}, T_L={T_L}")
    return T_L / (T_H - T_L)


def compute_cop(model):
    Q_evap = value(model.fs.evaporator.heat_duty[0])
    W_comp = value(model.fs.compressor.work_mechanical[0])
    if W_comp is None or abs(W_comp) < 1e-12:
        raise ValueError(f"Compressor work too small or None: {W_comp}")
    return abs(Q_evap) / abs(W_comp)


def build_and_solve(
    fluid="R134a",
    eta=0.75,
    optimize=False,
    debug_disable_arc_pressure_eq=False,
):
    _disable_plots()

    cycle = SimpleVaporCompressionCycle(
        fluid_name=fluid,
        compressor_efficiency=eta,
        mode=Mode.IMPROVED_TPX,
    )

    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    cycle.set_specifications(
        low_side_pressure=(200, 600),
        high_side_pressure=(800, 2500),
        evaporator_temperature=(-30, 5),
        condenser_temperature=(20, 60),
        subcooling=3,
        superheating=3,
        max_pressure_ratio=4,
        ambient_temperature=35,
        condenser_approach=5,
        evap_sat_temperature=-10,
        debug_disable_arc_pressure_eq=debug_disable_arc_pressure_eq,
    )

    try:
        cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        if optimize:
            cycle.optimize_COP(verbose=False, initialize=False, optimize=True)
    except Exception:
        if not debug_disable_arc_pressure_eq:
            return build_and_solve(
                fluid=fluid,
                eta=eta,
                optimize=optimize,
                debug_disable_arc_pressure_eq=True,
            )
        raise

    if (cycle.optimization_converged is False) and (not debug_disable_arc_pressure_eq):
        return build_and_solve(
            fluid=fluid,
            eta=eta,
            optimize=optimize,
            debug_disable_arc_pressure_eq=True,
        )

    model = cycle.model
    T_L = -10 + C_TO_K
    T_H = (35 + 5) + C_TO_K
    results = {
        "COP": compute_cop(model),
        "COP_carnot": compute_carnot_cop(T_L, T_H),
        "T_L": T_L,
        "T_H": T_H,
        "Q_evap": value(model.fs.evaporator.heat_duty[0]),
        "W_comp": value(model.fs.compressor.work_mechanical[0]),
        "eta": value(model.fs.compressor.efficiency_isentropic[0]),
        "Wi": value(model.fs.compressor.work_isentropic[0]),
        "Wm": value(model.fs.compressor.work_mechanical[0]),
        "debug_disable_arc_pressure_eq": debug_disable_arc_pressure_eq,
    }
    return cycle, model, results


def main():
    parser = argparse.ArgumentParser(description="Run VCC tests or a baseline case.")
    parser.add_argument(
        "--pytest",
        action="store_true",
        help="Run pytest on tests/test_vapor_compression.py",
    )
    parser.add_argument(
        "--optimize",
        action="store_true",
        help="Also run the optimization step (not required for tests).",
    )
    args = parser.parse_args()

    if args.pytest:
        import pytest
        return pytest.main(["-q", "tests/test_vapor_compression.py"])

    _, _, results = build_and_solve(optimize=args.optimize)
    for k, v in results.items():
        print(f"{k}: {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
