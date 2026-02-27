import csv
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pytest
from idaes.core.solvers import get_solver
from pyomo.environ import ConcreteModel, Objective, Var

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import simple_vcc_codex as vcc


def _skip_if_solver_unavailable():
    solver = get_solver()
    if solver is None or not solver.available(exception_flag=False):
        pytest.skip("IDAES solver is not available via get_solver().")
    return solver


def _skip_if_ma57_unavailable():
    solver = _skip_if_solver_unavailable()
    solver_name = getattr(solver, "name", "").lower()
    if "ipopt" not in solver_name:
        return

    test_model = ConcreteModel()
    test_model.x = Var(initialize=1.0)
    test_model.obj = Objective(expr=(test_model.x - 1.0) ** 2)

    try:
        solver.options = {"linear_solver": "ma57"}
        result = solver.solve(test_model, tee=False)
    except Exception as exc:  # pragma: no cover - depends on solver install
        message = str(exc).lower()
        if "ma57" in message or "hsl" in message:
            pytest.skip(f"MA57 linear solver not available: {exc}")
        raise

    term = getattr(result.solver, "termination_condition", None)
    if term is None:
        return
    term_str = str(term).lower()
    if term_str not in {"optimal", "locallyoptimal", "feasible", "locally optimal"}:
        if "ma57" in term_str or "hsl" in term_str:
            pytest.skip(f"MA57 linear solver not available: {term}")


def test_cop_vs_ambient_r134a(tmp_path):
    _skip_if_ma57_unavailable()

    vcc.plt.show = lambda *args, **kwargs: None

    cycle = vcc.SimpleVaporCompressionCycle(
        fluid_name="R134a",
        compressor_efficiency=0.75,
        mode=vcc.Mode.IMPROVED_TPX,
    )
    cycle.specify_initial_conditions()
    cycle.initialize(verbose=False)

    ambient_temps = np.linspace(20, 30, 21)
    cops = []

    for ambient_c in ambient_temps:
        cycle.set_specifications(
            low_side_pressure=(20, 300),  # kPa
            high_side_pressure=(500, 2000),  # kPa
            evaporator_temperature=(-20, 0),  # degC
            condenser_temperature=(float(ambient_c) + 5, float(ambient_c) + 15),  # degC
            subcooling=3,
            superheating=3,
        )

        cop, converged = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)

        if not converged or not math.isfinite(cop):
            term = cycle.last_termination_condition
            pytest.fail(
                "COP optimization failed at ambient={}C; termination_condition={}".format(
                    ambient_c, term
                )
            )

        cops.append(float(cop))

    csv_path = tmp_path / "cop_vs_ambient_r134a.csv"
    repo_csv_path = REPO_ROOT / "cop_vs_ambient_r134a.csv"
    with csv_path.open("w", newline="") as handle, repo_csv_path.open(
        "w", newline=""
    ) as repo_handle:
        writer = csv.writer(handle)
        repo_writer = csv.writer(repo_handle)
        writer.writerow(["ambient_C", "cop"])
        repo_writer.writerow(["ambient_C", "cop"])
        for ambient_c, cop in zip(ambient_temps, cops):
            row = [float(ambient_c), float(cop)]
            writer.writerow(row)
            repo_writer.writerow(row)

    fig, ax = plt.subplots()
    ax.plot(ambient_temps, cops, marker="o")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP")
    ax.set_title("COP vs Ambient Temperature (R134a)")
    fig.tight_layout()
    fig.savefig(tmp_path / "cop_vs_ambient_r134a.png")
    fig.savefig(REPO_ROOT / "cop_vs_ambient_r134a.png")
    plt.close(fig)

    tol = 1e-2
    diffs = np.diff(cops)
    assert np.all(diffs <= tol), "COP should be non-increasing with ambient temperature"

    for ambient_c, cop in zip(ambient_temps, cops):
        assert 0.5 < cop < 20, "COP out of expected range at ambient={}C".format(ambient_c)
