"""Regression tests for PLR compressor flow scaling.

Run:
  pytest -q
  pytest -q -m slow
"""

import pytest
import warnings
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

idaes = pytest.importorskip("idaes")
pyomo = pytest.importorskip("pyomo")

from idaes.core.solvers import get_solver
from vapor_compression_plr import SimpleVaporCompressionCyclePLR, Mode

# Suppress matplotlib non-interactive show warnings
warnings.filterwarnings(
    "ignore",
    message="FigureCanvasAgg is non-interactive*",
    category=UserWarning,
)


def _skip_if_solver_unavailable():
    solver = get_solver()
    if solver is None or not solver.available(exception_flag=False):
        pytest.skip("IDAES solver is not available via get_solver().")
    return solver


def _design_kwargs():
    return dict(
        ambient_temperature=35,
        condenser_approach=5,
        evap_sat_temperature=-10,
        superheating=3,
        subcooling=3,
        max_pressure_ratio=4,
    )


@pytest.mark.slow
def test_plr_flow_scaling():
    _skip_if_solver_unavailable()
    plt.show = lambda *args, **kwargs: None

    cycle = SimpleVaporCompressionCyclePLR(
        "R134a", compressor_efficiency=0.75, mode=Mode.IMPROVED_TPX, plr=0.75
    )
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    cycle.set_specifications(**_design_kwargs(), plr=0.5)

    # Full-load flow is fixed to 1 kg/s, all loop inlet flows should be PLR*1
    assert abs(cycle.model.fs.evaporator.inlet.flow_mass[0].value - 0.5) < 1e-8
    assert abs(cycle.model.fs.compressor.inlet.flow_mass[0].value - 0.5) < 1e-8
    assert abs(cycle.model.fs.condenser.inlet.flow_mass[0].value - 0.5) < 1e-8
    assert abs(cycle.model.fs.expansion_valve.inlet.flow_mass[0].value - 0.5) < 1e-8
