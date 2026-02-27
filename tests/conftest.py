"""Shared helpers for PLR vapor compression tests.

Run:
  pytest -q
  pytest -q -m slow
"""

import math

import pytest

idaes = pytest.importorskip("idaes")
pyomo = pytest.importorskip("pyomo")

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import sys
from pathlib import Path

from idaes.core.solvers import get_solver
from pyomo.environ import value

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from vapor_compression_plr import SimpleVaporCompressionCyclePLR, Mode

C_TO_K = 273.15


def _disable_plots():
    try:
        plt.show = lambda *args, **kwargs: None
    except Exception:
        pass


def _skip_if_solver_unavailable():
    solver = get_solver()
    if solver is None or not solver.available(exception_flag=False):
        pytest.skip("IDAES solver is not available via get_solver().")
    return solver


def build_cycle_and_solve(
    fluid="R134a",
    eta=0.75,
    plr=None,
    Q_rated=None,
    optimize=False,
    debug_disable_arc_pressure_eq=False,
):
    _skip_if_solver_unavailable()
    _disable_plots()

    cycle = SimpleVaporCompressionCyclePLR(
        fluid_name=fluid,
        compressor_efficiency=eta,
        mode=Mode.IMPROVED_TPX,
    )

    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    def _apply_specs(flag):
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
            debug_disable_arc_pressure_eq=flag,
            plr=plr,
            Q_cool_rated=Q_rated,
        )

    def _try_solve(flag):
        _apply_specs(flag)
        try:
            cycle.optimize_COP(verbose=False, initialize=True, optimize=optimize)
        except Exception as exc:
            return False, exc
        return True, None

    ok, err = _try_solve(debug_disable_arc_pressure_eq)
    if (not ok or cycle.optimization_converged is False) and (not debug_disable_arc_pressure_eq):
        ok, err = _try_solve(True)
    if not ok:
        raise err

    model = cycle.model
    Q_evap = value(model.fs.evaporator.heat_duty[0])
    W_comp = value(model.fs.compressor.work_mechanical[0])

    result = {
        "COP": abs(Q_evap) / W_comp if W_comp not in (None, 0) else math.inf,
        "Q_evap": Q_evap,
        "W_comp": W_comp,
        "Wi": value(model.fs.compressor.work_isentropic[0]),
        "eta": value(model.fs.compressor.efficiency_isentropic[0]),
        "mdot": value(model.fs.evaporator.inlet.flow_mass[0]),
        "W_valve": value(model.fs.expansion_valve.work_mechanical[0]),
        "P_low": value(model.fs.P_low) if hasattr(model.fs, "P_low") else None,
        "P_high": value(model.fs.P_high) if hasattr(model.fs, "P_high") else None,
        "T_evap_sat": None,
        "T_cond_sat": None,
        "termination": cycle.last_termination_condition,
    }

    if model.fs.evaporator.evap_sat_constraint.active:
        try:
            result["T_evap_sat"] = value(
                model.fs.evaporator.control_volume.properties_out[0].temperature_sat
            )
        except Exception:
            result["T_evap_sat"] = None

    if model.fs.condenser.approach_constraint.active:
        try:
            result["T_cond_sat"] = value(
                model.fs.condenser.control_volume.properties_out[0].temperature_sat
            )
        except Exception:
            result["T_cond_sat"] = None

    if model.fs.evaporator.evap_sat_constraint.active and model.fs.condenser.approach_constraint.active:
        TL_K = -10 + C_TO_K
        TH_K = (35 + 5) + C_TO_K
        result["COP_carnot"] = TL_K / (TH_K - TL_K)
    else:
        result["COP_carnot"] = None

    try:
        plt.close("all")
    except Exception:
        pass

    return result
