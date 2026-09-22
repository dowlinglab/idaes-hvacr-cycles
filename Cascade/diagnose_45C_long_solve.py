import os
import sys
import re

from idaes.core.solvers import get_solver
from pyomo.environ import value

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT_DIR = os.path.dirname(_THIS_DIR)
if _REPO_ROOT_DIR not in sys.path:
    sys.path.insert(0, _REPO_ROOT_DIR)

from vapor_compression_cascade import CascadeCycle


AMBIENT_C = 45.0
EVAP_SAT_TARGET_C = -29.0
COND_APPROACH_C = 9.0
SUPERHEAT_C = 3.0
SUBCOOL_C = 3.0
SUPERHEAT_MAX_C = 8.0
SUBCOOL_MAX_C = 8.0
COMPRESSOR_EFFICIENCY = 0.8
MAX_PRESSURE_RATIO = 20.0
CASCADE_APPROACH_DT_C = 3.0
PLR_VALUE = 0.75
CD_VALUE = 0.13

MAX_ITER_CAP = 50000

def build_r1234yf_cycle():
    cold_fluid = {
        "gh_component": "co2", "coolprop_name": "CO2", "custom_json": None,
        "efficiency": "COMPRESSOR_EFFICIENCY",
        "superheat_max": SUPERHEAT_MAX_C, "subcool_max": SUBCOOL_MAX_C,
        "Tmax_C": None,
    }

    hot_fluid ={}