import os
import sys
from pyomo.environ import value

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
import diagnose_45C as d

def build(fluid, ambient = 45.0):
    cycle, _ = d._build_cycle(fluid)
    cycle.specify_initial_conditions(hot_ambient_C= ambient, cold_evap_C= d.EVAP_SAT_TARGET_C)
    cycle.initialize(verbose = False)
    cold = dict(evap_sat_temperature = d.EVAP_SAT_TARGET_C,
                superheating = d.SUPERHEAT_C, superheating_max = d.SUPERHEAT_MAX_C,
                subcooling = d.SUBCOOL_C, subcooling_max = d.SUBCOOL_MAX_C,
                max_pressure_ratio = d.MAX_PRESSURE_RATIO)
    hot = dict(ambient_temperature = ambient, condenser_approach = d.COND_APPROACH_C,
               superheating = d.SUPERHEAT_C, superheating_max = d.SUPERHEAT_MAX_C,
               subcooling = d.SUBCOOL_C, subcooling_max = d.SUBCOOL_MAX_C,
               max_pressure_ratio = d.MAX_PRESSURE_RATIO)

    cycle.set_specifications(hot = hot, cold = cold,
                             cascade_approach_dT= d.CASCADE_APPROACH_DT_C,
                             plr = d.PLR_VALUE, cd = d.CD_VALUE)

    return cycle

cycle = build("R1234yf")
cop, ok = cycle.optimize_COP(verbose = True, initialize= True, optimize = True)
print("COP:", cop, "ok:", ok)