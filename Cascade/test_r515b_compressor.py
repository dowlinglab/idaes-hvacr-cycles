import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent/ "R515B_idaes_package"))


from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Compressor
from idaes.core.solvers import get_solver
from r515b_property_package import R515BParameterBlock

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic = False)
m.fs.props = R515BParameterBlock()
m.fs.comp = Compressor(property_package = m.fs.props)

# --- Fix inlet, outlet pressure, efficiency (values from the cascade run) ---
m.fs.comp.inlet.flow_mass[0].fix(1.0)        # kg/s
m.fs.comp.inlet.pressure[0].fix(5.0296e5)    # Pa
m.fs.comp.inlet.enth_mass[0].fix(3.9784e5)   # J/kg, superheated vapor
m.fs.comp.outlet.pressure[0].fix(7.5436e5)   # Pa
m.fs.comp.efficiency_isentropic[0].fix(0.80)

m.fs.comp.initialize(optarg={"tol": 1e-4})
res = get_solver().solve(m)
print("Solver:", res.solver.termination_condition)


# --- Check: outlet h should be ABOVE inlet h, and vapor (x = 1) ---
s_in = m.fs.comp.control_volume.properties_in[0]
s_out = m.fs.comp.control_volume.properties_out[0]
s_isen = m.fs.comp.properties_isentropic[0]
print("h  in / out / isentropic [kJ/kg]:",
      value(s_in.enth_mass)/1e3, value(s_out.enth_mass)/1e3, value(s_isen.enth_mass)/1e3)
print("x  in / out / isentropic:",
      value(s_in.vapor_frac), value(s_out.vapor_frac), value(s_isen.vapor_frac))
print("Work [kJ/kg]:", value(m.fs.comp.work_mechanical[0])/1e3)