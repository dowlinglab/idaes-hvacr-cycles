"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/10/2026

What this does
    Checks the R515B compressor at the 30 C cascade point, three ways:
      Part 1  The right answer, from the validated R515B core (plain Python).
      Part 2  One IDAES R515B state block: P_out + inlet entropy -> where
              does it land?
      Part 3  The IDAES compressor, started from Part 1's values instead
              of IDAES's own initialize().

Why
    In the cascade the R515B compressor outlet came out impossible
    (h dropped 398 -> 292 kJ/kg, two-phase). This splits the compressor's
    main calculation ("same entropy, higher pressure") into pieces to find
    which piece is wrong.

Result (10/10, sandbox)
    Part 1: h_isen 405.35, h_out 407.23 kJ/kg, vapor, work +9.39 kJ/kg.
    Part 2: state block lands on 405.35 kJ/kg, vapor -- correct.
    Part 3: compressor solves correctly (407.23 / 405.35 / +9.39).
    IDAES's own compressor initialize() fails (isentropic vapor_frac
    stuck at 1.01). So the bug is the starting guess, not the physics.

Run
    python check_r515b_isentrope.py 2>&1 | grep -v "DEBUG\\|Missing scaling" | tail -12
"""

import sys
from pathlib import Path

# R515B package + its validated mixture model live at the repo root
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "R515B_idaes_package"))
sys.path.insert(0, str(REPO / "R515B_props_validated"))

from scipy.optimize import root
from pyomo.environ import ConcreteModel, Constraint, value, units as pyunits
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Compressor
import r515b_helmholtz_core as core
from r515b_property_package import R515BParameterBlock

# Cascade 30 C point (hot loop, from the stream table)
P_IN, H_IN, P_OUT, ETA = 5.0296e5, 3.9784e5, 7.5436e5, 0.80   # Pa, J/kg, Pa, -


# --- Part 1: right answer from the validated core --------------------------
d1 = core.load_idaes_helmholtz_json(core.FLUID1)
d2 = core.load_idaes_helmholtz_json(core.FLUID2)
x1 = core.r515b_x1()
mw = x1 * core.mw_from_json(d1) + (1 - x1) * core.mw_from_json(d2)    # kg/mol

inlet = core.solve_singlephase_tr_at_ph_pseudopure(d1, d2, x1, P_IN, H_IN * mw, 301.0, 27.05 / mw)
s_in = core.mix_entropy_direct(d1, d2, inlet["T_K"], inlet["rho_molm3"], x1)        # J/mol/K

sol = root(core._isentrope_2eq_residual, x0=[315.0, P_OUT / (8.314 * 315.0)],
           args=(d1, d2, x1, P_OUT, s_in), method="hybr", tol=1e-12)
h_isen = core.mix_state(d1, d2, sol.x[0], sol.x[1], x1).h_jmol / mw                 # J/kg
h_out = H_IN + (h_isen - H_IN) / ETA
t_sat = core.solve_t_sat_at_p_pseudopure(d1, d2, x1, P_OUT, 312.0, 9000.0, 330.0)["T_K"]

print("PART 1 (core):  h_isen = %.2f  h_out = %.2f kJ/kg  work = %.2f kJ/kg  T_isen = %.2f K  T_sat = %.2f K"
      % (h_isen / 1e3, h_out / 1e3, (h_out - H_IN) / 1e3, sol.x[0], t_sat))


# --- Part 2: one IDAES state block at P_out with the inlet entropy -----------
m = ConcreteModel(); m.fs = FlowsheetBlock(dynamic=False)
m.fs.props = R515BParameterBlock()
m.fs.st = m.fs.props.build_state_block([0], defined_state=True)
st = m.fs.st[0]
st.flow_mass.fix(1.0); st.pressure.fix(P_OUT); st.enth_mass.fix(3.0e5)    # start inside the dome on purpose
m.fs.st.initialize()
st.enth_mass.unfix()
m.fs.s_spec = Constraint(expr=st.entr_mol == s_in * pyunits.J / pyunits.mol / pyunits.K)
r = get_solver().solve(m).solver.termination_condition
print("PART 2 (state): %s  h = %.2f kJ/kg  x = %.3f" % (r, value(st.enth_mass) / 1e3, value(st.vapor_frac)))


# --- Part 3: IDAES compressor, started from Part 1's values -----------------
m = ConcreteModel(); m.fs = FlowsheetBlock(dynamic=False)
m.fs.props = R515BParameterBlock()
c = m.fs.comp = Compressor(property_package=m.fs.props)
c.inlet.flow_mass[0].fix(1.0); c.inlet.pressure[0].fix(P_IN); c.inlet.enth_mass[0].fix(H_IN)
c.outlet.pressure[0].fix(P_OUT); c.efficiency_isentropic[0].fix(ETA)

# seed each state block at the right (P, h) instead of calling c.initialize()
for blk, P, h in ((c.control_volume.properties_in, P_IN, H_IN),
                  (c.properties_isentropic, P_OUT, h_isen),
                  (c.control_volume.properties_out, P_OUT, h_out)):
    s = blk[0]
    was_fixed = [(v, v.fixed) for v in (s.flow_mass, s.pressure, s.enth_mass)]
    s.flow_mass.fix(1.0); s.pressure.fix(P); s.enth_mass.fix(h)
    blk.initialize()
    for v, f in was_fixed:
        if not f:
            v.unfix()
c.work_mechanical[0].set_value(h_out - H_IN)
c.deltaP[0].set_value(P_OUT - P_IN)

r = get_solver().solve(m).solver.termination_condition
o, i = c.control_volume.properties_out[0], c.properties_isentropic[0]
print("PART 3 (comp):  %s  h_out = %.2f  h_isen = %.2f kJ/kg  x_out = %.3f  work = %.2f kJ/kg"
      % (r, value(o.enth_mass) / 1e3, value(i.enth_mass) / 1e3, value(o.vapor_frac),
         value(c.work_mechanical[0]) / 1e3))
