"""
Phase 3b flash-seeding test (isolation).

Reproduces the two-phase evaporator-inlet flash-init failure on a SINGLE state
block, then tests a single-phase-seed-via-Ambrose-Walton recipe to make the
generic (modular cubic-PR) flash converge onto the two-phase state.

Condition = evaporator inlet: R-32, T = -29 C, P = Psat(-29) (AW) -> on the dome.

Once a recipe converges here, apply the same seeding inside
vapor_compression_cubic.initialize() for the evaporator's two-phase inlet.

Author: Shilpa Narasimhan   Support: Claude AI   Date: 2026-07-20
"""
import matplotlib; matplotlib.use("Agg")
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from phase_1_cubic_eos_validation import make_config, METHODS
from phase_2_saturation_dome import ambrose_walton_psat   # validated AW Psat

p = METHODS["NIST"]
T = -29.0 + 273.15                                   # K
Psat = ambrose_walton_psat(T, p["Tc"], p["Pc"], p["omega"])   # Pa
print(f"T = {T:.2f} K, AW Psat = {Psat/1e5:.3f} bar\n")


def new_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**make_config(p))
    m.fs.state = m.fs.props.build_state_block([0], defined_state=True)
    sb = m.fs.state[0]
    sb.flow_mol.fix(1.0)
    sb.mole_frac_comp["R32"].fix(1.0)
    sb.temperature.fix(T)
    return m, sb


# ---- Attempt 1: plain init at the two-phase point (expected to fail) ----
print("=== Attempt 1: plain init at (T, Psat) [on the dome] ===")
m, sb = new_block()
sb.pressure.fix(Psat)
try:
    m.fs.state.initialize(outlvl=idaeslog.WARNING)
    print(f"  OK  phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")


# ---- Attempt 2: single-phase seed (P below Psat -> superheated vapor), then
#      move P to Psat and solve into two-phase ----
print("\n=== Attempt 2: single-phase seed below Psat, then solve up to Psat ===")
m, sb = new_block()
sb.pressure.fix(0.70 * Psat)          # below Psat at this T -> single-phase vapor
try:
    m.fs.state.initialize(outlvl=idaeslog.WARNING)
    print(f"  seed init OK (single-phase), phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
    sb.pressure.fix(Psat)             # move onto the dome
    res = get_solver().solve(m)
    print(f"  solve to Psat: {res.solver.termination_condition}")
    print(f"  phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}, "
          f"hL = {value(sb.enth_mol_phase['Liq']):.1f}, hV = {value(sb.enth_mol_phase['Vap']):.1f} J/mol")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")


# ---- Attempt 3: WELL-POSED two-phase spec. Seed single-phase, then fix quality
#      and FREE pressure (Gibbs F=1: fix T + quality, let P solve to Psat).
#      Do NOT fix both T and P on the dome (that is the degenerate case). ----
print("\n=== Attempt 3: seed single-phase, then fix quality + free P (F=1) ===")
m, sb = new_block()
sb.pressure.fix(0.5* Psat)          # single-phase vapor seed
m.fs.state.initialize(outlvl=idaeslog.WARNING)
print(f"  seed init OK, phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
sb.pressure.unfix()                    # +1 DOF
sb.phase_frac["Vap"].fix(0.2)          # -1 DOF: quality = 0.2 (well-posed with T fixed)
# assert degrees_of_freedom(m) == 0, "DOF != 0"
res = get_solver().solve(m)
print(f"  solve: {res.solver.termination_condition}")
print(f"  P = {value(sb.pressure)/1e5:.3f} bar (AW Psat = {Psat/1e5:.3f}), "
      f"phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
print(f"  hL = {value(sb.enth_mol_phase['Liq']):.1f}, hV = {value(sb.enth_mol_phase['Vap']):.1f} J/mol")


# ---- Attempt 4: seed single-phase, then PIN P = Psat (two distinct roots exist
#      there) and FREE T, with quality fixed. Freeing P (Attempt 3) let the solver
#      escape to subcooled liquid (single root -> trivial hL=hV). Pinning P keeps
#      it where two phases genuinely exist. ----
print("\n=== Attempt 4: seed single-phase, then fix P=Psat + quality, free T ===")
m, sb = new_block()
sb.pressure.fix(0.90 * Psat)
m.fs.state.initialize(outlvl=idaeslog.WARNING)
print(f"  seed init OK, phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
sb.pressure.fix(Psat)                  # pin at the REAL Psat (two roots exist)
sb.temperature.unfix()                 # +1 DOF (T solves to Tsat)
sb.phase_frac["Vap"].fix(0.2)          # -1 DOF: quality
assert degrees_of_freedom(m) == 0, "DOF != 0"
res = get_solver().solve(m)
print(f"  solve: {res.solver.termination_condition}")
print(f"  T = {value(sb.temperature)-273.15:.2f} C (want -29), P = {value(sb.pressure)/1e5:.3f} bar")
print(f"  phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
print(f"  hL = {value(sb.enth_mol_phase['Liq']):.1f}, hV = {value(sb.enth_mol_phase['Vap']):.1f} J/mol "
      f"(distinct => real two-phase)")


# ---- Attempt 5: same well-posed spec as Attempt 4 (fix P=Psat + quality, free
#      T), but ANCHOR T near the true value with a tight temporary bound so
#      Newton cannot escape to the trivial (single-root) branch. Attempts 3 & 4
#      both converged "optimal" but to the TRIVIAL solution (hL=hV, T or P far
#      from the dome) -- the trivial solution always satisfies equal-fugacity
#      once both roots collapse to one, so fixing quality alone did not prevent
#      the solver from wandering there. A tight bound on the freed variable
#      keeps Newton in the correct basin. ----
print("\n=== Attempt 5: like Attempt 4, but T bounded tightly around Tsat ===")
m, sb = new_block()
sb.pressure.fix(0.90 * Psat)
m.fs.state.initialize(outlvl=idaeslog.WARNING)
print(f"  seed init OK, phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
sb.pressure.fix(Psat)
sb.temperature.unfix()
sb.temperature.setlb(T - 5.0)          # anchor: Tsat +/- 5 K, not the full domain
sb.temperature.setub(T + 5.0)
sb.temperature.set_value(T)            # warm-start AT the expected Tsat
sb.phase_frac["Vap"].fix(0.2)
assert degrees_of_freedom(m) == 0, "DOF != 0"
res = get_solver().solve(m)
print(f"  solve: {res.solver.termination_condition}")
print(f"  T = {value(sb.temperature)-273.15:.2f} C (want -29), P = {value(sb.pressure)/1e5:.3f} bar")
print(f"  phase_frac[Vap] = {value(sb.phase_frac['Vap']):.4f}")
hL, hV = value(sb.enth_mol_phase['Liq']), value(sb.enth_mol_phase['Vap'])
print(f"  hL = {hL:.1f}, hV = {hV:.1f} J/mol, latent = {hV-hL:.1f} J/mol "
      f"(~{(hV-hL)/0.052024:.0f} J/kg) -- distinct & physical => real two-phase" if abs(hV-hL) > 1
      else f"  hL = {hL:.1f}, hV = {hV:.1f} J/mol -- STILL TRIVIAL")
