"""
entropy_dome_proximity_check.py -- directly test the "dome proximity"
theory for the GCGP T_amb=20 compressor anomaly (task #37): did the
isentropic block's own flash, in the ALREADY-CONVERGED main solve, land
on a genuine single-phase-vapor point, or did it land inside the
two-phase region (a mixed vapor/liquid root instead of the intended
superheated-vapor one)?

REWRITE (2026-07-27): the original version of this script tried to
force-fix the isentropic block to (T=T_dew, P=P_high) and re-solve it,
using the on-demand `temperature_dew` property, to read off a dew-point
entropy s_dew and compare it to s_in. That approach turned out to be
fundamentally unreliable for two compounding reasons, discovered by
actually running it:
  1. Accessing `comp_isen.temperature_dew[...]` for the FIRST time AFTER
     the main solve triggers IDAES to construct a brand-new Var+Constraint
     pair on the fly (on-demand property) -- it does NOT automatically
     solve them. The printed "T_dew=307.150K" was consistently the SAME
     unsolved default value, not a converged physical answer.
  2. Force-fixing (T, P, flow, comp) and re-solving that freshly
     constructed sub-piece repeatedly hit "Too few degrees of freedom
     (rethrown)!" -- structural degeneracy right at the phase boundary
     that fixing/unfixing individual variables couldn't cleanly resolve
     (confirmed by phase_frac_Vap printing exactly 0.5000 -- a generic
     default, meaning that solve never actually moved).

SECOND REWRITE (2026-07-27): tried reading `comp_isen.phase_frac["Vap"]`
straight off the converged solve (no re-fixing/re-solving) as the
"trustworthy" signal. Still unreliable: it printed exactly 0.5000 even
for NIST T_amb=25 -- a fully validated, known-good case whose COP
(2.8005) matches the trusted baseline exactly. A variable frozen at a
generic default even in the GOOD case cannot be trusted as a live
indicator for the ANOMALOUS case either -- this isentropic StateBlock's
`phase_frac` var evidently isn't the thing actually driving/reflecting
the flash outcome here (possibly inert/unlinked in this particular
construction). Abandoned.

**Current approach**: compare the ALREADY-SOLVED `entr_mol` (this
block's total molar entropy, which by the isentropic constraint equals
the compressor inlet's own entropy) against the ALREADY-SOLVED
`entr_mol_phase["Vap"]` (the entropy the EoS computes for a pure-vapor
state at this same converged T, P) -- both are properties of the SAME
already-converged block, read with zero additional fixing or solving.
  - entr_mol ~= entr_mol_phase["Vap"]: the solved state IS (for all
    practical purposes) pure vapor -- clean branch, regardless of what
    `phase_frac` happened to read.
  - entr_mol measurably BELOW entr_mol_phase["Vap"] (and if available,
    closer to entr_mol_phase["Liq"]): direct, first-hand proof the
    solved state is a genuine two-phase blend -- the wrong-branch
    symptom, confirmed with no re-solve at all.
Also reports T_isen and, purely as approximate context (NOT this
method's own dew point -- that anchor is known-invalid for GCGP/SPGP per
the pass-6 finding), the real-R32 CoolProp saturated-vapor temperature
at P_high.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-27
"""
import sys
import CoolProp.CoolProp as CP
from pyomo.environ import value
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

METHODS = sys.argv[1].split(",") if len(sys.argv) > 1 else ["GCGP", "NIST"]
AMBIENTS = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [10, 15, 20, 25]


def check_one(method, T_amb):
    T_cond_sat = T_amb + 9
    vc = SimpleVaporCompressionCycle("R32", compressor_efficiency=0.9999,
                                      mode=Mode.IMPROVED_TPX, method=method)
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=T_cond_sat)
    vc.initialize(verbose=False)
    vc.set_specifications(
        ambient_temperature=T_amb, condenser_approach=9, evap_sat_temperature=-29,
        superheating=0, subcooling=0, max_pressure_ratio=10,
    )
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)

    m = vc.model
    comp_isen = m.fs.compressor.properties_isentropic[0]

    s_in = value(comp_isen.entr_mol)  # == compressor inlet entropy, by the isentropic constraint
    T_isen = value(comp_isen.temperature)
    P_high = value(comp_isen.pressure)
    s_vap = value(comp_isen.entr_mol_phase["Vap"])
    try:
        s_liq = value(comp_isen.entr_mol_phase["Liq"])
    except Exception:
        s_liq = float("nan")

    T_sat_high = CP.PropsSI('T', 'P', P_high, 'Q', 1, vc.fluid_name)  # real-R32 CoolProp reference only
    gap_vs_vap = s_in - s_vap  # ~0 => pure vapor; negative & sizable => blend

    print(f"{method:>6}  T_amb={T_amb:>3}  converged={str(converged):>5}  COP={cop:>9.4f}  "
          f"s_in={s_in:>10.4f}  s_vap={s_vap:>10.4f}  s_liq={s_liq:>10.4f}  "
          f"s_in-s_vap={gap_vs_vap:>+9.4f}  T_isen={T_isen:>8.3f}K  "
          f"T_sat_high(CoolProp,real-R32)={T_sat_high:>8.3f}K  "
          f"T_isen-T_sat_high={T_isen - T_sat_high:>+8.3f}K")
    return gap_vs_vap


print(f"{'method':>6}  {'T_amb':>8}  {'converged':>10}  {'COP':>10}  {'s_in':>10}  "
      f"{'s_vap':>10}  {'s_liq':>10}  {'s_in-s_vap':>12}  {'T_isen':>10}  "
      f"{'T_sat_high':>20}  {'T_isen-T_sat':>15}")
for method in METHODS:
    for T_amb in AMBIENTS:
        try:
            check_one(method, T_amb)
        except Exception as e:
            print(f"{method:>6}  T_amb={T_amb:>3}  CRASHED: {type(e).__name__}: {e}")
