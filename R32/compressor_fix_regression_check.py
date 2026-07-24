"""
compressor_fix_regression_check.py -- diagnose why NIST at T_amb=20 now
gives COP=3.4225 (fallback=False, converged=True per phase4_cop_NIST.csv)
instead of the PREVIOUSLY CONFIRMED COP=3.1624 (also converged=True,
"EXIT: Optimal Solution Found", clean stream table: evap out 247.33K,
cond out 304.51K, valve-to-evap 245.44K) from right before the compressor
initialization fix (_initialize_compressor_with_retry rewritten to bypass
IDAES's built-in init_isentropic()).

Both results report converged=True with no fallback needed, so this is
NOT the same kind of bug as before (dome-boundary degeneracy, composition
contradiction, etc.) -- it looks like the SAME system of equations has
multiple valid roots, and changing the compressor's warm-start path (via
the new, more robust initialization) shifted which root Newton finds.

This script re-runs the exact NIST/T_amb=20 case and prints the full
stream table PLUS actual superheat/subcool at the evaporator/condenser
outlets (should be ~0 for the intended SH=SC=0 ideal-cycle spec) and the
COP/Carnot ratio (Phase 3a's established band for this cycle: ~0.75-0.78,
see PHASE3_NOTES.md). If SH/SC have drifted away from ~0 and/or COP/Carnot
is well outside that band, this new "converged" point is a DIFFERENT,
less-intended branch, not a genuine improvement.

Extended (2026-07-24): parametrized by method (argv[1], default NIST) to
diagnose the SAME class of issue for GCGP, whose T_amb=20 case regressed
from an earlier validated -0.00% (COP=3.1899) to +7.65% (COP=3.4340) after
the compressor fixes (passes 1-3) -- a similar branch-selection anomaly to
NIST's, but NOT fixed by the same later attempt (isentropic phase_frac
fix, pass 4) that helped NIST somewhat -- pass 4 made the FULL grid worse
and was reverted. GCGP has different Tc/Pc/omega than NIST, so it lands at
a different pressure ratio/discharge Tsat for the same ambient spec --
this script lets us see exactly where GCGP's solution differs.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
import sys
from pyomo.environ import value
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

METHOD = sys.argv[1] if len(sys.argv) > 1 else "NIST"
T_EVAP_SAT = -29.0
T_AMB = 20
T_COND_SAT = T_AMB + 9


def carnot_cop(t_evap_c, t_cond_c):
    T_L = t_evap_c + 273.15
    T_H = t_cond_c + 273.15
    return T_L / (T_H - T_L)


vc = SimpleVaporCompressionCycle("R32", compressor_efficiency=0.9999,
                                  mode=Mode.IMPROVED_TPX, method=METHOD)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=T_COND_SAT)
vc.initialize(verbose=False)
vc.set_specifications(
    ambient_temperature=T_AMB, condenser_approach=9, evap_sat_temperature=-29,
    superheating=0, subcooling=0, max_pressure_ratio=10,
)
cop, converged = vc.optimize_COP(verbose=True, initialize=True, optimize=False)

m = vc.model
e_out = m.fs.evaporator.control_volume.properties_out[0]
c_out = m.fs.condenser.control_volume.properties_out[0]
comp_in = m.fs.compressor.control_volume.properties_in[0]
comp_out = m.fs.compressor.control_volume.properties_out[0]
comp_isen = m.fs.compressor.properties_isentropic[0]

# This package exposes temperature_bubble/temperature_dew (NOT
# temperature_sat, which is a Helmholtz-only attribute -- for a pure fluid
# bubble == dew == Tsat, the convention used throughout this session).
e_tsat = value(e_out.temperature_bubble["Vap", "Liq"])
c_tsat = value(c_out.temperature_bubble["Vap", "Liq"])
SH = value(e_out.temperature) - e_tsat
SC = c_tsat - value(c_out.temperature)
carnot = carnot_cop(T_EVAP_SAT, T_COND_SAT)

entropy_residual = value(comp_isen.entr_mol) - value(comp_in.entr_mol)

print(f"\n{'='*70}\n  REGRESSION CHECK: {METHOD}, T_amb=20\n{'='*70}")
print(f"COP = {cop:.4f}, converged = {converged}")
print(f"Carnot COP = {carnot:.4f}, COP/Carnot = {cop/carnot:.4f}  "
      f"(Phase 3a established band: ~0.75-0.78)")
print(f"Evaporator outlet: T={value(e_out.temperature):.3f} K, "
      f"Tsat={e_tsat:.3f} K, superheat={SH:+.3f} K")
print(f"Condenser outlet:  T={value(c_out.temperature):.3f} K, "
      f"Tsat={c_tsat:.3f} K, subcool={SC:+.3f} K")
print(f"Compressor inlet:      T={value(comp_in.temperature):.3f} K, "
      f"P={value(comp_in.pressure):.1f} Pa, h={value(comp_in.enth_mol):.2f} J/mol")
print(f"Compressor ISENTROPIC:  T={value(comp_isen.temperature):.3f} K, "
      f"P={value(comp_isen.pressure):.1f} Pa, h={value(comp_isen.enth_mol):.2f} J/mol, "
      f"entropy match residual={entropy_residual:.6e}")
print(f"Compressor REAL outlet: T={value(comp_out.temperature):.3f} K, "
      f"P={value(comp_out.pressure):.1f} Pa, h={value(comp_out.enth_mol):.2f} J/mol")
print(f"Compressor pressure ratio: "
      f"{value(comp_out.pressure)/value(comp_in.pressure):.3f}")
print(f"Efficiency: {value(m.fs.compressor.efficiency_isentropic[0]):.4f}")
print(f"Ideal (isentropic) enthalpy rise:  "
      f"{value(comp_isen.enth_mol) - value(comp_in.enth_mol):.2f} J/mol")
print(f"Actual enthalpy rise (should be ideal/efficiency): "
      f"{value(comp_out.enth_mol) - value(comp_in.enth_mol):.2f} J/mol")
print(f"Evaporator heat duty: {value(m.fs.evaporator.heat_duty[0]):.2f} W")
print(f"Compressor work: {value(m.fs.compressor.work_mechanical[0]):.2f} W")

print("\nFull stream table:")
m.fs.report()
