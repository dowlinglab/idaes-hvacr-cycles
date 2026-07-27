"""
compressor_fix_regression_check.py -- detailed per-cell diagnostic for the
R32 cubic-PR cycle: full stream table, evap/cond superheat-subcool, and
the compressor's real-outlet vs isentropic T/h side by side, plus COP and
COP/Carnot.

Extended (2026-07-27): now loops over the WHOLE ambient sweep (10/15/20/25
C) in one run instead of a single hardcoded T_amb=20, per the request to
"step back and print the entire COP vs ambient loop" rather than keep
zooming in on one point. Still parametrized by method (argv[1], default
GCGP -- this is the one under active investigation, task #37). Pass a
second argument to check a single ambient instead of the full sweep, e.g.
`python3 compressor_fix_regression_check.py GCGP 20`.

Original purpose (2026-07-24): diagnose why NIST/GCGP's compressor
initialization rewrite shifted which root Newton finds at T_amb=20 --
COP/Carnot's established band (Phase 3a, PHASE3_NOTES.md: ~0.75-0.78) and
SH/SC near 0 are the tells for whether a "converged=True" result is on the
intended branch or a different, less-intended one. The isentropic-vs-real-
outlet T/h relationship (should track together for efficiency~1) is the
tell for the specific "wrong root near the two-phase boundary" bug this
whole session has been chasing.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
import sys
from pyomo.environ import value
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

METHOD = sys.argv[1] if len(sys.argv) > 1 else "GCGP"
AMBIENTS = [int(sys.argv[2])] if len(sys.argv) > 2 else [10, 15, 20, 25]
T_EVAP_SAT = -29.0


def carnot_cop(t_evap_c, t_cond_c):
    T_L = t_evap_c + 273.15
    T_H = t_cond_c + 273.15
    return T_L / (T_H - T_L)


def run_one(method, T_amb):
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
    e_out = m.fs.evaporator.control_volume.properties_out[0]
    c_out = m.fs.condenser.control_volume.properties_out[0]
    comp_in = m.fs.compressor.control_volume.properties_in[0]
    comp_out = m.fs.compressor.control_volume.properties_out[0]
    comp_isen = m.fs.compressor.properties_isentropic[0]

    e_tsat = value(e_out.temperature_bubble["Vap", "Liq"])
    c_tsat = value(c_out.temperature_bubble["Vap", "Liq"])
    SH = value(e_out.temperature) - e_tsat
    SC = c_tsat - value(c_out.temperature)
    carnot = carnot_cop(T_EVAP_SAT, T_cond_sat)
    entropy_residual = value(comp_isen.entr_mol) - value(comp_in.entr_mol)

    print(f"\n{'='*70}\n  {method}, T_amb={T_amb} C (T_cond_sat={T_cond_sat} C)\n{'='*70}")
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
    print(f"  T_out - T_isen = {value(comp_out.temperature)-value(comp_isen.temperature):+.3f} K, "
          f"h_out - h_isen = {value(comp_out.enth_mol)-value(comp_isen.enth_mol):+.3f} J/mol")
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

    return {"T_amb": T_amb, "cop": cop, "converged": converged,
            "T_isen": value(comp_isen.temperature), "T_out": value(comp_out.temperature)}


rows = [run_one(METHOD, T_amb) for T_amb in AMBIENTS]

print(f"\n{'='*70}\n  SUMMARY: {METHOD}\n{'='*70}")
print(f"{'T_amb':>7}{'COP':>9}{'T_isen':>10}{'T_out':>10}{'T_out-T_isen':>14}")
for r in rows:
    print(f"{r['T_amb']:>7}{r['cop']:>9.4f}{r['T_isen']:>10.3f}{r['T_out']:>10.3f}"
          f"{r['T_out']-r['T_isen']:>14.3f}")
