"""
compressor_branch_timing_debug.py -- pass 4/5/6 all tried to PATCH the GCGP
T_amb=20 branch problem with a bound/constraint tweak, and all three made
things worse elsewhere (see BREADCRUMB_07-20.md, 2026-07-24 pass 4/5/6
entries). Before trying a pass 7 patch, this script instead answers a more
basic question: WHEN does the compressor real outlet drift onto the bad
(near-Tsat) branch -- is it already wrong right after `vc.initialize()`
(a warm-start problem), or does it only go wrong once `set_specifications()`
and the real coupled solve (`optimize_COP()`) run (a solve-path problem)?

`_initialize_compressor_with_retry()` computes a deliberately good warm
start for the real outlet (anchored to the isentropic estimate via the
efficiency ratio -- see its Step 3). If that good value is still in place
right after `vc.initialize()` but drifts away during/after
`set_specifications()`+`optimize_COP()`, the bug is in the COUPLED solve
(something in the rest of the flowsheet is pulling the compressor outlet
off its good warm start) -- not in initialization itself, and not
something a compressor-local bound can necessarily fix in isolation.

Prints compressor real-outlet and isentropic T/h at three checkpoints:
1. Right after vc.initialize()
2. Right after vc.set_specifications() (bounds/constraints applied, but
   the real coupled solve inside optimize_COP() has NOT run yet)
3. Right after optimize_COP()'s coupled solve

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-24
"""
import sys
from pyomo.environ import value
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

METHOD = sys.argv[1] if len(sys.argv) > 1 else "GCGP"
T_AMB = 20
T_COND_SAT = T_AMB + 9


def snapshot(label, m):
    comp_out = m.fs.compressor.control_volume.properties_out[0]
    comp_isen = m.fs.compressor.properties_isentropic[0]
    print(f"\n--- {label} ---")
    print(f"  Real outlet:  T={value(comp_out.temperature):.3f} K, "
          f"h={value(comp_out.enth_mol):.2f} J/mol, "
          f"lb={comp_out.temperature.lb}, ub={comp_out.temperature.ub}")
    print(f"  Isentropic:   T={value(comp_isen.temperature):.3f} K, "
          f"h={value(comp_isen.enth_mol):.2f} J/mol, "
          f"lb={comp_isen.temperature.lb}, ub={comp_isen.temperature.ub}")
    print(f"  T_out - T_isen = {value(comp_out.temperature) - value(comp_isen.temperature):+.3f} K, "
          f"h_out - h_isen = {value(comp_out.enth_mol) - value(comp_isen.enth_mol):+.3f} J/mol")


vc = SimpleVaporCompressionCycle("R32", compressor_efficiency=0.9999,
                                  mode=Mode.IMPROVED_TPX, method=METHOD)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=T_COND_SAT)
vc.initialize(verbose=False)
snapshot("CHECKPOINT 1: right after vc.initialize()", vc.model)

vc.set_specifications(
    ambient_temperature=T_AMB, condenser_approach=9, evap_sat_temperature=-29,
    superheating=0, subcooling=0, max_pressure_ratio=10,
)
snapshot("CHECKPOINT 2: right after set_specifications() (pre coupled solve)", vc.model)

cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
snapshot("CHECKPOINT 3: right after optimize_COP() (post coupled solve)", vc.model)

print(f"\nCOP = {cop:.4f}, converged = {converged}")
