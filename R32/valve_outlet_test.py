"""
valve_outlet_test.py -- test whether giving the expansion valve's OUTLET an
explicit, correctly-placed temperature guess fixes its wrong-branch
convergence (T=283.32K, 39K above the model's own dome at this pressure,
yet phase_frac reporting ~99% liquid -- an internal contradiction: that far
above the boiling point, it should be overwhelmingly vapor).

Hypothesis: with no state_args given at all for the outlet, IDAES defaults
to copying the INLET's temperature (284.15K) as the starting guess -- the
same "default guess copies inlet unchanged" pattern already found and fixed
for the compressor. The real answer should be near what the evaporator's
inlet already converged to earlier in this same run (T=248.737K,
phase_frac[Vap]=0.064) -- since the valve's outlet and the evaporator's
inlet are literally the same physical point.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from idaes.core.util.exceptions import InitializationError
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode
import logging as pylog

vc = SimpleVaporCompressionCycle(
    "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX
)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)

# Let this run through evaporator/compressor/condenser/valve-inlet (all
# confirmed working) and fail at the valve's outlet, as we already know it
# does. We don't care about that failed attempt -- we're about to redo the
# outlet step ourselves with a better guess.
try:
    vc.initialize(verbose=False)
    print("vc.initialize() completed without raising (unexpected at this point)\n")
except InitializationError as e:
    print(f"vc.initialize() failed as expected: {e}\n")

m = vc.model
evap_in = m.fs.evaporator.control_volume.properties_in[0]
valve_in = m.fs.expansion_valve.control_volume.properties_in[0]
valve_out = m.fs.expansion_valve.control_volume.properties_out[0]

T_evap_in = value(evap_in.temperature)
print(f"Evaporator inlet (already solved earlier this run): T = {T_evap_in:.3f} K, "
      f"phase_frac[Vap] = {value(evap_in.phase_frac['Vap']):.6f}")

print("\n=== Before the retry ===")
print(f"  valve outlet T = {value(valve_out.temperature):.3f} K "
      f"(tbub = {value(valve_out.temperature_bubble['Vap','Liq']):.3f} K)")
print(f"  valve outlet phase_frac[Vap] = {value(valve_out.phase_frac['Vap']):.6f}")

# --- The fix under test: explicit state_args using the evaporator inlet's
# already-known-good temperature, instead of letting the outlet default to
# copying the valve's own inlet temperature (284.15K, wrong direction). ---
print(f"\nRetrying expansion_valve.initialize() with an explicit outlet guess "
      f"(T={T_evap_in:.3f} K, from the evaporator inlet)...")
try:
    m.fs.expansion_valve.initialize(
        outlvl=pylog.WARNING,
        state_args={
            "flow_mol": 1.0,
            "temperature": T_evap_in,
            "pressure": value(m.fs.expansion_valve.outlet.pressure[0]),
            "mole_frac_comp": {"R32": 1.0},
        },
        optarg={"tol": 1e-4, "constr_viol_tol": 1e-4, "acceptable_tol": 1e-3},
    )
    print("expansion_valve.initialize() SUCCEEDED with explicit outlet guess")
except Exception as e:
    print(f"expansion_valve.initialize() still failed: {type(e).__name__}: {e}")

print("\n=== After the retry ===")
T_out = value(valve_out.temperature)
tbub_out = value(valve_out.temperature_bubble["Vap", "Liq"])
vap_frac_out = value(valve_out.phase_frac["Vap"])
print(f"  valve outlet T = {T_out:.3f} K, tbub = {tbub_out:.3f} K "
      f"(T is {'BELOW' if T_out < tbub_out else 'ABOVE'} tbub)")
print(f"  valve outlet phase_frac[Vap] = {vap_frac_out:.6f}")
print(f"  valve outlet entr_mol = {value(valve_out.entr_mol):.4f} J/mol/K")

print("\n=== Verdict ===")
consistent = (T_out < tbub_out and vap_frac_out < 0.5) or (T_out > tbub_out and vap_frac_out > 0.5)
if consistent:
    print("  CONSISTENT -- T-vs-tbub relationship now agrees with phase_frac.")
    print("  (Compare against the evaporator inlet's own T/phase_frac above --")
    print("   the two should now be close, since they're the same physical point.)")
else:
    print("  STILL INCONSISTENT -- T-vs-tbub and phase_frac disagree, same")
    print("  contradiction as before. Need a different fix (e.g. a bound on")
    print("  phase_frac here too, or freeing/bounding T explicitly).")
