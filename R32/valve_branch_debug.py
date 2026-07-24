"""
valve_branch_debug.py -- targeted diagnostic for the expansion valve inlet's
apparent wrong-branch convergence (phase_frac[Vap] ~1.0 when it should be
~0.0, entropy mismatched from the condenser outlet feeding it).

Idea: propagate_state() only copies the four canonical state variables
(flow, T, P, composition) between connected blocks -- NOT derived variables
like phase_frac. Each block builds phase_frac from scratch via IDAES's own
FTPx.state_initialization() function, which uses a simple rule: if the fixed
temperature is below the bubble-point temperature at this pressure, guess
"almost all liquid" (vap_frac ~ 1e-5). Given our T=284.15K is comfortably
below tbub=301.856K here, that rule SHOULD produce a correct starting guess.

This script runs the real (already-fixed) vc.initialize(), then calls
state_initialization() again, by itself, directly on the valve's inlet
block -- not to re-solve anything, just to see what guess that function
computes for THIS exact state. Comparing that to what the real solver
actually converged to tells us whether:
  (a) the guess itself is correct, and the numerical solve is what moves it
      to the wrong branch afterward, or
  (b) the guess itself is already wrong, meaning the problem is upstream of
      the solve entirely.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 2026-07-23
"""

from pyomo.environ import value
from idaes.models.properties.modular_properties.state_definitions.FTPx import (
    state_initialization,
)
from vapor_compression_cubic import SimpleVaporCompressionCycle, Mode

vc = SimpleVaporCompressionCycle(
    "R32", compressor_efficiency=0.9999, mode=Mode.IMPROVED_TPX
)
vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=29)

try:
    vc.initialize(verbose=False)
    print("vc.initialize() completed without raising\n")
except Exception as e:
    print(f"vc.initialize() failed: {type(e).__name__}: {e}\n")

valve_in = vc.model.fs.expansion_valve.control_volume.properties_in[0]

T = value(valve_in.temperature)
tbub = value(valve_in.temperature_bubble["Vap", "Liq"])
vap_frac_actual = value(valve_in.phase_frac["Vap"])
entr_actual = value(valve_in.entr_mol)

print("=== What the real solver actually converged to ===")
print(f"  T = {T:.3f} K, tbub = {tbub:.3f} K  (T is {'BELOW' if T < tbub else 'ABOVE'} tbub)")
print(f"  phase_frac[Vap] = {vap_frac_actual:.6f}")
print(f"  entr_mol = {entr_actual:.4f} J/mol/K")

print("\n=== What FTPx.state_initialization()'s own rule says the guess SHOULD be ===")
print("(calling it directly on this same block -- this only recomputes the")
print(" guess variables, it does not re-run the solver)")
state_initialization(valve_in)
vap_frac_guess = value(valve_in.phase_frac["Vap"])
print(f"  recomputed phase_frac[Vap] guess = {vap_frac_guess:.6f}")

print("\n=== Interpretation ===")
if vap_frac_guess < 0.01 and vap_frac_actual > 0.5:
    print("  The GUESS is correct (near-zero, i.e. 'almost all liquid', as expected")
    print("  for T well below tbub) but the ACTUAL converged answer is nowhere near")
    print("  it (near-1.0, 'almost all vapor'). This means the starting point was")
    print("  fine -- the numerical solve itself is what pulls the state over to the")
    print("  wrong branch. The fix needs to target the SOLVE (e.g. bounding phase_frac")
    print("  near the correct branch during initialization, like the evaporator recipe),")
    print("  not the initial guess.")
elif vap_frac_guess > 0.5:
    print("  The GUESS ITSELF is already wrong (near-1.0) before any solve runs.")
    print("  That would mean something upstream of the solve -- e.g. tbub/tdew being")
    print("  computed inconsistently for this specific block -- is the root cause,")
    print("  not the numerical solve. Worth checking tbub/tdew values directly next.")
else:
    print(f"  Guess = {vap_frac_guess:.6f}, actual = {vap_frac_actual:.6f} -- neither")
    print("  case above matched cleanly; look at the raw numbers above directly.")
