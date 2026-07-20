# Phase 0 — Cubic-EoS Solvability Test

**Script:** `phase_0_cubic_roots_test.py`
**Author:** Shilpa Narasimhan · Support: Claude AI
**Date:** 2026-07-20

## Purpose

Prove that IDAES can build and solve the Peng-Robinson cubic equation of state,
at **0 degrees of freedom**, for *every* R-32 parameter set in the Colon-group
collaboration (NIST, GCGP, SPGP). This is the entry gate for the whole pipeline:
if a parameter set cannot even be solved as a single fixed state, it cannot be
run in a vapor-compression cycle.

Phase 0 tests **solvability only** — not accuracy. Property validation is Phase 1
(Z, Cp vs `PR_EOS.py`) and Phase 2 (saturation dome vs Linde).

## What the test does

1. `METHODS` holds one dict per parameter set: `Pc`, `Tc`, `omega`, and Shomate
   `A–E`. `omega` is the Pitzer acentric factor from the Linde vapor pressure at
   Tr = 0.7 (`omega = -1 - log10(Psat/Pc)`), using each method's own Tc, Pc.
2. `make_config(p)` returns the IDAES generic-property configuration for one set:
   both phases use the cubic PR EoS; ideal-gas Cp/h/s use the NIST Shomate method
   fed the coefficients; saturation uses an Antoine init guess fit to Linde; VLE
   via SmoothVLE + LogBubbleDew + log-fugacity.
3. For each method: build the package, create one state block, fix R-32 at
   **20 °C, 10 bar** (flow = 1 mol/s, x = 1), confirm DOF = 0, initialize, solve,
   and print T, P, h, s, and phase fractions.

No `try/except` (per SOP): a failure stops the run loudly at the offending method.

## Result — GATE PASSED

All three methods built at DOF = 0 and reported "Optimal Solution Found" (SPGP:
"Found feasible point for square problem"). SPGP solved despite its negative
acentric factor, which we had flagged as a possible failure point.

| method | solved | h [J/mol] | s [J/mol·K] |
|--------|--------|-----------|-------------|
| NIST   | yes    | 1002.5    | 26.05       |
| GCGP   | yes    | 1966.7    | -42.98      |
| SPGP   | yes    | 45130.6   | -125.53     |

## Caveats (read before trusting numbers)

- **Absolute h/s are not physical yet.** The package uses a zero datum
  (`F = G = H = 0`, `include_enthalpy_of_formation = False`), so the absolute
  values are arbitrary. Only datum-independent quantities (Z, Cp, and h/s
  *differences*) are meaningful — and COP uses only differences. To get
  Linde-comparable absolute h/s, anchor F, G to the IIR datum.
- **SPGP's h is ~20–40× the others** — consistent with its poor Cp/critical
  properties from the method comparison. It converges, but expect it to be the
  problem child in the cycle.
- **Phase fractions (~9% liq / ~91% vap at 20 °C/10 bar) are a SmoothVLE
  artifact, not physics.** A pure fluid off its saturation curve is single-phase
  (at 20 °C, R-32's Psat ≈ 14.6 bar, so at 10 bar it is superheated vapor).
  SmoothVLE smears the transition to keep the solver stable and returns a
  fractional split near the boundary. Harmless here; note it when reading
  single-phase properties later.
- **`W1002` init warnings are floating-point dust.** `log_mole_frac_tbub/tdew`
  have upper bound 0 (log of a mole fraction ≤ 1). For a pure fluid the mole
  fraction is exactly 1, so the log sits on 0; the solver lands ~1e-9 above it
  and IDAES warns. Physically zero, no effect on the solve.

## Gate summary

- Prerequisite: `idaes get-extensions` in `myidaesenv` (cubic-root functions +
  ipopt). Confirmed via `idaes environment-info`.
- **Gate: all three methods solve at DOF = 0.** ✅ Passed.

## Next: Phase 1

Pull `Z` (and `cp_mol`) from the solved state block and compare to `PR_EOS.py` at
20 °C, 10 bar. Target: vapor **Z ≈ 0.876**. If Z and Cp match, the IDAES package
is provably the same physics as the standalone model, datum aside.
