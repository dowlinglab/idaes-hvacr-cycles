# Phase 1 — Cubic-EoS Property Validation (IDAES vs vanilla PR)

**Scripts:** `phase_1_cubic_eos_validation.py`, `pr_eos_lib.py`
**Author:** Shilpa Narasimhan · Support: Claude AI
**Date:** 2026-07-20

## Purpose

Prove that the IDAES generic Peng-Robinson property package computes the *same*
EoS physics as the standalone hand-written model (`PR_EOS.py`), for every R-32
parameter set in the Colon-group collaboration (NIST, GCGP, SPGP). Phase 0
proved *solvability*; Phase 1 proves *equivalence* of the underlying physics —
using datum-independent quantities so the arbitrary enthalpy/entropy reference
does not confound the comparison.

## Method

At a single test state (**20 °C, 10 bar**), for each parameter set:

1. **IDAES side** — build the generic PR package via `make_config(p)`, solve at
   DOF = 0, read the vapor-root compressibility `compress_fact_phase["Vap"]` and
   the total heat capacity `cp_mol_phase["Vap"]`.
2. **Vanilla side** — `vanilla_z_vap()` and `vanilla_cp_ideal()`, parameterized
   PR/Shomate functions (same physics as `pr_eos_lib.py`) taking Tc, Pc, omega,
   and the Shomate coefficients as arguments so one function serves all methods.
3. **Gate** — `|Z_idaes − Z_vanilla| < 1e-3` per method. IDAES's compiled cubic
   solver vs an independent numpy-`roots` solver of the same cubic.
4. **Anchor** — confirm `vanilla_z_vap` reproduces the *actual* standalone code
   `pr_eos_lib.z_roots` for NIST (proves the helper is faithful, so the GCGP/SPGP
   gates — which have no standalone reference — are trustworthy by transitivity).

**Why compressibility Z is the gate:** Z is datum-independent and is the direct
solution of the cubic — exactly "did IDAES solve the same equation of state."

## Result — ALL GATES PASSED

| method | Z_idaes | Z_vanilla | \|dZ\| | Cp_idaes (total) | Cp_ideal |
|--------|---------|-----------|--------|------------------|----------|
| NIST   | 0.8762  | 0.8762    | 2.0e-8 | 48.43            | 42.46    |
| GCGP   | 0.8562  | 0.8562    | 2.4e-8 | 20.76            | 14.20    |
| SPGP   | 0.8380  | 0.8380    | 2.8e-8 | 198.16           | 194.03   |

**Helper anchor (NIST):** vanilla 0.8762 vs `pr_eos_lib` 0.8758, |dZ| = 4.75e-4
(the small residual is only because `pr_eos_lib` derives omega slightly
differently, 0.277 vs 0.2769).

IDAES and the vanilla PR agree on Z to **~1e-8 (eight decimals)** for all three
methods, including SPGP with its negative acentric factor. The two solvers are
computing the identical root of the identical equation.

## Caveats

- **The Cp columns intentionally differ.** `Cp_idaes` is the *total* Cp
  (ideal-gas + PR departure); `Cp_ideal` is ideal-gas only. The gap (~4–7
  J/mol·K) is the real-gas departure at 10 bar — positive, as expected for a
  vapor below its inversion. Not an error. An exact Cp gate would add the PR
  Cp-departure term to the vanilla side; Z is the decisive check and it passed,
  so this is optional.
- **Absolute h/s still arbitrary** (zero datum) — not tested here; that is
  Phase 2's job, where the datum must be handled (anchor to IIR or compare
  differences).
- **W1002 warnings** are the same pure-fluid floating-point dust as Phase 0
  (log of a mole fraction sitting ~1e-9 above its upper bound of 0). They come
  from Pyomo, not the IDAES logger, so `outlvl` does not suppress them. Harmless.

## Gate summary

- **Gate: IDAES vapor Z matches vanilla PR to < 1e-3 for all three methods.**
  ✅ Passed at ~1e-8. Conclusion: the IDAES generic PR package is provably the
  same EoS as the standalone model, datum aside.

## Next: Phase 2

Build the saturation dome from the IDAES package and compare hf/hg/sf/sg to the
Linde datasheet — the first *accuracy* check. Expect NIST to track Linde well,
GCGP moderately, SPGP poorly (consistent with the Cp-method comparison). This is
where the datum must be resolved: anchor the IDAES package to the IIR reference,
or compare enthalpy/entropy *differences*.
