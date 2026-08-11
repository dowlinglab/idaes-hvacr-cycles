# R32 vapor-compression cycle -- file guide

R32 (difluoromethane) vapor-compression cycle models, comparing four
thermodynamic property sources -- a reference Helmholtz-energy equation
of state, and three Peng-Robinson cubic-EoS property fits (NIST, GCGP,
SPGP) -- against each other via COP (coefficient of performance).

## Core files

**`phase6_final.py`** -- the full COP sweep across all four methods
(Helmholtz, NIST, GCGP, SPGP), at four ambient temperatures (10/15/20/25
C). Each method is run through an iterative warmstart chase: if a
method's COP-vs-ambient trend isn't monotonically decreasing on the
first attempt, it's re-solved warmstarted from its own prior attempt's
converged points, up to 10 attempts, until the trend is fully converged
and physically monotonic. Reports a combined COP table plus %-difference
vs. Helmholtz, and how many attempts each method needed.

**`phase6_final_sensitivity_Tc.py`** -- sensitivity analysis on NIST's
critical constants. Perturbs NIST's critical temperature (Tc) and/or
critical pressure (Pc) by +/-{0.01%, 0.1%, 0.5%, 1%} (plus a 0%
baseline) relative to NIST's fitted values, in three sweeps: Tc only
(Pc held at baseline), Pc only (Tc held at baseline), and the full 9x9
grid of both varied together. Every perturbed combination goes through
the same warmstart chase as phase6_final, across the same four
ambients. Reports to `phase6_sensitivity_NIST.xlsx` (one sheet per
sweep).

**Property-package files** -- define the actual thermodynamic models
both scripts above build their cycles from:

- `phase_1_cubic_eos_validation_refstate.py` -- NIST/GCGP/SPGP's
  critical constants (Tc, Pc, omega), Shomate ideal-gas Cp
  coefficients, and IIR-calibrated reference-state offsets (F/G),
  plus `make_config()`, which turns those into an IDAES
  `GenericParameterBlock` configuration.
- `vapor_compression_cubic_refstate.py` -- the cubic-PR cycle model
  (`SimpleVaporCompressionCycle`/`CubicCycle`), used for the NIST/GCGP/
  SPGP methods.
- `vapor_compression.py` -- the Helmholtz-EoS cycle model
  (`SimpleVaporCompressionCycle`/`HelmCycle`), molar basis, used for
  the Helmholtz reference method.

## Everything else

All other files in this folder were used for debugging and/or
validation during development -- earlier/superseded versions of the
cycle and property-package files, per-phase validation and sweep
scripts (phase_0 through phase5), unit-level debug scripts (valve/
compressor/composition-loop initialization issues), and property-fit
comparison scripts. They aren't needed to reproduce `phase6_final.py`'s
or `phase6_final_sensitivity_Tc.py`'s results, but are kept for
reference on how the property packages and cycle model were built up
and validated.
