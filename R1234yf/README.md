# R1234yf Custom IDAES Helmholtz Property Package

Standalone Python/Pyomo implementation of a Helmholtz equation-of-state property
evaluator for pure R1234yf (2,3,3,3-tetrafluoropropene, HFO-1234yf), based on
Lemmon & Akasaka (2022). This is a **separate, unrelated project** from the R515B
mixture-model work elsewhere in this repo (`R515B_props_validated/`,
`PROJECT_CONTEXT.md`) -- different refrigerant, different approach (this one hand
writes the Helmholtz expressions directly in Pyomo rather than reusing IDAES's
compiled `general_helmholtz` functions, because R1234yf is not one of IDAES's
built-in registered pure fluids).

Detailed, chronological build/debug history lives in `BREADCRUMB.md`
(gitignored, not pushed -- see "Repo conventions" below). This file is the
current-state summary: what each file is, what works right now, and what doesn't.

## Status as of 2026-08-14

**The critical-point validation test passes.** Running `test_critical_point.py`
evaluates all seven core thermodynamic properties (P, h, s, cv, cp, w, rho) at
delta=tau=1 (the critical point) with no errors, no NaN/inf, and pressure/density
matching the expected critical constants to within numerical precision:

```
Pressure:  3.3843 MPa   (expected 3.3844 MPa -- 0.003% off)
Density:   476.69 kg/m^3 (expected 476.69 kg/m^3 -- exact)
```

This is a real result, not a placeholder -- it required fixing four separate bugs
(two in the EOS coefficient data, two in how the code is structured/called). Full
details in `BREADCRUMB.md`'s "Session 3" entry.

**Not yet done:** validation away from the critical point (NIST saturation-pressure
and pVT data are already identified as the right source, see `BREADCRUMB.md`, but
never actually run against the code), integration into an actual IDAES flowsheet/
cycle model, and the IDAES-proper-base-class rework needed before this can
participate in a real Pyomo NLP solve rather than standalone point evaluation (see
"Known limitations" below).

## Files

### `r1234yf.json`
The EOS parameter file: critical constants (Tc, Pc, rhoc, MW, R), ideal-gas
Helmholtz coefficients (`n0`, `g0` -- IDAES ideal Type 01 structure: log term +
linear terms + 3 Planck-Einstein oscillators), and residual Helmholtz coefficients
(`n`, `d`, `t`, `c`, `a`, `b`, `e`, `g` -- IDAES residual Type 02 structure: 5 plain
polynomial terms + 5 exponentially-damped terms + 7 Gaussian bell terms, 17 total).
Coefficients are Lemmon & Akasaka (2022) values as published in CoolProp's own
R1234yf fluid definition. Two real bugs in this file were found and fixed
2026-08-14: a missing damping-exponent (`c`) array and wrong term-boundary count
(the file previously described the residual as 2 categories, "10 power-law + 7
Gaussian," when it's actually 3: 5 undamped + 5 damped + 7 Gaussian), and a
three-way cyclic mislabeling of the Gaussian terms' `e`/`b`/`g` values (traced back
to a wrong CoolProp-parameter-name mapping written into `r1234yf_property_package.py`'s
comments, which is presumably where the original error came from). Also contains
`aux.delta_l_sat_approx`/`delta_v_sat_approx` (saturated-density correlations used
only to locate the T=273.15K reference state for the enthalpy/entropy offset) and
`transport.surface_tension` (unused by the current property calculations).

### `R1234yf.py`
**The current, working implementation.** Two classes:

- `R1234yfPropertyParameterBlock` -- loads `r1234yf.json`, stores all constants and
  EOS coefficients as plain Python floats/dicts/lists (not Pyomo `Param`/`Set` --
  see "Known limitations" below for why), and implements the actual math: ideal and
  residual dimensionless Helmholtz energy (`alpha_ideal`, `alpha_residual`,
  `alpha_total`) and their tau/delta partial derivatives up to second order
  (`alpha0_delta`, `alphar_delta_delta`, etc.), then the derived thermodynamic
  property formulas (`pressure`, `enthalpy`, `entropy`, `heat_capacity_v`,
  `heat_capacity_p`, `speed_of_sound`, `density`, `compressibility_factor`) built
  from those derivatives via the standard Helmholtz-EOS identities. Enthalpy and
  entropy are rebased to a NIST-style reference state (h=200 kJ/kg, s=1.00 kJ/(kg K)
  at saturated liquid, T=273.15K) via `h_offset`/`s_offset`, computed once during
  `build()`.
- `R1234yfPropertyStateBlock` -- a Pyomo `Block` with real `Var`s for state
  variables (P, h primary; T, rho, delta, tau auxiliary) and `Constraint`s linking
  them, intended for eventual use inside a solved flowsheet. Not yet exercised by
  any test.

Instantiation pattern (unusual, see "Known limitations"): `params =
R1234yfPropertyParameterBlock(); params.build()` -- `build()` is called manually,
not via Pyomo's normal automatic construction.

### `test_critical_point.py`
The validation script. Instantiates `R1234yfPropertyParameterBlock`, evaluates all
seven properties at delta=tau=1, and checks for errors/NaN/inf plus a sanity check
against the known critical pressure and density. **Currently passing.** Run with:
```
python3 test_critical_point.py
```

### `r1234yf_property_package.py`
**Superseded / stale -- not the current implementation.** An earlier, smaller
(175-line) attempt at the same thing, written against IDAES's real
`PropertyParameterBlock`/`PropertyStateBlock` base classes (imported but never
actually subclassed -- the file is just a flat top-level script that loads the JSON
and prints a parameter summary, "Block 1" only, no property math at all). Notable
because its inline comments contain the WRONG CoolProp-to-IDAES parameter mapping
for the Gaussian terms (`e` from beta, `b` from gamma -- backwards), which is
almost certainly where the `r1234yf.json` mislabeling bug originated. **This file
will now also error if run**, since `r1234yf.json`'s `last_term_residual` changed
from a 2-element to a 3-element list and this script's print logic assumes the old
2-element format. Kept for history; not fixed, since `R1234yf.py` has fully
superseded it. If IDAES's real base-class architecture is revisited later (see
"Known limitations"), this file is a starting point, not `R1234yf.py`.

### `test_output.txt`
**Stale.** Captured output from a test run before the 2026-08-14 bug fixes -- shows
the old Param-construction error, not the current passing state. Safe to delete or
regenerate (`python3 test_critical_point.py > test_output.txt`); not automatically
kept in sync.

### `BREADCRUMB.md`
Full chronological build/debug log, gitignored (see below). The authoritative
history if you need to know *how* a decision was reached, not just what the current
state is.

## Known limitations / open design questions

**Constants and coefficients are plain Python, not Pyomo objects.** This was a
deliberate fix (2026-08-14) for a real bug: `R1234yfPropertyParameterBlock`
subclasses bare `pyo.Block` and calls `build()` manually, which is not a pattern
plain Pyomo supports for `Param`/`Set` sub-components (they never get
auto-constructed, causing "cannot iterate/evaluate before constructed" errors the
moment anything tries to read them). Storing constants as plain floats/dicts
sidesteps this entirely and is confirmed working for point evaluation. It does
**not** solve the general problem: if this package is ever meant to provide state
variables (T, rho, P, h, etc.) that a Pyomo NLP solver actually solves for inside a
real flowsheet, those specific variables need to be genuine Pyomo `Var` objects
(as `R1234yfPropertyStateBlock` already attempts), and the whole class hierarchy
likely needs rebuilding on IDAES's own `PhysicalParameterBlock`/`StateBlockData`
base classes (via `declare_process_block_class`) rather than bare `pyo.Block`, to
get proper automatic construction. `r1234yf_property_package.py`'s unfinished
import of `PropertyParameterBlock`/`PropertyStateBlock` shows this was already
recognized once before and not completed.

**Cp/speed-of-sound diverge at the critical point.** This is correct physics, not a
bug -- both formulas involve dividing by `(dP/drho)_T`, which is genuinely zero
exactly at a true critical point. Confirmed the test output's huge Cp value
(~1.24e15 J/(kg K)) is expected at delta=tau=1 specifically; not yet checked
whether Cp behaves sanely at nearby, non-critical states.

**No phase-equilibrium (VLE) handling at all.** Everything here is single-phase
point evaluation given (delta, tau). Saturated-liquid/vapor density correlations
exist in the JSON (`aux.delta_l_sat_approx`/`delta_v_sat_approx`) but are only used
internally to locate the reference state -- there's no bubble/dew-point solver, no
two-phase region handling, nothing analogous to the R515B mixture work's VLE
machinery.

## Repo conventions

- `BREADCRUMB.md` is gitignored (see `.gitignore`) and stays local -- it's the
  working log, not meant to be pushed.
- This file (`README.md`) and the actual code/data files (`R1234yf.py`,
  `r1234yf.json`, `r1234yf_property_package.py`, `test_critical_point.py`) are
  meant to be pushed.
- `__pycache__/` and `*.nl` files are gitignored per the existing `.gitignore`.
