# R1234yf Custom IDAES Helmholtz Property Package

Standalone Python/Pyomo implementation of a Helmholtz equation-of-state property
evaluator for pure R1234yf (2,3,3,3-tetrafluoropropene, HFO-1234yf), based on
Lemmon & Akasaka (2022). This is a **separate, unrelated project** from the R515B
mixture-model work elsewhere in this repo (`R515B_props_validated/`,
`PROJECT_CONTEXT.md`) -- different refrigerant, different approach (this one hand
writes the Helmholtz expressions directly in Pyomo rather than reusing IDAES's
compiled `general_helmholtz` functions, because R1234yf is not one of IDAES's
built-in registered pure fluids).

This file is the current-state summary: what each file is, what works right
now, and what doesn't.

## Status as of 2026-08-17

**The full saturated-dome / p-H diagram machinery is built and validated.**
Beyond the 2026-08-14 critical-point result (still valid, see below),
`R1234yf_validation.py` now builds a complete p-H diagram: the two-phase
saturation dome (bubble/dew density solved point-by-point), constant-quality
lines inside the dome, isotherms, and isentropes (including correct handling
of isentropes that cross the two-phase dome). This has been validated two
independent ways:

1. **Against the Danfoss R1234yf p-H chart** (`R1234yf SI Units.pdf`) -- the
   dome shape, isotherm pattern, and all 38 of the chart's real isentrope
   values (775-1575 J/(kg K) in steps of 50, 1625-2125 in steps of 25) are
   reproduced and visually match.
2. **Against independent NIST experimental data** -- Richter, McLinden &
   Lemmon (2011, J. Chem. Eng. Data 56, 3254-3264), 135 raw experimental
   points (30 vapor-pressure, 105 p-rho-T) compared point-by-point against
   this code's own computed values (not against the paper's own, different,
   EOS fit). Vapor pressure: mean deviation 0.056%, worst 0.541% (at
   T=250.002K, the coldest point tested). p-rho-T, normal region (93 pts):
   mean 0.045%, max 0.277%. p-rho-T, near-critical region (12 pts, compared
   on a pressure basis per the source paper's own methodology since density
   is hypersensitive near Tc): mean 0.169%, max 0.473%.

Full detail, methodology, and caveats (what is and isn't covered) are in
`VALIDATION_REPORT.md`.

**Older result, still valid:** the critical-point test (`test_critical_point.py`)
evaluates all seven core thermodynamic properties (P, h, s, cv, cp, w, rho) at
delta=tau=1 (the critical point) with no errors, no NaN/inf, and pressure/density
matching the expected critical constants to within numerical precision:

```
Pressure:  3.3843 MPa   (expected 3.3844 MPa -- 0.003% off)
Density:   476.69 kg/m^3 (expected 476.69 kg/m^3 -- exact)
```

This required fixing four separate bugs (two in the EOS coefficient data, two
in how the code is structured/called).

**Not yet done:** integration into an actual IDAES flowsheet/cycle model, and
the IDAES-proper-base-class rework needed before this can participate in a
real Pyomo NLP solve rather than standalone point evaluation (see "Known
limitations" below). Also not yet validated: Table 7 speed-of-sound/cp data
from the NIST paper, and pixel-level overlay against the Danfoss chart (only
visual comparison so far). See `VALIDATION_REPORT.md` section 7 for the full
list of what remains open.

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
  entropy are rebased to a NIST/IIR-style reference state (h=200 kJ/kg, s=1.00
  kJ/(kg K) at saturated liquid, T=273.15K) via `h_offset`/`s_offset`, computed
  once during `build()`. Confirmed computationally (2026-08-17) that both offsets
  are correctly 0.0 -- the underlying EOS's own ideal-gas constants already put
  raw h/s at this state within 0.0001% of the IIR target (200000.23 J/kg and
  1000.0011 J/(kg K) respectively), so no additional rebasing is needed.
- `R1234yfPropertyStateBlock` -- a Pyomo `Block` with real `Var`s for state
  variables (P, h primary; T, rho, delta, tau auxiliary) and `Constraint`s linking
  them, intended for eventual use inside a solved flowsheet. Not yet exercised by
  any test.

Instantiation pattern (unusual, see "Known limitations"): `params =
R1234yfPropertyParameterBlock(); params.build()` -- `build()` is called manually,
not via Pyomo's normal automatic construction.

### `R1234yf_validation.py`
**The p-H diagram driver, built and validated 2026-08-17.** Given a built
`R1234yfPropertyParameterBlock`, this script:

- Solves the two-phase saturation dome point-by-point (per-temperature
  bubble/dew density pairs from T just above the triple point up to Tc).
- Draws constant-quality lines inside the dome, and isotherms across the
  full diagram.
- Constructs isentropes matching the Danfoss chart's real entropy values
  (38 values total), including correctly stitching together the single-phase
  branches with the two-phase dome-crossing segment where an isentrope
  passes through the two-phase region (`compute_full_isentrope`, which
  internally calls `isentrope_point_residuals`, `solve_isentrope_point`,
  `compute_isentrope`, and `isentrope_two_phase_segment`). Two real bugs were
  found and fixed here: a reversed two-phase segment that caused isentropes
  to visibly zigzag/cross over themselves, and a silent near-critical solver
  failure where `least_squares` reported success but had actually converged
  to the wrong root (fixed by checking the actual entropy residual, not just
  the solver's own success flag, and trying multiple seed points).
- Renders the full diagram to `dome.png`.

Also contains the reference-state verification and the NIST experimental
point-by-point comparison logic. Run with:
```
python3 R1234yf_validation.py
```

### `VALIDATION_REPORT.md`
**Written 2026-08-17.** The validation report: critical point (solved vs. the
EOS's own circular fit vs. Tanaka & Higashi 2010's independent measurement),
reference-state check, saturation dome solver convergence, the p-H diagram
vs. the Danfoss chart, the two bugs found and fixed in isentrope construction,
the full 135-point NIST experimental comparison (vapor pressure + p-rho-T,
both described above), and an honest list of what is still not validated
(Table 7 speed-of-sound/cp data, pixel-level diagram overlay, full-flowsheet
behavior).

### `dome.png`
The rendered p-H diagram: saturation dome, quality lines, isotherms, and
isentropes, produced by `R1234yf_validation.py`.

### `R1234yf SI Units.pdf`
The reference Danfoss R1234yf p-H chart (SI units). Used as the visual
target for the diagram and as the source of the 38 real isentrope entropy
values used in `R1234yf_validation.py`.

### `test_critical_point.py`
The critical-point validation script. Instantiates `R1234yfPropertyParameterBlock`,
evaluates all seven properties at delta=tau=1, and checks for errors/NaN/inf plus a
sanity check against the known critical pressure and density. **Currently passing.**
Run with:
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

**Phase-equilibrium (VLE) handling is now point-solver-level, not
flowsheet-level.** As of 2026-08-17, `R1234yf_validation.py` has a working
saturation dome solver (per-temperature bubble/dew density pairs), quality
lines, and isentropes/isotherms that correctly cross the two-phase region --
this is genuine bubble/dew-point VLE machinery, not just internal reference-state
plumbing. What's still missing is integration of this into
`R1234yfPropertyStateBlock`/a real flowsheet: there's no flash calculation
exposed as a Pyomo constraint set, nothing analogous to the R515B mixture
work's in-flowsheet VLE machinery.

## Repo conventions

- This file (`README.md`) and the actual code/data/report files --
  `R1234yf.py`, `r1234yf.json`, `r1234yf_property_package.py`,
  `test_critical_point.py`, `R1234yf_validation.py`, `VALIDATION_REPORT.md`,
  `dome.png`, `R1234yf SI Units.pdf` -- are meant to be pushed.
- `__pycache__/` and `*.nl` files are gitignored per the existing `.gitignore`.
