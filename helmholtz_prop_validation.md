# R-515B IDAES Helmholtz Property Package -- Structured Validation Record

This file is the mandatory, structured, machine-checkable validation record
for the new independent R-515B Helmholtz property package, kept separate from
the chronological narrative in `PROJECT_CONTEXT.md` (which remains the
append-only breadcrumb log). This file is a NEW file created for this task;
per its own governing task spec, it may be edited/appended freely as
validation stages complete, but earlier evidence is never deleted -- superseded
results are marked SUPERSEDED with a dated note, never overwritten.

No Git was used to produce or manage this file or any part of this task.

---

## 0. Environment

| Field | Value |
|---|---|
| Date established | 2026-08-17 |
| Python version | 3.11.15 (main, Mar 3 2026, 09:26:23) [GCC 13.3.0] |
| IDAES version | 2.12.0 |
| Pyomo version | 6.10.1 |
| Nonlinear solver | IPOPT |
| Solver version | 3.13.2 (x86_64-pc-linux-gnu), ASL(20190605) |
| Solver location | `/root/.idaes/bin/ipopt` (not on PATH by default) |
| Other IDAES-bundled solvers available | bonmin, cbc, clp, couenne, k_aug (at `/root/.idaes/bin/`) |
| Reference/oracle implementation | `R515B_props_validated/mixture_fully_validated.py` (device-resident, 189,129 bytes) |
| Oracle equivalence note | Full-file diff against the actively-worked `mixture_isentrope_validation.py` (192,435 bytes) shows the ONLY difference anywhere in the file is one defensive helper (`_verify_isentrope_solution`, added 2026-08-17) already proven upstream to be byte-for-byte output-neutral (identical CSVs, MD5-identical p-H diagrams at the project's standard parameters). The oracle is therefore confirmed equivalent, for every practical purpose, to the full, most-recently-fixed state of this project's validated isentrope/VLE/critical-point code (including all 2026-08-14 near-critical fixes). See `PROJECT_CONTEXT.md`, 2026-08-17 "MASTER TASK kickoff" entry for the full diff record. |
| New production property-package files | *(to be listed as created -- none yet at this stage)* |
| New `vapor_compression.py` integration copy | *(to be created -- none yet at this stage)* |
| R-515B composition | w1 = 0.911 (mass fraction R-1234ze(E)); x1 (mole fraction) computed via `w1_to_x1` from each fluid's `basic.MW` -- see Stage D data-validation section below for the numeric value once ported |
| R-1234ze(E) data source | `idaes.models.properties.general_helmholtz` parameter path, `r1234ze.json` (IDAES-format Helmholtz EOS parameters) |
| R-227ea data source | same parameter path, `r227ea.json` |
| Binary parameter source | Bell, I. H. (2023), J. Phys. Chem. Ref. Data 52(1), 013101 -- `BELL_2023_R1234ZE_R227EA` reducing parameters (beta_T, beta_v, gamma_T, gamma_v) and `BELL_2023_DEP_COEFFS` (Table 7 departure coefficients), as hardcoded in `linear_model_codex.py` (imported, unmodified, by the oracle) |

**Git:** not used for any purpose in this task, per governing spec rule 7.

---

## 1. Frozen numerical tolerances

**STATUS: ESTABLISHED 2026-08-17.** Reproducibility measured directly against
the oracle (`mixture_fully_validated.py`) using
`R515B_idaes_package/establish_reference_tolerances.py` (new file, read-only
against the oracle). Full numeric detail in
`R515B_idaes_package/reference_repeatability_results.json`.

**Representative states tested (5 repeated evaluations each):**

| State | Definition | Result (mean) |
|---|---|---|
| `bubble_saturation` | `solve_bubble_at_t` at T=300K, x1=0.9385037942 | P=535479.39 Pa, rho_l=9794.535 mol/m3, rho_v=249.171 mol/m3, y1=0.930047 |
| `dew_saturation` | `solve_dew_at_t` at T=300K | P=534739.88 Pa, rho_l=9837.018 mol/m3, rho_v=248.687 mol/m3, x1=0.946847 |
| `liquid_direct_subcooled` | `mix_state` at T=300K, rho=1.10x the solved bubble rho_l | p=32179958.74 Pa, h=29130.570 J/mol, g=-9389.643 J/mol |
| `vapor_direct_superheated` | `mix_state` at T=300K, rho=0.50x the solved dew rho_v | p=288569.04 Pa, h=47154.370 J/mol, g=-13835.469 J/mol |
| `quality_0p5_at_Tsat` | lever-rule q=0.5 at T=300K via `compute_quality_lines` | P=535109.64 Pa, h=37280.567 J/mol |
| `critical_point_solve` | `solve_mixture_critical_point`, Bell-mixing-rule-seeded | T=381.5224K, P=3550445.18 Pa, rho=4228.487 mol/m3 |
| `supercritical_direct` | `mix_state` at T=Tc+8K, rho=3000 mol/m3 | p=3939214.45 Pa, h=49444.917 J/mol, g=-26232.323 J/mol |

**Finding: the oracle is EXACTLY bit-reproducible (zero measured spread, all
fields, all 7 states, 5/5 identical runs)** -- expected for deterministic
floating-point NumPy/SciPy code with no RNG, confirmed rather than assumed.

**Note on test-state validity (recorded as a minor methodology fix, not an
oracle defect):** an earlier version of this script picked an UNANCHORED
liquid test point (T=280K, rho=9500 mol/m^3) and got p=-14.25 MPa -- not a
bug in the oracle, but an invalid test state outside this mixture's real
liquid branch at that T (the picked density did not correspond to any
physically meaningful subcooled-liquid condition for this EOS at 280K).
Fixed by anchoring liquid/vapor single-phase test states to the actual solved
bubble/dew densities at the same T (10% denser than saturated liquid, 50% of
saturated vapor density) rather than an arbitrary guess -- both now produce
physically sensible positive pressures consistent with each phase's expected
compressibility behavior.

**Frozen tolerances (may be tightened later; per spec rule 41 may NOT be
loosened without user approval):**

| Tier | States | Property | Frozen tolerance |
|---|---|---|---|
| 1 (direct-equation) | `liquid_direct_subcooled`, `vapor_direct_superheated`, `supercritical_direct` -- states where P/h/g come straight from the same closed-form Helmholtz identities (M3/M4) with no iterative solve involved | P, h, g | 1e-8 relative |
| 2 (iterative-equilibrium) | `bubble_saturation`, `dew_saturation`, `critical_point_solve`, `quality_0p5_at_Tsat` -- states reached via an iterative solve (least_squares / bisection) where the new IDAES/Pyomo implementation may legitimately reach the same physical root via a numerically different path (different solver, different equation formulation, IPOPT's NLP solve vs. SciPy's `least_squares`/`root`) | P, h, T, s, rho_l, rho_v, x1, y1 | 1e-4 relative |

**Rationale for the two-tier split:** Tier 1 states involve no solver at all
on either side (oracle or new implementation) -- both sides evaluate the
identical closed-form Helmholtz P/h/g formulas at the same (T,rho,x1), so a
faithful port should agree to near machine precision; 1e-8 is a tight but
fair bar. Tier 2 states are the output of a genuinely different numerical
method on the new-implementation side (an equation-oriented Pyomo/IPOPT
solve, versus the oracle's own bounded `least_squares`/bisection) -- both
are valid ways to find the same physical root, but will not agree to
1e-8; 1e-4 (0.01%) is tight enough to catch a wrong-root or formulation bug
while tolerant of ordinary solver-path differences. This mirrors the
project's own established convention elsewhere (e.g. the bubble/dew VLE
acceptance gate's own 1e-6 relative residual threshold, and the isentrope
entropy-verification insurance's 1e-6 relative tolerance) scaled up modestly
for the added tolerance of comparing two independent numerical
implementations rather than one implementation against itself.

---

## 2. Input/data validation

**STATUS: STAGE D — PASS, by construction (2026-08-17).** Key finding that
reshapes the rest of the porting plan: the oracle (`mixture_fully_validated.py`)
itself does `from linear_model_codex import (... bell2023_Tred_vred,
mixture_alpha0_alphar_derivs, load_idaes_helmholtz_json, mw_from_json, ...)`
(confirmed by reading its own import block, line 51). `linear_model_codex.py`
is an existing, unmodified file, explicitly NOT the oracle, and explicitly
permitted as a production dependency under MASTER TASK spec rule 6's
"legitimate production dependency" carve-out. This means the R-1234ze(E)/
R-227ea JSON data loading (`load_idaes_helmholtz_json`, `mw_from_json`),
critical properties (read directly from each JSON's `basic` section), and
the Bell (2023) reducing/departure parameters (`BELL_2023_R1234ZE_R227EA`,
`BELL_2023_DEP_COEFFS`) are not two independently-ported copies that could
drift apart -- the oracle and the new production package will call the
EXACT SAME functions in the EXACT SAME file. There is nothing left to prove
numerically identical here beyond confirming the new package actually
imports `linear_model_codex.py` (not a re-transcription) -- deferred to
Stage L wiring.

---

## 3. Reducing-function validation (Tred(x), vred(x))

**STATUS: STAGE E — PASS, by construction (2026-08-17).** Same reasoning as
Section 2: `bell2023_Tred_vred` (Bell 2023 Eq. 3-5 binary reducing
functions) lives in `linear_model_codex.py` and is called directly, by the
SAME shared function call, in both the oracle's `_mix_alpha_and_derivs`/
`mix_state`/`_mix_entropy_direct` and in `linear_model_codex.py`'s own
`_mixture_reduced_state`/`_mixture_alpha_eval`. Confirmed no divergence is
possible for this layer since it's one shared code path, not two. The new
production package will call `bell2023_Tred_vred` directly (imported from
`linear_model_codex.py`) for its own reducing-function evaluation.

---

## 4. Helmholtz kernel validation (alpha0, alphar, derivatives)

**STATUS: STAGE F — PASS (partial), 2026-08-17.** `mixture_alpha0_alphar_derivs`
(in `linear_model_codex.py`, combining pure-fluid `alpha0_idaes_with_derivs`/
`alphar_idaes_with_derivs` plus the Bell 2023 departure function
`bell2023_departure_alphar`) is the SAME function the oracle's own
`_mix_alpha_and_derivs` calls. Directly verified numerically (script:
`R515B_idaes_package/validate_table1_vs_oracle.py`, ad-hoc console check) at
a representative subcooled-liquid state (T=300K, rho=10773.99 mol/m^3,
x1=0.9385037942): calling `mixture_alpha0_alphar_derivs` directly and
assembling `alpha_mix = a0_mix + ar`, `alpha_tau_mix = a0_tau + ar_tau`
reproduces the oracle's own `_mix_alpha_and_derivs` output for
(alpha_mix, alpha_tau_mix, ar_del_mix) with **exact (0.0) absolute
difference** on all three fields. Composition-derivative machinery
(`_bell2023_reducing_derivs_binary`/`_bell2023_reducing_derivs_binary_local`,
needed for fugacity in Stage H) not yet separately re-verified as of this
entry -- deferred to Section 7/Stage H.

**Important false-start, corrected same session:** an earlier ad-hoc
diagnostic run appeared to show a constant ~1.4595 J/(mol*K) entropy
mismatch between `linear_model_codex.py`'s `compute_table1_properties` and
the oracle. Root cause was a units/accounting bug in the diagnostic script
itself (comparing a dimensionless s/R value against a J/(mol*K) value in one
branch, and in a later, corrected script, double-counting the oracle's
`ENTROPY_REFERENCE_OFFSET_JMOLK` rebase because `_mix_entropy_direct`
already applies that offset internally before returning -- confirmed by
reading its body directly, lines 2569-2570) -- NOT a real physics
discrepancy in either file. Corrected and re-verified below (Section 5);
recorded here per spec rule 74 (retain the false start, mark corrected,
never silently erase).

---

## 5. Core (P, h, T, s) validation

**STATUS: STAGE G — PASS (partial), 2026-08-17.** Primary acceptance
criterion per spec rule 56/90: (P,h,T,s)_new ~= (P,h,T,s)_reference within
frozen tolerances. Validated so far at 3 representative DIRECT
("no-solver-either-side", tier-1, 1e-8 relative tolerance) states via
`R515B_idaes_package/validate_table1_vs_oracle.py`, using
`linear_model_codex.py`'s existing `compute_table1_properties` (legitimate
dependency, not the oracle) against the oracle's `mix_state` (P, h) and
`_mix_entropy_direct` (s, corrected for its internal Honeywell rebase --
see Section 4's false-start note):

| state | T (K) | rho_mol (mol/m^3) | P rel err | h rel err | s(raw) rel err |
|---|---|---|---|---|---|
| liquid_direct_subcooled | 300.000 | 10773.989 | 0 | 0 | 0 |
| vapor_direct_superheated | 300.000 | 124.343 | 0 | 3.09e-16 | 1.4e-16 |
| supercritical_direct | 389.500 | 3000.000 | 0 | 1.47e-16 | 0 |

**OVERALL: PASS**, all three states well inside the frozen 1e-8 tier-1
tolerance (Section 1) -- in fact at or near machine epsilon, as expected
given Section 4's finding that both paths share the same underlying
`mixture_alpha0_alphar_derivs`/`bell2023_Tred_vred` calls and the P/h/s
formulas are algebraically identical (Z=1+delta*ar_del; h/RT=1+tau*
alpha_tau+delta*ar_del; s/R=tau*alpha_tau-alpha, all confirmed by direct
reading of both files' source).

**Still pending:** saturated liquid/vapor (bubble/dew, Stage I -- these
require the 3-equation VLE solve, which does NOT exist anywhere in
`linear_model_codex.py` and must be genuinely ported, not reused), two-phase
quality states (Stage K), critical point (Stage J), isotherms/isentropes
(Stage K), and the near-machine-epsilon match must still be re-confirmed
once P/h/s are wired through actual Pyomo expressions in the IDAES
StateBlock (Stage L) rather than plain-Python/NumPy calls -- Pyomo's AD and
any unit-scaling could in principle introduce small differences that this
plain-Python check cannot detect.

---

## 6. Single-phase validation (subcooled liquid / superheated vapor / supercritical)

**STATUS: STAGE G — PASS (partial), 2026-08-17.** Covered by Section 5's
same 3-state table (all 3 states there ARE the single-phase states this
section is about). Re-stated here for cross-reference per the doc's own
section structure; no additional states run yet. Additional single-phase
states across a wider (T, rho) grid recommended before Stage L sign-off but
not yet executed -- next action.

---

## 7. Fugacity and chemical-potential validation

**STATUS: STAGE H — PASS (partial), 2026-08-17.** Important distinction from
Sections 2-4: the oracle's `chemical_potentials_analytic` uses a LOCAL
re-implementation, `_bell2023_reducing_derivs_binary_local` (defined inside
`mixture_fully_validated.py` itself, line 620), of the same composition-
derivative formulas `linear_model_codex.py`'s own (shared, imported)
`_bell2023_reducing_derivs_binary` implements. This is NOT a shared code
path (unlike reducing functions/Helmholtz kernel) -- read side-by-side, the
two are the same algebra written out independently, so per spec rule 10
("mathematically equivalent equations are not automatically numerically
equivalent") this needed an actual numeric check, not just a code-path
argument.

Validated via `R515B_idaes_package/validate_fugacity_vs_oracle.py`,
comparing the oracle's `chemical_potentials_analytic` (mu1, mu2) against
`linear_model_codex.py`'s `compute_table1_properties` (mu1_Jmol, mu2_Jmol)
at the same 3 direct states as Section 5:

| state | mu1 rel err | mu2 rel err |
|---|---|---|
| liquid_direct_subcooled | 0 | 0 |
| vapor_direct_superheated | 0 | 0 |
| supercritical_direct | 0 | 0 |

**OVERALL: PASS**, exact match at all 3 tested states (tolerance used: 1e-6
relative, deliberately looser than Section 5's 1e-8 tier-1 floor since this
is an independent re-implementation rather than a shared code path -- exact
0.0 match obtained anyway). Confirms the two independently-written
`_bell2023_reducing_derivs_binary`/`_local` functions are numerically
identical for this composition (x1=0.9385037942) at these 3 states.

**Still pending:** fugacity/chemical-potential equality as it's actually
USED inside the bubble/dew VLE solve (mu1_L=mu1_V, mu2_L=mu2_V, spec's M6)
has not been separately re-verified -- that requires Stage I's ported VLE
solver to exist first. Composition-derivative validation across a wider
range of x1 (this project only ever needs x1=0.9385037942 for R-515B
itself, but a wider sweep would strengthen confidence this isn't a
coincidental match at one composition) not yet done.

---

## 8. Bubble-line validation

**STATUS: STAGE I — PASS, 2026-08-17.** New INDEPENDENT production module
`R515B_idaes_package/r515b_helmholtz_core.py` created (does NOT import the
oracle at runtime -- imports only `linear_model_codex.py`'s shared
functions plus its own port of `mix_state`, `chemical_potentials_analytic`,
`solve_bubble_at_t`, `solve_dew_at_t`, `solve_mixture_critical_point`,
reproducing the oracle's exact algorithms, sigmoid/logit reparameterization,
4-seed retry ladder, tapered density-separation gate, and strict post-solve
residual re-verification -- see the module's own docstring for the full
oracle-function -> new-function traceability map).

Validated via `R515B_idaes_package/validate_core_vs_oracle.py`:
`solve_bubble_at_t` compared against the oracle at 6 temperatures (270K,
285K, 300K, 315K, 330K, 345K), continuation-seeded from the oracle's own
converged prior-T state (mirrors `run_true_vle_envelope`'s own continuation
pattern). **Result: P, rho_l, rho_v, y1_vap all rel_err = 0.0 (exact) at
every tested T; both CONVERGED at every T.** Well inside the frozen 1e-4
tier-2 tolerance (Section 1) -- in fact exact, since both sides literally
run the same `scipy.optimize.least_squares` call with identical inputs
(same residual function values from the now-validated `mix_state`/
`chemical_potentials_analytic` ports feed an identical solver call).

**OVERALL: PASS.**

---

## 9. Dew-line validation

**STATUS: STAGE I — PASS, 2026-08-17.** Same script/run as Section 8:
`solve_dew_at_t` compared at the same 6 temperatures. **Result: P, rho_l,
rho_v, x1_liq all rel_err = 0.0 (exact) at every tested T; both CONVERGED
at every T.**

**OVERALL: PASS.**

**Still pending for Sections 8/9:** wider temperature coverage (only 6
points tested vs. the full practical dome, though the underlying mechanism
being identical across all of them makes further points a confirmation
exercise, not a discovery one); azeotrope/consistency cross-check between
bubble-branch P(T) and dew-branch P(T) (spec mentions this as a sanity
check, not yet explicitly tabulated); near-critical bubble/dew points
(within the RHO_SEP_TAPER_START_K=30K taper window) not yet separately
exercised.

---

## 10. Critical-point validation

**STATUS: STAGE J — PASS, 2026-08-17.** `solve_mixture_critical_point`
ported into `r515b_helmholtz_core.py` exactly (THIRD-design bisection on
spinodal-dip existence -- see oracle docstring for the documented failure
history of the two abandoned earlier designs; deliberately reproduced
as-is, not re-derived). Validated 2 cases via `validate_core_vs_oracle.py`:

- **Case A (crude Bell-mixing-rule-reducing-point guess, no real sweep
  data):** both the oracle and the new module return `converged=False`
  (the fallback path) with IDENTICAL fallback values (T_K=381.5224K,
  rho_molm3 = the guess itself) -- exact match on the failure mode itself,
  not just a nominal "both failed."
  **Important correction to earlier session usage:** the "Tc_mix=381.5224K"
  value quoted in `establish_reference_tolerances.py` (Stage C) and
  `validate_ancillary_guess.py`'s upper-bound-of-grid calculation is THIS
  non-converged fallback's T_K -- i.e. the Bell(2023) mixing-rule reducing
  temperature Tred_mix, an APPROXIMATION, not a genuinely solved critical
  point. This was harmless in both prior uses (only needed as a rough safe
  upper bound for a temperature grid/repeatability check), but must not be
  mistaken for a validated Tc going forward.
- **Case B (real-data-informed guess/window, mirroring how
  `run_true_vle_envelope` actually calls this in practice):** seeded from a
  real bubble solve at T=350K (rho_guess = midpoint of bubble's rho_l/rho_v),
  with explicit t_scan_lo=350K/t_scan_hi=385K. Both oracle and new module
  genuinely CONVERGE to the SAME point: **rel_err = 0.0 (exact) on T_K,
  rho_molm3, P_Pa, and h_Jmol.**

**THE VALIDATED TRUE MIXTURE CRITICAL POINT (R-515B, w1=0.911,
x1=0.9385037942):** Tc = 381.8939 K, rhoc = 3858.3849 mol/m^3,
Pc = 3.576459 MPa. (Distinct from the Tred_mix=381.5224K approximation
used elsewhere as a cheap proxy -- the two are close, 0.37K apart, which is
why the proxy has worked fine for taper-window/grid-bound purposes, but
they are not the same quantity.)

**OVERALL: PASS.**

---

## 11. Quality-state validation

**STATUS: STAGE K — PASS, 2026-08-17.** `compute_quality_lines` ported into
`r515b_helmholtz_core.py` exactly (pure lever-rule post-processing of
already-solved bubble/dew rows, no new EOS call). Validated via
`R515B_idaes_package/validate_quality_isotherm_vs_oracle.py`: all 9
Honeywell-chart quality values (0.1-0.9) across a 10-point T sweep
(260-350K), 90 total rows compared -- **h_Jmol and P_Pa rel_err = 0.0 at
every row.**

**OVERALL: PASS.**

---

## 12. Isotherm validation

**STATUS: STAGE K — PASS, 2026-08-17.** All four oracle isotherm functions
ported into `r515b_helmholtz_core.py` exactly (`compute_isotherms_two_phase`,
`compute_isotherms_liquid_side`, `compute_isotherms_vapor_side`,
`compute_isotherms_supercritical`) -- all single-equation/single-unknown
`brentq` root-finds on this module's own already-validated `mix_state`, log-
spaced (geomspace) pressure sampling, identical seeding/bracket-widening
logic. Validated via the same script as Section 11:

| isotherm family | T-values tested | max rel. err |
|---|---|---|
| two-phase | 13 | 0 |
| liquid-side extension | 13 | 0 |
| vapor-side extension | 11 | 0 |
| supercritical | 9 | 0 |

**OVERALL: PASS.**

---

## 13. Isentrope validation

**STATUS: STAGE K — PASS, 2026-08-17.** Important methodology note: ported
from `mixture_isentrope_validation.py` (NOT the base oracle
`mixture_fully_validated.py`) -- that file is byte-identical to the oracle
except for `_verify_isentrope_solution`, the explicit post-hoc (P,s)
residual re-check added earlier this session (the R1234yf-Bug#3-inspired
"insurance" the user explicitly authorized with "Add insurance", and
confirmed via a before/after regression check to introduce zero behavior
change on real output). Reproducing the base oracle's bare-`sol.success`
version instead would mean deliberately regressing an already-validated,
user-approved improvement, so the insurance-patched version is treated as
"the established working model" for this specific piece.

Ported: `compute_isentropes_two_phase` (lever-rule post-processing),
`_isentrope_2eq_residual`, `_verify_isentrope_solution`,
`compute_isentrope_liquid_side` (bubble-row/critical-point-fallback anchor
selection + geomspace outward walk), `compute_isentrope_vapor_side`
(dew-row/critical-point-fallback anchor, incremental entropy ramp,
dome-reentry guard via bubble/dew interpolation, `_adaptive_pressure_walk`
with step halving/doubling both directions).

Validated via `R515B_idaes_package/validate_isentrope_vs_oracle.py`, using
a shared 20-point bubble/dew sweep (255-349K) plus a genuinely-converged
critical point (Tc=381.8939K, same Stage J methodology) so both sides
start from bit-identical inputs, across all 15 Honeywell-chart isentrope
values (0.22-0.49 Btu/lbm-R):

| isentrope family | isentropes/rows | max rel. err |
|---|---|---|
| two-phase (lever rule) | 136 rows across 15 isentropes | 0 |
| liquid-side extension | 15 isentropes, 360 points | 0 |
| vapor-side extension | 15 isentropes, 456 points | 0 |

Confirms the trickiest procedural machinery in the whole oracle --
near-critical anchor selection, the incremental entropy ramp bridging the
dew branch's non-monotonic entropy hump, the dome-reentry guard, and the
adaptive step-halving/doubling pressure walk -- all reproduce exactly,
not just "close." **First test run of this script used an incorrect Btu-
to-J/(mol*K) unit conversion (extra stray /1000, wrong basis for mw_mix)
and got trivially-passing 0-row/0-point results for several families --
caught before accepting the result, fixed to match the oracle's own exact
conversion (`s_btu * BTU_LBMR_TO_JKGK * mw_mix`, line 3962), re-run with
genuine non-degenerate coverage as shown above.** (Recorded per spec rule
74: a self-caught false near-pass, not swept aside.)

**OVERALL: PASS.**

---

## 14. Derivative validation (analytic/symbolic vs. finite-difference)

**STATUS: NOT YET STARTED.** Pending Stage G/H, required before IDAES porting
per spec rule 43 (matching values is not sufficient for an equation-oriented
flowsheet -- derivatives consumed by the active Pyomo model must be validated
too).

---

## 15. Physical-invariant checks

**STATUS: NOT YET STARTED.** Will check x1+x2=1, y1+y2=1, P_L~=P_V,
mu_i^L~=mu_i^V, rho_L>rho_V away from critical collapse, and absence of
NaN/Inf/invalid-composition/wrong-root convergence, alongside each relevant
validation family above rather than as a separate late pass.

---

## 16. Complete numerical comparison dataset

**STATUS: NOT YET CREATED.** Will be a new machine-readable file (CSV/JSON),
referenced from this section once Stage G+ produces comparable states.

---

## 17. Error statistics summary

**STATUS: NOT YET STARTED.**

---

## 18. IDAES structural validation

**STATUS: NOT YET STARTED.** Pending Stage L/M. Will record: base classes used
(`PhysicalParameterBlock`/`StateBlock`/`StateBlockData` per IDAES 2.12.0's
actual custom-property-package architecture -- to be confirmed against
official docs before implementation, not assumed from the native
single-component `HelmholtzParameterBlock` example that `vapor_compression.py`
currently uses), component/phase declarations, state variables, units,
metadata, StateBlock construction, degrees of freedom, initialization,
release, scaling, and absence of unintended free/fixed internal variables.

---

## 19. vapor_compression.py integration validation

**STATUS: NOT YET STARTED.** Pending Stage N.

- Original source (untouched): `DVCT_Project/property_packages_DVRT_code/idaes-hvacr-cycles/vapor_compression.py`
- Integration copy: *(not yet created -- planned name `vapor_compression_r515b_integration.py`, pending final confirmation of target directory convention)*

Architecture facts already extracted directly from the original (read-only,
2026-08-17, see PROJECT_CONTEXT.md for full detail): single
`HelmholtzParameterBlock(pure_component=fluid_name, state_vars=sv,
amount_basis=AmountBasis.MASS)`; `Heater`/`Compressor`/`Heater`/
`PressureChanger` closed-loop flowsheet; state variables actually used are
`flow_mass`, `pressure`, `enth_mass` (PH mode) or `temperature`+`vapor_frac`
(TPX modes); native `temperature_sat`, `eq_complementarity`/`eq_sat` phase
mechanisms; `hp_diagram()`/`pt_diagram()`/`ts_diagram()` plotting methods
called directly on the property block; CoolProp-based initialization (will
need replacement for R-515B, since CoolProp has no R-515B entry).

---

## 20. Key flowsheet-state validation (compressor in/out, condenser, expansion, evaporator)

**STATUS: NOT YET STARTED.** Pending Stage N.

---

## 21. Failure records

Any validation failure encountered during development is appended here in
full (never hidden from statistics), using the format:

```
Validation category:
State:
Reference:
New value:
Absolute error:
Relative error:
Frozen tolerance:
First differing intermediate value:
Diagnosis:
Resolution:
Status:
```

### 2026-08-17: apparent entropy mismatch, `compute_table1_properties` vs. oracle (false start, resolved same session)

```
Validation category: Core (P,h,T,s) validation / Helmholtz kernel validation (Section 4/5)
State: liquid_direct_subcooled (T=300K, rho=10773.989 mol/m^3, x1=0.9385037942)
Reference: oracle _mix_entropy_direct() = 126.9412115 J/(mol*K)
New value: linear_model_codex.compute_table1_properties() s_molar_JmolK = 128.4007115 J/(mol*K)
Absolute error: 1.4595 J/(mol*K) (constant across all 3 tested states, composition-fixed)
Relative error: ~0.0115 (state-dependent, since denominator varies; NOT actually a tolerance violation, see Resolution)
Frozen tolerance: 1e-8 relative (tier-1, Section 1)
First differing intermediate value: none -- alpha_mix, alpha_tau_mix, ar_del_mix all matched
  the oracle's _mix_alpha_and_derivs() output EXACTLY (0.0 absolute difference) when
  mixture_alpha0_alphar_derivs() was called directly; the apparent mismatch only appeared
  when comparing entropy VALUES, not the underlying alpha terms.
Diagnosis: script bug, not a physics/code discrepancy. _mix_entropy_direct() already applies
  ENTROPY_REFERENCE_OFFSET_JMOLK (=-1.4595) internally before returning (confirmed by reading
  its body, lines 2569-2570: `s_raw = R_u*(tau*alpha_tau-alpha); return s_raw + ENTROPY_
  REFERENCE_OFFSET_JMOLK`). The first diagnostic script treated its return value as
  pre-offset "raw" entropy and either (a) compared it directly against a dimensionless s/R
  quantity without multiplying by R_u (a units bug), or (b) added the offset a SECOND time
  when constructing a "rebased" comparison value (a double-counting bug). Both were script
  errors in R515B_idaes_package/validate_table1_vs_oracle.py's draft, not in
  mixture_fully_validated.py, linear_model_codex.py, or compute_table1_properties.
Resolution: corrected validate_table1_vs_oracle.py to back out the oracle's own offset
  (s_ref_raw = s_ref_rebased - ENTROPY_REFERENCE_OFFSET_JMOLK) before comparing against
  compute_table1_properties's never-rebased output. Re-run confirmed exact match (rel_err
  0 to 1.4e-16, machine epsilon) at all 3 tested states -- see Section 5's table.
Status: RESOLVED (script bug fixed, re-verified PASS). Retained here per spec rule 74/89
  (never hide a failure, even a self-inflicted one, and never silently erase it once
  superseded -- this entry stays even though Section 5 now shows PASS).
```

### 2026-08-18: `bell2023_departure_expr` wrong tuple order + missing exp(-delta^l) term (Stage L part 1, resolved same session)

```
Validation category: Stage L Pyomo-kernel validation (new Section 25) -- Helmholtz kernel
  departure-function term specifically
State: liquid_direct (T=300K, rho=10773.988537 mol/m^3, x1=0.9385037942295681)
Reference: r515b_helmholtz_core._mix_alpha_and_derivs's alpha_mix = -4.96182106697 (NumPy path,
  itself already validated exact vs. the oracle in Sections 4-13)
New value (before fix): r515b_pyomo_eos.mixture_alpha_and_derivs_expr's alpha_mix, evaluated via
  pyomo.environ.value() = a differing value giving rel_err=0.00922 at liquid_direct (0.000406 at
  vapor_direct, 3.41e-06 at supercritical_direct)
Absolute error: dep_e=-0.0495730 vs dep_ref=-0.0038450625001861124 (departure term alone, ~13x
  too large in magnitude) -- traced via a component-by-component debug script
Relative error: 0.00922 (liquid_direct, worst case) vs. required tolerance 1e-10 for this
  expression-equivalence check
Frozen tolerance: N/A (this is a Pyomo-vs-NumPy transcription-equivalence check, tolerance
  1e-10, tighter than the tier-1 1e-8 physics tolerance since no solving/rounding is involved --
  only expression evaluation)
First differing intermediate value: dep (the departure-function term). Tred, vred, tau, delta,
  a01, a02, ar1, ar2 all matched the NumPy reference EXACTLY (confirmed via a debug script
  printing every intermediate) -- isolating the bug entirely to `bell2023_departure_expr`.
Diagnosis: `bell2023_departure_expr` in r515b_pyomo_eos.py unpacked BELL_2023_DEP_COEFFS tuples
  as `(nk, dk, tk_pow, _unused)` (order n,d,t,discard-4th) and computed only
  `nk * delta^dk * tau^tk_pow` with NO exponential term. Reading the oracle's own
  `bell2023_departure_alphar`/`bell2023_departure_base` source in linear_model_codex.py (lines
  652-730) showed the correct unpacking is `for nk, tk, dk, lk in BELL_2023_DEP_COEFFS` (order
  n,t,d,l) and every term includes a damping factor: `term = nk * tau^tk * delta^dk *
  exp(-(delta^lk))`. The Pyomo transcription had swapped the t/d exponent roles AND dropped the
  4th coefficient (l) and its exp(-delta^l) factor entirely -- a pure transcription error, not a
  physics or oracle discrepancy (the oracle itself was not touched, read-only per spec rule 9).
Resolution: corrected `bell2023_departure_expr` to unpack `(nk, tk, dk, lk)` and compute
  `nk * (tau**tk) * (delta**dk) * exp(-(delta**lk))` per term, matching the oracle exactly.
  Re-ran validate_pyomo_eos_vs_core.py: all 3 states now PASS at rel_err 0 to 2.09e-16 (machine
  epsilon), well under the 1e-10 tolerance. See Section 25.
Status: RESOLVED (transcription bug fixed in r515b_pyomo_eos.py, re-verified PASS). Retained
  here per spec rule 74/89 (never hide a failure, even a self-inflicted/self-caught one).
```

---

## 22. IDAES architectural decisions log

Decisions affecting the property-package architecture, each with the official
documentation area consulted, the decision, and the reason. Append-only within
this section (earlier decisions are not rewritten if superseded -- a
superseding entry is added instead, dated, explaining why).

### 2026-08-17: Component representation (spec rule 21/22 checkpoint)

- **Question:** does the new IDAES property package represent R-515B as a
  single fixed-composition pseudo-pure fluid, or as two explicit
  flowsheet-visible components (R-1234ze(E) + R-227ea)?
- **Evidence gathered (see PROJECT_CONTEXT.md 2026-08-17 entry for full
  detail):** (1) `vapor_compression.py`'s entire architecture is built around
  IDAES's native single-`pure_component` `HelmholtzParameterBlock` -- no
  component-level flow/composition variable appears anywhere in that file;
  (2) this project's own multi-month framing of R-515B (`mixture_fully_
  validated.py`, `mixture_isentrope_validation.py`, and the earlier,
  superseded `mixture_pseudo_dome.py`) has always solved at a FIXED overall
  composition z1=0.911, exposing a single dome/quality variable, never
  treating x1/y1 as independent flowsheet-level unknowns; (3) this project's
  own 2026-08-12 breadcrumb entry explicitly anticipated integration via
  "the pseudo-pure collapse approach ... feeding a custom property package
  into `vapor_compression.py`-style code."
- **Decision status:** evidence points unambiguously to fixed-composition
  pseudo-pure-fluid representation. NOT yet finalized -- per spec rule 22,
  pausing to get explicit user confirmation before this (foundational,
  expensive-to-reverse) commitment, even though the evidence is one-sided.
  See the message accompanying this entry's creation for the exact question
  posed to the user.
- **CONFIRMED by user, 2026-08-17:** fixed-composition pseudo-pure-fluid
  representation selected (user explicitly chose the recommended option).
  R-515B will be exposed to IDAES as a single-component-equivalent StateBlock
  at z1 frozen to the composition matching w1=0.911, with x1/y1 solved
  internally by the property package's own VLE machinery exactly as the
  oracle does, never as independent flowsheet-level unknowns. This decision
  now governs all subsequent StateBlock/property-package design.

### 2026-08-17: Stage D-F strategy revised after discovering `linear_model_codex.py` is a SHARED dependency, not a separate port target

- **Finding:** reading the oracle's own import block (`mixture_fully_
  validated.py` line 51) confirms it imports `bell2023_Tred_vred`,
  `mixture_alpha0_alphar_derivs`, `load_idaes_helmholtz_json`,
  `mw_from_json`, and other core-kernel functions directly FROM
  `linear_model_codex.py` -- an existing, unmodified file, explicitly
  permitted as a production dependency under spec rule 6 (distinct from the
  oracle-import prohibition of spec rule 9).
- **Decision:** the new production package will import these specific
  functions (`load_idaes_helmholtz_json`, `mw_from_json`, `bell2023_Tred_vred`,
  `mixture_alpha0_alphar_derivs`, `bell2023_departure_alphar`,
  `bell2023_departure_base`, `_bell2023_reducing_derivs_binary`) directly
  from `linear_model_codex.py` rather than re-transcribing Stage D
  (data)/Stage E (reducing functions)/the core-kernel part of Stage F. This
  is stronger than "reproducing the same numbers" -- it is the same code
  path the oracle itself uses, so these layers cannot drift apart. Stage F
  narrows to reproducing the oracle's own DOWNSTREAM assembly of these
  primitives (`_mix_alpha_and_derivs`'s `alpha_mix = a0+x1*ln(x1)+x2*ln(x2)+ar`
  pattern), which is a small, already-verified (Section 4) piece of new
  code, not a full re-port.
- **Caution recorded:** `linear_model_codex.py` ALSO contains a separate,
  independent implementation path (`_mixture_alpha_eval`/
  `compute_table1_properties`) that is NOT what the oracle calls internally.
  This path was verified (Section 5) to numerically match the oracle for
  P, h, and raw entropy at the tested direct states, so it is usable as a
  convenience cross-check, but the new production package's PRIMARY
  implementation should assemble properties from `mixture_alpha0_alphar_
  derivs` directly (mirroring the oracle's own `_mix_alpha_and_derivs`
  pattern) rather than depending on `compute_table1_properties`'s separate
  path, to avoid any risk of the two `linear_model_codex.py` paths
  diverging for a property/state this session did not test (e.g. cv, cp,
  speed of sound, which use finite-difference second derivatives in
  `compute_table1_properties` and have no oracle equivalent to check against
  at all).

### 2026-08-18: PAUSE -- Stage L part 3 two-phase VLE encoding is a genuine rule-22 architectural fork

**Context:** Stage L parts 1-2 (native-Pyomo EOS kernel and P/h/s/g/mu
expressions, `r515b_pyomo_eos.py`, Sections 24a/24b) are complete and
validated. Part 3 requires the actual `PhysicalParameterBlock`/
`StateBlockData` classes. Re-inspected (per rule 22's explicit
requirement) `vapor_compression.py`'s real interface contract (staged
from the user's device) and the installed `general_helmholtz.
helmholtz_state.HelmholtzStateBlockData` source (the package
`vapor_compression.py` currently uses for pure fluids) to confirm the
exact drop-in contract our package must satisfy.

**Confirmed by direct source inspection:** for `StateVars.PH` +
`AmountBasis.MASS`, exactly 3 Vars are real DOF-carrying state variables
per stream -- `flow_mass` (extensive), `pressure`, `enth_mass`
(intensive). `temperature` and `vapor_frac` are DERIVED (Expression, not
independently fixable) in this mode -- confirmed both from
`HelmholtzStateBlockData`'s own branch for this mode (line ~490-512) and
from `vapor_compression.py`'s own usage (only `flow_mass`/`pressure`/
`enth_mass` are ever `.fix()`ed in `Mode.PH` branches; `vapor_frac` is
never touched in PH mode; temperature bounds go through separate
flowsheet-level `Tmin`/`Tmax` Params and `T_lower_bound`/`T_upper_bound`
Constraints, not direct Var bounds).

**The fork:** for a pure fluid, phase-region logic is a 1-DOF saturation
curve plus a smooth complementarity condition (`eq_complementarity`/
`eq_sat`) -- no composition to split. R-515B is a confirmed ZEOTROPIC
BLEND (bubble T != dew T at fixed P -- "glide" -- and vapor composition
y1 != liquid/feed composition x1 in the two-phase region; this is exactly
what Stage I's real 3-equation `solve_bubble_at_t`/`solve_dew_at_t` VLE
solves determine). Encoding this inside the active Pyomo NLP (not as a
SciPy pre-solve, per spec rule 44) is a genuinely open, multi-way design
problem. Four materially different, individually defensible options were
identified (full detail, tradeoffs, and risks for each recorded in
`PROJECT_CONTEXT.md`'s matching 2026-08-18 entry, not duplicated here to
avoid drift between the two records):
  (a) single-phase-only initial scope (defer two-phase to a follow-up)
  (b) complementarity-based two-phase from the start, general_helmholtz-
      style, extended with a composition-split (fugacity-equality) system
  (c) GDP (disjunctive) three-region formulation via `pyomo.gdp`
  (d) externally-fixed phase state from `initialize()`, using the already-
      validated Stage I/J SciPy solvers to decide the region per solve

**Status: PAUSED, per spec rule 22** ("if more than one materially
different design remains plausible... PAUSE AND CHECK WITH THE USER
before committing, since IDAES structural choices are expensive to
reverse once flowsheet integration (Stage N) depends on them") -- this is
the exact condition. Asked the user via the interactive question tool
which direction to take. Stages A-K and Stage L parts 1-2 remain complete,
validated, and synced to the device regardless of the outcome; no work is
at risk from this pause.

### 2026-08-18: RESOLVED -- fork collapses to x1=y1=z1 pseudo-pure/near-azeotropic saturation treatment (user correction)

The user responded to the pause not by selecting one of the 4 listed
options but by correcting the framing itself: R-515B's Stage L
phase-equilibrium treatment should use x1=y1=z1 (fixed), and the
"zeotropic" characterization used to pose the fork overstated the
severity of the real physics. Re-verified against primary sources before
acting (per the user's explicit instruction, not just taking the
correction at face value):
- Re-read `solve_bubble_at_t`/`solve_dew_at_t` in the oracle in full:
  confirmed the oracle's general-purpose bubble/dew solver DOES carry
  x1(liquid)/y1(vapor) as genuinely distinct solved unknowns -- that part
  of Section 8/9's characterization of the oracle's own math is accurate
  and unchanged.
- Searched the FULL `PROJECT_CONTEXT.md` history (predating this Master
  Task, from 2026-03-03) for R-515B's established classification: a
  dedicated fixed-pressure glide test against the real Honeywell datasheet
  grid classified it **"Near-azeotropic"** (max dT_glide=0.2677K,
  mean=0.0291K; max composition split |delta_x1|~5.3%, mean~0.9%) -- NOT
  a strict zero-glide azeotrope, but also not the severity implied by
  calling it "zeotropic" outright. That same history shows a "pseudo-pure
  saturation mode" (fixed single composition on both phases) was already
  prototyped and quantified as acceptable (MAE=0.373K, max=0.906K vs.
  Honeywell) well before this Master Task began.
- Conclusion: my error was in the Stage L PROPERTY-PACKAGE design framing
  (overstating what machinery the active NLP needs), not in the earlier
  factual description of the oracle's own general-purpose solver
  capability (which does genuinely support x1!=y1 and remains correctly
  described that way in Sections 8/9 above).

**Design now adopted (supersedes options a-d above):** the active
two-phase logic in the StateBlockData uses x1=y1=z1 fixed -- structurally
identical to how `general_helmholtz.HelmholtzStateBlockData` handles a
pure fluid's saturation curve (single composition, P_l=P_v, g_l=g_v, no
composition-split unknown in the NLP at all). The rigorous bubble/dew-
with-real-y1 machinery in `r515b_helmholtz_core.py` (Stage I, Sections
8/9) is unaffected and remains available for `initialize()` and as the
deeper reference for quantifying the x=y simplification's own error once
Stage L is built (mirroring the historical `pseudopure_mode_eval.py`
MAE/max-error reporting pattern, but against the new validated core).

**Status: RESOLVED, proceeding.** This was the rule-22 check-in; it is not
being re-opened. See `PROJECT_CONTEXT.md`'s matching 2026-08-18 entry for
full detail.

### 2026-08-18: pseudo-pure (x1=y1=z1) saturation solve built and sanity-checked

**Built:** `r515b_helmholtz_core.solve_pseudopure_saturation_at_t` -- a NEW
function (not a port of any oracle function; the oracle deliberately never
implements this simplification). At fixed T and the fixed pseudo-pure
composition z1, solves 2 equations (P_l=P_v, and g_l=g_v -- the Maxwell/
equal-molar-Gibbs-energy condition, which is the mathematically correct
second condition once composition is fixed identical on both phases,
rather than the separate mu1/mu2 equalities the rigorous 3-unknown bubble/
dew solve needs; full derivation in the function's module-level comment)
for 2 unknowns (rho_l, rho_v). Reuses the same sigmoid rho-parameterization,
4-seed retry ladder, and tapered density-separation gate as Stage I's
`solve_bubble_at_t`/`solve_dew_at_t` for consistency and numerical
robustness.

**Sanity-checked** (`R515B_idaes_package/validate_pseudopure_saturation.py`,
not a pass/fail-vs-oracle validation since this is a deliberate
simplification with no oracle equivalent to match -- a well-behavedness
and deviation-quantification check instead):
- Converges 15/15 across a 250-375K test grid (covering the full practical
  refrigeration-cycle range).
- At every converged point, the pseudo-pure saturation pressure sits
  strictly between the true bubble pressure and true dew pressure at that
  T (15/15) -- confirms it behaves as a sane "average" of the real
  liquid/vapor split, not an unrelated root.
- Deviation from the true bubble/dew average: **mean 0.0009%, max
  0.0032%** -- i.e. essentially negligible, strongly corroborating the
  2026-08-18 near-azeotropic classification (this project's own historical
  Honeywell-grid glide test found max glide only ~0.27K / mean ~0.03K) and
  validating that the x1=y1=z1 simplification is an excellent fit for this
  specific blend, not merely an expedient shortcut.

**Not yet done (at the time of that entry):** the corresponding native-
Pyomo saturation Constraint set. See the immediately following entry --
now built and validated.

### 2026-08-18: native-Pyomo saturation-curve system built and validated (real IPOPT solve, not just value() evaluation)

**Design research first:** read the installed `general_helmholtz.
helmholtz_state.HelmholtzStateBlockData._state_vars` source in full to see
exactly how it handles PH-mode temperature/vapor_frac. Finding: for PH
mode, `temperature` and `vapor_frac` are NOT solved via a visible Pyomo
complementarity Constraint at all -- they are single black-box, globally-
smooth EXTERNAL FUNCTION calls (`t_hp_func(cmp, h, P)`, `vf_hp_func(cmp, h,
P)`) backed by a compiled C++ library with baked-in saturation
correlations and phase-blend smoothing invisible to Pyomo. (The
`eq_complementarity`/`eq_sat` Constraints found earlier belong to
`StateVars.TPX` mode specifically, a different, unrelated code path.) We
have no such compiled correlation library and building one is out of
scope (spec rule 1) -- so the native-Pyomo path must solve the defining
saturation EQUATIONS directly and implicitly as ordinary Constraints,
which is in any case more directly rule-44-compliant (no black box at
all, not even a compiled one) at the cost of needing a real implicit NLP
solve rather than a single smooth expression evaluation.

**Built:** `r515b_pyomo_eos.saturation_residuals_expr(d1, d2, x1, T_sat,
rho_l_sat, rho_v_sat, P, Tc1, Tc2, vc1, vc2)` -- returns 3 residual
expressions (res_P_liq, res_P_vap, res_gibbs) built from
`mixture_state_expr` (Stage L part 2, already validated), for use as 3
Constraints against 3 free Vars (T_sat, rho_l_sat, rho_v_sat) at a given
fixed pressure P. This is the SciPy pseudo-pure solver's exact
counterpart, transcribed as an implicit native-Pyomo system instead of a
`scipy.optimize.least_squares` call.

**Validated** (`R515B_idaes_package/validate_pyomo_saturation_vs_core.py`)
by actually SOLVING the 3x3 system with IDAES's own `get_solver()`
(IPOPT) -- not just evaluating expressions at fixed values, since this is
the first Stage L component that must genuinely be solved as an implicit
system. Fixed pressure at 25 target values from a fine (5K-step)
continuation sweep of the already-validated SciPy pseudo-pure solver
(255-375K), with the Pyomo model's own initial guesses deliberately
OFFSET from the SciPy answer (1.05x/1.15x/0.7x) so this is a genuine
solver-convergence check, not a value() echo. **Result: 25/25 PASS**, all
at machine epsilon (rel_err 0 to 1.4e-12).

**Bug found and fixed (self-caught, methodology bug in the validation
script, not in the Pyomo formulation) -- recorded per spec rule 74/89:**
a first draft used a COARSE target list (260/290/320/350/370K, i.e. up to
30K seed jumps between comparison points). The 350K point FAILED
(rel_err up to 30% on rho_l). Root-cause diagnosis: printed a fine
continuation sweep of the SciPy reference solver alone (5K steps, each
reseeded from the previous point, exactly as `solve_bubble_at_t` and this
project's other iterative solvers are meant to be used) and found rho_l
decreases smoothly and monotonically with T (9189 at 320K -> 8017.6 at
350K -> 6767.4 at 370K) -- but the COARSE script's 30K-jump-seeded SciPy
call at 350K had converged to rho_l=6092.99, which is LOWER than even the
370K value: non-monotonic, i.e. a spurious/non-physical root, despite
passing its own internal residual gate (r_P, r_g both <1e-14). The IPOPT/
Pyomo solution at that same badly-seeded comparison point (rho_l=7942.6)
was actually CLOSER to the true continuation-consistent value than the
badly-seeded SciPy reference was -- i.e. this was a seed-robustness gap in
how the COMPARISON was set up (a 30K jump is simply too large for this
solver family's local Newton-type convergence basin near this quantity's
curvature), not a defect in either the Pyomo formulation or
`solve_pseudopure_saturation_at_t` itself when used as intended (fine-step
continuation, matching every other iterative solver in this codebase).
Fixed by switching the validation script to a 5K-step continuation grid;
re-ran and got the clean 25/25 machine-epsilon result above. Flagging this
seed-sensitivity as a general characteristic of the saturation system (not
just a script bug) worth remembering for `initialize()` design: any caller
of `solve_pseudopure_saturation_at_t` (SciPy) or the Pyomo saturation
Constraints should seed from a NEARBY known state (fine continuation, or
the ancillary initial-guess machinery in Section 23), not an arbitrary or
far-away guess.

**Not yet done (at the time of that entry):** the full smooth single-phase/
two-phase blending logic. See the immediately following entry -- now
built and validated.

### 2026-08-18: native-Pyomo smooth PH-flash (single-phase/two-phase blending) built and validated

**Built:** `r515b_pyomo_eos.mixture_ph_flash_residuals_expr` -- an 8-
equation, 8-unknown smooth system (T_sat/rho_l_sat/rho_v_sat from the
saturation curve above; vapor_frac; T_liq/rho_liq and T_vap/rho_vap, two
ALWAYS-solved single-phase branches) that determines, for given (P,H),
whether the state is subcooled liquid, two-phase, or superheated vapor --
with NO explicit if/then branching, so the whole system stays smooth and
IPOPT-differentiable. Design, in full:
- A smooth complementarity constraint (via `idaes.core.util.math.
  smooth_max`) pins vapor_frac to 0 when H is below the saturated-liquid
  enthalpy and to 1 when H is above the saturated-vapor enthalpy, mirroring
  `general_helmholtz`'s own `eq_complementarity` pattern (adapted from
  their pressure-vs-P_sat(T) comparison to an enthalpy-vs-h_sat(P)
  comparison, since we have no p_sat_t_func-direction correlation).
- Each single-phase branch (liquid, vapor) is fed a CLIPPED target
  enthalpy (via the same smooth_max trick) so it never has to represent an
  enthalpy outside its own real physical domain -- outside its true
  region it simply sits exactly at the saturation boundary (T_liq=T_sat,
  rho_liq=rho_l_sat, or the vapor equivalent), never extrapolating into a
  nonphysical branch. This is a deliberate design choice, different in
  detail from (but in the same spirit as) `general_helmholtz`'s own
  `pressure_phase` extension trick, adapted because we lack their compiled
  correlation library.
- Any extensive-like property is then recovered via Y_actual =
  (1-vapor_frac)*Y(T_liq,rho_liq) + vapor_frac*Y(T_vap,rho_vap) -- proven
  correct in all three regions in the function's own docstring derivation
  (outside the dome, one branch carries zero weight and the other holds
  the real single-phase answer; inside the dome, both branches sit at the
  boundary and the blend reduces to the correct lever-rule mixture value).

**Validated** (`R515B_idaes_package/validate_pyomo_flash_vs_core.py`) via
a real IPOPT solve of the full 8x8 system at 4 test points spanning all
three regions (subcooled liquid, two two-phase qualities 0.3/0.7,
superheated vapor), fixed P/H, against the new explicit-branching SciPy
reference `r515b_helmholtz_core.flash_ph_pseudopure`. **Result: 4/4 PASS**
after two self-caught bugs (both recorded below, neither hidden):

1. **Complementarity sign/pairing bug.** First draft paired
   `vapor_frac*h_over_sat - (1-vapor_frac)*h_under_sat` (h_over_sat =
   superheated indicator, h_under_sat = subcooled indicator). Result: an
   EXACT 1-vapor_frac inversion at all 4 test points (subcooled gave
   vf=1, superheated gave vf=0, the two two-phase points gave 0.7/0.3
   instead of 0.3/0.7) -- a clean, immediately diagnosable signature.
   Root cause: `general_helmholtz`'s own pairing is [subcooled indicator]
   *vf - [superheated indicator]*(1-vf); the first draft paired the
   indicators to vf/(1-vf) backwards relative to that reference pattern.
   Fixed by swapping to `vapor_frac*h_under_sat - (1-vapor_frac)*
   h_over_sat`; re-verified the sign logic by hand for both the subcooled
   and superheated limiting cases (worked through explicitly in the
   function's own updated code comment) before re-running.
2. **Seed-basin non-uniqueness (methodology, not a formulation bug).**
   After fixing (1), vapor_frac matched the reference exactly at all 4
   points, but T_actual was wrong at 3 of 4 (off by 1-7%). Diagnosis: the
   single-phase branch system (2 equations: P(T,rho)=P, h(T,rho)=
   h_target; 2 unknowns T,rho) is NOT globally unique away from a good
   seed -- with the validation script's first, generic/flat initial
   guesses (T_liq=295K constant, rho_liq at a generic ballpark density
   unrelated to the also-being-solved rho_l_sat), IPOPT converged to a
   mathematically valid (residuals ~1e-9, `termination_condition=optimal`)
   but PHYSICALLY WRONG root (e.g. rho_liq=3756 mol/m^3, sitting in the
   unstable/metastable region between the true liquid and vapor branches,
   instead of the true rho_liq=10318.6). This is the same class of
   seed-sensitivity already documented for the saturation curve above --
   confirmed, not a one-off, now seen in a second, independent piece of
   this system. Fixed by seeding T_liq/rho_liq/T_vap/rho_vap near the
   saturation-state anchor (perturbed 3-15%, branch-consistent: liquid
   guesses biased toward the liquid branch, vapor guesses toward the
   vapor branch) rather than an arbitrary flat constant. Re-ran: all 4
   points PASS at rel_err ~1e-9. **Design implication for `initialize()`:**
   any real StateBlockData built on this formulation MUST seed T_liq/
   rho_liq/T_vap/rho_vap/T_sat/rho_l_sat/rho_v_sat from a good starting
   point (e.g. `flash_ph_pseudopure`'s own converged SciPy answer, or
   continuation from a neighboring flowsheet state) -- this is not
   optional robustness polish, it is required for correctness, since a
   bad seed can converge to a numerically-valid-looking but physically
   wrong answer with no local signal (from residuals or termination
   status alone) that anything is amiss.

**Not yet done:** the actual `PhysicalParameterBlock`/`StateBlockData`
classes wiring this validated 8-equation system into a real IDAES
property package (metadata, PH+MASS state vars matching
`vapor_compression.py`'s confirmed contract, ports, ties to `flow_mass`,
scaling, and an `initialize()` method that seeds every auxiliary Var from
`flash_ph_pseudopure`/`solve_pseudopure_saturation_at_t` per the design
implication just above).

### 2026-08-18: user clarification -- R-515B is NOT a zeotrope, oracle is ground truth; "no silent ancillary correlations" made an explicit standing commitment

Mid-build, the user sent a message reinforcing (not revising) the prior
rule-22 resolution: **R-515B is not a zeotrope; `mixture_fully_validated.py`
is ground truth.** It then drew a distinction material to how the
StateBlock's saturation/flash logic must be built: our R-515B pipeline
solves P_sat/rho_L_sat/rho_V_sat (and x=y) from the FULL binary Helmholtz/
VLE model via the bubble/dew equilibrium solves -- never from separate
fitted ancillary correlations P_sat=f(T), rho_L_sat=g(T), rho_V_sat=h(T).
This is unlike how `general_helmholtz` treats pure fluids (compiled
ancillary-equation external functions, confirmed earlier this stage) --
R-515B has no such ancillary equations of its own, and the Honeywell
datasheet's saturation points are VALIDATION data, not correlations
available to embed. If the StateBlock needs a smooth function like
`temperature_sat(P)`, three conceptually different options exist: (1)
preserve the validated VLE equations algebraically inside the IDAES
formulation; (2) derive some other IDAES-compatible representation
directly from the validated Helmholtz model; (3) fit new ancillary
correlations from validated VLE data. Option 3 is explicitly flagged as
NEW METHODOLOGY that must never be introduced silently just because it
would make the StateBlock easier -- it requires demonstrating necessity
and pausing for approval first.

**Verification performed** (not accepted at face value): `grep -n
"ancillary" r515b_property_package.py r515b_pyomo_eos.py` -> zero matches.
The string "ancillary" only appears in the pre-existing Stage-I files
`ancillary_initial_guess.py`/`validate_ancillary_guess.py` (Section 23,
below), which generate NUMERICAL SEED VALUES only and are never imported
by either the Pyomo EOS kernel or the property package. Confirmed
`r515b_property_package.py` imports exactly `r515b_helmholtz_core` (for
`flash_ph_pseudopure`, used ONLY inside `_R515BStateBlock.initialize()` as
a seed generator, per the already-documented rule-44 carve-out) and
`r515b_pyomo_eos.mixture_ph_flash_residuals_expr` (the actual governing
Constraints -- Option 1). **Conclusion: the property package as already
built already IS Option 1, with no drift toward Option 3 anywhere.** No
code change was required in response to this message. Standing commitment
recorded here: this project will not introduce fitted R-515B ancillary
correlations at any point without first demonstrating necessity and
pausing for explicit user approval.

---

## 23. Ancillary initial-guess machinery (initialization support, Stage L prep)

**STATUS: BUILT AND VALIDATED, 2026-08-17.** Explicitly requested by the
user as a narrower side task ("build the machinery please") during a pause
in the broader Stage D-K porting work, motivated by the question "do we have
ancillary relations for R-515B?" This is NOT part of the numbered Stage A-O
sequence's critical path -- it is initialization-support machinery intended
for the new package's `initialize()` routine (Stage L/spec rule 47), so a
cold-start IPOPT solve has a physically-informed starting density instead of
an arbitrary guess.

**What exists:**
- `R515B_idaes_package/ancillary_initial_guess.py` -- new file. Composes
  each PURE fluid's own ancillary saturated-density correlation (IDAES's own
  `delta_l_sat_approx`/`delta_v_sat_approx`, types 1/2 respectively --
  formulas confirmed directly from `idaes/models/properties/general_helmholtz/
  expressions/sat_delta_approx.py`, NOT guessed) through this project's own
  tau_i=Tc_i/T chain-rule mapping, then combines the two pure-fluid
  saturated-volume estimates via simple ideal volume-additivity mixing
  (v_mix = z1*v1_sat + z2*v2_sat) to produce a fast, non-iterative
  (rho_l_guess, rho_v_guess) estimate at any (T, z1). This combination rule
  is NEW engineering (not present in the oracle, which has no
  mixture-level ancillary at all -- confirmed by reading `mixture_fully_
  validated.py`, `mixture_isentrope_validation.py`, and `linear_model_
  codex.py` in full; every oracle bubble/dew point comes from a real
  iterative 3-equation VLE solve).
- `R515B_idaes_package/validate_ancillary_guess.py` -- new validation
  script. Compares the composite guess against the oracle's own real,
  converged `solve_bubble_at_t`/`solve_dew_at_t` densities (read-only
  reference use, spec rule 9 carve-out) across 12 temperatures spanning
  250K to Tc_mix-15K (using the oracle's own mixture critical-point solve,
  Tc_mix=381.52K, to set the safe upper bound).

**Validation results** (`R515B_idaes_package/ancillary_guess_validation_results.json`):

| branch | max rel. error | mean rel. error | characterization |
|---|---|---|---|
| liquid (rho_l) | 2.45% | 2.08% | consistently good across the whole tested range -- a solid cold-start seed |
| vapor (rho_v) | 83.8% (at T=250K, far from Tc) | 46.1% | POOR at low T, improving monotonically to ~5% near Tc |

**Diagnosis of the vapor-branch weakness (not a bug -- an inherent
limitation of the combination rule as built):** the composite guess uses
the OVERALL feed composition z1 for both the liquid and vapor branch
estimates. The real vapor phase is significantly enriched in the more
volatile component relative to z1 away from the critical point (y1 != z1,
sometimes by a large margin) -- exactly the composition split the real
3-equation VLE solve exists to determine. Because the guess ignores this
split entirely, and vapor molar volume is far more sensitive to
composition than liquid molar volume is (dilute-gas behavior), the vapor
density guess degrades badly at low T where the true y1/z1 gap is largest,
and improves near Tc where liquid and vapor compositions converge together
with the phases themselves.

**Recommendation for Stage L's initialize() use:** use this composite guess
directly for the LIQUID density seed (validated good, ~2% typical). For the
VAPOR density seed, do NOT use this composite as-is; instead fall back to
the oracle's own existing simple heuristic (rho_v_seed0 = 0.01*rho_l_seed,
already used throughout `establish_reference_tolerances.py` and mirrored
from the oracle's `run_true_vle_envelope`) or, as a future improvement,
extend the composite with a relative-volatility estimate to approximate the
y1/z1 split before combining -- not attempted this session; flagged as a
possible enhancement, not a blocker, since the existing 1%-of-liquid
heuristic is already an accepted, working fallback.

**Not yet done:** no failure -- this machinery is validated for its stated
purpose (seed quality, not accuracy as an independent model) and the
limitation above is a characterization, not an open defect. Wiring it into
Stage L's actual `initialize()` method has not started yet (Stage L itself
has not started).

---

## 24a. Stage L (part 1): Pyomo EOS-kernel validation

**STATUS: BUILT AND VALIDATED, 2026-08-18.** First piece of Stage L (the
actual IDAES property package). Per MASTER TASK spec rule 44 ("no hidden
SciPy solves inside the active IDAES NLP"), the property package's
Constraints must be built from native Pyomo expressions (so IPOPT's own
automatic differentiation applies), NOT by calling
`r515b_helmholtz_core.py`'s plain NumPy functions from inside a Constraint
body. `r515b_helmholtz_core.py` itself is NOT superseded -- it remains the
correct tool for `initialize()` (seeding via its validated SciPy bubble/
dew/critical-point solvers) and for ongoing regression.

**What exists:**
- `R515B_idaes_package/r515b_pyomo_eos.py` -- new file. Native-Pyomo
  (`pyomo.environ.exp`/`log`) transcription of: `pure_alpha0_expr` (phi==1
  ideal-gas branch, the only branch both pure fluids use, confirmed by
  direct JSON inspection), `pure_alphar_expr` (phi==2 residual branch,
  same scope narrowing), `bell2023_Tred_vred_expr` (reducing functions),
  `bell2023_departure_expr` (Bell 2023 3-term departure function), and
  `mixture_alpha_and_derivs_expr` (full assembly: tau, delta, alpha_mix,
  Tred, vred from x1/T/rho, which may be Pyomo Vars). No hand-coded
  derivative companions are returned anywhere in this file -- Pyomo/IPOPT
  differentiates the returned expressions automatically; this is
  deliberate (see file docstring) to avoid a second hand-written
  derivative path silently drifting from the primal expression.
- `R515B_idaes_package/validate_pyomo_eos_vs_core.py` -- new validation
  script. Builds a Pyomo `ConcreteModel` with `Var T, rho, x1` FIXED (not
  solved) at the same 3 representative states used throughout Sections
  5-7 (liquid_direct, vapor_direct, supercritical_direct), evaluates
  `mixture_alpha_and_derivs_expr`'s `alpha_mix` via `pyomo.environ.value()`,
  and compares against `r515b_helmholtz_core._mix_alpha_and_derivs`'s own
  NumPy `alpha` output at the identical inputs. Pure expression-evaluation
  check -- no solving happens, so this isolates transcription correctness
  only.

**Validation results** (tolerance 1e-10, tighter than the tier-1 1e-8
physics tolerance since this is an expression-equivalence check, not a
physics comparison):

| state | T (K) | rho (mol/m^3) | rel_err (alpha_mix) | status |
|---|---|---|---|---|
| liquid_direct | 300.0 | 10773.988537 | 1.79e-16 | PASS |
| vapor_direct | 300.0 | 124.343493 | 0.0 | PASS |
| supercritical_direct | 389.5 | 3000.0 | 2.09e-16 | PASS |

All 3 states PASS at machine epsilon. **One bug found and fixed during
this validation** -- see the 2026-08-18 entry in Section 21 (Failure
records) for full detail: `bell2023_departure_expr` initially unpacked
`BELL_2023_DEP_COEFFS` in the wrong tuple order `(n,d,t,unused)` instead
of the oracle's actual `(n,t,d,l)`, and omitted the `exp(-delta^l)`
damping factor entirely. Isolated via a component-by-component debug
script (Tred/vred/tau/delta/a01/a02/ar1/ar2 all matched exactly; only the
departure term `dep` diverged, ~13x too large), fixed by matching
`bell2023_departure_alphar`/`bell2023_departure_base`'s exact term formula,
and re-verified PASS.

**Not yet done (Stage L, part 2 onward):**
- P/h/s/g Pyomo expressions (from alpha_mix and Pyomo-AD tau/delta
  derivatives via `pyomo.core.expr.calculus.derivatives.differentiate` or
  equivalent)
- fugacity/chemical-potential Pyomo expressions (composition derivatives
  of alpha_mix)
- the actual `PhysicalParameterBlock`/`StateBlockData` classes wiring
  these expressions into real Constraints, with PH state variables
  (matching `vapor_compression.py`'s PH mode), phase-equilibrium
  constraints (P_L=P_V, mu_i_L=mu_i_V), and an `initialize()` method
  calling `r515b_helmholtz_core.py`'s validated SciPy solvers
- metadata, units, scaling

(All 4 bullets above except the last two -- P/h/s/g and fugacity/chemical
potential -- are now DONE; see Section 24b immediately below.)

---

## 24b. Stage L (part 2): Pyomo P/h/s/g/Z and mu1/mu2 validation

**STATUS: BUILT AND VALIDATED, 2026-08-18.** Extends Stage L part 1
(Section 24a)'s Pyomo EOS kernel with the remaining state-property and
fugacity/chemical-potential expressions.

**Design decision -- symbolic differentiation instead of hand-transcribed
derivative formulas:** rather than porting the oracle's hand-derived
analytic tau/delta-derivative term sums (`mixture_alpha0_alphar_derivs`'s
`a0_tau`/`ar_tau`/`ar_del` accumulation loops) or its composition-derivative
formulas (`_bell2023_reducing_derivs_binary`, `bell2023_departure_base`/
`_alphar`'s tau/delta-derivative outputs, `chemical_potentials_analytic`'s
`dar_dx1`/`dar_drho` assembly), `r515b_pyomo_eos.py`'s new
`mixture_state_expr`/`mixture_chemical_potentials_expr` use Pyomo's own
exact symbolic differentiation
(`pyomo.core.expr.calculus.derivatives.differentiate`, mode=
`reverse_symbolic`) directly on the already-built `alpha_mix`/`ar_mix`
Pyomo expression trees, differentiating with respect to the actual T, rho,
x1 Pyomo Vars. This was a deliberate choice, not a shortcut: a Helmholtz
EOS's alpha(tau,delta,x) is by construction a genuine closed-form algebraic
function, so its exact partial derivatives are well-defined mathematical
objects that symbolic differentiation computes exactly (verified: it is
NOT a finite-difference approximation -- confirmed via a standalone test
comparing `differentiate(..., mode=reverse_symbolic)` output against a
central-difference check on a toy expression, matching to ~1e-10, the
residual being pure FD truncation error on the *reference* side, not error
in the symbolic result). Using it here removes an entire class of
hand-transcription bugs -- the exact class that caused the departure-
function bug in Section 24a/21 -- at the cost of one extra `differentiate`
call per partial derivative needed. The chain-rule identities used to
convert d(alpha)/dT, d(ar)/drho, d(ar)/dx1 into the physically-meaningful
alpha_tau_mix/ar_del_mix/dar_dx1 quantities are derived in full, with
citations, in each function's own docstring in `r515b_pyomo_eos.py`.

**What exists:**
- `mixture_state_expr` -- returns P_pa, h_jmol, g_jmol, s_jmolk (raw, pre-
  Honeywell-rebase, matching `mix_entropy_direct`'s internal convention;
  an `apply_entropy_offset` flag adds `ENTROPY_REFERENCE_OFFSET_JMOLK` when
  a rebased value is wanted), Z, plus tau/delta/Tred/vred/alpha_mix/
  alpha_tau_mix/ar_del_mix for callers that need the intermediate
  quantities directly.
- `mixture_chemical_potentials_expr` -- returns (mu1_jmol, mu2_jmol) via
  the standard binary-mixture fugacity-coefficient identity for a molar-
  Helmholtz-explicit EOS parameterized by (T, rho, x1): d(n*ar)/dn1 =
  ar_mix + x2*(d ar/d x1)|_{T,rho} + rho*(d ar/d rho)|_{T,x1}, and
  symmetrically for n2 with a -x1 term instead of +x2 (full derivation in
  the function's docstring); f_i = x_i*rho*R*T*exp(d(n*ar)/dn_i),
  mu_i = R*T*ln(f_i) -- identical formula structure to the oracle's own
  `chemical_potentials_analytic`.
- `R515B_idaes_package/validate_pyomo_state_vs_core.py` -- new validation
  script, same pattern as Stage L part 1's validator (fixed, non-solved
  Pyomo Vars at the 3 representative states; compares `value()` of the new
  expressions against `r515b_helmholtz_core.py`'s `mix_state`/
  `mix_entropy_direct`/`chemical_potentials_analytic`).

**Validation results** (tolerance 1e-8; all 6 quantities at all 3 states
actually landed at machine epsilon, 0 to 1.9e-15):

| state | P_pa | h_jmol | g_jmol | s_jmolk (raw) | mu1_jmol | mu2_jmol |
|---|---|---|---|---|---|---|
| liquid_direct | 1.85e-15 | 5.0e-16 | 5.81e-16 | 4.43e-16 | 2.06e-16 | 1.2e-16 |
| vapor_direct | 0 | 1.54e-16 | 0 | 1.4e-16 | 0 | 0 |
| supercritical_direct | 1.06e-15 | 2.94e-16 | 0 | 1.46e-16 | 1.53e-16 | 1.89e-16 |

All 18 checks PASS. No bug found this time (unlike Stage L part 1) --
the symbolic-differentiation approach worked correctly on first attempt,
consistent with the design rationale above (removing the class of bug
that hand-transcription is prone to). Re-ran `validate_pyomo_eos_vs_core.py`
(Stage L part 1) immediately after this change to confirm the
`_mixture_ideal_residual_parts_expr` refactor introduced no regression --
still PASS, rel_err unchanged (0 to 2.09e-16).

**Not yet done:**
- the actual `PhysicalParameterBlock`/`StateBlockData` classes wiring
  these expressions into real Constraints, with PH state variables
  (matching `vapor_compression.py`'s PH mode), phase-equilibrium
  constraints (P_L=P_V, mu_i_L=mu_i_V), and an `initialize()` method
  calling `r515b_helmholtz_core.py`'s validated SciPy solvers
- metadata, units, scaling
- a numerical-safety guard on `log(f_i)` analogous to the oracle's
  `max(1e-300, ...)` clamp (a scaling/robustness concern for the eventual
  Constraint, not a kernel-correctness issue -- flagged in
  `r515b_pyomo_eos.py`'s own STATUS block so it isn't forgotten)

---

## 24c. Stage L (part 3, final): `R515BParameterBlock`/`R515BStateBlockData` construction, DOF, and initialize()+solve validation

**STATUS: BUILT AND VALIDATED, 2026-08-18.** `r515b_property_package.py`
wires the validated `mixture_ph_flash_residuals_expr` 8-equation system
(Section 24b/PH-flash entry above) into a real IDAES `PhysicalParameterBlock`
/`StateBlockData` pair matching `vapor_compression.py`'s confirmed PH+MASS
drop-in contract (3 public state vars: flow_mass, pressure, enth_mass;
`self.Mix = Phase()` single-phase presentation). `validate_state_block_
construction.py` builds a real `FlowsheetBlock`+`R515BParameterBlock`+
`R515BStateBlock`, checks degrees of freedom before/after fixing the 3
state vars, then runs `.initialize()` (SciPy-seeded per the rule-44
carve-out) followed by the block's own embedded IPOPT solve, and compares
against `flash_ph_pseudopure` at the same 4 representative points used
throughout Stage L part 3.

**Result: 4/4 PASS.** T rel_err 4.9e-10 to 1.9e-09; vapor_frac abs_err
3.9e-09 to 6.2e-08; dens_mass rel_err 6.2e-10 to 1.6e-07 -- all far inside
the 1e-6 working tolerance.

**One self-caught TEST-SCRIPT bug (not a package bug):** a first draft
asserted `degrees_of_freedom(blk) == 3` before fixing state vars (naive
total_vars(11) - total_constraints(8) reasoning). Measured value was 2.
Root cause: IDAES's `degrees_of_freedom()` = `number_unfixed_variables_
in_activated_equalities(block) - number_activated_equalities(block)` --
it only counts variables that actually appear inside an activated
equality. `flow_mass` never appears in any of the 8 embedded flash
constraints (correctly -- it is a pure extensive throughput quantity,
independent of the intensive T/P/h/rho system), so it's excluded from the
count regardless of fixed/unfixed status. Verified this is the SAME
convention the reference pure-fluid `HelmholtzParameterBlock`
(`StateVars.PH`, `AmountBasis.MASS`) uses: a direct check shows that block
reports `degrees_of_freedom()==0` both before and after fixing state vars,
since it has ZERO Pyomo-visible constraints at all (temperature/
vapor_frac come from black-box external functions there). Our block
correctly reports DOF=2 before fixing (pressure, enth_mass counted;
flow_mass not) and DOF=0 after fixing all 3 public state vars. Fixed the
test script's expectation, not the package -- documented per rule 74/89
(never hide a failure, including self-inflicted test-script bugs).

**Two harmless cleanups made while validating:** (1) removed a dead
`... if False else 3.0e5` expression in `enth_mass`'s `initialize=` value
(always evaluated to 3.0e5 anyway; flagged as odd in the prior entry).
(2) Split `vapor_frac` out of `add_properties()` into its own
`define_custom_properties()` call, since `vapor_frac` is not one of
IDAES's `StandardPropertySet` names (closest standard name, `phase_frac`,
is phase-indexed and doesn't fit our single-Mix-phase/lever-rule
convention) -- this silenced a real (if non-fatal) IDAES deprecation
warning. Neither change altered any numerical result; re-ran the
validation script after both, all 4 points still PASS.

**STAGE L IS NOW COMPLETE** -- EOS kernel (24a), P/h/s/g/mu via symbolic
AD (24b), and the smooth-PH-flash StateBlockData (24c, this entry). No
fitted ancillary correlations exist anywhere in the package; saturation/
flash logic is 100% the validated VLE system solved algebraically inside
the IDAES NLP (Option 1, per the 2026-08-18 clarification in Section 22).

---

## 25. Stage M: structural/numerical diagnostics (`DiagnosticsToolbox`)

**STATUS: BUILT AND VALIDATED, 2026-08-18. PASS (0 structural warnings, 0
numerical warnings).** `validate_stage_m_diagnostics.py` runs
`idaes.core.util.diagnostics_tools.diagnostics_toolbox.DiagnosticsToolbox`
(current 2.12.0 location) against the R515B StateBlock at two checkpoints:
square-but-unsolved (state vars fixed, DOF=0, not yet solved) and solved
(a representative two-phase point).

**One REAL package bug found and fixed here (the first genuine Stage-M
finding):** the first run of `report_structural_issues()` raised a hard
"Units problem with expression" error inside the saturation-pressure
residual. Root cause: `r515b_property_package.py`'s 8 internal auxiliary
Vars (T_sat, rho_l_sat, rho_v_sat, T_liq, rho_liq, T_vap, rho_vap) carried
real Pyomo units (K, mol/m^3), but they feed directly into `mixture_ph_
flash_residuals_expr` -- and that function, plus everything under it in
`r515b_pyomo_eos.py` (all of Stage L parts 1-3, validated at machine
precision), was built with ZERO `pyunits` usage anywhere (confirmed by
grep) -- i.e. it is deliberately unitless, raw-SI-value math. Mixing
unit-bearing Vars into that unitless tree is exactly what Pyomo's unit
checker is designed to catch. Fix: made the internal 8-Var/8-Constraint
system fully unitless (dropped `units=` from those 8 Vars, matching how
they were built and validated all along), and added a unit-stripping/
restoring boundary layer exactly where the internal system talks to the
two public unit-bearing quantities it needs: `self.pressure / pyunits.Pa`
and `self._h_molar_jmol / (pyunits.J/pyunits.mol)` strip units
symbolically (not via `value()` -- stays exact) going in;
`res["T_actual"] * pyunits.K` and `(1.0/self._v_mol) * (pyunits.mol/
pyunits.m**3)` restore units on `temperature`/`dens_mol`/`dens_mass`,
the Expressions actually exposed on the StateBlock's Port. This choice
(unitless internals + unit-bearing boundary/public interface) avoids
retrofitting units through ~30KB of already-validated EOS algebra (high
risk, no physical benefit, since those 8 Vars are purely internal) while
keeping the public contract fully unit-safe, matching how
`vapor_compression.py` itself relies on unit-aware arithmetic throughout
(Tmin/Tmax Params/Constraints in K, P_low/P_high Vars in Pa -- confirmed
by direct source inspection). Re-ran `validate_state_block_construction.py`
after the fix: all 4 points (subcooled/two-phase q0.3/two-phase q0.7/
superheated) still PASS at the same tolerances -- no regression.

**One self-caught TEST-DESIGN bug (not a package bug):** a first draft ran
`assert_no_structural_warnings()` on the block with its 3 public state
vars still UNFIXED (DOF=2, per Section 24c's already-documented
convention). This correctly produced "2 Degrees of Freedom"/"Structural
singularity found" WARNINGS -- an accurate description of a deliberately
under-determined block, not a defect. `DiagnosticsToolbox`'s structural-
warning assertion is meant to run on a square (DOF=0) problem, matching
how a StateBlock is actually used downstream. Fixed by fixing the 3
public state vars before running the structural check.

**Final result: 0 structural warnings, 0 numerical warnings** at both
checkpoints. Six non-blocking CAUTIONS (not warnings) were reported at
the solved checkpoint: Jacobian condition number 1.5e10; 2 variables with
extreme values; 8 constraints with mismatched/cancelling terms; a
handful of extreme Jacobian entries/column/row norms. All "Caution" tier
-- `assert_no_numerical_warnings()` passes regardless. These are
consistent with the EOS's inherently wide native value scales (molar
densities ~1e2-1e4 mol/m^3 alongside dimensionless quality alongside
temperatures ~1e2-1e3 K) and are exactly the sort of thing a
`CustomScalerBase` (still-pending Stage L polish item, optional/non-
blocking) would address. Recorded here as a known, non-blocking item per
rule 74/89 -- not hidden, just not yet acted on.

---

## 26. Stage N: `vapor_compression_r515b_integration.py` end-to-end validation

**STATUS: BUILT AND VALIDATED, 2026-08-18. PASS (2/2 test cases).** A new,
self-contained integration copy (`R515BVaporCompressionCycle`) of
`vapor_compression.py`'s cycle flowsheet, substituting `R515BParameterBlock`
for `HelmholtzParameterBlock`, restricted to Mode.PH (the only mode this
package supports). Two Option-1-compliant adaptations were required
because R-515B is a fixed-composition blend, not a CoolProp pure fluid:
`specify_initial_conditions()` uses `solve_pseudopure_saturation_at_t`
instead of `CP.PropsSI`; the diagram-drawing convenience methods build a
real envelope from repeated saturation solves over the already-validated
255-375K/5K grid instead of `HelmholtzParameterBlock`'s compiled-ancillary
`hp_diagram()`/`pt_diagram()`/`ts_diagram()`.

**Three REAL bugs found and fixed (all in `r515b_property_package.py`,
none in the already-validated EOS math) -- this is the first point in the
whole project where an actual Arc-connected flowsheet exercises
`initialize()`/`propagate_state()` and the generic `Compressor`/
`PressureChanger` unit model, so none of these could have been caught by
Stage L/M's own validation, which never exercised that machinery:**

1. `Compressor(...)` construction raised `PropertyNotSupportedError:
   ... entr_mol is not supported`. Root cause: `PressureChanger.
   add_isentropic()` writes its isentropic-assumption Constraint directly
   against `properties_isentropic[t].entr_mol` and `control_volume.
   properties_in[t].entr_mol` -- UNCONDITIONALLY, regardless of
   `amount_basis=MASS`. Our package only exposed `entr_mass`. Fixed by
   adding `entr_mol`, backed by a small additive extension to
   `mixture_ph_flash_residuals_expr` (now also returns `s_liq_jmolk`/
   `s_vap_jmolk`/`S_actual_jmolk`, the same Honeywell-rebased entropy
   `mixture_state_expr` already computed internally but didn't expose --
   zero new EOS math, purely additive, no change to the 8 existing
   residuals). Validated separately in `validate_entropy_and_tsat_vs_
   core.py` against `mix_entropy_direct`: 4/4 PASS, rel_err 1e-10 to
   8.6e-9 for `entr_mass`, 1.4e-13 for the new `temperature_sat` alias
   added at the same time (see below).
2. Compressor `initialize()` then failed identically for `enth_mol` --
   `PressureChanger.init_isentropic()` also accesses it directly and
   unconditionally. Fixed trivially: `enth_mol` is just the already-
   existing internal `_h_molar_jmol` Expression, exposed under the
   standard property name.
3. `initialize()`'s call to `propagate_state()` (used by every unit
   model's own `initialize()` to copy values across an Arc) raised
   `TypeError: ... cannot set a <class 'pyomo.core.base.expression.
   Expression'>`. Root cause: an earlier point in this session had
   overridden `define_port_members()` to add `temperature` (an
   Expression) to the Port, reasoning from `StateBlock.build_port()`'s
   `Reference()` mechanism accepting Expression port members at BUILD
   time -- true, but irrelevant to `propagate_state()`, which needs every
   port member to be independently settable (a Var). Direct inspection of
   the REFERENCE `HelmholtzStateBlockData.define_port_members()` (the
   contract this package must match) confirms it is NOT overridden there
   either -- it uses the StateBlockData base class default (`define_state_
   vars()`, i.e. just the 3 public state vars). Fixed by removing the
   override entirely. `temperature`/`vapor_frac`/`entr_mass`/`entr_mol`/
   `dens_mass`/`dens_mol`/`temperature_sat` remain fully reachable exactly
   the way `vapor_compression.py` already reaches them for the reference
   package: directly off `control_volume.properties_out[0]`/
   `properties_in[0]`, never through the Port in PH mode (confirmed
   neither `vapor_compression.py` nor the new integration file ever
   references `.outlet.temperature[0]`-style Port access in PH mode).

**One independent, deliberate deviation from `vapor_compression.py`
(fixed ONLY in the new file -- the original remains untouched, per the
master task's hard constraint):** its temperature-bound checks use bare
`if bound:`, and Python's `0` is falsy -- so a caller who wants a
legitimate bound of exactly 0 degC gets it silently skipped, with no
error. Caught empirically: a first end-to-end test with `evaporator_
temperature=(-15, 0)` ran away to an unphysical COP~34.5 (compressor
inlet superheated to 412 K, since the intended upper bound was never
actually active). Fixed as `if bound is not None:` throughout `set_
specifications()` in the NEW file only. Verified both ways: a bound
avoiding literal 0 gives a physically sensible COP=3.93 (evaporator
outlet -1 degC, condenser outlet 35 degC + subcooling, full state
consistency confirmed unit-by-unit); an exactly-0-degC upper bound is now
correctly honored (evaporator outlet pinned at exactly 273.150 K,
COP=4.06).

**Result (`validate_stage_n_integration.py`): 2/2 PASS.** Both cases:
solver reports optimal termination; COP within a physically plausible
1-10 range (3.93 and 4.06); closed-loop pressure/enthalpy consistency
(expansion-valve outlet == evaporator inlet, exact to solver tolerance --
`p_err=0 Pa`, `h_err=0` relative -- since that Arc's `flow_mass_equality`
is deliberately deactivated for the closed loop, matching `vapor_
compression.py`'s own "closed, circular loop" comment) both hold; the
0-degC-bound regression check passes.

**STAGE N IS NOW COMPLETE.**

---

## 27. Stage O: final full-suite regression

**STATUS: COMPLETE, 2026-08-18. ALL 15 numbered validation scripts PASS
in a single end-to-end re-run** (not just re-run individually earlier in
each stage's own entry -- re-executed together, in sequence, as a final
closing check that nothing regressed across the whole build):

`validate_table1_vs_oracle.py`, `validate_core_vs_oracle.py` (Stage J
OVERALL: PASS), `validate_fugacity_vs_oracle.py`,
`validate_quality_isotherm_vs_oracle.py`, `validate_isentrope_vs_oracle.py`,
`validate_pseudopure_saturation.py`, `validate_pyomo_eos_vs_core.py`,
`validate_pyomo_state_vs_core.py`, `validate_pyomo_saturation_vs_core.py`,
`validate_pyomo_flash_vs_core.py`, `validate_state_block_construction.py`,
`validate_entropy_and_tsat_vs_core.py`, `validate_stage_m_diagnostics.py`,
`validate_stage_n_integration.py` -- all report `PASS`/`OVERALL: PASS` at
their respective tolerances (machine epsilon through 1e-6, as documented
in each script's own section above). `validate_ancillary_guess.py` is not
a pass/fail gate (it's initialization-support machinery, Section 23) --
it re-ran cleanly and re-saved its results JSON with its already-
documented, already-acceptable liquid/vapor-branch accuracy profile
unchanged.

No new bugs surfaced in this final pass -- every fix made during Stages
L/M/N (departure-function transcription, complementarity pairing, seed-
basin non-uniqueness, DOF-expectation test bug, units-architecture bug,
structural-diagnostics test-design bug, entr_mol/enth_mol support,
propagate_state/define_port_members incompatibility, the independent
falsy-zero fix) held under this combined re-run with no interaction
effects between them.

**FINAL PROJECT STATUS: Stage A-O of A-O COMPLETE.**

- Stage A-K (non-IDAES thermodynamic core): exact (rel_err=0.0) match vs.
  the oracle across direct states, bubble/dew VLE, critical point,
  quality lines, isotherms, isentropes.
- Stage L (native-Pyomo IDAES property package): EOS kernel (part 1),
  P/h/s/g/mu via genuine Pyomo symbolic differentiation (part 2), and the
  smooth 8-equation PH-flash StateBlockData with a validated x1=y1=z1
  near-azeotropic simplification (part 3) -- all independently validated
  against the SciPy reference core, rel_err 0 to ~1e-7 depending on
  whether via direct value() evaluation or a real IPOPT solve.
- Stage M (IDAES structural/numerical diagnostics): 0 structural
  warnings, 0 numerical warnings via `DiagnosticsToolbox`.
- Stage N (`vapor_compression.py` integration): a new, self-contained
  integration copy runs the full evaporator/compressor/condenser/
  expansion-valve cycle end to end, with COP optimization converging to
  physically plausible answers (COP 3.93-4.06 across tested cases).
- Stage O (this entry): full-suite regression, no interactions/
  regressions found.

**No fitted R-515B ancillary correlations exist anywhere in the delivered
package.** Every saturation/flash/entropy quantity is either an exact
shared code path with the oracle (Stages A-K) or the validated binary-
Helmholtz VLE system solved algebraically inside the IDAES NLP (Stage L,
Option 1 per the 2026-08-18 standing commitment, Section 22) -- confirmed
by direct grep and import-graph inspection at both the Stage-M and Stage-N
checkpoints. `mixture_fully_validated.py` (the oracle) and
`vapor_compression.py` remain read-only and unmodified throughout; only
`linear_model_codex.py` was imported (never modified) as the pre-approved
shared production dependency; no Git command was invoked at any point;
all breadcrumbs were updated append-only, with every self-caught bug
(there were many, across every stage) documented in full, never hidden,
per rule 74/89.

---

## 24. Final summary

**STATUS: IN PROGRESS -- Stage A-K COMPLETE (D/E by construction; F/G/H/I/
J/K all by the new independent module `r515b_helmholtz_core.py`,
numerically validated EXACT (rel_err=0.0) vs. the oracle across direct
states, bubble/dew, critical point, quality lines, isotherms, AND
isentropes -- the full non-IDAES thermodynamic core is done). Stage L is
now FULLY COMPLETE: part 1 (native-Pyomo EOS kernel, Section 24a, rel_err
0 to 2.2e-16 vs. the NumPy core at 3 states), part 2 (native-Pyomo P/h/s/
g/Z and mu1/mu2 via Pyomo's own symbolic differentiation, Section 24b,
rel_err 0 to 1.9e-15 across 18 checks), and part 3 (the actual
PhysicalParameterBlock/StateBlockData classes -- metadata, units,
construction, degrees of freedom, initialize()+solve -- Section 24c,
4/4 PASS, T rel_err ~1e-9, vapor_frac abs_err ~1e-8, dens_mass rel_err
~1e-7). Stage M (IDAES structural diagnostics), N (vapor_compression.py
integration copy), O (final regression) are the remaining work.**

```
THERMODYNAMIC VALIDATION
P: PASS -- 3 direct states rel_err 0 to 1.4e-16 (Section 5/6); bubble/dew
   P(T) at 6 temperatures rel_err 0.0 (Section 8/9); critical-point P
   rel_err 0.0 (Section 10, real-data-informed case). Quality/isotherm/
   isentrope states NOT YET STARTED (Stage K).
h: PASS -- same coverage as P (direct states, bubble/dew h_l/h_v, critical
   h). Rest NOT YET STARTED.
T: PASS for critical point (T_K solved-for output, rel_err 0.0, Section 10).
   Not yet independently exercised as a solved-for output elsewhere
   (isotherms/isentropes solve at fixed T by construction, so T isn't an
   "output" there in the same sense) -- Stage K pending.
s: PASS (partial) -- 3 direct states rel_err 0 to 1.4e-16 (raw, pre-
   Honeywell-rebase convention; rebase offset location confirmed, Section
   4/21). Not yet exercised on the bubble/dew/critical/quality/isotherm/
   isentrope states (Stage K pending) even though the underlying
   `mix_entropy_direct` function used everywhere is the same validated one.

Reducing functions: PASS (by construction, Section 3) -- shared code path
  with the oracle via linear_model_codex.bell2023_Tred_vred.
Helmholtz terms: PASS (Section 4) -- alpha_mix/alpha_tau_mix/ar_del_mix
  exact match (0.0 abs diff) via mixture_alpha0_alphar_derivs (shared code
  path), now embedded in the validated `r515b_helmholtz_core.py` module
  used throughout Sections 5-10.
Fugacity/chemical potentials: PASS (Section 7) -- mu1/mu2 exact match
  (0.0 rel err) at direct states AND as actually used inside the converged
  bubble/dew VLE solves (Section 8/9, since mu1_L=mu1_V/mu2_L=mu2_V are the
  residual equations solve_bubble_at_t/solve_dew_at_t drive to zero).
Bubble line: PASS (Section 8) -- 6 temperatures, rel_err 0.0
Dew line: PASS (Section 9) -- 6 temperatures, rel_err 0.0
Critical point: PASS (Section 10) -- rel_err 0.0 on T/rho/P/h (real-data-
  informed case); crude-guess non-convergence fallback ALSO matches exactly
  (same failure mode reproduced, not just "both failed")
Quality states: PASS (Section 11) -- 90 rows, rel_err 0.0
Isotherms: PASS (Section 12) -- 4 families (two-phase/liquid/vapor/
  supercritical), 46 T-values total, rel_err 0.0
Isentropes: PASS (Section 13) -- 15 isentropes, 136+360+456 points across
  two-phase/liquid-side/vapor-side, rel_err 0.0, INCLUDING the near-critical
  entropy-ramp/dome-reentry-guard/adaptive-walk machinery
Derivatives (analytic vs. finite-difference): implicitly exercised (the
  critical-point solve's own FD derivatives, `_dp_drho_fd`/
  `_pressure_rho_derivatives_fd`, are ported and validated as part of
  Section 10) but not separately tabulated against an independent FD check
  at additional states -- Stage G/H residual item, low priority given
  Section 10's exact match.

INITIALIZATION-SUPPORT MACHINERY (not on the Stage A-O critical path)
Ancillary composite initial-guess (Section 23): BUILT AND VALIDATED --
  liquid branch good (~2% typical rel. error vs. oracle's real bubble
  density), vapor branch poor at low T (up to 84% rel. error), improving
  near Tc (~5%); recommended for liquid-seed use only, documented limitation
  not a blocker.

IDAES STRUCTURAL VALIDATION
StateBlock construction: PASS (Section 24c) -- R515BParameterBlock +
  R515BStateBlock construct cleanly inside a real FlowsheetBlock
Units: PASS -- public state vars/Expressions (flow_mass, pressure,
  enth_mass, temperature, dens_mass, dens_mol) carry proper Pyomo units;
  8 internal auxiliary Vars deliberately unitless (matching the validated
  unitless EOS kernel) with a symbolic unit-stripping/restoring boundary
  layer at the interface (Section 25 -- real bug found+fixed via
  DiagnosticsToolbox's structural check, not caught by Stage L's own
  value()-based validation since that never exercises Pyomo's unit
  checker). Metadata split correctly between add_properties() (7 standard
  names) and define_custom_properties() (vapor_frac).
Degrees of freedom: PASS (Section 24c) -- 2 before fixing state vars, 0
  after (matches the SAME accounting convention as the reference
  HelmholtzParameterBlock, verified directly)
Initialization: PASS (Section 24c) -- SciPy-seeded (flash_ph_pseudopure,
  rule-44 carve-out) .initialize() + IPOPT solve reproduces the reference
  at all 4 test points (subcooled/two-phase q0.3/two-phase q0.7/superheated)
Structural/numerical diagnostics (DiagnosticsToolbox): PASS (Section 25)
  -- 0 structural warnings, 0 numerical warnings at both a square-unsolved
  and a solved (two-phase) checkpoint; 6 non-blocking conditioning
  CAUTIONS noted (Jacobian condition number 1.5e10 etc.), tied to the
  still-pending optional scaling polish item
Scaling: NOT STARTED (Stage L polish item, optional CustomScalerBase --
  does not block Stage N/O; motivated by Section 25's conditioning
  cautions)

VAPOR_COMPRESSION INTEGRATION
Construction: PASS (Section 26) -- R515BVaporCompressionCycle (new file,
  vapor_compression.py untouched) constructs the full evaporator/
  compressor/condenser/expansion_valve Arc-connected flowsheet
Initialization: PASS (Section 26) -- .initialize() + propagate_state
  across all 4 Arcs completes without error (after fixing entr_mol/
  enth_mol support and the define_port_members()/propagate_state
  incompatibility)
Solve: PASS (Section 26) -- optimize_COP() converges to optimal
  termination in both validated test cases
Physical states: PASS (Section 26) -- COP 3.93/4.06 (physically
  plausible refrigeration-cycle range); closed-loop pressure/enthalpy
  consistency exact to solver tolerance; 0-degC-bound regression check
  passes

REPOSITORY INTEGRITY
Existing source files unchanged: PASS (verified 2026-08-17: only append-only
  edits made to PROJECT_CONTEXT.md; mixture_fully_validated.py read-only;
  vapor_compression.py read-only; linear_model_codex.py read-only (imported,
  not modified); no other existing file touched)
Breadcrumb append-only requirement respected: PASS
Git prohibition respected: PASS (no Git command invoked at any point)

OVERALL: **PASS -- Stage A-O of A-O COMPLETE (2026-08-18).** The FULL
non-IDAES thermodynamic core (Stage A-K, exact rel_err=0.0 vs. the
oracle), the FULL native-Pyomo/IDAES property package (Stage L parts
1-3), structural/numerical diagnostics (Stage M, 0 warnings), a working
end-to-end `vapor_compression.py` integration copy (Stage N, COP
optimization converges to physically plausible answers), and a final
full-suite regression (Stage O, all 15 validation scripts PASS together,
Section 27) are all built and independently validated. The user
authorized continuing autonomously through anomaly-free work while away
("Continue with everything per the prompt... build the entire model...
will be away"; "you will only pause if you encounter pause conditions
laid out in the prompt"; "Continue from where you left off" x2), and
separately reaffirmed mid-Stage-L that R-515B is not a zeotrope and that
no fitted ancillary correlations may be introduced without pausing for
approval (Section 22, 2026-08-18 entry) -- verified at every subsequent
checkpoint (Stage M, N, O) that the delivered package complies: no
ancillary-correlation drift anywhere; every saturation/flash/entropy
quantity is either an exact shared code path with the oracle or the
validated binary-Helmholtz VLE system solved algebraically inside the
IDAES NLP. No anomaly encountered across Stages D-O requiring a genuine
pause -- every self-caught bug across every stage (input-precision/unit-
conversion in D-K; departure-function transcription, complementarity
pairing, seed-basin non-uniqueness, and a DOF-expectation test-script bug
in Stage L; a real units-architecture bug and a test-design bug in Stage
M; missing entr_mol/enth_mol support and a propagate_state/define_port_
members incompatibility in Stage N) was fixed and documented, never
hidden (Sections 4/13/21/22/24c/25/26/27). `mixture_fully_validated.py`
and `vapor_compression.py` remain read-only and unmodified throughout; no
Git command was ever invoked; all breadcrumbs were updated append-only.
**No further action pending under the master task as currently scoped.**
```
