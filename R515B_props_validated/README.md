# R515B_props_validated

Snapshot of the R-515B (R-1234ze(E)/R-227ea) mixture-model code, taken 2026-08-12 as the reference version for the paper's SI. This folder freezes the code state discussed in that debugging session; the live/working copies of these files remain at the repo root and may continue to evolve.

## Files

### `pressure_validated_model.py`
Independent implementation of the Bell (2023) binary Helmholtz mixture model for R-1234ze(E)/R-227ea (the R-515B surrogate pair), including chart-basis enthalpy alignment (`ChartReference`, `compute_props_with_chart_h`) and a density-bracketing root solver (`solve_rho_mass_for_P`) for pressure-to-density inversion. Three bugs were found and fixed in this file on 2026-08-12:

1. `beta_v` transcription error: `0.99290` corrected to `0.999290` (Bell 2023 Table 2).
2. Missing `/8.0` divisor in the `nu12` reducing-volume mixing rule (paper Eq. 5).
3. `compute_alpha_res_pure` (phi_residual_type 2/3/4 branches) was reusing the first term's exponent (`t1`) for every term in the sum instead of each term's own per-index exponent. This was the dominant bug (~92x effect on its own).

After all three fixes, this file's output matches `linear_model_codex.py` and `mixture_model_one_point.py` to floating-point/finite-difference precision at the test state T=298.15 K, rho=1179.8 kg/m3, w1=0.911 (p = 4936.44 kPa).

### `mixture_model_one_point.py`
The primary hand-maintained R-515B property-query implementation (author-edited). Contains the same Bell (2023) binary Helmholtz math (reducing functions, departure function, IDAES-Helmholtz-JSON pure-fluid evaluators) as a single-state p-h query interface. The same `beta_v` transcription fix (`0.99290` -> `0.999290`) was applied independently in this file on the same date.

## Known limitation carried by both files

Both implementations are validated to ~1.0-1.1% MAPE against the Honeywell Solstice N15 (R-515B) chart bubble/dew pressure data *when the model's own VLE solver determines the equilibrium liquid density self-consistently*. They are **not** validated for a liquid density supplied from an independent/measured source: at R-515B's real composition (x1 ~ 0.9385 mole fraction R-1234ze(E)), which falls outside the composition range Bell (2023) actually fit (x1 = 0.33-0.68, their Table 12), the model's own equilibrium liquid density comes out ~1.95% low relative to the real Honeywell value. Because of an extremely steep local dP/drho sensitivity on the liquid branch (a near-total ideal/residual cancellation in the compressibility factor), that ~2% density gap amplifies into a ~10x pressure error if a real/measured density is ever fed into the pressure model directly instead of being solved for self-consistently.

Practical implication: pressure/enthalpy/entropy outputs from a fully self-consistent solve (VLE bubble/dew, cycle COP, p-H / s-T diagrams, etc. -- where P, T, H, S, and rho are all computed jointly from the model, never mixing in an externally-measured density) stay reliable at the ~1-2% level. Liquid density itself, or any pressure computed from an externally supplied density, should not be trusted at R-515B's real composition without further correction.

### Independent confirmation (NIST)

This exact density-extrapolation gap is independently confirmed on NIST's own REFPROP issue tracker: [`usnistgov/REFPROP-issues#750`, "R-515B Density Uncertainty"](https://github.com/usnistgov/REFPROP-issues/issues/750). Another user ran the same composition-range argument (R-515B's real composition is ~0.939 mole fraction R-1234ze(E); Bell (2023)'s validation only covers 0.33-0.67) and asked NIST directly whether the paper's density-accuracy claims hold outside that range. A NIST REFPROP maintainer (`marciahuber`) responded by pointing to an updated `HMX.BNC` file and a dedicated `R515B.MIX` file with improved interaction parameters, available on request from REFPROP@NIST.GOV but not yet in the public release. Nobody in that thread confirmed the updated file's actual accuracy at the real R-515B composition, so it is the best available lead, not a confirmed fix.

## 2026-08-13 additions: Honeywell p-h chart reproduction (saturation dome, quality lines, isotherms, isentropes)

Feature-by-feature reproduction of the Honeywell Solstice N15 (R-515B) Technical Data Sheet's p-h chart (page 2), built as a linear sequence of file copies, each adding one chart element on top of the last. All isotherm/isentrope/quality-line work is presentation-layer only, added in separate validation-only functions -- the underlying SI bubble/dew solve (`run_true_vle_envelope`, `mix_state`) is untouched by any of it.

### `mixture_true_vle_copy.py`
The original ancestor file all of today's copies descend from. Only change today: docstring/header metadata (credited Claude AI alongside Codex, added an edit date, and a short note framing this file as the "vanilla" baseline that later files validate against).

### `mixture_dome_validation.py`
Independently-validated saturation-dome (bubble/dew) implementation, confirmed to 1.84% MAPE against the Honeywell P-T table by a separate cross-check.

### `mixture_dome_validation_pseudo_pure.py`
The file where the bubble-branch loss function was switched `cauchy` -> `huber` (fixing a non-monotonic enthalpy "kink"/S-turn), and where the true mixture critical-point solver was developed (bisection on spinodal-dip existence, converges to within 0.04% of Honeywell's Tc and 0.56% of Pc). Superseded by `Honeywell_T_P_revalidation.py` (see below), which copied this file's fixed/validated state forward as the base for all subsequent chart-reproduction work. Kept here for history; not further edited after the copy.

### `diagnose_spinodal.py`
Diagnostic-only script (not part of the production CLI) supporting the critical-point solver work above -- prints spinodal density pairs at trial temperatures to visualize the bisection search.

### `Honeywell_T_P_revalidation.py`
Copy of the fixed `mixture_dome_validation_pseudo_pure.py`. Adds: the 86-point Honeywell TDS P-T table (0-170F, page 3) revalidation (86/86 converged, 1.84% MAPE, matches the independent `mixture_dome_validation.py` cross-check), and the IP-unit (psia, Btu/lbm) conversion architecture (`to_honeywell_units`, `plot_envelope_honeywell_units`) -- conversion applied only at the final plotting step, never touching the SI core solve.

### `mixture_quality_line_validation.py`
Copy of `Honeywell_T_P_revalidation.py`. Adds quality lines (x=0.1-0.9, lever-rule interior two-phase states: h(x)=h_l+x*(h_v-h_l), pure post-processing of already-solved bubble/dew rows). Validated against the reference chart; all quality lines converge to the exact solved critical point at the dome apex.

### `mixture_isotherm_validation.py`
Copy of `mixture_quality_line_validation.py`. Adds the full isotherm family, all in red matching Honeywell's own color: two-phase segments (-20F to 220F, horizontal since T and P are linked on the saturation curve), subcooled-liquid extensions (up to the chart's 1350 psia ceiling), superheated-vapor extensions (down to the chart's 15 psia floor), and supercritical isotherms (240F-400F, entirely above the solved Tc=228.0F, no bubble/dew anchor). Also adds chart polish matching the Honeywell convention exactly: box border, axis label text/tick values (psia, Btu/lbm, non-uniform gridlines read off the reference chart), and `annotation_clip=False` on all labels (matplotlib silently drops an annotation whose anchor falls outside the current axis view, which is what caused a quality-line-label regression mid-session).

### `mixture_isentrope_validation.py`
Copy of `mixture_isotherm_validation.py`. Adds isentropes (constant specific entropy, blue matching Honeywell's own color, the 14 values printed on the TDS chart: 0.22-0.49 Btu/(lbm-R)) -- a two-phase segment (algebraic lever rule on entropy, no new EOS solve) plus subcooled-liquid/superheated-vapor single-phase segments (2-equation/2-unknown `scipy.optimize.root` solve for (T,rho) at fixed P and target entropy, since unlike an isotherm neither T nor rho is fixed here). Entropy itself is computed directly from the Helmholtz mixture model (`s/R = tau*alpha_tau - alpha`), matching the established direct-entropy pattern in `mixture_model_one_point.py`'s `compute_table1_properties`, not derived via `s=(h-g)/T`.

**Known issue as of 2026-08-13, since resolved/superseded -- see the 2026-08-14 section below:** the rendered isentropes originally visually mismatched the Honeywell chart. Root cause was eventually isolated to (a) a genuine vapor-side solver bug (anchor collision, fixed) and (b) a real entropy reference-state offset that grows toward the critical point (confirmed, but the correction-curve fix attempt for it was ultimately abandoned in favor of a structural explanation -- see below). Neither of the two candidate causes originally guessed here (unbounded-solver wrong-branch convergence; a flat reference-state offset) turned out to be the full picture.

## 2026-08-14 additions: isentrope solver fixes, entropy investigation, and a consolidated validation report

Continuation of the chart-reproduction work above. Two separate threads: (1) fixing real solver bugs in the vapor-side isentrope extension, and (2) a full investigation into why entropy specifically -- not pressure or enthalpy -- was the property most out of line with Honeywell's chart, ending in a structural explanation rather than a numerical correction.

### `mixture_isentrope_validation.py` (updated, this is the active/current file)
Two real bugs fixed in `compute_isentrope_vapor_side()`:

1. **Anchor-collision bug (isentropes 0.41-0.49 Btu/lb-R):** these five isentropes were all being seeded from the same single dew-branch point (the branch's own interior entropy peak), then asked to jump a large, increasingly extreme entropy gap (5.5 up to 44.8 J/(mol*K)) in one `scipy.optimize.root` call -- producing implausible shapes. Fixed with an incremental entropy ramp (small isobaric substeps, each capped and warm-started from the previous one) plus an upward pressure walk (these isentropes never re-enter the two-phase dome at any pressure, so walking up to the 1350 psia ceiling is physically correct, not just cosmetic), gated to only fire when the target entropy exceeds the dew branch's own maximum -- otherwise it can walk a line straight through the dome (a real regression that was caught, diagnosed, and fixed mid-session).
2. **0.39 Btu/lb-R partial-to-full fix:** this isentrope's target entropy sits just below the dew branch's cold-end value, so extending it downward/outward from its anchor always re-enters genuinely two-phase territory (confirmed three separate ways: adaptive step-halving, a temperature-stepped walk, and ideal-gas-biased seeding all independently hit the same wall). The fix was to walk the OTHER direction instead -- upward through the critical point, where no two-phase region exists at any entropy -- extending it from a 7-point stub to a full 47-point curve (477.6-1350 psia).
3. **0.37 Btu/lb-R -- reported as a known model limitation, not a bug.** Its liquid-side anchor (via the critical-point fallback) is numerically accurate to ~0.02% against an independently-interpolated fine grid, but the true near-critical point's own enthalpy is still ~2.9% below the model's solved critical enthalpy -- consistent with genuine steep near-critical curvature (dh/ds diverges approaching Tc for any real EOS), not a computational error.

All isentropes are validated against a real bubble/dew-line dome-membership check (not just point counts) and confirmed to never cross the two-phase interior anywhere on the current chart.

### The entropy investigation: correction curve attempted, then abandoned for a structural finding
A T-dependent entropy correction curve was built (`fit_entropy_correction.py`, 15 hand-digitized (s,P,H) points from the Honeywell chart, fit to a quadratic residual vs. temperature), wired into a working copy, and initially verified matching the fit exactly. It was then found to break solver convergence for three other isentropes via a step-discontinuity at the liquid/vapor branch switch temperature, and after evaluating an "apply after convergence" alternative (mathematically identical, doesn't remove the discontinuity), the correction-curve approach was abandoned entirely. **Both isentrope files currently use only the original flat reference-state offset (`ENTROPY_REFERENCE_OFFSET_JMOLK = -1.4595`), no T-dependent correction.**

In its place, a structural, code-grounded explanation was found for why entropy specifically is affected: of the three governing identities, `P = rho*R*T*(1 + delta*alphar_del)` and `H/(RT) = 1 + tau*alpha_tau + delta*alphar_del` are both built only from *derivatives* of the departure function (`alphar_del`, `alpha_tau`), while `S/R = tau*alpha_tau - alpha` is the only one built from the raw, undifferentiated `alpha` value. A constant/level-shift bias in the composition-extrapolated departure function (Bell 2023's fit range is x1=0.33-0.68; real R-515B is x1~0.9385) contributes zero to P and H's own derivatives-based formulas but passes directly into S with coefficient exactly -R. This is consistent with everything observed across this whole validation effort: P and H validate well almost everywhere, S is the one property persistently off, worst near the critical point.

### `mixture_isentrope_correction_validation.py` -- abandoned experiment, kept for history
The working copy where the T-dependent entropy correction above was built, wired in, found to regress solver convergence, and reverted. No longer kept in sync with `mixture_isentrope_validation.py` (diverged after the revert as further vapor-side fixes were made only to the active file). Not part of the current validated result.

### `mixture_fully_validated.py`
A frozen snapshot of `mixture_isentrope_validation.py`'s current state, byte-identical as of this entry, designated as the reference "fully validated vanilla" copy -- i.e. the version whose bubble/dew points are validated against Honeywell data and from which the rest of the codebase's correctness is assumed to follow.

### `mixture_isotherm_validation_correction.py`
A copy of `mixture_isotherm_validation.py`, created as a potential home for a density/pressure correction layer (see `honeywell_saturation_correction.py` below). Not currently wired to anything; status unresolved.

### `honeywell_saturation_correction.py`
A standalone density/pressure correction layer fit against digitized Honeywell saturation data (PCHIP, shape-preserving, fixed an earlier overshoot problem). Built and verified working, but superseded before being wired into anything: the isentrope investigation above concluded the real gap was an entropy-formula bias, not a density/pressure one, so this file's approach -- while functional -- isn't the fix for the isentrope mismatch. Kept for potential future use (e.g. if a density/pressure correction is separately needed for the isotherm liquid extension).

### `linear_model_codex.py`
Shared EOS module (Bell 2023 reducing functions, departure function, pure-fluid Helmholtz evaluators) imported by the other files in this folder. Infrastructure dependency, not a standalone deliverable.

### Diagnostic scripts (read-only, not part of the production CLI)
`check_entropy_reference.py`, `check_dew_entropy_range.py`, `diagnose_isentrope_026.py`, `diagnose_isentrope_vapor_037.py`, `diagnose_vapor_isentrope_anchor_collision.py`, `diagnose_037_critical_approach.py`, `diagnose_039_upward_from_anchor.py`, `plot_bubble_entropy_shape.py`, `plot_dew_entropy_shape.py`, and `fit_entropy_correction.py` (the abandoned correction's fitting script, see above). Each supports one specific finding described above; none modify the production files.

### `verification/`
Supporting plots and metadata (`r515b_true_vle_envelope*.png`, `current_model_honeywell_units_FIXED*.png`, `r515b_true_vle_metadata.json`) generated during this session's diagnostic work.

## Current overall status (2026-08-14)

The model is being treated as validated, with every known error and its magnitude recorded in `VALIDATION_REPORT.md` (this folder) -- that file is the authoritative, audit-ready summary of what this model gets right and wrong; this README describes what each file is and does. In short: pressure and enthalpy validate well almost everywhere (the ~1-2% density-extrapolation caveat above still applies to externally-supplied densities). Entropy is the weakest property, with a bias that grows approaching the critical point, for the structural reason described above -- not yet corrected, but now well understood and documented rather than papered over. The two near-critical isentropes (0.37, 0.39 Btu/lb-R) each have their own documented, non-bug explanation. None of this reaches the repository's cycle-simulation models today (they use IDAES's native pure-component property package, not this mixture model) -- `VALIDATION_REPORT.md` documents exactly how that would change if this mixture model is ever wired into a cycle's compressor block, since compressor duty is computed from an isentropic (constant-entropy) constraint.
