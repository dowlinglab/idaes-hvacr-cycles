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

**Known issue, NOT yet fixed:** the rendered isentropes visually mismatch the Honeywell chart (lines that should fan out across the liquid region instead crush together; several terminate before reaching the chart's pressure limits). Two candidate causes, neither resolved: (1) the unbounded `scipy.optimize.root` solve may be converging to the wrong phase branch when the seed-to-target entropy gap is large -- candidate fix is switching to bounds-constrained `scipy.optimize.least_squares`, matching the bubble/dew solve's existing sigmoid-bounded pattern; (2) a possible entropy reference-state offset between this EOS and Honeywell's own convention ("h=200kJ/kg, s=1.00kJ/kg-K, sat. liq. at 0C" footnote) -- Honeywell's real chart shows NONE of these 14 isentropes crossing into the two-phase region at all, which our current output contradicts, suggesting the absolute entropy scale (not just the tracing logic) may need reconciling.
