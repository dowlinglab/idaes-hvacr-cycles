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
