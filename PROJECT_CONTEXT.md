# Project Context

## Project Name
idaes-hvacr-cycles

## Last Updated
2026-03-03

## Current Objective
Establish and maintain a persistent project research memory backbone for disciplined development and reproducible engineering decisions.
Provide a stable query interface for pressure/enthalpy at a single state point while preserving strict model behavior.
Implement a pure-fluid Helmholtz-based saturation solver and P-H dome pipeline for phase-aware validation workflows.

## System Architecture Overview
- Pending detailed capture.
- This document will track separation of thermodynamic property calls, numerical grid generation, plotting, and validation routines as implementation evolves.
- Added pure-fluid saturation module (`helmholtz_saturation.py`) and wrapper CLI (`scripts/plot_ph_dome.py`) for dome CSV/figure generation.

## Key Technical Decisions
- 2026-03-02: Adopted `PROJECT_CONTEXT.md` as mandatory session-start/session-end context memory at repository root.
- 2026-03-02: Standardized required section schema for traceability across sessions.
- 2026-03-02: Defined "point" for p,h query API as `(T [K], rho_mass [kg/m^3], composition x1 [mol/mol], x2=1-x1)`.
- 2026-03-02: Standardized caller return units as `p [kPa]` and `h [kJ/kg]`.
- 2026-03-02: Chose strict-failure behavior for missing derivative mapping (no non-strict fallback added without explicit approval).
- 2026-03-02: Implemented strict fixed-composition chain-rule mapping in `mixture_alpha0_alphar_derivs` using existing `alpha0_idaes_with_derivs` and `alphar_idaes_with_derivs` outputs.
- 2026-03-02: Rejected a derivative-free or approximate fallback pathway; retained strict derivative-based formulation.
- 2026-03-02: Caller-layer composition input is mass fraction (`x1` in `query_ph_codex.py`), converted to mole fraction before calling `compute_pressure_enthalpy`.
- 2026-03-02: Renamed caller interface field to `w1` (mass fraction) in API and CLI for explicit basis clarity.
- 2026-03-02: Updated `linear_model_codex.py` CLI to also accept `--w1` (mass fraction) and convert internally to mole fraction `x1` before solving.
- 2026-03-02: Added a dedicated two-phase dome plotting utility that consumes externally provided saturation data (`P_bar`, `h_l_kJkg`, `h_v_kJkg`) instead of embedding a new flash solver.
- 2026-03-02: Implemented pure-fluid saturation equations with residuals `r1=P_l-P_v` and `r2=g_l-g_v` using Helmholtz identities and analytic derivative mappings from reduced variables.
- 2026-03-02: Chose Newton with damping as primary saturation solver, with SciPy root fallback and auxiliary-density fallback for robustness when strict solve stalls.
- 2026-03-02: Enforced strict saturation acceptance gates (no silent fallback acceptance):
  `r_P = |P_l-P_v|/max(1,0.5*(P_l+P_v)) <= 1e-6`,
  `r_mu = |g_l-g_v|/(R_u*T) <= 1e-6`,
  and `rho_l > rho_v*(1+1e-8)`.
- 2026-03-02: Auxiliary/surrogate density estimates are now recovery seeds and failure diagnostics only; they are never treated as converged saturation points.
- 2026-03-02: Saturation run output schema switched to per-temperature status rows with explicit `status ∈ {CONVERGED, FAILED_CHECKS, DIVERGED}` plus residual metrics.
- 2026-03-02: Added a separate fixed-composition pseudo-dome path (`mixture_pseudo_dome.py`) for exploratory R515B approximation without altering strict pure-fluid saturation logic.
- 2026-03-02: Per user request, executed pseudo-dome run from a copied file (`mixture_pseudo_dome_copy.py`) to isolate experimentation from the original implementation file.
- 2026-03-02: Added a separate true-VLE experimental solver module (`mixture_true_vle_copy.py`) using equations
  `P_l=P_v`, `mu1_l=mu1_v`, `mu2_l=mu2_v` with finite-difference `mu_i` from `A(T,V,n1,n2)`.
- 2026-03-02: Bubble/dew formulation choice for fixed overall composition `z`: bubble branch solved at `x_liq=z`, dew branch solved at `y_vap=z`.
- 2026-03-02: Changed core p,h API in `linear_model_codex.py` to accept `w1` mass fraction directly; mole fraction conversion now occurs internally (`x1_from_w1`) to eliminate interface ambiguity.
- 2026-03-02: Updated core enthalpy assembly in `compute_pressure_enthalpy` to explicit Lemmon-style mixture form:
  `h0/(RT) = 1 + sum_i x_i*tau_i*(d alpha0_i / d tau_i)`,
  `h/(RT) = h0/(RT) + tau*(d alpha_mix^r/d tau) + delta*(d alpha_mix^r/d delta)`.
- 2026-03-02: Implemented a new Table-1 property evaluator in `linear_model_codex.py` (`compute_table1_properties`) including
  `Z, p, h, u, s, Cv, Cp, w, f_i, phi_i, mu_i` for the binary mixture state.
- 2026-03-02: Explicit approximation choice for Table-1 closure:
  second residual derivatives (`a_tt^r, a_dd^r, a_dt^r`) and composition derivative terms
  (`∂(n a^r)/∂n_i` for fugacity) are currently computed by finite differences.
- 2026-03-02: Corrected Helmholtz residual expression mapping to IDAES form by using term-specific exponents `t[i]` in
  `alphar_idaes_with_derivs` for `phi_residual_type` 2/3/4 (previous code incorrectly reused `t[1]` across many terms).
- 2026-03-02: Replaced finite-difference composition derivatives for fugacity in `compute_table1_properties` with analytic chain-rule
  expressions for binary `∂(n α^r)/∂n_i` at fixed `(T,V,n_j)`.
- 2026-03-03: Migrated `linear_model_codex.py` Table-1 second-derivative path to analytic derivatives by adding:
  `alpha0_idaes_with_second_derivs`, `alphar_idaes_with_second_derivs`,
  `bell2023_departure_alphar_second_derivs`, and `mixture_alpha0_alphar_second_derivs`.
- 2026-03-03: Removed finite-difference second-derivative fallback from `compute_table1_properties` for `Cv/Cp/speed_of_sound`; retained `fd_rel` argument only for backward API compatibility metadata.

## Assumptions
- Date format uses ISO (`YYYY-MM-DD`) unless a different format is explicitly required.
- This repository is the authoritative workspace for context tracking.
- Fluid names map to IDAES Helmholtz JSON files available on the configured parameter path.
- During differentiation in `(tau, delta)`, composition is treated as constant.
- Mixture reducing functions are treated as composition-only for this derivative path (no explicit `∂/∂x` coupling terms yet).
- External API composition basis is mass fraction for usability; internal EOS/reducing rules continue to use mole fraction.
- Saturation solver currently targets pure fluids only and operates in molar density basis internally.
- Auxiliary saturation density fits in JSON are used as initial guesses and as final fallback if nonlinear root solves fail.

## Known Limitations
- Existing files may not yet conform to the new header/docstring/breadcrumb standards.
- Historical architectural decisions made before this context file may be incomplete until backfilled.
- Explicit composition-derivative coupling terms are not yet implemented in mixture derivative treatment.
- Two-phase region plotting currently depends on externally prepared saturation data; in-repo equilibrium flash generation is not yet implemented.
- Robust fallback paths remain for recovery and diagnostics, but those points are explicitly marked and excluded from converged dome branches.
- Near-critical continuation can still produce occasional `DIVERGED` states depending on step history; these are now logged and excluded from clean dome plots.
- Fixed-composition pseudo-dome is not true mixture VLE; phase compositions are not split and fugacity-equality per component is not enforced.
- True-VLE experimental solver currently has weak dew-branch robustness over broad temperature ranges; many dew states diverge without additional continuation/regularization.
- True-VLE chemical potentials are finite-difference estimates (not closed-form analytic composition derivatives), which can affect conditioning.
- `compute_table1_properties` still carries `fd_rel` for API compatibility, but it is no longer used in analytic derivative calculations.
- True-VLE branch robustness remains the main unresolved numerical limitation; low-temperature pressure bias vs Honeywell references persists despite improved residual closure.

## Validation Status
- Context file bootstrap complete.
- Full repository compliance audit pending.
- `query_ph_codex.py` API/CLI wiring is implemented.
- Strict derivative mapping path is implemented with fixed-composition chain-rule transforms; benchmark validation remains pending.
- Saturation solver artifacts generated:
  - `verification/dome.csv`
- `verification/ph_dome.pdf`
- `verification/saturation_run_metadata.json`
- `verification/saturation_review.md`
- Analytic derivative audit (first and second derivatives) complete for residual Helmholtz pathways in `linear_model_codex.py` with CSV artifacts:
  - `verification/residual_helmholtz_audit_details.csv`
  - `verification/residual_helmholtz_audit_summary.csv`
  and max relative errors up to `1.146e-08` for audited second-derivative mappings.
- Unit test status for saturation module:
  - `tests/test_saturation_solver.py`: 3 passed (strict residual acceptance in single-point check).
- Strict-gate dome run (UTC 2026-03-02T21:35:28.932838+00:00):
  - attempted: 200
  - converged: 199
  - failed/diverged: 1
  - output CSV: `verification/saturation_dome_run.csv`
  - clean plot: `verification/ph_dome_clean.png`
  - diagnostic plot: `verification/ph_dome_with_failures.png`
  - solver trace log: `verification/log_saturation_solver.txt`
- True-VLE copy run (UTC 2026-03-02T21:54:10.006931+00:00):
  - fluid pair: `r1234ze/r227ea`, `w1=0.911`, `z1=0.9385037942 mol/mol`
  - grid: `T=250..360 K`, `n=40`
  - bubble converged: `38/40`
  - dew converged: `19/40`
  - outputs:
    `verification/r515b_true_vle_bubble_copy.csv`,
    `verification/r515b_true_vle_dew_copy.csv`,
    `verification/r515b_true_vle_envelope_copy.png`,
    `verification/r515b_true_vle_metadata_copy.json`.

## Open Questions
- Which existing modules should be prioritized first for compliance retrofitting?
- Should unit conventions be consolidated into a single shared module/spec document?
- When should strict chain-rule mapping be implemented for pure-fluid-to-mixture derivative transfer?
- Should composition-coupling derivative terms (`∂/∂x` effects through reducing functions) be added in the next model extension?

## Next Concrete Steps
1. Audit existing source files for header/docstring/breadcrumb compliance.
2. Backfill missing architectural decisions and rejected alternatives into this file.
3. Add explicit unit-system declarations in relevant computational modules.
4. Enforce modular separation for properties, grids, plotting, and validation during upcoming edits.
5. Add regression tests for known state points to validate the implemented chain-rule mapping path.
6. Evaluate and, if needed, implement composition-coupling derivative terms in a controlled extension.
7. Add a phase-aware saturation/flash pipeline if direct in-code generation of dome boundaries is required.
8. Investigate the remaining near-critical diverged point(s) with adaptive local temperature-step refinement and improved critical-region conditioning.
9. Add regression checks that enforce no failed statuses in target operating temperature windows used for downstream P-H work.
10. Improve true-VLE dew-branch continuation with pressure-stepping or homotopy in composition and add failure-marked diagnostics before using outputs for signoff.
11. Replace finite-difference composition derivatives used for fugacity/chemical potentials with analytic mixture composition derivatives for production VLE stability.
12. Promote mixture flash/VLE solver from exploratory copy modules to a production module with explicit equations:
    `P_l=P_v`, `mu_i^l=mu_i^v` (all components), and material balance constraints.
13. Implement robust bubble/dew continuation (adaptive temperature step, predictor-corrector, and near-critical stop criteria) to obtain continuous converged envelope branches.
14. Produce standard VLE/saturation deliverables for each run:
    converged-only envelope plot, failure-marked diagnostic plot, and per-T status CSV with residuals.
15. Run formal validation against Honeywell/NIST datasets and record per-state pressure/enthalpy errors plus bubble/dew pressure-curve errors.
16. Add dedicated automated tests for VLE path:
    endpoint consistency (`x->0/1`), fugacity-equality residual checks, and regression baselines for agreed reference states.
17. Continue strict session discipline: read/update `PROJECT_CONTEXT.md` at start/end of each significant VLE/dome iteration with solver settings, tolerances, and failure diagnostics.

## Change Log (Chronological)
- 2026-03-02: Created `PROJECT_CONTEXT.md` with required persistent structure and initial baseline entries for project discipline adoption.
- 2026-03-02: Added point-definition breadcrumb `(T, rho_mass, x1)` and standardized returned units `(kPa, kJ/kg)` for p,h query interface.
- 2026-03-02: Recorded strict blocker status: end-to-end execution currently fails by design at `mixture_alpha0_alphar_derivs` `NotImplementedError`.
- 2026-03-02: Added `query_ph_codex.py` API/CLI wrapper with explicit NotImplementedError reporting and no non-strict fallback.
- 2026-03-02: Expanded docstrings in `linear_model_codex.py` with units, assumptions, failure modes, JSON schema expectations, and implementation-status breadcrumbs.
- 2026-03-02: Implemented `mixture_alpha0_alphar_derivs` using chain-rule mapping from pure-fluid derivatives to mixture reduced variables:
  `tau_i = (Tc_i/Tred)*tau`, `delta_i = (rho_red/rho_c_i)*delta`,
  `d(alpha)/d(tau) = d(alpha)/d(tau_i)*(Tc_i/Tred)`,
  `d(alpha)/d(delta) = d(alpha)/d(delta_i)*(rho_red/rho_c_i)`.
- 2026-03-02: Updated enthalpy ideal contribution in `compute_pressure_enthalpy` to use mapped mixture ideal derivative (`H0/RT = 1 + tau*a0_tau`), equivalent to pure-fluid weighted form under fixed composition.
- 2026-03-02: Documented current derivative scope: fixed composition, reducing functions treated as composition-only, no explicit `∂/∂x` coupling yet.
- 2026-03-02: Changed `query_ph_codex.py` interface interpretation of `x1` to mass fraction [kg/kg] and added conversion to mole fraction via
  `x1_mole = (w1/MW1) / ((w1/MW1) + (w2/MW2))` prior to `compute_pressure_enthalpy`.
- 2026-03-02: Updated `query_ph_codex.py` caller interface naming from `x1` to `w1` (`--w1` on CLI) to make mass-fraction basis explicit to users.
- 2026-03-02: Updated `linear_model_codex.py` CLI to accept either `--x1` (mole fraction) or `--w1` (mass fraction), with explicit `w1 -> x1` conversion for consistency with wrapper behavior.
- 2026-03-02: Added `plot_two_phase_dome_codex.py` to visualize two-phase P-H regions from external saturation datasets and optional quality lines.
- 2026-03-02: Removed `plot_two_phase_dome_codex.py` per user direction to scrap plotting code additions.
- 2026-03-02: Added `helmholtz_saturation.py` with pure-fluid saturation solver (`saturation_point_at_T`), dome builder (`compute_saturation_dome`), CSV writer, and P-H dome plotting helper.
- 2026-03-02: Added `scripts/plot_ph_dome.py` CLI wrapper and `tests/test_saturation_solver.py` pytest coverage for derivative consistency and solver behavior.
- 2026-03-02: Generated verification artifacts under `verification/` and recorded run metadata for fluid `r1234ze` over `T=240..372 K` (`n=100`).
- 2026-03-02: Review commit reference: `c8d0fdc`.
- 2026-03-02: Refactored `helmholtz_saturation.py` solver acceptance to strict gates (`r_P`, `r_mu`, density ordering) and disabled fallback-as-converged behavior.
- 2026-03-02: Added explicit solver statuses (`CONVERGED`, `FAILED_CHECKS`, `DIVERGED`), per-temperature CSV diagnostics (`verification/saturation_dome_run.csv`), and failed-point trace log (`verification/log_saturation_solver.txt`).
- 2026-03-02: Added clean/diagnostic dome plots:
  `verification/ph_dome_clean.png`,
  `verification/ph_dome_with_failures.png`.
- 2026-03-02: Strict run summary (`r1234ze`, `T=240..372 K`, `n=200`, `maxiter=50`, gates at `1e-6`):
  converged `199/200`; worst failed point `T=370.010050 K` with
  `r_P=3.1735e-01`, `r_mu=9.4479e-02`, `jac_cond=1.2813e+16`, `step_norm=9.5649e-02`.
- 2026-03-02: Added `mixture_pseudo_dome.py` and copied it to `mixture_pseudo_dome_copy.py` for approximation-only R515B pseudo-dome testing.
- 2026-03-02: Ran copy-based R515B pseudo-dome (`w1=0.911`, `T=240..372 K`, `n=180`): `177` converged, `3` diverged near upper-temperature end.
  Outputs:
  `verification/r515b_pseudo_dome_copy_run.csv`,
  `verification/r515b_pseudo_dome_copy_clean.png`,
  `verification/r515b_pseudo_dome_copy_with_failures.png`,
  `verification/r515b_pseudo_dome_copy_metadata.json`.
- 2026-03-02: Added `mixture_true_vle_copy.py` and executed true-VLE exploratory run (`w1=0.911`, `T=250..360 K`, `n=40`).
  Bubble branch: `38` converged / `2` failed.
  Dew branch: `19` converged / `21` failed.
  Generated:
  `verification/r515b_true_vle_bubble_copy.csv`,
  `verification/r515b_true_vle_dew_copy.csv`,
  `verification/r515b_true_vle_envelope_copy.png`,
  `verification/r515b_true_vle_metadata_copy.json`.
- 2026-03-02: Updated `linear_model_codex.compute_pressure_enthalpy` signature from mole-fraction input (`x1`) to mass-fraction input (`w1`) and added internal conversion helper `x1_from_w1`.
- 2026-03-02: Updated dependent call sites in `query_ph_codex.py` and `third_model_codex.py` to pass `w1` directly (removed duplicate external conversion logic).
- 2026-03-02: Smoke checks after interface update:
  `linear_model_codex.py --w1 0.911 --T 433.15 --rho 608.702` -> `p_kPa=5146.19`, `h_kJkg=496.164`;
  `query_ph_codex.py` returned matching values.
- 2026-03-02: Implemented strict Lemmon-style mixture enthalpy expression in `linear_model_codex.compute_pressure_enthalpy`
  using explicit component-wise `tau_i` ideal derivatives and explicit residual mixture derivative terms.
  Post-change smoke check at `w1=0.911`, `T=433.15 K`, `rho=608.702 kg/m^3` remained:
  `p_kPa=5146.1876`, `h_kJkg=496.1639`.
- 2026-03-02: Added `compute_table1_properties` and supporting helper functions in `linear_model_codex.py` for Lemmon/Tillner-Roth Table-1 style property reporting, including fugacity/chemical-potential terms.
- 2026-03-02: Verified syntax and smoke run for Table-1 path at
  `w1=0.911`, `T=433.15 K`, `rho=608.702 kg/m^3`; representative outputs:
  `Z=0.27580`, `p_kPa=5146.19`, `h_kJkg=496.164`, `phi1=0.6003`, `phi2=0.5108`.
- 2026-03-02: Corrected Table-1 expressions per review:
  `s0/R = h0/(RT) - alpha0_mix - 1` and
  `cv0/R = -sum_i x_i*tau_i^2*(d2 alpha_i^0 / d tau_i^2)`.
  Implemented component-wise finite-difference `alpha_i^0` second derivatives for `cv0`.
- 2026-03-02: Fixed residual Helmholtz mapping bug (`phi_residual_type` 2/3/4) by switching from constant `t[1]` usage to term-wise `t[i]`.
  This changed reference-state outputs at `w1=0.911`, `T=433.15 K`, `rho=608.702 kg/m^3` from:
  previous `p_kPa=5146.1876`, `h_kJkg=496.1639`
  to corrected mapping `p_kPa=8849.8207`, `h_kJkg=446.1515`.
- 2026-03-02: Implemented analytic `∂(n α^r)/∂n_i` replacement for fugacity terms in `compute_table1_properties`.
  Validation against prior finite-difference reference at `w1=0.911`, `T=433.15 K`, `rho=608.702 kg/m^3`:
  analytic `d_nares_dn1=-1.3515584912843746` vs FD `-1.35155849117985` (relative diff `7.7e-11`);
  analytic `d_nares_dn2=-1.388274778060891` vs FD `-1.3882747739473285` (relative diff `3.0e-9`).
- 2026-03-03: Investigated R515A true-VLE root-solver failures in `mixture_true_vle_copy.py`.
  Observed `scipy.optimize.root(method='hybr')` divergence into nonphysical states, causing EOS overflow in `alphar_idaes_with_derivs` (`delta**di` overflow) during vapor-density iterations.
  Actionable next step: replace unconstrained HYBR solve with bounded/damped least-squares (`least_squares` with variable bounds in transformed space) and strict non-finite-state rejection before EOS evaluation.
- 2026-03-03: Added `scripts/validate_r515a_saturation.py` to compare model predictions against NIST TN 2063 Table-7 R515A reference point at `T=277.6 K`, including SI and IP outputs (kPa/psia and kJ/kg/Btu-lbm).
- 2026-03-03: Executed R515A validation using `scripts/validate_r515a_saturation.py` against NIST TN 2063 Table-7 point (`T=277.6 K`, `w1=0.88`).
  Solver-based VLE results: bubble `P=262.32 kPa` (`+3.97%` vs `252.31 kPa`), dew `P=261.49 kPa` (`+3.64%`), and latent heat `171.89..172.46 kJ/kg` (`-2.02%..-1.69%` vs `175.43 kJ/kg`).
  Reference-density pressure cross-check still shows large liquid-side inconsistency (`P_liq≈8471.93 kPa` at reference `rho_l`), indicating remaining model/state-definition inconsistency on dense-liquid branch despite converged bubble/dew residual gates.
- 2026-03-03: Added `scripts/plot_r515a_validation_errorbars.py` and generated
  `verification/r515a_validation_error_bars.png` to visualize model-vs-reference
  error bars at R515A `T=277.6 K` for pressure and latent heat (bubble/dew solver outputs vs NIST TN 2063 reference).
- 2026-03-03: Added consolidated reference-suite validation script `scripts/validate_r515a_reference_suite.py` and executed it.
  Output: `verification/r515a_reference_suite_validation.csv` with 7 checks across Bell 2023 Table 13 (pure and binary alpha_r checks) and NIST TN 2063 Table-7 R515A saturation metrics.
  Key findings: pure alpha_r checks for R1234ze and R227ea match to machine precision; binary alpha_r check at `z1=0.4` has `-1.014%` relative deviation; R515A saturation pressure mismatch remains `+3.64%..+3.97%`, latent heat mismatch `-1.69%..-2.02%` at `T=277.6 K`.
- 2026-03-03: Generated updated R515B true-VLE dome outputs using `mixture_true_vle_copy.py` at `w1=0.911` over `T=250..360 K` (`n=120`).
  Raw outputs: `verification/r515b_true_vle_bubble_20260303.csv`, `verification/r515b_true_vle_dew_20260303.csv`, `verification/r515b_true_vle_envelope_raw_20260303.png`, `verification/r515b_true_vle_metadata_20260303.json`.
- 2026-03-03: Postprocessed run with `scripts/postprocess_true_vle_outputs.py` to produce clean converged and diagnostic artifacts:
  `verification/r515b_true_vle_diagnostics_20260303.csv`,
  `verification/r515b_true_vle_converged_only_20260303.csv`,
  `verification/r515b_true_vle_envelope_clean_20260303.png`,
  `verification/r515b_true_vle_envelope_with_failures_20260303.png`.
  Convergence summary: bubble `118/120` converged, dew `91/120` converged (overall `209/240`).
- 2026-03-03: Ran direct comparison against Honeywell Solstice N15 (R-515B) PT chart values (49 points, `100..2500 kPa`) at chart temperatures.
  Generated `verification/r515b_honeywell_pt_comparison_20260303.csv` with bubble/dew statuses and model pressures.
  Result summary: fully converged at 42/49 chart points (`17.98..90.58 °C`); notable positive pressure bias at lower converged temperatures (e.g., ~`+41.8%` at `17.98 °C`) decreasing toward high temperatures (~`+2.1%` at `90.58 °C`).
- 2026-03-03: Created visual diagnostics for PT comparison:
  `verification/r515b_honeywell_pt_overlay_20260303.png` and
  `verification/r515b_honeywell_pt_error_20260303.png`.
- 2026-03-03: Implemented robust bounded nonlinear solve in `mixture_true_vle_copy.py` using `scipy.optimize.least_squares(method='trf')` replacing unconstrained `root(hybr)` for bubble/dew solves.
  Architectural change details:
  - Reparameterized densities as `(rho_v, rho_l-rho_v)` to enforce `rho_l > rho_v`.
  - Added explicit bounds on state variables (`rho` bounded to `1e-9..2e4 mol/m^3`).
  - Added non-finite/exception residual penalty handling to avoid EOS overflow crashes.
- 2026-03-03: Trial outcome for Honeywell PT sweep (49 chart points) with new least-squares solver:
  `verification/r515b_honeywell_pt_comparison_leastsq_20260303.csv`.
  Fully converged points dropped to `13/49` under strict gates (`r_P<=1e-6`, `r_mu<=1e-6`), compared with prior run `42/49`.
  Interpretation: bounded solver improves numerical stability but currently underperforms on strict chemical-potential convergence with present residual scaling/FD chemical-potential noise; further tuning required before adopting as default.
- 2026-03-03: Added error-figure diagnostics for Honeywell PT comparison runs:
  `verification/r515b_honeywell_pt_error_original_20260303.png`,
  `verification/r515b_honeywell_pt_error_leastsq_20260303.png`, and
  `verification/r515b_honeywell_pt_error_compare_20260303.png`.
- 2026-03-03: Generated R515A p-H envelope artifacts with current solver (`w1=0.88`, `T=250..360 K`, `n=120`):
  `verification/r515a_true_vle_bubble_20260303.csv`,
  `verification/r515a_true_vle_dew_20260303.csv`,
  `verification/r515a_true_vle_envelope_raw_20260303.png`,
  `verification/r515a_true_vle_metadata_20260303.json`.
- 2026-03-03: Postprocessed R515A outputs:
  `verification/r515a_true_vle_diagnostics_20260303.csv`,
  `verification/r515a_true_vle_converged_only_20260303.csv`,
  `verification/r515a_true_vle_envelope_clean_20260303.png`,
  `verification/r515a_true_vle_envelope_with_failures_20260303.png`.
  Convergence summary: bubble `64/120`, dew `110/120` (overall `174/240` rows).
- 2026-03-03: Re-ran consolidated R515A validation suite (`verification/r515a_reference_suite_validation.csv`):
  Bell pure check values match to machine precision; Bell binary `alpha_r` check remains `-1.014%` relative deviation;
  NIST TN 2063 Table-7 saturation at `277.6 K` gives pressure mismatch `+3.64%..+3.97%` and latent-heat mismatch `-1.69%..-2.02%`.
- 2026-03-03: Added R515A trust-map artifacts for interpretability:
  `scripts/plot_r515a_trust_map.py`,
  `verification/r515a_trust_map_20260303.png`, and
  `verification/r515a_trust_map_20260303.csv`.
  Trust classes are currently heuristic:
  GREEN = both branches converged with strict residuals,
  YELLOW = partial/looser convergence,
  RED = no branch convergence.
  External check anchor at `T=277.6 K` includes available NIST saturation errors in plot annotation.
- 2026-03-03: Implemented analytic chemical-potential replacement in `mixture_true_vle_copy.py`:
  replaced finite-difference `mu_i` from `A(T,V,n)` with fugacity-based analytic expressions using
  `d(n*alpha^r)/dn_i` and Bell reducing/departure composition derivatives (aligned with Table-1 style formulation).
  Updated solver residual evaluations to use `chemical_potentials_analytic`.
- 2026-03-03: Post-change checks:
  - R515A reference suite (`verification/r515a_reference_suite_validation.csv`) retained prior pressure/latent-heat mismatch at 277.6 K, but achieved tighter residual closure (`r_mu` near machine precision).
  - Honeywell R515B PT sweep with analytic `mu` (`verification/r515b_honeywell_pt_comparison_analyticmu_20260303.csv`) improved strict convergence from `13/49` (least-squares + FD mu) and `42/49` (older run) to `48/49` points.
  - Remaining issue: strong positive pressure bias at low temperatures persists (model/inputs consistency issue, not solver-equilibrium closure issue).
- 2026-03-03: Completed full docstring compliance pass for VLE path and supporting model structures.
  Updated all top-level functions/classes in `mixture_true_vle_copy.py` to include Purpose, Inputs/Outputs (with units), Assumptions, Failure modes, References, and numerical stability notes.
  Also added missing class docstring for `Bell2023PairParams` in `linear_model_codex.py`.
  Verification: `python3 -m py_compile mixture_true_vle_copy.py linear_model_codex.py` and AST docstring audit report no missing top-level docstrings.
- 2026-03-03: Evaluated denominator-rescaling experiment in VLE residuals (pressure and chemical-potential normalization tweaks) to test low-temperature behavior sensitivity.
  Outcome: modified denominators degraded low-temperature solver robustness and produced broad divergence in a targeted Honeywell low-T check.
  Action taken: reverted denominator scaling to prior form while retaining analytic chemical-potential formulation.
  Verified post-revert baseline on Honeywell PT sweep: `48/49` converged, matching pre-experiment analytic-mu behavior.
- 2026-03-03: Design-gap review (stepwise walkthrough) in progress.
  D1 disposition approved by reviewer: keep model hardwired to R1234ze(E)/R227ea for now (intentional scope constraint).
  D2 approved: identified branch-coupled same-temperature continuation as a numerical design gap.
  Current behavior updates dew seeds from bubble convergence at the same T; agreed direction is to split continuation state by branch (separate bubble/dew seed tracks) in future implementation.
- 2026-03-03: D4 action completed per reviewer instruction: removed hardcoded plot titles from VLE plotting paths to avoid mislabeling across blends.
  Updated files:
  `mixture_true_vle_copy.py` (main envelope plot title removed),
  `scripts/postprocess_true_vle_outputs.py` (clean/failure plot titles removed).
- 2026-03-03: Executed residual Helmholtz audit with reproducible script `scripts/audit_residual_helmholtz.py`.
  Outputs:
  `verification/residual_helmholtz_audit_details.csv`,
  `verification/residual_helmholtz_audit_summary.csv`.
  Audit scope: pure `alphar` derivatives, Bell departure derivatives, mixture `ar_tau/ar_del` mapping, and composition derivatives `d(n*ar)/dn_i`.
  Results (max relative error):
  - pure_alphar: `1.406e-09`
  - departure_alphar: `1.117e-09`
  - mixture_ar_tau_delta: `2.998e-09`
  - composition_d_nar_dn: `1.026e-09`
  Interpretation: residual Helmholtz implementation and first-derivative mappings are internally consistent with finite differences across audited state ranges; remaining low-temperature PT mismatch is unlikely to be caused by first-derivative coding errors in residual Helmholtz terms.
- 2026-03-03: Completed analytic second-derivative implementation in `linear_model_codex.py` for Table-1 closure.
  Added analytic functions:
  `bell2023_departure_alphar_second_derivs` and `mixture_alpha0_alphar_second_derivs`,
  and switched `compute_table1_properties` to analytic `a_tt^r`, `a_dd^r`, `a_td^r` paths.
  Smoke check at `w1=0.911`, `T=433.15 K`, `rho=608.702 kg/m^3` remains numerically stable:
  `p_kPa=8849.8207`, `h_kJkg=446.1515`.
- 2026-03-03: Extended `scripts/audit_residual_helmholtz.py` to include second-derivative audits for:
  pure residual (`alphar_tautau`, `alphar_taudel`, `alphar_deldel`),
  Bell departure second derivatives, and
  mixture mapped second derivatives.
  Re-ran audit and regenerated:
  `verification/residual_helmholtz_audit_details.csv` and
  `verification/residual_helmholtz_audit_summary.csv`.
  Summary maxima (relative):
  `pure_alphar_second=6.509e-09`,
  `departure_alphar_second=2.460e-09`,
  `mixture_ar_second=1.146e-08`.
  Interpretation: analytic first/second derivative implementations are internally consistent with finite-difference checks across sampled states.
- 2026-03-03: Re-read `PROJECT_CONTEXT.md` before Honeywell PT reproduction run (standing instruction compliance).
- 2026-03-03: Reproduced R515B Honeywell PT comparison using current code state and the existing 49-point chart reference set from
  `verification/r515b_honeywell_pt_comparison_analyticmu_20260303.csv` (`T_C`, `P_chart_kPa` columns).
  Ran `mixture_true_vle_copy.run_true_vle_envelope('r1234ze','r227ea', w1=0.911, T_vals)` and saved refreshed output:
  `verification/r515b_honeywell_pt_comparison_20260303_refresh.csv`.
  Result summary (both branches converged): `42/49` points; converged-point MAPE `19.02%`, mean bias `-18.45%`, max absolute percent error `42.77%`.
  Interpretation: current branch-level VLE run exhibits significantly lower modeled pressure than the Honeywell chart over many converged points; this is materially worse than earlier analytic-mu comparison artifacts and indicates a regression/path inconsistency in the current reproduction workflow.
- 2026-03-03: Read `PROJECT_CONTEXT.md` before plotting (standing instruction compliance).
- 2026-03-03: Regenerated Honeywell PT error visualization from refreshed comparison CSV:
  input `verification/r515b_honeywell_pt_comparison_20260303_refresh.csv`;
  outputs:
  `verification/r515b_honeywell_pt_error_refresh_20260303.png` and
  `verification/r515b_honeywell_pt_overlay_refresh_20260303.png`.
  Plot summary: `42/49` fully converged points; converged-point MAPE `19.0188%`; bias `-18.4532%`.
- 2026-03-03: Rolled back the latest analytic-derivative migration in `linear_model_codex.py` per user instruction.
  Restored pre-migration behavior:
  - `alpha0_idaes_with_derivs` and `alphar_idaes_with_derivs` as first-derivative evaluators.
  - `bell2023_departure_alphar` first-derivative form.
  - finite-difference helper `_fd_2d` reinstated.
  - `compute_table1_properties` reverted to finite-difference path for second derivatives (`a_tt^r`, `a_dd^r`, `a_td^r`) and finite-difference fugacity composition derivatives (`d(n*alpha^r)/dn_i`).
  Smoke check after rollback at `w1=0.911`, `T=433.15 K`, `rho=608.702 kg/m^3`:
  `p_kPa=8849.8207`, `h_kJkg=446.1515`.
- 2026-03-03: Reproduced Honeywell R515B PT comparison after derivative rollback using
  `mixture_true_vle_copy.run_true_vle_envelope` and the 49-point chart reference set from
  `verification/r515b_honeywell_pt_comparison_analyticmu_20260303.csv`.
  Output: `verification/r515b_honeywell_pt_comparison_rollback_20260303.csv`.
  Metrics: `42/49` fully converged, converged-point MAPE `19.0188%`, max absolute percent error `42.7712%`, bias `-18.4532%`.
  Conclusion: rollback did not recover Honeywell PT alignment; mismatch pattern is unchanged relative to the immediate pre-rollback refresh run.
- 2026-03-03: Experimented with pressure residual scaling rollback in `mixture_true_vle_copy.py` nonlinear solves:
  changed bubble/dew residual equation from
  `r1 = (P_l - P_v)/max(1,0.5*(P_l+P_v))`
  to unscaled `r1 = (P_l - P_v)` (kept strict acceptance metric `r_P` unchanged for reporting).
- 2026-03-03: Post-change Honeywell R515B PT run (49 chart points) written to
  `verification/r515b_honeywell_pt_comparison_pressure_noscale_20260303.csv`.
  Outcome degraded materially: `24/49` fully converged, converged-point MAPE `99.974%`,
  mean bias `-99.974%`, pressures collapsing near ~100 kPa over broad T range.
  Conclusion: removing pressure normalization in residual equations is not viable for this solver setup.
- 2026-03-03: Undid the D5 pressure-residual rollback experiment in `mixture_true_vle_copy.py`.
  Restored residual normalization in bubble/dew solves to:
  `r1 = (P_l - P_v)/max(1, 0.5*(P_l + P_v))`.
- 2026-03-03: Re-ran Honeywell PT comparison after restoring scaling:
  `verification/r515b_honeywell_pt_comparison_restored_scaling_20260303.csv`.
  Metrics returned to prior baseline: `42/49` fully converged, converged-point MAPE `19.0188%`, bias `-18.4532%`.
- 2026-03-03: Attempted hard rollback to checkpoint commit `69f444a` for
  `linear_model_codex.py` and `mixture_true_vle_copy.py` to recover prior `48/49` Honeywell PT behavior.
  Reproduction run output: `verification/r515b_honeywell_pt_comparison_reverted_to_69f444a_20260303.csv`.
  Result did not recover target: `35/49` both-branch converged (`bubble 45/49`, `dew 37/49`).
  Interpretation: the previously observed `48/49` state was not equivalent to this checkpoint snapshot alone (likely depended on later uncommitted solver-state/code path not captured by this commit).
- 2026-03-03: User direction updated: no thermodynamic/science changes; only convergence-focused numerical adjustments are allowed for current iteration (solver behavior, FD step sizing, continuation strategy).
- 2026-03-03: Attempted restore to prior checkpoint `69f444a` for
  `linear_model_codex.py` and `mixture_true_vle_copy.py` and re-ran Honeywell PT comparison.
  Output: `verification/r515b_honeywell_pt_comparison_reverted_to_69f444a_20260303.csv`.
  Result: `35/49` both-branch converged (`bubble 45/49`, `dew 37/49`), which did not recover previously observed `48/49`.
- 2026-03-03: Current tuning objective set to recover/improve PT convergence from this restored baseline using convergence-only edits (no EOS/science formula changes).
- 2026-03-03: Restored pre-`69f444a` 42/49 baseline file states from local git object snapshots:
  `mixture_true_vle_copy.py` <- blob `3b9735f02d84fc83e79b5a16b41926950d737a55`
  `linear_model_codex.py` <- blob `6cf521a0353523044991ccd38b9b73575754dce0`.
  Verification run output: `verification/r515b_honeywell_pt_comparison_restored_42baseline_20260303.csv`.
  Result: `42/49` both-branch converged (`bubble 48/49`, `dew 42/49`), matching the previously observed 42-baseline pattern.
- 2026-03-03: Started convergence-only FD/Jacobian tuning from confirmed 42/49 baseline (no thermodynamic equation changes).
  Active path note: VLE solver currently uses analytic `chemical_potentials_analytic`; FD tuning targeted nonlinear solver Jacobian differencing and scaling only.
- 2026-03-03: Added solver numeric knobs in `mixture_true_vle_copy.py`:
  `LSQ_DIFF_STEP`, `LSQ_MAX_NFEV`, `LSQ_METHOD`, `LSQ_X_SCALE`.
  Performed sweeps:
  - `verification/r515b_fdstep_sweep_20260303.csv`
  - `verification/r515b_fdstep_sweep_small_20260303.csv`
  - `verification/r515b_nfev_sweep_20260303.csv`
  - `verification/r515b_lsq_numerics_sweep_20260303.csv`
- 2026-03-03: Best convergence-only setting selected from sweep:
  `LSQ_METHOD='trf'`, `LSQ_X_SCALE=1.0`, `LSQ_DIFF_STEP=1e-7`, `LSQ_MAX_NFEV=800`.
  Confirmation run output: `verification/r515b_honeywell_pt_comparison_tuned_solver_20260303.csv`.
  Result: `49/49` both-branch converged (`bubble 49/49`, `dew 49/49`), with converged-point MAPE `16.1491%`.
- 2026-03-03: Generated tuned Honeywell PT diagnostic plots from
  `verification/r515b_honeywell_pt_comparison_tuned_solver_20260303.csv`:
  `verification/r515b_honeywell_pt_overlay_tuned_20260303.png` and
  `verification/r515b_honeywell_pt_error_tuned_20260303.png`.
  Converged points: `49/49`; MAPE `16.1491%`.
- 2026-03-03: Added repository hygiene ignore rules to `.gitignore` for generated diagnostics/plot artifacts (verification images and top-level generated isotherm/isochor/PT mismatch plot files). Existing `*.log` ignore remains in place for log files.
- 2026-03-03: Updated `.gitignore` to also ignore CSV outputs (`verification/*.csv` and `*.csv`) per user request for excluding generated tabular artifacts from pushes.
- 2026-03-03: Updated `.gitignore` to ignore JSON outputs as requested (`verification/*.json` and `*.json`) alongside CSV/plot artifacts.
- 2026-03-03: Reverted JSON ignore addition per user correction; `.gitignore` now excludes plots and CSV artifacts, but does not ignore JSON files.
