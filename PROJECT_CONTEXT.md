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
- 2026-03-03: Ran a focused Honeywell PT solver-parameter sweep to check if MAPE can be reduced further without changing thermodynamic equations.
  Output: `verification/r515b_solver_mape_sweep_focused_20260303.csv`.
  Findings:
  - Best raw MAPE observed: `15.0987%` using `LSQ_METHOD='dogbox'`, `LSQ_DIFF_STEP=5e-8` (or `1e-7`), `LSQ_MAX_NFEV=1200`, `LSQ_X_SCALE=1.0`.
  - Tradeoff: this setting only converged `33/49` both-branch points, so it fails the strict full-coverage objective.
  - Best full-coverage setting remains unchanged: `LSQ_METHOD='trf'`, `LSQ_X_SCALE=1.0`, `LSQ_DIFF_STEP=1e-7`, `LSQ_MAX_NFEV=800`, with `49/49` converged and MAPE `16.1491%`.
  Decision: keep `trf` settings for production runs unless user approves reduced-coverage calibration mode.
- 2026-03-03: Implemented convergence-only guess-strategy upgrades in `mixture_true_vle_copy.py`:
  (1) separated continuation seeds by branch (bubble and dew no longer cross-share `rho_l/rho_v`),
  (2) added bounded adaptive retry seeds per point (`rho_l`/`rho_v` perturbations and `x/y` perturbation ±0.01).
  Verification run output: `verification/r515b_honeywell_pt_comparison_guess_tuned_20260303.csv`.
  Result: strict convergence stayed `49/49` (`bubble 49/49`, `dew 49/49`), MAPE remained `16.1491%` (no net change vs tuned baseline).
  Retry usage diagnostics on the 49-point set: bubble used primary seed on all points (`retry=0` for `49/49`), dew used one fallback-retry point (`retry=1` for `1/49`).
  Interpretation: current MAPE floor appears model-dominant under the existing equation set; guess strategy improved robustness margin but not pressure-bias magnitude.
- 2026-03-03: Generated SI cross-check handoff figure for interaction-parameter source review:
  `verification/r515b_bell_si_crosscheck_figure_20260303.png`.
  Panel A shows Bell 2023 Table-13 anchor relative errors from `verification/r515a_reference_suite_validation.csv`;
  Panel B shows current Honeywell R515B PT relative-error trend vs temperature from
  `verification/r515b_honeywell_pt_comparison_guess_tuned_20260303.csv` (`w1=0.911`, `49/49` converged, MAPE `16.15%`).
- 2026-03-03: Generated fresh R515B true-VLE dome artifacts with current tuned solver (no science changes):
  command settings: `w1=0.911`, `T=250..360 K`, `n=200`.
  Outputs:
  - `verification/r515b_true_vle_bubble_dome_20260303.csv`
  - `verification/r515b_true_vle_dew_dome_20260303.csv`
  - `verification/r515b_true_vle_metadata_dome_20260303.json`
  - `verification/r515b_true_vle_envelope_dome_raw_20260303.png`
  Postprocessed diagnostics:
  - `verification/r515b_true_vle_diagnostics_dome_20260303.csv`
  - `verification/r515b_true_vle_converged_only_dome_20260303.csv`
  - `verification/r515b_true_vle_envelope_clean_dome_20260303.png`
  - `verification/r515b_true_vle_envelope_with_failures_dome_20260303.png`
  Convergence summary: bubble `200/200`, dew `200/200`.
- 2026-03-03: Implemented crossover-aware dome postprocessing/visualization updates in
  `scripts/postprocess_true_vle_outputs.py` (no thermodynamic solver logic changes):
  - added `detect_crossover_indices(h_liq, h_vap)` and `get_crossovers(dome_csv_path)`,
  - added crossover diagnostics CSV output and warning path with optional `--strict` non-zero exit,
  - split bubble/dew plotted lines at crossover indices and highlighted inverted tie-lines in red,
  - added companion reviewer plots: P-T overlay, h-T overlay, and P-h composition-colorized view,
  - added optional `--continuous-envelope` diagnostic interpolation proxy (visualization-only; not TP-flash).
  New artifacts from the 2026-03-03 dome run:
  - `verification/r515b_dome_crossovers_20260303.csv`
  - `verification/r515b_true_vle_pt_overlay_dome_20260303.png`
  - `verification/r515b_true_vle_ht_overlay_dome_20260303.png`
  - `verification/r515b_true_vle_ph_colored_dome_20260303.png`
  Crossover count detected: `17`.
  Representative crossover: `T = 306.381910 K`,
  `h_liq = 426.975227 kJ/kg`, `h_vap = 396.888137 kJ/kg`,
  `P_bubble = 491.356425 kPa`, `P_dew = 777.490842 kPa`,
  `x_liq = 0.938504`, `y_vap = 0.938504`,
  action: tie-line plotted in red and row recorded in `verification/r515b_dome_crossovers_20260303.csv`.

## 2026-03-02 — Honeywell comparison run (automated diagnostics)

- Honeywell source: Solstice-N15-TDS_EN.pdf, page 1 (composition) & page 3 (PT table).
- Honeywell composition used for comparison: w1=0.911 (mass fraction), w2=0.089. Source: Solstice-N15-TDS_EN.pdf, p.1.
- Honeywell mass fractions: w1=0.911 (R1234ze), w2=0.089 (R227ea). Converted mole fraction used: x1=0.93850379423, x2=0.0614962057704. MWs used from model JSONs: MW1=0.1140416 kg/mol, MW2=0.17002886 kg/mol.
- Honeywell mixture MW note: datasheet mixture MW (117.48 kg/kmol) logged as reference only; composition conversion used component MWs from IDAES parameter files.
- Bubble comparison CSV: `diagnostics/honeywell_vs_model_bubble_20260302.csv`
- Dew comparison CSV: `diagnostics/honeywell_vs_model_dew_20260302.csv`
- REFPROP/CoolProp comparison CSV: `diagnostics/ref_compare_20260302.csv` (status: not_available for usable saturation points in this environment; CoolProp import succeeded but returned NaN for queried blend saturation pressures at selected temperatures)
- Parameter dump path: `diagnostics/params_used_20260302.txt`
- Single point deep dive JSON: `diagnostics/single_point_deepdive_T_253.93_20260302.json`
- Honeywell-marker overlay figure: `diagnostics/r515b_ph_overlay_honeywell_markers_20260302.png`
- MAPE (bubble vs Honeywell) = 16.418081620111% ; convergence count = 49/49
- MAPE (dew vs Honeywell) = 15.880078702164% ; convergence count = 49/49
- Solver note: all Honeywell-grid points returned `solver_status=CONVERGED`; fugacity residual max values are near machine precision (order 1e-15 in logged CSVs).
- Worst bubble absolute-percent-error deep dive (auto-selected):
  T = 253.93 K (-19.22 C), P_ref = 100.0 kPa,
  P_bubble_model = 180.269170622232 kPa, P_dew_model = 176.393221333737 kPa,
  P_at_liq_density(x_overall) = 180.269170622259 kPa,
  P_at_vap_density(x_overall) = 180.166638373148 kPa,
  ln(f_l/f_v) = [1.7231058903339965e-15, 3.446211780667993e-15].
- Bell parameter sanity result: `diagnostics/params_used_20260302.txt` reports
  `max_abs_diff_reducing=0`, `max_abs_diff_departure=0`, `flag_error_gt_1e-10=False` against Bell (2023) Table 7 literals currently encoded in code.
- Quick conclusion: under current no-science-change model, Honeywell PT bubble/dew MAPE remains ~16%; no parameter-load mismatch was detected in the Bell pair constants/coefficient dump.
- Next action: keep EOS parameters unchanged; if lower MAPE is required, proceed via approved scientific path (reference-tool baseline on an environment with working mixture saturation calls and/or model-parameter recalibration).

- 2026-03-03: Per user-directed sensitivity check, excluded the lowest Honeywell PT point (`T=-19.22 C`, `P_ref=100 kPa`) and recomputed comparison metrics on the remaining 48 points (no solver/model changes).
  Filtered outputs:
  - `diagnostics/honeywell_vs_model_bubble_excl_minus19_20260303.csv`
  - `diagnostics/honeywell_vs_model_dew_excl_minus19_20260303.csv`
  Results:
  - bubble: convergence `48/48`, MAPE `15.087850599234%` (from `16.418081620111%` on 49-point set)
  - dew: convergence `48/48`, MAPE `14.619388230673%`
  Interpretation: excluding `-19.22 C` improves MAPE modestly, but error remains well above single-digit threshold.
- 2026-03-03: User-directed solver-only tuning attempt (explicitly no science edits) expanded least-squares controls in `mixture_true_vle_copy.py`:
  added configurable `LSQ_FTOL`, `LSQ_XTOL`, `LSQ_GTOL`, `LSQ_LOSS`, `LSQ_F_SCALE` (equations/parameters unchanged).
- 2026-03-03: Randomized solver hyperparameter search on Honeywell bubble PT objective (`w1=0.911`, full 49-point set):
  results file: `diagnostics/solver_tuning_bubble_randomsearch_20260303.csv`.
  Best candidate found:
  `method='trf'`, `diff_step=1e-6`, `max_nfev=800`, `x_scale=1.0`,
  `loss='cauchy'`, `f_scale=10.0`, `ftol=xtol=gtol=1e-14`.
  Bubble validation CSV: `diagnostics/honeywell_vs_model_bubble_solver_tuned_20260303.csv`.
  Bubble outcome: convergence `49/49`, MAPE `1.122370545966%` (and `1.050552283380%` excluding `-19.22 C`).
- 2026-03-03: Same tuned solver profile tested on dew branch:
  CSV `diagnostics/honeywell_vs_model_dew_solver_tuned_20260303.csv`,
  convergence `49/49` but MAPE `63.178603412185%`.
  Interpretation: this solver profile is bubble-objective-specific and degrades dew consistency; defaults were not switched globally pending user decision.
- 2026-03-03: Updated density variable mapping in `mixture_true_vle_copy.py` from exponential transform to bounded scaled-sigmoid transform (user-requested stability exploration; no EOS parameter/formula changes).
  New parameterization:
  - `rho_v = lo_v + (hi_v-lo_v)*sigmoid(u0)`
  - `drho  = lo_d + (hi_d-lo_d)*sigmoid(u1)`
  - `rho_l = rho_v + drho` with existing physical ordering/caps.
  Implemented bounds correspond to user-selected upper density target near `2 g/cc` (conservative molar bounds:
  `RHO_MAP_MIN_MOLM3=5.0`, `RHO_MAP_MAX_MOLM3=2.0e4`, `DRHO_MAP_MIN_MOLM3=1.0`, `DRHO_MAP_MAX_MOLM3=2.0e4`).
  Validation run with bubble-tuned solver profile:
  `diagnostics/honeywell_vs_model_bubble_sigmoidmap_20260303.csv` =>
  convergence `49/49`, MAPE `1.122370545963%` (effectively unchanged vs prior bubble-tuned mapping).
- 2026-03-03: Dew validation with scaled-sigmoid mapping and same tuned solver profile:
  `diagnostics/honeywell_vs_model_dew_sigmoidmap_20260303.csv` =>
  convergence `49/49`, MAPE `15.880078702165%` on 49 points and
  `14.619388230673%` excluding `-19.22 C`.
  Interpretation: scaled-sigmoid mapping itself does not fix dew-vs-Honeywell mismatch; dew branch remains materially higher-error than bubble under current no-science-change model.
- 2026-03-03: User requested dew hyperparameter search excluding the `-19.22 C` point; bubble settings kept unchanged.
  Dew-only search artifact:
  `diagnostics/solver_tuning_dewonly_excl_minus19_randomsearch_20260303.csv`.
  Best dew profile found:
  `LSQ_METHOD='trf'`, `LSQ_DIFF_STEP=1e-8`, `LSQ_MAX_NFEV=800`, `LSQ_X_SCALE=2.0`,
  `LSQ_LOSS='huber'`, `LSQ_F_SCALE=3.0`, `LSQ_FTOL=LSQ_XTOL=LSQ_GTOL=1e-10`.
  Validation CSV using that profile on 48-point set:
  `diagnostics/honeywell_vs_model_dew_tuned_excl_minus19_20260303.csv`.
  Dew result: convergence `48/48`, MAPE `0.971571003176%`.
  Plot artifact:
  `diagnostics/honeywell_vs_model_dew_tuned_excl_minus19_20260303_plot.png`.
- 2026-03-03: Tested the same dew-tuned profile with `-19.22 C` added back (full 49-point Honeywell set):
  - CSV: `diagnostics/honeywell_vs_model_dew_tuned_full49_20260303.csv`
  - plot: `diagnostics/honeywell_vs_model_dew_tuned_full49_20260303_plot.png`
  - convergence `49/49`, MAPE `1.021220878991%` (vs `0.971571003176%` on 48-point subset; delta `+0.049649875815` points).
  - `-19.22 C` row under this profile: `P_model=103.404415 kPa` vs `P_ref=100 kPa` (`+3.404415%`).
- 2026-03-03: Generated combined both-branch tuned comparison package on full 49-point Honeywell grid
  (bubble uses bubble-tuned profile; dew uses dew-tuned profile):
  - CSV: `diagnostics/honeywell_vs_model_both_tuned_full49_20260303.csv`
  - plot: `diagnostics/honeywell_vs_model_both_tuned_full49_20260303_plot.png`
  - convergence: bubble `49/49`, dew `49/49`
  - MAPE: bubble `1.122370545963%`, dew `1.021220878991%`, simple average-pressure MAPE `1.071795712477%`.
- 2026-03-03: Generated new R515B saturation-dome package using branch-specific tuned solver profiles
  (bubble profile for bubble branch; dew profile for dew branch) with no EOS/science changes:
  - raw branch CSVs:
    - `verification/r515b_true_vle_bubble_dome_branch_tuned_20260303.csv`
    - `verification/r515b_true_vle_dew_dome_branch_tuned_20260303.csv`
  - metadata:
    - `verification/r515b_true_vle_dome_branch_tuned_metadata_20260303.json`
  - postprocessed diagnostics:
    - `verification/r515b_true_vle_diagnostics_dome_branch_tuned_20260303.csv`
    - `verification/r515b_true_vle_converged_only_dome_branch_tuned_20260303.csv`
    - `verification/r515b_dome_crossovers_branch_tuned_20260303.csv`
  - figures:
    - `verification/r515b_true_vle_envelope_clean_dome_branch_tuned_20260303.png`
    - `verification/r515b_true_vle_envelope_with_failures_dome_branch_tuned_20260303.png`
    - `verification/r515b_true_vle_pt_overlay_dome_branch_tuned_20260303.png`
    - `verification/r515b_true_vle_ht_overlay_dome_branch_tuned_20260303.png`
    - `verification/r515b_true_vle_ph_colored_dome_branch_tuned_20260303.png`
  Run summary: bubble `200/200` converged, dew `200/200` converged, crossover count `0`.
- 2026-03-03: User requested no postprocessing dependency for dome interpretation.
  Generated direct raw dome visualization from branch CSV outputs only (no postprocess transformations):
  `verification/r515b_true_vle_envelope_raw_nopostprocess_20260303.png`.
  Raw-plot source branches: `verification/r515b_true_vle_bubble_dome_branch_tuned_20260303.csv`,
  `verification/r515b_true_vle_dew_dome_branch_tuned_20260303.csv` (both `200/200` converged).

## 2026-03-03 — β-parameterized two-phase band added

- Two-phase mixture enthalpy now computed via lever rule:
    h = (1 - β) h_l + β h_v
- References:
    - Smith, Van Ness & Abbott
    - Prausnitz et al.
    - MIT Unified Thermodynamics Notes (Node 69)
- Purpose:
    Correct representation of mixture VLE envelope in P–h space.
- 2026-03-03: Re-generated saturation dome using branch-tuned solver profiles and in-module β-band plotting:
  - `verification/r515b_true_vle_bubble_dome_beta_20260303.csv`
  - `verification/r515b_true_vle_dew_dome_beta_20260303.csv`
  - `verification/r515b_true_vle_envelope_beta_band_20260303.png`
  - `verification/r515b_true_vle_dome_beta_metadata_20260303.json`
  Convergence: bubble `200/200`, dew `200/200`.

## 2026-03-03 — Frozen true VLE reference module

- Created: `mixture_vle_true_reference.py`
- Copied from: `mixture_true_vle_copy.py`
- Purpose: preserve μ-equality mixture VLE solver for future non-azeotropic mixtures
- Guard: `tests/test_true_vle_reference_frozen.py` (3-point regression)
- Source commit hash at copy time: `b3ec3a5b110deb36eb47e109259ea4a8a6ed6eae`
- 2026-03-03: Ran targeted Honeywell-grid azeotrope logic test at each PT-table temperature:
  bubble solve with `x=z`, dew solve with `y=z`, then logged
  `delta_x1 = y1_from_bubble - x1_from_dew`.
  Output CSV: `diagnostics/honeywell_azeotrope_dx_test_20260303.csv`.
  Summary:
  - points: `49`, converged bubble/dew: `49/49` each
  - `max |delta_x1| = 0.052089715988` (at `T=-19.22 C`)
  - `mean |delta_x1| = 0.008993712315`
  - `p95 |delta_x1| = 0.033902202016`
  Interpretation: `delta_x1` is not numerically tiny across the full table, so the tuned branch solutions are not behaving as an azeotrope-like `x≈y` path over the full range.

## 2026-03-03 — Fixed-Pressure Glide Test

- CSV: `diagnostics/honeywell_glide_test_20260303.csv`
- Max dT_glide: `0.267650315562 K`
- Mean dT_glide: `0.029085725643 K`
- p95 dT_glide: `0.142774483082 K`
- Convergence: `49/49` bubble fix-P solves, `49/49` dew fix-P solves
- Worst pressure point: `P=100.0 kPa` with
  `T_bubble=252.890667853297 K`,
  `T_dew=253.158318168859 K`,
  `dT_glide=0.267650315562 K`.
- Failures: none (`bubble_fixP_status='OK'` for all rows, `dew_fixP_status='OK'` for all rows).

## 2026-03-03 — Mixture Behavior Classification (from fixed-P glide)

- Classification: **Near-azeotropic**
- Rule used:
  - Strict azeotrope if `Max dT_glide < 0.1 K`
  - Near-azeotropic if `0.1 K <= Max dT_glide <= 1 K`
  - Zeotropic if `Max dT_glide > 1 K`
- Computed metrics from `diagnostics/honeywell_glide_test_20260303.csv`:
  - `Max dT_glide = 0.267650315562 K`
  - `Mean dT_glide = 0.029085725643 K`
  - `p95 dT_glide = 0.142774483082 K`

## 2026-03-03 — Fixed-Pressure Composition Split Check

- CSV: `diagnostics/honeywell_glide_dx_check_20260303.csv`
- Definition used: `delta_x1_abs = |y1_bubble(P) - x1_dew(P)|`
- Summary metrics:
  - `max |delta_x1| = 0.053102925668`
  - `mean |delta_x1| = 0.009162151968`
  - `p95 |delta_x1| = 0.034530793824`
- Negligibility assessment:
  - Composition split is **not negligible** over the full pressure range (cold/low-pressure end shows largest separation).
- Correlation with glide:
  - Pearson correlation between `|delta_x1|` and `|dT_glide|` is `0.977504486911` (strong positive correlation).

## 2026-03-03 — Step 4 Charting Mode Decision and Outputs

- Classification input from Step 2: **Near-azeotropic**.
- Applied protocol Case A (optional pseudo-pure saturation mode) in a separate diagnostic script:
  `scripts/pseudopure_mode_eval.py` (keeps frozen reference module untouched).
- Generated artifacts:
  - `diagnostics/r515b_pseudopure_dome_20260303.csv`
  - `diagnostics/r515b_pseudopure_vs_honeywell_20260303.png`
  - `diagnostics/r515b_pseudopure_summary_20260303.json`
- Quantitative deviation metrics (from summary JSON):
  - initial run had `n_total=49`, `n_ok=42`, `n_failed=7`
  - `MAE(T_mid vs Honeywell T) = 0.373356946118 K`
  - `max |T_mid - T_honeywell| = 0.905506988922 K`
  - fixed-pressure error metric in this mode is structurally `0` (pressure matched by construction at each row).
- 2026-03-03: Bracketing refinement applied per user direction in `scripts/pseudopure_mode_eval.py`:
  - first-pass windows in mass-density space:
    vapor `0.001..50 kg/m^3`, liquid `50..max(1500, 3*max(rhoc)) kg/m^3`
  - adaptive fallback full-range sweep (`1e-3..max(liq_hi, 2500) kg/m^3`) selecting first/last sign-change brackets.
  Re-run result:
  - `n_ok=49`, `n_failed=0` (all prior 7 `FAILED_NO_RHO_BRACKET` points at `2200..2500 kPa` resolved).
  - 2400 kPa P(rho) diagnostic written to:
    `diagnostics/pseudopure_prho_sweep_2400kPa_20260303.csv`.
- 2026-03-03: Bug/fix standardization applied with shared utility:
  “Bug: density bracketing bounds were set to ρ∈[0.001,2.0] kg/m³, which excludes liquid-like roots at multi-MPa pressures. Fix: unify density bracketer with vapor+liquid windows and adaptive upper bound based on critical density scale; apply to pseudo-pure dome and any true VLE routines that bracket densities.”
  Implementation:
  - shared utility module: `density_bracketing.py`
  - integrated in: `scripts/pseudopure_mode_eval.py` and `pressure_validated_model.py`
  Defaults now explicit in mass-density units:
  - `rho_min = 1e-4 kg/m^3`
  - `rho_split = 50 kg/m^3`
  - `rho_max = min(2000, 3*max(rhoc_i)) kg/m^3`
  Failure logging now includes final bounds via `format_bounds_log(...)` at bracket failures.
- 2026-03-03: Generated updated Honeywell-vs-model comparison figure from tuned full-49 branch results:
  `diagnostics/honeywell_vs_model_comparison_tuned_full49_20260303.png`.
- 2026-03-03: Regenerated saturation dome after comparison-review approval using current tuned branch profiles and in-module β-band plotting:
  - `verification/r515b_true_vle_bubble_dome_regen_20260303.csv`
  - `verification/r515b_true_vle_dew_dome_regen_20260303.csv`
  - `verification/r515b_true_vle_envelope_regen_20260303.png`
  - `verification/r515b_true_vle_dome_regen_metadata_20260303.json`
  Convergence: bubble `200/200`, dew `200/200`.
- 2026-03-03: Generated clipped-domain P-h view from regenerated dome with axis limits requested for review:
  `h = 150..500 kJ/kg`, `P = 1..100 bar`.
  Output: `verification/r515b_true_vle_envelope_regen_20260303_h150_500_p1_100.png`.
- 2026-03-03: Layer-by-layer plotting protocol STEP 0 executed.
  Created plotting driver scaffold module: `plot_ph_r515b_layers.py`.
  Artifacts:
  - `diagnostics/plots/ph_layer0_scaffold_20260303.png`
  - `diagnostics/plots/ph_layer_status_20260303.json`
  Status JSON records `completed_layers=[0]`, `inputs_used=[]`.
- 2026-03-03: Layer-by-layer plotting protocol STEP 1 executed.
  Honeywell reference data available in current repo is PT-only (no direct Honeywell enthalpy points).
  Artifacts:
  - `diagnostics/plots/ph_layer1_honeywell_reference_20260303.png` (PT reference plot)
  - `diagnostics/plots/ph_layer1_inputs_20260303.csv` (exact columns: `T_C`, `P_ref_kPa`)
  - `diagnostics/plots/ph_layer_status_20260303.json` updated to `completed_layers=[0,1]`.
- 2026-03-03: Layer-by-layer plotting protocol STEP 2 executed (true VLE boundaries only; no new solves).
  Inputs:
  - `verification/r515b_true_vle_bubble_dome_regen_20260303.csv`
  - `verification/r515b_true_vle_dew_dome_regen_20260303.csv`
  - `verification/r515b_true_vle_dome_regen_metadata_20260303.json`
  Artifacts:
  - `diagnostics/plots/ph_layer2_true_vle_boundary_20260303.png`
  - `diagnostics/plots/ph_layer2_true_vle_boundary_20260303.csv`
  - `diagnostics/plots/ph_layer_status_20260303.json` updated to `completed_layers=[0,1,2]`.
  Boundary checks on aligned points (`n=200`):
  - crossings (`h_dew < h_bubble`): `0`
  - minimum enthalpy separation: `18.621116283849 kJ/kg`
  - maximum enthalpy separation: `30.418098616050 kJ/kg`
  - missing/gap indicators in branch T arrays: `0` for bubble, `0` for dew.
- 2026-03-03: Layer-by-layer plotting protocol STEP 3 executed (β-band fill from conjugate pairs).
  Artifacts:
  - `diagnostics/plots/ph_layer3_true_vle_beta_band_20260303.png`
  - `diagnostics/plots/ph_layer3_beta_grid_20260303.csv`
  - `diagnostics/plots/ph_layer_status_20260303.json` updated to `completed_layers=[0,1,2,3]`.
  Implementation notes:
  - β levels used: `0.0, 0.1, ..., 1.0`
  - per pair, `h(beta)=(1-beta)h_l+beta h_v`, `P(beta)=P_pair` (horizontal tie-lines)
  - geometric validity check: `invalid_pairs(h_dew < h_bubble)=0`, `min_h_span=18.621116283849 kJ/kg`.
- 2026-03-03: Layer-by-layer plotting protocol STEP 4 executed.
  Chosen family: **Isotherms** (Option A, first-pass context lines).
  Artifacts:
  - `diagnostics/plots/ph_layer4_isotherms_20260303.png`
  - `diagnostics/plots/ph_layer4_isotherms_20260303.csv`
  - `diagnostics/plots/ph_layer_status_20260303.json` updated to `completed_layers=[0,1,2,3,4]`.
  Construction details:
  - temperature set [C]: `[-19, -5, 10, 25, 40, 55, 70, 85]`
  - vapor-like density window: `0.001..40 kg/m^3`
  - liquid-like density window: `120..2000 kg/m^3`
  - solved points: `2240 OK`, `0 FAILED`
  Two-phase artifact control:
  - no flash was invoked in this layer;
  - curves were generated only from separated vapor-like / liquid-like density windows for context overlays.
- 2026-03-03: Layer-by-layer plotting protocol STEP 5 executed (optional pseudo-pure overlay).
  Artifacts:
  - `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.png`
  - `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.csv`
  - `diagnostics/plots/ph_layer_status_20260303.json` updated to `completed_layers=[0,1,2,3,4,5]`.
  Status summary:
  - pseudo-pure points `OK=49`, `FAILED=0` (status codes included in CSV).
  - pseudo-pure dotted overlay tracks near-azeotropic datasheet-style boundary context over the solved range.
- 2026-03-03: Step-5 pseudo-pure liquid-branch jump discontinuity diagnosed (for ongoing work with TARS).
  Observation from `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.csv`:
  - largest liquid enthalpy drop occurs between `900 -> 950 kPa`:
    `h_l: 356.846010 -> 267.713337 kJ/kg` (`Δh_l = -89.132673 kJ/kg`).
  Root-selection evidence from `diagnostics/r515b_pseudopure_dome_20260303.csv`:
  - at `<=900 kPa`: `bracket_note=WINDOW_OK`, `rho_l ~ 1708..1731 mol/m^3`
  - at `>=950 kPa`: `bracket_note=FALLBACK_FIRST_LAST_SIGN_CHANGE`, `rho_l ~ 9159..8961 mol/m^3`
  Interpretation:
  - discontinuity is caused by liquid-root branch switching during fallback bracketing,
    not by EOS parameter changes.
  Next numeric-only direction (no science edit):
  - apply branch-continuity root selection (prefer liquid root nearest previous
    liquid state in `rho_l`/`h_l`) when multiple sign-change candidates exist.
- 2026-03-03: STEP 5.1 instrumentation completed in `scripts/pseudopure_mode_eval.py` (no root-selection behavior change).
  Added full candidate logging for all sign-change brackets in both windows per pressure:
  - `diagnostics/pseudopure_root_candidates_20260303.csv`
  Fields include:
  - `rho_bracket_lo`, `rho_bracket_hi`, `rho_root`, `h_root_kJkg`, `root_rank`, `phase_window`, `root_status`
  Candidate-count check (problematic region):
  - `900 kPa`: vapor candidates `1`, liquid candidates `4`
  - `950 kPa`: vapor candidates `0`, liquid candidates `5`
  - `1000 kPa`: vapor candidates `0`, liquid candidates `5`
  Confirmation:
  - multiple liquid roots exist in this pressure range; current jump is consistent with branch-selection ambiguity.
- 2026-03-03: STEP 5.2 continuity selector implemented in `scripts/pseudopure_mode_eval.py` (numerical-only selection logic).
  Selection rule:
  - if previous liquid state exists, select liquid root minimizing `|rho_l-rho_l_prev|`
  - first accepted point keeps standard liquid-window bracket root (fallback to highest-density liquid candidate only if no bracket-match candidate is found).
  Guardrail:
  - if `|h_l-h_l_prev| > 20 kJ/kg`, emit warning and reference candidate CSV.
  Re-run outcome (`diagnostics/r515b_pseudopure_dome_20260303.csv`):
  - liquid-side jump discontinuity in 900–1000 kPa region removed.
  - previous jump (`900->950 kPa: -89.13 kJ/kg`) replaced by smooth progression:
    `h_l(900)=356.846`, `h_l(950)=358.118`, `h_l(1000)=359.245 kJ/kg`.
  - largest remaining adjacent liquid-step magnitude across curve: `5.981 kJ/kg`.
  - jump warnings emitted: `0`.
- 2026-03-03: STEP 5.3 completed in `scripts/pseudopure_mode_eval.py`.
  Fallback policy update:
  - fallback remains bracket-discovery only (`bracket_note` retained for diagnostics)
  - root selection is continuity-driven over explicit full-grid root candidates
    (`vapor = lowest-density root`, `liquid = continuity-selected from remaining roots`).
  Re-run results:
  - `diagnostics/r515b_pseudopure_dome_20260303.csv`: `OK=49`, `FAILED=0`
  - `diagnostics/r515b_pseudopure_vs_honeywell_20260303.png` regenerated
  - `diagnostics/pseudopure_root_candidates_20260303.csv` regenerated
  Continuity checks:
  - 900/950/1000 kPa liquid branch remains continuous after selection:
    `h_l = 264.688, 267.713, 270.642 kJ/kg`
  - maximum adjacent `|Δh_l|` across run: `12.269 kJ/kg` (below `20 kJ/kg` guardrail)
  - jump warnings emitted: `0`.
- 2026-03-03: Layer-by-layer plotting protocol STEP 6 executed (final Honeywell-style composite).
  Code update:
  - `plot_ph_r515b_layers.py` enhanced with `--step 6` composite builder.
  Artifacts:
  - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.png`
  - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.pdf`
  - `diagnostics/plots/ph_layer6_manifest_20260303.json`
  - `diagnostics/plots/ph_layer_status_20260303.json` updated to `completed_layers=[0,1,2,3,4,5,6]`
  STEP-6 visualization settings:
  - pressure axis: log scale
  - plot window: `h=150..500 kJ/kg`, `P=1..100 bar`
  - legend entries:
    `True VLE bubble (liq)`, `True VLE dew (vap)`,
    `Two-phase beta-band`, `Isotherms`,
    `Pseudo-pure overlay (datasheet mode)`.
  Manifest records short git hash `af2ae3d` and exact input CSV dependencies.
- 2026-03-03: Corrected STEP-6 pseudo-pure overlay input chain to use continuity-fixed pseudo-pure data.
  Issue:
  - STEP-6 initially consumed stale `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.csv`
    generated before continuity selector updates, which retained a visible liquid-side jump.
  Action:
  - refreshed `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.csv` and
    `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.png` from
    `diagnostics/r515b_pseudopure_dome_20260303.csv` (continuity-fixed source),
  - regenerated STEP-6 outputs:
    `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.png`,
    `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.pdf`,
    `diagnostics/plots/ph_layer6_manifest_20260303.json`.
  Verification:
  - updated Layer-5 plotting CSV now has `max adjacent |Δh_l| = 12.269 kJ/kg`,
    and `h_l(900,950,1000 kPa) = 264.688, 267.713, 270.642 kJ/kg`.
- 2026-03-03: Isotherm crossover diagnostics protocol — DIAG STEP 1 completed (reproducibility/provenance export).
  Source inputs:
  - `diagnostics/plots/ph_layer4_isotherms_20260303.csv`
  - `diagnostics/plots/ph_layer2_true_vle_boundary_20260303.csv`
  Output files (one per isotherm):
  - `diagnostics/isotherms/T_m19p0C_raw_points_20260303.csv`
  - `diagnostics/isotherms/T_m5p0C_raw_points_20260303.csv`
  - `diagnostics/isotherms/T_10p0C_raw_points_20260303.csv`
  - `diagnostics/isotherms/T_25p0C_raw_points_20260303.csv`
  - `diagnostics/isotherms/T_40p0C_raw_points_20260303.csv`
  - `diagnostics/isotherms/T_55p0C_raw_points_20260303.csv`
  - `diagnostics/isotherms/T_70p0C_raw_points_20260303.csv`
  - `diagnostics/isotherms/T_85p0C_raw_points_20260303.csv`
  Exported fields include:
  - `T_K`, `rho_basis`, `rho_value`, `P_kPa`, `h_kJkg`, `phase_tag_guess`,
    `Z`, `status`, `hb_interp_kJkg`, `hd_interp_kJkg`, `inside_2phase_flag`.
  Counts per isotherm (`n_total`, `n_OK`, `n_flagged_2phase`):
  - `-19 C`: `280`, `280`, `0`
  - `-5 C`: `280`, `280`, `0`
  - `10 C`: `280`, `280`, `2`
  - `25 C`: `280`, `280`, `5`
  - `40 C`: `280`, `280`, `9`
  - `55 C`: `280`, `280`, `14`
  - `70 C`: `280`, `280`, `20`
  - `85 C`: `280`, `280`, `4`
  Two-phase interpolation overlap window used:
  - `P in [5.902676, 23.206615] bar`.
- 2026-03-03: Isotherm crossover diagnostics protocol — STEP D1 completed (raw vs filtered overlays).
  Generated per-isotherm artifacts:
  - `diagnostics/isotherms/plots/T_-19C_raw_vs_filtered_20260303.png`
  - `diagnostics/isotherms/plots/T_-5C_raw_vs_filtered_20260303.png`
  - `diagnostics/isotherms/plots/T_10C_raw_vs_filtered_20260303.png`
  - `diagnostics/isotherms/plots/T_25C_raw_vs_filtered_20260303.png`
  - `diagnostics/isotherms/plots/T_40C_raw_vs_filtered_20260303.png`
  - `diagnostics/isotherms/plots/T_55C_raw_vs_filtered_20260303.png`
  - `diagnostics/isotherms/plots/T_70C_raw_vs_filtered_20260303.png`
  - `diagnostics/isotherms/plots/T_85C_raw_vs_filtered_20260303.png`
  plus flagged-point CSVs:
  - `diagnostics/isotherms/T_<Tc>_points_with_flags_20260303.csv` for each listed temperature.
  Summary file:
  - `diagnostics/isotherms/isotherm_d1_summary_20260303.csv`
  Focus-isotherm metrics (`n_total`, `n_flagged_2phase`):
  - `40 C`: `280`, `9`
  - `55 C`: `280`, `14`
  - `70 C`: `280`, `20`
  Observation from D1:
  - filtering `inside_2phase_flag=True` points alone did **not** remove loop behavior in order-connected polylines,
    indicating ordering/branch-mixing effects remain and D2 boundary-stop alone may be insufficient.
- 2026-03-03: Isotherm crossover diagnostics protocol — STEP D3.1 completed (ordering vs branch-mixing).
  Focus temperatures: `40 C`, `55 C`, `70 C`.
  For each isotherm (filtered to `inside_2phase_flag=False`), generated:
  - A: current-order polyline:
    `diagnostics/isotherms/plots/T_<Tc>_polyline_current_order_20260303.png`
  - B: enthalpy-sorted polyline:
    `diagnostics/isotherms/plots/T_<Tc>_polyline_h_sorted_20260303.png`
  - C: branch-split (`rho_split = median(rho_filtered)`) and each branch h-sorted:
    `diagnostics/isotherms/plots/T_<Tc>_polyline_branch_split_20260303.png`
  Summary CSV:
  - `diagnostics/isotherms/isotherm_d31_summary_20260303.csv`
  Results:
  - `40 C`: `rho_split=29.486730`, branch counts `(vapor=135, liquid=136)`,
    loops: `A=True`, `B=False`, `C=False`.
  - `55 C`: `rho_split=24.387718`, branch counts `(vapor=133, liquid=133)`,
    loops: `A=True`, `B=False`, `C=False`.
  - `70 C`: `rho_split=19.401984`, branch counts `(vapor=130, liquid=130)`,
    loops: `A=True`, `B=False`, `C=False`.
  Interpretation:
  - current isotherm crossover issue is polyline parameterization/order driven.
  - both h-sorting and branch-split+h-sorting remove loops in tested focus temperatures.
- 2026-03-03: Isotherm crossover diagnostics protocol — STEP D4 (70 C only, pre-plot stability filter check) completed.
  Input:
  - `diagnostics/isotherms/T_70C_filtered_OK_only_20260303.csv`
  Method:
  - computed numerical `dP/drho` per point using symmetric perturbation on rho-sorted `(rho, P)` data
  - set `mechanically_stable = (dPdrho > 0)`
  Output:
  - `diagnostics/isotherms/T_70C_filtered_OK_only_with_stability_20260303.csv`
  Counts (before plotting):
  - total points: `260`
  - mechanically stable: `194`
  - mechanically unstable: `66`
  Branch breakdown (`rho_split=1.4142135623730951 kg/m^3`):
  - vapor-like: stable `87`, unstable `9`
  - liquid-like: stable `107`, unstable `57`
- 2026-03-03: Isotherm crossover diagnostics protocol — STEP D4 plotting (70 C stable-only) completed.
  Stable-only branch-split artifacts:
  - `diagnostics/isotherms/T_70C_branch_split_stable_only_20260303.csv`
  - `diagnostics/isotherms/plots/T_70C_branch_split_polyline_stable_only_20260303.png`
  Outcome:
  - vertical spikes: removed
  - remaining loops in vapor-like branch: none
  - remaining loops in liquid-like branch: none
  Stable branch counts:
  - vapor-like `87`
  - liquid-like `107`
- 2026-03-03 — Mechanical Stability Filter Added for Isotherms

- Implemented filter: retain only states where `(∂P/∂ρ)_T > 0`
- Purpose: remove spinodal/unstable states from single-phase isotherm plots
- No changes to EOS or thermodynamic formulation
- References: Callen; Smith, Van Ness & Abbott; Prausnitz et al.
- 2026-03-03: Final isotherm plotting protocol (F1-F5) executed for `T=70 C` with rho-parameterized piecewise branch logic.
  Inputs:
  - `diagnostics/isotherms/T_70C_filtered_OK_only_with_stability_20260303.csv`
  Rules applied:
  - no enthalpy sorting
  - filter: `inside_2phase_flag=False`, `status=OK`, `mechanically_stable=True`
  - fixed split `rho_split=1.4142135623730951 kg/m^3`
  - piecewise segmentation in original `plot_idx` order when `sign(Δh)` flips
  Output:
  - `diagnostics/isotherms/plots/T_70C_final_stable_rho_parametrized_20260303.png`
  Counts:
  - filtered stable points: `194` (`vapor=87`, `liquid=107`)
  - vapor segments: `1`
  - liquid segments: `5`
  Reported checks:
  - vertical spikes gone: `No`
  - vapor branch smooth: `Yes`
  - liquid branch smooth: `Yes`
- 2026-03-03: Implemented block-based stable rho-parameterized isotherm plotting utility (plotting-only, no EOS edits).
  New module:
  - `scripts/plot_isotherms_stable_blocks.py`
  Algorithm implemented:
  - `rho_grid = logspace(1e-3, 1200, 800)`
  - evaluate `P(T,rho,z)`, `h(T,rho,z)`, numerical `dPdrho`
  - strict stability mask `dPdrho > 1e-6`, with automatic retry at `1e-5` when `n_blocks > 2`
  - contiguous stable-block segmentation over rho-grid indices
  - branch choice: lowest-density block as vapor, highest-density block as liquid
  - no enthalpy sorting, no cross-gap stitching, separate branch plotting
  - monotonicity assertions on branch pressure sequences.
  70 C run artifacts:
  - `diagnostics/isotherms/plots/T_70C_final_stable_block_rho_parametrized_20260303.png`
  - `diagnostics/isotherms/stable_block_segmentation_summary_20260303.csv`
  - `diagnostics/isotherms/stable_block_segmentation_summary_20260303.json`
  70 C diagnostics:
  - stable blocks: `3` (tolerance raised to `1e-5` per protocol)
  - vertical segments disappear: `True`
  - monotonicity assertions pass: `True`
- 2026-03-03: Updated stable-block isotherm branch identification to recover true vapor basin and endpoint-consistent labeling (plotting-only).
  File updated:
  - `scripts/plot_isotherms_stable_blocks.py`
  Numerical/labeling changes (no physics edits):
  - soft stability gate: `stable_mask = dPdrho > -tol_soft`, `tol_soft=1e-4` (auto-tighten to `1e-5` when block count > 2)
  - monotonic pressure block-pruning with jitter allowance `min(diff(P_block)) >= -P_tol`, `P_tol=1e-3*max(P_block)`
  - branch labeling by VLE endpoint proximity in scaled `(P,h)` coordinates (dew for vapor high-P endpoint, bubble for liquid low-P endpoint)
  - branch endpoint trimming to nearest VLE endpoint (rho-order, no h-sorting)
  - mandatory per-isotherm diagnostic print line implemented.
  70 C rerun diagnostics (post-trim endpoints):
  - `rho_grid_min=0.001`, `rho_vapor_candidate_min=0.001`
  - stable blocks before/after pruning: `3 / 3`
  - vapor endpoint vs dew: `(1650.81 kPa, 432.351 kJ/kg)` vs `(1616.81 kPa, 418.657 kJ/kg)`, distance `0.138543`
  - liquid endpoint vs bubble: `(3132.18 kPa, 395.255 kJ/kg)` vs `(2032.27 kPa, 396.131 kJ/kg)`, distance `0.541291`
  Acceptance checks:
  - vertical segments disappear: `True`
  - monotonicity assertions pass: `True`
- 2026-03-03: Docstring-only thermodynamic reference update for isotherm stabilization logic.
  Updated file:
  - `scripts/plot_isotherms_stable_blocks.py`
  Scope:
  - documentation only (no numerical or physics behavior changes)
  - added structured docstrings for isotherm sweep, stability filtering, spinodal screening,
    block segmentation, branch labeling, and plotting CLI
  - added explicit thermodynamic basis and mathematical conditions including
    `(∂P/∂ρ)_T` criteria
  - added required literature references (Span 2000; Callen 1985; Prausnitz et al. 1999;
    Lemmon et al. 2018; Bell & Lemmon mixture EOS literature)
  Compliance note:
  - EOS, mixing rules, departure terms, fugacity, and VLE solver logic unchanged.
- 2026-03-03: Reverted the large thermodynamic-reference docstring expansion in
  `scripts/plot_isotherms_stable_blocks.py` per user request.
  Scope:
  - documentation text reduced to concise docstrings only
  - no numerical or plotting behavior changes.

## 2026-03-03 — Layer Memory (P-h Build Protocol Snapshot)

- Layer 0 — Scaffold only
  - Purpose: initialize axes/plot status scaffolding.
  - Artifacts:
    - `diagnostics/plots/ph_layer0_scaffold_20260303.png`
    - `diagnostics/plots/ph_layer_status_20260303.json`

- Layer 1 — Honeywell reference anchor
  - Purpose: reference input sanity anchor from available Honeywell PT data in repo.
  - Artifacts:
    - `diagnostics/plots/ph_layer1_honeywell_reference_20260303.png`
    - `diagnostics/plots/ph_layer1_inputs_20260303.csv`
  - Note: current in-repo Honeywell input used here is PT-only.

- Layer 2 — True VLE boundaries
  - Purpose: plot bubble and dew branch boundaries from true-VLE outputs.
  - Artifacts:
    - `diagnostics/plots/ph_layer2_true_vle_boundary_20260303.png`
    - `diagnostics/plots/ph_layer2_true_vle_boundary_20260303.csv`

- Layer 3 — Two-phase beta band
  - Purpose: fill/interpolate interior two-phase region using beta parameterization.
  - Artifacts:
    - `diagnostics/plots/ph_layer3_true_vle_beta_band_20260303.png`
    - `diagnostics/plots/ph_layer3_beta_grid_20260303.csv`

- Layer 4 — Isotherm context overlays
  - Status: removed from layered pipeline.
  - Reason: scrapped rho-sweep isotherm overlay layer to eliminate
    metastable/two-phase interior tracing artifacts in final layered plot.
  - Artifacts: deleted (`ph_layer4_isotherms_*` removed).

- Layer 5 — Pseudo-pure overlay
  - Purpose: optional near-azeotropic/datasheet-style pseudo-pure overlay.
  - Artifacts:
    - `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.png`
    - `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.csv`
  - Current source alignment: regenerated from continuity-fixed
    `diagnostics/r515b_pseudopure_dome_20260303.csv`.

- Layer 6 — Final composite (Honeywell-style)
  - Purpose: combine Layers 2–5 with presentation settings.
  - Artifacts:
    - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.png`
    - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.pdf`
    - `diagnostics/plots/ph_layer6_manifest_20260303.json`
  - Isotherm overlay excluded (no isotherm legend entry, no Layer-4 dependency).

- Layer 7 — Isentropes (planned, not yet implemented)
  - Purpose: add constant-entropy guide curves for chart interpretation.
  - Status: planned only; no Layer-7 artifacts generated yet in `diagnostics/plots/`.
  - Next expected artifacts (when implemented):
    - `diagnostics/plots/ph_layer7_isentropes_<date>.png`
    - `diagnostics/plots/ph_layer7_isentropes_<date>.csv`

- Layer 8 — Isochors (planned, not yet implemented)
  - Purpose: add constant-density guide curves for chart interpretation.
  - Status: planned only; no Layer-8 artifacts generated yet in `diagnostics/plots/`.
  - Next expected artifacts (when implemented):
    - `diagnostics/plots/ph_layer8_isochors_<date>.png`
    - `diagnostics/plots/ph_layer8_isochors_<date>.csv`

## 2026-03-03 — Isotherm Layer Removal + REFPROP-Style Replacement

- Removed from layered pipeline:
  - isotherm overlay logic in `plot_ph_r515b_layers.py` (Layer-6 no longer reads Layer-4 CSV)
  - isotherm-only helper functions tied to layered overlay
  - old rho-parameterized isotherm helper module: `scripts/plot_isotherms_stable_blocks.py` deleted
- Deleted isotherm-layer artifacts:
  - `diagnostics/plots/ph_layer4_isotherms_20260303.csv`
  - `diagnostics/plots/ph_layer4_isotherms_20260303.png`
  - prior `diagnostics/isotherms/` diagnostic folder removed
- Regenerated Layer-6 outputs (without isotherm legend entries):
  - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.png`
  - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.pdf`
  - `diagnostics/plots/ph_layer6_manifest_20260303.json`
  - `diagnostics/plots/ph_layer_status_20260303.json`

- Added standalone REFPROP-style isotherm algorithm module:
  - `scripts/plot_refprop_style_isotherm.py`
  - Algorithm:
    - single-phase vapor branch from low-density root over `P in [Pmin, Pdew(T)-eps]`
    - single-phase liquid branch from high-density root over `P in [Pbub(T)+eps, Pmax]`
    - straight two-phase connector between `(h_bub, P_bub)` and `(h_dew, P_dew)`
    - no interior two-phase rho-sweep points plotted.
  - References embedded in docstrings:
    - NIST REFPROP docs (v10, 2018)
    - REFPROP-docs metastable notes
    - EES/REFPROP straight-line connector plotting convention
    - Span (2000), Bell et al. (2014)

- Generated 70 C REFPROP-style deliverables:
  - figure: `diagnostics/isotherms_refprop_style/T_70C_REFPROP_style_isotherm.png`
  - csv: `diagnostics/isotherms_refprop_style/T_70C_REFPROP_style_isotherm.csv`
    with columns `[P_bar, h_kJkg, phase_flag]`, `phase_flag ∈ {vapor, two_phase_connector, liquid}`.
- Printed diagnostics:
  - `P_bub(70 C)=20.322731 bar`, `h_bub=396.130635 kJ/kg`
  - `P_dew(70 C)=16.168085 bar`, `h_dew=418.657380 kJ/kg`
  - vapor points: `80`, liquid points: `80`
- 2026-03-03: Debug overlay run executed for dome + REFPROP-style isotherm at `70 C` only.
  Command scope:
  - `scripts/overlay_refprop_isotherms_on_dome.py --temps-c=70 --nv 30 --nl 30`
  Artifacts:
  - `diagnostics/plots/ph_dome_with_refprop_style_isotherms_20260303.png`
  - `diagnostics/plots/ph_dome_with_refprop_style_isotherms_20260303.csv`
  Note:
  - This specific run is single-temperature (70 C) for debugging speed.

## 2026-03-03 — Pseudo-Pure Iso-Diagram Switch (Honeywell-Style)

- Scope honored:
  - no EOS/mixing/departure/fugacity/VLE solver edits
  - plotting/overlay logic only.

- Layer-6 plotting logic switched to pseudo-pure saturation representation:
  - file updated: `plot_ph_r515b_layers.py`
  - removed true-VLE boundary dependence for Layer-6 rendering
  - Layer-6 now uses:
    - pseudo-pure saturation boundaries from `diagnostics/plots/ph_layer5_pseudopure_overlay_20260303.csv`
    - two-phase fill between `h_f(T)` and `h_g(T)` at `P_sat(T)`
    - quality lines `x=0.1..0.9`
  - chart window set to `1..35 bar`.
  - updated Layer-6 legend:
    - `Sat. liquid boundary`
    - `Sat. vapor boundary`
    - `Two-phase region`
    - `Quality lines x=0.1..0.9`

- Old isotherm attempt artifacts removed:
  - `diagnostics/isotherms_refprop_style/` (deleted)
  - `diagnostics/plots/ph_dome_with_refprop_style_isotherms_20260303.png` (deleted)
  - `diagnostics/plots/ph_dome_with_refprop_style_isotherms_20260303.csv` (deleted)

- Added pseudo-pure overlay generator:
  - `scripts/plot_pseudopure_isodiagram.py`
  - outputs:
    - `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`
    - `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm.csv`
  - `phase_flag` values: `{vapor, two_phase_connector, liquid}`.

- 70 C validation prints from pseudo-pure logic:
  - `P_sat(70 C)=16.169297312 bar`
  - `h_f(70 C)=301.269002697 kJ/kg`
  - `h_g(70 C)=418.655032863 kJ/kg`
  - vapor segment `P_max=16.159297311 bar` (target `P_sat-eps=16.159297312`)
  - liquid segment `P_min=16.179297274 bar` (target `P_sat+eps=16.179297312`)
  - connector pressure `=16.169297312 bar` (exact `P_sat`).

- Regenerated Layer-6 pseudo-pure outputs:
  - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.png`
  - `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.pdf`
  - `diagnostics/plots/ph_layer6_manifest_20260303.json`
- 2026-03-03: Implemented continuation-based liquid branch tracing guards in
  `scripts/plot_pseudopure_isodiagram.py` (plotting-only), including:
  - anchor at interpolated `rho_f(T)` from pseudo-pure saturation table
  - first-point solve at `P_sat+eps` with initial/bracket near `rho_f`
  - continuation with `rho_previous` for subsequent liquid pressures
  - hard guards: `rho > rho_f(T)` and `(dP/drho)_T > 0`.
  Observed 70 C diagnostic blocker under strict guards:
  - interpolated `rho_f(70 C) = 978.452499443 kg/m^3`
  - at `P_sat+1 kPa`, scan over `rho in [rho_f, 2000] kg/m^3` found
    `0` sign-change brackets for `P(T,rho)-P_target`, i.e., no valid
    liquid-root candidate satisfying `rho > rho_f`.
  Consequence:
  - strict-guard run produced connector + vapor segment but no accepted liquid
    segment points for 70 C (`phase_flag` counts: vapor 20, connector 2, liquid 0
    in latest debug run with reduced point count).

## 2026-03-03 — Pseudo-Pure Isotherm Continuation Anchors (Root Selection Update)

- File updated (plotting/numerics only): `scripts/plot_pseudopure_isodiagram.py`
- Change implemented:
  - Replaced pseudo-table density anchors with direct TP root anchors at saturation pressure:
    - `rho_f(T) := solve_rho_mass_for_P(T, P_sat(T), phase_hint="liquid")`
    - `rho_g(T) := solve_rho_mass_for_P(T, P_sat(T), phase_hint="vapor")`
  - Kept continuation and monotonic guards:
    - liquid branch monotone increasing in rho with pressure
    - vapor branch monotone decreasing in rho with pressure
    - mechanical-stability gate `(dP/drho)_T > 0`
  - Kept first-point bracket diagnostics for 70 C:
    - liquid: `[rho_f, 1.05*rho_f]` then `rho_hi *= 1.2`
    - vapor: `[0.95*rho_g, rho_g]` then `rho_lo *= 0.8`.

- 70 C run diagnostics after anchor fix (`python3 scripts/plot_pseudopure_isodiagram.py --temps-c=70 --Pmin-bar=1 --Pmax-bar=35 --eps-kpa=1 --nv=80 --nl=80`):
  - `P_sat(70 C) = 1,616,929.731223707 Pa = 16.169297312 bar`
  - `rho_f(70 C) = 463.069935688 kg/m^3`
  - `rho_g(70 C) = 87.204087463 kg/m^3`
  - first liquid bracket at `P_sat + 1 kPa`:
    - `[463.069935688, 472.331334402] kg/m^3`
    - `F(lo) = -0.999996974 kPa`, `F(hi) = 654.324985057 kPa`
    - solved first liquid root: `rho_first_liquid = 463.084906212 kg/m^3`
    - `dP/drho_first_liquid = 66.789800933 kPa/(kg/m^3)`
  - solved first vapor root: `rho_first_vapor = 87.129201411 kg/m^3`
    - `dP/drho_first_vapor = 13.362471418 kPa/(kg/m^3)`
  - monotonicity violations:
    - liquid `0`
    - vapor `0`

- Current known limitation (explicit):
  - Pseudo-saturation table enthalpies (`h_l_kJkg`, `h_v_kJkg`) used for dome drawing are not on the same numerical enthalpy basis as single-phase values returned in the current TP branch path of `pressure_validated_model`.
  - Effect: branch-start `h` continuity checks against table `h_f/h_g` can fail even when pressure-root continuation is correct.
  - This is a plotting-basis consistency issue, not an EOS-equation change.

- Updated artifacts:
  - `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`
  - `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm.csv`
  - `diagnostics/pseudopure_iso/pseudopure_isotherm_overlays_20260303.csv`

## 2026-03-03 — Rollback for 70C Isotherm Debug (User-directed)

- User-directed rollback applied:
  - Reverted recent forced root-anchor/continuation edits in
    `scripts/plot_pseudopure_isodiagram.py`.
  - Intent: debug baseline pseudo-pure isotherm behavior at `T=70 C`.

- 70 C debug run executed (single-temperature only):
  - Command:
    `python3 scripts/plot_pseudopure_isodiagram.py --temps-c=70 --Pmin-bar=1 --Pmax-bar=35 --eps-kpa=1 --nv=80 --nl=80`
  - Outputs regenerated:
    - `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`
    - `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm.csv`

- Printed diagnostics:
  - `P_sat(70C)=16.169297312 bar`
  - `h_f(70C)=301.269002697 kJ/kg`
  - `h_g(70C)=418.655032863 kJ/kg`
  - `vapor_Pmax_bar=16.159297311` (target `P_sat-eps`)
  - `connector_pressure_bar=16.169297312` (exact `P_sat`).

- Current debug observation:
  - CSV phase counts: `vapor=80`, `two_phase_connector=2`, `liquid=0`.
  - This confirms current baseline failure mode is absence of liquid branch
    in the 70 C pseudo-pure isotherm path.

## 2026-03-03 — 70C Isotherm Debug: Removed Liquid Guard

- User-directed change (plotting/root-selection only):
  - Removed liquid-branch rejection guard in `scripts/plot_pseudopure_isodiagram.py`:
    - deleted block that rejected/retried whenever `rho < rho_f_mass`.

- Motivation:
  - At 70 C, pseudo-table anchor `rho_f_mass` (~978.45 kg/m^3) was inconsistent
    with the TP root returned by current pressure solve at `P_sat+1 kPa`
    (~463.08 kg/m^3), causing all liquid points to be skipped.

- Re-run (70 C only) results:
  - Command:
    `python3 scripts/plot_pseudopure_isodiagram.py --temps-c=70 --Pmin-bar=1 --Pmax-bar=35 --eps-kpa=1 --nv=80 --nl=80`
  - Phase counts in CSV:
    - vapor = 80
    - liquid = 80
    - two_phase_connector = 2
  - Artifacts regenerated:
    - `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`
    - `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm.csv`

- Current note:
  - Liquid branch now exists for debugging, but enthalpy-basis consistency against
    pseudo-saturation endpoint values remains an open issue.

## 2026-03-03 — A/B Engine Diagnostic at 70C (Same T,P grid)

- Diagnostic generated to isolate vapor-branch mismatch source:
  - `diagnostics/ab_compare_70C_pvm_vs_mtvle_20260303.csv`
  - comparison setup:
    - pressure grid includes low-P points and near-saturation points
    - solve roots with `pressure_validated_model` (PVM)
    - evaluate the same `(T,rho)` in `mixture_true_vle_copy` (MTVLE)
    - compare `p` and `h` deltas at identical state points.

- Key results:
  - Vapor branch at low pressure (`P=1..35 kPa`) is nearly consistent between engines:
    - `dp_vapor ~ O(1e-5..1e-2) kPa`
    - `dh_vapor ~ O(1e-3..1e-2) kJ/kg`.
  - Near saturation (`P~1616.93 kPa`), vapor state diverges strongly at same rho:
    - `dp_vapor (MTVLE - PVM) ~ -89.5 kPa`
    - `dh_vapor (MTVLE - PVM) ~ -11.76 kJ/kg`.
  - Liquid-state divergence is much larger at same rho:
    - near saturation: `dp_liquid ~ -1981 kPa`, `dh_liquid ~ -57.75 kJ/kg`
    - low pressure (liquid root branch): large nonphysical pressure mismatch confirms branch/path inconsistency.

- Interpretation for current debugging:
  - Missing alignment at isotherm endpoint is not numerical precision.
  - The plotted dome endpoints and isotherm single-phase points are sourced from thermodynamic paths/engines that are not mutually consistent near saturation.

## 2026-03-03 — Pseudo-Pure Isotherm Pressure Window Default Updated

- User-directed plotting-window correction:
  - Updated `scripts/plot_pseudopure_isodiagram.py` defaults:
    - `P_plot_max_bar` in `build_pseudopure_isotherm_segments`: `35.0 -> 100.0`
    - CLI default `--Pmax-bar`: `35.0 -> 100.0`

- Verification run (70 C, defaults):
  - command: `python3 scripts/plot_pseudopure_isodiagram.py --temps-c=70 --nv=80 --nl=80`
  - output CSV check:
    - `liquid_count=80`
    - `liquid_P_max_bar=100.000000187`

- Artifacts regenerated:
  - `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`
  - `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm.csv`

## 2026-03-03 — Root Audit 70C (Bracket vs Fallback)

- Added diagnostics script:
  - `scripts/audit_isotherm_roots_70c.py`
  - Purpose: audit all sign-change intervals and selected roots for 70 C isotherm branches.

- Run outputs:
  - `diagnostics/root_audit_70C_summary_20260303.csv`
  - `diagnostics/root_audit_70C_intervals_20260303.csv`

- Headline findings (70 C):
  - `steps_total = 160` (vapor 80 + liquid 80)
  - `multi_interval_steps = 150`
  - `fallback_steps = 157`
    - `fallback_vapor = 79/80`
    - `fallback_liquid = 78/80`
  - Interval-count distribution:
    - vapor: all `80` steps had `3` sign-change intervals
    - liquid: `70` steps had `3` intervals, `10` steps had `1` interval.

- Interpretation:
  - Current branch tracing is dominated by fallback root selection rather than
    stable local continuation bracket picks.
  - This supports the hypothesis that root-jumping/missing-branch behavior is
    primarily driven by bracket selection and fallback basin choice, not plotting.

## 2026-03-03 — Bracket Enumeration + Continuation Selector (Pseudo-Pure Isotherm)

- Implemented in `scripts/plot_pseudopure_isodiagram.py` (no EOS/reference-state edits):
  - Added full bracket enumeration helper:
    - `enumerate_pressure_root_brackets(...)`
  - Added solve-all-candidates helper:
    - `_solve_candidates_from_brackets(...)`
  - Replaced branch inversion with stateful continuation tracer:
    - `_trace_branch_with_continuation(...)`

- Selection logic now enforces hard gates per candidate:
  - pressure residual `|P-P_target| <= tolP` (default `1e-2 kPa`)
  - mechanical stability `(dP/drho)_T > 0`
  - basin guard via `rho_split = sqrt(rho_f * rho_g)`
    - vapor: `rho < rho_split`
    - liquid: `rho > rho_split`
  - monotonic continuation (upward pressure march): `rho_k >= rho_{k-1}` after first accepted point.

- Fallback policy changed:
  - fallback solver can be invoked only as seed for local bracket retry;
  - fallback rho is never accepted directly.
  - reported metric: `fallback_accept_count` (must remain 0).

- Added adaptive retry:
  - if multiple brackets and no eligible candidate, retry once at half-step pressure.

- Required 70 C instrumentation added:
  - selection trace CSV:
    - `diagnostics/pseudopure_iso/T_70C_selection_trace_20260303.csv`
  - per-step fields include:
    - number of brackets
    - candidate rhos
    - rejection reasons
    - chosen rho/cost
    - fallback seed usage
    - stop reason if branch terminates.

- 70 C run results after selector replacement:
  - output CSV phase counts:
    - vapor = 80
    - liquid = 80
    - connector = 2
  - fallback usage:
    - `fallback_invocations = 0`
    - `fallback_accept_count = 0`
  - branch stop reasons: empty (both branches completed)
  - liquid enthalpy minimum: `266.34 kJ/kg` (>100 kJ/kg).

- Regenerated artifacts:
  - `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm.csv`
  - `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`
  - `diagnostics/pseudopure_iso/T_70C_selection_trace_20260303.csv`

## 2026-03-03 — Dome Cap Added (Visualization)

- Updated `scripts/plot_pseudopure_isodiagram.py` to add an explicit cap line at
  the highest available saturation pressure by connecting the corresponding
  `h_f` and `h_g` endpoints.
- Scope: plotting-only visual closure of dome top; no EOS/solver changes.
- Regenerated artifact:
  - `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`

## 2026-03-03 — Dome Cap Converted to Computed Cap (No Visual Arc)

- Removed the prior visual-only rounded cap in `scripts/plot_pseudopure_isodiagram.py`.
- Added computed cap routine:
  - `_compute_dome_cap_from_saturation(p_bar, h_f, h_g, n_fit=8)`
  - Method: fit top-range `Delta h(P) = h_g - h_f` linearly and solve for
    `Delta h = 0` to estimate closure pressure, then infer cap enthalpy.
- Plot now closes dome top using two computed connector segments from the highest
  available saturation endpoints to the computed cap point.

- Current computed cap (from `diagnostics/r515b_pseudopure_dome_20260303.csv`):
  - `p_cap_bar = 48.9219296004868`
  - `h_cap_kJkg = 426.3612968931469`

- Scope remains plotting-only; EOS and saturation-table generation unchanged.
