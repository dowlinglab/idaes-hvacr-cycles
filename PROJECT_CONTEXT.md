# Project Context

## Project Name
idaes-hvacr-cycles

## Last Updated
2026-03-02

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
- In `compute_table1_properties`, Cp/speed-of-sound/fugacity terms currently depend on finite-difference numerical derivatives rather than analytic second derivatives/composition partials.

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
