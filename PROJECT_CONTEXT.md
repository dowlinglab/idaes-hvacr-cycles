# Project Context

## Project Name
idaes-hvacr-cycles

## Last Updated
2026-03-10

## Current Objective
Establish and maintain a persistent project research memory backbone for disciplined development and reproducible engineering decisions.
Provide a stable query interface for pressure/enthalpy at a single state point while preserving strict model behavior.
Implement a pure-fluid Helmholtz-based saturation solver and P-H dome pipeline for phase-aware validation workflows.
Stabilize an IDAES-based PLR+HX implementation that reproduces the frozen non-IDAES lumped baseline behavior for R134a and R1234ze(E) ambient sweeps.

## System Architecture Overview
- Pending detailed capture.
- This document will track separation of thermodynamic property calls, numerical grid generation, plotting, and validation routines as implementation evolves.
- Added pure-fluid saturation module (`helmholtz_saturation.py`) and wrapper CLI (`scripts/plot_ph_dome.py`) for dome CSV/figure generation.

## Key Technical Decisions
- 2026-03-24: Added a dedicated runner for the plain IDAES PLR-only cycle path (`run_plr_only_cop_vs_ambient.py`) to generate COP-vs-ambient curves without the HX-copy models. Current runner policy uses `PLR=0.75`, `CD=0.13`, fixed evaporator saturation target `-29 C`, condenser approach `9 C`, and ambient sweep `10..45 C` for `R134a` and `R1234ze(E)`.
- 2026-03-24: Executed the new PLR-only runner. `R134a` curve artifacts were generated successfully (`cop_vs_ambient_plr_only_r134a.csv/.png/.pdf`) with `6/8` converged points over `10..45 C`; `10 C` and `15 C` failed with condenser subcooling residuals in the current base IDAES PLR path. `R1234ze(E)` is not yet runnable through the plain `vapor_compression_plr.py` path because that model lacks a working combined IDAES/CoolProp fluid-name mapping for the refrigerant.
- 2026-03-23: For lumped-HX sizing research, pulled an official Carrier evaporator catalogue with model-level capacities and airflows (`Evaporator MT-LT catalogue`, Carrier India) plus Carrier Tenor supermarket condenser family data. Current status: evaporator side has enough published rating information for a first back-calculated `UA` estimate at a catalog point; condenser side still needs a model-level rating sheet or equivalent detailed product data before a comparably specific `UA` back-calculation can be trusted.
- 2026-03-23: For industry-meeting preparation, clarified the recommended IDAES `Heater`-block lumped-HX closure strategy: retain the refrigerant-side cycle topology from `vapor_compression_plr.py`, deactivate SH/SC/approach target constraints for the first lumped build, and replace those closures with evaporator/condenser UA-driven duty equations of the form `Q = epsilon * C_air * DeltaT_drive` using the phase-change-dominant approximation `epsilon = 1 - exp(-UA/C_air)`; compute actual SH/SC only as post-solve diagnostics.
- 2026-03-10: Froze the non-IDAES PLR+HX epsilon-NTU baseline as the reference behavior target for subsequent IDAES replication work.
- 2026-03-10: Active boundary policy for cold-storage runs:
  `T_evap_sat` derived from `T_cold_sp-10` to `T_cold_sp-8` bounds and
  `T_cond_sat` derived from `T_amb+8` to `T_amb+10` bounds.
- 2026-03-10: PLR treatment remains post-correction only (`COP_part = PLF * COP_full`) with normalized reference refrigerant flow in cycle equations.
- 2026-03-10: Repository organization update: active model/runners kept in root; exploratory artifacts moved to `for_review/` and `for_review/` added to `.gitignore`.
- 2026-03-10: IDAES `HeatExchangerNTU`, `HeatExchanger1D`, `HeatExchangerLumpedCapacitance`, and `HeatExchanger` (0-D) copy-based attempts were run against the frozen closure policy but are currently non-convergent on ambient sweeps.
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
- 2026-03-04: For PLR-cycle brainstorming, agreed not to hard-fix evaporator saturation temperature to a single value for cold-storage runs.
- 2026-03-04: Proposed evaporator saturation band relative to cold storage setpoint:
  `T_evap_sat = T_cold_storage - (8 to 10 C)` (example at `-20 C` storage gives `-30 to -28 C`).
- 2026-03-04: Proposed condenser saturation band relative to ambient:
  `T_cond_sat = T_ambient + (8 to 10 C)` instead of fixed approach.
- 2026-03-04: Ambient sweep range for current cold-storage PLR studies updated to `10 to 25 C`.
- 2026-03-04: Confirmed compressor isentropic-efficiency default alignment target across base and PLR code: `0.75`.
- 2026-03-04: Confirmed compressor vapor-only outlet guard is intentionally retained for physical validity and solver robustness.

## Assumptions
- Date format uses ISO (`YYYY-MM-DD`) unless a different format is explicitly required.
- This repository is the authoritative workspace for context tracking.
- Fluid names map to IDAES Helmholtz JSON files available on the configured parameter path.
- During differentiation in `(tau, delta)`, composition is treated as constant.
- Mixture reducing functions are treated as composition-only for this derivative path (no explicit `∂/∂x` coupling terms yet).
- External API composition basis is mass fraction for usability; internal EOS/reducing rules continue to use mole fraction.
- Cold-storage PLR studies use condenser-side ambient coupling; evaporator is coupled via refrigerant-side constraints (no explicit evaporator air-side HX in base model).
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
- Current IDAES PLR+HX copy variants (NTU/1D/LC/0-D) have not yet achieved converged COP-vs-ambient sweeps under the frozen cold-storage boundary policy (latest runs: `0/8` converged points for R134a and R1234ze(E)).

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
- For PLR cold-storage runs, should evaporator/condenser approach bands be implemented as hard bounds, soft penalties, or explicit HX UA-driven equations?
- Which minimal IDAES HX formulation/closure set should be the stabilization target before reintroducing higher-fidelity spatial models?

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
18. For PLR cold-storage modeling, draft and review a non-hard setpoint constraint spec with:
    `T_evap_sat` band, `T_cond_sat` band, and per-point reporting of achieved approaches.
19. Add manual-audit outputs for refrigeration first-law checks:
    `Q_evap`, `W_comp`, `Q_cond`, and residual `Q_cond - (Q_evap + W_comp)`.
20. Stabilize one single-point IDAES PLR+HX solve (`R134a`, `Tamb=20 C`) under frozen bounds before reattempting ambient sweeps.
21. Add explicit scaling coverage for HX `U`, `A`, heat terms, and compressor/valve work in active IDAES PLR+HX copy.
22. Once one-point convergence is obtained, perform warm-start continuation across ambient (`10..45 C`) and then validate `% of Carnot` against the frozen non-IDAES reference trend.

## Change Log (Chronological)
- 2026-03-10: Added copy-ready markdown table for the current 70 C pseudo-pure isotherm values:
  `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm_table_20260310.md`
  with all vapor, connector, and liquid points formatted for slide use.
- 2026-03-10: Extended `plot_ph_r515b_layers.py` STEP-6 renderer with optional single-isotherm overlay support via `--isotherm-csv` and `--isotherm-label`.
- 2026-03-10: Regenerated the pseudo-pure Honeywell-style p-H chart with the current 70 C pseudo-pure isotherm overlaid:
  `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.png`,
  `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.pdf`,
  `diagnostics/plots/ph_layer6_manifest_20260303.json`.
- 2026-03-10: Re-generated the previously approved pseudo-pure Honeywell-style p-H chart figure from the existing Layer-5 pseudo-pure overlay artifact using `plot_ph_r515b_layers.py --step 6 --stamp 20260303`.
- 2026-03-10: Refreshed chart-style artifacts with pseudo-pure saturation boundaries, two-phase fill, and quality lines `x=0.1..0.9`:
  `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.png`,
  `diagnostics/plots/ph_layer6_final_honeywell_style_20260303.pdf`,
  `diagnostics/plots/ph_layer6_manifest_20260303.json`.
- 2026-03-10: Updated `scripts/plot_pseudopure_isodiagram.py` to user-directed debug mode that removes all dome plotting elements (fill, sat boundaries, quality lines, computed cap) and renders only the requested isotherm overlay segments plus connector.
- 2026-03-10: Regenerated overlay-only artifacts with the project conda interpreter:
  `diagnostics/pseudopure_iso/ph_dome_pseudopure_iso_overlays.png`,
  `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm.csv`,
  `diagnostics/pseudopure_iso/T_70C_selection_trace_20260310.csv`.
- 2026-03-10: Added root-level frozen non-IDAES PLR+HX alias `vapor_compression_plr_hx.py` from epsilon-NTU copy path.
- 2026-03-10: Added `for_review/` directory for non-active scripts/artifacts and updated `.gitignore` to ignore `for_review/`.
- 2026-03-10: Added and tested copy-based IDAES PLR+HX variants:
  `vapor_compression_plr_hx_idaes_ntu.py`, `vapor_compression_plr_hx1d.py`,
  `vapor_compression_plr_hx_lc.py`, and `vapor_compression_plr_hx_0d.py`
  with paired R134a/R1234ze runners and SH/SC CSV reporting.
- 2026-03-10: Current convergence status for latest ambient sweeps (`10..45 C`, frozen bounds):
  all IDAES copy variants above report `0/8` converged points for both R134a and R1234ze(E).
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
- 2026-03-04: Added PLR/cold-storage brainstorming notes and manual-check equations in `verification/plr_brainstorm_manual_checks_2026-03-04.md`.
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

## 2026-03-09 — PLR Model Breadcrumbs (`vapor_compression_plr.py`)

- Scope and intent:
  - `vapor_compression_plr.py` is the active IDAES PLR cycle implementation for
    full-load COP solve plus PLF-based post-correction.
  - Current soft-constraint direction is to enforce evaporator/condenser
    approach bands with nonnegative slack variables and objective penalty.

- PLR/CD and COP handling:
  - `self.model.fs.cop` remains the solved full-load COP variable.
  - PLF is computed by helper: `PLF = 1 - CD * (1 - PLR)`.
  - Part-load COP is computed as post-processing:
    `COP_part = PLF * COP_full`.
  - `set_specifications(...)` includes `plr` and `cd` optional overrides and
    updates model parameters through `.set_value(...)`.

- Soft approach constraints (target band 8-10 K):
  - Parameters: `approach_min`, `approach_max`, `soft_approach_weight`.
  - Slacks: `eps_evap_low`, `eps_evap_high`, `eps_cond_low`, `eps_cond_high`.
  - Temperature references:
    - `cold_storage_T` (K)
    - `ambient_T_soft` (K)
  - Equations represent:
    - evaporator approach: `cold_storage_T - T_sat_evap in [8, 10]` (soft)
    - condenser approach: `T_sat_cond - ambient_T_soft in [8, 10]` (soft)
  - Activation is controlled in `set_specifications(...)` via
    `use_soft_approach`.

- Known implementation corrections already made:
  - Fixed PLR parameter initialization typo (`self.plr` instead of undefined
    `self.plf`).
  - Corrected soft ambient param assignment to `self.model.fs...`.
  - Corrected condenser soft constraint references from invalid
    `approach.min/max` to `approach_min/max`.
  - Added getters:
    - `get_full_load_cop()`
    - `get_part_load_cop()`

- Documentation updates:
  - Module header updated with author/codex/QA breadcrumb template.
  - Detailed method docstrings added with governing math and intent.

## Session Breadcrumb (2026-03-09)

### Scope
- Focus shifted to PLR-vs-base overlay reproducibility for cold-storage ambient sweeps.
- User requested strict separation: keep `vapor_compression.py` unchanged and apply PLR in a copied model only.

### Files Added/Updated
- Added: `vapor_compression_plr_only.py` (copy of `vapor_compression.py` + PLR post-correction only).
- Updated: `vapor_compression_plr.py` overwritten with `vapor_compression_plr_only.py` contents per user request.
- Updated run scripts:
  - `run_plr_cold_storage_r134a_copy.py`
  - `run_plr_cold_storage_r1234zee_copy.py`

### PLR Logic (Current)
- Implemented as post-correction only:
  - `PLF = 1 - CD * (1 - PLR)`
  - `COP_part = PLF * COP_full`
- Added guards:
  - `PLR in [0, 1]`
  - `CD in [0, 1]`
- No soft/hard HX approach constraints in the active PLR model path.

### Active Overlay Setup (Current)
- Overlays now include 3 curves on same ambient grid:
  - `Carnot COP`
  - `vapor_compression COP`
  - `vapor_compression_plr COP`
- Ambient grid: `10..45 C` in `5 C` steps.
- Ambient coupling active in both base and PLR runs via hard condenser relation:
  - `T_cond_sat = T_amb + condenser_approach`
- Current `condenser_approach`: `20 C`.
- Current evap outlet temperature bounds passed in runs:
  - `(-55 C, -20 C)` (interpreted as `(Tsp-35, Tsp)` with `Tsp=-20 C`).
- Pressure bounds used in overlays:
  - R134a: low `(60,120) kPa`, high `(500,1000) kPa`
  - R1234ze: low `(40,80) kPa`, high `(350,800) kPa`

### Latest Run Status
- With current settings above:
  - R134a overlay convergence: `8/8` for base and PLR
  - R1234ze overlay convergence: `8/8` for base and PLR
- Figure outputs:
  - `cop_vs_ambient_plr_cold_storage_r134a_copy.png/.pdf/.csv`
  - `cop_vs_ambient_plr_cold_storage_r1234zee_copy.png/.pdf/.csv`

### Notes For Next Session
- User asked to potentially delete `vapor_compression_plr_only.py` after explicit confirmation.
- If overlays look unexpectedly flat, first check whether condenser ambient coupling (`condenser_approach`) is being passed/activated.

## Session Breadcrumb (2026-03-09, Late)

### Current Agreed Direction
- Brainstorm mode only for next phase (no HX code edits yet).
- User-selected future HX method: epsilon-NTU with 3 zones.
- Must remain in IDAES model ecosystem.

### Current Codebase State
- `vapor_compression.py` remains unchanged.
- `vapor_compression_plr_only.py` was deleted.
- Kept only `vapor_compression_plr.py` as PLR model copy.
- Removed other PLR variant files (`_SN`, `_ambient_bounds`, `_condsoft`).

### Active Overlay Scripts and Outputs
- Scripts:
  - `run_plr_cold_storage_r134a_copy.py`
  - `run_plr_cold_storage_r1234zee_copy.py`
- Current overlays include 3 curves on same grid:
  - Carnot COP
  - vapor_compression COP
  - vapor_compression_plr COP
- Current run settings in scripts:
  - ambient grid: `10..45 C` step `5 C`
  - condenser approach coupling active: `T_cond_sat = T_amb + 20 C`
  - evaporator temperature bounds: `(-55 C, -20 C)`
  - PLR settings: `PLR=0.75`, `CD=0.13`
- Latest rerun status: both fluids `8/8` converged for base and PLR.

### PLR Logic Status
- PLR applied as post-correction only in `vapor_compression_plr.py`:
  - `PLF = 1 - CD*(1-PLR)`
  - `COP_part = PLF * COP_full`
- Guards currently enforced:
  - `PLR in [0,1]`
  - `CD in [0,1]`

### IDAES HX Note (Confirmed)
- Current evap/condenser are modeled as `Heater` units, not explicit HX units.
- Confirmed IDAES epsilon-NTU unit model to use in next phase:
  - `HeatExchangerNTU` (`idaes.models.unit_models`).

### User Preference Reminders
- Do not make unsolicited suggestions.
- Do not modify code unless explicitly requested.

## Session Breadcrumb (2026-03-09, Overlay-4 on Copy)

### Change Scope
- Kept original files untouched:
  - `vapor_compression.py`
  - `vapor_compression_plr.py`
  - existing `run_plr_cold_storage_*_copy.py`
- Added a PLR copy module:
  - `vapor_compression_plr_copy.py` (direct copy of `_plr` at this step)
- Added new overlay runner copies:
  - `run_plr_cold_storage_r134a_plr_copy_overlay.py`
  - `run_plr_cold_storage_r1234zee_plr_copy_overlay.py`

### Overlay Curves Requested
- New overlay runners produce 4 curves:
  - Carnot COP
  - vapor_compression COP
  - PLR only COP (`PLR * COP_full`)
  - PLR + COP (`PLF * COP_full`, where `PLF = 1 - CD*(1-PLR)`)

### Output Filenames
- R134a:
  - `cop_vs_ambient_plr_cold_storage_r134a_plr_copy_overlay.csv`
  - `cop_vs_ambient_plr_cold_storage_r134a_plr_copy_overlay.png`
  - `cop_vs_ambient_plr_cold_storage_r134a_plr_copy_overlay.pdf`
- R1234ze:
  - `cop_vs_ambient_plr_cold_storage_r1234zee_plr_copy_overlay.csv`
  - `cop_vs_ambient_plr_cold_storage_r1234zee_plr_copy_overlay.png`
  - `cop_vs_ambient_plr_cold_storage_r1234zee_plr_copy_overlay.pdf`

### Execution Result (Overlay-4 on Copy)
- Executed:
  - `python3 run_plr_cold_storage_r134a_plr_copy_overlay.py`
  - `python3 run_plr_cold_storage_r1234zee_plr_copy_overlay.py`
- Convergence:
  - R134a: vapor_compression `8/8`, PLR copy `8/8`
  - R1234ze: vapor_compression `8/8`, PLR copy `8/8`
- Figures generated:
  - `cop_vs_ambient_plr_cold_storage_r134a_plr_copy_overlay.png`
  - `cop_vs_ambient_plr_cold_storage_r1234zee_plr_copy_overlay.png`

## Session Breadcrumb (2026-03-09, NTU Copy Start)

### User Clarification Captured
- "PLR + COP" means: keep current cycle + PLR structure, replace heater coils with `HeatExchangerNTU`, start single NTU per coil, keep PLR as post-correction, calibrate UA to baseline, then move to zoned NTU.
- User requested counterflow selection.

### Copy-Only Implementation Added
- New file created (no edits to originals):
  - `vapor_compression_plr_hx_ntu_copy.py`
- Model structure in copy:
  - Evaporator: `HeatExchangerNTU` (air hot side, refrigerant cold side)
  - Condenser: `HeatExchangerNTU` (refrigerant hot side, air cold side)
  - Compressor + expansion valve retained
  - PLR retained as post-correction (`COP_part = PLF * COP_full`)

### Counterflow Status
- Counterflow epsilon-NTU relation implemented explicitly for both HX blocks:
  - `epsilon = (1-exp(-NTU*(1-Cr))) / (1-Cr*exp(-NTU*(1-Cr))+eps_reg)`

### Current Solve Status
- Initial NTU copy was overconstrained; corrected to DOF = 0.
- Current blocker: one-point smoke solve still fails to converge robustly (Ipopt max-iterations / bad status in this first pass).
- No changes were made to `vapor_compression.py` or `vapor_compression_plr.py`.

## Session Breadcrumb (2026-03-09, Overlay Refresh)
- Re-ran 4-curve overlay scripts:
  - `run_plr_cold_storage_r134a_plr_copy_overlay.py`
  - `run_plr_cold_storage_r1234zee_plr_copy_overlay.py`
- Curves included in each figure:
  - Carnot COP
  - vapor_compression COP
  - PLR only COP
  - PLR + COP
- Convergence: both fluids 8/8 for vapor_compression and PLR copy runs.
- Refreshed outputs:
  - `cop_vs_ambient_plr_cold_storage_r134a_plr_copy_overlay.png/.pdf/.csv`
  - `cop_vs_ambient_plr_cold_storage_r1234zee_plr_copy_overlay.png/.pdf/.csv`

## Session Breadcrumb (2026-03-09, NTU Stabilization Attempt)
- User requested overlays with: Carnot, vapor_compression, PLR only (cycle degradation), and NTU.
- Attempted NTU stabilization in `vapor_compression_plr_hx_ntu_copy.py`:
  - Replaced effectiveness relation with phase-change dominant form: `epsilon = 1 - exp(-NTU)`.
  - Deactivated fragile constraints for this copy solve path:
    - superheating_constraint
    - subcooling_constraint
    - vapor_constraint
  - Increased Ipopt iteration limit and made solve failure handling non-fatal.
- Status:
  - DOF confirmed 0 after specification.
  - Still encountering non-convergence (`infeasible` and `maxIterations`) on smoke solve.
  - Stable overlays remain available from non-NTU copy scripts.

## Session Breadcrumb (2026-03-09, Zoned NTU Copy)
- Added new copy model:
  - `vapor_compression_plr_hx_ntu_zoned_copy.py`
- Implemented zones:
  - Evaporator: `evap_tp`, `evap_sh`
  - Condenser: `cond_ds`, `cond_tp`, `cond_sc`
- Preserved PLR post-correction (`get_part_load_cop`).
- Added per-zone NTU getter (`get_zone_ntu`).
- Smoke-test status (R134a, ambient=20 C case):
  - DOF corrected from -1 to 0.
  - Still not converging (`infeasible` / `maxIterations`).
  - Current reported NTU values are available but solve is not successful.

## Session Breadcrumb (2026-03-09, Air-Side Boundary Update)
- Updated NTU copy models to enforce user-requested air-side conditions:
  - Evaporator air inlet fixed at `-20 C`
  - Condenser air inlet set from ambient sweep input (`10..45 C` in runners)
- Files patched:
  - `vapor_compression_plr_hx_ntu_copy.py`
  - `vapor_compression_plr_hx_ntu_zoned_copy.py`
- Verification checks:
  - Single-NTU copy: evap air `-20 C`, condenser air tracks ambient (tested at `15 C`)
  - Zoned NTU copy: both evap zones air `-20 C`, all condenser zones air track ambient (tested at `45 C`)

## Session Breadcrumb (2026-03-09, Convergence Debug Progress)
- Debug status for zoned NTU copy (`vapor_compression_plr_hx_ntu_zoned_copy.py`):
  - Direct solve at target bounds still unstable.
  - Relaxed -> tightened continuation now works through bound tightening to target at fixed ambient 20 C.
- Ambient sweep test with continuation (R134a, Tamb 10..45 C):
  - Converged with positive COP at 10, 15, 20 C.
  - 25, 30 C: non-converged (maxIterations).
  - 35, 40, 45 C: solver returned but with negative COP (physically invalid cooling point).
- Interpretation:
  - Model can find mathematically feasible states at higher ambient that do not represent cooling operation.

### Additional Debug Result (Cooling-Mode Constraint)
- Added `Q_evap_total >= 0` constraint in zoned NTU copy to prevent heating-mode solutions.
- Continuation sweep outcome (R134a, Tamb 10..45 C, current bounds):
  - Converged + positive NTU COP at 10, 15, 20 C.
  - Non-converged or infeasible for 25 C and above.

### Evaporator Bound Test (Tsp-50 to Tsp)
- Tested NTU continuation with evaporator bounds expanded to `(-70 C, -20 C)`.
- Result unchanged versus prior run:
  - Converged/positive at Tamb: 10, 15, 20 C
  - Failed at Tamb >= 25 C (maxIterations/infeasible)
- Conclusion: high-ambient non-convergence is not resolved by widening evaporator temperature bounds alone.

## Session Breadcrumb (2026-03-09, Evap-vs-Setpoint Coding)
- Updated active overlay scripts to compute evaporator bounds from cold storage setpoint:
  - `evaporator_temperature = (Tsp-35, Tsp)` from `cold_storage_setpoint_c`.
  - Files:
    - `run_plr_cold_storage_r134a_plr_copy_overlay.py`
    - `run_plr_cold_storage_r1234zee_plr_copy_overlay.py`
- Updated NTU copy model APIs to support setpoint-based evaporator bounds directly:
  - New args in `set_specifications`:
    - `cold_storage_setpoint`
    - `evap_offset_bounds`
  - Files:
    - `vapor_compression_plr_hx_ntu_copy.py`
    - `vapor_compression_plr_hx_ntu_zoned_copy.py`
- Verified mapping example:
  - `cold_storage_setpoint=-20`, `evap_offset_bounds=(-50,0)` => bounds `(-70, -20) C`.

## Session Breadcrumb (2026-03-09, HX1D Copy Created)
- Added new PLR copy with 1D heat exchangers:
  - `vapor_compression_plr_hx1d_copy.py`
- Core model changes on copy only:
  - Replaced evaporator and condenser with `HeatExchanger1D`
  - Countercurrent flow pattern
  - Explicit finite elements (constructor arg, default 8)
  - Fixed geometry and distributed `heat_transfer_coefficient[t,x]`
  - Refrigerant property package retained as Helmholtz two-phase-capable package
  - `EnergyBalanceType.enthalpyTotal` on both HX sides
  - PLR remains post-correction (`COP_part = PLF * COP_full`)
- Boundary conventions in copy:
  - Evaporator air inlet fixed to cold-storage setpoint (default `-20 C`)
  - Condenser air inlet follows ambient input
- Smoke build/spec check:
  - Model constructs and accepts specs, with finite elements applied
  - Reported `DOF = -1` in current quick check (not yet tuned to converged run)

## Session Breadcrumb (2026-03-09, HX1D Datasheet Anchors Updated)
- User-specified equipment references are now the explicit anchors for HX1D copy defaults:
  - Condenser reference:
    - Bitzer/Buffalo Trident `HX-518-1` (`HVB-2-4R-4P`)
    - Source: `https://www.bitzer.de/shared_media/documentation/hx-518-1-au.pdf`
  - Evaporator reference:
    - Norlake Split-Pak sheet (`NASJ125RL4` / `WL6A094SDAS`)
    - Source: `https://norlake.com/wp-content/uploads/2020/07/nasj125rl4.pdf`
- Active default anchors in `vapor_compression_plr_hx1d_copy.py`:
  - `finite_elements=4`
  - Evaporator: `area=6.0 m^2`, `length=1.11 m`, `U=80 W/m^2-K`, `air_flow=31 mol/s`
  - Condenser: `area=100.0 m^2`, `length=1.914 m`, `U=90 W/m^2-K`, `air_flow=356 mol/s`
- Notes:
  - Values are first-pass engineering anchors for initialization/stability and can be refined
    after steady-state convergence is robust across ambient sweep points.

## Session Breadcrumb (2026-03-09, HX1D Explicit Scaling Applied)
- Applied explicit HX1D scaling in:
  - `vapor_compression_plr_hx1d_copy.py`
- Added `_apply_recommended_scaling()` and called it:
  - after flowsheet/constraint build
  - in `set_specifications()` before `calculate_scaling_factors`
- Scaling targets implemented:
  - Pressure vars/constraints: `1e-6`
  - Enthalpy vars: `1e-5`
  - Temperature vars/constraints: `1e-2`
  - Heat/work vars/constraints: `1e-4`
  - Mass flow vars: `1`
  - COP var/constraint: `1`
- Spot-check results after scaling (R134a, high-side `(500, 2000)` kPa):
  - `Tamb=15 C`: non-converged (`infeasible`)
  - `Tamb=45 C`: non-converged (`infeasible`)
- Interpretation:
  - Scaling is active and reduces generic scaling warnings, but convergence still requires
    manifold relaxation/continuation around pressure and approach constraints.

## Session Breadcrumb (2026-03-09, New IDAES Lumped-Capacitance Copy)
- Added new copy model:
  - `vapor_compression_plr_hx_lc_copy.py`
- Implemented from PLR structure with coil-only swap:
  - evaporator -> `HeatExchangerLumpedCapacitance`
  - condenser -> `HeatExchangerLumpedCapacitance`
  - compressor/valve/PLR post-correction kept aligned with PLR logic.
- Added LC/HX config parameters:
  - `UA_evap_hot`, `UA_evap_cold`, `UA_cond_hot`, `UA_cond_cold`
  - `wall_heat_capacity_evap`, `wall_heat_capacity_cond`
  - wall temperature initialization params
- Applied resistance-addition style default UA split:
  - evap: hot=640 W/K, cold=1920 W/K
  - cond: hot=36000 W/K, cold=12000 W/K
- Good-news checkpoint:
  - Model now constructs successfully with `HeatExchangerLumpedCapacitance`.
  - `flow_pattern` integration fixed.
  - DOF after specifications is now `0` (closed model).
- Current status:
  - Ambient smoke sweep (`10,15,25,35,45 C`) still non-converged.
  - Repeated dominant residuals are at condenser LC equations:
    - `delta_temperature_in_equation`
    - `heat_transfer_equation`
    - `wall_temperature_eq`
    - `condenser.cold_side.material_balances[CO2]`

## 2026-03-09: IDAES 0D Condenser-3-Zone Copy (Debug State)

- Added copy module: `vapor_compression_plr_hx_0d_cond3.py`.
- Implemented condenser as three serial `HeatExchanger` (0-D) blocks:
  - `desuperheater` -> `condenser` (Underwood delta-T callback) -> `subcooler`.
- Refrigerant property setup switched to Helmholtz `PhaseType.LG` and `StateVars.PH`.
- Added copy runners:
  - `run_plr_cold_storage_r134a_hx_0d_cond3.py`
  - `run_plr_cold_storage_r1234zee_hx_0d_cond3.py`
- Generated zoned non-IDAES debug references on ambient 10..45 C:
  - `debug_zoned_ref_r134a_10_45.csv`
  - `debug_zoned_ref_r1234zee_10_45.csv`
- Generated pointwise comparison files:
  - `debug_compare_zoned_vs_idaes_cond3_r134a.csv`
  - `debug_compare_zoned_vs_idaes_cond3_r1234zee.csv`
- Current result: IDAES cond3 runs converge at `0/8` points for both fluids; zoned non-IDAES references are converged at all points.

## Session Breadcrumb (2026-03-10, Zoned Warm-Start / Air-Chain Debug)

### Active File Under Debug
- `/Users/snarasi2/idaes-hvacr-cycles/vapor_compression_plr_hx_0d_cond3.py`

### User-Requested Immediate Focus
- Pause broad sweeps and debug one ambient point first.
- Print PFD + stream table for that single point to verify continuity.
- Stabilize warm-start bridge from non-IDAES zoned model without nonphysical air temperatures.

### Current Observed Failure Signature
- Solver remains non-converged at the one-point debug run.
- Air-side seeded states can jump to nonphysical values (~500 K / 226.85 C) before solve.
- Residual pattern indicates inconsistent initialization between manual inlet/outlet fixes and arc-coupled stream equalities.

### Root-Cause Hypothesis (Current)
1. Air-side over-specification during initialization:
   - A downstream HX inlet can be seeded/fixed inconsistently while also arc-linked.
2. Warm-start bridge inconsistency:
   - Zoned-duty-derived air temperature reconstruction (`T += Q/C_air`) can produce unrealistic seeds on mismatched basis.
3. Refrigerant enthalpy continuity conflict:
   - Hardcoded inlet seeds in `initialize()` can disagree with zoned-profile handoff values.

### Next Debug Actions (Single-Point Only)
1. Apply strict "one fixed condenser-air inlet" policy (all downstream air inlet temperatures unfixed).
2. Keep serial condenser-air arcs physically consistent and use ambient-range temperature seeds only.
3. Enforce refrigerant seed continuity at DS->COND->SC interfaces.
4. Remove/disable hardcoded condenser/subcooler inlet enthalpy anchors when zoned warm-start is active.
5. Re-run one-point case and report only PFD + stream table before any ambient sweep.

### Guardrail Reminder
- Do not alter `vapor_compression.py` baseline.
- Keep all HX experiments isolated in copy files.

## Session Breadcrumb (2026-03-10, End-of-Day Cond3 Re-Chain Debug)

### Active Debug File
- `/Users/snarasi2/idaes-hvacr-cycles/vapor_compression_plr_hx_0d_cond3.py`

### User-Requested Structural Patches Applied
1. Re-chained condenser-train air path to:
   - `Ambient -> SC_cold_in -> COND_cold_in -> DS_cold_in -> Exhaust`
2. Enforced single fixed condenser-train air inlet:
   - only `subcooler.cold_side_inlet` (`T` and `P`) fixed to ambient/1 atm.
3. Added high-pressure manifold constraints:
   - `P_comp_out = P_DS_in = P_COND_in = P_SC_in` (explicit equalities).
4. Added subcooler refrigerant-side zero-dP closure:
   - `P_SC_hot_out = P_SC_hot_in`.
5. Added air-side temperature caps:
   - `T_air_out <= 350 K` on evaporator hot-side out and condenser-train cold-side outs.
6. Added subcooler outlet enthalpy initialization kickstart in fallback seed:
   - `h_SC_out_seed = h_f - 5000 J/kg`.

### Additional Fixes Applied During Debug
- Prevented warm-start routine from reactivating arc `pressure_equality` constraints (to avoid pressure over-closure with explicit manifold constraints).
- Added deterministic fallback warm-start seed when zoned bridge fails, using target `P_low/P_high` and staged refrigerant enthalpy seeds.

### Latest One-Point Test Case
- Fluid: `R134a`
- Ambient: `20 C`
- Cold storage setpoint: `-20 C`
- Bounds used:
  - `P_low: 60..200 kPa`
  - `P_high: 500..4000 kPa`
  - `T_evap_sat band: [Tsp-10, Tsp-8] C`
  - `T_cond_sat band: [Tamb+8, Tamb+10] C`

### Latest Solver Status
- Final one-point solve result: **not converged**.
- Reported termination in solve path: `other` with message `Too few degrees of freedom (rethrown)` during initialization solve attempts.
- Post-patch model-level closure check before solve can report `DOF=0`, but initialization-phase solve still hits the above termination.

### Last Printed State Snapshot (Non-Converged)
- Refrigerant states remained on physically plausible pressure levels for S1..S4 except SC outlet pressure anomaly persists in snapshot (`S3a` high).
- Air-side outlet temperatures in snapshot still pinned at `226.85 C` (nonphysical), indicating the current initialization path is not reaching a valid physical branch.

### Next-Step Reminder for Restart
- First action tomorrow should be to diagnose initialization-specific DOF/over-closure separately from steady-state DOF:
  - inspect which constraints are active inside `initialize()` solve calls,
  - avoid solving full model while temporary fixed anchors are still active,
  - then rerun one-point table/PFD before any ambient sweep.

## Session Breadcrumb (2026-03-10, Late-Night Residual Diagnosis)

### Requested Focus
- Explain persistent nonphysical air temperatures and SC pressure jump in `vapor_compression_plr_hx_0d_cond3.py` after closure patches.

### Structural Status (Current)
- `subcooler_hot_dp0` is active.
- Built-in `subcooler.hot_side.pressure_balance[0]` is deactivated to avoid duplicate zero-dP equations.
- Condenser-zone areas are fixed; `cond_area_sum` deactivated to avoid redundancy.
- Current model-level DOF at set-spec is `0`.
- IDAES structural diagnostics no longer show structural singularity; remaining warning is potential evaluation errors in HX LMTD equations.

### One-Point Debug Outcome (R134a, Tamb=20 C)
- Solve status: not converged (`solver_exception`).
- Key residuals at failed iterate:
  - `subcooler_hot_dp0` residual: `+1.281e6 Pa`
  - `P_high_comp_out` residual: `-4.517e5 Pa`
  - `subcooler.heat_transfer_equation[0]` residual: `-175.87`
  - `subcooler.delta_temperature_out_equation[0]` residual: `+21.05`
- SC duty at failed iterate:
  - `Q_hot = 0 W`, `Q_cold = 0 W` (effectively no SC heat transfer occurring in failed state)

### Interpretation Captured
- The displayed air temperatures (e.g., `226.85 C`) are not physical predictions.
- They are unconverged iterate values after NLP failure; constraints tying pressure and HX delta-T are not satisfied.
- Pressure jump `SC_in -> SC_out` is also an unconverged artifact, not a solved thermodynamic result.

### Next Debug Target
- Focus on numerical robustness of HX equations (delta-T/LMTD behavior and initialization path), not additional structural closure edits.

## Session Breadcrumb (2026-03-10, SC Balance Snapshot)

### User Request
- Provide explicit mass and energy balances around Subcooler (SC) hot and cold sides using latest one-point run.

### One-Point Context
- Fluid: `R134a`
- Ambient: `20 C`
- Solver status: non-converged (`solver_exception`)

### SC State Values at Failed Iterate
- Hot side (refrigerant):
  - `m_dot_hot = 1.0 kg/s`
  - `h_hot_in = 241722.392 J/kg`
  - `h_hot_out = 236722.392 J/kg`
  - `Q_hot = 0.0 W`
- Cold side (air/flue gas package):
  - `N_dot_cold = 356.0 mol/s`
  - `h_cold_in = -250.0 J/mol`
  - `h_cold_out = -250.0 J/mol`
  - `Q_cold = 0.0 W`
  - `T_cold_in = 20.0 C`
  - `T_cold_out = 226.85 C`

### Balance Interpretation
- Hot-side mass continuity holds (`m_in = m_out`).
- Hot-side energy balance is not satisfied at this failed iterate:
  `m_dot*(h_out-h_in)+Q_hot = -5000 W`.
- Cold-side energy equation is algebraically zero with current values, but the paired temperature/enthalpy state is nonphysical due to failed convergence.
- Cross-side duty closure (`Q_hot + Q_cold = 0`) is trivially zero here, but not representative of a physical solved SC duty.

### Practical Note
- These are failed-iterate diagnostics, not converged performance outputs.

## Session Breadcrumb (2026-08-12, `mixture_model_one_point.py` Derivative Audit + Bell (2023) Citation Check)

### Scope
- Read-only walkthrough of `mixture_model_one_point.py` (no edits made to that
  file this session; user is writing/editing its functions directly).
- User is deciding whether to keep central finite differences for residual
  second derivatives or switch to IDAES's compiled exact-derivative external
  functions; this session gathered the facts needed for that decision.

### Confirmed: Which Derivatives Are Analytic vs Finite-Difference (current file)
- First derivatives of alpha^r (pure-fluid `alphar_idaes_with_derivs` and
  departure `bell2023_departure_alphar`), and their mixture chain-rule
  combination (`ar_tau_mix`, `ar_del_mix`): analytic, no FD.
- Second derivatives of the mixture alpha^r (`ar_tautau`, `ar_deldel`,
  `ar_taudel`, feeding `cv`/`cp`/speed of sound in `compute_table1_properties`):
  finite difference, via `_fd_2d` wrapping a local `ar_func(tau,delta)` closure.
- Ideal term's second tau-derivative (`a01_tautau_i`/`a02_tautau_i`, feeding
  `cv0`): also finite difference, via `_a0_tautau_i`.
- Composition derivative for fugacity (`dar_dx1` via
  `_bell2023_reducing_derivs_binary`): analytic in this file — note this
  differs from `linear_model_codex.py`, where `PROJECT_CONTEXT.md` history
  above records that the analytic-composition-derivative migration was rolled
  back per user instruction and that file still uses finite-difference
  fugacity composition derivatives. The two mixture modules currently disagree
  on this point.
- Confirmed `_fd_2d`'s three formulas (`f_tautau`, `f_deldel`, `f_taudel`) are
  all standard second-order **central** (symmetric ±h) finite differences, not
  one-sided; user has decided to keep central finite differences rather than
  switch to IDAES's compiled exact-derivative external functions for now.

### IDAES External Function Architecture (researched, not implemented)
- Confirmed via IDAES docs + a maintainer's GitHub Discussion #428: IDAES's
  general Helmholtz EOS implements the *entire* EOS (ideal + residual terms,
  all derived properties) inside a compiled C++ shared library per fluid
  (e.g. `iapws95_external.so`, `swco2_external.so` style naming), wired into
  Pyomo via `ExternalFunction`, explicitly to give the NLP solver exact
  analytic first *and* second derivatives (quoted from the parameter-file
  docs: "External functions allow for exact first and second derivatives to
  be calculated and used by the solver").
- Could NOT confirm the literal internal C++ symbol names: `github.com`
  directory/API browsing was blocked (`robots.txt` disallow, then a 403 from
  the REST API), and neither readthedocs page enumerates them. Did not
  attempt to route around this with curl per the session's web-fetch policy.
- Practical implication if this path is revisited: IDAES's built-in machinery
  only covers pure fluids. The mixture layer in this repo (reducing
  functions, Bell departure term, composition derivatives) is custom Python
  and would still need its own derivative treatment regardless of whether the
  pure-fluid pieces come from IDAES's exact external functions or from this
  file's own analytic formulas.

### Bell (2023) Citation — Correction/Clarification
- The paper cited as ref [2] in `Helmholtz_R515B_mixutre_math.pdf` (Bell,
  Jaeger, Breitkopf, *Fluid Phase Equilibria* 463, 87-108, 2018) is the
  general theoretical departure-function framework paper — it is **not** the
  source of the specific fitted parameters hardcoded in
  `BELL_2023_R1234ZE_R227EA` / `BELL_2023_DEP_COEFFS`.
- The correct source for those specific numbers: Ian H. Bell, "Mixture Model
  for Refrigerant Pairs R-32/1234yf, R-32/1234ze(E), R-1234ze(E)/227ea,
  R-1234yf/152a, and R-125/1234yf," *J. Phys. Chem. Ref. Data* 52(1), 013101
  (2023), DOI 10.1063/5.0135368. Single-author, NIST. Matches the citation
  already present in this file's docstrings. Freely available NIST PDF:
  https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=935321
- Confirmed via that PDF: Table II (interaction parameters `beta_T`,
  `gamma_T`, `beta_v`, `gamma_v` per pair) and Table VII (3-row departure
  coefficients `n_k,t_k,d_k,l_k` for R-1234ze(E)/227ea specifically) match
  the structure implemented in code.

### CONFIRMED: `beta_v` Transcription Error in `BELL_2023_R1234ZE_R227EA`
- User supplied the actual paper PDF (`013101_1_online_1.pdf`, the JPCRD
  52(1):013101 (2023) article itself). Verified directly via `pdftotext`
  extraction of that file (not an AI-mediated read) — Table 2, row
  "R-1234ze(E)/227ea": `betaT,12=1.001247`, `gammaT,12=0.989180`,
  `betav,12=0.999290`, `gammav,12=1.001581`. Table 7 (departure coefficients
  for this pair) matches `BELL_2023_DEP_COEFFS` exactly, digit for digit,
  for all three rows.
- **Confirmed discrepancy:** code's `BELL_2023_R1234ZE_R227EA.beta_v =
  0.99290` (5 digits after the decimal) is missing a "9" versus the paper's
  printed `0.999290` (6 digits) — a ~0.65% relative difference in that one
  constant. The other three interaction parameters (`beta_T`, `gamma_T`,
  `gamma_v`) and all three departure-function coefficients are correct as
  currently hardcoded.
- Not yet fixed in code (user is editing `mixture_model_one_point.py`
  directly this session; this breadcrumb documents the finding for
  whenever that edit is made). This value feeds `vc12` -> `vred` -> `delta`
  for every mixture-level calculation in the file, so correcting it is a
  small science/parameter fix, not a documentation-only edit.
- Source: Ian H. Bell, "Mixture Model for Refrigerant Pairs R-32/1234yf,
  R-32/1234ze(E), R-1234ze(E)/227ea, R-1234yf/152a, and R-125/1234yf,"
  *J. Phys. Chem. Ref. Data* 52(1), 013101 (2023), DOI 10.1063/5.0135368,
  Table 2.

### Other Findings From This Session
- `mixture_alpha0_alphar_derivs` is dead code: defined but never called
  anywhere else in `mixture_model_one_point.py`. `compute_pressure_enthalpy`
  re-derives the same a0/ar mixture assembly inline instead of calling it.
  Its `a0_mix` correctly omits the ideal-mixing-entropy term (per PDF Eq. 19's
  qualifier that the term is only added "whenever full mixture thermodynamic
  functions are assembled"), consistent with `_mixture_alpha_eval` which does
  include that term for the fuller Table-1 path.
- In `_mixture_alpha_eval`, the `_Tred`/`_vred` values unpacked from
  `_mixture_reduced_state` are reused a few lines later (to recompute `c1`,
  `c2`, `rho_red_mol`) despite the leading-underscore "unused" naming
  convention that IS honored at this function's other two call sites in the
  file (inside `ar_func` and in `compute_table1_properties`'s composition-
  derivative block). Not a bug (values are consistent), but redundant
  computation and a misleading name at that one call site.
- Reconfirmed `Tred`/`vred` from `bell2023_Tred_vred` are mixing-rule
  "reducing" quantities, not the mixture's true critical temperature/volume;
  they only equal `Tc,i`/`vc,i` at the pure-component composition limits
  (verified numerically: `x1=1 -> Tred=Tc1` exactly, `x1=0 -> Tred=Tc2`
  exactly, cross term vanishes at both edges).

## Session Breadcrumb — 2026-08-12 (later same day): Docstring/Comment Pass on Bell (2023) Math Code

User's request: "Please write detailed docstring tying back to the paper,
next to each parameter, define what it is in a comment. This code is very
difficult to follow." Scope: documentation only, no logic changes.

### What was done
Expanded docstrings and added per-parameter/per-variable inline comments
(each tied to a specific PDF equation number from `Helmholtz_R515B_mixutre_math.pdf`
and/or a table number from the actual Bell 2023 paper) throughout:
- `Bell2023PairParams` dataclass — per-field comments tying `beta_T`/`beta_v`/
  `gamma_T`/`gamma_v` to PDF Eqs. 5-8 and paper Table 2 columns.
- `BELL_2023_R1234ZE_R227EA` / `BELL_2023_DEP_COEFFS` — comment block
  documenting the verified paper Table 2 / Table 7 values and flagging the
  `beta_v` discrepancy inline next to the value itself (value left
  unchanged: `0.99290`, not corrected to `0.999290`).
- `bell2023_Tred_vred` — expanded docstring including an explicit
  "Tred/vred are NOT the true mixture critical point" caveat, plus inline
  comments on `theta_T`/`theta_v`, `Tc12`/`vc12`, `Tred`/`vred`.
- `bell2023_departure_alphar` and `bell2023_departure_base` — expanded
  docstrings (PDF Eq. 14/16) and per-line comments on `pref`/`val`/
  `val_tau`/`val_del` and the per-term derivative logic in the loop.
- `mixture_alpha0_alphar_derivs` — docstring now explicitly flags this as
  **dead code** (no call sites in the file; `compute_pressure_enthalpy`
  re-derives the same assembly inline) and notes it differs subtly from
  `_mixture_alpha_eval` (its `a0_mix` omits the ideal entropy-of-mixing
  term). Every local variable commented (`Tc1`/`Tc2`/`MW1`/`MW2`/
  `rhoc1_mol`/`rhoc2_mol`, `c1`/`c2`/`k1`/`k2`, `tau1`/`tau2`/`delta1`/
  `delta2`, the four IDAES evaluator calls, `a0_mix`/`a0_tau_mix`/`ar_mix`/
  `ar_tau_mix`/`ar_del_mix`).
- `_mixture_reduced_state` — expanded docstring (Purpose/Inputs/Outputs per
  variable/Assumptions/Failure modes/References) and full inline comments.
- `_mixture_alpha_eval` — expanded docstring explicitly stating this
  function computes NO entropy (s/R is assembled downstream in
  `compute_table1_properties`), and corrected the user's own truncated/
  imprecise inline comment on `_Tred` (previously: "Reducing (critical
  property for pure fluid) temperature for mixtur" — cut off and wrong,
  since Tred is not a critical property for 0<x1<1). New comment makes the
  non-critical-point caveat explicit. All other locals commented.
- `_bell2023_reducing_derivs_binary` — was untouched before this pass;
  now has a full docstring (explains this is the x1-derivative of
  `bell2023_Tred_vred`'s Tred/vred, in closed form) plus inline comments on
  `DT`/`DV`/`theta_T`/`theta_v`/`dtheta_T`/`dtheta_v`/`Tc12`/`vc12`/
  `dxx_theta_T`/`dxx_theta_v`/`dTred_dx1`/`dvred_dx1`.

### Process note (near-miss, self-corrected)
Mid-pass, a `device_stage_files` call (done to double check the live device
mtime before writing back) landed the device's copy of
`mixture_model_one_point.py` directly into the same container path that was
being hand-edited, silently overwriting the in-progress edits in this
session's working copy. Caught it immediately (grepped for a known-added
comment string; got zero matches) and re-verified: the device file's
content/size (47726 bytes) matched the documented pre-pass baseline exactly,
confirming the user had not made new edits in the interim — so this was safe
to treat as "redo the same edits on the correct base," not a real conflict.
All edits above were then re-applied from scratch onto that fresh base.
**Takeaway for future sessions:** `device_stage_files` overwrites the
container-local file at that path; never call it on a file with unsaved
in-progress `Edit` tool changes still sitting only in the container copy —
stage BEFORE starting a round of edits, not mid-way through.

### Verification before delivery
- `python3 -m py_compile mixture_model_one_point.py` — passed.
- Grepped all previously-verified numeric constants (`beta_T`, `beta_v`,
  `gamma_T`, `gamma_v`, all three `BELL_2023_DEP_COEFFS` rows) post-edit —
  all unchanged from their pre-edit values, including the still-unfixed
  `beta_v=0.99290` bug (intentionally left as-is; only documented).
- Device file mtime checked immediately before delivery/commit
  (1786558391355, 47726 bytes) — unchanged from the pre-pass baseline,
  confirming no unsaved user edits were at risk of being clobbered.
- Delivered via `SendUserFile` and committed to disk via
  `device_commit_files` with `expectedMtimeMs` set to that same value (not
  `force`) — written successfully, zero rejections.

### Still not done in this file (not part of this request's scope)
- The `beta_v` numeric fix itself (`0.99290` -> `0.999290`) — confirmed bug,
  documented in comments, not applied. User said "I will write the
  functions," so this is intentionally left for the user.
- No functional/behavioral changes were made anywhere in this pass; every
  edit was additive (docstrings, `##`-prefixed inline comments) or a
  correction to an existing comment's wording, never to executable code.

## Session Breadcrumb — 2026-08-12 (later same day): Commented Out Dead Code

User's request: "Please comment out dead code." User also independently
added two small inline comments of their own (`## coefficient of fugacity
for fluid-1` / `-fluid-2` on `phi1`/`phi2`, compute_table1_properties)
between the previous breadcrumb and this request -- noted, not reverted.

### What was done
- `mixture_alpha0_alphar_derivs` (the function flagged as dead code in the
  earlier docstring pass) was commented out in place -- every line prefixed
  with `#`, preceded by a header block explaining why, rather than deleted.
  Chose comment-out over delete since (a) the user's request was literally
  "comment out," and (b) the earlier docstring already called this a
  "possible refactor target" worth keeping visible.
- **Before commenting it out, verified system-wide (not just within this
  file) that nothing depends on the mixture_model_one_point.py copy of this
  function.** `grep -rn "mixture_alpha0_alphar_derivs"` across the whole
  repo (via device_bash) turned up 6 other files that import/call a
  function of this exact name: `mixture_pseudo_dome.py`,
  `mixture_pseudo_dome_copy.py`, `mixture_true_vle_copy.py`,
  `mixture_vle_true_reference.py`, `scripts/validate_r515a_reference_suite.py`,
  `scripts/audit_residual_helmholtz.py`. Checked each one's import
  statement: **all six import from `linear_model_codex.py`, not from
  `mixture_model_one_point.py`.** `linear_model_codex.py` has its own
  independent (and apparently actively-used) definition of a function with
  this same name at its own line 733 -- a different file, unaffected by
  this edit. Confirmed via `ast.parse` that `mixture_alpha0_alphar_derivs`
  no longer appears as a live top-level function in
  `mixture_model_one_point.py` after the edit, and `py_compile` passes.

### New finding surfaced by this check: THREE MORE unused functions in this file
While checking call counts for every function to decide what counts as
"dead code," found that these three are *also* never called anywhere in
`mixture_model_one_point.py`, and (like the case above) other repo files
that share their names import from `linear_model_codex.py` instead, not
from here:
- `pure_tau_delta_from_T_rhomol` (line ~206) -- 0 call sites in this file.
- `_nares_from_n` (line ~1553) -- 0 call sites in this file. (Used
  elsewhere in the repo, e.g. `scripts/audit_residual_helmholtz.py`, but
  that script imports its own `_nares_from_n` from `linear_model_codex.py`.)
- `compute_table1_properties` (line ~1585, ~175 lines) -- 0 call sites in
  this file. This is the "fuller Table-1 property evaluator" -- the
  function containing cv/cp/s/speed-of-sound/fugacity/phi/mu, including the
  `phi1`/`phi2` fugacity-coefficient lines discussed with the user just
  before this request. `_cli()` (the module's only entry point) calls
  `compute_pressure_enthalpy` exclusively, never `compute_table1_properties`.

**Did NOT comment these three out in this pass.** Reasons: (1) the user's
request most directly matched the one function already labeled
"CURRENTLY UNUSED" in the prior docstring pass -- treating "dead code" as
referring to that, not yet confirmed to mean "everything with zero call
sites"; (2) `compute_table1_properties` in particular is large (~175
lines), is the function the user was just asking detailed questions about
(`phi1`/`phi2`), and other files in the repo may plausibly still import the
*mixture_model_one_point.py* copies of these three (not yet exhaustively
re-checked with the same rigor as `mixture_alpha0_alphar_derivs` above --
only confirmed that files sharing their names import from
`linear_model_codex.py`, not that literally zero importers exist for the
`mixture_model_one_point.py` copies). Flagged to the user directly in-chat
to confirm before touching these three.

### Verification before delivery
- `python3 -m py_compile mixture_model_one_point.py` -- passed.
- `ast.parse` + walk confirms `mixture_alpha0_alphar_derivs` is absent from
  the set of top-level `FunctionDef` names post-edit (21 functions remain).
- Device file mtime checked immediately before delivery/commit
  (1786559369971, 71880 bytes) -- matched the last-known state, confirming
  the small `phi1`/`phi2` comment the user added themselves was accounted
  for and nothing else had changed underneath this edit.
- Delivered via `SendUserFile` and committed via `device_commit_files` with
  `expectedMtimeMs` (not `force`) -- written successfully, zero rejections.

## Standing Instruction — 2026-08-12: Response Format

User request: "Too much at once. Standing instructions, try to summarize
points in bullets not more than 5, commit this to your memory, update your
context breadcrumbs."

**Going forward, in this project, keep chat responses short: bullet points,
max 5 per response, instead of long prose explanations.** This overrides
the general prose-preference default for this specific project/user going
forward -- applies to technical explanations, diagnostic summaries, etc.
Still fine to include code blocks/commands inline where needed; the
constraint is on the surrounding explanation length/structure, not on
omitting necessary specifics (numbers, file/line refs) entirely -- just
say them tersely, one bullet each, rather than paragraphs.

Also applies to file-edit summaries and breadcrumb-update confirmations
after this point.

## R-515B Liquid-Branch Pressure Validation — 2026-08-12

### Test point used
Real, independently-published Honeywell Solstice N15/R-515B data point:
T=298.15 K (25 C / 77 F), saturated-liquid density rho=1179.8 kg/m3,
w1=0.911 (mass fraction R-1234ze(E), matches Climalife/Honeywell TDS
composition 91.1%/8.9%). Expected Psat approx 497.4 kPa from the
Honeywell P-T table/P-h chart (Reference State on chart: h=200 kJ/kg,
s=1.00 kJ/kg-K, sat. liq. at 0C).

### Finding: ~9.9x pressure over-prediction, root-caused to Z-cancellation
Running `mixture_model_one_point.py`'s CLI at this exact point:
`p_kPa = 4936.55`, `h_kJkg = 236.63` -- about 9.9x the expected ~497.4 kPa.
- Added full diagnostic `print()` block inside `compute_pressure_enthalpy`
  (the only function `_cli()` actually calls) printing every intermediate:
  x1/x2, MWmix, rho_mol, Tc1/Tc2/vc1/vc2, Tred/vred, tau/delta, k1/k2/c1/c2,
  tau1/tau2/delta1/delta2, ar1_del_i/ar2_del_i, dar_del, ar_del_mix, Z.
- Confirmed Z = 1 + delta*ar_del_mix = 0.1983 at this state -- i.e. `1` and
  `delta*ar_del_mix` (approx -0.80) are both large and nearly cancel. This
  liquid-region near-total-cancellation means any small absolute error in
  ar_del_mix gets amplified hugely in Z, and thus in p.
- Independently reimplemented R-227ea's `phi_residual_type=2` formula from
  scratch (own script, actual JSON coefficients from
  `~/Desktop/DVCT_Project/property_package/r227ea.json`) at the printed
  state tau2=1.2574207613617305, delta2=2.873295245082979: reproduced
  ar2_del_i = 1.3389118299412626 to 16 digits, plus a finite-difference
  cross-check confirming the analytic derivative is correct. **Ruled out a
  coding bug in `alphar_idaes_with_derivs`** -- it faithfully implements
  R-227ea's own published EOS at this mapped state.
- k2 = 1.20987 (mixing-rule chain-rule factor mapping mixture delta onto
  R-227ea's own reduced state) maps to an "equivalent" R-227ea density of
  approx 1707 kg/m3 at this T -- flagged as the likely site of the real
  issue (R-227ea's own critical density is only 594.25 kg/m3), but not yet
  confirmed against R-227ea's own true saturated-liquid density at 25C
  (web search inconclusive; recommended next step is
  `CoolProp.PropsSI('D','T',298.15,'Q',0,'R227EA')` in the user's own env).

### Ruled out: NOT a regression in mixture_model_one_point.py
User pushed back hard that nothing material had changed and a prior
`linear_model_codex.py` should already have been validated. Verified
empirically, not by argument:
- Staged `linear_model_codex.py` (47349 bytes, untouched) and
  `r1234ze.json`/`r227ea.json` (from
  `~/Desktop/DVCT_Project/property_package/`) into sandbox.
- Ran `linear_model_codex.py`'s own `compute_pressure_enthalpy("r1234ze",
  "r227ea", 0.911, 298.15, 1179.8)` directly (via idaes-module-stubbing
  technique -- see Technical Notes below): returned
  `p_kPa=4936.545765922068`, `h_kJkg=236.63017973453339` -- **identical**
  to `mixture_model_one_point.py`'s result.
- Diffed code: `alphar_idaes_with_derivs` byte-identical between the two
  files; `BELL_2023_R1234ZE_R227EA` / `bell2023_Tred_vred` /
  `bell2023_departure_alphar` are the same formulas, same still-present
  `beta_v=0.99290` transcription bug (should be 0.999290 per Bell 2023
  Table 2 -- confirmed earlier project, NOT fixed since user owns this
  file: "I will write the functions").
- Conclusion: this is a real, pre-existing property of the Bell (2023)
  mixture model as already implemented -- not something introduced by
  today's rewrite/edits.

### Reconciled with prior breadcrumb history (2026-03-02/03)
Cross-checked against `scripts/audit_residual_helmholtz.py` results
already in these breadcrumbs (lines ~367-398, validated against
`linear_model_codex.py`):
- Pure alphar vs Bell (2023) Table 13: max rel error 1.406e-09 (1st deriv),
  6.509e-09 (2nd deriv) -- machine precision. Consistent with today's
  independent reimplementation check above; pure-fluid layer is doubly
  confirmed correct.
- Binary alpha_r at composition z1=0.4 (near 50/50): known persistent
  -1.014% deviation -- small, pre-existing, unresolved mixing-layer
  imperfection, not zero even in the "validated" module.
- R-515B vs Honeywell PT chart via `mixture_true_vle_copy.py`'s bubble/dew
  VLE solver: mean bias -18.45%, max abs error 42.77% (one run), range
  approx +41.8% to +2.1% (another run) -- much smaller than today's ~890%
  (9.9x) error, and in some runs the OPPOSITE direction (under- vs
  over-prediction).
- Why today's test exposed a much larger error than any prior check:
  (1) composition mismatch -- Table 13 check used z1=0.4, today's is
  x1 approx 0.938, heavily lopsided, and departure-fit accuracy is not
  guaranteed uniform across composition; (2) the VLE solver in the PT-chart
  check finds its own self-consistent (rho_l, rho_v) rather than being
  forced to use the literal real density -- if the EOS surface is bad at
  the true density, the solver can drift to a different (rho_l, rho_v)
  where its own equations still balance, masking the full severity;
  (3) today, for the first time, the literal published real saturated-liquid
  density (1179.8 kg/m3) was plugged in directly with no solver free to
  compensate, landing on a state very close to worst-case for the
  Z-cancellation amplification described above.
- Bias-direction discrepancy (today over-predicts; historical PT sweep
  under-predicts) is NOT yet fully reconciled -- noted as a real open
  question, not resolved.

### Still open / not yet done
- Confirm R-227ea's own real saturated-liquid density near 298 K (via
  CoolProp/REFPROP in user's env) to check whether k2's ~1707 kg/m3 mapped
  state is unrealistically dense for R-227ea specifically.
- Whether Bell (2023)'s departure/mixing fit was ever validated against
  deep-liquid states at all (vs. mainly VLE/near-critical data) -- would
  explain why error is small near VLE/Table-13 conditions but huge here.
- The three other zero-call-site functions flagged 2026-08-12
  (`pure_tau_delta_from_T_rhomol`, `_nares_from_n`,
  `compute_table1_properties`) -- still not commented out, user never
  confirmed.
- `beta_v=0.99290` transcription bug -- still not fixed (user owns the
  file).
- Reconciling over- vs under-prediction direction between today's fixed-
  density test and the historical VLE-solver PT sweep.

### Technical notes for future sessions
- To run code that imports `idaes` (not installed in sandbox): create fake
  `types.ModuleType` objects for `idaes`, `idaes.models`,
  `idaes.models.properties`, `idaes.models.properties.general_helmholtz`
  (with a stub `get_parameter_path`), insert into `sys.modules`, then load
  the target file via `importlib.util.spec_from_file_location` +
  `module_from_spec`. **Must register the module in `sys.modules` BEFORE
  calling `spec.loader.exec_module()`** (e.g.
  `sys.modules["linear_model_codex"] = lmc`) -- otherwise
  `@dataclass(frozen=True)` introspection throws
  `AttributeError: 'NoneType' object has no attribute '__dict__'` because
  `dataclasses._is_type` looks up `sys.modules.get(cls.__module__)`.
- Actual IDAES Helmholtz JSON parameter files live at
  `~/Desktop/DVCT_Project/property_package/r1234ze.json` and `r227ea.json`
  (NOT `r1234ze_parameters.json`/`r227ea_parameters.json`, which are
  Pyomo external-function files: nl_file/expr_map/var_map/param -- not
  consumed by this Python-level EOS code path).
  r1234ze.json basic: R=0.07290727, MW=114.0416, Tc=382.513, rhoc=489.238,
  Pc=3634.86; eos: phi_residual_type=2, last_term_residual=[5,10,16].
  r227ea.json basic: R=0.04890029904335064, MW=170.02886, Tc=374.9,
  rhoc=594.2508657, Pc=2925; eos: phi_residual_type=2,
  last_term_residual=[5,11,18].
- Watch for adding `print()` to `_mixture_reduced_state`,
  `_mixture_alpha_eval`, or `compute_table1_properties` -- none of these
  are called by `_cli()`; only `compute_pressure_enthalpy` is on the
  CLI's call path, so debug prints elsewhere silently never fire. Hit this
  trap twice this session.

## Reconciling the "1.01-1.1% validated" slide vs today's ~9.9x finding — 2026-08-12

User showed a 5-month-old slide ("R515A/B Thermodynamic Package") claiming
R-515B validated against Honeywell datasheet with average bubble/dew error
1.01-1.1%. Confirmed this is real and already in these breadcrumbs, NOT
contradicted by today's finding -- it's a different code path:
- Exact match found at 2026-03-03 entries: bubble-branch MAPE
  `1.122370545966%` (`diagnostics/honeywell_vs_model_bubble_solver_tuned_20260303.csv`),
  dew-branch MAPE `1.021220878991%` (full 49-point,
  `diagnostics/honeywell_vs_model_dew_tuned_full49_20260303.csv`), combined
  package average `1.071795712477%`
  (`diagnostics/honeywell_vs_model_both_tuned_full49_20260303.csv`). This
  is the slide's "1.01-1.1%" number.
- Critical mechanism: this ran through `mixture_true_vle_copy.py`'s
  bubble/dew VLE solver, which SOLVES for its own self-consistent
  `(rho_l, rho_v)` rather than being handed the literal published density
  -- exactly the "solver free to compensate" masking mechanism already
  flagged when explaining today's finding to the user.
- Also required BRANCH-SPECIFIC tuned nonlinear-least-squares
  hyperparameters to reach that number: bubble-tuned profile
  (`method='trf', diff_step=1e-6, loss='cauchy', f_scale=10.0,
  ftol=xtol=gtol=1e-14`) gives bubble MAPE 1.12% but dew MAPE 63.18% on
  the SAME profile; a separate dew-tuned profile
  (`diff_step=1e-8, x_scale=2.0, loss='huber', f_scale=3.0`) was needed to
  get dew down to ~1%. So "1.01-1.1%" is two different tuned solves, not
  one universal evaluation with fixed solver settings.
- Today's test used a different function entirely --
  `compute_pressure_enthalpy` in `mixture_model_one_point.py` -- fed the
  literal real saturated-liquid density (1179.8 kg/m3) directly, with no
  solver and no per-branch tuning free to compensate. That's what exposed
  the ~9.9x liquid-branch error the solver-based validated pipeline never
  surfaces.
- Conclusion: both results are correct and NOT in conflict. The
  underlying EOS/departure-function surface has a real, large error at
  this deep-liquid density; the solver's freedom to pick its own density
  (plus branch-specific tuning) was masking that error in the previously
  validated/presented pipeline. The slide's "model agreement" is real for
  the bubble/dew-solver pathway specifically, not a general certification
  of the direct fixed-density evaluation path tested today.

## Root cause nailed down: extreme density sensitivity, not a code mismatch — 2026-08-12

User asked directly: "shouldn't mixture_model_one_point.py predict the same
values as mixture_true_vle_copy.py?" Answered empirically by running
`mixture_true_vle_copy.py`'s `solve_bubble_at_t` directly (bubble-tuned LSQ
profile from the 2026-03-03 breadcrumbs: `method='trf', diff_step=1e-6,
loss='cauchy', f_scale=10.0, ftol=xtol=gtol=1e-14`) at T=298.15K,
z1=0.938504 (from w1=0.911):

- Solver converges (`status=CONVERGED`, r_P=7.4e-14, r_mu=1.5e-15) to its
  OWN self-consistent equilibrium liquid density:
  `rho_l_molm3=9846.44` -> `rho_l_kgm3=1156.80 kg/m3`, with
  `P_model=506.89 kPa` -- closely matching the CSV row already on file
  (`509.47 kPa @ 25.17C`, `pct_err=1.89%`) and the real Honeywell value
  (~497-500 kPa).
- **Critically: 1156.80 kg/m3 is only ~1.95% below the real experimental
  saturated-liquid density (1179.8 kg/m3)** -- NOT wildly different. The
  earlier hypothesis that the solver "drifts to a very different density"
  is WRONG; the solver's density is very close to the real one.
- Fed the solver's own equilibrium density (1156.80 kg/m3) directly into
  `compute_pressure_enthalpy("r1234ze","r227ea",0.911,298.15,1156.80)`:
  returns `p_kPa=506.89` -- matches the VLE solver's own P_model to 5
  significant figures. **`compute_pressure_enthalpy` and
  `mixture_true_vle_copy.py` fully agree when given the same density.**
  There is NO code/formula discrepancy between the two files.
- The entire ~9.9x gap comes from feeding in the literal real density
  (1179.8) instead of the model's own equilibrium density (1156.80) --
  a mere **1.95% density difference produces a ~9.7x (870%) pressure
  difference** (506.89 kPa -> 4936.55 kPa). This is a direct, now
  numerically pinned-down demonstration of the Z-cancellation
  amplification mechanism described earlier: `dp/drho` is enormous in this
  liquid region, so the EOS is extremely sensitive to small density errors
  right at this state.
- Reframes the open question: it is no longer "why does the model
  disagree with itself" (it doesn't) -- it is "why does the Bell (2023)
  departure/mixing fit's own equilibrium liquid density come out ~2% low
  relative to the real Honeywell/experimental value at this composition
  and T," since that small density error is what the extreme local
  pressure-sensitivity blows up into the ~10x error. Still open: whether
  this 2% density error itself is within/beyond Bell's expected fit
  tolerance, and whether it grows or shrinks at other T/composition.

Script used for this check (not saved to repo, sandbox-only):
loaded `linear_model_codex.py` + `mixture_true_vle_copy.py` via the
idaes-stub + `parameter_path` override technique (see earlier entry),
patched `mtv.LSQ_METHOD/LSQ_DIFF_STEP/.../LSQ_F_SCALE` module globals to
the tuned bubble profile before calling `solve_bubble_at_t` directly, with
seed `rho_l_seed = 0.8*(z1*rhoc1/mw1 + (1-z1)*rhoc2/mw2)` (mol/m3, matching
`run_true_vle_envelope`'s own seed formula).

## Diagnosing the 2% density gap — 2026-08-12 (continued)

Ran two more empirical checks in the same sandbox to isolate the root
cause of the ~2% model-vs-real liquid-density gap identified above.

### Check 1: beta_v transcription bug RULED OUT as the cause
Patched `BELL_2023_R1234ZE_R227EA.beta_v` from `0.99290` to the Bell 2023
Table 2 value `0.999290` in a scratch copy of `linear_model_codex.py` and
re-ran the identical bubble solve (T=298.15K, z1=0.938504, tuned LSQ
profile):
- `beta_v=0.99290` (current): `rho_l=1156.8049 kg/m3`
- `beta_v=0.999290` (corrected): `rho_l=1156.8054 kg/m3`
- Difference: `0.0005 kg/m3` -- **negligible**. This known transcription
  bug is NOT a meaningful contributor to the ~2% density gap. Still worth
  fixing for correctness/Table-2 fidelity, but it will not close this gap.

### Check 2: Solver is NOT seed-robust across composition -- new finding
Ran the same single-shot bubble solve (no continuation, generic seed
`0.8*(z1*rhoc1+(1-z1)*rhoc2)`) across a composition sweep at T=298.15K:

| z1 | rho_l (kg/m3) | P (kPa) | status |
|---|---|---|---|
| 0.9999 | 112.93 | 1175.68 | CONVERGED (spurious) |
| 0.99 | 112.67 | 1171.46 | CONVERGED (spurious) |
| 0.9385 (real R515B) | 1156.80 | 506.89 | CONVERGED (physical) |
| 0.90 | 455.85 | 536.42 | CONVERGED (spurious) |
| 0.70 | 112.58 | 1076.35 | CONVERGED (spurious) |
| 0.50 | 115.34 | 1028.59 | CONVERGED (spurious) |
| 0.30 | 118.77 | 989.06 | CONVERGED (spurious) |
| 0.10 | 122.54 | 955.57 | CONVERGED (spurious) |
| 0.0001 | 769.62 | 124508.72 | CONVERGED (spurious) |

All rows report `status=CONVERGED` (residual gates satisfied,
`r_P, r_mu` tiny) EXCEPT the real R515B composition landed, apparently by
luck of the generic seed being close enough, on the physically sensible
liquid root. Every other composition converged to a clearly nonphysical
root (liquid density ~110-120 kg/m3 is vapor-like, or P>100,000 kPa at
z1~0). **The bubble/dew residual equations (P_l=P_v, mu1_l=mu1_v,
mu2_l=mu2_v) have multiple roots, and `status=CONVERGED` alone does not
guarantee the physical liquid branch was found** -- this is a real gap in
the solver's convergence gating, independent of the EOS accuracy question.
- Important caveat: the actual validated 2026-03-03 pipeline
  (`run_true_vle_envelope`) uses temperature CONTINUATION (each converged
  point seeds the next), which is far more robust than these single-shot
  naive-seed calls -- so this does NOT invalidate the existing 1.01-1.1%
  MAPE results, which were generated via continuation. It does mean any
  *new* single-point solver call (like this diagnostic work, or future
  one-off checks) needs a physically-informed seed, not the generic
  `0.8*rhoc` formula, or it can silently converge to a spurious root while
  still reporting CONVERGED.

### Check 3: Near-pure R-1234ze(E), physically-seeded, for comparison
With a manually physically-informed seed (`rho_l~1170 kg/m3`,
`rho_v~25 kg/m3`, both realistic ballpark guesses for an HFO near 25C),
z1=0.9999 converges to: `rho_l=1163.07 kg/m3`, `rho_v=26.32 kg/m3`,
`P=498.54 kPa`. This is NOT yet cross-checked against a trusted
independent reference (REFPROP/NIST/CoolProp) for pure R-1234ze(E)'s real
saturated-liquid density at 25C -- web search this turn returned candidate
sources (NIST ThermoML, refrigerants.com data sheet, TEGA p-T tables) but
no fetched/verified number yet. **This is the necessary next step**: get a
trusted real rho_l for PURE R-1234ze(E) at 298.15K and compare directly to
this 1163.07 kg/m3 model prediction, to determine whether the ~2% gap
already exists in the PURE fluid EOS/fit (inherited from Bell's own
pure-fluid correlation, not a mixing-rule issue) or is introduced by the
Bell (2023) departure/mixing layer specifically when going from pure to
the real R515B blend composition.

### Updated diagnostic plan (still open)
1. Fetch a trusted real pure R-1234ze(E) saturated-liquid density at
   298.15K (CoolProp/REFPROP in user's own env, or a cited NIST/vendor
   p-T-rho table) and compare to the 1163.07 kg/m3 model value above.
2. If pure-fluid density already shows ~2% gap: root cause is in the pure
   EOS fit itself (upstream of mixing), not fixable by touching the Bell
   mixing-rule coefficients.
3. If pure-fluid density matches real data well (<0.5% say): root cause is
   specific to the Bell (2023) departure/mixing layer at this composition
   -- would point back at the departure coefficients (Table 7) or the
   reducing-function mixing rule (Tred/vred, beta_T/gamma_T/gamma_v -- beta_v
   already ruled out) as the place introducing the ~2% error.
4. Separately (lower priority, solver-robustness hygiene): any future
   single-point (non-continuation) solver diagnostic call should seed with
   a real physically-plausible density, not the generic `0.8*rhoc`
   formula, given the multiple-roots behavior found in Check 2.

## Pure R-1234ze(E) check completed: pure fluid EOS is CLEAN — 2026-08-12

Completed diagnostic-plan step 1 (fetch trusted real pure-fluid data,
compare to model). Source: refrigerants.com R-1234ze Reference Guide PDF
(imperial-unit saturation table).

- Table rows fetched: 75F -> P=69.9 psia, rho_l=72.83 lb/ft3;
  80F -> P=76.0 psia, rho_l=72.27 lb/ft3.
- Linearly interpolated to 77F (=25.00C=298.15K): `rho_l=1163.04 kg/m3`,
  `P=498.77 kPa` (unit conversions: 1 lb/ft3=16.01846 kg/m3,
  1 psia=6.894757 kPa).
- Compared directly to this session's physically-seeded near-pure
  R-1234ze(E) model solve (z1=0.9999, Check 3 above):
  `rho_l=1163.07 kg/m3` (diff **+0.0029%**), `P=498.54 kPa`
  (diff **-0.0455%**).
- **Conclusion: the pure R-1234ze(E) EOS/fit matches this real reference
  table to <0.05% in BOTH density and pressure at 298.15K.** The pure
  fluid layer is now confirmed clean by three independent checks (Bell
  Table 13 machine-precision alphar match, this session's independent
  from-scratch alphar reimplementation, and now a real external reference
  table match).
- **This conclusively narrows the ~2% liquid-density gap (and the ~10x
  pressure amplification it causes) to the Bell (2023) departure/mixing
  layer specifically at the R515B blend composition (x1=0.9385,
  8.9 mass% R227ea)** -- not the pure-fluid EOS. `beta_v` already ruled
  out (negligible effect). Remaining suspects: the Table 7 departure
  coefficients, or the `beta_T/gamma_T/gamma_v` reducing-function mixing
  parameters (only `beta_v` tested so far).
- Natural next check: verify pure R-227ea's own real saturated-liquid
  density near 298K the same way (still not done -- flagged as open since
  the 2026-08-12 diagnostic session began), then test the mixing layer in
  isolation (e.g., turn off/zero the Bell departure term and see how far
  a simple linear/ideal mixing rule alone lands from the real 1179.8 kg/m3
  blend density, to bound how much the departure term should be
  contributing).

Source: https://refrigerants.com/wp-content/uploads/2020/03/R1234ze_RefGuid.pdf

## R-227ea real-density lookup attempt — 2026-08-12 (inconclusive so far)

Tried to mirror the successful R-1234ze(E) real-table check for pure
R-227ea near 298.15K, to test whether ITS pure-fluid layer is also clean
(same diagnostic logic: if both pure fluids are clean, the ~2% blend
density error is conclusively in the Bell mixing/departure layer).

- `refrigerants.com/wp-content/uploads/2020/03/R227ea_RefGuid.pdf` --
  404, guessed URL pattern from the R1234ze/R134a guides does not exist
  for R227ea on that site.
- FM-200 (HFC-227ea fire-suppressant) datasheet
  (`suppressionsystems.com/.../FM-200-Clean-Agent.pdf`) -- WebFetch failed
  on an http/https redirect loop, not retrieved.
- CoolProp fluid page (`coolprop.org/fluid_properties/fluids/R227EA.html`)
  -- has critical point only (Tc=374.9K, rhoc=594.2 kg/m3, matching the
  JSON already on file), no saturation-table density value at 25C.
- Web search surfaced several paywalled/abstract-only academic sources
  (ResearchGate "Saturated densities and critical properties of HFC-227ea",
  ACS JCED papers, NIST ThermoML) but no directly retrievable number yet.

**Not yet resolved.** Most likely path forward: user's own
CoolProp/REFPROP environment (`CoolProp.PropsSI('D','T',298.15,'Q',0,
'R227EA')`) to get this number directly, since public web sources have not
yielded a usable free-text table for R227ea the way they did for
R1234ze(E).

## Session summary delivered — 2026-08-12

Delivered a full session-summary CSV (`r515b_session_summary.csv`) to the
user via SendUserFile for import into their own Google Sheet tracking
record. Contains one row per major finding/action across this entire
debugging arc (docstrings, dead-code cleanup, the ~10x pressure-error
investigation and its full root-cause chain, the beta_v test, the solver
seed-robustness finding, the R1234ze(E) external-table confirmation, and
the still-open R227ea lookup). Breadcrumbs remain the authoritative
detailed record; the CSV is a condensed index of it for external tracking.

## User fixed the beta_v transcription bug — 2026-08-12

User hand-edited `mixture_model_one_point.py` directly and fixed the
long-flagged `beta_v` transcription error at line ~610:
`beta_v=0.99290` -> `beta_v=0.999290` (matches Bell 2023 Table 2,
column "beta_v,12"), with an inline comment marking it as the confirmed
fix. Verified via `device_stage_files` (mtime 1786563865853) +
`python3 -m py_compile` -- compiles cleanly, no syntax issues introduced.

Reminder already logged earlier this session (Check 1 under "Diagnosing
the 2% density gap"): this fix is correct and worth having, but it is
**not expected to close the ~10x liquid-branch pressure gap** -- the
sandboxed test of this exact same change showed the model's own
equilibrium liquid density shifting by only 0.0005 kg/m3 (1156.8049 ->
1156.8054 kg/m3), i.e. negligible. The real ~2% density gap driving the
~10x pressure amplification is still believed to live in the Table 7
departure coefficients or the beta_T/gamma_T/gamma_v mixing parameters,
not beta_v.

## beta_v re-verified against paper text directly — 2026-08-12

User pushed back: "Checked against Table 2. No transcription error."
Re-extracted Table 2 directly from `/tmp/paper.txt` (Bell 2023 paper text)
rather than relying on memory/prior breadcrumb claims:

```
Pair (1/2)              betaT,12    gammaT,12   betav,12    gammav,12   Fij
R-1234ze(E)/227ea       1.001 247   0.989 180   0.999 290   1.001 581   1.0
```

The paper's table formatting inserts a space after the 3rd decimal digit
for every entry in every row (confirmed across all 5 rows of Table 2,
e.g. "0.999 032", "1.007 852", "0.999 710" for other pairs) -- this is a
typographic digit-grouping convention, not two separate numbers. Read
correctly, betav,12 for R-1234ze(E)/227ea = **0.999290** (six decimal
digits), matching exactly what the user already fixed in
`mixture_model_one_point.py` (line 610). This re-confirms, from the primary
source directly (not secondhand from prior breadcrumbs), that `0.99290`
(five digits, missing a trailing 9) was the transcription error, and the
fix applied is correct. Flagged to user in case they were reading a
different PDF edition/printing or a rendering artifact caused the
"290" to look separate/dropped on their end.

## Clarified: delta is NOT an independent solver variable

User asked how the true-VLE solver (`mixture_true_vle_copy.py`) produces
`delta`, and whether it's an independent variable. Traced through
`mix_state()` (line ~285-345): the solver's actual independent/unknown
variables are `rho_l`, `rho_v` (molar densities, parameterized via
sigmoid transforms `u0,u1` for bounded/stable optimization) and `y1` (or
`x1` for dew) -- NOT delta. Inside `mix_state`, at every residual
evaluation: `tred, vred = bell2023_Tred_vred(x1, x2, ...)` is computed
from the (fixed, for that phase) composition, then
**`delta = rho_mol * vred`** (line 335) -- a fully deterministic algebraic
function of the solved density and the composition-dependent reducing
volume. So `delta` is a *derived* quantity recomputed at every solver
iteration, never solved for directly.

## Full equation-by-equation code-vs-paper check — 2026-08-12

User confirmed (independently) that the Bell 2023 pair parameters and
Table 7 departure coefficients are correct. Went further this turn and
checked the CODE'S FORMULAS (not just the coefficient values) against the
paper's actual equations 1-7, term by term, using `/tmp/paper.txt`.

### Reducing function (paper Eq 3-5) vs `bell2023_Tred_vred`
- Eq 3: `Yred(x) = x1^2*Ycrit1 + x2^2*Ycrit2 + 2*x1*x2*[(x1+x2)/(betaY^2*x1+x2)]*Yij`
  Code: `Tred = x1**2*Tc1 + x2**2*Tc2 + 2*x1*x2*theta_T*Tc12` with
  `theta_T = (x1+x2)/((beta_T**2)*x1+x2)` -- **exact match** (same for vred).
- Eq 4: `Tij = betaT,ij*gammaT,ij*sqrt(Tcrit_i*Tcrit_j)`.
  Code: `Tc12 = beta_T*gamma_T*np.sqrt(Tc1*Tc2)` -- **exact match**.
- Eq 5: `vij = betav,ij*gammav,ij*(1/8)*(vcrit_i^(1/3)+vcrit_j^(1/3))^3`.
  Code: `vc12 = beta_v*gamma_v*((vc1**(1/3)+vc2**(1/3))**3)/8.0` -- **exact
  match**.

### Departure function (paper Eq 6-7) vs `bell2023_departure_alphar`
- Eq 6: `alphar_dep = x1*x2*Fij*alphar_ij(tau,delta)` (binary, F_ij=1.0 per
  Table 2). Code: `pref=x1*x2`, returns `pref*val` -- **matches** (F_ij=1.0
  implicit, correct per Table 2).
- Eq 7: `alphar_ij = sum_k n_k*tau^tk*delta^dk*exp(-sgn(lk)*delta^lk)`.
  Code: `expv = np.exp(-(delta**lk))` -- this omits the explicit `sgn(lk)`
  function, but **all three lk values in Table 7 are 1.0** (confirmed
  against paper: k=0,1,2 all have l=1), and `sgn(1)=1`, so
  `exp(-sgn(1)*delta^1) = exp(-delta)` is exactly what the code computes.
  **No discrepancy given this table's actual l_k values** (would only
  matter if some l_k were 0, which none are here).
- Table 7 coefficients cross-checked directly against paper text: k=0
  `(n=-0.057178, t=1.290298, d=1, l=1)`, k=1 `(0.031318, 0.038796, 2, 1)`,
  k=2 `(-0.027496, 2.640532, 3, 1)` -- **exact match** to
  `BELL_2023_DEP_COEFFS` in the code, digit for digit.

### Conclusion: no coding bug found anywhere in the mixing/departure layer
Combined with the already-verified corresponding-states chain-rule mapping
(`tau_i = tau*c_i = Tc_i/T`, `delta_i = delta*k_i = rho/rhoc_i`, both
dimensionally/formulaically consistent with paper Eq 2's definition of
alpha^rCS) and the already-ruled-out `beta_v` bug, **every piece of the
Bell (2023) mixing/departure math has now been verified correct against
the primary source**: pair parameters (Table 2), departure coefficients
(Table 7), the reducing-function formulas (Eq 3-5), the departure-function
formula (Eq 6-7), and the CS chain-rule mapping (Eq 2).

**This rules out the leading hypothesis from earlier this session** (that
the ~2% liquid-density gap lives in the Table 7 coefficients or the
beta_T/gamma_T/gamma_v parameters as a transcription/formula bug). No bug
has been found there after this level of scrutiny.

**Updated leading hypothesis**: the ~2% density discrepancy at this
specific state (x1=0.9385, T=298.15K, deep liquid) may be a genuine, small
residual inaccuracy inherent to the published Bell (2023) correlation
itself in this particular composition/density regime -- not a coding
error in this codebase. Published mixture-model fits are not guaranteed
to have uniform accuracy across the whole composition/density space; nothing
in the paper explicitly claims deep-liquid, highly-asymmetric-composition
density accuracy at the 0.1% level. This would need to be checked against
Bell's own paper for any explicitly reported density-fit residuals/error
bars in this region (not yet checked) to distinguish "known/expected
model limitation" from "still worth investigating further."

## ROOT CAUSE FOUND: composition extrapolation, not a bug — 2026-08-12

Checked the paper's own reported fit statistics and, critically, the
actual composition range of the experimental data used to fit the
R-1234ze(E)/227ea departure function (Section 4.5, Table 12).

### Paper's reported accuracy (Section 4.5)
> "The density measurements cover the full density range and are
> reproduced with an AARD of 0.03% and a maximum deviation of 0.5%."
> "...the new model fits the bubble-point measurements with an AARD of
> 0.05%."

These numbers looked hard to reconcile with our observed ~2% density
error -- until checking WHERE that accuracy was measured.

### Table 12: composition range of the fitting/validation data
| Kind | Source | N | x1 (mole frac. R-1234ze(E)) | T (K) |
|---|---|---|---|---|
| PVT | Fortin (2023) | 164 | **0.33-0.67** | 230.0-400.0 |
| SOS | Rowane (2022) | 313 | **0.33-0.67** | 230.0-345.0 |
| VLE | Outcalt (2021) | 29 | **0.33-0.68** | 270.0-360.0 |

**Every single dataset used to develop and validate this binary departure
function covers only x1 = 0.33-0.68 (roughly 1:2 to 2:1 molar ratio).**

### The real R-515B composition is far outside this range
Our test composition, from the actual Honeywell R-515B blend
(w1=0.911 mass -> x1=0.9385 mole fraction R-1234ze(E)), is nowhere near
the 0.33-0.68 range the model was ever fit or validated against. This is
a **significant extrapolation** -- x1=0.9385 is about 38% of the way
beyond the upper edge of the fitted range, well outside where the
paper's own reported 0.03%/0.5% density accuracy applies.

### Conclusion
**This is the root cause.** Not a coding bug, not a transcription error --
every coefficient, parameter, and formula in the mixing/departure layer
has now been verified correct against the primary source (see prior
entries this session). The ~2% liquid-density error (and the ~10x
pressure amplification it causes via the Z-cancellation mechanism) is a
genuine model-extrapolation limitation: the Bell (2023) R-1234ze(E)/227ea
departure function was never fit or validated anywhere near the real
R-515B composition, so its accuracy there is unconstrained by the
underlying experimental data. This fully closes the investigation loop
started at the beginning of this session (the ~9.9x pressure
over-prediction at the real Honeywell R-515B state point).

### Implication for downstream use
Any property calculation for the *real* R-515B blend (x1~0.9385) using
this Bell (2023) pair correlation is operating in extrapolated territory,
not validated territory -- worth flagging in any paper/report that cites
this model's accuracy for R-515B specifically, since the headline 0.03-
0.05% accuracy numbers do not apply at this composition.

## Why does the VLE solver still fit well despite the extrapolation? — 2026-08-12

User asked the natural follow-up: if x1=0.9385 is a genuine extrapolation
outside the fitted range (x1=0.33-0.68), why did the 2026-03 VLE-solver
pipeline still achieve ~1-2% pressure MAPE at this exact composition?
Answered with a numeric sensitivity comparison at the converged bubble
point (T=298.15K, rho_l=9846.44 mol/m3, rho_v=235.92 mol/m3):

- Computed finite-difference dP/drho on both branches at the same
  converged state:
  - Liquid: `dP/drho = 20420.76 Pa per (mol/m3)`
  - Vapor: `dP/drho = 1823.62 Pa per (mol/m3)`
- Converted to dimensionless elasticity (%change in P per %change in
  rho): **liquid elasticity = 396.7, vapor elasticity = 0.85 -- the
  liquid branch is ~467x more sensitive (in relative terms) than the
  vapor branch.** (Vapor elasticity ~0.85 ~ 1 makes sense: vapor is close
  to ideal-gas-like at this state, p ~ rho*R*T.)
- **Mechanism**: the VLE solver never needs the "correct" liquid density
  as an input -- it solves jointly for (rho_l, rho_v, y1) such that
  P_liquid=P_vapor and both components' chemical potentials match. Because
  the vapor phase is thermodynamically well-behaved (gentle, near-linear
  p-vs-rho relationship) and the phase-equilibrium condition anchors the
  solution, the resulting equilibrium PRESSURE stays close to correct even
  when the same underlying departure-function extrapolation error causes
  the LIQUID density to land ~2% away from the true experimental value.
  The ~2% density error gets "absorbed" almost invisibly on the liquid
  side (needs only a ~0.005% density change to shift P by 2%, per the
  elasticity above) while the vapor side and the coupling constraints do
  the real work of pinning down the correct pressure.
- **This is why pressure MAPE (~1-2%, the metric the 2026-03 pipeline and
  the paper's own AARD numbers use) and density accuracy are NOT
  interchangeable metrics for this pair.** A model can extrapolate poorly
  in liquid density while still reporting excellent bubble/dew pressure
  agreement, precisely because of this ~467x asymmetry in sensitivity
  between the two phases. Directly imposing the literal real liquid
  density (as this session's original diagnostic test did) bypasses the
  solver's natural "vapor-anchored" pressure-pinning mechanism entirely,
  which is why it exposed a completely different, much larger error than
  any prior pressure-based validation metric ever could.

This fully answers the apparent paradox: extrapolation degrades density
accuracy substantially at this composition, but pressure/bubble-point
accuracy is comparatively protected by the physics of vapor-liquid
equilibrium itself, not because the underlying departure-function fit is
actually good at this composition.

## pressure_validated_model.py: fixed beta_v, added entry point — 2026-08-12

Per user request, made the following edits to `pressure_validated_model.py`
(committed to device, mtime verified before each write, no conflicts):

1. **beta_v fix**: `default_interaction()` (line 341) had
   `betarho12=0.99290` -- same bug pattern as `mixture_model_one_point.py`.
   Fixed to `betarho12=0.999290` per Bell 2023 Table 2 (re-verified
   directly from `/tmp/paper.txt` this session -- see earlier entry).
2. **Docstring header rewritten**: Author: Shilpa Narasimhan; Technical
   support: Claude AI and Codex; QA/Testing Responsibility: Shilpa;
   Creation date: (original, unrecorded); re-created 2026-08-12; plus a
   Purpose-of-file paragraph and a breadcrumb note explaining why the
   `__main__` block was added.
3. **Added a `__main__` entry point**: this file previously had zero
   `if __name__=="__main__"` blocks and only one `print()` (an
   error-path debug message inside `solve_rho_mass_for_P`) -- running it
   directly produced no output at all. Added a demo block that calls
   `compute_props_with_chart_h()` at the real Honeywell R-515B reference
   state used throughout this session (T=298.15K, rho=1179.8 kg/m3,
   w1=0.911) and prints T/rho/w1/w2, tau/delta/alpha, p, raw model h,
   the chart offset, and the chart-aligned h.
4. **Verified it runs**: staged a copy + `density_bracketing.py` into a
   sandbox, stubbed `idaes` the same way as prior checks, ran it with
   `runpy.run_path(..., run_name="__main__")` -- confirmed it now prints
   output end-to-end. `python3 -m py_compile` also passes on the live file.

### New finding surfaced while verifying the demo output (NOT yet fixed)
The demo run printed `p = 45972.93 kPa` at the real R-515B state -- even
more wrong than the ~4936.55 kPa seen in the other files at the same
state. Tracing why: this file's own reducing-function implementation
(`reducing_functions()`, line ~266) computes the cross critical volume as

```python
nu12 = betarho12 * gammanu12 * ((nu1 ** (1 / 3) + nu2 ** (1 / 3)) ** 3)
```

but the paper's Eq 5 is
`vij = betav,ij * gammav,ij * (1/8) * (vcrit_i^(1/3) + vcrit_j^(1/3))^3`
(confirmed directly from `/tmp/paper.txt`; `linear_model_codex.py`'s
equivalent line correctly divides by `8.0`). **This file's `nu12` is
missing the `/8.0` divisor entirely** -- an independent, real bug in this
specific file's formula (separate from the beta_v number typo), not yet
fixed since it was outside the scope of what was asked this turn. Flagged
to the user for a decision on whether to fix it.

## Fixed the missing /8.0 divisor in pressure_validated_model.py — 2026-08-12

Applied the fix flagged in the previous entry:
`nu12 = betarho12*gammanu12*((nu1**(1/3)+nu2**(1/3))**3)` ->
added `/ 8.0` divisor, matching Bell 2023 paper Eq 5 and
`linear_model_codex.py`'s equivalent (already-correct) `vc12` line.
Compiles cleanly (`python3 -m py_compile`), committed to device.

**Re-ran the `__main__` demo after the fix -- pressure is still very
wrong.** Before: `p=45972.93 kPa`. After: `p=45830.71 kPa` (delta shifted
from 4.477 to 2.374, so the fix clearly changed the reducing-volume
calculation as expected, but the resulting pressure barely moved and is
still ~92x the real ~497 kPa -- even further off than the ~4936.55 kPa
already-wrong result from `linear_model_codex.py`/
`mixture_model_one_point.py` at the same state).

**This means there is at least one more, separate issue in this file**
beyond the /8.0 divisor -- the fix was correct (verified against the
paper and against the already-validated sibling file) but did not resolve
the order-of-magnitude pressure error. Leading suspect not yet
investigated: `p_h_from_T_rho`'s method of isolating the residual
delta-derivative -- it finite-differences the FULL mixture alpha
(ideal+residual+departure combined) via `alpha_and_partials_tau_delta`,
then subtracts the analytic ideal-gas contribution `1.0/delta` to recover
`a_del_res`. This algebraic isolation trick was reasoned through
analytically this session and looks structurally sound (since each pure
component's `ln(delta_i)` term chain-rules back to exactly `1/delta_mix`
summed over x1+x2=1), but has NOT been numerically cross-checked against
an independent direct calculation of the residual derivative the way the
other files' `ar_del_mix` was verified earlier this session. This file has
not been given the same level of scrutiny as `linear_model_codex.py` /
`mixture_model_one_point.py` and should not be treated as validated.
Flagged to user rather than continuing to fix without being asked, since
this looks like it may have multiple independent issues.

## Scrutinizing pressure_validated_model.py (in progress) — 2026-08-12

Per user request ("Scrutinize it"), doing a full structural + numeric
review of `pressure_validated_model.py`, matching the rigor already
applied to `linear_model_codex.py`/`mixture_model_one_point.py` this
session. Working notes so far (analysis before full numeric verification):

- `reducing_functions()`: now matches Bell 2023 Eq 3-5 exactly after the
  /8.0 fix (already verified in a prior entry).
- `compute_alpha_mix()`'s ideal term includes an explicit ideal
  mixing-entropy contribution: `alpha0_mix = x1*alpha0_pure_1 +
  x2*alpha0_pure_2 + x1*ln(x1) + x2*ln(x2)`. Reasoned through analytically:
  this term is delta-independent (constant at fixed composition), so it
  does NOT affect pressure at all -- p only depends on the delta-derivative
  of alpha. Not the source of the pressure bug, but worth noting as a
  structural difference from the other files (which don't carry this term
  explicitly since they only ever compute the residual delta-derivative
  directly, never round-tripping through a full ideal+residual alpha).
- `p_h_from_T_rho()`'s core trick: it finite-differences the FULL
  alpha_total(tau,delta) = alpha0_mix + alphares_mix + alpha_dep via
  central differences (`alpha_and_partials_tau_delta`), then computes
  `a_del_res = a_del - 1.0/delta` to isolate the residual-only
  delta-derivative (relying on the analytic fact that
  d(alpha0_mix)/d(delta_mix) = 1/delta_mix exactly, derived via chain rule
  through each pure component's ln(delta_i) term with delta_i=delta_mix*k_i
  and x1+x2=1). This reasoning checks out analytically on paper.
- Two remaining candidate bug locations, not yet numerically tested:
  (1) `compute_alpha_res_pure()` -- a from-scratch reimplementation of the
  IDAES phi_residual_type formulas (generalized for types 1-4), independent
  of the already-validated `alphar_idaes_with_derivs` in the other files.
  If this has an error for phi_residual_type=2 (used by both R1234ze(E)
  and R227ea), it would corrupt both the pure terms AND leak into the
  finite-differenced `a_del` used for pressure.
  (2) The finite-difference step size/central-difference implementation
  itself in `alpha_and_partials_tau_delta` (dtau=ddel=1e-6 defaults) --
  unlikely to cause a ~90x error on its own, but not yet ruled out.
- Next step: directly numerically compare `compute_alpha_res_pure`'s
  output against the already-validated `alphar_idaes_with_derivs` at the
  same (tau_i, delta_i) for each pure fluid, and compare this file's
  `a_del_res` against the other files' analytically-computed `ar_del_mix`
  at the identical mixture state, to localize the remaining bug precisely.

## ROOT CAUSE OF THE ~92x PRESSURE ERROR FOUND — 2026-08-12

Numerically compared `compute_alpha_res_pure()`'s output against the
already-validated `alphar_idaes_with_derivs()` at the identical
(tau_i, delta_i) for both pure fluids at the real R-515B test state
(T=298.15K, rho=1179.8 kg/m3, w1=0.911):

```
R1234ze: pvm=-1.9498204099  lmc=-3.0504755055  diff=1.101
R227ea:  pvm=-1.0498301869  lmc=-2.7817490669  diff=1.732
```

Large, unambiguous mismatch on BOTH pure fluids -- confirms the bug is in
`compute_alpha_res_pure` itself, not in the mixing/reducing layer.

### Exact bug, located and confirmed
In `compute_alpha_res_pure`, the function reads `t1 = float(data["t"]["1"])`
once at the top, then:
- **phi_residual_type == 1** (correct): reads each term's own exponent,
  `ti = float(data["t"][str(index)])`, inside every loop.
- **phi_residual_type == 2, 3, and 4** (buggy): every single term in every
  loop uses `tau ** t1` -- the FIRST term's exponent reused for ALL terms,
  instead of each term's own `data["t"][str(index)]`.

Confirmed this matters with real numbers: R-1234ze(E)'s JSON `t` exponents
range from 0.223 to 2.2 across its 16 residual terms (`t1=1.0` is just one
of many distinct values) -- so reusing `t1` for every term completely
corrupts the sum for any fluid using phi_residual_type 2/3/4. Both
R-1234ze(E) and R-227ea use phi_residual_type=2, so this bug hits every
single evaluation in this file for this pair.

**This is the actual root cause of the ~92x pressure error**, not the
/8.0 divisor (that was a real, separate, correctly-fixed bug, but a much
smaller contributor) and not the finite-difference isolation trick (which
was reasoned through and is structurally sound). The
`p_h_from_T_rho`/`alpha_and_partials_tau_delta` finite-difference machinery
is fine; it was just differentiating an already-wrong `alpha_total`
because `compute_alpha_res_pure` was feeding it wrong pure-fluid residual
values for phi_residual_type 2/3/4.

**Not yet fixed** -- reported to user for a decision before editing,
per this session's now-established for pressure_validated_model.py.

## Bug fixed and verified — pressure_validated_model.py now matches the validated model — 2026-08-12

Fixed `compute_alpha_res_pure`'s phi_residual_type 2/3/4 branches: every
`tau**t1` replaced with `tau**ti` reading each term's own exponent from
`data["t"][str(index)]`, matching the already-correct `phi==1` branch and
the validated `alphar_idaes_with_derivs`. Re-ran both verification checks:

- Pure-fluid alphar now matches `alphar_idaes_with_derivs` **exactly**
  (`diff=0.000e+00` for both R-1234ze(E) and R-227ea) at the test state.
- Full `__main__` demo at T=298.15K, rho=1179.8 kg/m3, w1=0.911:
  `p=4936.44 kPa`, `h=236.63 kJ/kg` (raw) -- matches
  `linear_model_codex.py`/`mixture_model_one_point.py`'s already-verified
  result (`p=4936.55 kPa`, `h=236.63 kJ/kg`) to within ~0.002% (residual
  difference is finite-difference truncation noise from
  `alpha_and_partials_tau_delta`'s `dtau=ddel=1e-6` central-difference
  step, not a remaining bug).

**`pressure_validated_model.py` is now confirmed consistent with the rest
of the validated codebase** -- three independent bugs found and fixed
this session (`beta_v` transcription: 0.99290->0.999290; missing `/8.0`
divisor in `nu12`; hardcoded `t1` instead of per-term `ti` in
`compute_alpha_res_pure`'s phi 2/3/4 branches), plus a working `__main__`
entry point and updated docstring header (Author: Shilpa Narasimhan,
Technical support: Claude AI and Codex, re-created 2026-08-12).
`python3 -m py_compile` passes. Committed to device with mtime guard,
zero rejections at each step.

## 2026-08-12 (cont'd): Impact on "the DVCT model" (the cycle models in DVCT_Project)

- User asked how the ~2% liquid-density / ~10x pressure-error finding affects "the DVCT model." Investigated `/Users/snarasi2/Desktop/DVCT_Project` (the parent project folder) and confirmed the repo we've been auditing (`/Users/snarasi2/idaes-hvacr-cycles`) also exists as a nested clone at `DVCT_Project/property_packages_DVRT_code/idaes-hvacr-cycles/`.
- `vapor_compression.py` (both the nested DVCT_Project copy, 31379 bytes, and the live top-level copy, 40971 bytes) is a `SimpleVaporCompressionCycle` built on IDAES's **native single-component** `HelmholtzParameterBlock(pure_component=fluid_name, ...)` plus `CoolProp.PropsSI(..., fluid_name)` for initialization -- i.e. it is a PURE-FLUID cycle model (R134a, R1234ze(E), etc.), not a mixture-aware one.
- Dispatched a subagent to grep every cycle file in the repo (`vapor_compression*.py`, `vapor_compression_plr*.py`, `simple_vcc_codex.py`, `run_plr_cold_storage*.py`, `run_vcc_tests.py`, `test_cop_vs_ambient_r134a_root.py`, `plot_ph_r515b_layers.py`) for any import/reference to `mixture_model_one_point`, `linear_model_codex`, `pressure_validated_model`, `mixture_pseudo_dome`, `mixture_true_vle`/`mixture_vle_true_reference`, or `R515B`. **Zero matches in any cycle file.** Every cycle model uses only IDAES's native pure-component property package (+ CoolProp for pure fluids); none build a custom JSON EOS/StateBlock for R-515B.
- `plot_ph_r515b_layers.py` is confirmed to be a standalone plotting script only -- it reads pre-generated CSV artifacts (e.g. `ph_layer5_pseudopure_overlay_<stamp>.csv`) and renders a P-h chart; it does not call the mixture model directly.
- `idaes_cycle_math.tex` explicitly states the cycle models use IDAES's pure-component Helmholtz package (no mixture mention); `helmholtz_mixture_math.tex` documents the Bell-mixture workflow as a self-contained single-state p-h query tool, with no described integration into a cycle model.
- **Conclusion / answer given to user**: today, the R-515B mixture-model bug/limitation has **zero direct impact** on any cycle model in this repo or in DVCT_Project -- the two tracks are currently disconnected; no cycle file calls the mixture code. The finding becomes relevant only if/when someone wires R-515B into a cycle model (e.g. via the pseudo-pure collapse approach in `mixture_pseudo_dome.py`, feeding a custom property package into `vapor_compression.py`-style code). At that point: COP/pressure/enthalpy-based cycle outputs would likely stay accurate to ~1-2% (same self-consistent-solve protection mechanism already established), but any density-derived downstream quantity computed from the condenser-outlet subcooled-liquid state (charge/mass inventory, receiver sizing) would inherit the ~2% density bias -- exactly the state (`vapor_frac` fixed to 0.0, i.e. the liquid branch) where the Z-cancellation sensitivity is steepest.

## 2026-08-12 (cont'd): DVCT self-consistency clarification + Git checkpoint

- Clarified for the user: for DVCT, P/T/H/S/rho should all be computed together from one self-consistent model solve (same T/composition, solved jointly) rather than substituting an externally-measured density into the pressure equation. This keeps P, H, and S inside the already-validated ~1-2% MAPE regime (same protection mechanism as the bubble/dew pressure checks).
- Explicitly flagged: self-consistency protects P/H/S from being wrong, but does NOT fix rho itself -- the model's own self-consistent liquid rho at R-515B's real composition is still the same ~1.95%-low value (1156.80 vs real 1179.8 kg/m3). Any DVCT output that reports/uses liquid density directly (charge/mass inventory, receiver sizing) still carries this caveat.
- Explained why P is uniquely vulnerable to amplification but H/S are not: P depends on Z = 1 + delta*alphar_delta, a near-total cancellation (~0.2) at this liquid state, so a small absolute error in alphar_delta becomes a huge relative error in Z/P. H and S depend on different combinations (tau-derivative terms, or delta*alphar_delta added to a much larger non-cancelling baseline for h) -- so the same ~2% density/departure-function error should not blow up similarly for H or S.
- Flagged as an open item: only pressure has been independently validated against real Honeywell chart data this session (PT comparison). Enthalpy and entropy have NOT been independently checked against real Honeywell (T,h) or (T,s) chart points the same way -- the "H/S shouldn't blow up" conclusion is a structural/EOS argument, not yet an empirical confirmation. Recommended next step: pull a real Honeywell chart h or s value and repeat the same direct-comparison test used for pressure/density.

### Git checkpoint (repo: /Users/snarasi2/idaes-hvacr-cycles, branch PLR_vanilla_prop)

- Committed (locally; NOT pushed -- this sandbox has no outbound network access, 403 from proxy on `git push`):
  1. `4e4bfb9` "Debugging runs for paper" -- the three pressure_validated_model.py bug fixes (beta_v, /8.0, t1->ti) plus new `R515B_final_validation/` snapshot folder (pressure_validated_model.py, mixture_model_one_point.py, README.md explaining each file and the density limitation).
  2. `fded32f` "Track mixture_model_one_point.py" -- started tracking the previously-untracked root file, noting it's cross-validated against linear_model_codex.py to floating-point precision.
  3. `4b44948` "Rename R515B_final_validation to R515B_props_validated" -- folder rename per user request.
- PROJECT_CONTEXT.md (this file) is intentionally excluded from all of the above commits per user instruction: "We should not be pushing any breadcrumb markdowns."
- New standing instruction from user (2026-08-12): going forward, Claude gives the full add/commit/push command sequence for the user to review and run themselves, with the commit message subject to user approval -- Claude does not execute `git commit`/`git push` directly anymore.
- Immediate follow-up request (pending, not yet executed by Claude per the new standing instruction): user wants the "Debugging runs for paper" commit message amended to include the NIST REFPROP-issues #750 link (https://github.com/usnistgov/REFPROP-issues/issues/750) that independently documents this exact density-extrapolation gap. Since commit 4e4bfb9 is not HEAD (two commits sit on top of it) and nothing has been pushed yet, this requires an interactive rebase (`git rebase -i HEAD~3`, mark it `reword`) run by the user themselves -- command + proposed new message text given to user for approval, not yet confirmed/run.

## 2026-08-12: SESSION PIN -- stopping point, resume diagnosis tomorrow

### Correction to earlier "P/H/S protected vs density exposed" framing
- User correctly pushed back on an oversimplification. Revised understanding: the departure-function extrapolation error (Bell 2023 fit only covers x1=0.33-0.68; real R-515B is x1~0.9385) flows into EVERY quantity built from alpha^r_mix / d(alpha^r)/d(delta) / d(n*alpha^r)/dn_i -- density, vapor composition y1, enthalpy, entropy, chemical potentials. Nothing is actually error-free.
- Pressure is not "protected" because it's more correct -- it's protected because of a structural coincidence: Z = 1 + delta*alphar_delta is a near-total cancellation (~0.2) at the deep-liquid state, so the *relative/percentage* error in P stays small even though the *absolute* departure-function error is the same size everywhere else.
- Enthalpy (h = RT*(1+tau*alpha_tau+delta*alphar_delta)) contains the exact same delta*alphar_delta term as Z, so it picks up the same absolute error as pressure -- it just doesn't get amplified into a large *relative* error because h's total magnitude is much larger than Z's near-zero value.
- Fugacity/chemical potential (f_i = x_i*rho*R*T*exp(d(n*alphar)/dn_i), from `chemical_potentials_analytic` in `mixture_true_vle_copy.py`) has rho entering benignly (via ln(rho), well-conditioned, no cancellation) -- but the residual term d(n*alphar)/dn_i is still built from the same biased alphar_mix/alphar_delta, so mu_i (and hence the solved y1, rho_l, rho_v) directly inherits the extrapolation error. This is literally the mechanism by which the already-documented ~2% density bias gets baked into the bubble/dew solve.
- **New open item**: y1 (equilibrium vapor composition) has never been independently checked against real Honeywell dew-point/vapor-composition data, same gap as h/s.

### p-H dome image review (screenshot provided by user)
- Dome shows saturated liquid/vapor boundaries, quality lines x=0.1..0.9, and a T=70C isotherm.
- Good consistency signal: the two-phase (dashed) segment of the 70C isotherm lands at approximately the same h as the sat-liquid curve (~300 kJ/kg) and sat-vapor curve (~420 kJ/kg) at the same P (~16-17 bar) -- bubble/dew and two-phase isotherm construction agree with each other internally.
- Not yet independently confirmed: whether the dome's h-axis values match real Honeywell chart data (same open item as above). Pressure axis is presumed ~1-2% accurate by extension of the already-validated bubble/dew MAPE, but this specific chart/dataset has not itself been directly checked against Honeywell data.
- **Bug confirmed still present** (matches the original slide's item 5, "Quality lines replicated per Honeywell. Isotherm computations show issues- debugging"): the single-phase branches of the 70C isotherm (subcooled liquid to the left, superheated vapor to the right of the dome) appear as disconnected "hook" curves that shoot up and away near the top of the plot instead of smoothly meeting the dome boundary at the bubble/dew points.
- Working hypothesis (not yet confirmed against code): this is a *different* failure mode from the ~2% composition-extrapolation density bias. Drawing a single-phase isotherm branch requires root-finding rho from P(T,rho,x)=P_target at fixed T -- below the critical point this P(rho) function is S-shaped (stable liquid root, unstable/metastable middle branch, stable vapor root), and near the saturation boundary these roots crowd together. If the single-phase root-finder isn't tightly bracketed to the physically stable branch there, it can jump to the wrong root or lose continuity right at the dome edge -- producing exactly this kind of disconnected-hook artifact. This is a root-selection/robustness issue, not a magnitude/accuracy issue like the density bias.
- **Not yet done**: haven't located or read the actual script that generates the isotherm CSV consumed by `plot_ph_r515b_layers.py` (that file only plots pre-made CSVs, it doesn't compute them) -- needed to confirm or refute the root-selection hypothesis above.

### Git state as of pin
- Repo `/Users/snarasi2/idaes-hvacr-cycles`, branch `PLR_vanilla_prop`. Local history is currently messy after a rebase mishap (two commits both ended up titled "Debugging runs for paper" -- see reflog). A safe recovery recipe (`git reset --soft 3a59890` + one clean recommit with the NIST-link-inclusive message) was given to the user but **not yet confirmed as run**.
- New standing instruction (this session): Claude gives full git command sequences for user approval/execution; Claude does not run `git commit`/`git push` itself going forward.
- Nothing has been pushed to `origin` yet (sandbox has no outbound network access; push must happen from the user's own machine).
- PROJECT_CONTEXT.md (this file) is intentionally excluded from all commits per user instruction.

### Open items heading into tomorrow (priority order suggested)
1. Confirm/finish the git history cleanup (reset --soft + clean recommit + push) -- currently unconfirmed.
2. Locate and read the isotherm-generating script; confirm or refute the single-phase root-selection hypothesis above.
3. Validate enthalpy and entropy against real Honeywell (T,h)/(T,s) chart points (same direct-comparison method already used for pressure/density).
4. Validate vapor composition y1 against real Honeywell dew-point/composition data (newly flagged today).
5. Still pending from earlier: draft REFPROP@NIST.GOV email requesting updated HMX.BNC/R515B.MIX, and/or start the empirical density correction anchored to real (T,rho) data as a stopgap.

## 2026-08-12 (cont'd): Analyzed NIST's newly-issued R515B.MIX + HMX.BNC files -- major finding

User obtained and uploaded the "latest" NIST REFPROP files (R515B.MIX, HMX.BNC) -- presumably the files referenced in the earlier GitHub issue (usnistgov/REFPROP-issues#750). Read both directly.

### R515B.MIX
- Predefined-mixture file. Mole fraction R-1234ze(E) = 0.9385037942295682 -- matches our x1~0.9385 (derived from w1=0.911) essentially exactly. Good independent confirmation of the composition we've been using all session.
- Line 2 numbers (MW=117.4846 g/mol, "Tc"=382.045 K, "Pc"=3601.40 kPa~36.0 bar, 4.2308) are almost certainly initial-guess/bookkeeping parameters for REFPROP's saturation solver, not EOS inputs. MW checks out exactly against a mole-fraction-weighted average of R1234ze(E) (114.04) and R227ea (170.03). Pc~36.0 bar matches the visual estimate from the p-H dome screenshot reviewed earlier.

### HMX.BNC -- KEY FINDING
- File's own changelog (line 30): `12-10-25 MLH, Add interaction model for mixture 1234ze(E)/R227ea to allow better R515B, based on https://doi.org/10.1063/5.0135368` (the Bell 2023 DOI). Confirms R515B support really was added to REFPROP recently (Dec 2025), consistent with the user's recollection and the earlier GitHub issue thread.
- The binary-pair entry (`?R1234ze(E)/R227ea [R1234ZEE/R227EA]`) contains a row explicitly named/documented elsewhere in the file as `B08  Bell (2023) model for R1234ze(E)/R227ea`. Its 5 parameters: **1.001247, 0.98918, 0.99929, 1.001581, 1.0** -- byte-for-byte identical to Bell (2023) Table 2, i.e. IDENTICAL to what we already have coded in linear_model_codex.py/mixture_model_one_point.py/pressure_validated_model.py and already re-verified against the paper this session.
- Also present in the same block: an `XR0` row (`1.0, 0.993523, 1.0, 1.0`) -- confirmed via the file's own `#MXM` model-definition section to be a *different* model type ("XR0  Reducing functions only" -- a plain Kunz-Wagner-style reducing-parameter rule with NO departure/excess function, used as a fallback/simpler model elsewhere in the file for many other pairs). Not clear from the text file alone which of XR0 vs B08 REFPROP actually dispatches to at runtime for this pair -- that's determined by REFPROP's compiled code, not visible in this file.
- **Conclusion**: NIST's official, newly-added REFPROP support for R515B is NOT a separate refit targeting R-515B's real composition -- it is the exact same general Bell (2023) R-1234ze(E)/R-227ea binary model we already implemented and validated, with the exact same Table 2 coefficients, fit over the exact same composition range (x1=0.33-0.68 per Table 12). This independently corroborates that the ~2% liquid-density extrapolation gap found this session is very likely NOT a bug specific to our code -- it is very likely present in the official NIST/REFPROP R515B calculation too, since it runs the identical math at the identical (out-of-range) composition.
- **Recommended next step (not yet done)**: if real REFPROP software (the executable/library, distinct from these text parameter files) is available anywhere, run this exact R515B.MIX at the Honeywell test state (T=298.15 K) and report the liquid density -- this would give an authoritative, independent number to cite in the paper's SI, and would resolve the XR0-vs-B08 dispatch ambiguity noted above.

## SESSION PIN (end of day, 2026-08-12)

Stopping here for the day per user request. Tomorrow's punch list (unchanged from the earlier pin, plus the new NIST-file finding above):
1. Confirm/finish the git history cleanup (`git reset --soft 3a59890` + one clean recommit with NIST-link-inclusive message + push) -- still unconfirmed as of this pin.
2. Locate and read the isotherm-generating script; confirm or refute the single-phase root-selection hypothesis for the disconnected-isotherm-hook bug.
3. Validate enthalpy and entropy against real Honeywell (T,h)/(T,s) chart points.
4. Validate vapor composition y1 against real Honeywell dew-point/composition data.
5. If possible, run the actual NIST R515B.MIX/HMX.BNC through real REFPROP software at the Honeywell test state, to get an independent density number and resolve whether NIST's official implementation reproduces the same ~2% gap.
6. Still pending: draft REFPROP@NIST.GOV email (may now be less urgent given the finding above -- NIST's shipped file appears to use the same model, not an improved one) and/or the empirical density-correction stopgap.

## 2026-08-13: New day, resumed session

### Morning check-in
- Confirmed the git history cleanup from yesterday went through cleanly: `git log` on branch `PLR_vanilla_prop` now shows a single, well-formed commit `125c976 "Debugging runs for paper"` sitting directly on `3a59890`, no more duplicate-message mess. Only `PROJECT_CONTEXT.md`/`.DS_Store` show as locally modified (expected). Item 1 from yesterday's punch list is CLOSED.
- Asked the user where to focus today (isotherm bug / h-s validation / y1 validation / REFPROP cross-check). User chose a different, larger priority: **build the saturation dome and match it against the Honeywell datasheet.**
- Created a 5-item task list for this: (1) extract Honeywell saturation table from TDS PDF, (2) regenerate model dome at w1=0.911 via the validated bubble/dew solver, (3) align reference states (chart-basis h/s), (4) compute per-quantity error metrics vs Honeywell, (5) plot overlay + deliver.

### Honeywell TDS extraction (task 1 -- DONE)
- Source: `eaf1cb70-8077001rastdssolsticen15ltren.pdf` (Solstice N15 / R-515B Technical Data Sheet), page 3.
- **Important scope correction**: this TDS only contains a Temperature(°F) vs Pressure(psig) saturation table -- there is NO enthalpy, entropy, or density data in this document. The "match against Honeywell datasheet" task is therefore a P-T (saturation curve) comparison, not a full h/s/rho validation. The TDS explicitly describes R-515B as "a nonflammable (A1) azeotropic blend" -- consistent with bubble ~ dew (negligible glide), which is why the 2026-03-03 bubble MAPE (1.12%) and dew MAPE (1.02%) came out so close to each other against this same style of data.
- Cleanly parsed via pdfplumber + regex into 86 unique (T_F, P_psig) pairs, T=0..170°F step 2°F, monotonically increasing pressure, no interleaving/parsing errors. This is a larger, freshly-extracted dataset than the 49-point set referenced in the 2026-03-03 breadcrumbs (which may have come from a chart image rather than this table, or was a curated subset).
- Saved to `/tmp/dome_check/honeywell_pt.csv` in this sandbox (T_F, P_psig columns; not yet unit-converted to SI or committed anywhere in the repo).
- Real h/s/rho validation against Honeywell (open items 3-4 from yesterday) will still need either a digitized P-h chart or another data source -- this TDS does not provide it.

### In progress (task 2, not yet complete)
- About to run the already-validated `solve_bubble_at_t` (mixture_true_vle_copy.py, tuned solver settings: `method='trf', diff_step=1e-6, max_nfev=800, x_scale=1.0, loss='cauchy', f_scale=10.0, ftol=xtol=gtol=1e-14`) at w1=0.911 across all 86 converted temperatures, then compare model P against the real Honeywell P at each point -- same methodology as the existing 2026-03-03 checks, but on this freshly-extracted, larger dataset, run fresh today rather than reused from old logs.

### New, larger project direction (stated by user this morning -- important for scoping future work)
User laid out the plan beyond today's dome check:
1. Building and validating this saturation dome against Honeywell completes "vanilla package" validation (i.e., the from-scratch Python implementations: linear_model_codex.py / mixture_model_one_point.py / pressure_validated_model.py / mixture_true_vle_copy.py).
2. **Next major phase**: port this into an actual IDAES property package -- build a custom Helmholtz MIXTURE property package (analogous to IDAES's existing native single-component `HelmholtzParameterBlock`, but for binary/multi-component mixtures).
3. **Long-term goal**: this custom mixture property package is intended to eventually be merged into IDAES-PSE upstream via a pull request. User will provide the PR details/target when we get to that stage.
4. **Design requirement**: the package must be general-purpose -- it should take (a) pure-fluid JSON files (the existing IDAES-Helmholtz-JSON format already used for r1234ze.json/r227ea.json), (b) departure-function terms, and (c) binary interaction parameters as inputs, and build out the full Helmholtz mixture property surface from those -- NOT hardcoded to just R-1234ze(E)/R-227ea/Bell-2023. Because this is aimed at an eventual upstream merge, the user explicitly said "we need to be thorough."
5. **Immediate scope**: for now, focus specifically on the pseudo-pure-fluid treatment of R-515B (i.e., treating the fixed-composition blend as if it were a single substance with its own effective saturation dome/properties) rather than the fully general N-component mixture package -- that generalization is the eventual IDAES-PSE-facing deliverable, but not today's task.
- This reframes the "isotherm bug," "h/s validation," and "y1 validation" open items from yesterday as work that will matter MORE once we're building the IDAES-mergeable package (correctness there needs to be bulletproof for a PR review), not just for the paper's SI.

---
## 2026-08-13 (afternoon) -- Dome validation physics + per-region correction-curve plan

### mixture_dome_validation.py comment/cleanup pass (completed, pushed to device)
- Added exhaustive math-to-code comments (tags M1-M6) throughout the file the
  user authored (`mixture_dome_validation.py`), verified zero logic change via
  py_compile + bit-for-bit output match vs `mixture_true_vle_copy.py`.
- Per user instruction, trimmed the file to saturation-dome-ONLY scope: removed
  the two-phase-region (lever-rule/quality-line) docstring text, the beta-band
  plotting code in plot_envelope, and the `mixture_two_phase_enthalpy` function
  entirely. That calculation belongs to a separate flash/cycle module, not dome
  validation.
- Added a "FULL REFERENCES" bibliography block (full citations, not vague
  name-drops) for every literature reference used in the file's docstrings:
  - [B23] Bell, I. H. (2023). "Mixture Model for Refrigerant Pairs
    R-32/1234yf, R-32/1234ze(E), R-1234ze(E)/227ea, R-1234yf/152a, and
    R-125/1234yf." J. Phys. Chem. Ref. Data, 52(1), 013101.
    https://doi.org/10.1063/5.0135368
  - [LT99] Lemmon, E. W., & Tillner-Roth, R. (1999). "A Helmholtz energy
    equation of state for calculating the thermodynamic properties of fluid
    mixtures." Fluid Phase Equilibria, 165(1), 1-21.
    https://doi.org/10.1016/S0378-3812(99)00262-9
  - [BCL99] Branch, M. A., Coleman, T. F., & Li, Y. (1999). "A Subspace,
    Interior, and Conjugate Gradient Method for Large-Scale Bound-Constrained
    Minimization Problems." SIAM J. Sci. Comput., 21(1), 1-23.
    https://doi.org/10.1137/S1064827595289108 (the scipy "trf" solver method)
  - Every other function's References section either cites one of these by
    tag with a specific note, or is explicitly marked "not an external
    literature citation" when it's really just an internal identity/algorithm
    choice -- no more vague/unclear references.
- Also fixed the M2a comment (x1*ln(x1)+x2*ln(x2)) to acknowledge the EPS_X
  clip happens before ln() is ever evaluated, closing a forward-reference gap
  the user caught (comment described the ideal math without mentioning the
  singularity guard that makes it safe).
- File is live on device at /Users/snarasi2/idaes-hvacr-cycles/mixture_dome_validation.py,
  NOT committed to git (per standing instruction: user commits, I only give
  command sequences on request).
- NOTE: this file changed under me multiple times during this session (user
  editing live: shebang/coding-declaration lines removed, a URL comment
  rejoined) -- always re-staged and matched the user's current formatting via
  expectedMtimeMs guards rather than force-overwriting. No git push performed
  (device_bash has no outbound network access in this sandbox -- confirmed
  again, user must push from their own terminal).

### Physics established this session: WHY dome pressure is (partially)
### protected from the density-extrapolation bias, and where that breaks down
- Root mechanism (recap + refined): Bell (2023) departure function alphar_dep
  fit only over x1 in [0.33, 0.68]; real R-515B is x1~=0.9385 -- extrapolated,
  ~1.95% liquid density bias at T=298.15K, w1=0.911 (established prior session,
  1156.80 vs 1179.8 kg/m3).
- Refined mechanism for why PRESSURE is less exposed than density:
  1. Departure-function terms scale with delta (reduced density) by
     construction -- so their ABSOLUTE size (and absolute bias) is small on
     the dilute vapor branch (small delta_v, near-ideal-gas) and large on the
     dense liquid branch (large delta_l).
  2. Z = 1 + delta*alphar_del is a near-total cancellation (~0.2) specifically
     on the liquid branch at R-515B's real composition -- amplifies whatever
     absolute bias exists there into a large RELATIVE error in Z (and hence
     in rho_l when inverted).
  3. The bubble/dew solve's reported dome pressure is a JOINT/coupled
     condition (P_liq=P_vap, mu_i matched) -- the well-conditioned vapor
     branch anchors this shared value, so the OUTPUT pressure ends up much
     closer to correct than the poorly-conditioned liquid branch alone would
     suggest. This is NOT immunity -- pressure still carries some residual
     bias (observed ~1-4% vs Honeywell) -- it's just far less exposed because
     it's pinned by the better-behaved side of the coupled system.
  4. Density (rho_l specifically) has no such protection -- it's a DIRECT
     readout evaluated exactly where the bias is largest (dense liquid) AND
     amplified by the Z~0.2 cancellation. Hence density shows the full ~2%
     bias while pressure shows only ~1-4%.
- Where the pressure protection BREAKS DOWN: single-phase subcooled liquid
  (off the dome, no vapor phase present at all). There, pressure is evaluated
  or solved purely on the liquid branch alone -- no vapor-branch anchor to
  temper the Z~0.2 cancellation. Expect single-phase-liquid pressure
  predictions to be MORE exposed to the bias than dome pressure is, possibly
  worse than the density error itself. This generalizes the existing M3 code
  comment ("never substitute a real/externally-measured rho into this formula
  directly") -- the self-consistent coupling is what buys pressure its
  protection on the dome; single-phase liquid doesn't have that coupling.
- Single-phase superheated vapor: expected to remain well-behaved (small
  delta, Z~1, small absolute departure bias) same as the vapor side of the
  dome -- no new exposure there.
- Two-phase interior (between bubble/dew at fixed T): NOT computed by this
  file (deliberately out of scope, see above). If built elsewhere, any
  interior state is a lever-rule interpolation of the bubble/dew ENDPOINT
  values (for extensive properties like h) -- so it inherits whatever bias is
  already baked into rho_l/rho_v at those endpoints; pressure inside the
  two-phase region for a MIXTURE is its own value between P_bubble(T) and
  P_dew(T) (bubble and dew are distinct loci for mixtures, unlike a pure
  fluid where they coincide) -- would need its own equilibrium solve, not a
  simple interpolation, if ever implemented.

### Correction-curve plan (discussed, not yet built -- user's own build)
- User proposed: fit an empirical correction curve against digitized
  Honeywell data, SEPARATELY for each region -- one curve for two-phase
  (dome), one for single-phase liquid, one for single-phase vapor -- rather
  than one blanket curve applied everywhere.
- This directly resolves the identifiability problem flagged: a single
  correction curve fit on dome P-T (which is damped/insensitive to the
  density bias per the mechanism above) would under-correct if blindly
  reapplied to density or single-phase-liquid pressure, which are fully
  exposed to the same underlying bias.
- CURRENT DATA GAP (blocking): we only have ONE Honeywell reference dataset
  in hand -- the 86-point saturation T(F)/P(psig) table parsed from the
  Solstice N15 TDS PDF (page 3; delivered to user as honeywell_pt.csv,
  file_uuid 0b3736c5-e708-46df-b3a3-2e6c2ee96e6e). That supports the
  TWO-PHASE/DOME correction curve directly. It does NOT support liquid or
  vapor single-phase curves -- no density, superheated-vapor, or
  subcooled-liquid tables were found on that page.
- OPEN QUESTIONS before liquid/vapor curves can be built:
  1. Does the Honeywell TDS PDF have additional pages (superheated vapor
     table, P-h chart, subcooled liquid data) beyond page 3 that haven't been
     extracted yet?
  2. What is the provenance/source of the 1179.8 kg/m3 reference liquid
     density figure used in the earlier 1.95% density-uncertainty
     calculation (from prior session day)? If independently measured or from
     a non-REFPROP Honeywell publication, it's a legitimate liquid-region
     anchor point. If it came from REFPROP itself, it's circular (same
     extrapolated Bell 2023 model) and can't serve as an independent
     correction target -- need to verify before using it.
  3. Is Honeywell's own P-T table independently measured, or itself generated
     from REFPROP/an EOS? Matters for whether even the two-phase curve is a
     truly independent anchor or partially circular. Not yet confirmed.
- NEXT STEP (deferred, DVCT integration scoping): user stated "I need the
  validated model for use within DVCT" -- i.e., the eventual goal is to make
  this validated R-515B pseudo-pure model usable as a working fluid inside
  the DVCT vapor-compression cycle code (`vapor_compression.py` and variants
  in this repo). User chose "scope the integration now" (conceptual only, no
  code) when asked via clarifying question, but then redirected to continue
  dome validation first -- DVCT-integration scoping is PARKED, not started.
  Recall from prior session: DVCT currently has ZERO exposure to this mixture
  model bug since it only uses IDAES's native PURE-COMPONENT Helmholtz
  package today -- R-515B as a mixture isn't wired in at all yet. Whatever
  interface DVCT expects (Pyomo/IDAES StateBlock vs. simpler standalone
  correlation lookup) has not yet been investigated in this session -- next
  time DVCT integration comes up, start by reading vapor_compression.py in
  /mnt/user-data/uploads/Desktop/DVCT_Project/property_packages_DVRT_code/idaes-hvacr-cycles/
  to see how it currently sources properties before scoping what would need
  to change.

### Immediate next step per user: proceed with dome validation (their own
### script work). No action items on my side right now beyond support as
### requested -- user is building/running the model themselves per standing
### instruction ("I will write the script and run the model. I want to build
### it myself").

---
## 2026-08-13 (afternoon, cont'd) -- MAJOR FINDING: multiple spurious roots, not just a 2% bias

- User ran the new dome CLI, got a kinked/S-shaped bubble line in P-h space.
  Verified via full 80-point sweep comparison: mixture_dome_validation.py and
  mixture_true_vle_copy.py give byte-for-byte identical output (P, rho_l,
  status) across the entire T range -- NOT a regression from today's comment
  edits. The underlying (already "validated") math produces this on its own.
- User uploaded two older bubble CSVs for comparison
  (r515b_true_vle_bubble_20260303.csv, r515b_true_vle_bubble_copy.csv,
  file_uuids 260e5ddc-28dd-42cc-b211-4d6cbd59d4a2 and
  b0b172ee-ece6-463e-85eb-77c9afe1e1c5). Found:
  - The 20260303 CSV has a HARD DISCONTINUOUS JUMP: rho_l goes 430.4 -> 684.1
    kg/m3 between T=276.8K and T=277.7K in a single step, with P actually
    DECREASING (4.46->3.97 bar) despite T increasing. This is a root-switch
    mid-continuation, not smooth physics -- the "old, correct-looking" dome
    already had this bug baked in; it just wasn't caught/inspected closely
    at the time.
  - The "_copy" CSV stays on a THIRD, even-lower, nearly-flat branch
    (415.9-473.1 kg/m3) across its whole T range -- likely an earlier/buggier
    run predating some fix, not a clean reference either.
  - Today's fresh continuation sweep (T=250-360K) is ALSO non-monotonic:
    rho_l goes 692->696(local max @258K)->679(local min @287K)->778(local
    max @323K)->650 kg/m3 -- a double-hump, not a physical liquid density
    curve.
  - Above ~T=280K, today's sweep and the old 20260303 CSV's "high branch"
    agree closely (both ~715-780 kg/m3 in the 290-330K range) -- so it's not
    a full regression, it's specifically the LOW-T region (and any
    non-continuation single-point solve) that's unstable.
  - CRITICAL: none of the three distinct roots found so far (430, 460, or
    680-780 kg/m3) are anywhere close to the real ~1180 kg/m3 Honeywell-
    adjacent liquid density. This is 40-60% off, not the ~2% composition-
    extrapolation bias quantified earlier in the session. The earlier 1.95%
    figure (1156.80 vs 1179.8 kg/m3, from a DIFFERENT code path/methodology,
    likely pressure_validated_model.py or mixture_model_one_point.py, not
    this mu-equality continuation solver) needs re-examination -- it's not
    yet clear that figure was even measuring the same thing as this solver
    produces.
  - Also recall from earlier in this session: a single-point smoke test
    (T=298.15K, w1=0.911) using a cruder seed heuristic
    (rho_l_seed=0.8*(z1*rhoc1+(1-z1)*rhoc2)) gave YET ANOTHER value: 460.3
    kg/m3 -- a fourth distinct self-consistent root at essentially the same
    T. Confirms extreme seed-sensitivity: which of several spurious roots
    you land on depends entirely on the starting guess, not just on T.
- CONCLUSION: this is a genuine solver root-selection / multiple-roots
  problem in the (unmodified, shared) mu-equality 3-equation system at
  R-515B's real (out-of-Bell-2023-fit-range) composition -- not something
  introduced by today's comment/trim edits, and not fully explained by the
  previously-quantified ~2% extrapolation bias. The extrapolated departure
  surface may be poorly-conditioned/multi-valued this far outside its fitted
  x1 range, producing multiple self-consistent-but-nonphysical minima that
  the continuation and single-point solves can each separately fall into.
- NEXT: user asked "why am I unable to reproduce the dome from before" --
  answered that the "before" dome likely either (a) came from a different
  script (pressure_validated_model.py/mixture_model_one_point.py, single
  fixed-point validated model, not this continuation-based mu-equality
  solver) or (b) had the same jump bug present but unnoticed. Root-selection
  robustness (better seeding/bracketing, or detecting and rejecting the
  low-density branch as nonphysical) is the actual next problem to solve
  before dome validation can proceed -- NOT yet started, awaiting user
  direction (this is their own script/build per standing instruction).

---
## 2026-08-13 (afternoon, decision) -- Pivoting to pseudo-pure (azeotrope) assumption

- User decided: adopt the pseudo-pure assumption (x1=y1, fixed composition
  both phases) for R-515B dome modeling, justified by R-515B being
  manufactured/marketed as a near-azeotropic blend (Honeywell TDS: "a
  nonflammable (A1) azeotropic blend").
- This means `mixture_pseudo_dome.py` (2-unknown solve: rho_l, rho_v; P_l=P_v
  and g_l=g_v at fixed x1) is now the INTENDED working model going forward,
  not a placeholder to be replaced -- its header TODO ("Replace with full
  mixture fugacity-equality flash formulation") is now stale and should be
  revisited/reworded to reflect this as a deliberate choice, not a shortcut
  to fix later.
- Supporting evidence (from old 20260303 true-VLE CSV, high-T tail only,
  since low-T rows there are contaminated by the root-jump bug): at
  T=358-360K, x1=0.9385 (fixed) vs solved y1~=0.9415-0.9418 -- a small
  (~0.3%) but nonzero gap. Have NOT yet confirmed this gap stays small
  across the full T range of interest (250-360K) -- flagged this to user,
  they chose to proceed with the pseudo-pure assumption without that check
  for now.
- Practical effect: pseudo-pure's 2-unknown system (rho_l, rho_v only) should
  sidestep the free-y1-driven multiple-roots/branch-jumping problem found
  in mixture_dome_validation.py's true VLE solver -- much simpler root
  topology, likely far more numerically robust for the dome sweep.
- NOTE: pseudo-pure's coexistence check (g_l=g_v, single scalar condition)
  is a slightly WEAKER condition than full component chemical-potential
  equality (mu1_l=mu1_v AND mu2_l=mu2_v) even at true azeotropic
  composition -- worth keeping in mind if any downstream check seems too
  easy to satisfy.
- STATUS: this is the user's own script/build (per standing instruction --
  they write/run it themselves). No code written by me. Awaiting user's
  next move (likely: run mixture_pseudo_dome.py's CLI over the Honeywell T
  range and compare against the Honeywell P-T table).

---
## 2026-08-13 (afternoon) -- Ran mixture_pseudo_dome.py, T=250-360K, w1=0.911, n=140

- 129/140 converged. Low-T end (250-259.5K, 11 of first 12 points) mostly
  DIVERGED -- likely a seeding/initialization issue at the start of
  continuation, not yet diagnosed further.
- From T=260.3K onward: ALL converged, and rho_l is smooth + monotonic-ish
  (dips slightly 260-272K, then rises smoothly to 517.5 kg/m3 at T=360K) --
  NO jumps, NO double-hump. Pseudo-pure (2-unknown) formulation is far more
  numerically well-behaved than the true-VLE (3-unknown, free-y1) solver, as
  expected.
- BUT: rho_l magnitude across the whole converged range is 429-517 kg/m3 --
  even FURTHER from the real ~1180 kg/m3 Honeywell-adjacent liquid density
  than the true-VLE solver's "high branch" (680-780 kg/m3) was. Pseudo-pure
  fixed the shape/robustness problem but made the magnitude problem worse,
  not better.
- P at T~298.3K = 4.71 bar (471 kPa) -- have not yet cross-checked this
  against the Honeywell P-T table's interpolated value at that T; should
  do that next.
- STATUS: reported to user, awaiting direction on whether to debug the
  low-T divergence, investigate the magnitude gap, or something else. Not
  yet clear whether the magnitude gap here is the same departure-function
  extrapolation bias mechanism discussed earlier, or a units/setup issue
  specific to this file (rho_l this far off, this consistently, across the
  ENTIRE converged range, is worth double-checking rather than assuming
  it's "just" the known extrapolation bias).

---
## 2026-08-13 (afternoon) -- Codex's tuned solver knobs are already in place; pivoted to pseudo-pure fix

- Checked mixture_dome_validation.py's LSQ_* constants against Codex's
  2026-03-03 "tuned_solver" settings (trf, diff_step=1e-7, max_nfev=800,
  x_scale=1.0) -- already identical, already hardcoded. Confirmed these
  knobs only ever bought CONVERGENCE (42/49 -> 49/49), never accuracy
  (MAPE stayed ~16-19% throughout that tuning). Re-running true-VLE with
  same knobs, no thermo changes, will NOT close the Honeywell gap.
- User decided (again, reconfirmed): proceed with pseudo-pure
  (mixture_pseudo_dome.py) as the working model, azeotrope-justified.
- Diagnosed the low-T divergence found earlier in mixture_pseudo_dome.py
  (11/12 points DIVERGED, T=250-259.5K): isolated a specific failing point
  (T=252.4K, no continuation seed) and found r_mu already tiny (1.86e-5)
  but r_P huge (1.33) when scipy.optimize.root(method="hybr") stalls
  ("not making good progress"). Root cause: badly scaled/mismatched
  residuals confusing hybr's dogleg step control.
- FIX IDENTIFIED (convergence-only, zero thermodynamic changes, mirrors
  Codex's approach exactly): swap `method="hybr"` -> `method="lm"`
  (Levenberg-Marquardt) in the `root(...)` call inside
  `pseudo_saturation_point_at_t` (~line 207). Tested in sandbox: same
  T=252.4K point goes from r_P=1.33 (fail) to r_P=-1.4e-13 (machine zero).
  Full 140-point sweep (T=250-360K): 140/140 CONVERGED with this one-line
  change, vs 129/140 before. Low-T branch smooth: 438.8 kg/m3 @250K down to
  430.3 kg/m3 @262K, no jumps.
  IMPORTANT CAVEAT (not yet resolved): even with 140/140 convergence, the
  absolute rho_l magnitude (~430-520 kg/m3 across the whole range) is still
  far from the real ~1180 kg/m3 Honeywell-adjacent density -- this fix is
  ONLY a convergence-robustness fix, exactly like Codex's tuning was. It
  does NOT address the magnitude/accuracy gap.
- User has NOT yet applied this edit themselves (per standing instruction,
  they write/edit their own script) -- I gave them the exact before/after
  line via chat, they were about to make the change when this breadcrumb
  was written.
- Also clarified for user: the `y` in `residual(y)` inside
  mixture_pseudo_dome.py is NOT vapor mole fraction y1 -- it's just the
  solver's raw 2-element unconstrained parameter vector
  [log(rho_v), log(gap)]. There is no y1 variable anywhere in this file;
  the pseudo-pure/azeotrope assumption (y1=x1) is enforced structurally by
  passing the SAME x1,x2 into _mixture_state_from_rhomol for both the
  liquid and vapor evaluations, never introducing a separate vapor
  composition term at all.
- NEXT STEPS (not yet done): (1) user applies the hybr->lm edit and reruns
  the CLI to confirm full convergence on their end; (2) once confirmed,
  compare the pseudo-dome's P-T output against the Honeywell table (same
  comparison methodology used for the true-VLE 20260303 runs); (3) the
  rho_l magnitude gap (~430-520 vs ~1180 kg/m3) still needs its own
  investigation -- not yet started, and NOT explained by the composition-
  extrapolation bias mechanism discussed earlier this session (that
  mechanism was quantified at ~2%, not ~55-65%).

---
## 2026-08-13 (evening) -- Diffed current file against Codex's actual smooth-dome commit

- User's explicit request: edit mixture_dome_validation.py to do what Codex did
  to get the right (smooth) answer -- NOT the pseudo-dome file. (Cleared up
  after repeated back-and-forth confusion earlier in the session.)
- Found Codex's real "smooth dome" commit via git history on
  mixture_true_vle_copy.py: af2ae3d ("Checkpoint, full working saturation
  dome for non-azeotropic blends") -- matches the manifest hash referenced
  in the 2026-03-03 STEP-6 breadcrumb.
- Diffed af2ae3d's mixture_true_vle_copy.py against the LIVE
  mixture_dome_validation.py. Current file ALREADY has every mechanism
  af2ae3d has: sigmoid-bounded density mapping (_scaled_sigmoid,
  _rho_param_to_states, guarantees rho_l>rho_v structurally), the 4-seed
  retry loop (_attempt), separated bubble/dew continuation chains, same
  LSQ_* knobs. No missing machinery.
- Ran current live file's own solver (T=250-360K, n=140, w1=0.911) with
  retry-index + jump + crossover instrumentation:
  - 140/140 CONVERGED, every single point on retry=0 (first/continuation
    seed), zero fallback attempts ever used.
  - rho_l(T): zero jumps >50 mol/m3 (checked); h_l(T): zero jumps >500 J/mol;
    bubble/dew crossovers (h_v<h_l): 0.
  - CONCLUSION: this file's bubble branch, run this way, is already smooth
    end-to-end. The S-shaped kink the user saw is NOT reproduced by this
    sweep.
- Two live hypotheses for the discrepancy (unresolved):
  (1) user's kinked plot used different CLI args (T range/n/w1) than tested
      here -- need their exact command to reproduce.
  (2) it's a plotting/ordering artifact, not a solver artifact -- 2026-03-03
      breadcrumbs explicitly flagged "current isotherm crossover issue is
      polyline parameterization/order driven" as a known failure mode
      distinct from root-jumping.
- STATUS: blocked on user providing the exact command/args that produced
  their kinked P-h plot. No edits made to mixture_dome_validation.py --
  diagnostics only, nothing to fix yet since the file's own solver output
  came back clean.

---
## 2026-08-13 (evening) -- CONFIRMED FIX: bubble-branch S-kink root cause found and tested

- User uploaded the real Solstice N15 (R-515B) Honeywell TDS PDF (primary
  source, dated 05/26). Confirms Tc=228F=382.0K -- matches this session's
  earlier Kay's-rule Tc_mix estimate (382.04K) almost exactly. Also
  confirms Honeywell's own published P-T table only goes to 170F (349.8K)
  -- explains why the user's CLI Tmax=349.8 was chosen (matches datasheet
  range exactly).
- Digitized the datasheet's page-3 P-T table (86 points, 0-170F, 2F steps)
  and ran the CURRENT live file's run_true_vle_envelope at those exact T's:
  86/86 converged, smooth monotonic error trend, MAPE 32.18% (rel err
  76.99% at T=0F down to 3.83% at T=170F). This IS the real accuracy
  baseline against a confirmed-primary-source table -- higher MAPE than
  Codex's old 16.15% only because this table includes many more very-low-P
  points near 0F where tiny absolute Pa error -> huge relative %. Same
  underlying bias shape as Codex's old result, just weighted differently.
- CORRECTED an earlier mistaken diagnosis: I had claimed the P-h "S-shape
  kink" was h_l(T) turning over (non-monotonic) near the critical point.
  Directly tested this -- FALSE. h_l(T) is monotonically increasing
  throughout T=250-360K under current settings; there is no turnover.
- ACTUAL root cause (confirmed by directly generating and viewing the PNG
  for the user's exact CLI command, `--w1 0.911 --Tmin 255.4 --Tmax 349.8
  --n 80`): h_l(T) is monotonic but has an extremely steep, narrow
  acceleration/deceleration around T~282-298K (P~6-11 bar) -- slow, then
  near-vertical, then slow again. This reads as a visible "S" kink in
  log(P) vs h space even though nothing is discontinuous or non-monotonic.
  Reproduced this exact kinked plot at both n=80/Tmax=349.8 AND
  n=200/Tmax=360 with the file's CURRENT (shared, single-profile) LSQ
  settings.
- Compared directly against Codex's actual smooth branch_tuned PNG
  (`verification/r515b_true_vle_envelope_clean_dome_branch_tuned_20260303.png`,
  metadata: `verification/r515b_true_vle_dome_branch_tuned_metadata_20260303.json`,
  Tmin=250K/Tmax=360K/n=200, 200/200 converged both branches, 0
  crossovers) -- Codex's dome at the SAME T range has NO kink.
- Diffed the two runs' solver settings. Codex used TWO DIFFERENT LSQ
  profiles, one per branch (current live file uses ONE shared profile for
  both):
  - bubble profile: method=trf, diff_step=1e-6, max_nfev=800, x_scale=1.0,
    ftol=xtol=gtol=1e-14, loss=cauchy, f_scale=10.0
  - dew profile: method=trf, diff_step=1e-8, max_nfev=800, x_scale=2.0,
    ftol=xtol=gtol=1e-10, loss=huber, f_scale=3.0
  (current file's single shared profile: diff_step=1e-7, ftol=xtol=gtol=
  1e-12, loss=linear, f_scale=1.0 -- this is the LSQ_* block already
  confirmed present weeks ago, but it was NEVER split per-branch.)
- TESTED AND CONFIRMED THE FIX: monkey-patched the live module's LSQ_*
  globals to Codex's bubble profile ONLY (cauchy/f_scale=10/diff_step=1e-6/
  tol=1e-14), reran both n=80 (Tmax=349.8) and n=200 (Tmax=360) sweeps and
  replotted -- S-kink is GONE in both, output matches Codex's smooth PNG
  almost exactly (enthalpy range ~334-405 kJ/kg bubble at n=200, same as
  Codex's 332-403 kJ/kg).
- FIX NEEDED (not yet applied to live file -- user's own edit per standing
  instruction): split the current single LSQ_* block into two named sets
  (bubble keeps the bare LSQ_* names since solve_bubble_at_t already
  reads those; add LSQ_*_DEW constants with Codex's dew values) and update
  solve_dew_at_t's _attempt()'s least_squares(...) call to reference the
  _DEW versions. No change needed to solve_bubble_at_t's call site itself
  (bare LSQ_METHOD/LSQ_DIFF_STEP/etc. just need their VALUES changed to
  Codex's bubble profile).
- STATUS: root cause fully confirmed via direct before/after plot
  comparison (not just CSV inspection). Walkthrough given to user; awaiting
  their edit + rerun to confirm on their end.

---
## 2026-08-13 (evening) -- FIX APPLIED: dew branch wired to LSQ_*_DEW profile

- User had already applied Step 1/2 (bubble LSQ_* values -> cauchy profile)
  and Step 3 (added LSQ_*_DEW constants block) themselves before this edit.
- At user's explicit request ("Can't find that function make these edits"),
  I made the one remaining edit myself: in solve_dew_at_t's _attempt()
  (~line 1270), changed the least_squares(...) call's 9 keyword args
  (method/ftol/xtol/gtol/max_nfev/x_scale/diff_step/loss/f_scale) from the
  bare LSQ_* names to the LSQ_*_DEW names. solve_bubble_at_t's call site
  untouched (still reads bare LSQ_* names, which the user already pointed
  at the bubble/cauchy profile).
- Verified via py_compile (clean) and a full smoke run of the user's exact
  CLI command (--w1 0.911 --Tmin 255.4 --Tmax 349.8 --n 80): 80/80 bubble
  converged, generated p-h plot -- S-kink confirmed GONE, bubble enthalpy
  range ~334-400 kJ/kg, matches Codex's smooth reference shape.
- Edited file sent back to device and committed to
  mixture_dome_validation.py (mtime-guarded, no force).
- STATUS: fix applied and verified end-to-end. Next natural step (not yet
  started): re-run the real accuracy check (P-T vs the digitized Honeywell
  TDS table, 86 points) with this fixed file to confirm the ~32% MAPE bias
  finding from the previous breadcrumb entry is unchanged (expected --
  this was a shape/kink fix only, not a thermodynamic/accuracy change).

---
## 2026-08-13 (evening) -- MAJOR FIX: bubble branch was converging to a spurious near-equal-density root all session

- Discovered while checking the pseudo-pure conversion's density output:
  bubble branch's "liquid" and "vapor" densities were nearly IDENTICAL
  (e.g. rho_l=602.4, rho_v=601.4 mol/m3 at T=255.4K) -- not real phase
  separation. Passed the old strict gate (`rho_l > rho_v*(1+1e-8)`)
  because that gate only rejects EXACT equality, not "suspiciously close."
  The normalized P/mu residuals are trivially near-zero for ANY two nearby
  single-phase states (continuity), so this degenerate root satisfies the
  gate to ~1e-10 -- six orders of magnitude inside the 1e-6 cutoff.
- Compared against the SAME file's dew branch at the same T: dew was
  finding a properly-separated pair all along (rho_l~11064 mol/m3, ~1300
  kg/m3) -- much closer to the real R-515B liquid density than anything
  the bubble branch had produced all session (the ~40-65% magnitude gap
  flagged repeatedly this session, never explained, was this bug).
- FIX applied to `mixture_dome_validation.py`:
  1. New constant `RHO_SEPARATION_MIN_RATIO = 3.0`; both solve_bubble_at_t
     and solve_dew_at_t's strict convergence gate now requires
     `rho_l > rho_v * RHO_SEPARATION_MIN_RATIO`, not just `rho_l > rho_v`.
  2. `run_true_vle_envelope` now solves DEW first at each T (was: bubble
     first), then passes that T's converged dew (rho_l, rho_v) into
     solve_bubble_at_t as a new `dew_rho_l0`/`dew_rho_v0` param.
  3. In solve_bubble_at_t, the dew-seeded attempt is now tried FIRST
     (originally tried it last, which still let bubble's own continuation
     lock onto a real-but-wrong branch before ever reaching the dew seed).
- VERIFIED (T=255.4-349.8K, n=80, w1=0.911):
  - 80/80 bubble converged, retry=0 on all 80 (dew seed accepted
    immediately every time).
  - rho_l/rho_v separation ratio now 8.2-198.8 (was ~1.002).
  - bubble and dew liquid densities now agree to ~1-2% at every T (e.g.
    T=298.4K: bubble 1155.8 kg/m3 vs dew's close match) -- consistent with
    near-azeotrope expectation, unlike before.
  - At T=298.4K (~77F): rho_l=1155.8 kg/m3 vs Honeywell-datasheet-quoted
    real value 1179.9 kg/m3 -- 2.0% off. This FINALLY matches the ~2%
    Bell-2023 composition-extrapolation bias theorized early this session,
    not the ~40-65% gap that had been unexplained since.
  - Re-ran the 86-point digitized Honeywell P-T comparison (same table as
    the earlier 32.18% MAPE finding): MAPE dropped from 32.18% -> 1.84%,
    max abs error 76.99% -> 4.19% (worst point, T=0F), best point 0.43%
    (T=170F). This is a real accuracy fix, not just a plotting/shape fix.
  - P-h plot: still smooth, no kink, enthalpy range now ~180-315 kJ/kg
    (bubble)/~365-420 kJ/kg (dew) -- much narrower/more physical than the
    prior wrong-density branch's ~334-405 kJ/kg bubble range.
- File sent back to device and committed (mtime-guarded).
- OPEN: mixture_dome_validation_pseudo_pure.py (the single-function 2-
  unknown pseudo-pure conversion made earlier this session) inherits the
  SAME degenerate-root risk but has no separate dew branch to borrow a
  seed from (bubble and dew are literally the same call in that file).
  Has NOT yet been fixed or re-verified against this same bug -- next
  step, not yet started.

---
## 2026-08-13 (evening) -- Catch-up entry: pseudo-pure conversion of mixture_dome_validation_pseudo_pure.py

(Logging this now -- happened chronologically BEFORE the "MAJOR FIX" entry
above, between the S-kink fix and discovering the degenerate-root bug.)

- User created `mixture_dome_validation_pseudo_pure.py` as a copy of the
  (S-kink-fixed) true-VLE file, intending to convert it to the pseudo-pure/
  azeotrope assumption. On inspection it was still byte-identical in
  substance to the true-VLE file -- still 3 unknowns (rho_l, rho_v, free
  y1/x1), nothing pseudo-pure yet despite the filename.
- Explained the conversion: drop the free composition unknown, force
  x1=y1=z1 in BOTH phases (not as a soft residual -- as a hard substitution,
  never introducing a separate x1/y1 variable at all), and replace the
  3-equation system (P equality + mu1 equality + mu2 equality) with a
  2-equation system (P equality + Gibbs-energy equality g_l=g_v). Flagged
  that g_l=g_v is a slightly WEAKER test than full mu1/mu2 equality, but is
  the only independent condition left once composition is shared.
- User asked "Shouldn't x1=y1 be a constraint?" -- clarified hard
  substitution (no separate variable to violate) is strictly better than a
  soft residual (which only gets close to zero, not exactly), and is
  justified here specifically because the Honeywell R-515B datasheet lists
  "Zero glide" as a stated property (bubble T == dew T at fixed P, i.e.
  x1=y1 exactly, not approximately).
- At user's explicit request ("Lets do this. Can you key in these
  changes?"), edited `mixture_dome_validation_pseudo_pure.py` directly:
  merged solve_bubble_at_t + solve_dew_at_t into one function,
  solve_pseudo_pure_at_t(d1,d2,t_k,z1,rho_l0,rho_v0) -- 2 unknowns (u0,u1
  decode to rho_l,rho_v via the existing sigmoid mapping), 2 residuals
  (r1=P equality, r2=g equality via st_l.g_jmol vs st_v.g_jmol). Updated
  run_true_vle_envelope to call this once per T and append the SAME row to
  both bubble_rows and dew_rows (h_l_Jmol = bubble point, h_v_Jmol = dew
  point, from one solve -- no line lost). Also updated the plot title and
  module docstring/header (previously stale "true VLE" language) and
  bumped version to v0.3.0-pseudopure.
- Verified via py_compile (clean) and a full CLI smoke run
  (--w1 0.911 --Tmin 255.4 --Tmax 349.8 --n 80): compiled and ran cleanly.
  THIS is the run whose output (rho_l~603.7, rho_v~602.6 mol/m3 at
  T=255.4K -- nearly identical) first surfaced the degenerate-root bug
  documented in the "MAJOR FIX" entry above. The pseudo-pure conversion
  itself (2-unknown/2-residual structure) is believed correct; the bad
  output is inherited from the same underlying solver/gate bug as the
  true-VLE file's bubble branch, not a mistake in the conversion logic.
- File sent back to device and committed (mtime-guarded) BEFORE the
  degenerate-root bug was discovered -- i.e. the version currently on disk
  still has the bug. NOT yet re-fixed with the RHO_SEPARATION_MIN_RATIO +
  seed-ordering changes applied to mixture_dome_validation.py. This file
  has no separate dew branch to borrow a corrective seed from (bubble and
  dew are literally the same function call here) -- fixing it will need a
  different approach, e.g. multi-start across a wide rho_v range and
  picking the most-separated converged candidate, rather than
  cross-seeding from a sibling branch. Not yet started.

## Session Breadcrumb — 2026-08-13: Pseudo-Pure Density Bug Fix Strategy

### Root Cause Analysis (Stage 5 walkthrough)
Degenerate-root bug in `mixture_pseudo_dome.py` occurs because:
- Convergence residuals use relative differences (normalized by pressure/energy scale)
- Two nearly-identical densities (e.g., 603 vs 602 mol/m³) trivially satisfy relative-difference equations
- **Two mathematical solutions exist near any starting guess**: fake (rho_l ≈ rho_v, physically wrong) and real (rho_l >> rho_v, true phase equilibrium)
- Solver method `hybr` with current settings accepts the fake solution (settles early on nearby good residuals)
- Continuation (carry forward previous T's answer) locks in fake solution for all subsequent temperatures
- Result: 55-65% liquid-density error (~65 kg/m³ vs real ~1180 kg/m³)

### Fix Strategy (Three-Part)
1. **Tighten density check** (line 143): change `rho_l > rho_v * (1+1e-8)` to `rho_l > rho_v * 3.0`
   - Real equilibrium has 30-100x separation; fake has only 1.0x
   - 3x is conservative, will not reject physical solutions
   - Stops accepting degenerate roots

2. **Solver method change** (line 207): swap `method="hybr"` to `method="lm"`
   - LM is damped least-squares, less aggressive, avoids early acceptance of fake answer
   - Targets low-T divergence issue mentioned in earlier work

3. **Optional: tighter residual gates** (line 213)
   - Current: `r_p <= 1e-6 and r_mu <= 1e-6`
   - Consider tightening if needed after first run

### Expected Outcome
- Target: 140/140 converged points (full sweep)
- Check: rho_l values smooth across T range (no jumps from fake-to-real transitions)
- Then: compare P-T output vs Honeywell table
- Then: investigate remaining ~2% bias vs real 1180 kg/m³ value

### Implementation Status
- Changes keyed in: density check (3x), solver method (hybr→lm)
- First test run: check convergence count and rho_l smoothness

### CORRECTION: Changes were keyed into the wrong file
- The 3x density-check edit was first mistakenly applied to
  `mixture_pseudo_dome.py` (single-density-pair pseudo-pure CLI), then
  reverted per user instruction. hybr->lm solver-method edit was also
  applied and reverted to that same wrong file (test run failed anyway:
  sandbox has no `scipy` installed, so it was never actually executed).
- Correct target file confirmed by user: `mixture_dome_validation_pseudo_pure.py`
  (the bubble/dew mu-equality VLE solver, NOT the simpler pseudo_dome.py).
  3x density-check fix (`rho_l > rho_v * 3.0`, replacing
  `rho_l > rho_v * (1.0 + 1e-8)`) applied at BOTH acceptance-gate sites in
  this file: `solve_bubble_at_t` (~line 1128) and `solve_dew_at_t`
  (~line 1292). This file uses `scipy.optimize.least_squares`
  (bubble profile: `LSQ_LOSS="cauchy"`; dew profile: `LSQ_LOSS_DEW="huber"`),
  not `scipy.optimize.root(method=...)`, so the earlier hybr->lm idea
  does not apply here at all -- that was specific to the (wrong) pseudo_dome.py
  file's different solver call.

### Honeywell Solstice N15 (R-515B) Datasheet — SI Reference Values (source: user-uploaded TDS, 2026 edition, pages 1-3)
- Composition: R-1234ze(E)/R-227ea, 91.1%/8.9% mass -- matches w1=0.911 already used.
- Explicitly states "Zero glide" as a stated product property -- reconfirms azeotrope/pseudo-pure justification.
- Critical temperature: 228.0 F = 382.04 K (matches earlier Kay's-rule estimate).
- Critical pressure: 507 psig = 521.7 psia = 3596 kPa.
- Critical density: 31.03 lbm/ft3 = 497.1 kg/m3.
- Liquid density @ 77 F (298.15 K): 73.65 lbm/ft3 = 1180.0 kg/m3 (matches prior 1179.8-1179.9 kg/m3 reference used in Stage 5 true-VLE fix).
- Vapor density @ 77 F: 1.69 lbm/ft3 = 27.07 kg/m3.
- Saturated pressure @ 77 F: 57.45 psig = 72.146 psia = 497.5 kPa.
- Molecular weight: 117.5 lbm/lb-mol = 117.5 g/mol.
- Page 3 contains a dense P-T table, 0-170 F in 2 F steps (86 points) -- same table already used for the true-VLE Honeywell comparison; available again once pseudo-pure branch is validated.

### Bubble-branch CSV inspection (`r515b_true_vle_bubble.csv`, user-provided run) — kink root-caused
- Dew branch (`r515b_true_vle_dew.csv`): 140/140 CONVERGED, `rho_v_molm3` smooth and monotonic across the full 250-360 K sweep (44 -> 1298 mol/m3). No issue on dew side.
- Bubble branch: also reports 140/140 CONVERGED (residual gates ok), BUT `h_l_Jmol` is NOT monotonic -- decreases smoothly from 41316 J/mol (T=250K) down to a minimum ~14186 J/mol at T=278.49K, then reverses and climbs to 14628 J/mol by T=280.86K before continuing to rise. This turnover is the S-shaped kink visible in the user's plotted P-h envelope.
- Root cause confirmed NOT a residual-tolerance problem: `r_P`/`r_mu` at the turnover rows (37-41) are already 1e-11 to 1e-13, far tighter than the 1e-6 gate -- the solver is finding two different, both cleanly-converged, physically-inconsistent branches, not failing to converge.
- Three candidate fixes evaluated (not yet implemented):
  1. Tighten residual gates further -- REJECTED, residuals already far below any plausible tighter gate on both branches; would not discriminate between them.
  2. Fix continuation seeding (e.g., extrapolate trend from last 2-3 converged points instead of raw carry-forward of the single previous point) -- assessed as the most likely real fix; targets the actual mechanism (seed lands in wrong solution basin near T~278K).
  3. Adjust solver loss/settings (e.g., bubble's `LSQ_LOSS="cauchy"` -> `huber`, matching dew's already-working profile) -- complementary, cheap to try alongside #2, not a standalone fix since basin selection depends on both seed and loss function together.
- Decision: implement seeding fix (#2) first; test loss-function change (#3) as a secondary tweak if turnover persists. Not yet implemented.
- Next validation step once fixed: re-run against Honeywell P-T table (page 3 of TDS, 86 points) and against the liquid-density spot check (1180.0 kg/m3 @ 298.15K) from this same datasheet.

### 2026-08-13 — Option #3 implemented: bubble loss function switched cauchy -> huber
- User chose to try #3 (solver loss/settings) before #2 (seeding fix).
- Edited `mixture_dome_validation_pseudo_pure.py` module-level constants (~line 197-198):
  `LSQ_LOSS = "cauchy"` -> `LSQ_LOSS = "huber"`
  `LSQ_F_SCALE = 10.0` -> `LSQ_F_SCALE = 3.0`
  New values now match the dew branch's already-working profile
  (`LSQ_LOSS_DEW="huber"`, `LSQ_F_SCALE_DEW=3.0`). All other bubble LSQ_*
  constants (DIFF_STEP, MAX_NFEV, METHOD, X_SCALE, FTOL, XTOL, GTOL)
  left unchanged.
- Rationale logged inline in code comment: cauchy's aggressive
  down-weighting of large residuals was hypothesized to let the bubble
  solver settle early on a nearby-but-wrong solution basin near T~278K,
  producing the non-monotonic h_l turnover (S-kink) confirmed in the
  prior CSV inspection entry above.
- NOT yet run/verified -- user needs to execute the CLI (sandbox lacks
  scipy) and report back whether the h_l turnover near T=278.49K is gone
  and whether 140/140 still converges.
- If turnover persists after this change, next step is Option #2
  (trend-extrapolated continuation seeding instead of raw single-point
  carry-forward), as originally recommended as the more fundamental fix.

### 2026-08-13 — VERIFIED: Option #3 (cauchy->huber) fixed the bubble-branch kink
- User ran the CLI and confirmed via plot: bubble line (liq) now smooth
  and monotonic across the full range, h_l ~175-355 kJ/kg; dew line (vap)
  also smooth, ~360-420 kJ/kg. No S-turn/kink, no loop. Both branches
  visually clean on the log-P vs h envelope plot.
- Option #2 (continuation-seeding fix) NOT needed -- #3 alone resolved it.
- Root cause confirmed in practice: the bubble branch's `cauchy` loss was
  sufficient by itself to cause the non-monotonic turnover; switching to
  `huber` (matching dew) eliminated it without any seeding changes.
- `mixture_dome_validation_pseudo_pure.py` current confirmed-good state:
  `LSQ_LOSS="huber"`, `LSQ_F_SCALE=3.0` (bubble), density-check gate
  `rho_l > rho_v * 3.0` (both branches). This combination is now the
  working baseline for this file.

### Open items (in order), post-kink-fix
1. Compare pseudo-pure P-T output against Honeywell datasheet P-T table
   (page 3 of TDS, 86 points, 0-170 F) -- not yet done for this
   (now-fixed) pseudo-pure file.
2. Re-check liquid density accuracy: does `rho_l` at T=298.15K
   (77 F) now land near the Honeywell reference 1180.0 kg/m3 (within the
   ~2% bias already characterized for the true-VLE file), or does the
   ~55-65% deep-liquid density gap (documented earlier in this file, see
   "R-515B Liquid-Branch Pressure Validation" and "Root cause nailed
   down" sections above) still apply to the pseudo-pure path specifically?
   Not yet checked.

### 2026-08-13 — Dome truncation at critical point identified and fixed (taper + overlay)
- User reported the fixed-kink plot still doesn't close: bubble line ends
  ~355 kJ/kg, dew line ends ~420 kJ/kg at similar P (~45 bar), with a gap
  between them instead of converging to one point.
- Root cause: the same 3x density-separation gate that fixed the S-kink
  is fundamentally incompatible with the critical point, where
  `rho_l -> rho_v` by definition. As T -> Tc, even the TRUE physical
  solution's density ratio shrinks below 3x and gets rejected/marked
  DIVERGED, truncating the dome before bubble/dew can meet. This is the
  tradeoff flagged earlier when the user asked "would we not lose real
  roots?" -- confirmed to bite right at the critical point specifically.
- User decision: taper the gate (not just relax it uniformly), AND
  overlay/mark the near-critical segment distinctly in the plot rather
  than presenting it with the same confidence as the rest of the dome.
- Implemented in `mixture_dome_validation_pseudo_pure.py`:
  1. New module-level constants (~line 210 area): `RHO_SEP_RATIO_FAR=3.0`,
     `RHO_SEP_RATIO_NEAR_TC=1.05`, `RHO_SEP_TAPER_START_K=30.0`.
  2. New helper `_rho_separation_min_ratio(t_k, tred_mix)`: linear taper
     from 3.0 (at/below `Tred_mix - 30K`) down to 1.05 (at/above
     `Tred_mix`). Floor stays > 1.0 so exact rho_l==rho_v degeneracy is
     still rejected even at Tc.
  3. `tred_mix` computed once at the top of both `solve_bubble_at_t` and
     `solve_dew_at_t` via the existing `bell2023_Tred_vred(z1, z2, tc1,
     tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)` call (same mixture reducing
     temperature already used inside `mix_state`) -- NOT a hardcoded
     literature value, so the taper stays self-consistent with the
     model's own EOS at whatever composition is passed in. This is also
     the Tred already confirmed elsewhere in this file's history to match
     the Honeywell datasheet's Tc=382.04K almost exactly at w1=0.911.
  4. Acceptance gate in both `_attempt()` closures changed from
     `rho_l > rho_v * 3.0` to `rho_l > rho_v * min_ratio` where
     `min_ratio = _rho_separation_min_ratio(t_k, tred_mix)`.
  5. Each returned row dict gains two new fields: `rho_sep_min_ratio`
     (the actual ratio required/used at that T) and `near_critical`
     (bool, True once inside the 30K taper window) -- lets downstream
     code (plotting, CSV consumers) distinguish these points without
     re-deriving the taper logic.
  6. `plot_envelope()` updated to split each branch into far-field
     (solid line) and near-critical (dashed, alpha=0.6, separate legend
     entry: "Bubble/Dew line (near-critical, tapered gate)") segments,
     directly implementing the user's "produce that part as an overlay"
     request -- the near-critical closing segment is now visually
     flagged as solved-under-relaxed-confidence rather than drawn
     identically to the rest of the dome.
- Verified: `python3 -m py_compile mixture_dome_validation_pseudo_pure.py`
  passes cleanly. NOT yet run end-to-end (sandbox lacks scipy) -- user
  needs to execute the CLI and confirm (a) the dome now closes to a
  single point near Tred_mix, (b) the near-critical segment is visually
  distinct (dashed) in the output figure, (c) full sweep convergence
  count.

### How the three taper constants were chosen (2026-08-13) -- explicitly judgment calls, NOT data-derived
User asked directly how the taper numbers were picked. Honest accounting,
logged since these are starting values that may need retuning once real
near-critical output exists:
- `RHO_SEP_RATIO_FAR=3.0`: NOT new -- reused from the earlier
  already-validated conservative floor (real separation 30-100x vs fake
  degenerate root ~1.0x; discussed and agreed in the "would we not lose
  real roots?" exchange above).
- `RHO_SEP_RATIO_NEAR_TC=1.05`: a judgment call. Must satisfy two
  constraints -- strictly >1.0 (else an exact rho_l==rho_v degenerate
  root would pass even at Tc, defeating the gate's purpose) but close to
  1.0 (so the TRUE near-critical solution, whose own ratio approaches
  1.0 as T->Tc, isn't rejected). 5% above exact equality was picked as a
  round, defensible buffer -- NOT derived from actual converged
  near-critical density data, which was not available (no run of this
  file's near-Tc region has been done yet).
- `RHO_SEP_TAPER_START_K=30.0`: also a judgment call. No literature
  critical-scaling data or converged near-critical output was consulted
  to determine exactly how far below Tred_mix=382K the true physical
  ratio actually drops below 3.0x. 30K (~23% of the 250-380K sweep
  range) was chosen as a plausible width without empirical verification.
- Flagged as the first things to retune if the next run's dome still
  doesn't fully close or the near-critical segment looks numerically
  noisy/unstable: widen `RHO_SEP_TAPER_START_K`, and/or move
  `RHO_SEP_RATIO_NEAR_TC` closer to 1.0.

### 2026-08-13 — Taper confirmed working; remaining gap traced to Tmax, not the gate
- User ran with the taper (`--Tmax 380`, `--w1 0.911`) and reported "kind
  of worked, but not quite closed" -- plot showed dashed near-critical
  overlay on both branches, but bubble stopped ~375 kJ/kg and dew
  stopped ~410 kJ/kg, still a ~35 kJ/kg gap between them.
- Inspected user-uploaded `r515b_true_vle_bubble.csv` and
  `r515b_true_vle_dew.csv` tails (last rows, T up to 380.0K):
  - BOTH branches: 140/140 CONVERGED, reaching all the way to T=380.0K
    (the requested Tmax) with no premature DIVERGED truncation. The
    taper itself is not cutting the sweep short.
  - Densities ARE converging toward each other across branches as
    designed: at T=380.0K, bubble gives (rho_l=5345.80, rho_v=2864.24
    mol/m3, ratio=1.867) vs dew gives (rho_l=5340.30, rho_v=2865.30
    mol/m3, ratio=1.864) -- nearly identical despite being two
    independent equation systems (bubble solves y1 free with x1=z1
    fixed; dew solves x1 free with y1=z1 fixed). This confirms the taper
    is doing its job: both branches are genuinely approaching the same
    physical critical state, not just being permitted to converge to
    unrelated nearby roots.
  - Root cause of the visible gap: Tmax=380K sits ~2.04K short of
    Tred_mix=382.04K (the mixture reducing/critical-temperature proxy
    used by the taper itself). Converting the final rows to kJ/kg
    (mw_mix~117.48 g/mol): bubble h_l(380K)=375.95 kJ/kg, dew
    h_v(380K)=410.70 kJ/kg -- matches the ~35 kJ/kg gap seen in the
    plot. P is also climbing steeply near the end (3.45 MPa at 380K vs
    Pc=3.596 MPa), so the P-h curve is nearly vertical in this last
    stretch -- a small remaining T gap produces a large visual h gap
    purely from curve steepness, not from any solver/gate failure.
- Conclusion: the taper mechanism itself is validated and working
  correctly (converging densities confirm this). The dome's apparent
  non-closure is fully explained by Tmax not being pushed close enough
  to Tred_mix, not by any remaining flaw in `_rho_separation_min_ratio`
  or the gate logic.
- Recommended next run: `--Tmax 381.95 --n 200` (denser sampling near
  the steep near-critical region) to pull both branches' endpoints much
  closer together visually. Exact closure to a single point is not
  expected -- `RHO_SEP_RATIO_NEAR_TC=1.05` still requires nonzero
  separation, so the solver can approach but never solve exactly at
  Tred_mix (consistent with the true critical point being a genuine
  mathematical singularity, rho_l=rho_v exactly). Not yet run by user.

### 2026-08-13 — User tried Tmax=390 (above Tc), then settled on Tmax=382
- User first tried `--Tmax 390`, reported it "did not fix" the closure.
  Flagged to user: 390K is ABOVE Tred_mix=382.04K -- there is no
  liquid/vapor phase split above the critical temperature at all (single
  supercritical phase only), so pushing Tmax higher was asking the
  solver to find an equilibrium that does not physically exist. Expected
  outcomes there are either correct DIVERGED status past ~382K, or (red
  flag) spurious convergence to a non-physical root, since
  `RHO_SEP_RATIO_NEAR_TC=1.05` is a flat floor for any T>=Tred_mix and
  doesn't get looser above Tc -- but a coincidental root satisfying it
  without being physically meaningful was flagged as possible. Full
  CSV/plot from the T=390 run was not reviewed (user moved on before
  sharing it).
- User decision: cap `--Tmax 382` instead (0.04K short of
  Tred_mix=382.04K) -- essentially as close to true Tc as the model's
  own reducing-temperature estimate allows.
- Expectation set for this run: even at 382K, the very last few points
  approaching Tred_mix may legitimately DIVERGE, since the TRUE physical
  rho_l/rho_v ratio there is expected to already be below the
  `RHO_SEP_RATIO_NEAR_TC=1.05` floor (real ratio -> 1.0 exactly at Tc).
  This is expected/correct behavior, not a bug -- the dome should look
  nearly closed but may not include the literal last handful of
  temperature steps right at the edge. Recommended pairing with a denser
  `--n` (e.g. 200) given how steep the P-h curve is in this region.
  Command not yet run by user for this Tmax value.

### 2026-08-13 — Rigorous critical-point solve implemented (user chose "solve the true critical point" over cosmetic closure)
- User tried Tmax=390 (above Tc, flagged as unphysical -- no phase split
  exists above critical), then Tmax=382 (0.04K short), then reported the
  SAME plot at Tmax=382.4 (0.36K ABOVE Tred_mix=382.04K, technically
  already supercritical) -- gap still visible in all cases, confirming
  no amount of Tmax-tuning alone can close it: bubble/dew are two
  independently-solved equation systems that can only coincide exactly
  AT the true critical point, and `RHO_SEP_RATIO_NEAR_TC=1.05` (a strict
  >1.0 floor by design) means the solver can never land exactly there.
- Asked user to choose: (a) cosmetic closure -- draw a dotted line
  connecting the last converged bubble/dew points, clearly labeled as
  extrapolated/not solver-verified, OR (b) solve the actual mixture
  critical point via the EOS's own criticality conditions and use it as
  a genuine shared closing vertex. User chose (b).
- Implemented in `mixture_dome_validation_pseudo_pure.py`:
  1. `_pressure_rho_derivatives_fd(d1,d2,t_k,rho_mol,x1,rel_step=1e-4)`:
     new helper, central finite-difference (dP/drho)_T and
     (d2P/drho2)_T at fixed (T,x1), built on top of the already-existing
     `mix_state` (no analytic third-derivative-of-alphar machinery
     exists in this file, so FD on P itself is the simplest correct
     path).
  2. `solve_mixture_critical_point(d1,d2,z1,t_guess,rho_guess)`: new
     function, 2-unknown/2-equation `least_squares` solve for
     (Tc_mix, rhoc_mix) satisfying (dP/drho)_T=0 and (d2P/drho2)_T=0 at
     FIXED composition z1 (consistent with this file's pseudo-pure
     treatment -- standard fixed-N criticality condition, not full
     multicomponent critical-locus tracing). Residuals nondimensionalized
     by fixed p_scale/rho_scale constants (from the initial guess) so
     both equations are comparable magnitude despite different physical
     units. New dedicated tuning constants: `LSQ_METHOD_CRIT="trf"`,
     `LSQ_MAX_NFEV_CRIT=400`, `LSQ_FTOL/XTOL/GTOL_CRIT=1e-13`,
     `CRIT_RESID_GATE=1e-4` (looser than VLE's 1e-6 -- FD-based second
     derivative is inherently noisier). Initial guess: Tred_mix and
     1/vred_mix from the EXISTING `bell2023_Tred_vred` call -- i.e. the
     same mixing-rule reducing point already used for the taper, since
     it's already confirmed close to the real Honeywell Tc (382.04K).
  3. `run_true_vle_envelope` now computes Tred_mix/vred_mix once up
     front (composition-only, doesn't depend on the T-sweep), calls
     `solve_mixture_critical_point` once, and returns a 4th output:
     `crit_point` (dict: T_K, rho_molm3, P_Pa, h_Jmol, x1, converged,
     dp_drho_resid, d2p_drho2_resid, iterations, notes). Signature change
     -- callers updated (`_cli()`).
  4. `plot_envelope` gains new optional param `crit_point`. If provided
     and `crit_point["converged"]` is True, its (h,P) coordinate is
     appended as the shared FINAL vertex of both the bubble-near and
     dew-near dashed arrays (so both dashed lines terminate at the
     literal same point -- genuine closure), plus a separate black-star
     marker with its own legend entry "Critical point (solved)". If not
     converged/None, dome is left exactly as before -- no silent
     fabricated closure.
  5. `_cli()` updated: unpacks 4-tuple from `run_true_vle_envelope`,
     passes `crit_point` into `plot_envelope`, and records the full
     `crit_point` dict under a new `"critical_point"` key in the
     metadata JSON output.
  6. `save_csv`/CSV schema for bubble_rows/dew_rows themselves
     UNCHANGED -- crit_point is kept as a separate structure (not
     smuggled into the per-branch CSVs with a fake status string),
     avoiding any DictWriter/schema fragility.
- Verified: `python3 -m py_compile mixture_dome_validation_pseudo_pure.py`
  passes cleanly. NOT yet run end-to-end (sandbox lacks scipy) -- user
  needs to execute the CLI and confirm (a) `critical_point.converged` is
  True in the output metadata JSON, (b) the dome now visually closes to
  the black-star point with both dashed lines meeting it, (c) the solved
  Tc_mix is close to Tred_mix=382.04K / the Honeywell reference 382.04K
  as a sanity check.

### 2026-08-13 — Nested-Brent solve run #1 also FAILED: window mis-centered
- User ran the new nested-bracketed-Brent solve (--Tmin 250 --Tmax 380
  --n 140). Result: `bubble_converged=140/140`, `dew_converged=140/140`
  (full sweep clean), but `critical_point.converged: false`:
  `T_K=377.82`, `rho_molm3=2959.94`, `dp_drho_resid=9.37e-13` (Stage 1 --
  the inner spinodal-density bracket -- nailed a genuine root to machine
  precision), `d2p_drho2_resid=1.1528` (~11,500x over the 1e-4 gate,
  WORSE than the earlier flat-solve attempt's 0.1175).
- Diagnosed: rho_molm3=2959.94 sits almost exactly at the inner scan
  window's LOWER edge (0.7*rho_guess=0.7*4228=2959.6, matching to 4
  significant figures) -- a boundary artifact, not a genuine interior
  root. Confirmed by cross-referencing the SAME run's own converged dew
  branch: dew's vapor density at T=380K is 2864.2 mol/m^3, which is
  BELOW the search window's lower bound of 2960 -- i.e. the window was
  clipping off part of the very branch it needed to resolve, before it
  could ever reach the true (lower, still-separating) vapor-side
  spinodal density.
- Also confirmed from this run's own data that the true critical
  temperature MUST be above 380K: bubble/dew both still converge to
  clearly-separated real phases at T=380K (rho_l=5345.8 vs rho_v=2864.2,
  ratio~1.87, comfortably above the tapered gate's floor there) -- so
  the previous run's T_c=377.82 (below 380K) was independently confirmed
  wrong, not just numerically suspicious.
- Fix implemented in `mixture_dome_validation_pseudo_pure.py`:
  1. `solve_mixture_critical_point` signature extended with four new
     OPTIONAL override parameters: `rho_scan_lo`, `rho_scan_hi`,
     `t_scan_lo`, `t_scan_hi` (absolute K/mol-m^3 values, not
     offsets/factors). When omitted, falls back to the original
     `CRIT_*_SCAN_*`-factor-based windows around `t_guess`/`rho_guess`
     exactly as before (so the function remains usable standalone with
     just a guess, e.g. as a fallback when no sweep data exists).
  2. `run_true_vle_envelope` restructured: the critical-point solve now
     runs AFTER the bubble/dew T-sweep loop (previously ran BEFORE it,
     using only the Bell-mixing-rule guess). Tred_mix/vred_mix are still
     computed early (needed as the fallback guess/upper T bound), but
     `solve_mixture_critical_point` itself is deferred.
  3. After the sweep, finds the last CONVERGED bubble row and last
     CONVERGED dew row (via `next(... for r in reversed(...))`). If both
     exist: builds `crit_rho_guess` = midpoint of the last converged
     bubble `rho_l` and dew `rho_v` (physically grounded in real solved
     EOS behavior, not just the mixing-rule reducing density); builds
     explicit window overrides `rho_scan_lo=0.85*rho_v_last`,
     `rho_scan_hi=1.15*rho_l_last` (margins outside the real equilibrium
     densities, since spinodal points lie INSIDE the equilibrium
     bubble/dew envelope, between rho_v_last and rho_l_last, not outside
     it -- so a window with modest margin beyond each side safely
     contains both spinodal branches without excluding either);
     `t_scan_lo=t_last-2.0` (small margin below the highest T the sweep
     itself reached), `t_scan_hi=Tred_mix+CRIT_T_SCAN_HI_K` (unchanged
     upper bound, tied to the mixing-rule estimate since the sweep never
     goes supercritical). If no usable converged sweep data exists
     (e.g. Tmax set far below critical), falls back to the original
     Tred_mix/vred_mix-only guess with default windows.
- Verified: `python3 -m py_compile mixture_dome_validation_pseudo_pure.py`
  passes cleanly. NOT yet run end-to-end (sandbox lacks scipy) -- user
  needs to re-run and check: `critical_point.converged` should now be
  True, `T_K` should land ABOVE 380K (confirmed necessary from this
  run's own bubble/dew data), and both residuals should be comfortably
  under 1e-4.

### 2026-08-13 — Nested-Brent run #2 also FAILED: window correct, resolution too coarse
- User re-ran with the sweep-data-derived window fix. Result: window WAS
  correctly placed this time (notes confirmed rho window
  [2435.5,6147.7] mol/m^3, T window [378.000,384.522] K -- matches the
  intended design), but STILL `converged: false` with
  `"no outer sign-change bracket found"` -- the fallback path fired
  again, this time for a genuinely different reason than run #1
  (mis-centered window is now ruled out).
- Diagnosed: the "hump" between the two spinodal densities (where
  dP/drho<0) narrows toward EXACTLY zero width as T approaches the true
  critical temperature, by definition of criticality. The scan
  resolution at the time (`CRIT_RHO_SCAN_N=40` inner, `CRIT_T_SCAN_N=24`
  outer) was likely too coarse to ever sample a point actually falling
  inside this narrowing feature near the top of the search window --
  i.e. dP/drho may never have registered negative at any scanned point
  close to Tc, so `_find_sign_change_bracket` had no sign flip to find
  even though the window correctly contained the true crossing.
- Fix: raised `CRIT_RHO_SCAN_N` 40->150 and `CRIT_T_SCAN_N` 24->30
  (~4-5x more evaluations). Since this solve runs ONCE per CLI
  invocation (not per sweep temperature point), the added cost is a
  flat one-time runtime increase, not a per-point multiplier -- traded
  for the ability to resolve a narrower hump closer to Tc. Flagged to
  user that the run will take noticeably longer this time.
- Verified: `python3 -m py_compile mixture_dome_validation_pseudo_pure.py`
  passes cleanly. NOT yet run end-to-end. If STILL no bracket found
  after this resolution increase, next diagnostic step would be to log
  the actual dP/drho values at each outer-scanned T (not just
  pass/fail) to see how close to zero the minimum |dP/drho| gets even
  without a clean sign flip -- would help distinguish "still too coarse"
  from "a deeper issue in the window/algorithm design."

### 2026-08-13 (earlier same day) — First critical-point solve attempt FAILED; diagnosed and rewritten using literature method
- User ran the flat-2D-least_squares critical-point solve. Metadata
  showed `converged: false`: T_K=381.52 (plausible location), but
  `dp_drho_resid=0.00088` (9x over the 1e-4 gate) and
  `d2p_drho2_resid=0.1175` (~1175x over gate). `notes` said "`xtol`
  termination condition is satisfied" -- diagnosed as the classic
  symptom of a solver STALLING (step size vanished) rather than actually
  reaching the root, which is expected behavior at a genuine critical
  point: the criticality conditions' Jacobian is singular there by
  definition (that's literally what makes it a "critical"/degenerate
  point), so a flat simultaneous 2D solver structurally loses traction
  exactly where it needs precision most.
- User uploaded `1-s2.0-S0378381216305349-main.pdf`: Bell, I.H., Jager,
  A. (2017), "Calculation of critical points from Helmholtz-energy-
  explicit mixture models," Fluid Phase Equilibria 433, 159-173,
  https://doi.org/10.1016/j.fluid.2016.10.030 -- same author (Ian H.
  Bell) as the Bell (2023) mixing-rule paper already used throughout
  this project. Read via `pdftotext -layout` (pdftoppm/poppler-utils
  page-render was unavailable in sandbox; text extraction worked fine).
- Key findings applied from the paper:
  1. Confirms the general multicomponent critical-point conditions are
     Legendre-transform matrix determinants L1=det(L*)=0 (first
     criticality condition) and M1=det(M*)=0 (second), built from
     (1/RT)*(d2A/dni/dnj)_{T,V} entries -- BUT for a FIXED-composition
     (N=1-effective) system, exactly this project's pseudo-pure/azeotrope
     treatment, these reduce to the classical single-component
     conditions (dP/dV)_T=0 and (d2P/dV2)_T=0 already being used here.
     Confirms the PHYSICS/equations in this file were already correct;
     only the numerical method needed fixing.
  2. Section 3.1 explicitly states Newton-Raphson "requires... quite
     good estimates" near critical, and that first-order (Newton)
     methods need a much better starting guess than second-order
     (Halley's) methods -- directly explains why the flat least_squares
     attempt stalled even with a good (Tred_mix-based) starting guess.
  3. Section 1.2's literature review cites Hoteit et al.'s approach
     "based on nested and bounded iterations of Brent's method,"
     reporting it "more than a thousand times faster" than an
     alternative while retaining similar reliability -- validates
     nested bracketed root-finding (not flat least_squares, not
     unconstrained minimization) as the literature-endorsed practical
     method for exactly this problem.
- Rewrote `solve_mixture_critical_point` in
  `mixture_dome_validation_pseudo_pure.py` using a NESTED BRACKETED
  BRENT'S METHOD design (fixed-composition reduction of the paper's
  general contour-tracing strategy):
  1. New `_dp_drho_fd(d1,d2,t_k,rho_mol,x1,rel_step)`: lightweight
     2-point central-difference first derivative only (cheaper than the
     existing 3-point `_pressure_rho_derivatives_fd`, since this is
     called many times inside the inner scan/bracket search).
  2. New `_find_sign_change_bracket(func,x_lo,x_hi,n_scan)`: generic
     scan-then-bracket helper -- scans n_scan points, returns the first
     verified sign-change bracket (or None), required because
     `scipy.optimize.brentq` needs a pre-verified bracket and does not
     search for one itself. Non-finite/exception points are skipped as
     gaps rather than treated as false sign changes.
  3. New `_spinodal_rho_at_t(d1,d2,t_k,x1,rho_scan_lo,rho_scan_hi,...)`:
     Stage 1 -- at FIXED T, brackets and roots (dP/drho)_T=0 in rho
     alone via `_find_sign_change_bracket` + `brentq`. Returns
     (None, None) if no bracket found (e.g. T is already supercritical
     for this density window) rather than raising, so the outer stage
     can treat it as a scan gap.
  4. New module constants: `CRIT_RHO_SCAN_LO_FACTOR=0.7`,
     `CRIT_RHO_SCAN_HI_FACTOR=1.6` (inner rho window, narrow since we
     already have a good density guess -- unlike the paper's from-scratch
     0.5x-to-wide search), `CRIT_RHO_SCAN_N=40`, `CRIT_T_SCAN_LO_K=6.0`,
     `CRIT_T_SCAN_HI_K=3.0` (outer T window around t_guess),
     `CRIT_T_SCAN_N=24`. `CRIT_RESID_GATE=1e-4` retained from the prior
     attempt.
  5. `solve_mixture_critical_point` REWRITTEN (same external signature
     and return dict schema, so `run_true_vle_envelope`/`plot_envelope`/
     `_cli()` wiring from the prior entry needs no changes): outer
     function `g(T)` = (d2P/drho2)_T evaluated at Stage 1's spinodal
     density for that T (raises if Stage 1 can't bracket, caught by the
     outer `_find_sign_change_bracket`'s scan). Finds outer sign-change
     bracket in T, refines via `brentq`, then re-derives rho_c from the
     final T_c via Stage 1 one more time. If no outer bracket is found,
     returns `converged=False` with a diagnostic note identifying which
     search window failed, rather than raising or silently returning a
     bad point. Old flat-least_squares constants (`LSQ_*_CRIT`) removed
     entirely (structurally replaced, not just superseded).
  6. Added `brentq` to the `scipy.optimize` import line (alongside
     existing `least_squares, root`).
- Verified: `python3 -m py_compile mixture_dome_validation_pseudo_pure.py`
  passes cleanly. NOT yet run end-to-end (sandbox lacks scipy) -- user
  needs to re-run the CLI and check the new `critical_point` block in
  metadata JSON: expect `converged: true` this time, with both residuals
  comfortably under `1e-4`, and `T_K` close to Tred_mix=382.04K /
  Honeywell's 382.04K.

### 2026-08-13 — THIRD design: bisection on spinodal-dip existence (avoids noisy d2P/drho2 entirely)
- Bumping scan resolution 4-5x (run #2 fix) made NO difference --
  identical failure, identical window. Ruled out "too coarse" as the
  explanation. Wrote a standalone diagnostic script
  (`diagnose_spinodal.py`, written directly into the repo) to print raw
  (dP/drho)_T values across the real search window at 11 temperatures,
  citing Bell & Jager (2017)'s explicit warning that "multi-fluid model"
  criticality contours "can be not smooth" for some mixture/EOS
  combinations as the reason to get real data before guessing again.
- User ran it. Results were highly informative:
  - The spinodal "dip" (region where dP/drho<0) DOES exist and shrinks
    CLEANLY and monotonically from T=376K (min dP/drho=-89.1) down to
    T=381.5K (min=-5.53), then is GONE entirely by T=382.0K (min=+1.52,
    fully positive) -- exactly the expected, well-behaved physical
    transition. The true Tc sits between 381.5K and 382.0K.
  - BUT the d2P/drho2 values printed alongside at each T's dip minimum
    were NOT monotonically shrinking toward zero as expected -- they
    flip sign essentially at random: -5.9e-3, -7.3e-3, +4.0e-3, -8.1e-4,
    +1.9e-3, -4.2e-3, -1.1e-3 across consecutive ~0.5-1K steps from
    376K to 381.5K. Diagnosed as finite-difference NOISE dominating a
    genuinely-shrinking-toward-zero physical signal: the true curvature
    IS getting small near critical (as expected), but the FD noise
    floor does not shrink with it, so noise wins exactly where the
    second attempt's design needed precision most. This explains BOTH
    prior failures: the outer search was root-finding on a quantity
    that was never going to cross zero cleanly.
- THIRD redesign of `solve_mixture_critical_point` in
  `mixture_dome_validation_pseudo_pure.py`, avoiding the second
  derivative entirely:
  1. New `_find_all_sign_changes(func,x_lo,x_hi,n_scan)`: like
     `_find_sign_change_bracket` but returns EVERY bracket found (not
     just the first) -- needed to locate BOTH spinodal roots (vapor-side
     AND liquid-side) from one scan.
  2. New `_spinodal_pair_at_t(d1,d2,t_k,x1,rho_scan_lo,rho_scan_hi,...)`:
     Stage 1 -- finds BOTH spinodal densities (rho_v_sp, rho_l_sp) at
     fixed T using only first-derivative brackets (already proven
     reliable to ~1e-13). Returns None if fewer than 2 sign changes
     found (no dip exists -- T is at/above Tc).
  3. Stage 2 completely redesigned: instead of root-finding a noisy
     quantity, BISECT on T using dip-existence itself (a robust,
     discrete, always-well-defined test: `_spinodal_pair_at_t` returns
     a pair or None, nothing in between to be noisy about). Explicit
     preconditions checked before bisecting: dip MUST exist at
     `t_scan_lo` and MUST NOT exist at `t_scan_hi` (returns
     `converged=False` with a specific diagnostic note identifying which
     precondition failed, rather than bisecting on a false assumption).
  4. New constants: `CRIT_BISECT_MAX_ITERS=50`, `CRIT_BISECT_T_TOL_K=1e-4`
     (bisection stops once the T-bracket is this narrow -- reached in
     ~log2(9K/1e-4K)~=17 steps, well under the cap), `CRIT_WIDTH_FRAC_GATE=0.02`
     (converged also requires the final spinodal-density gap, as a
     fraction of rho_scale, to be below 2%).
  5. `converged` now determined ONLY by the T-bracket tolerance being
     reached AND the width-fraction gate -- `dp_drho_resid`/
     `d2p_drho2_resid` are still computed and reported in the output
     dict for diagnostic continuity, but are explicitly NOT part of the
     convergence decision anymore (documented inline as to why).
  6. Old `_spinodal_rho_at_t` (single-root Stage 1 from design #2) and
     the old `CRIT_RESID_GATE`/`CRIT_T_SCAN_N` outer-scan-count constant
     removed entirely (structurally superseded, not just unused).
     `_find_sign_change_bracket` (design #2's helper) kept as-is
     (harmless, still a valid generic utility, not currently called but
     left in case future single-root use arises).
- Verified: `python3 -m py_compile mixture_dome_validation_pseudo_pure.py`
  passes cleanly. `run_true_vle_envelope`'s call site is unaffected
  (same keyword-argument names: rho_scan_lo/hi, t_scan_lo/hi, t_guess,
  rho_guess). `plot_envelope`'s usage (`crit_point["converged"]`,
  `["h_Jmol"]`, `["P_Pa"]`) also unaffected -- schema-compatible. NOT yet
  run end-to-end. Expect: `converged: true`, `T_K` between 381.5 and
  382.0 (per the diagnostic script's own findings), `spinodal_width_molm3`
  very small, and a `notes` field showing the final T-bracket width.

### 2026-08-13 — CONFIRMED: bisection design (#3) converged, dome closes, matches Honeywell
- User initially pasted a STALE metadata.json (byte-identical
  `generated_at_utc` timestamp to a message from two turns earlier,
  missing the new `spinodal_width_molm3` field, old design-#2 notes
  text) -- flagged as a copy-paste mixup, not a regression, since the
  closed-dome plot already shown could only have come from a converged
  design-#3 run. Asked user to `cat` the file fresh.
- Fresh read confirmed SUCCESS:
  `T_K=381.8939`, `rho_molm3=3858.39`, `P_Pa=3576456.65` (3576.5 kPa),
  `h_Jmol=46587.87`, `converged: true`,
  `dp_drho_resid=1.64e-6`, `d2p_drho2_resid=6.41e-6` (both tiny --
  reported for diagnostic continuity only, not the actual convergence
  gate per the design-#3 rationale, but reassuringly small anyway),
  `spinodal_width_molm3=9.888` (0.241% of rho_scale, comfortably under
  the 2% `CRIT_WIDTH_FRAC_GATE`), `iterations=17` (matches the
  ~log2(9K/1e-4K)~=17 predicted bisection step count).
- Cross-checked against the Honeywell Solstice N15 TDS reference values
  (already logged earlier this session):
  - T_K=381.894 vs Honeywell Tc=382.04K -- **0.146K off (0.038%)**.
  - P_Pa=3576.5 kPa vs Honeywell Pc=3596 kPa (507 psig) -- **0.54% off**.
  Both comfortably within the ~2% extrapolation bias already
  characterized for this Bell-2023-based model elsewhere in this
  project -- strong independent validation that the solved critical
  point is physically correct, not just numerically self-consistent.
- Dome now visually closes: bubble and dew near-critical dashed
  segments both terminate at the identical (h,P) coordinate, marked with
  a black-star "Critical point (solved)" legend entry, confirmed via
  user-shared plot.
- STATUS: the "close the gap" task is COMPLETE. Three solver-design
  iterations were required (flat least_squares -> noisy
  d2P/drho2-based nested Brent -> width-based bisection on dip
  existence), each diagnosed empirically from real run output rather
  than guessed, with the final design validated both internally
  (residuals/width all comfortably converged) and externally (matches
  independent Honeywell reference data to within known model bias).

### 2026-08-13 — Conceptual aside: what is a spinodal (no code changes)
User asked for a plain explanation of "spinodal," used throughout the
critical-point work above. Answered conversationally, not logged in
detail here since no code/decisions changed -- summary for future
reference: the spinodal is where (dP/drho)_T=0, the mechanical-stability
limit of a single phase (distinct from the bubble/dew equilibrium
densities, which sit OUTSIDE the spinodal pair -- real phase separation
happens before the fluid ever reaches actual mechanical instability, so
ordering is vapor_eq < vapor_spinodal < liquid_spinodal < liquid_eq).
Two spinodal points (vapor-side, liquid-side) exist below Tc, merging
into one point (where d2P/drho2=0 too) exactly at the critical
temperature -- the physical basis for this file's Stage-1 (find
spinodal at fixed T)/Stage-2 (find T where spinodal's curvature -> 0)
nested search strategy.

### 2026-08-13 (later) — File-identity clarification + Honeywell 86-point P-T revalidation built

- Clarified three separate, similarly-named files so future sessions don't
  conflate them: `mixture_pseudo_dome.py` (old 2-unknown pseudo-pure file,
  the hybr->lm thread from before the bubble/dew architecture split --
  superseded, not an open bug); `mixture_dome_validation.py` (the
  user-authored true-VLE file, already independently fixed/validated to
  1.84% MAPE vs the same Honeywell 86-point table, ~2% liquid-density bias
  at 77F); and the file actively worked on this whole session, now saved
  by the user as `Honeywell_T_P_revalidation.py` (a renamed copy of
  `mixture_dome_validation_pseudo_pure.py`'s fully-fixed state -- huber
  loss, tapered rho-separation gate, bisection critical-point solver all
  already present). A pasted "1.01-1.1% average error" claim was traced to
  an older Codex-authored version predating today's fixes, not a
  verified number for the current file.
- User directed: build the Honeywell revalidation in
  `Honeywell_T_P_revalidation.py` directly (this file's own docstring
  already states this exact purpose: "if the dew point and bubble points
  are validated against Honeywell data, the code is assumed to be
  validated").
- Re-extracted the ground-truth 86-point P-T table (0-170F, 2F steps,
  psig) directly from page 3 of the Honeywell Solstice N15 TDS PDF via
  `pdftotext -layout` (not manual retyping/memory -- the sandbox still
  lacks `pdftoppm`/poppler for page-image rendering, same workaround as
  the Bell & Jager paper earlier this session). Verified programmatically:
  86 points, strictly monotonic 2F steps, 0 to 170F.
- Added to `Honeywell_T_P_revalidation.py`:
  - `HONEYWELL_TDS_PT_TABLE_F_PSIG`: the embedded 86-point reference table.
  - `_f_to_k`, `_psig_to_pa`: unit conversions (psig -> psia via +14.696,
    then exact psi->Pa factor 6894.757293168361).
  - `run_honeywell_pt_comparison(fluid1, fluid2, w1, out_csv,
    out_summary_json)`: calls the EXISTING `run_true_vle_envelope` with
    `t_vals` = the table's 86 T_K values directly (already ascending, so
    continuation seeding works unmodified), then compares model pressure
    to the chart pressure at each point. x1=y1 simplification implemented
    as: average bubble and dew P into one P_model per T (falling back to
    whichever single branch converged if only one did), justified by the
    datasheet's own "Zero glide" characterization of this blend (page 1)
    -- Honeywell itself publishes only one P-T curve, not separate
    bubble/dew curves, for this composition. Writes a per-point CSV and a
    summary JSON (MAPE_pct, bias_pct, max_abs_err_pct, convergence counts,
    critical_point passthrough).
  - Wired into `_cli()` via a new `--honeywell-compare` flag (plus
    `--honeywell-csv`/`--honeywell-summary` path args); when set, this
    replaces the linspace dome sweep entirely (Tmin/Tmax/n ignored).
  - Table's T range (0-170F = 255.37-349.82K) is comfortably subcritical
    relative to Tred_mix (~382K), so this comparison never touches the
    near-critical taper/bisection code paths -- no interaction risk with
    today's earlier critical-point fix.
- Verified via `python3 -m py_compile Honeywell_T_P_revalidation.py`
  (clean) and a standalone `ast.literal_eval` check confirming the
  embedded table has exactly 86 points in strict 2F steps. NOT yet run
  end-to-end (sandbox lacks scipy, per this project's standing
  limitation) -- next step: user runs
  `python Honeywell_T_P_revalidation.py --w1 0.911 --Tmin 250 --Tmax 380
  --n 5 --honeywell-compare` (Tmin/Tmax/n are ignored in this mode but
  currently still required by argparse) and reports back the printed
  MAPE/bias/max-error summary and/or the two output files
  (`verification/r515b_honeywell_pt_revalidation.csv`,
  `verification/r515b_honeywell_pt_revalidation_summary.json`).

### 2026-08-13 (later still) — CONFIRMED: Honeywell_T_P_revalidation.py validated, 86/86, MAPE 1.84%

- User ran `--honeywell-compare` and pasted back both the summary JSON and
  the full 86-row CSV. Results:
  - 86/86 points CONVERGED (0 diverged), all via `xtol`/`gtol` termination,
    retry=0 on every single point (first seed always sufficient -- no
    reliance on the 4-seed retry fallback anywhere across the whole range).
  - MAPE_pct = 1.8396%, bias_pct = 1.8396% -- EXACTLY equal, confirming
    every one of the 86 points has the same sign of error (model
    consistently OVERpredicts chart pressure, never under). This is a
    systematic bias, not scattered noise/cancellation.
  - Error is largest at the cold end: 4.19% at T=0F (chart P=0.7psig=
    ~106kPa absolute -- a small absolute Pa error is amplified in percent
    terms at such low absolute pressure), and shrinks monotonically down
    to 0.43% at T=170F (chart P=255.7psig, the table's warmest/highest-P
    point). No sign flips, no local error spikes -- clean monotonic decay.
  - Critical point (recomputed as part of this same run, via
    `run_true_vle_envelope`'s post-sweep solve): T_K=381.893 vs
    Honeywell's 382.04K (0.04% off), P_Pa=3,576,404 vs Honeywell's Pc=
    507psig=3,596,432Pa (0.56% off) -- consistent with the standalone
    critical-point check earlier this session.
  - Cross-check: this 1.84% MAPE is numerically indistinguishable from the
    separately-debugged `mixture_dome_validation.py`'s own 1.84% MAPE
    result against the SAME 86-point table (different file, different bug
    history, same underlying Bell-2023 EOS/composition) -- two
    independent codebases landing on the same figure is strong evidence
    this ~1.8% is a real characteristic of the model (the ~2%
    Bell-2023-composition-extrapolation bias theorized/documented
    repeatedly throughout this project), not a coincidental artifact of
    either file's specific fix history.
- STATUS: `Honeywell_T_P_revalidation.py` is now validated end-to-end
  against the Honeywell Solstice N15 TDS -- P-T curve (86/86 points,
  1.84% MAPE) AND critical point (Tc/Pc both within ~0.5%). This closes
  the "compare against Honeywell" task opened earlier this session.
  Output files: `verification/r515b_honeywell_pt_revalidation.csv`,
  `verification/r515b_honeywell_pt_revalidation_summary.json`.

### 2026-08-13 (milestone) — MILESTONE CLOSED: R-515B saturation curve (bubble/dew dome) validation

- Formally marking this milestone DONE. Scope covered: the true-VLE
  bubble/dew solve for R-515B (R-1234ze(E)/R-227ea, w1=0.911) in
  `Honeywell_T_P_revalidation.py`, validated two ways against the
  Honeywell Solstice N15 TDS:
  1. Full 86-point P-T table (0-170F, 2F steps): 86/86 converged,
     MAPE=1.84%, bias=1.84% (systematic, one-directional, monotonically
     shrinking from 4.19% at 0F to 0.43% at 170F).
  2. True critical point (rigorous bisection-on-spinodal-dip solve):
     Tc within 0.04%, Pc within 0.56% of Honeywell's published values.
  Both results cross-checked against the independently-debugged
  `mixture_dome_validation.py` (same 1.84% MAPE), and the residual ~1.8%
  error is understood as the model's known Bell-2023-composition-
  extrapolation bias, not an unresolved defect.
- OUT OF SCOPE for this milestone (unchanged from this file's own
  docstring, written back on 2026-03-02): interior two-phase
  (quality/lever-rule) states. The saturation dome (bubble line + dew
  line only) is what was validated; no quality-line calculation exists
  yet anywhere in this file.
- NEXT MILESTONE (opening now): quality lines (constant vapor-fraction x
  contours, e.g. x=0.1...0.9) inside the two-phase region, referencing
  the Honeywell TDS's "PRESSURE AND ENTHALPY" p-h chart on page 2. Not
  yet scoped -- see the next entry below once scoping is settled with
  the user (which file this lives in, and what "validated against the
  Honeywell chart" means given that chart is an image, not extractable
  text data like the page-3 P-T table was).

### 2026-08-13 (later) — New file mixture_quality_validation.py created; ported taper + critical-point solver from Honeywell_T_P_revalidation.py

- User created `mixture_quality_validation.py` as their own copy of
  `mixture_dome_validation.py`'s post-"MAJOR FIX" state (dew-first
  reordering, flat `RHO_SEPARATION_MIN_RATIO=3.0` gate, NO taper, NO
  critical-point solver -- i.e. the pre-taper state that hit 1.84% MAPE).
  This is the file for the upcoming quality-line work (implicitly answers
  the still-open "which file" scoping question with "a new file", per the
  earlier AskUserQuestion that got interrupted/declined).
- Two explicit edits requested and made:
  1. `plot_envelope`'s bubble/dew lines changed from default-color-cycle
     to solid black, `lw=2.0` (both `ax.plot` calls).
  2. Ported the FULL taper + true-critical-point solver stack from
     `Honeywell_T_P_revalidation.py` into this file verbatim (same
     physics/constants/functions: `RHO_SEP_RATIO_FAR/NEAR_TC/TAPER_START_K`,
     `_rho_separation_min_ratio`, `_pressure_rho_derivatives_fd`,
     `_dp_drho_fd`, `_find_sign_change_bracket`, `_find_all_sign_changes`,
     `_spinodal_pair_at_t`, `solve_mixture_critical_point`,
     `CRIT_*` constants), replacing the flat `RHO_SEPARATION_MIN_RATIO`
     gate in both `solve_bubble_at_t`/`solve_dew_at_t` (added `tred_mix`
     computation to each, mirroring the Honeywell file). `run_true_vle_envelope`
     now returns a 4-tuple `(bubble_rows, dew_rows, z1, crit_point)`,
     solving the critical point AFTER the sweep using real converged
     bubble/dew data to build the search window (same design as the
     Honeywell file). `_cli()` and metadata JSON updated to match.
  3. IMPORTANT DEVIATION from the Honeywell file, per explicit user
     instruction ("No demarcation between tapering and critical point"):
     `plot_envelope` here does NOT split converged points into
     far/near-critical subsets and does NOT use dashed/alpha styling for
     the near-critical segment. Both the tapered acceptance gate (which
     rows count as CONVERGED) and the critical-point closure (appending
     the solved (h,P) as each line's shared final vertex) are still
     present and functionally identical to the Honeywell file -- only
     the VISUAL demarcation between "far" and "near-critical" points was
     removed, per user preference for a single clean black line.
- Verified via `python3 -m py_compile mixture_quality_validation.py`
  (clean). NOT yet run end-to-end (sandbox lacks scipy) -- not yet
  re-validated against Honeywell after this port (the file's docstring
  still says its purpose is Honeywell-based validation of "vanilla
  code"; the port should reproduce the same ~1.84% MAPE / critical-point
  match as Honeywell_T_P_revalidation.py once run, since the underlying
  physics is identical, but this has not been explicitly re-confirmed on
  this specific file/copy).
- STATUS: user said they will make the next set of changes themselves
  (quality-line work). Task #3 ("Scope quality-line work") still open --
  scoping question (validation approach: lever-rule-only vs digitized
  chart-image overlay) was never answered since the AskUserQuestion call
  was interrupted; revisit if/when needed.

### 2026-08-13 (later still) — mixture_quality_validation.py abandoned; fresh mixture_quality_line_validation.py created as a direct copy of Honeywell_T_P_revalidation.py

- Root cause of user confusion this stretch: `mixture_quality_validation.py`
  (the file ported taper+critical-point INTO, above) was a copy of
  `mixture_dome_validation.py`'s state, which STILL had the pre-fix
  `LSQ_LOSS="cauchy"`/`LSQ_F_SCALE=10.0` bubble-branch config -- the exact
  S-kink bug fixed earlier this session, but only ever in
  `Honeywell_T_P_revalidation.py`/`mixture_dome_validation_pseudo_pure.py`,
  never in `mixture_dome_validation.py` itself. So the "not smooth" dome
  the user was seeing was a REAL, separate bug (stale loss config), not a
  cosmetic marker issue -- flagged this explicitly but user said "No" to
  that fix and instead chose to reset.
- User decided to abandon `mixture_quality_validation.py` and start over
  by copying `Honeywell_T_P_revalidation.py` directly (already has the
  correct `LSQ_LOSS="huber"` on BOTH branches, confirmed by grep after
  the copy) into a new file: `mixture_quality_line_validation.py`.
- Edited `mixture_quality_line_validation.py`'s `plot_envelope` (ported
  as-is from Honeywell_T_P_revalidation.py, including the far/near-critical
  dashed-overlay split and star marker at Tc): removed the far/near split
  entirely, removed the dashed near-critical styling, removed the
  `marker="*"` critical-point star. Both bubble and dew branches are now
  ONE continuous solid black line (`lw=2.0`), closing at the solved
  critical point coordinate (still appended as each line's shared final
  vertex -- only the visual demarcation/marker was removed, not the
  underlying closure logic). Verified via `python3 -m py_compile` (clean).
- STATUS: `mixture_quality_line_validation.py` is now the active file for
  quality-line work. Not yet run end-to-end (sandbox lacks scipy) --
  should reproduce Honeywell_T_P_revalidation.py's smooth dome (correct
  loss config already in place) with the requested minimal-marker
  presentation. `mixture_quality_validation.py` is considered abandoned/
  superseded as of this entry; do not continue building on it.

### 2026-08-13 (later still) — Honeywell p-h chart (page 2) provided; IP-unit conversion added as a final-step-only, validation-only path

- User shared the actual Honeywell Solstice N15 TDS p-h chart image
  (page 2, "PRESSURE AND ENTHALPY") -- this is the chart with the real
  printed quality lines (x=0.1...0.9) that the next milestone (quality
  lines) will eventually need to reference. Confirmed its axes are IP
  units: Pressure (psia, log scale, 15-1350), Enthalpy (Btu/lbm, 70-230).
  Footnote states its reference state explicitly: "h = 200 kJ/kg,
  s = 1.00 kJ/kg-K; sat. liq. at 0C" (IIR/ASHRAE-style convention) --
  flagged that this file's EOS-native enthalpy zero-point is NOT
  necessarily the same, so absolute h values may stay offset even after
  unit conversion; only shape/width comparisons are guaranteed valid
  until that's separately reconciled (not yet done).
- IMPORTANT ARCHITECTURE DECISION from user: "The conversion should be
  applied on final values only. We will use the code as is for DVCT, but
  validation should be against Honeywell." I.e. the core SI solve and
  the DVCT-facing `plot_envelope` must NEVER be touched by unit
  conversion -- IP conversion belongs in a separate, final-step-only,
  validation-only path.
- Implemented in `mixture_quality_line_validation.py`:
  - `BTU_LBM_TO_KJ_KG = 2.326` (exact) constant, next to the existing
    `PSI_TO_PA`.
  - `to_honeywell_units(h_kj_kg, p_pa) -> (h_btu_lbm, p_psia)`: pure
    final-value unit conversion, explicitly documented as never to be
    called inside the solve or inside `plot_envelope`.
  - `plot_envelope_honeywell_units(...)`: a validation-only TWIN of
    `plot_envelope` -- identical SI data prep, identical
    CONVERGED-only/critical-point-closure/single-black-line-no-marker
    presentation, but every (h,P) point is run through
    `to_honeywell_units` as the LAST step before plotting. `plot_envelope`
    itself is completely unmodified.
  - Wired into `_cli()`: new `--fig-honeywell-units` arg (default
    `verification/r515b_true_vle_envelope_honeywell_units.png`), called
    right after the existing SI `plot_envelope` call so both figures are
    always produced together from the same solve.
- Verified via `python3 -m py_compile mixture_quality_line_validation.py`
  (clean). NOT yet run end-to-end (sandbox lacks scipy).
- Command to run: `python mixture_quality_line_validation.py --w1 0.911
  --Tmin 250 --Tmax 380 --n 140` -- now produces BOTH
  `verification/r515b_true_vle_envelope.png` (SI, DVCT) and
  `verification/r515b_true_vle_envelope_honeywell_units.png` (Btu/lbm,
  psia, for visual comparison against the Honeywell TDS page-2 chart).

### 2026-08-13 (later still) — Quality lines (lever-rule interior states) added; user self-implemented plan, I found and fixed 5 real bugs

- Explained the underlying physics conversationally first (user asked
  "how are we computing everything"): bubble/dew lines are genuine EOS
  solves; quality lines are NOT a new solve at all, just the lever rule
  h(x) = h_l + x*(h_v - h_l) applied to already-solved bubble/dew
  enthalpies at each T, using P = avg(P_bubble, P_dew) as the shared
  pressure (same zero-glide simplification as run_honeywell_pt_comparison).
  User confirmed understanding ("Oh so we are not re-solving anything").
- Gave a full 6-change implementation plan for `compute_quality_lines`
  (new function), `QUALITY_LINE_VALUES` (new constant, x=0.1...0.9),
  and wiring `quality_lines` through both `plot_envelope` and
  `plot_envelope_honeywell_units` (thin gray lines, `lw=0.75`, labeled
  "x=0.1" etc., no change to the solver).
- User applied the plan themselves in `mixture_quality_line_validation.py`.
  User then asked me to "check the code now" -- `py_compile` passed but
  that only checks syntax, not runtime name resolution, so I read the
  diff by hand and found 5 real bugs that would have crashed or
  silently produced wrong output the moment the file was actually run:
  1. `qualities: List[float] = Quality_line_values` -- wrong-case name,
     and the constant was never defined under that name anywhere (would
     `NameError` at module-load time, since default-arg values are
     evaluated when Python reads the `def` statement).
  2. `quality_lines = compute_quality_lines(bubble_rows, dew_rows)`
     placed at MODULE level (right after `run_true_vle_envelope`'s
     `return`), not inside `_cli()` -- `bubble_rows`/`dew_rows` don't
     exist in module scope, and `compute_quality_lines` isn't even
     defined yet at that point in the file (declared further down).
  3. `if abs[b["T_K"] - d["T_K"]]>1.0e-6:` -- square brackets instead of
     parens; indexes the `abs` builtin instead of calling it
     (`TypeError: not subscriptable`).
  4. `return lines` indented one level too deep, inside the outer
     `for b,d in zip(...)` loop -- would have returned after only the
     FIRST converged (bubble,dew) pair, giving every quality line just
     1 point instead of the full T sweep.
  5. The real `plot_envelope(...)`/`plot_envelope_honeywell_units(...)`
     calls in `_cli()` never passed `quality_lines=quality_lines` at
     all -- so even with 1-4 fixed, no quality lines would ever reach
     either plot.
  Also found and removed a 6th issue: an orphaned duplicate constant
  `Quality_line_values = [...]` (wrong case, unused dead code) sitting
  near `BTU_LBM_TO_KJ_KG` -- likely a first attempt at change #1 that
  was never wired up and left behind.
- Fixed all 6 directly (user said "Fix em"): moved the constant
  definition to `QUALITY_LINE_VALUES` right above `compute_quality_lines`,
  fixed the default-arg reference, fixed `abs(...)`, fixed `return lines`
  indentation, moved `quality_lines = compute_quality_lines(...)` into
  `_cli()` right after the real `run_true_vle_envelope(...)` call, added
  `quality_lines=quality_lines` to both plot calls, and deleted the
  orphaned duplicate constant.
- Verified via `python3 -m py_compile mixture_quality_line_validation.py`
  (clean) plus targeted grep confirming no remaining
  `Quality_line_values`/`abs[` patterns and that `QUALITY_LINE_VALUES`/
  the `quality_lines=quality_lines` wiring are present exactly once each
  in the right places. NOT yet run end-to-end (sandbox lacks scipy) --
  next step: user runs the file and checks the resulting plots for thin
  gray quality lines (x=0.1...0.9) inside the dome on both the SI figure
  and the Honeywell-units figure.

### 2026-08-13 (later still) — Quality lines: color fix + critical-point convergence fix

- User: quality lines "should be black in color, difficult to see." Both
  `plot_envelope` and `plot_envelope_honeywell_units`: changed
  `color="gray", alpha=0.7` -> solid `color="black"` (no alpha) for both
  the quality-line strokes and their "x=0.1" etc. labels. `lw=0.75` kept
  unchanged so they stay visually distinct from the bold `lw=2.0` dome
  boundary despite now being the same color.
- User then shared a screenshot showing quality lines fanning out near
  the dome apex WITHOUT converging to a single point, asked "shouldn't
  they meet at the critical point?" -- correct catch, real bug: physically
  at Tc, h_l=h_v (that IS the definition of critical), so h(x)=h_l+x*(h_v-h_l)
  collapses to the same value for every x there -- all quality lines
  should terminate at exactly the solved critical point, same as the
  bubble/dew boundary already does. Root cause: `compute_quality_lines`
  only sees `bubble_rows`/`dew_rows` from the T-sweep, has no knowledge of
  `crit_point` at all, so each line stopped at whatever the last
  converged near-critical T happened to be -- close to, but not exactly,
  Tc, so lines terminated at slightly different points instead of one.
  User then confirmed via a Honeywell chart crop ("This is how Honeywell's
  is") showing the real chart's quality lines DO visibly converge tightly
  at the apex, confirming the fix direction.
- FIX (both plot functions, not `compute_quality_lines` itself -- kept
  the closure logic in the plotting layer, mirroring exactly how the
  bubble/dew boundary's own critical-point closure already works): inside
  the `if quality_lines:` loop, when `have_crit` is True, append
  `(hc,pc)` (SI: `hc_si,pc_si`, converted at the final step in the
  Honeywell-units version) as each quality line's own final point, so
  every line -- boundary AND all 9 quality lines -- shares the exact same
  terminal vertex.
- Verified via `python3 -m py_compile mixture_quality_line_validation.py`
  (clean). NOT yet run end-to-end (sandbox lacks scipy) -- next step:
  user re-runs and confirms quality lines now converge cleanly to the
  black star-free apex, matching the Honeywell reference shape.

### 2026-08-13 (later still) — Legend removed; MILESTONE CLOSED: quality lines validated

- Removed `ax.legend(fontsize=8)` from both `plot_envelope` and
  `plot_envelope_honeywell_units` (2 occurrences) -- plots are now bare
  lines/labels with no legend box. Verified via `python3 -m py_compile`
  (clean).
- MILESTONE CLOSED (user instruction): quality lines in
  `mixture_quality_line_validation.py` are now marked VALIDATED. Recap of
  what "validated" covers here: lever-rule interior states
  (h(x)=h_l+x*(h_v-h_l) at shared P=avg(P_bubble,P_dew)), computed as
  pure post-processing of the already-solved (and previously validated)
  bubble/dew boundary -- no new EOS calls. Rendered as solid black
  `lw=0.75` lines (fixed from an initial hard-to-see gray), each labeled
  "x=0.1" through "x=0.9", all 9 lines closing at the exact solved
  critical-point coordinate (fixed from an initial bug where they
  terminated at scattered near-critical points instead of one) --
  confirmed against a crop of the actual Honeywell TDS p-h chart (page 2)
  showing the same convergent-at-apex behavior. 5 implementation bugs
  (undefined constant reference, stray module-level call, `abs[...]` vs
  `abs(...)`, mis-indented return, missing `quality_lines=` wiring) were
  found and fixed along the way -- see the "user self-implemented plan, I
  found and fixed 5 real bugs" breadcrumb entry above for the full list.
  No numeric Honeywell ground-truth exists for quality-line PLACEMENT
  specifically (the TDS only publishes the P-T table numerically; the
  p-h chart with printed quality lines is an image) -- "validated" here
  means lever-rule correctness + critical-point convergence + visual
  shape match against that chart image, not a MAPE-style numeric
  cross-check like the P-T table got.
- NEXT MILESTONE (not started): isotherms (constant-temperature lines)
  through the two-phase region, referencing the Honeywell chart's red
  temperature lines (e.g. "220 F", "200 F", "180 F" visible in the
  crop the user shared). User will make their OWN copy of
  `mixture_quality_line_validation.py` for this work (mirroring how
  `mixture_quality_line_validation.py` itself was a copy of
  `Honeywell_T_P_revalidation.py`) -- explicit instruction: "Don't do
  anything now." No code changes made or planned until the user directs
  otherwise on the new copy.

### 2026-08-13 (later still) — Isotherms milestone started: user copy created, physics explained, plan given (not yet built)

- User created `mixture_isotherm_validation.py` (their own copy, per the
  pattern above). Explained the physics before any code: unlike quality
  lines (diagonal, lever-rule interpolation between bubble/dew), an
  isotherm INSIDE the two-phase dome is a horizontal segment, because T
  and P are not independent inside the dome -- at fixed T there is
  exactly one saturation P, shared by both bubble and dew ends. So an
  isotherm at T is just the segment from (h_l,P) to (h_v,P) at that T --
  the same P-equality fact already enforced inside the 3-equation solve,
  no lever rule needed.
- Flagged the one real design choice: Honeywell draws isotherms at round
  Fahrenheit values (0,20,...,220F -- matches the crop the user shared
  earlier), which won't land on the existing T-sweep's Kelvin grid.
  Offered two options: interpolate from the existing sweep (cheap, small
  error) vs. solve FRESH at the exact round-F temperatures (seeded from
  the nearest converged sweep point). Recommended the fresh-solve option
  as consistent with this project's general preference for solving
  exactly over approximating (e.g. the critical point being solved
  rather than interpolated). User did not object to that recommendation.
- User explicit scope decision: "Let's plot inside two-phase first" --
  i.e. isotherms OUTSIDE the dome (subcooled liquid, superheated vapor)
  are explicitly deferred; those need a genuinely different calculation
  (solving the single-phase EOS for density at an arbitrary (T,P), not
  just referencing the bubble/dew boundary) and are NOT part of this
  milestone yet.
- Gave a full 7-change plan for `mixture_isotherm_validation.py`:
  - `ISOTHERM_VALUES_F = [0,20,...,220]` constant.
  - `compute_isotherms_two_phase(d1, d2, z1, bubble_rows, dew_rows,
    t_values_f=ISOTHERM_VALUES_F) -> Dict[float, Dict]`: for each target
    T_F, converts to K via the existing `_f_to_k` helper, finds the
    nearest already-converged sweep point in bubble_rows/dew_rows for
    seeding, then calls `solve_bubble_at_t`/`solve_dew_at_t` FRESH at the
    exact target T (real EOS solve, not interpolation), records
    `{T_F, T_K, P_Pa (avg), h_l_Jmol, h_v_Jmol}` if both converge.
  - `plot_envelope` and `plot_envelope_honeywell_units`: add an
    `isotherms: Optional[Dict[float, Dict]] = None` parameter to each,
    and a plotting loop drawing the 2-point horizontal segment per
    isotherm (`lw=1.0`, black, labeled e.g. "220 F"); the Honeywell-units
    version converts via `to_honeywell_units` at the final step only,
    same architecture as quality lines.
  - `_cli()`: compute `isotherms = compute_isotherms_two_phase(d1, d2,
    z1, bubble_rows, dew_rows)` after `d1`/`d2` are loaded, pass
    `isotherms=isotherms` into both plot calls.
- STATUS: plan given, NOT yet applied to `mixture_isotherm_validation.py`
  -- no edits made to that file. Next step: user applies (or asks me to
  apply) the 7 changes above.

### 2026-08-13 (later still) — Isotherms: user partially self-implemented, I completed + wired + verified, scoped to ONE test isotherm

- User had started applying the plan themselves: `ISOTHERM_VALUES_F`
  (full 12-value list) and the start of `compute_isotherms_two_phase`
  were present, but the function was cut off mid-statement with a real
  syntax error -- `nearest_b = min(converged_bubble, key=lambda
  r: abs(r["T_K"]-t_k), default=)` (empty `default=`), and the rest of
  the function body (finding `nearest_d`, the two solve calls, the
  CONVERGED check, building/returning the `isotherms` dict) was entirely
  missing -- the `def` fell straight through into `def plot_envelope(`.
- User asked me to "add the plot envelope functions, check everything,
  let us plot one isotherm first." Did all three:
  1. Scoped `ISOTHERM_VALUES_F` down to `[100.0]` (one test value) for
     this first validation pass -- comment left in place noting to
     expand back to the full `[0,20,...,220]` set once confirmed working.
  2. Completed `compute_isotherms_two_phase`: replaced the broken
     `default=` pattern with an explicit `if not converged_bubble or not
     converged_dew: return isotherms` early-out guard (avoids needing
     `min()`'s `default` kwarg at all), then finds `nearest_d`
     independently, calls `solve_bubble_at_t`/`solve_dew_at_t` FRESH at
     the exact target `t_k` (verified both signatures match the actual
     calls by reading them directly: `(d1,d2,t_k,z1,rho_l0,rho_v0,y10)`
     and `(...,x10)` respectively -- no `dew_rho_l0`/`dew_rho_v0` extra
     param in this file's lineage, unlike the abandoned
     `mixture_quality_validation.py`), checks both CONVERGED, and builds/
     returns the `isotherms` dict.
  3. Added `isotherms: Optional[Dict[float, Dict]] = None` to both
     `plot_envelope` and `plot_envelope_honeywell_units`, plus a plotting
     block in each drawing the two-point horizontal segment (`lw=1.0`,
     black, labeled e.g. "100 F") -- Honeywell-units version converts via
     `to_honeywell_units` at the final step only, matching the standing
     architecture.
  4. Wired into `_cli()`: `isotherms = compute_isotherms_two_phase(d1,
     d2, z1, bubble_rows, dew_rows)` added right after `quality_lines =
     ...`, and `isotherms=isotherms` added to both plot calls.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean) plus targeted grep confirming no `default =)`/`abs[` patterns
  remain and that the constant/function/params/wiring each appear
  exactly once in the right places.
- NOT yet run end-to-end (sandbox lacks scipy). Next step: user runs the
  file and checks for a single horizontal "100 F" line inside the dome
  on both figures; once confirmed, expand `ISOTHERM_VALUES_F` back to
  the full `[0,20,...,220]` Honeywell set.

### 2026-08-13 (later still) — Isotherms expanded to full Honeywell set

- User confirmed the single "100 F" test isotherm and said "let's plot
  all lines." `ISOTHERM_VALUES_F` in `mixture_isotherm_validation.py`
  expanded from `[100.0]` back to the full `[0.0, 20.0, 40.0, ...,
  220.0]` (12 values, matching the Honeywell TDS p-h chart's printed
  isotherm labels, page 2). No other changes needed -- everything else
  (`compute_isotherms_two_phase`, the two plotting hooks, the `_cli()`
  wiring) was already generic over the list length.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean). NOT yet run end-to-end (sandbox lacks scipy) -- next step:
  user runs the file and checks for all 12 horizontal isotherm lines
  inside the dome on both figures.

### 2026-08-13 (later still) — Isotherms colored red (matches Honeywell); labels hard to read against black quality lines/boundary

- User shared a zoomed screenshot of the dome apex: with everything
  (boundary, quality lines, isotherms) drawn black, the isotherm labels
  were hard to distinguish from the diagonal quality lines crossing
  them. Requested: "Make isotherms red like in Honeywell."
  (Honeywell's own p-h chart does use red for its isotherm family,
  confirmed from the full-chart image shared earlier this session --
  green=density, blue=entropy, red=temperature, black=quality/saturation
  boundary.)
- Changed both isotherm plotting blocks (`plot_envelope` and
  `plot_envelope_honeywell_units`) from `color="black"` to `color="red"`
  for both the horizontal segment (`ax.plot`) and its "{T} F" label
  (`ax.annotate`) -- 2 lines changed in each function, 4 total. Quality
  lines and the bubble/dew boundary remain black, unchanged.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean) plus grep confirming both isotherm annotate calls now say
  `color="red"`. NOT yet run end-to-end (sandbox lacks scipy) -- next
  step: user re-runs and confirms isotherm lines/labels are now clearly
  legible in red against the black boundary and quality lines.

### 2026-08-13 (later still) — Accidental overwrite lost downstream isotherm wiring; recovered + fixed 2 new bugs

- User: "Next, we tackle the liquid side" (subcooled-liquid isotherm
  extension, deferred milestone) -- before I could scope that, user sent
  a rapid sequence: "Something was not saved" / "In the file" / "I am
  going to overwrite" / "I will need you to go in and fix." User had
  overwritten `mixture_isotherm_validation.py` with their own version of
  `compute_isotherms_two_phase`, which wiped out everything downstream
  of it that I'd added: the `isotherms` parameter on both `plot_envelope`
  and `plot_envelope_honeywell_units`, both red-isotherm plotting blocks,
  and the `_cli()` wiring (`isotherms = compute_isotherms_two_phase(...)`
  + `isotherms=isotherms` on both plot calls) -- confirmed via grep
  finding zero matches for any of that.
- User's rewritten `compute_isotherms_two_phase` also introduced 2 new
  bugs:
  1. `nearest_d = min(converged_dew, key=lambda r: abs(r["T_k"]-t_k),
     default=None)` -- lowercase `T_k`, but every row dict uses `T_K`
     (capital K). Would `KeyError` the instant this ran.
  2. `"h_v_Jmol": b["h_v_Jmol"]` -- pulled the vapor enthalpy from `b`
     (bubble result, vapor state at bubble's own solved `y1`) instead of
     `d` (dew result, vapor state at the fixed feed composition z1) --
     would silently plot a physically wrong isotherm vapor endpoint even
     though it wouldn't crash. Also normalized `"T_f"` -> `"T_F"` (key
     casing) for consistency with the rest of the file, though this one
     wasn't a functional bug since the plotting code keys off the outer
     dict key, not this inner one.
- Fixed both bugs in place, then re-added everything that was lost: the
  `isotherms` parameter on both plot functions, both red-line plotting
  blocks (unchanged from before -- same `lw=1.0`, `color="red"` styling
  from the prior breadcrumb entry), and the `_cli()` wiring.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean), plus grep confirming zero remaining `"T_k"`/`b["h_v_Jmol"]`
  patterns and that `isotherms=isotherms`/`color="red"` etc. are all
  present again in exactly the right places.
- LESSON for future sessions: if the user reports "something wasn't
  saved" after they've edited a file locally, always re-grep the actual
  on-disk state before assuming prior work is still there -- local
  edits/overwrites from the user's own editor are NOT visible to this
  session until read fresh.
- Isotherms (two-phase region) milestone marked completed in the task
  list. Next up (not yet scoped/started): subcooled liquid side --
  isotherm segments extending OUTSIDE the dome, which need a genuinely
  different single-equation/single-unknown (rho_l at fixed T,P) solve
  rather than the bubble/dew boundary lookup used here.

### 2026-08-13 (later still) — Liquid-side (subcooled) isotherm extension: explained, scoped, implemented, wired, verified

- User: "Now lets do the other isotherms" (the subcooled-liquid side,
  deferred from the entry above). Explained the physics first, no code,
  per standing pattern:
  - Outside the dome there is only ONE phase at the fixed feed
    composition z1 -- no equilibrium condition, so unlike bubble/dew this
    is a single equation in a single unknown: find rho such that
    `mix_state(T,rho,z1).p_pa == P_target`, solved via bracketed
    `brentq` (already imported), not `least_squares`.
  - Physically, liquid is nearly incompressible, so h rises only slightly
    with P at fixed T -- each isotherm continues from its bubble point
    almost straight up with a slight rightward bend, matching the
    near-vertical lines left of the dome on a real p-h chart.
  - Open question flagged: how far up in P to extend. Initially proposed
    a placeholder -- `LIQUID_EXT_P_MAX_FACTOR = 1.2` (1.2x solved
    critical pressure) with `LIQUID_EXT_N_POINTS = 5` -- pending a better
    answer.
- User then shared the full-resolution Honeywell p-h chart image:
  "Objective is to reproduce this." This settled the open question: the
  chart's own y-axis ceiling is 1350 psia, and the subcooled-liquid
  isotherm extensions are drawn in BLACK (distinct from the red in-dome
  segments). Superseded the 1.2x-Pc placeholder with the chart's actual
  ceiling.
- Implemented in `mixture_isotherm_validation.py`:
  - Added `"rho_l_molm3": b["rho_l_molm3"]` to `compute_isotherms_two_phase`'s
    returned per-T dict -- needed to seed the liquid-side root-find.
  - `LIQUID_EXT_P_MAX_PA = 1350.0 * 6894.757293168361` (hardcoded PSI->Pa
    factor, not the later-defined `PSI_TO_PA` constant, to avoid a
    module-load-order `NameError`) and `LIQUID_EXT_N_POINTS = 6`.
  - `compute_isotherms_liquid_side(d1, d2, z1, isotherms, p_max_pa=LIQUID_EXT_P_MAX_PA,
    n_points=LIQUID_EXT_N_POINTS) -> Dict[float, List[Dict]]`: for each
    isotherm, starts at `(P_Pa=p_sat, h=h_l_Jmol)`, then for each
    subsequent target pressure in `np.linspace(p_sat, p_max_pa,
    n_points)[1:]` brackets and solves `P(T,rho)=P_target` via `brentq`,
    seeded from the previous point's density (`rho_hi = min(rho_seed*1.5,
    RHO_MAX_MOLM3)`), stopping early (rather than guessing a wider
    bracket) if a `ValueError` (no sign change) is raised.
  - Wired `isotherms_liquid` parameter + a black-line plotting block into
    both `plot_envelope` (SI) and `plot_envelope_honeywell_units` (IP,
    converted point-by-point via the existing `_convert_all`/
    `to_honeywell_units` helpers, final-step-only, consistent with the
    standing unit-conversion architecture).
  - Wired into `_cli()`: `isotherms_liquid = compute_isotherms_liquid_side(d1,
    d2, z1, isotherms)`, passed to both `plot_envelope(...)` and
    `plot_envelope_honeywell_units(...)` calls.
- Mid-session bookkeeping snag (self-inflicted, not a user overwrite this
  time): my own edit calls landed the file in a temporarily messy state
  -- an early, superseded copy of the constants (`LIQUID_EXT_P_MAX_FACTOR
  = 1.2`, `LIQUID_EXT_N_POINTS = 5`, from the placeholder plan above) was
  still sitting in the file alongside my later, chart-based replacement
  (`LIQUID_EXT_P_MAX_PA`, `LIQUID_EXT_N_POINTS = 6`), plus a leftover
  dead inner function `f` in `compute_isotherms_liquid_side` (copy-paste
  debris from when the working closure was renamed to `g`). Caught via a
  targeted grep (`^LIQUID_EXT|^def compute_isotherms_liquid_side`) before
  trusting the file was clean. Fixed: removed the orphaned
  `LIQUID_EXT_P_MAX_FACTOR`/first `LIQUID_EXT_N_POINTS=5` pair (unused --
  nothing referenced them) and the dead `f` function.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean) and a follow-up grep confirming exactly one `LIQUID_EXT_P_MAX_PA`
  and one `LIQUID_EXT_N_POINTS` definition remain. NOT yet run end-to-end
  (sandbox lacks scipy) -- next step: user runs the file and checks the
  black subcooled-liquid isotherm extensions render correctly (near-
  vertical, all reaching the 1350 psia ceiling) alongside the existing
  red in-dome isotherms, black quality lines, and black boundary.
- Liquid-side isotherm extension milestone marked completed in the task
  list. Remaining unscoped pieces of the full Honeywell p-h chart
  reproduction: blue entropy lines, green specific-volume/density lines,
  and the red isotherm continuations into the superheated-vapor region
  (diagonal lines above/right of the dome) -- none of these discussed yet.

### 2026-08-13 (later still) — Liquid-side isotherms recolored black->red

- User: "Isotherms should all be in red." The liquid-side extension had
  been colored black (distinguishing it from the in-dome red segments);
  user wants one continuous red line per isotherm instead. Changed
  `color="black"` -> `color="red"` for the `isotherms_liquid` plotting
  block in both `plot_envelope` and `plot_envelope_honeywell_units` (1
  line each). In-dome segments were already red; no change there.
  Verified via `python3 -m py_compile` (clean).

### 2026-08-13 (later still) — Added -20F isotherm (Honeywell chart has 2 more than our 0-220F set)

- User shared a crop of the actual Honeywell chart's low-enthalpy corner
  showing more red isotherm lines than our set produces there; asked to
  identify the 2 missing ones. After a couple rounds of back-and-forth
  (the crop resolution made exact labels hard to read), user confirmed:
  0F (which has a two-phase/in-dome portion) and -20F.
  - 0.0 was already in `ISOTHERM_VALUES_F` -- no change needed for that
    one; flagged as worth confirming its in-dome red segment actually
    renders (depends on whether the envelope sweep's --Tmin is low
    enough that a fresh bubble/dew solve at T=0F has a nearby seed to
    converge from -- not yet verified end-to-end, sandbox lacks scipy).
  - -20.0 was NOT in the list. User self-added it directly in their
    editor: `ISOTHERM_VALUES_F = [-20.0, 0.0, 20.0, 40.0, ..., 220.0]`
    (13 values total now, still 20F steps above 0F, with -20F as a new
    low-end extra).
- Re-grepped the file per standing lesson (always verify on-disk state
  after user edits) -- confirmed `ISOTHERM_VALUES_F` is the only change,
  all downstream wiring (`compute_isotherms_two_phase`,
  `compute_isotherms_liquid_side`, both plot functions, `_cli()`) is
  generic over the list and untouched/still correct. Verified via
  `python3 -m py_compile mixture_isotherm_validation.py` (clean).
- NOT yet run end-to-end -- next step: user reruns and confirms both the
  -20F isotherm's two-phase segment (if the sweep range/solver converges
  that low) and its liquid-side extension render correctly, and that the
  0F isotherm's in-dome portion is now visible too.

### 2026-08-13 (later still) — Liquid-side isotherm was kinked, not smooth; switched to log-spaced pressure sampling

- User shared a crop of the rendered liquid-side extension: a visible
  sharp elbow partway up, shallow segment near the dome then abruptly
  steep for the rest, vs. Honeywell's smooth curve ("Our isotherm single
  phase is a kinky line, Honeywell is smooth").
- Diagnosis: not a physics bug -- an under-sampling/plotting artifact.
  `compute_isotherms_liquid_side` only had `LIQUID_EXT_N_POINTS = 6`
  points, spaced via `np.linspace` (linear in Pa) across a huge range
  (as low as ~15 psia up to the 1350 psia ceiling), while the plot's
  P-axis is log-scaled (`ax.set_yscale("log")`). Liquid compressibility
  (dh/dP) is largest right at the saturation line and flattens out
  moving away from it, so nearly all the curve's real bend landed inside
  the first linear-spaced segment -- drawn as one straight, visibly
  kinked line instead of a smooth curve.
- Fix, both in `mixture_isotherm_validation.py`:
  - `LIQUID_EXT_N_POINTS` raised 6 -> 25.
  - `compute_isotherms_liquid_side`: `np.linspace(p_sat, p_max_pa,
    n_points)` -> `np.geomspace(p_sat, p_max_pa, n_points)` -- log-spaced
    pressure targets land at visually even intervals on the log P-axis,
    and concentrate more points near p_sat where the curve actually
    bends, instead of wasting resolution on the nearly-linear high-P
    tail. Root-finding logic (brentq, bracket seeded from the previous
    point's density) unchanged.
  - Updated the function's docstring to record this reasoning.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean). NOT yet run end-to-end -- next step: user reruns and confirms
  the liquid-side lines now render as smooth curves.

### 2026-08-13 (later still) — Honeywell reference: superheated-region isotherm continuation (not yet built)

- User shared a crop of the actual Honeywell chart's upper-right region
  (outside/above the dew line): a red isotherm labeled "240" curving
  over the top, alongside a green (density) line and a blue
  "0.35 Btu/lb-R" (entropy) line. This is the superheated-vapor-region
  isotherm continuation already flagged as unscoped in the prior
  "Liquid-side isotherm extension" entry's closing note.
- No code changes made -- logging this reference for when that milestone
  is scoped. Still-unbuilt pieces of the full Honeywell p-h chart
  reproduction, per that same note: blue entropy lines, green
  specific-volume/density lines, and red isotherm continuations into the
  superheated-vapor region (this crop). Liquid-side (subcooled) isotherm
  work is otherwise considered done pending the smoothness re-check
  above.

### 2026-08-13 (later still) — Liquid-side smoothness fix confirmed working ("Works"); vapor-side (superheated) milestone scoped, not yet built

- User confirmed the geomspace/25-point fix rendered smoothly ("Works").
  Liquid-side (subcooled) isotherm extension milestone fully closed.
- User: "Now, what do we do for vapor side, explain first." Explained
  conceptually, no code yet, per standing pattern:
  - Same single-equation/single-unknown architecture as the liquid side
    (still only one phase outside the dome, no VLE machinery), but
    mirrored in direction: starts at the DEW point and pushes P DOWN
    (toward superheat) instead of starting at the bubble point and
    pushing P up (toward subcooling). Same `mix_state(T,rho,z1).p_pa ==
    P_target` equation, same bracketed `brentq`, just walking the
    density seed downward each step instead of upward.
  - Needs `rho_v_molm3` added to `compute_isotherms_two_phase`'s stored
    dict (from the dew-row solve), mirroring how `rho_l_molm3` (from the
    bubble-row solve) was added to seed the liquid-side extension.
  - Pressure floor mirrors the liquid ceiling: the chart's own bottom
    axis value, 15 psia (converted to Pa the same hardcoded way as
    `LIQUID_EXT_P_MAX_PA`), instead of the 1350 psia top.
  - Sampling: applying the just-learned lesson immediately -- log-spaced
    (`np.geomspace`) pressure targets from the start, not linear, since
    the P-axis is log-scaled and linear spacing already proved to
    produce a visibly kinked line on the liquid side.
  - Expected shape: at fixed T, h should rise slightly as P drops toward
    the ideal-gas limit (departure function shrinks), so these should
    slope up-and-right away from the dew line -- consistent with the
    diagonal "240" isotherm crop shared earlier.
  - Planned new function `compute_isotherms_vapor_side` (mirrors
    `compute_isotherms_liquid_side`'s structure/signature), plotted red
    (same color as the rest of each isotherm), via a new
    `isotherms_vapor` param on both plot functions + `_cli()` wiring.
  - NOT yet built -- awaiting explicit go-ahead before writing code, per
    standing pattern (explain conceptually first, build only after
    "go ahead"/"build it"/equivalent).

### 2026-08-13 (later still) — Correction: requested 240-400F isotherms are supercritical, not dew-point vapor extensions

- User: "We need 240,260,280,300,320,340,360,380,400" (from the same
  superheated-region crop), then clarified: "These are vapor-only lines
  in addition to the existing range from -20F to 220F in 20F intervals."
- Checked our own solved critical temperature before scoping further:
  Tc = 382.045 K (from the Honeywell revalidation, 0.04% off Honeywell's
  published value) = 228.0F. ALL nine requested values (240-400F) are
  above Tc.
- This invalidates the just-planned `compute_isotherms_vapor_side`
  approach (extend from an existing dew point) for these specific
  values -- above Tc there is no phase separation at all, so there's no
  bubble/dew point to anchor from. These are genuinely different: a
  single continuous phase across the ENTIRE chart pressure range
  (15-1350 psia), not an extension of a two-phase point.
- Revised plan, explained to user, not yet built:
  - New constant `VAPOR_ONLY_ISOTHERM_VALUES_F = [240, 260, 280, 300,
    320, 340, 360, 380, 400]`, kept separate from `ISOTHERM_VALUES_F`
    (-20 to 220F, the set with real bubble/dew points).
  - New function `compute_isotherms_supercritical`: same governing
    equation/brentq approach as the liquid/vapor extensions, but with no
    saturation point to seed from -- starts at the low-P end (15 psia)
    using an ideal-gas density estimate (rho ~ P/(R*T)) as the initial
    bracket guess, then walks up to 1350 psia in log-spaced
    (`np.geomspace`) steps, seeding each bracket from the previous
    point's solved density (mirrors the walking pattern already used in
    `compute_isotherms_liquid_side`, just spanning the full range in one
    pass with no saturation anchor).
  - Plotted red (consistent with every other isotherm), wired as a third
    isotherm dict (`isotherms_supercritical` or similar) alongside
    `isotherms`/`isotherms_liquid` in both plot functions + `_cli()`.
  - Whether the originally-planned `compute_isotherms_vapor_side` (dew-
    point-anchored superheated extension for the EXISTING -20 to 220F
    subcritical set) is still needed is a separate, still-open question
    -- not addressed by this correction, not yet built either.
  - NOT yet built -- awaiting explicit go-ahead.

### 2026-08-13 (later still) — Supercritical isotherm function built + wired; explicit axis limits set

- User asked "How did we guess the liquid one?" -- explained that
  `compute_isotherms_liquid_side` never guesses a density from nothing:
  its first point reuses the real solved `rho_l_molm3` from the bubble
  solve, and every point after that seeds from the PREVIOUS point's own
  solved density (bracket `[rho_seed, min(rho_seed*1.5, RHO_MAX_MOLM3)]`)
  -- the only judgment call is that 1.5x multiplier, never the density
  itself. This is what breaks for the supercritical case (no bubble
  point at all to start from).
- User: "Ok, please write the function." Implemented in
  `mixture_isotherm_validation.py`, placed after
  `compute_isotherms_liquid_side`:
  - `VAPOR_ONLY_ISOTHERM_VALUES_F = [240, 260, 280, 300, 320, 340, 360,
    380, 400]`, `SUPERCRIT_P_MIN_PA = 15.0 * 6894.757293168361` (chart's
    own y-axis floor, same hardcoded-conversion pattern as
    `LIQUID_EXT_P_MAX_PA` for module-load-order reasons),
    `SUPERCRIT_N_POINTS = 25`.
  - `compute_isotherms_supercritical(d1, d2, z1, t_values_f=...,
    p_min_pa=SUPERCRIT_P_MIN_PA, p_max_pa=LIQUID_EXT_P_MAX_PA,
    n_points=SUPERCRIT_N_POINTS) -> Dict[float, List[Dict]]`: same
    governing `mix_state(T,rho,z1).p_pa==P_target` / bracketed `brentq`
    approach as the liquid extension, log-spaced (`np.geomspace`) from
    the start. Walks LOW-P to HIGH-P deliberately -- the very first
    point has no real neighboring density, so its bracket comes from an
    ideal-gas estimate `rho ~ P/(R*T)` (R_GAS=8.314462618 J/mol-K, used
    ONLY to size that one bracket, never in the actual EOS call), which
    is most reliable at low P where real-gas departure is smallest.
    Every subsequent point reuses the previous point's actual solved
    density (wider bracket than the liquid case, 0.5x-2x vs 1.0x-1.5x,
    since supercritical density can shift more per pressure step near
    the Widom line), with one widening retry (0.1x-5x) before giving up
    on a given point.
  - Wired `isotherms_supercritical` parameter + a red plotting block
    (same color as every other isotherm segment) into both
    `plot_envelope` and `plot_envelope_honeywell_units`; wired
    `compute_isotherms_supercritical(d1, d2, z1)` into `_cli()`, passed
    to both plot calls alongside `isotherms`/`isotherms_liquid`.
- Separately, user specified explicit axis extents: "Axes limits:
  Pressure in [15,1400] psia, [enthalpy] in [70,230] Btu/lbm]." Added
  `ax.set_xlim(70.0, 230.0)` / `ax.set_ylim(15.0, 1400.0)` to
  `plot_envelope_honeywell_units` only (the IP-units validation plot --
  `plot_envelope` stays in SI/bar/kJ-kg, not addressed by this
  instruction). Note: 1400 psia is the axis view limit; the liquid/
  supercritical extension ceiling itself (`LIQUID_EXT_P_MAX_PA`) is
  still 1350 psia, unchanged -- left as a small margin between the data
  and the axis edge unless told otherwise.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean) plus grep confirming `isotherms_supercritical` and the new
  axis-limit calls appear in exactly the expected places. NOT yet run
  end-to-end (sandbox lacks scipy) -- next step: user runs the file and
  checks the 9 new supercritical isotherms render as smooth curves
  across the full pressure range, and that both IP-unit figure axes now
  match the Honeywell chart's own extents.

### 2026-08-13 (later still) — Built the missing vapor-side (superheated) extension for the -20 to 220F set

- User shared a screenshot of the rendered dome and said "it did not
  plot older connections." Code review of the two-phase isotherm block
  (`if isotherms:`, drawing `ax.plot([h_l,h_v],[p,p],color="red")`)
  showed no wiring bug -- labels and lines share the same loop/dict
  entry, so if the "220 F"..."20 F" labels were rendering, the in-dome
  segments should be too. Asked for clarification.
- User clarified: "The vapor side of 0-220F." This was the originally-
  planned `compute_isotherms_vapor_side` (dew-point-anchored superheated
  extension for the EXISTING -20 to 220F subcritical set) -- explicitly
  left as an open question in the 2026-08-13 "Correction: requested
  240-400F isotherms are supercritical" entry above, and never actually
  built (only `compute_isotherms_supercritical`, for the DIFFERENT
  240-400F set with no dome anchor at all, was built). Confirmed now
  needed.
- Implemented in `mixture_isotherm_validation.py`, mirroring
  `compute_isotherms_liquid_side` with the direction reversed:
  - Added `"rho_v_molm3": d["rho_v_molm3"]` to
    `compute_isotherms_two_phase`'s returned per-T dict (from the dew-
    row solve) -- needed to seed the vapor-side root-find, same role
    `rho_l_molm3` (from the bubble row) plays for the liquid side.
  - `VAPOR_EXT_P_MIN_PA = 15.0 * 6894.757293168361` (chart's own y-axis
    floor, same hardcoded-conversion pattern as `LIQUID_EXT_P_MAX_PA`)
    and `VAPOR_EXT_N_POINTS = 25`.
  - `compute_isotherms_vapor_side(d1, d2, z1, isotherms,
    p_min_pa=VAPOR_EXT_P_MIN_PA, n_points=VAPOR_EXT_N_POINTS) ->
    Dict[float, List[Dict]]`: starts at the dew point `(P_Pa=p_sat,
    h=h_v_Jmol)`, walks P DOWN via `np.geomspace(p_sat, p_min_pa,
    n_points)[1:]` (log-spaced from the start, applying the liquid-
    side's kink-avoidance lesson immediately rather than relearning it).
    Each step's bracket searches BELOW the previous point's density
    (`rho_hi = rho_seed`, `rho_lo = max(rho_seed/1.5, 1e-3)`) since
    vapor density falls as pressure drops, mirroring the liquid side's
    upward-searching bracket exactly in reverse. Same brentq/ValueError-
    break pattern otherwise.
  - Wired `isotherms_vapor` parameter + a red plotting block (same color
    as every other isotherm segment) into both `plot_envelope` and
    `plot_envelope_honeywell_units`; wired
    `compute_isotherms_vapor_side(d1, d2, z1, isotherms)` into `_cli()`,
    passed to both plot calls alongside the other three isotherm dicts.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean) plus grep confirming `isotherms_vapor` appears in all expected
  places (dict-building, both function signatures, both plotting blocks,
  `_cli()` computation + both plot calls). NOT yet run end-to-end
  (sandbox lacks scipy) -- next step: user reruns and confirms the
  superheated-vapor continuations now render for the -20 to 220F set,
  connecting smoothly to each in-dome red segment at the dew point.

### 2026-08-13 (later still) — Chart polish: missing supercritical labels, smoothness, box border, exact Honeywell axis convention

- User shared a rendered crop of our supercritical-region output: the
  240-400F lines were plotted but completely unlabeled, and looked less
  smooth than Honeywell's. Then shared the FULL, high-resolution
  Honeywell Solstice N15 (R-515B) reference chart (previously only seen
  in partial crops) -- this pinned down exact conventions that had been
  inferred/approximated up to now.
- Root cause of the missing labels: `compute_isotherms_supercritical`'s
  plotting blocks only ever had `ax.plot(...)`, no `ax.annotate(...)` --
  unlike the subcritical isotherms, which get their single label for
  free from the original two-phase block (anchored at the dew point),
  the supercritical lines never had ANY code path that could label them
  (no dew point to anchor to). Confirmed by re-grepping both plotting
  blocks before editing.
- Fixes made in `mixture_isotherm_validation.py`:
  1. Added `ax.annotate(f"{t_f:.0f} F", (h_sc[-1], p_sc[-1]), ...)` to
     both `isotherms_supercritical` plotting blocks (`plot_envelope` and
     `plot_envelope_honeywell_units`) -- labeled at each line's last
     (highest-pressure) computed point, mirroring how Honeywell's own
     chart places these labels near the top of the chart where the fan
     of lines is most spread out and legible (confirmed directly against
     the reference image: "240 F" through "380 F" sit in a row around
     1050-1200 psia).
  2. `VAPOR_EXT_N_POINTS` 25->40 and `SUPERCRIT_N_POINTS` 25->40 (both
     already log-spaced via `np.geomspace` per the earlier kink-avoidance
     fix; more points was the same lever that fixed the liquid side's
     smoothness).
  3. Added an explicit box border to BOTH plot functions: `for spine in
     ax.spines.values(): spine.set_visible(True);
     spine.set_linewidth(1.0)` -- matches Honeywell's own fully-boxed
     chart style.
  4. `plot_envelope_honeywell_units` axis convention now matches the
     reference chart exactly (read directly off the full image, not
     approximated): axis label text changed from square brackets to
     parentheses (`"Enthalpy (Btu/lbm)"`, `"Pressure (psia)"`), and the
     y-axis now uses Honeywell's own explicit, non-uniform gridline
     values -- `[15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 300, 450,
     600, 750, 900, 1050, 1200, 1350]` -- via `ax.set_yticks(...)` +
     `ax.set_yticklabels(...)`, with matplotlib's default log-axis minor
     ticks suppressed (`ax.set_yticks([], minor=True)`) since they don't
     match Honeywell's set. X-axis ticks set to Honeywell's own 20-unit
     steps: `[70, 90, 110, 130, 150, 170, 190, 210, 230]`. `plot_envelope`
     (SI, bar/kJ-kg) was NOT touched by the label/tick changes -- only
     the box border was added there, since the bracket-vs-parenthesis
     and custom-tick instructions were specifically about matching
     Honeywell's own (IP-unit) chart.
- Verified via `python3 -m py_compile mixture_isotherm_validation.py`
  (clean) plus grep confirming the annotate calls, spine block, and new
  axis label/tick calls all appear in the expected places in both
  functions. NOT yet run end-to-end (sandbox lacks scipy) -- next step:
  user reruns and visually compares the output directly against the
  full reference chart now on hand.

### 2026-08-13 (later still) — Isochores confirmed skipped (density-bias reasoning); isentropes scoped next; quality-line-label report under investigation

- User asked what the blue ("Btu/lb-R") and green ("lbm/ft^3") lines on
  the Honeywell chart are. Identified: blue = isentropes (constant
  specific entropy, Btu/(lbm-R) IS the entropy unit); green = isochores
  (constant density, equivalently constant specific volume).
- User: "I think we skip isochors since we get a huge error if we do,"
  then pointed at the existing breadcrumb section "Physics established
  this session: WHY dome pressure is (partially) protected from the
  density-extrapolation bias, and where that breaks down" (2026-08-13
  afternoon entry, well above this session's isotherm work) as the
  reasoning to revisit.
- This CORRECTED something said earlier in this session (previous turn):
  isochores were characterized as "easiest, most accurate" since they
  need no brentq root-find (rho is already a native state variable,
  direct evaluation via mix_state(T,rho,z1)). That's true only in the
  narrow numerical-solver sense. The physics point from the referenced
  breadcrumb section is different and controls here: dome PRESSURE is
  only partially protected from the Bell-2023 composition-extrapolation
  bias (fit range x1 in [0.33,0.68]; real R515B x1~=0.9385) because the
  bubble/dew solve is a COUPLED condition -- the well-behaved vapor
  branch anchors the shared pressure value, tempering the ~2% bias that
  otherwise shows up raw in liquid density (Z~0.2 near-cancellation on
  the liquid branch amplifies absolute bias into large relative error).
  An isochore has NO such protection: it's built by picking density
  values directly and sweeping T, with no coupled physical constraint
  pinning the result -- exactly the case the pre-existing M3 code
  comment warns about ("never substitute a real/externally-measured rho
  into this formula directly"). Confirmed: isochores stay skipped, for
  this physics reason, not a numerical-difficulty reason.
- User: "I think we should plot isentropes" -- confirmed as next scope
  target, but "I will create a new copy for isentropes" -- user is
  making their own copy of the file for this work; I am NOT to create
  that file.
- User also reported: "Quality lines have lost their labels. Add
  labels." Checked both `ax.annotate(f"x={q:.1f}", ...)` calls (SI
  `plot_envelope` and `plot_envelope_honeywell_units`) via fresh grep +
  read -- both present and correctly wired, unchanged from the
  already-validated quality-line milestone. NOT a code regression found
  on inspection. Two live hypotheses, neither confirmed: (a) visual
  crowding now that the supercritical/vapor-side isotherm additions put
  many more red lines/labels in the same region, or (b) something in
  actual rendering not visible from source review alone. Asked user for
  a screenshot before making any speculative edit -- NOT yet fixed,
  NOT yet diagnosed with certainty.

### 2026-08-13 (later still) — Quality-line labels: root cause found + fixed (matplotlib annotation clipping); vapor_compression.py checked for density exposure

- User shared the rendered `plot_envelope_honeywell_units` output. Chart
  looked otherwise correct (dome, in-dome red segments + labels 20F-220F,
  supercritical 240F-400F fan with labels, box border, Honeywell axis
  convention all present) but confirmed: zero visible "x=0.1"..."x=0.9"
  quality-line labels anywhere.
- Root cause: matplotlib's `ax.annotate(..., annotation_clip=None)`
  (the default) silently DROPS the entire annotation if its anchor xy
  point falls outside the axes' current view limits -- unlike a plotted
  line, which just clips visually and still partly renders. Quality
  lines' label anchor is `(hq[0], pq[0])`, the FIRST row of the T-sweep
  (lowest T => lowest P). The immediately-preceding change in this same
  session added `ax.set_ylim(15.0, 1400.0)` to this exact plot -- the
  sweep's coldest point's saturation pressure is very likely below the
  15 psia floor, so every quality-line label's anchor fell outside the
  view and got silently dropped, while the LINES themselves still
  rendered fine (line clipping doesn't remove data, just clips
  visually). Isotherm labels survived because their anchors sit inside
  the dome or near the chart top, well within [15,1400].
  Confirmed this is specific to `plot_envelope_honeywell_units`: the SI
  `plot_envelope` never got a `set_ylim`/`set_xlim` call, so it wasn't
  exposed to this failure mode.
- Fix: added `annotation_clip=False` to ALL six `ax.annotate(...)` calls
  in the file (both quality-line label calls, both in-dome isotherm
  label calls, both supercritical isotherm label calls, across both
  plot functions) -- not just the two that were confirmed broken,
  since any of them could hit the identical failure mode if axis limits
  are ever tightened further. Verified via `python3 -m py_compile`
  (clean).
- Separately, user asked: "This should not affect DVCT computations
  right- inability to predict density?" (following the isochore/
  density-bias discussion above) then said "Check vapor_compression.py."
  Read `/Users/snarasi2/idaes-hvacr-cycles/vapor_compression.py`
  (the `SimpleVaporCompressionCycle` DVCT flowsheet) in full. Findings:
  - It uses `idaes.models.properties.general_helmholtz.HelmholtzParameterBlock`
    with `pure_component=fluid_name` -- IDAES's own BUILT-IN single-fluid
    Helmholtz property package, NOT our R1234ze(E)/R227ea custom mixture
    EOS. This file is currently NOT wired to any of today's (or this
    session's) mixture/dome/isotherm work at all -- two separate tracks.
  - Every unit operation (evaporator/compressor/condenser/expansion
    valve) is driven off `pressure`, `temperature`, `enth_mass`,
    `entr_mass`, `vapor_frac` -- property-package OUTPUTS for a given
    state. Raw density is never referenced anywhere in this file.
  - Conclusion given to user: today's changes don't touch DVCT (already
    true architecturally, `plot_envelope` untouched). The density-bias
    exposure question is currently moot for THIS file specifically,
    since it doesn't consume our mixture EOS at all yet. Once the
    planned custom mixture Helmholtz property package (see "Next major
    phase" in the 2026-08-13-afternoon entry, upstream IDAES-PSE PR
    goal) is built and wired into a cycle model like this one, the
    relevant question becomes whether that package exposes T/P/h/s
    (inheriting the coupled dome's partial protection) or raw rho
    (inheriting the full, unprotected bias) to the flowsheet -- not yet
    determined, since that integration hasn't happened.

### 2026-08-13 (later still) — Isentrope validation started (mixture_isentrope_validation.py, user's own copy); entropy calc corrected to direct Helmholtz identity

- User made their own copy `mixture_isentrope_validation.py` (from
  `mixture_isotherm_validation.py`) per their earlier stated intent.
  User: "Ok. Lets do isentrope validation. The points on Honeywell
  datasheet are: 0.22,0.26,0.28,0.3,0.32,0.34,0.35,0.37,0.39,0.41,0.43,
  0.45,0.47,0.49" (Btu/(lbm-R)).
- Mid-scoping, user asked: "No isentropes in 2-phase region because it
  is not an isentropic process?" -- clarified this is a real, standard
  distinction: an isentrope is a line of constant STATE entropy (like an
  isotherm or quality line), not a claim about which process the fluid
  is undergoing. Two-phase mixture entropy is well-defined and
  extensive (lever-rule computable, same as enthalpy), so the two-phase
  segment stays in scope -- confirmed against standard refrigerant p-h
  chart convention and the Honeywell reference image's own blue lines,
  which appear continuous through the dome. User accepted this and
  confirmed proceeding.
- Checked `mix_state` in the new file and found it returns `p_pa,
  h_jmol, g_jmol, rho_mol, rho_mass, x1` -- no entropy field.
- Presented a concrete plan (constants, entropy helper, 3 compute
  functions, plotting wiring, CLI wiring) -- see prior "Tell me what
  needs to be done" exchange -- user confirmed: "Yes, should be in blue
  as in Honeywell. x-axis limits are [70,230] F" (the x-limits restate
  the already-set enthalpy-axis bounds on `plot_envelope_honeywell_units`,
  no new change needed there).
- Implemented in `mixture_isentrope_validation.py`, inserted right
  before `plot_envelope`:
  - `BTU_LBMR_TO_JKGK = 4186.8` (exact, Btu_IT/(lbm-R) -> J/(kg-K), same
    Btu_IT basis as the existing `BTU_LBM_TO_KJ_KG=2.326`, scaled by the
    exact 9/5 R-to-K degree-size ratio), `ISENTROPE_VALUES_BTU_LBMR`
    (the 14 target values), `ISENTROPE_P_MAX_PA`/`ISENTROPE_P_MIN_PA`
    (reuse `LIQUID_EXT_P_MAX_PA`/`VAPOR_EXT_P_MIN_PA`), `ISENTROPE_N_POINTS=40`.
  - `compute_isentropes_two_phase(d1,d2,bubble_rows,dew_rows,s_values_jmolK)`:
    walks the EXISTING bubble/dew T-grid; at each row, s_l(T)/s_v(T) come
    from the entropy helper at the already-solved states; if a target s
    falls within [s_l(T),s_v(T)], solves algebraically for the quality
    x(T) and gets h via the same lever rule as `compute_quality_lines`.
    Pure post-processing, NO root-find -- a genuine simplification versus
    the original plan (which had proposed a P-anchored crossing search).
  - `compute_isentrope_liquid_side` / `compute_isentrope_vapor_side`:
    genuinely new solver pattern for this file -- 2 equations
    ([P(T,rho,z1)=P_target, s(T,rho,z1)=s_target]), 2 unknowns (T,rho),
    solved via `scipy.optimize.root` (method="hybr"), NOT `brentq` (which
    only handles 1 unknown). Anchor: among rows where the isentrope is
    reachable on that branch (s_l(T)>s_target for liquid, s_v(T)<s_target
    for vapor), picks the closest-matching row, warm-starts by solving AT
    that row's own saturation P first, then walks P outward (log-spaced)
    to the 1350 psia ceiling / 15 psia floor, seeding each step from the
    previous point. This single anchor-and-walk approach naturally
    handles BOTH the "crosses the dome" and "never touches the dome"
    cases (an isentrope that's colder/denser or hotter/thinner than
    saturation across the ENTIRE T-grid just produces a full-range
    single-phase line) -- eliminated the need for a separate
    "supercritical-style" 3rd function that had been flagged as the
    least-certain part of the original plan.
  - Wired `isentropes_two_phase`/`isentropes_liquid`/`isentropes_vapor`
    params into `plot_envelope`, plotted blue (matching Honeywell,
    2026-08-13: "in blue as in Honeywell"), `annotation_clip=False` from
    the start (applying the earlier quality-line-label lesson
    immediately rather than relearning it). Each isentrope's 3 possible
    segments share one color/line style so they read as one continuous
    curve; ONE label per isentrope, anchored preferentially at the
    liquid branch's topmost point (falls back to the two-phase segment's
    top, then the vapor branch's bottom), matching Honeywell's own
    near-the-top label placement.
- **Correction (important):** user caught that entropy was computed via
  the Gibbs-energy rearrangement `s=(h-g)/T` (using `mix_state`'s
  already-returned `h_jmol`/`g_jmol`) rather than directly from the
  Helmholtz model: "No. S is directly computed. Go back and look at the
  linar_model/other codes." Investigation (via subagent search across
  the whole codebase) found the established DIRECT pattern in
  `mixture_model_one_point.py`'s `compute_table1_properties()`
  (line ~1601, also present in `linear_model_codex.py`): entropy via the
  standard identity `s/R = tau*alpha_tau - alpha`, there assembled from
  split ideal (`s0_over_r = h0_over_rt - a0_mix - 1`) + residual
  (`s_over_r = s0_over_r + tau*ar_tau - ar`) pieces. Algebraically
  proved these are the SAME formula (h and g are both already
  Helmholtz-derived via the same alpha/alpha_tau, so (h-g)/T reduces to
  exactly `R*(tau*alpha_tau-alpha)`) -- but per the user's explicit
  correction, replaced the h-g shortcut with a NEW function,
  `_mix_entropy_direct(d1,d2,t_k,rho_mol,x1)`, computing
  `R_u*(tau*alpha_tau-alpha)` DIRECTLY from `_mix_alpha_and_derivs`'s
  own alpha/alpha_tau output (recomputing tau independently, same
  self-contained-identity pattern `mix_state` itself already uses) --
  matching the codebase's established "gold" direct-entropy reference
  rather than deriving it as a byproduct of h and g. Replaced ALL 5 call
  sites (`compute_isentropes_two_phase` x2,
  `_isentrope_2eq_residual`, both anchor-finding loops in
  `compute_isentrope_liquid_side`/`compute_isentrope_vapor_side`).
  Verified via `python3 -m py_compile` (clean) and grep confirming zero
  remaining references to the old `_entropy_from_state` name.
- STILL IN PROGRESS at time of writing: `plot_envelope`'s isentrope
  plotting block is wired; `plot_envelope_honeywell_units`'s isentrope
  block and the `_cli()` wiring (converting `ISENTROPE_VALUES_BTU_LBMR`
  to J/(mol-K) via `mw_mix`, calling the three compute functions, passing
  results to both plot calls) are NOT yet done. Next step: finish that
  wiring, `python3 -m py_compile`, then hand off the run command.

### 2026-08-13 (later still) — Isentrope wiring completed

- Finished the remaining pieces: added the identical 3-segment blue
  plotting block (mirroring `plot_envelope`'s, but converting per-point
  via the existing `_convert_all` helper, final-step-only) to
  `plot_envelope_honeywell_units`, right after its supercritical-isotherm
  block.
- Wired `_cli()`: `s_values_jmolK = [s_btu * BTU_LBMR_TO_JKGK * mw_mix for
  s_btu in ISENTROPE_VALUES_BTU_LBMR]` (Btu/(lbm-R) -> J/(mol-K), the
  only place this conversion happens), then
  `compute_isentropes_two_phase(d1,d2,bubble_rows,dew_rows,s_values_jmolK)`,
  `compute_isentrope_liquid_side(d1,d2,z1,bubble_rows,s_values_jmolK)`,
  `compute_isentrope_vapor_side(d1,d2,z1,dew_rows,s_values_jmolK)`, all
  three results passed into both `plot_envelope(...)` and
  `plot_envelope_honeywell_units(...)` calls alongside the existing
  isotherm dicts.
- Verified via `python3 -m py_compile mixture_isentrope_validation.py`
  (clean) and grep cross-checking every compute-function definition
  against its `_cli()` call site for matching argument order/count --
  all three consistent. NOT yet run end-to-end (sandbox lacks scipy) --
  next step: user runs the file and checks the 14 isentropes render in
  blue, roughly matching the diagonal shape/slope-changes in the
  Honeywell reference image (steeper in subcooled liquid, bending
  through the two-phase interior, shallower in superheated vapor), with
  the corrected direct-entropy calculation now in place.
- Isentrope validation milestone: implementation complete, empirically
  UNVERIFIED (no scipy in this sandbox). Isentrope work marked
  in-progress -> ready-for-user-run in the task list.

### 2026-08-13 (later still) — Isentrope run results: visually off vs Honeywell; two candidate causes identified, NEITHER fixed yet

- User ran the file and shared both figures. Assessment: "We are
  cimplerely off I think" (completely off). Symptoms: 0.22/0.26 Btu/lb-R
  land roughly where Honeywell's do (far left, near-vertical, well
  separated), but 0.28 through 0.37 are crushed together in a narrow
  band instead of fanned out across the liquid region like the
  reference, and several labels sit deep inside the dome interior
  instead of near the chart top -- meaning those lines terminate early.
- First diagnosis (solver-side): `compute_isentrope_liquid_side` /
  `compute_isentrope_vapor_side` use `scipy.optimize.root` (method=
  "hybr") with NO bounds on either unknown (T,rho). Nothing constrains
  the solve to stay on the liquid (or vapor) branch -- when the anchor
  row's own entropy isn't close to the target, the solver can converge
  to a completely different, unphysical, or wrong-phase root instead of
  the intended one. Proposed fix (NOT yet applied): switch to
  `scipy.optimize.least_squares` with explicit bounds (rho > rho_l_sat(T)
  for liquid, rho < rho_v_sat(T) for vapor), same style as the existing
  bubble/dew solve's sigmoid-bounded density parameterization.
- Second, more fundamentally important correction from the user:
  "None of the isentropes in Honeywell plot cross into two phase. Ours
  seem to land in the two-phase region." This directly contradicts the
  earlier "isentropes DO cross the two-phase interior, standard chart
  convention" answer given a few exchanges prior -- for THIS specific
  set of 14 target values on the REAL Honeywell chart, apparently none
  of them ever intersect the dome at all (all fully single-phase across
  their entire span). Diagnosis: this points at a REFERENCE-STATE
  OFFSET between our EOS's entropy and Honeywell's own convention
  (Honeywell's chart footnote: "h=200kJ/kg, s=1.00kJ/kg-K, sat. liq. at
  0degC") -- the SAME category of caveat already flagged for enthalpy
  (`to_honeywell_units`'s docstring: "Does NOT reconcile the Honeywell
  chart's reference-state footnote... shape/width comparisons are
  valid, absolute enthalpy position may be offset") but never actually
  checked for entropy. If our EOS's native entropy scale has an additive
  offset relative to Honeywell's, using Honeywell's LITERAL isentrope
  values (0.22-0.49 Btu/lb-R) as our targets will land some of them
  inside our computed dome even if the true physical states they
  represent never touch Honeywell's own dome -- not necessarily a bug
  in the two-phase crossing logic itself.
  Proposed next step (NOT yet done): compute our EOS's saturated-liquid
  entropy at a convenient reference point (e.g. sat. liquid at 32F/0C)
  and compare against Honeywell's reference anchor (s=1.00 kJ/kg-K at
  sat. liq., 0C) to check for/quantify an offset, before deciding
  whether an entropy reference-shift belongs in the final-step-only
  conversion (alongside `to_honeywell_units`) or whether the solver
  bounding fix above is the primary issue.
- NEITHER fix has been applied yet -- both are diagnosed but open.
  Isentrope work status: BLOCKED pending a decision on which to tackle
  first (or both). This is the direct next task for a future session.

### 2026-08-13 (end of session) — Wrapping up: breadcrumbs updated, user will push to DowlingLab git themselves

- User: "I think We will update the breadcrumbs and push the code to
  dowlinglab git. I will push them." Per the long-standing instruction
  already on record in this file ("We should not be pushing any
  breadcrumb markdowns"), `PROJECT_CONTEXT.md` itself is NOT part of
  what gets committed/pushed -- only the actual code files
  (`mixture_isentrope_validation.py`, `mixture_isotherm_validation.py`,
  `Honeywell_T_P_revalidation.py`, etc.) go to git, and the user pushes
  from their own terminal (this sandbox has no outbound git/network
  access, confirmed multiple times in earlier sessions).
- User also asked about "a markdown file that we uploaded" -- checked
  this session's uploads folder directly: no `.md` file present (only
  PNGs -- chart screenshots and the Honeywell reference image -- two
  PDFs -- the Honeywell Solstice N15 TDS and a journal article
  `1-s2.0-S0378381216305349-main.pdf` -- and a few VLE-run CSVs/JSON).
  Flagged the candidate `.md` files that DO exist in the repo itself
  (`README.md`, `R515B_props_validated/README.md`,
  `verification/saturation_review.md`,
  `verification/saturation_strict_fix_summary.md`,
  `verification/plr_brainstorm_manual_checks_2026-03-04.md`,
  `diagnostics/pseudopure_iso/T_70C_pseudopure_isotherm_table_20260310.md`)
  in case one of those is what the user meant -- user confirmed: "It is
  called readme.md." Initially resolved to the top-level `README.md`
  (repo root, confirms "dowlinglab" = Prof. Alexander Dowling's group at
  Notre Dame, NSF ERC EARTH hub) -- but user then shared a git file
  listing screenshot showing the ACTUAL intended file: the `README.md`
  inside `R515B_props_validated/` (last commit "Rename
  R515B_final_validation to R515B_props_validated," yesterday --
  "We uploaded this yesterday").
- Read `R515B_props_validated/README.md` in full. Contents: documents
  the same ~2% liquid-density extrapolation bias already established
  earlier this session (Bell 2023 departure function fit range
  x1=0.33-0.68 vs R515B's real x1~=0.9385), the same 3 bug fixes in
  `pressure_validated_model.py` (beta_v 0.99290->0.999290, missing /8.0
  in the nu12 mixing rule, the t1-vs-per-term-exponent bug), and the
  same "self-consistent solve is protected, externally-supplied density
  is not" conclusion. NEW information not previously in this breadcrumb:
  an independent confirmation of this exact density-extrapolation gap on
  NIST's own REFPROP issue tracker --
  [`usnistgov/REFPROP-issues#750`, "R-515B Density Uncertainty"](https://github.com/usnistgov/REFPROP-issues/issues/750).
  Another user raised the same composition-range argument there; a NIST
  REFPROP maintainer (`marciahuber`) responded pointing to an updated
  `HMX.BNC` file and a dedicated `R515B.MIX` file with improved
  interaction parameters, available on request from REFPROP@NIST.GOV but
  NOT yet in the public release, and NOT independently confirmed accurate
  by anyone in that thread -- best available lead, not a confirmed fix.
  Worth revisiting if the density-bias issue (and its downstream
  isentrope reference-state-offset hypothesis, see above) needs a deeper
  fix than what's achievable with the current Bell (2023) departure
  function alone.

### SESSION SUMMARY (2026-08-13) -- for quick orientation next session

This entire long session progressively rebuilt the Honeywell Solstice
N15 (R-515B) p-h chart feature-by-feature, all in
`mixture_isotherm_validation.py` and (most recently) its copy
`mixture_isentrope_validation.py`:
1. Honeywell 86-point P-T revalidation (`Honeywell_T_P_revalidation.py`) -- CLOSED, 1.84% MAPE.
2. IP-unit (psia/Btu-lbm) conversion architecture, final-step-only, validation-only -- CLOSED.
3. Quality lines (x=0.1-0.9) -- CLOSED, validated against the reference chart.
4. Isotherms: two-phase (red, horizontal segments) -- CLOSED. Liquid-side subcooled extension (red, up to 1350 psia) -- CLOSED, smoothed via geomspace. Vapor-side superheated extension (red, down to 15 psia) -- CLOSED. Supercritical isotherms 240-400F (red, no dome anchor, ideal-gas-seeded) -- CLOSED.
5. Chart polish -- CLOSED: box border, Honeywell-exact axis labels/tick values (psia/Btu-lbm plot only), annotation_clip=False fix (matplotlib silently drops off-view-limit labels).
6. Isentropes (blue, 14 Honeywell values 0.22-0.49 Btu/lb-R) -- IN PROGRESS / BLOCKED. Two-phase segment (algebraic lever-rule) + liquid/vapor single-phase segments (2D scipy.optimize.root walk) implemented and wired, entropy calc corrected to the DIRECT Helmholtz identity (`_mix_entropy_direct`, matching `mixture_model_one_point.py`'s `compute_table1_properties` gold reference) per explicit user correction -- but the RENDERED result is visually wrong vs Honeywell (lines crushed together, wrong dome-crossing behavior). Two open, unfixed hypotheses: (a) unbounded `scipy.optimize.root` jumping to wrong-branch roots -- fix: switch to bounded `least_squares`; (b) an entropy reference-state offset vs Honeywell's own convention (same unresolved caveat as enthalpy) -- fix: check/quantify offset at a reference point (sat. liquid, 32F) before deciding whether to add an entropy reference shift to the final-step-only conversion. NEXT SESSION SHOULD START HERE.
7. Still fully unscoped: green density/specific-volume lines, and the isentrope superheated-region continuation's visual correctness (blocked on the above).
LESSON reconfirmed this session: always re-grep/re-read on-disk state before trusting prior edits persisted, since the user edits these files concurrently in VS Code.
