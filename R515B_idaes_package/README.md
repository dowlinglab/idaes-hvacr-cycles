# R515B_idaes_package

A new, independent IDAES `PhysicalParameterBlock`/`StateBlockData` property
package for R-515B (R-1234ze(E)/R-227ea, w1=0.911), built to be a drop-in
replacement for `HelmholtzParameterBlock` inside `vapor_compression.py`-style
flowsheets. This is a **different architecture** from the reference/validation
work in `R515B_props_validated/`: that folder holds the plain-Python/SciPy
"oracle" model (point queries only, no Pyomo NLP participation); this folder
reproduces the same validated Bell (2023) physics as native Pyomo
expressions/constraints so R-515B can actually be solved inside a real IDAES
flowsheet by IPOPT, the same way R134a/R1234ze(E)/R227ea already are via
`general_helmholtz`.

The build was carried out as a staged "MASTER TASK" (Stages A through O),
with every numeric result, tolerance, and design decision recorded in
`helmholtz_prop_validation.md` (repo root) -- that file, not this README, is
the authoritative validation record. This README is the file-by-file map.

## Status

All 27 stages (0-27) in `helmholtz_prop_validation.md` are complete: kernel
validation, single-phase/bubble/dew/critical-point validation, native-Pyomo
EOS/flash/saturation validation, `DiagnosticsToolbox` structural checks
(Stage M), and a full end-to-end `vapor_compression_r515b_integration.py`
flowsheet solve (Stage N) all pass against the oracle
(`R515B_props_validated/mixture_fully_validated.py`, confirmed functionally
equivalent to the actively-worked `mixture_isentrope_validation.py`). A
downstream COP-vs-ambient comparison against R134a was also run, both as an
ad hoc sweep (`cop_sweep_r515b_vs_r134a.py`) and, later, against the
project's established PLR benchmark conditions (`vapor_compression_plr_
r515b.py`).

## Files

### Core production modules (meant to be pushed)

**`r515b_helmholtz_core.py`** -- Stage F-K production core: plain NumPy/SciPy
reimplementation of the oracle's validated physics (reducing functions,
Helmholtz kernel derivatives, bubble/dew/critical-point/isentrope solvers,
pseudo-pure saturation), independent of the oracle at runtime (imports only
from the pre-existing `linear_model_codex.py`, not from any oracle file).
Every function documents which oracle function it reproduces and how.

**`r515b_pyomo_eos.py`** -- Stage L (part 1): the same EOS kernel translated
into native Pyomo expressions (`pyomo.environ.exp`/`log`, real Pyomo `Var`s)
so IPOPT's automatic differentiation can solve it as part of the flowsheet's
equation-oriented NLP, rather than calling a black-box SciPy function from
inside a Constraint body. Includes the native-Pyomo smooth PH-flash system
(`mixture_ph_flash_residuals_expr`) and saturation-curve system
(`saturation_residuals_expr`).

**`r515b_property_package.py`** -- Stage L (part 3, final): the actual IDAES
`R515BParameterBlock`/`R515BStateBlockData` classes wiring the Pyomo EOS
kernel and PH-flash system into a real, drop-in-compatible IDAES property
package (`StateVars.PH` + `AmountBasis.MASS`, matching the exact 3-state-Var
contract -- `flow_mass`, `pressure`, `enth_mass` -- that `general_helmholtz`'s
`HelmholtzStateBlockData` already exposes, confirmed by direct inspection of
the installed IDAES 2.12.0 source).

**`ancillary_initial_guess.py`** -- Fast, non-iterative initial-guess
machinery for R-515B saturated liquid/vapor molar density, built by composing
the pure-component ancillary saturation-density correlations already shipped
in IDAES's own `r1234ze.json`/`r227ea.json` through the project's existing
Bell (2023) reducing-function mapping. Used only to seed `initialize()`'s
cold solve (no "previous point" to continue from inside a flowsheet); never
used in place of a real iterative solve.

**`vapor_compression_r515b_integration.py`** -- Stage N: a new, self-contained
integration copy of the vapor-compression-cycle flowsheet
(`vapor_compression.py`, kept strictly read-only/unmodified per the task
spec), substituting `R515BParameterBlock` for `HelmholtzParameterBlock`.
Hardcodes `Mode.PH` behavior only (the only mode `R515BParameterBlock`
implements) and adapts `specify_initial_conditions()`/state-var fixing for
R-515B's fixed-composition-blend nature.

**`vapor_compression_plr_r515b.py`** -- PLR-enabled R-515B cycle, built the
same way as `R1234yf/vapor_compression_plr_r1234yf.py`: a copy of
`vapor_compression_r515b_integration.py` with the same additive PLR/CD diff
(`_compute_plf`, `fs.plr`/`fs.cd` Params, `get_full_load_cop()`/
`get_part_load_cop()`, plus the superheat/subcool floor+cap fix) grafted on,
so R-515B can run through the project's established PLR COP-vs-ambient
benchmark (PLR=0.75, CD=0.13, evap_sat=-29C, condenser approach=9C,
superheat/subcool=3C, compressor_efficiency=0.75) alongside
R134a/R1234ze(E)/R1234yf. Adds `_sat_at_t_robust`, a continuation-seeded
wrapper needed because -29C (244.15K) sits below Stage L's originally
validated 255-375K saturation grid.

**`establish_reference_tolerances.py`** -- Stage C: measures the oracle's own
run-to-run reproducibility at deterministic representative states (liquid,
vapor, two-phase) before any new code is judged against it, and uses the
measured spread to set frozen regression tolerances with a safety margin.

**`cop_sweep_r515b_vs_r134a.py`** -- Downstream analysis (not part of the
Stage A-O build): sweeps R-515B and R134a through the same simplified
vapor-compression cycle over evaporating temperature at a fixed 40degC
condensing temperature, to compare converged COP directly. Produces
`cop_sweep_r515b_vs_r134a.json`/`_results.json`/`.html` (all generated
outputs, gitignored -- see below). Superseded for benchmark purposes by
`vapor_compression_plr_r515b.py`'s PLR-condition sweep, but kept as the
original ad hoc comparison.

### Validation scripts (`validate_*.py`, meant to be pushed)

Each validates one stage of the build against either the oracle directly or
the previous, already-validated layer beneath it (kernel -> Pyomo EOS -> Pyomo
flash/saturation -> full state block -> full flowsheet):

- `validate_table1_vs_oracle.py` -- Stage D/F/G: `linear_model_codex.py`'s
  `compute_table1_properties` vs. the oracle's `mix_state`/entropy.
- `validate_fugacity_vs_oracle.py` -- Stage H: fugacity/chemical-potential
  outputs vs. the oracle's `chemical_potentials_analytic`.
- `validate_core_vs_oracle.py` -- Stage F/G/H: `r515b_helmholtz_core.py` vs.
  the oracle at bit-identical input states.
- `validate_quality_isotherm_vs_oracle.py` -- Stage K (part 1): quality-line
  and isotherm machinery vs. the oracle.
- `validate_isentrope_vs_oracle.py` -- Stage K (part 2): isentrope machinery
  vs. `mixture_isentrope_validation.py` specifically.
- `validate_pseudopure_saturation.py` -- Stage L support: sanity-checks the
  new pseudo-pure saturation solve against the rigorous bubble/dew branches
  (not a pass/fail oracle comparison -- a deliberate simplification check).
- `validate_ancillary_guess.py` -- validates `ancillary_initial_guess.py`'s
  composite guess is seed-quality only, not a substitute saturation model.
- `validate_pyomo_eos_vs_core.py` -- Stage L (part 1): native Pyomo EOS-kernel
  expressions vs. the already-validated NumPy core.
- `validate_pyomo_state_vs_core.py` -- Stage L (part 2): native Pyomo P/h/s/g/Z
  and mu1/mu2 expressions vs. the NumPy core.
- `validate_pyomo_saturation_vs_core.py` -- Stage L (part 3): native Pyomo
  saturation-curve system, solved as a genuine implicit NLP, vs. the SciPy
  pseudo-pure saturation solve.
- `validate_pyomo_flash_vs_core.py` -- Stage L (part 3): native-Pyomo smooth
  PH-flash system vs. the explicit-branching SciPy reference.
- `validate_entropy_and_tsat_vs_core.py` -- validates the entropy/T_sat
  additions to the PH-flash system and `R515BStateBlockData` against the
  SciPy reference.
- `validate_state_block_construction.py` -- Stage L (part 3, final):
  `R515BParameterBlock`/`R515BStateBlock` construction, degrees-of-freedom,
  and a real IPOPT solve inside an actual `FlowsheetBlock`.
- `validate_stage_m_diagnostics.py` -- Stage M: `DiagnosticsToolbox`
  structural/numerical checks on the property package.
- `validate_stage_n_integration.py` -- Stage N: end-to-end exercise of
  `R515BVaporCompressionCycle` (construction, initialization, solve, COP
  sanity/consistency checks).

## Generated outputs (gitignored, not pushed)

`cop_sweep_r515b_vs_r134a.py` is the only generated-output-producing script
in this folder whose script itself is pushed; its outputs are not:
`ancillary_guess_validation_results.json`, `cop_sweep_r515b_vs_r134a.json`,
`cop_sweep_r515b_vs_r134a_results.json`, `cop_sweep_r515b_vs_r134a.html`,
`cop_vs_ambient_plr_r515b_vs_r134a.html`, and `reference_repeatability_
results.json` (from `establish_reference_tolerances.py`) are all
run-generated result/plot files, reproducible by re-running the corresponding
script, and are excluded via `.gitignore`. `__pycache__/` is likewise
excluded.

## Known limitations

- **No TPX state-var support.** `R515BParameterBlock` implements only
  `StateVars.PH` + `AmountBasis.MASS` (matching `vapor_compression.py`'s own
  usage); the ORIGINAL_TPX/IMPROVED_TPX branches are omitted entirely.
- **-29C benchmark point required a continuation-seeded workaround**
  (`_sat_at_t_robust` in `vapor_compression_plr_r515b.py`) because it sits
  below the originally validated 255-375K saturation grid; the underlying
  solver is fragile in that range with its default retry seeds.
- Inherits the same entropy caveat as the oracle it's built from (see
  `R515B_props_validated/README.md`): entropy is the weakest validated
  property, with a bias that grows approaching the critical point.

## Full validation record

See `helmholtz_prop_validation.md` (repo root) for the complete, staged,
numeric validation record (Sections 0-27) underlying every file above.
