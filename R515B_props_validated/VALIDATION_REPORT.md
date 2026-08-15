# R-515B Mixture Model — Validation Report: Known Errors and Limitations

**Scope:** Helmholtz-EOS mixture model for R-1234ze(E)/R-227ea (Bell 2023 departure
function, corresponding-states mixing rules), evaluated at mass fraction w1=0.911
(mole fraction x1≈0.9385) to approximate the commercial blend R-515B (Honeywell
Solstice N15). Validated against Honeywell's own Solstice N15 Technical Data Sheet
(TDS): its physical-properties table and its p-H chart.

**Status as of 2026-08-14:** the model is being called validated. This document
records every known error and limitation found during that validation, so the model
is used with its actual accuracy envelope in mind rather than an assumed one. Full
session-by-session diagnostic trace lives in `PROJECT_CONTEXT.md`; this document is
the distilled, audit-ready summary.

**Important scope boundary:** none of the findings below reach the cycle simulation
models in this repository (`vapor_compression*.py` family). Those use IDAES's native
pure-component Helmholtz package plus CoolProp and do not import any mixture-model
module discussed here. Everything in this report is confined to the standalone
mixture-property validation track (`R515B_props_validated/` and related scripts).

**Checked directly against `DVCT_Project` (2026-08-14):** `vapor_compression.py`
there builds its property package via `HelmholtzParameterBlock(pure_component=
fluid_name, ...)` -- IDAES's own native pure-component engine -- and does not import
`R515B_props_validated/`, `linear_model_codex.py`, or any file this report covers.
So today, none of the above reaches DVCT. However,
`DVCT_Project/property_packages_DVRT_code/code_linear_model.py` is an early,
incomplete draft explicitly intended to build R515A/R515B mixture properties for
future use in a cycle (it currently errors if run -- references an undefined
`mixture_data` and a misspelled helper function -- and nothing imports it yet). If
that file is finished and wired into a cycle later, everything in this report
becomes directly relevant at that point, especially Section 2.2 (entropy), since
isentropic compressor/valve modeling depends heavily on entropy. Treat that
integration as a trigger to revisit this report.

**Specific mechanism for COP, confirmed against IDAES's own compressor code
(`idaes/models/unit_models/pressure_changer.py`, ~line 798-817):** the compressor's
ideal discharge state is found by constraining `properties_isentropic.entr_mol ==
properties_in.entr_mol` at the discharge pressure -- literally an isentrope point
solve, the same computation validated throughout this report. Actual compressor
work (and therefore COP = `evaporator.heat_duty / compressor.work_mechanical`) is
derived directly from that isentropic enthalpy rise. So once R515B mixture
properties feed this block, COP inherits Section 2.2's entropy error directly and
structurally, not incidentally -- and specifically at the DISCHARGE-pressure/
condenser-side state, which for many operating envelopes sits closer to critical
(where the entropy error is largest) than the evaporator-side inlet state does.
Today, with `pure_component=fluid_name` feeding IDAES's own native engine, this
constraint is unaffected.

---

## 1. Root cause underlying nearly everything below

**The Bell (2023) R-1234ze(E)/R-227ea departure function was fit and validated only
for mole fraction x1 = 0.33–0.68.** Real R-515B sits at x1≈0.9385 — well outside that
range. This was confirmed to be a genuine extrapolation limitation, not a coding
error: every candidate coding bug (parameter transcription, reducing-function
formulas, departure-function formula, corresponding-states chain rule) was
independently re-verified correct against the primary source before this conclusion
was reached.

Everything in Section 2 is a downstream consequence of this one fact, at varying
magnitude depending on how directly each property depends on the raw departure
surface versus its derivatives (see 2.4 for why entropy is hit hardest).

---

## 2. Current open errors / limitations (not fixed — the model's real accuracy envelope)

### 2.1 Deep-liquid and critical density error
- Deep-liquid density (T=298.15K, away from critical): **~2% error**, amplified by
  extreme dP/drho sensitivity — force-feeding the literal biased density into a
  pressure calculation produces up to a ~9.7x pressure error. The production VLE
  solver does NOT do this (it self-consistently solves for its own density), which
  is why bubble/dew **pressure** predictions still validate well (~1–1.1% MAPE)
  despite this density error underneath.
- Critical density: **-8.8%** (3858.4 vs Honeywell's 4230.8 mol/m³) — much larger than
  critical temperature (-0.04%, essentially exact) or critical pressure (-0.57%,
  small). The critical-point inversion (dP/drho=0, d²P/drho²=0) is far more sensitive
  to departure-surface error than the simple algebraic Tc combining rule.
- Away from critical, liquid density is consistently **under**-predicted and vapor
  density consistently **over**-predicted (opposite-sign, structured bias, not
  noise) — both in the same ~2–4% band as pressure at those conditions (i.e. no
  special "protection" for pressure away from the near-critical region).

### 2.2 Entropy — the most affected property, growing error approaching critical
- At Honeywell's stated reference state (T=273.15K, sat. liquid, 0°C): raw model
  entropy is **+2.5% high** (1.012 vs Honeywell's defined 1.00 kJ/kg·K). A flat,
  T-independent offset (`ENTROPY_REFERENCE_OFFSET_JMOLK = -1.4595` J/mol·K) is
  applied and corrects this exactly at that one calibration point.
- Away from the calibration temperature, the residual entropy error **grows**:
  liquid-side isentropes 0.32/0.34/0.35 Btu/lb-R show a positional drift of roughly
  **50 to 100+ psia low** relative to Honeywell's chart (visual/pixel estimate, not
  exact), increasing with entropy value (i.e., worsening toward the critical
  region).
- **An explicit T-dependent (quadratic) correction was developed, fitted against 15
  hand-digitized chart points, and wired in — then abandoned entirely** after it
  introduced a real step-discontinuity (-1.84 J/mol·K) at the liquid/vapor branch
  switch, breaking solver convergence for 3 isentropes. The codebase currently
  retains only the flat offset above; **the growing near-critical entropy drift is
  therefore a known, currently-uncorrected limitation.**
- Structural explanation for why entropy specifically is hit hardest: pressure and
  enthalpy are built from *derivatives* of the departure Helmholtz function
  (`ar_del`, `alpha_tau`); entropy is the only one of the three built directly from
  the raw, undifferentiated departure value itself (`S/R = tau*alpha_tau - alpha`).
  A roughly constant "level" bias in the departure surface is structurally
  insulated in P and H (a constant offset contributes ~zero to its own derivative)
  but passes directly into S with coefficient -R. At the reference state this
  splits roughly 50/50 between a density-mediated share and an independent
  formula-only share; that split was not re-checked near critical, where it may
  shift further toward entropy.

### 2.3 0.37 Btu/lb-R isentrope does not visually touch the critical point
- On Honeywell's chart, 0.37 Btu/lb-R (liquid-side only — below the model's own
  critical entropy) visually touches the dome tip exactly at the critical point. In
  the model's rendering it approaches but stops visibly short.
- Verified this is **not** a fallback-anchor bug: the anchor point matches an
  independently-interpolated, much-finer-grid true bubble-line point to ~0.02% in
  enthalpy. But that true point's own enthalpy is still **~2.9% below the critical
  point's own enthalpy**, despite being within ~0.5% of critical in pressure and
  ~0.1% in temperature (~0.25K from Tc) — consistent with genuine, steep near-critical
  divergence of dh/ds as T→Tc, not a computational bug.
- **Decision (2026-08-14, this session): reported to the user as a known model
  limitation rather than pursued further.**

### 2.4 0.39 Btu/lb-R isentrope — substantially improved, still incomplete
- Fixed this session: the vapor side now gets a genuine 47-point extension
  (477.6–1350 psia, all independently verified outside the two-phase dome) by
  walking upward from its near-critical anchor toward and past the critical point,
  instead of the downward direction every earlier attempt used.
- **Still missing:** any downward/low-pressure vapor segment, and any liquid-side
  segment. Both were checked multiple independent ways (adaptive pressure stepping,
  temperature-stepped continuation, biased Newton seeding) and confirmed to be
  **genuinely two-phase** in that direction for this specific entropy value, not a
  numerical artifact — i.e., not expected to be fixable with the current
  formulation without a fundamentally different approach (e.g. solving the interior
  two-phase isentrope path directly, which is currently disabled by design, see 2.5).

### 2.5 No isentrope segments are drawn inside the two-phase dome
- By deliberate design (matching Honeywell's own chart convention), the model never
  draws an isentrope through the interior of the two-phase region — only the
  liquid-side and vapor-side single-phase segments. An earlier attempt to compute
  these interior segments (`compute_isentropes_two_phase()`) was found to trace
  constant-*quality* lines, not constant-*entropy* lines (a conceptual bug, not
  fixed), and the function has been disabled ever since. This is a scope boundary
  more than an "error," but it means isentropes very close to critical entropy
  (0.37, 0.39) can only be drawn where they leave the dome, never through it —
  directly relevant to 2.3 and 2.4 above.

### 2.6 Critical-point solver has mild grid-dependence
- The solved critical temperature varies by **~0.37K (~0.1%)** depending on the
  input temperature grid's range/resolution (381.894K on the standard grid vs.
  381.522K on a grid pushed closer to Tc). Not deeply investigated; a secondary
  contributor to the residual uncertainty in 2.3's near-critical numbers.

### 2.7 Solver seed robustness (usage caveat, not a model-accuracy error)
- The production VLE pipeline (temperature continuation, each point seeding the
  next) is robust and unaffected by this. But any *new, single-shot* solver call
  outside that established continuation — e.g. a quick one-off diagnostic — uses a
  generic composition-agnostic seed that was shown to converge to spurious
  (non-physical) roots at most compositions away from the real R-515B composition,
  while still reporting `status=CONVERGED`. Anyone extending this codebase with new
  single-point solves should seed physically, not assume convergence implies the
  physical root.

---

## 3. Errors found and fixed during validation (for the audit trail)

These no longer affect the current model but are recorded because they were real
bugs, not just extrapolation limitations, and because their fixes are part of what
makes the current state "validated":

| # | Bug | Fix | Magnitude |
|---|-----|-----|-----------|
| 1 | `beta_v` transcription typo (0.99290 vs. correct 0.999290, Bell 2023 Table 2) | Corrected to 0.999290 | Negligible (~0.0005 kg/m³ shift) |
| 2 | Missing `/8.0` divisor in reducing-volume mixing formula (`pressure_validated_model.py` only) | Added divisor | Contributed to a ~92x pressure error in that file |
| 3 | Hardcoded first-term temperature exponent (`t1`) reused for every residual term (`pressure_validated_model.py`) | Index exponent per-term | Root cause of the ~92x pressure error above; fixed to ~0.002% agreement with the validated pipeline |
| 4 | Two-phase isentrope function traced constant-*quality* lines, not constant-*entropy* lines | Disabled the function entirely (deliberate, permanent) | Every 2-phase isentrope segment was structurally wrong in shape |
| 5 | Vapor-side isentrope anchor collision: 5 isentropes (0.41–0.49) all shared one anchor (the dew branch's own entropy peak), each requiring one large, ill-conditioned Newton jump | Incremental entropy ramp (small isobaric substeps) + upward pressure walk + dome-reentry guard | Entropy gaps of 5.5–44.8 J/mol·K bridged in one bad step; now 79 clean points each, full 15–1350 psia range |
| 6 | Isentrope Newton residual mixed Pa-scale and J/(mol·K)-scale terms unnormalized, breaking the solver's Jacobian conditioning | Normalized to relative/dimensionless residuals | Was a prerequisite fix underlying essentially all later isentrope work |
| 7 | Entropy reference-state offset vs. Honeywell's defined convention | Flat correction, `ENTROPY_REFERENCE_OFFSET_JMOLK = -1.4595` J/mol·K | +2.5% at reference state, corrected to ~exact there (see 2.2 for what remains uncorrected away from that point) |
| 8 | Dome kink/elbow near critical in rendered plots | Not a code bug — caused by choosing `--Tmax` well below the mixture's real Tc, leaving ~30K of unsolved grid near the dome tip | Raise `--Tmax` close to (but below) the real Tc for the given composition |
| 9 | Honeywell TDS pressure-unit inconsistency (table uses psig, chart figures use psia) | Documented; chart-read values are NOT given the +14.696 psi gauge→absolute conversion (table values still are) | Would be a 14.696 psi error if mishandled |

---

## 4. Summary judgment

Pressure and enthalpy validate well nearly everywhere (~1–2% typical, worse only very
close to critical). Entropy is the weakest property throughout, with a real and
currently uncorrected drift that grows approaching the critical point — this is the
most important caveat for anyone using this model's entropy output, isentrope
placement, or anything derived from S near critical. The two isentropes closest to
critical entropy (0.37 liquid-only, 0.39 straddling) are each individually
documented above as partial or visibly imperfect near the dome tip specifically
because that is where the underlying composition-extrapolation limitation is most
exposed. Away from the critical region and away from deep-liquid conditions, the
model performs close to Bell (2023)'s own stated accuracy for its fitted composition
range.
