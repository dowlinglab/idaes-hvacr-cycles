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

**Do we actually anticipate operating near critical? Tested directly (2026-08-14)
against realistic R-515B operating conditions, not just argued abstractly:**
evaporator-exit saturated-vapor entropy across the entire practical evaporating
range (-15C to +10C, i.e. essentially all of R-515B's target medium-temperature
refrigeration/chiller/AC envelope) sits nearly flat at **0.3955-0.3959 Btu/lb-R** --
almost exactly the band straddling isentropes 0.39 (partially fixed this session)
and 0.41 (fully fixed). This is NOT a rare edge case; it is the default compressor-
inlet entropy for essentially any realistic evaporating temperature for this fluid.

Ran the actual isentropic-discharge calculation (the same one IDAES's compressor
model performs) for 4 realistic (T_evap, T_cond) pairs: (-10C,40C), (0C,40C),
(5C,45C), (-10C,50C). **With saturated vapor (no superheat) at the compressor
inlet, all 4 cases land inside the two-phase dome at the condenser pressure instead
of reaching superheated vapor** -- confirmed not a fluke: the saturated-vapor
entropy at the evaporating temperature (194.5-194.65 J/mol*K) is consistently LOWER
than the saturated-vapor entropy at the condensing temperature (195.6-195.9
J/mol*K) in this model, the classic signature of "wet compression." This may be a
genuine physical trait of this refrigerant blend (HFO-based blends like R-1234ze
are known to have flatter/wetter vapor domes than older HFCs) rather than purely a
model artifact -- but it is happening in exactly the region flagged as least
trustworthy in this report (Section 2.2), so the magnitude should not be taken at
face value without independent confirmation.

**Practically resolved by realistic superheat:** re-ran the same 4 cases with 5K of
compressor-inlet superheat (a small, standard amount -- real systems typically run
5-15K, specifically to prevent liquid slugging). All 4 succeed cleanly, staying
outside the dome entirely, with entropy shifting to ~0.399 Btu/lb-R and higher --
into the range already confirmed fully fixed. **Conclusion: the near-critical/
wet-compression fragility is real and would be hit by any cycle model that assumes
idealized saturated vapor at the compressor inlet, but is avoided by including even
modest, realistic superheat.** Any DVCT cycle representation of this refrigerant
should model compressor-inlet superheat explicitly rather than assuming saturated
vapor, both because real systems do this and because it sidesteps this exact
fragility.

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
- **Independently confirmed against real experimental data, not just a single
  spot-check** — see Section 3 below: 67 literature (T,P,ρ) points spanning
  T=254–362K, P=0.87–12.27MPa give a strikingly consistent -1.7% to -2.5% model
  bias (AAD%=1.87%), essentially the same systematic ~2% underprediction found
  earlier from a single point, now backed by a much larger independent dataset.
- **Root cause of the ~2% mixture bias narrowed down further (Section 3.2):** the
  SAME literature source's pure-component data shows both R-1234ze(E) and R227ea
  individually accurate to a few tenths of a percent (0.22% and 0.32% AAD
  respectively) — 5-8x smaller than the mixture's 1.87%. This directly confirms
  (not just infers) that the bias lives in the Bell (2023) mixing/departure-function
  layer, not the underlying pure-fluid EOS.

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

## 3. Independent literature validation: Kang et al. (2024) experimental density

### 3.1 Mixture density (R-515B, Table 4)

**First comparison in this project against data that is not the Honeywell TDS chart.**
Kang, K., Yang, S., Cui, J., Gu, Y. (2024), "Theoretical study and experimental
verification of the viscosities of azeotropic refrigerant R515B," *International
Journal of Refrigeration* 168, 59–69 (https://doi.org/10.1016/j.ijrefrig.2024.08.012)
reports real vibrating-wire measurements of R-515B liquid density and viscosity at
w1=91.1/8.9% R-1234ze(E)/R227ea (matching this project's composition exactly),
T=254.13–362.34K, P=0.87–12.27MPa, with declared combined expanded uncertainty of
0.2% (k=2) for density — a much stronger, quantitative check than reading points off
a chart by eye.

**Method:** the paper's Table 4 (68 points stated; 67 discrete rows could be cleanly
transcribed from the printed table image — see `kang2024_r515b_density_table4.csv`)
was compared point-by-point against this model's own predicted liquid density at
each exact (T,P), solved independently via `validate_against_kang2024.py` (new
script, this repo). Solving density at fixed (T,P) directly (rather than via the
model's own saturation/VLE solve) turned out to need real care: a first attempt used
`brentq` with a wide density bracket and got wildly wrong results (35–67% "error")
for many points — traced to the bracket spanning an unstable/unphysical root of the
Helmholtz EOS's P(ρ) curve (a Van der Waals-loop fold below the spinodal), not a
model error at all. A second attempt narrowed the bracket to ±40–60% of the
experimental density and still picked up the wrong root for several points nearer
the critical region. The fix that worked: a damped Newton iteration in molar density,
seeded exactly at the experimental value (already known accurate to 0.2%) and using
the model's own finite-difference `dP/dρ` (`_dp_drho_fd`, already used elsewhere in
this repo for the critical-point solve) — this follows the local slope from a
trusted seed and cannot jump to a distant, unrelated root the way a blind bracket
search can. All three attempts are preserved in the script's own comments/docstring
as a documented lesson, consistent with this report's practice of recording
methodology failures, not just final numbers.

**Result, after the fix — clean and highly consistent across all 67 points:**

| Metric | Value |
|---|---|
| N compared | 67 of 67 loaded (0 solver failures) |
| AAD% (mean absolute deviation) | **1.87%** |
| Mean signed bias | **-1.87%** (model reads low, i.e. underpredicts density) |
| Deviation range | -2.46% to -1.68% |
| Scatter | very low — essentially a flat systematic bias, not noise (see deviation plot) |

This independently confirms, with 67 real experimental points across the model's
full liquid T/P range (not a single hand-picked state), the same ~2% deep-liquid
density underprediction already documented in Section 2.1 from one spot-check at
T=298.15K. The consistency is notable: deviation stays within a ~0.8 percentage-point
band (-1.68% to -2.46%) across a 108K temperature range and an ~11.4 MPa pressure
range, with a mild trend toward smaller |deviation| at higher pressure. This is
strong evidence the ~2% liquid-density bias is a stable, systematic property of the
composition-extrapolation limitation (Section 1) — not scatter, not a fluke of one
test condition, and not something that gets worse within the ordinary liquid range
sampled here (it is specifically the near-critical region, Section 2.1's other
finding, where the density error grows much larger, to ~8.8%).

Full point-by-point results: `verification/kang2024_r515b_density_comparison.csv`.
Deviation plot (vs. both P and T): `verification/kang2024_r515b_density_deviation.png`.

### 3.2 Pure-component density (R-1234ze(E) and R227ea, Tables 3 and 2) — isolates WHERE the mixture bias comes from

The paper measured both pure components with the same lab, same vibrating-wire
method, same T/P grid as the R-515B mixture data above. This directly answers a
question flagged as open since 2026-08-12 ("Pure R-227ea's real saturated-liquid
density near 25°C still unconfirmed... isolate the Bell departure term directly")
— WHERE does the ~1.87% mixture bias actually come from: the pure-fluid EOS layer,
or the mixing/departure-function layer built on top of it?

Compared via `validate_pure_components_against_kang2024.py` (new script, same
damped-Newton method as 3.1, run at the pure-component composition limits z1=1.0
and z1=0.0 — confirmed beforehand that the mixing formulas evaluate cleanly at both
pure limits, no divide-by-zero):

| Fluid | N | AAD% | Mean bias | Range |
|---|---|---|---|---|
| R-1234ze(E) (pure, z1=1.0) | 67/67 | **0.225%** | -0.02% (~unbiased) | -1.24% to +0.30% |
| R227ea (pure, z1=0.0) | 68/68 | **0.323%** | +0.32% | -0.10% to +0.73% |
| *(for comparison) R-515B mixture* | *67/67* | *1.87%* | *-1.87%* | *-2.46% to -1.68%* |

**This settles the open question: the pure-component layer is NOT where the bias
comes from.** Both pure fluids are accurate to a few tenths of a percent — 5-8x
smaller than the mixture's ~1.87% — confirming what was only hypothesized before
(2026-08-12: "the remaining ~2% density error must be in the Bell (2023)
departure/mixing layer"). This is now a direct, quantitative confirmation rather
than an inference from ruling other things out. R-1234ze(E)'s result also
corroborates the project's very first pure-fluid check (a single point at 25°C,
matched to ~0.003%) — now backed by 67 points across the full range, with the
caveat that R-1234ze(E) shows more scatter than R227ea (a cluster of points in the
303-342K range reach -0.5% to -1.24%, worth a closer look if this matters for a
specific application, though the average remains essentially unbiased).

Full results: `verification/kang2024_r1234ze_density_comparison.csv`,
`verification/kang2024_r227ea_density_comparison.csv`. Deviation plot:
`verification/kang2024_pure_components_density_deviation.png`.

### 3.3 Critical point cross-check (Table 5, independently REFPROP-sourced)

The paper's Table 5 gives R-515B's critical properties as used in its own modeling,
explicitly sourced from REFPROP (Lemmon et al. 2018) — a second, independent source
beyond Honeywell's own TDS table, letting the earlier critical-point comparison
(Section 2.1) be cross-checked against two independent references instead of one.

| Quantity | Model | Paper (REFPROP) | Diff | Honeywell TDS (earlier finding) | Model vs. Honeywell |
|---|---|---|---|---|---|
| Molecular weight | 117.485 g/mol | 117.480 g/mol | **+0.004%** | — | — |
| Tc | 381.894 K | 382.030 K | **-0.036%** | 382.039 K | -0.038% |
| Pc | 3.5765 MPa | 3.5839 MPa | **-0.208%** | 3.5970 MPa | -0.570% |

Molecular weight matches almost exactly (trivial to get right, but a clean sanity
check that the mass-fraction-to-mole-fraction conversion is correct). Tc agrees with
both independent sources to within 0.04% either way — Honeywell's own table and
REFPROP agree with each other to within 0.002%, and the model sits right between
them. **Pc is a more interesting result: against REFPROP specifically, the model is
only -0.21% off, notably better than the -0.57% found earlier against Honeywell's
own stated value.** Honeywell's own table and REFPROP disagree with each other by
about 0.36% on Pc — so part of the originally-reported -0.57% gap reflects
disagreement between the two reference sources themselves, not solely the model's
own error. The model's true Pc accuracy is better characterized as ~0.2-0.6%
depending on which reference is treated as ground truth, not a single number.

Critical density (the -8.8% finding in Section 2.1) could NOT be cross-checked here
— Table 5 does not report it, and REFPROP's own value wasn't independently pulled
in this pass. Acentric factor (ω, also in Table 5) is a cubic/corresponding-states
EOS input; this model is a full multiparameter Helmholtz EOS and does not use or
report an acentric factor, so that column isn't comparable and is noted as not
applicable rather than forced into a comparison.

### 3.4 What else the paper has that isn't checkable against this model

The paper's Tables 6-8 (PC-SAFT parameters, residual-entropy-scaling viscosity
model coefficients, and viscosity AAD% comparisons) are all specific to the
authors' own separate viscosity-prediction model. This project's model computes
thermodynamic state properties (p, h, s, ρ) via a Helmholtz EOS and does not
compute viscosity at all — genuinely out of scope, not a gap in this validation.

---

## 4. Errors found and fixed during validation (for the audit trail)

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

## 5. Summary judgment

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
range. Liquid density now has a genuine independent literature check (Section 3,
Kang et al. 2024): the mixture's ~1.9% systematic underprediction is confirmed
across 67 real experimental points, and — importantly — the SAME data shows both
pure components individually accurate to a few tenths of a percent, pinning the
~2% mixture bias squarely on the mixing/departure-function layer rather than
leaving it as an inference. The same source's independently-REFPROP-sourced
critical properties also refine the critical-pressure comparison: against REFPROP
the model is only -0.21% off on Pc, better than the -0.57% found earlier against
Honeywell's own table, because Honeywell's and REFPROP's own stated values differ
from each other by ~0.36%. Critical density and the entropy-side findings remain
without an independent literature check as of this pass — the next-highest-value
target if further verification is wanted, since entropy is the property this
report already flags as least trustworthy.
