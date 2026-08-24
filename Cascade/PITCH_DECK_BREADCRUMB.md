# Pitch Deck Breadcrumb — EARTH-ERC Paper Pitch

Separate from `Cascade/BREADCRUMB.md` (which tracks the cascade *model*
development). This one tracks the paper-pitch deck build, so the model
breadcrumb doesn't get diluted with deck/figure work.

Gitignored, internal notes only — never reference this file in
user-facing docs or the deck itself.

---

## 2026-08-21 — Prior session data loss + recovery

A prior session that hit the org's monthly spend limit turned out to be a
fully separate, now-inaccessible cloud sandbox. Nothing it produced
survived (its own breadcrumb, figures, GWP/TEWI analysis, PFD schematic) —
none of it had been pushed to the connected device. Confirmed via
`device_list_dir` (recursive) and a grep of this session's own transcript
for pitch-related file paths — no matches, nothing recoverable.

Told the user plainly and asked for re-upload. They re-sent:
- `Dowling_Han_EARTH_Seed_Project_Proposal__07182025.pdf` (EARTH seed
  proposal — PI Dowling/Co-PI Han, defines funded scope: R515B/CO2
  cascade in DVRT/IDAES, Modelica dynamic sim, benchmark vs
  R134a-or-R513a/CO2)
- `Shilpa__August_2026_Meeting__Google_Docs.pdf` ("outline doc" — GWP
  table, Figure 1-4 layout plan, COP formula definitions, CAPEX/OPEX
  placeholder index numbers, explicit note that placeholder Fig 3b
  curves were NOT internally energy-consistent)
- `06_26_2026_Shilpa_group_meeting.pptx` ("example deck" — 16-slide
  reference for a different Pyomo.DoE paper pitch, used purely as the
  template/theme basis)

Rebuilt everything from these plus already-validated repo data (cascade
physics itself was intact — nothing needed re-deriving).

## 2026-08-21 — `pitch_deck/` working dir + figures

Working dir: `/home/claude/work/pitch_deck/` (cloud sandbox only until
pushed via SendUserFile + device_commit_files to
`Cascade/pitch_deck/` on the device).

- **`make_cascade_pfd.py`** → `fig_cascade_pfd.png`. PFD schematic,
  matplotlib `FancyBboxPatch`/`FancyArrowPatch`, populated with REAL
  stream data from `cascade_full_stream_comparison.csv` at ambient=25C
  (COP_full=2.4919, hot_mass_flow=1.898 kg/s). Hit a
  `ValueError: At least one value in the dash list must be positive`
  from a broken `ax.annotate(arrowprops=dict(lw=0, ls='--'))` call —
  removed it (was an unnecessary leftover connector).

- **`make_projected_cascade.py`** → `fig_projected_cascade.png`.
  Projected (not simulated) COP-vs-ambient curves for R515B/CO2 and
  R1234yf/CO2 cascades: validated R134a/CO2 cascade COP scaled by each
  fluid's own single-stage COP ratio to R134a at the same ambient.
  Explicitly labeled projected throughout — R515B/CO2 and R1234yf/CO2
  cascades have NOT been built/solved yet.

- **`make_capex_opex.py`** → `fig_capex_opex_cascade.png` +
  `capex_opex_results.csv`. CAPEX index = outline doc's placeholder
  single-stage table (R134a=1.00, R515B=1.16, R1234yf=1.34) × 1.82
  cascade uplift (outline doc's own R134a cascade/single-stage ratio,
  applied uniformly). OPEX index = 1/COP (real solved values,
  normalized to R134a single-stage=1.00) + stated O&M adder (R515B +5%,
  R1234yf +15%, qualitative from outline doc's zeotropic-blend/A2L
  maintenance notes).
  Results: CAPEX single-stage/cascade — R134a 1.00/1.82, R515B
  1.16/2.11*, R1234yf 1.34/2.44*. OPEX single-stage/cascade — R134a
  1.00/0.98, R515B 1.07/1.05*, R1234yf 1.21/1.19* (* = projected cascade).

- **`make_tewi_figure.py`** → `fig_tewi_breakdown.png` +
  `tewi_results.csv`. TEWI direct/indirect breakdown, 6 bars (3
  single-stage + 3 cascade, 2 of 6 hatched/projected). Went through two
  revisions on user request:
  1. Charge basis: started with an illustrative 5 kg charge, not made
     sufficiently prominent — user asked "we computed for unit charge?"
     then "unit charge please" → changed `CHARGE_KG` 5.0 → 1.0.
  2. Guideline source: started with AIRAH 2012 (Australian) for
     lifetime/leak-rate/EOL params, inconsistent with the deck's
     EPA/US regulatory framing — user asked for "an american guideline"
     → replaced with ASHRAE Journal LCCP (Nov 2018, worked example:
     15-yr life, 85% recovery/15% EOL loss) for the direct-term
     structure, and EPA GreenChill (2022 webinar,
     `gc-webinar-data-driven-leak-reduction-2022-04-12`) for the
     25%/yr US average leak rate (5%/yr shown as best-practice
     reference line).
  Final params: N_YEARS=15, L_ANNUAL=0.25, L_ANNUAL_BEST_PRACTICE=0.05,
  RECOVERY=0.85, CHARGE_KG=1.0, GRID_FACTOR=0.367 (EIA 2023 US avg).
  Hit an `IndexError` hatching bars 5/6 on a 6-bar container (should be
  4/5) — fixed.
  Final numbers (unit charge, 15yr, 25%/yr leak, 85% recovery):
  R134a 5.58t direct / 68.25t indirect / 73.83t total (7.6% direct);
  R515B 1.14 / 69.23 / 70.38 (1.6%); R1234yf 0.02 / 71.53 / 71.54
  (0.0%); Cascade R134a/CO2 5.58/67.14/72.72 (validated); Cascade
  R515B/CO2 1.14/68.11/69.25 (projected); Cascade R1234yf/CO2
  0.02/70.37/70.38 (projected).

## 2026-08-21 — Deck assembly

**`build_deck.py`** → `Cascade_Paper_Pitch.pptx`. Copied
`06_26_2026_Shilpa_group_meeting.pptx` as `template.pptx`, stripped all
slides (`_sldIdLst` manipulation + `drop_rel`) while keeping the
slide master/theme, rebuilt ~15 slides using `slide_layouts[0]`
(TITLE) and `[1]` (TITLE_ONLY) with custom `add_slide`/`add_bullets`/
`add_table`/`add_picture_slide` helpers. Verified visually via
`libreoffice --headless --convert-to pdf` + `pdf2image` render to PNG
before delivery.

**Known stale content, not yet fixed:** slide 12 ("TEWI Methodology")
and slide 13 caption still reference the OLD AIRAH-based parameters
(10-yr life, 15%/yr leak, 5kg charge) — offered to update, not yet
actioned as of this entry.

All figures + CSVs + the deck pushed to device at
`Cascade/pitch_deck/` via SendUserFile + device_commit_files.

## 2026-08-21 — Fact-check pass (user-built slides)

User built their own slides (in `08_21_2026_Shilpa_group_meeting.pptx`,
NOT `build_deck.py`'s output) combining my figures with their own bullet
text. Two rounds of fact-checking against current script/data:

- **Slide 16 ("Figure 9 Environmental Impact")**: bullets 1 and 3 were
  stale (pre-unit-charge, pre-American-guideline numbers); bullet 4
  actively contradicted the chart's own title ("90% recovery... 15%/yr
  leak" vs chart's "25%/yr leak, 15% EOL loss" = 85% recovery). Gave
  corrected bullet-4 text (25%/yr EPA GreenChill leak rate, 85% recovery
  / 15% EOL per ASHRAE Journal LCCP, 15-yr life). Bullets 1 and 3
  correction offered, not yet requested.
- **Slide 20 ("Environmental Computations")**: formula and parameter
  values all correct, but citation still said "AIRAH 2012 Best Practice
  Guideline" in two places — same stale-attribution issue. Gave the
  fix: swap to "ASHRAE Journal LCCP, Nov 2018" (formula line) and split
  attribution on the parameters line (ASHRAE Journal for n=15yr/
  recovery=85%, EPA GreenChill for L_annual=25%/yr).

Also gave: CAPEX/OPEX 4-bullet explainer + references for slide 15,
and two rounds of "4-point conclusions" — first generic, then remapped
to the paper's actual Results subsections (3.1 single-stage screening,
3.2 cascade performance, 3.3 techno-economic, 3.4 environmental impact)
after seeing the real deck's table of contents in
`08_21_2026_Shilpa_group_meeting.pptx`.

## Outstanding as of this entry

- Slide 12/13 in `build_deck.py`'s own output (`Cascade_Paper_Pitch.pptx`)
  still has stale AIRAH/5kg/10yr text — not yet fixed.
- Slide 20's citation fix (AIRAH → ASHRAE Journal/EPA GreenChill split)
  given as text, not yet applied to any file.
- Slide 16 bullets 1 and 3 correction offered, not yet applied.
- No CAPEX/OPEX slide exists in `build_deck.py`'s deck yet — only in the
  user's own `08_21_2026` deck (slide 15, "AI Generated Place holder").
- R515B/CO2 and R1234yf/CO2 cascades remain projected only — direct
  simulation is the next real modeling task, not a deck task.
