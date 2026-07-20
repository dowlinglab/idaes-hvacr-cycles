# Phase 3b — Change Log: Helmholtz -> cubic PR in the VC cycle

**File edited:** `vapor_compression_cubic.py` (a copy of `vapor_compression_plr.py`,
the working Phase-3a Helmholtz cycle).
**Goal:** swap the Helmholtz EoS for the cubic PR generic property package
(Phases 0-2) so COP is computed *inside* the IDAES flowsheet with our model.
**Author:** Shilpa Narasimhan · Support: Claude AI · **Date:** 2026-07-20

Line numbers are AFTER the Step-2 edits (they shift as edits are made; treat as
approximate anchors).

---

## Step 1 — property package (`__init__`)   [STATUS: DONE]

- Added imports (~lines 27-30):
  `from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock`
  `import os, sys; sys.path.append(os.path.dirname(os.path.abspath(__file__)))`
  `from phase_1_cubic_eos_validation import make_config, METHODS`
- Line 66: `self.model.fs.properties = GenericParameterBlock(**make_config(METHODS["NIST"]))`
  (replaced the `HelmholtzParameterBlock(...)` block).
- Removed the `if mode == Mode.PH: sv = StateVars.PH else StateVars.TPX` block
  (generic package is FTPx; no PH/TPX choice). `state_vars=sv` / `AmountBasis`
  no longer used.
- Old `from idaes.models.properties.general_helmholtz import (...)` left in place
  (unused, harmless).
- Reason: FTPx is molar; there is no mass-basis / PH state option.
- Verified: syntax OK; no `self.model.fs.properties = HelmholtzParameterBlock`,
  no `sv = StateVars` remaining.

---

## Step 2 — mechanical basis/phase renames (Find & Replace)   [STATUS: DONE]

Rationale: Helmholtz ran in MASS basis; the generic PR package runs in MOLAR
basis (FTPx). Same physical quantities, different variable names/units. `vapor_frac`
(quality, time-indexed) becomes `phase_frac["Vap"]` (phase-name-indexed) in the
modular framework.

| Find | Replace | # | Lines (approx, post-edit) |
|------|---------|---|---------------------------|
| `flow_mass` | `flow_mol` | 5 | 115, 382, 485, 486, 548 |
| `enth_mass` | `enth_mol` | 8 | 388, 389, 406, 428, 429, 489, 490, 900 |
| `entr_mass` | `entr_mol` | 1 | 903 |
| `vapor_frac[0]` | `phase_frac["Vap"]` | 14 | 411, 434, 493, 495, 505, 506, 578, 580, 632, 634, 688, 690, 795, 796 |

Notes:
- One `vapor_frac` intentionally NOT changed: line 478, inside a comment
  ("...bound_vapor_frac..."). Left as-is.
- Units change with the rename: flow kg/s -> mol/s; enthalpy J/kg -> J/mol;
  entropy J/kg.K -> J/mol.K. Guess values and bounds fed in later steps must be
  molar accordingly.
- `flow_mol_equality` (line 115) is the arc-expansion constraint name after the
  rename (was `flow_mass_equality`); it is deactivated to leave the loop open
  (closed-cycle flow set once).

---

## Step 3 — hand-edit remaining Helmholtz-only refs   [STATUS: DONE]

### 3b — eq_complementarity removed   [DONE]
Deleted the 4 `.eq_complementarity.deactivate()` calls (were lines ~579, 632,
688, 795) for evaporator / compressor / condenser / expansion_valve; replaced
with `# Phase 3b:` comments. For the evaporator the guarding
`if self.mode == Mode.IMPROVED_TPX:` was removed with it (it guarded only that
call). For the expansion valve the `elif Mode.IMPROVED_TPX:` block became empty,
so replaced its body with `pass`. Reason: generic SmoothVLE has no
complementarity variable.

### 3c — eq_sat removed   [DONE]
Deleted the `.eq_sat.activate()` call (was line ~798), inside the same expansion
valve `elif` now bodied with `pass`.

### 3d — diagrams removed   [DONE]
Commented out all Helmholtz-only diagram calls + their `plt.plot/plt.show/
add_warning` overlays in three methods:
- `draw_thermodynamic_diagrams` (was ~266-274) -> body is `pass`.
- `specify_initial_conditions` (was ~351-367) -> removed the hp/pt/ts overlays
  and the now-unused `self.S_init` loop; h_init/p_init/T_init guesses kept.
- `report_solution` (was ~911-924) -> removed the three diagram overlays.
Verified: 0 actual `eq_complementarity`/`eq_sat` calls, 0 `*_diagram()` calls,
syntax OK.

### 3a — temperature_sat -> temperature_bubble["Vap","Liq"]   [DONE]
Find & Replace all `temperature_sat` -> `temperature_bubble["Vap", "Liq"]`
(6 refs: constraint lines 140/146/156/162/168 and the comment ~235). Why: the
generic PR package has no `temperature_sat`; for a pure fluid the Vap-Liq bubble
temperature equals the saturation temperature, so the superheat/subcool/vapor
constraints reference `temperature_bubble["Vap","Liq"]` instead. Verified: 0
`temperature_sat` remain, 6 `temperature_bubble["Vap","Liq"]` present.

### 3e — CoolProp guesses to molar   [NO EDIT NEEDED]
Investigated: every `enth_mol.fix(self.h_init[...])` (lines 364, 365, 382, 404,
405) sits inside an `if self.mode == Mode.PH:` branch (362/381/403). We run
`Mode.IMPROVED_TPX`, so the `else` path runs and fixes TEMPERATURE (`T_init`, in
K) and pressure -- both basis-independent. `self.h_init` (the CoolProp mass-unit
H guesses) is therefore never used for fixing in our mode, so no unit conversion
is required. Bonus: this avoids a CoolProp-vs-PR reference-frame mismatch (their
absolute enthalpy datums differ) that would otherwise corrupt the guess. The
`PropsSI('P'/'T')` calls that ARE used return pressure/temperature, which are
datum-/basis-independent. (Line 317's `PropsSI('T','H',...)` is inside a
triple-quoted comment block -- already dead.)

### 3h — phase_frac: port -> control-volume property   [DONE]
After the `vapor_frac[0]` -> `phase_frac["Vap"]` rename (Step 2), all 12 active
uses were on PORTS (`unit.outlet.phase_frac["Vap"]`, `unit.inlet.phase_frac`).
Helmholtz exposed `vapor_frac` on the port, but the generic FTPx port carries
only flow_mol / mole_frac_comp / T / P -- `phase_frac` lives inside the state
block. Replaced (Find & Replace all):
  `.outlet.phase_frac["Vap"]` -> `.control_volume.properties_out[0].phase_frac["Vap"]`
  `.inlet.phase_frac["Vap"]`  -> `.control_volume.properties_in[0].phase_frac["Vap"]`
Now 0 port refs, 14 control-volume refs. Why: otherwise every `.fix/.setlb/.setub/
.unfix` on `port.phase_frac` would AttributeError at build.
CAVEAT (for Step 4): fixing `phase_frac` interacts with the flash degrees of
freedom (Phase-2 lesson -- fixing a phase fraction while other things are fixed
can over-/under-determine the flash). Expect to revisit these fixes when the
model is initialized/solved.

### 3f — guard against Mode.PH   [DONE]
Added in `__init__` (after `self.mode = mode`):
`assert mode != Mode.PH, "Mode.PH is not supported with the generic PR (FTPx) package"`.
Why: the generic package uses the FTPx (molar T-P-flow-x) state definition, which
has no pressure-enthalpy state. Rather than leave ~13 dead `if Mode.PH:` branches
that would silently mis-specify the model if PH were ever passed, we fail fast
with a clear message. (The dead branches are left in place, now unreachable.)

### 3g — remove unused Helmholtz imports   [DONE]
Deleted the `from idaes.models.properties.general_helmholtz import
(HelmholtzParameterBlock, PhaseType, StateVars, AmountBasis)` block (was lines
3-8), replaced with a `# Phase 3b` comment. Why: after Steps 1-3, none of those
four names are used anywhere in the file (each appeared only in the import). The
cycle now depends only on the generic cubic PR package (via make_config), so
importing general_helmholtz was dead weight and misleading about the file's EoS.
Verified: no `HelmholtzParameterBlock`/`StateVars`/`AmountBasis`/`PhaseType`
usages remain (only the explanatory comment mentions "Helmholtz").

---

## Step 4 — build / initialize / solve / converge   [STATUS: TODO]

Build the object, fix construction errors, then initialize (expect molar-unit and
enth bound issues -- `relax_enth_bounds` in molar, R-32 vapor ~26 kJ/mol), solve
the ideal cycle (eta ~ 1, SH=SC=0) at T_amb = 20, and compare COP to Phase 3a
(3.19). That comparison is Phase 4.
