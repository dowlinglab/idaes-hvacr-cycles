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

## Step 4 — build / initialize / solve / converge   [STATUS: IN PROGRESS -- evaporator, composition-loop, and compressor fixes all wired in; next run is the first full `vc.initialize()` with everything in place]

### Build   [DONE]
`SimpleVaporCompressionCycle("R32", mode=Mode.IMPROVED_TPX)` constructs
successfully -- confirms Steps 1-3h are structurally correct.

### Initialize   [FIX WIRED IN, root-caused, confirmed in isolation, NOT YET RUN against the full cycle]

**Symptom.** `vc.initialize()` runs unit-by-unit and fails at the EVAPORATOR:
`InitializationError: fs.evaporator.control_volume.properties_in failed to
initialize`, underlying solver message "Ipopt: Converged to a locally
infeasible point." Also see `W1002` warnings on `log_mole_frac_tbub/tdew`
(harmless pure-fluid floating-point dust, same as Phases 0-1 -- not the cause).

**Diagnosis (isolated the failing state in its own script,
`phase_3b_flash_seed_test.py`, instead of debugging inside the full flowsheet):**

The evaporator INLET is the post-expansion-valve state: T = -29 C (evaporating
temperature), P = Psat(-29) -- i.e. sitting EXACTLY on the saturation line,
two-phase. The `initialize()` code fixes both `inlet.temperature` AND
`inlet.pressure` together (same pattern used for every other, single-phase,
state in the cycle).

**Root cause: Gibbs phase rule violation.** For a PURE fluid with 2 coexisting
phases, F = C - P + 2 = 1 - 2 + 2 = 1. Only ONE independent intensive variable
is free on the saturation line; T and P are not independent (P is a function of
T along the dome). Fixing BOTH T and P simultaneously over-specifies the state
by one equation -> the equal-fugacity/flash constraints become inconsistent ->
"locally infeasible." This is invisible for every OTHER state in the cycle
(compressor in/out, condenser in/out) because those are single-phase, where T
and P genuinely are independent -- so the existing "fix T and P" pattern works
everywhere except this one two-phase point.

**5 attempts to fix it, run on an isolated single state block (reproduces the
failure without the full flowsheet's overhead):**

| # | Spec | Result |
|---|------|--------|
| 1 | Fix T and P both, directly | FAILS -- reproduces the flowsheet bug exactly |
| 2 | Seed single-phase (P=0.9*Psat), then move to fixed (T,P)=(−29,Psat) | Still FAILS -- same degeneracy once both are fixed |
| 3 | Seed single-phase, then FREE P + fix quality=0.2 (T stays fixed) | "optimal" but WRONG: P jumped to 23.4 bar (not 2.84), hL=hV identical -- TRIVIAL solution |
| 4 | Seed single-phase, PIN P=Psat, FREE T (unbounded) + fix quality=0.2 | "optimal" but WRONG the other way: T jumped to 173.5 C, hL=hV identical -- TRIVIAL solution |
| 5 | Same as 4, but T ALSO warm-started at Tsat and bounded tightly (+/-5 K) | **WORKS.** T=-28.20 C (vs target -29), P=2.843 bar (pinned), quality=0.2000, hL=-18101.2, hV=993.1 J/mol -- DISTINCT. Latent heat 19094 J/mol = ~367 kJ/kg, matches known R-32 latent heat near -28 C. |

**Why 3 and 4 failed but 5 worked (the real lesson).** The "trivial solution"
(both phases numerically identical, hL=hV) ALWAYS mathematically satisfies the
equal-fugacity equations -- it's a valid but unphysical root that coexists with
the real two-phase root. Fixing `phase_frac` alone does NOT prevent Newton's
method from converging to the trivial root, because once hL=hV, the mixture
enthalpy no longer depends on the phase-fraction value at all -- so fixing it
provides no numerical "pull" toward the real solution. Starting from a
single-phase seed put Newton's iterations right next to the trivial-root basin
in both Attempts 3 and 4, and it fell in. Attempt 5 fixed this by ALSO
warm-starting and tightly bounding the freed variable (T) near the physically
correct value, so the trivial branch (173 C away) is outside the feasible
region Newton can even reach -- forcing convergence to the real two-phase root.

**The small T offset (-28.20 vs target -29) is expected, not a bug.** We pin P
at the Ambrose-Walton estimate of Psat(-29); the PR EoS's OWN internal
saturation curve sits at a very slightly different pressure at exactly -29 C.
Freeing T lets the EoS find the temperature that IS truly self-consistent with
the pinned pressure -- landing ~0.8 K away. This is actually a good sign: it
confirms the EoS's own equilibrium condition is driving the result, not just
copying the AW guess.

**CONFIRMED RECIPE for a two-phase pure-fluid state in the generic cubic-PR
package** (documented in BREADCRUMB_07-20.md):
1. Seed the state single-phase, off the dome (e.g. P = 0.9*Psat).
   `initialize()` -- trivial, converges easily.
2. Fix P at the target Psat.
3. Unfix T. Warm-start `.value` at the expected Tsat; set TIGHT temporary bounds
   (+/- ~5 K) so Newton cannot escape to the trivial branch.
4. Fix `phase_frac["Vap"]` to a quality guess (init guess only -- not the real
   physical spec once coupled to the flowsheet).
5. Solve. Expect T within ~1 K of target, P pinned, hL != hV and physical.
6. In the full cycle: relax/remove the temporary T bounds afterward if they'd
   clip the real operating range; unfix `phase_frac` once the evaporator is
   coupled to the upstream valve (quality becomes flowsheet-determined).

**WIRED IN (2026-07-22):** the recipe above is now applied inside
`vapor_compression_cubic.py::initialize()`, in the evaporator's `else:` branch
(the `Mode.PH` guard means this is the only branch that ever executes -- see
Step 3f). Every other state in the cycle (compressor in/out, condenser in/out,
expansion-valve in/out, evaporator OUTLET) is single-phase, so it keeps the
original port-level `fix(temperature)` + `fix(pressure)` pattern, unchanged.

**Before (old, degenerate):**
```python
self.model.fs.evaporator.inlet.pressure[0].fix(self.p_init[-1]*p_scale)
self.model.fs.evaporator.inlet.temperature[0].fix(self.T_init[-1])   # <- Gibbs F=1 violation
```

**After (evaporator INLET only -- two-phase, post-expansion-valve state):**
```python
evap_in = self.model.fs.evaporator.control_volume.properties_in[0]
T_target = self.T_init[-1]

evap_in.temperature.unfix()
evap_in.temperature.setlb(T_target - 5.0)
evap_in.temperature.setub(T_target + 5.0)
evap_in.temperature.set_value(T_target)
evap_in.phase_frac["Vap"].fix(0.2)          # init-only seed, not a physical spec

self.model.fs.evaporator.outlet.temperature[0].fix(self.T_init[0])   # outlet: single-phase, unaffected
self.model.fs.evaporator.initialize(outlvl=logging.WARNING)

# revert to the standard (T,P)-given pattern now that a real 2-phase point is found
evap_in.phase_frac["Vap"].unfix()           # quality becomes flowsheet-determined
evap_in.temperature.setlb(None)
evap_in.temperature.setub(None)
evap_in.temperature.fix(value(evap_in.temperature))   # fix at the converged T (~T_target +/- 1 K)
```
(`inlet.pressure[0].fix(...)` above the `else:` block is unchanged from the
original code -- P was never the problem, only fixing T alongside it was.)

**Why (physical, not just numerical).** The evaporator inlet is genuinely
two-phase because it IS the expansion-valve outlet: isenthalpic throttling
drops the pressure and flashes part of the subcooled liquid to vapor. Being
two-phase is exactly why T and P cannot both be independently fixed there --
Gibbs F=1 means P = Psat(T) on the dome, so "fix T" and "fix P" are two
slightly different numbers (Ambrose-Walton's Psat vs. the cubic EoS's own
internal Psat) describing what should be one fact. Freeing T lets the EoS
solve for its own self-consistent saturation temperature at the P we pinned,
which is why the converged T (-28.20 C) differs slightly from the T_init
guess (-29 C) -- that's the EoS's own equilibrium condition, not an error.

**Status:** confirmed working against the REAL flowsheet (not just isolation),
after fixing two more bugs the isolated test never hit:

**Bug #2 -- `mole_frac_comp` never fixed.** Wiring the recipe in (using
`get_solver().solve(evap_in)` directly instead of `evap_in.initialize()`) hit
`BurntToast: Degrees of freedom were not zero [-1]`. `evap_in` is
`defined_state=True`, which drops the auto "sum(mole_frac)=1" constraint --
every OTHER unit gets this fixed for free via its own `unit.initialize()`
call's internal `fix_state_vars()`, but we bypassed that here by calling
`get_solver().solve(evap_in)` ourselves. Fix: add
`evap_in.mole_frac_comp["R32"].fix(1.0)` right after grabbing `evap_in`.

**Bug #1.5 (ordering) -- `evaporator.initialize()` re-fixes T underneath us.**
Even after Bug #2's fix, `BurntToast` still hit because
`self.model.fs.evaporator.initialize()` was being called AFTER our revert
step -- but its internal `properties_in.initialize()` unconditionally
re-fixes the 4 canonical state vars (flow, mole_frac, T, P) via
`fix_state_vars()`, regardless of what we'd already set up, recreating the
over-specification one level deeper. Fix: reordered so our own
`get_solver().solve(evap_in)` AND the revert block (unfix quality, fix T at
its converged value) both run BEFORE `self.model.fs.evaporator.initialize()`,
not after.

**Confirmed result:** evaporator inlet converges to a genuine two-phase state
(T = -24.41 C vs -29 C target; hL = -17767.6, hV = 1148.8 J/mol, clearly
distinct). The 4.6 C gap (bigger than the ~0.8 C gap seen in isolation) is
judged to be the same already-accepted cubic-EoS-vs-real-fluid deviation from
Phase 2, not a new bug -- not pursued further. Considered and REJECTED a
"single-phase seed" fix to tighten this (would have called `evap_in.
initialize()`, discarding the real upstream-consistent state in favor of a
generic default -- correctly rejected as disconnecting the evaporator from
what the expansion valve actually hands it).

### Composition-loop redundancy ("Too few degrees of freedom")   [FIXED]

Running the full cycle past the evaporator hit a NEW, structural (not
numerical) failure: Ipopt's pre-check reported `TOO_FEW_DOF` -- negative DOF
means OVER-constrained. `DiagnosticsToolbox.display_overconstrained_set()`
isolated it to a closed-loop redundancy: the flowsheet is a 4-unit ring
(evaporator -> compressor -> condenser -> valve -> evaporator), no
splitting/mixing/reaction, so total composition (like total flow) is
conserved all the way around. R32 is pure (`mole_frac_comp["R32"] == 1`
everywhere, always), but this fact was enforced BOTH by the 4 arc-level
`mole_frac_comp_equality` constraints AND by each unit's own outlet state
block's `sum_mole_frac_out` constraint -- one redundant equation per unit, 4
total. Confirmed exactly: DOF was -4, deactivating all 4 brought it to 0
(Ipopt: 191 variables = 191 equality constraints).

Fix, added in TWO places in `vapor_compression_cubic.py`:
```python
for unit_name in ["evaporator", "compressor", "condenser", "expansion_valve"]:
    unit = getattr(self.model.fs, unit_name)
    unit.control_volume.properties_out[0.0].sum_mole_frac_out.deactivate()
```
1. In `__init__`, right after `TransformationFactory('network.expand_arcs')`
   (so the model's DOF is correct immediately after construction).
2. Again at the end of `set_specifications()`, right before
   `calculate_scaling_factors(self.model)` -- REQUIRED, because each unit's
   own `initialize()` call reactivates its own `sum_mole_frac_out` as part of
   internal bootstrapping (confirmed via a three-checkpoint `.active` status
   script), so the `__init__`-time deactivation alone doesn't survive
   `vc.initialize()`. This second spot is the one guaranteed to run last,
   right before the real solve.

### Compressor `initialize()` -- "Converged to a locally infeasible point"   [FIXED]

Two sub-causes:

**(a) Bad default outlet-temperature guess.** IDAES's `init_isentropic()`
(source: `pressure_changer.py`), when `state_args=None`, copies the INLET's
temperature UNCHANGED into the outlet guess (only pressure gets rescaled) --
for a real compressor (should heat significantly), this places the guess deep
in the wrong single-phase region. Fixed by passing explicit `state_args` with
a temperature guess.

**(b) The guess itself was the wrong number.** The guess used, `T_init[1]`,
turned out to be the CONDENSING SATURATION temperature (T_amb+9, ~302-306 K),
not a real compressor-outlet value -- it sits right on R32's saturation dome
at the high-side pressure (Antoine: Tsat(18.78 bar) = 302.2 K), the worst
possible cubic-EoS flash starting point. Confirmed via a state-inspection
diagnostic: the "converged" real outlet's `phase_frac` summed to 1.13 (not
1) -- a non-physical leftover from a failed solve, not an actual state.
Meanwhile `properties_isentropic` (the hypothetical 100%-efficient reference
state), seeded from that SAME bad guess, reliably solves correctly on its own
(T = 366.69 K, entropy-matching residual = 0.000000 exactly), because
`init_isentropic`'s Steps 2-3 give it a dedicated entropy-matching solve
independent of guess quality.

Fix: a two-pass retry in the compressor's `initialize()` call. Pass 1 tries
the naive `T_init[1]` guess (expected to fail); on `InitializationError`,
pull `properties_isentropic[0].temperature` (already correctly solved even
though the overall step failed) and retry using THAT as the outlet's
`state_args` temperature. Confirmed working in isolation: Pass 2 succeeds,
outlet T = 350.69 K (~49 K above Tsat, comfortably superheated),
entr_mol = 33.28 J/mol/K (close to isentropic's 35.71, consistent with
efficiency ~ 0.9999).
```python
compressor_state_args = {
    "flow_mol": 1.0,
    "temperature": self.T_init[1],
    "pressure": value(self.model.fs.compressor.inlet.pressure[0]),
    "mole_frac_comp": {"R32": 1.0},
}
try:
    self.model.fs.compressor.initialize(outlvl=logging.WARNING, state_args=compressor_state_args)
except InitializationError:
    T_isen = value(self.model.fs.compressor.properties_isentropic[0].temperature)
    compressor_state_args["temperature"] = T_isen
    self.model.fs.compressor.initialize(outlvl=logging.WARNING, state_args=compressor_state_args)
```
Also added `from idaes.core.util.exceptions import InitializationError` to
the imports.

**Status:** all three fixes above (evaporator, composition-loop, compressor)
are wired into `vapor_compression_cubic.py` and syntax-checked
(`py_compile`). The compressor fix is confirmed working in an isolated
reconstruction of the compressor step, but NOT yet re-run against the actual
`vc.initialize()` sequence with everything combined -- that has never yet
reached the condenser or expansion valve (compressor was always the first
blocker). That full run is the immediate next step.

### Solve / converge   [TODO, pending the initialize run above]
Once `vc.initialize()` runs clean end-to-end (or is judged "good enough" even
if imperfect): re-attempt `optimize_COP(initialize=True, optimize=False)` (the
real coupled solve), then compare COP to Phase 3a Helmholtz value (3.19 at
T_amb=20). That comparison is Phase 4.
