# Breadcrumb — 2026-07-20 end of session

Read this first before resuming. Full history: PHASE0-3_NOTES.md, PHASE3B_CHANGES.md.

## Where things stand

**Phases 0-2: DONE, gates passed.** Cubic PR package (NIST/GCGP/SPGP) solves at
DOF=0 (Phase 0), matches vanilla PR_EOS.py on Z to ~1e-8 (Phase 1), matches Linde
dome to cubic-EoS accuracy (Phase 2).

**Phase 3a: DONE.** Helmholtz ideal cycle (eta~1, SH=SC=0), R-32, Te=-29C.
COP vs T_amb (10/15/20/25 C) = 3.95/3.53/3.19/2.91, all below Carnot at a steady
~0.75-0.78 fraction. Validated against Suhengki 2026 (cold storage, COP~2.8) and
Daikin 2016 (AC, COP 5.01, REFPROP). Pushed to origin/Colon_group_collaboration.

**Phase 3b: IN PROGRESS (2026-07-22, evening) — evaporator two-phase fix,
the closed-loop composition-redundancy fix, AND the compressor's
`initialize()` failure are all now RESOLVED. Next: re-run `vc.initialize()`
end-to-end (should now reach the condenser/expansion valve, never yet
exercised) and then the full `optimize_COP()` coupled solve. See "Update
(2026-07-22, evening)" below for the two newest fixes.**

## Update (2026-07-22, later): two more findings past the original recipe

**Bug #2 found: mole_frac_comp was never fixed.** Wiring the Attempt-5 recipe
into `vapor_compression_cubic.py` (unfix T, bound/warm-start, fix
`phase_frac["Vap"]=0.2`, then `get_solver().solve(evap_in)` directly instead of
calling `evaporator.initialize()` on this state) hit a NEW error:
`BurntToast: Degrees of freedom were not zero [-1]`. Diagnosed via a targeted
script printing `degrees_of_freedom(evap_in)` and every unfixed Var inside it
-- `DOF(evap_in) = 1`, and `mole_frac_comp[R32]` was in the unfixed list. R32
is pure, but this state block is `defined_state=True`, which drops the
"sum(mole_frac)=1" constraint (it expects the constraint or an explicit fix
from the caller) -- and nobody ever fixes `mole_frac_comp` anywhere in this
file. It was never a problem elsewhere because every OTHER unit's
`unit.initialize()` call auto-fixes it via `fix_state_vars()` internally; we
bypassed that automatic step for the evaporator by solving `evap_in`
ourselves. Fix: `evap_in.mole_frac_comp["R32"].fix(1.0)` -- same line the
isolated test script already had (`phase_3b_flash_seed_test.py`), just
forgotten in translation. Confirmed: DOF -> 0, solve -> optimal.

**Bug #1.5 (ordering): `evaporator.initialize()` itself re-fixes T.** Even
after the DOF fix, the FIRST attempt still hit `BurntToast` because the
sequence was: free T -> fix quality -> call `self.model.fs.evaporator.
initialize()`. That call's internal `properties_in.initialize()` ALWAYS
force-fixes the 4 standard state vars (flow, mole_frac, T, P) via
`fix_state_vars()`, regardless of what we'd already set up -- so it re-fixed T
right on top of our already-fixed quality, recreating the same
over-specification one level deeper. Fix: reordered so our own direct solve
(`get_solver().solve(evap_in)`) AND the revert (unfix quality, fix T at its
converged value) both happen BEFORE `evaporator.initialize()` is called, not
after. By the time `evaporator.initialize()` runs, T is already fixed (its
internal re-fix is then a no-op) and quality is already free -- standard
pattern, no conflict.

**Result once both fixes applied:** evaporator no longer crashes. But the
converged inlet T = -24.41 C, not -29 C (target) -- a 4.6 C miss, bigger than
the ~0.8 C miss seen in isolation. Confirmed NOT a trivial-root collapse:
hL = -17767.6, hV = 1148.8 J/mol, diff = 18916.4 J/mol (~364 kJ/kg latent
heat, physically sensible, matches R-32). So it's a REAL two-phase point, just
not at the exact T we guessed.

**Rejected fix (important -- don't redo this):** proposed adding a
"single-phase seed" step before freeing T (fix P at 0.9x target, call
`evap_in.initialize()`, then move P back) to mirror the isolated test's exact
procedure and hopefully tighten the T match. User correctly rejected this:
`evap_in.initialize()` would forcibly overwrite/discard whatever's already in
`evap_in` (including any real upstream-consistent guess) with a generic
single-phase default -- disconnecting the state from the rest of the cycle for
no good reason. NOT implemented. Left as-is.

**Reframe (the real answer on the 4.6 C offset):** -29 C was never a target
the cubic EoS is obligated to hit exactly -- it's OUR spec, converted to a
pressure via CoolProp (a high-accuracy real-fluid correlation), and then the
cubic PR EoS solves for ITS OWN self-consistent T at that pressure. Phase 2
already established the cubic model's saturation curve doesn't exactly match
the real fluid's ("matches Linde dome to cubic-EoS accuracy," not exact). The
4.6 C gap is very plausibly just that already-known, already-accepted model
deviation showing up again here, not a new bug. Not yet cross-checked
numerically against Phase 2's actual dome-comparison error near -28/-29 C --
worth doing if this becomes relevant again, but not blocking.

## Update (2026-07-22, evening): closed-loop DOF fix + compressor fix, both resolved

**Bug #3: closed-loop composition redundancy ("Too few degrees of freedom",
DOF=-4).** Ran the full-cycle try/except test (see "New plan" below) and hit a
NEW structural error, not a numerical one: Ipopt's own pre-check reported
`TOO_FEW_DOF` (negative DOF = OVER-constrained, not under -- corrected myself
on this mid-session). Diagnosed with `DiagnosticsToolbox.
display_overconstrained_set()`: the flowsheet is a closed 4-unit ring
(evaporator -> compressor -> condenser -> valve -> evaporator) with no
splitting/mixing/reaction, so composition is conserved all the way around,
just like flow. R32 is pure, so `mole_frac_comp["R32"] == 1` everywhere,
always -- but this trivial fact was being enforced BOTH by the 4 arc-level
`mole_frac_comp_equality` constraints AND by each unit's own OUTLET state
block's local `sum_mole_frac_out` constraint. That's one redundant equation
per unit, 4 total -- confirmed exactly: DOF was -4 before, deactivating all 4
`sum_mole_frac_out` constraints brought it to exactly 0 (Ipopt header:
191 variables = 191 equality constraints). This is the same *class* of bug as
the earlier single-arc `flow_mol_equality` deactivation, but shows up locally
at every unit (4 places) rather than needing just one arc broken -- and it's a
general closed-loop artifact, not pure-fluid-specific (would apply to mixtures
too).

Fix (now in `vapor_compression_cubic.py`, both in `__init__` right after arc
expansion, AND re-asserted at the end of `set_specifications()` right before
`calculate_scaling_factors`):
```python
for unit_name in ["evaporator", "compressor", "condenser", "expansion_valve"]:
    unit = getattr(self.model.fs, unit_name)
    unit.control_volume.properties_out[0.0].sum_mole_frac_out.deactivate()
```
Needed in BOTH places: each unit's own `initialize()` call reactivates its own
`sum_mole_frac_out` as part of internal bootstrapping (confirmed via a
three-checkpoint `.active` status script), so the `__init__`-time deactivation
alone doesn't survive `vc.initialize()`. Re-deactivating at the end of
`set_specifications()` -- the one place guaranteed to run last, right before
the real solve -- is what actually sticks.

**Bug #4: compressor `initialize()` -- "Converged to a locally infeasible
point."** Two sub-causes, both fixed:

1. *(Already partially fixed earlier today)* IDAES's `init_isentropic()`
   defaults (`state_args=None`) copy the INLET's temperature UNCHANGED into
   the outlet guess, only rescaling pressure -- placing the guess deep in the
   wrong single-phase branch. Fixed by passing explicit `state_args`.
2. *(New, root cause of the remaining failure)* The temperature we were
   passing, `T_init[1]`, turned out to be the CONDENSING SATURATION
   temperature (T_amb+9 ~ 302-306 K), not a real compressor-outlet guess -- it
   sits right on R32's saturation dome at the high-side pressure (Antoine:
   Tsat(18.78 bar) = 302.2 K), the worst possible starting point for a cubic
   EoS flash. Confirmed via a state-inspection diagnostic: the real outlet's
   `phase_frac` summed to 1.13 (not 1) -- a non-physical leftover from a
   failed solve, not a converged state. Meanwhile `properties_isentropic`
   (the hypothetical 100%-efficient reference state), seeded from the SAME
   bad guess, reliably solves to the correct answer on its own (T = 366.69 K,
   entropy-matching residual = 0.000000 exactly) because it gets its own
   dedicated entropy-matching solve inside `init_isentropic`'s Steps 2-3,
   independent of how bad the initial guess is.

   Fix: a two-pass retry, now in `vapor_compression_cubic.py`'s compressor
   init block. Pass 1 tries the naive `T_init[1]` guess (expected to fail);
   on `InitializationError`, pull `properties_isentropic[0].temperature`
   (already correctly solved even though the overall step failed) and retry
   with THAT as the outlet's `state_args` temperature. Confirmed working:
   Pass 2 succeeds, outlet T = 350.69 K (comfortably superheated, ~49 K above
   Tsat), entr_mol = 33.28 J/mol/K (close to isentropic's 35.71, consistent
   with efficiency ~ 0.9999).
   ```python
   try:
       self.model.fs.compressor.initialize(outlvl=logging.WARNING, state_args=compressor_state_args)
   except InitializationError:
       T_isen = value(self.model.fs.compressor.properties_isentropic[0].temperature)
       compressor_state_args["temperature"] = T_isen
       self.model.fs.compressor.initialize(outlvl=logging.WARNING, state_args=compressor_state_args)
   ```
   Needed `from idaes.core.util.exceptions import InitializationError` added
   to the imports.

**Not yet done:** re-run `vc.initialize()` end-to-end with both fixes in
place -- it has never yet reached the condenser or expansion valve (compressor
was always the first blocker). Then re-attempt the full `optimize_COP()`
coupled solve (testing H3: that a better sequential warm-start resolves the
earlier local-infeasibility/Restoration-Failed numerical issues there too).

## New plan: stop perfecting `initialize()`, try the real coupled solve instead

`vc.initialize()` (sequential, unit-by-unit) is now failing one step further
downstream, at the COMPRESSOR: `fs.compressor.control_volume.properties_out`
-> "Ipopt: Converged to a locally infeasible point." Compressor in/out are
single-phase -- NOT the Gibbs-rule issue, cause not yet diagnosed (possibly an
enth_mol bound issue, analogous to Phase 3a's Helmholtz enth_mass bound fix,
never yet checked/relaxed for the molar cubic package).

Key realization (from re-reading `optimize_COP`, lines 855+): `vc.initialize()`
is NOT the thing that actually needs to succeed for a real answer.
`optimize_COP(initialize=True, optimize=False)` does something different --
`solver.solve(self.model, ...)` on the WHOLE model at once (every unit, every
constraint, coupled), with the COP objective deactivated. That's the actual
full-system solve. AND `set_specifications()` unfixes everything `initialize()`
touched anyway (lines 536-551) before applying the real cycle spec. So the
compressor failure inside the fragile sequential `initialize()` may not
actually block reaching a real answer -- it's just a warm-start step, and a
partial/imperfect warm-start might still be enough for the full coupled solve
to converge.

**Test in progress (not yet run/confirmed):** wrap `vc.initialize()` in
try/except so we push forward regardless of the compressor error, then call
`set_specifications(ambient_temperature=20, condenser_approach=9,
evap_sat_temperature=-29, superheating=0, subcooling=0, max_pressure_ratio=10)`
and `optimize_COP(verbose=True, initialize=True, optimize=False)` directly, to
see if the full coupled solve succeeds despite the imperfect warm-start.
Awaiting output from the user's terminal.

## The core problem (fully diagnosed today)

Swapped Helmholtz -> generic cubic PR package inside `vapor_compression_cubic.py`
(all mechanical/structural edits DONE, logged line-by-line in
PHASE3B_CHANGES.md: property block, molar renames, phase_frac port->CV fix,
temperature_bubble, dropped eq_complementarity/eq_sat/diagrams). Builds fine.
Fails at the EVAPORATOR INLET during init: "Ipopt: locally infeasible."

**Root cause (confirmed via isolated single-state-block tests,
`phase_3b_flash_seed_test.py`):** the evaporator inlet is a two-phase point
(T=-29C, P=Psat). For a PURE fluid at 2 phases, Gibbs F=1 -- T and P are NOT
independent. The flowsheet init fixes BOTH T and P simultaneously -> degenerate
-> infeasible. This is NOT a scaling issue (Helmholtz had the same missing
scaling-factor warnings and still worked fine -- scaling is a red herring).

**Attempts on the isolated test block (all in `phase_3b_flash_seed_test.py`):**
- Attempt 1: fix (T,P) both on the dome -> FAILS (confirms the bug).
- Attempt 2: seed single-phase (P below Psat), then move to P=Psat with T still
  fixed -> still FAILS (same degeneracy once both are fixed again).
- Attempt 3: seed single-phase, FREE P + fix quality=0.2 (T fixed) -> "optimal"
  but WRONG: P jumped to 23.4 bar, hL=hV (identical) -> TRIVIAL ROOT collapse.
- Attempt 4: seed single-phase, PIN P=Psat, FREE T (unbounded) + fix quality=0.2
  -> "optimal" but WRONG again: T jumped to 173.5 C, hL=hV -> trivial collapse
  in the other direction.
- **Attempt 5 (CONFIRMED WORKING, 2026-07-22):** same as Attempt 4, but the
  freed T is ALSO (a) warm-started at the expected Tsat value and (b) bounded
  tightly (+/- 5 K) around it, so Newton cannot wander to the trivial branch.
  Result: solve optimal, T = -28.20 C (target -29, small AW-vs-PR-internal-Psat
  offset, expected/fine), P = 2.843 bar (pinned), phase_frac[Vap] = 0.2000,
  **hL = -18101.2, hV = 993.1 J/mol -- DISTINCT.** Latent heat = 19094 J/mol =
  ~367 kJ/kg, matches known R-32 latent heat near -28 C. REAL two-phase state.

**Root lesson (why 3 & 4 failed but 5 worked):** the trivial solution
(hL=hV, both "phases" identical) ALWAYS satisfies equal-fugacity, so simply
fixing quality does not prevent Newton from wandering to it from a single-phase
seed -- once hL=hV, phase_frac's fixed value doesn't affect any other residual,
so it provides no pull toward the real two-phase branch. The fix is to ALSO
anchor the freed variable (warm-start + tight bound) near the true saturation
value so Newton's basin of attraction is the real solution, not the trivial one.

## CONFIRMED RECIPE for a two-phase state in the generic cubic-PR package

1. Seed the state single-phase, off the dome (e.g. P = 0.9*Psat -> superheated
   vapor). Call `initialize()` -- trivial, converges easily.
2. Fix P at the target Psat (from Ambrose-Walton or wherever the cycle gets it).
3. Unfix T. Warm-start its `.value` at the expected Tsat, and set TIGHT temporary
   bounds around it (+/- ~5 K) so it cannot escape to a trivial single-phase
   region far from the dome.
4. Fix `phase_frac["Vap"]` to a quality guess (e.g. 0.2 for the evaporator inlet
   -- this is just an INIT guess; the real quality is set later by energy
   balance once the valve/evaporator are coupled).
5. Solve. Expect T to land within ~1 K of the AW/target Tsat (small AW-vs-
   internal-PR-Psat mismatch is normal), P pinned, quality as fixed, hL != hV
   and physically sensible.
6. (For the full cycle) after init, relax/remove the temporary T bounds if they
   would clip the real operating range, and unfix `phase_frac` once the
   evaporator is coupled to the upstream valve (quality becomes flowsheet-
   determined, not an independent spec) -- check DOF at each step.

## Next actions (in order, superseded/updated 2026-07-22 evening)

1. **[DONE]** Evaporator two-phase recipe wired in and confirmed against the
   real flowsheet (mole_frac_comp DOF bug + fix_state_vars reordering bug,
   both fixed).
2. **[DONE]** Closed-loop composition-redundancy fix (Bug #3): 4x
   `sum_mole_frac_out` deactivated in both `__init__` and
   `set_specifications()`. Confirmed DOF -> 0 (191=191).
3. **[DONE]** Compressor `initialize()` fix (Bug #4): explicit `state_args` +
   two-pass retry using the isentropic block's own solved temperature.
   Confirmed working in isolation (Pass 2 succeeds, T=350.69K). Just wired
   into `vapor_compression_cubic.py` -- NOT yet re-tested against the full
   `vc.initialize()` sequence.
4. **[NEXT]** Re-run `vc.initialize()` end-to-end. Should now get past the
   compressor for the first time and reach the condenser + expansion valve
   (never yet exercised -- may surface their own analogous issues).
5. Then re-attempt `optimize_COP(initialize=True, optimize=False)` (the real
   coupled solve) -- testing whether a fully-successful sequential warm start
   resolves the earlier local-infeasibility / Restoration-Failed numerical
   issues seen there (hypothesis H3, not yet confirmed).
6. If the full solve succeeds: proceed straight to comparing COP vs Phase 3a
   Helmholtz (3.19 at T_amb=20). That comparison is Phase 4 (trust gate).
7. Once a COP number exists: update PHASE3B_CHANGES.md Step 4 status from
   "wired in, awaiting run" to actual results, then commit/push (Task #21).

## Key conceptual points to remember (asked about this session, now resolved)

- "Generic flash" = the VLE/phase-equilibrium solver component INSIDE our
  generic cubic PR package (SmoothVLE + LogBubbleDew + log_fugacity) — not a
  separate package. We are and have always been using the generic cubic EoS;
  the flash is one sub-piece of it that struggles for a pure fluid on the dome.
- Why it struggles: built for MIXTURES (composition differences anchor the
  solve); a pure fluid has no composition difference between phases, bubble=dew
  degenerate, F=1 means T,P aren't independent -> singular/trivial-root prone
  exactly on the saturation line. Helmholtz avoids this because it has a
  dedicated pure-fluid flash (quality as a direct state var, direct Psat).
- Phase 2's dome computation NEVER exercised this: it evaluated the two EDGES
  (saturated liquid root, saturated vapor root) as separate SINGLE-PHASE
  packages — never asked for coexistence in one state block. The evaporator
  inlet is a genuine INTERIOR two-phase mixture — coexistence is unavoidable
  there, unlike in Phase 2.
- Ambrose-Walton's role here: NOT a flash replacement. It supplies the Psat
  guess/pin used in the seeding recipe (Attempts 2-4). The actual robustness
  comes from the seeding STRATEGY (single-phase seed -> well-posed two-phase
  spec), not from AW itself.
- Considered and set aside: a custom pure-fluid cubic-PR property package
  (mirrors Helmholtz's quality-based approach exactly) — Route B, a real build,
  fallback only if seeding ultimately fails. Also considered and REJECTED:
  computing 4 state points externally and skipping the flowsheet — rejected
  because the whole point of using IDAES is solving the coupled system, not
  hand-computed state points (user was explicit about this).
- Molar-basis strategy for de-risking (Helmholtz MASS->MOLE first, prove COP
  unchanged, THEN swap EoS) is pinned as task #22 but NOT used — we did the
  rename and EoS swap together instead, per user's later call ("No. Put a pin
  on it for later" -> then proceeded directly). Still on the task list if a
  future basis issue comes up.

## File inventory (R32/)

- `phase_0_cubic_roots_test.py`, `phase_1_cubic_eos_validation.py` (has
  `make_config`, `METHODS` — imported everywhere), `phase_2_saturation_dome.py`
  (has `ambrose_walton_psat` — imported by the flash test) — all DONE, gates
  passed, pushed.
- `vapor_compression_plr.py` — Phase 3a Helmholtz cycle (working, unchanged,
  pushed). Use this as the Phase-4 comparison baseline.
- `vapor_compression_cubic.py` — Phase 3b cubic-PR cycle (structural edits done,
  init blocked at evaporator inlet as above). Pushed (WIP commit).
- `phase_3a_helmholtz_cop.py`, `phase_3a_validation_literature.py` — Phase 3a
  drivers, done, pushed.
- `phase_3b_flash_seed_test.py` — isolated flash test, Attempts 1-4 written;
  **Attempt 4 needs to be run and confirmed next session.**
- `PHASE3_NOTES.md`, `PHASE3B_CHANGES.md` — living docs, keep appending.
- `cop_vs_ambient_r32.csv/.png` — Phase 3a results, pushed.

## Task list state (as of this session)

#19 Phase 3 (cycle->COP) still in_progress — Phase 3b is the open sub-item.
#20 Phase 4 (COP vs Helmholtz trust gate) pending on 3b.
#21 verify/document/push each phase — ongoing, keep doing per-phase.
#22 molar-basis strategy — pinned/deferred, not currently the active path.

## Immediate next command to run

```
cd ~/Desktop/DVCT_Project/property_packages_DVRT_code/idaes-hvacr-cycles/R32
python phase_3b_flash_seed_test.py
```
Read the Attempt 4 block output. If T~-29C, P~2.843 bar, phase_frac[Vap]=0.2000,
and hL != hV with sensible magnitudes -> the fix is confirmed, go implement it
in `vapor_compression_cubic.py`'s evaporator-inlet initialization.
