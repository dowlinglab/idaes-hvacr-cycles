# Note — Getting saturated properties out of IDAES (Gibbs phase rule)

**Date:** 2026-07-20 · Shilpa Narasimhan

At saturation a pure fluid sits on a single P = Psat(T) curve where liquid and
vapor coexist. You cannot just fix (T, P) and expect coexistence -- if P != Psat
you land in a single phase (20 C / 10 bar, for example, is superheated vapor).

## Why fixing T alone pins the state

Gibbs phase rule: **F = C - P + 2**. A pure fluid has C = 1 component; two
coexisting phases means P = 2, so **F = 1 - 2 + 2 = 1**. On the saturation curve
there is exactly ONE free variable. Choose T, and Psat, the saturated liquid,
and the saturated vapor are all determined.

## How to impose it in IDAES (FTPx DOF swap)

FTPx state variables are flow, T, P, x. Normally all four are fixed -> 0 DOF ->
a single flash at the named (T, P). To sit on the saturation curve instead:

1. Fix flow, x, T; fix P at a guess; initialize.
2. Free the pressure (+1 DOF).
3. Fix the vapor fraction `phase_frac["Vap"] = 0.5` (-1 DOF) -- the equation
   "two phases coexist."
4. Net DOF still 0. The solver finds the pressure at which two phases coexist at
   that T (= Psat) and returns both phase states in one solve:
   `enth_mol_phase["Liq"]` = hf, `enth_mol_phase["Vap"]` = hg (same for entropy).

## Why the 0.5 is arbitrary

For a pure fluid the per-phase properties at coexistence depend only on T (from
the phase rule). The vapor fraction sets only how much of each phase exists by
amount, not what the saturated liquid/vapor is or the pressure. So 0.5, 0.3, 0.9
all give the identical Psat, hf, hg. Pick 0.5 because a solidly two-phase target
is easiest for the solver. (A mixture would differ -- composition shifts with the
fraction -- but R-32 is pure.)

---

# Robust extraction method (final) and results

## Why the flash approach was abandoned

Fixing (T, Psat) and running the full VLE `initialize()` (bubble/dew + equal
fugacity) is numerically fragile: the flash sits exactly on the saturation line
where liquid and vapor are on the verge of swapping, and ipopt intermittently
"converges to a locally infeasible point." It failed unpredictably point to
point across NIST/GCGP/SPGP.

## The robust method

We do not need equilibrium -- only each phase's cubic ROOT at (T, Psat): the
liquid root gives hf/sf, the vapor root gives hg/sg. So we build TWO single-phase
packages (liquid-only, vapor-only) with the cubic EoS and NO VLE machinery
(`make_phase_config`, `eval_phase`). A single-phase block at fixed (T, P) has
nothing to flash -- it just evaluates that phase's root. Psat per method from
Ambrose-Walton (`ambrose_walton_psat`, closed form from Tc/Pc/omega). Every point
solved "Optimal Solution Found" for all three methods, full -50..74 C window.

## Results (IIR-anchored, vs Linde)

| series | P_MAPE% | hf_MAE [kJ/kg] | hg_MAE | sf_MAE [kJ/kg.K] | sg_MAE |
|--------|---------|----------------|--------|------------------|--------|
| NIST   | 1.53    | 12.31          | 17.70  | 0.0385           | 0.0589 |
| GCGP   | 11.91   | 21.62          | 22.13  | 0.0754           | 0.0758 |
| SPGP   | 43.66   | 80.57          | 119.33 | 0.2792           | 0.4399 |

- **NIST**: ~12-18 kJ/kg on h -- cubic-EoS accuracy, matches the standalone
  `compare_cp_methods.py` (14 / 15.7). Gate PASSED.
- Ranking **NIST << GCGP < SPGP** reproduced, consistent with the standalone
  comparison. The IDAES package reproduces both the dome and the relative quality
  of the three parameter sets.

## Gate

- **Gate: NIST dome matches Linde to cubic-EoS accuracy (~15 kJ/kg on h).**
  ✅ Passed (12.3 / 17.7). GCGP moderate, SPGP poor -- as expected.

---

# Debugging log — how the "trivially infeasible pressure" error was diagnosed

**Symptom.** During the full dome sweep, a long Pyomo/nl_writer traceback ending
in:
`InfeasibleConstraintException: model contains a trivially infeasible variable
'fs.state[0].pressure' (fixed value 153.15... outside bounds [10000.0, 1e7]).`

**How the root cause was found (read the traceback bottom-up):**

1. **The last line is the signal; the stack above it is noise.** The 30 lines of
   Pyomo internals (nl_writer, visitor, ampl) just say "while writing the model
   for the solver." The final exception line names the exact fault: pressure was
   *fixed* to 153.15 Pa, which is outside its bounds [1e4, 1e7].

2. **Trace the fixed value to its source.** Pressure is fixed by
   `sat_point(p, T)` to `antoine_psat(T)`. The sweep loop starts at the first
   Linde row, -130 C = 143.15 K. Evaluating `antoine_psat(143.15)` gives
   10^(4.60123 - 959.89766/(143.15-13.71589)) bar = 1.53e-3 bar = 153 Pa --
   exactly the offending value.

3. **Cross-check against physical data.** Linde Psat at -130 C = 0.001312 bar =
   131 Pa. So 153 Pa is physically correct -- the *bounds* are wrong, not the
   calculation. (Same class of problem for temperature: 143 K < the 200 K lower
   bound.)

4. **Root cause.** `state_bounds` in `make_config` were set for the Phase 0/1
   single test state (20 C, 10 bar) -- pressure floor 1e4 Pa, temperature floor
   200 K. Phase 2 sweeps the *full* Linde range down to -130 C / 131 Pa, far
   below those floors.

**Fix.** Widen the bounds for the sweep: temperature floor 200 -> 140 K,
pressure floor 1e4 -> 1.0 Pa. Upper bounds already cover 78 C / 57.7 bar.

**Takeaways.**
- Read tracebacks bottom-up: the final exception line is the fact; the stack is
  the path.
- "Trivially infeasible ... outside bounds" = a fixed variable violates its
  declared bounds -> check what value it was fixed to and whether the bounds fit
  the range you're now running over.
- Bounds chosen for one operating point silently break when the same package is
  reused over a wider sweep.
