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
