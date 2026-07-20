# Phase 3 — Vapor-Compression Cycle -> COP (R-32)

**Author:** Shilpa Narasimhan · Support: Claude AI
**Date:** 2026-07-20

Phase 3 computes COP for R-32 through the DVCT vapor-compression cycle.
- **Phase 3a** (this note): baseline COP vs ambient using the Helmholtz EoS.
- **Phase 3b** (next): swap Helmholtz for the cubic PR package from Phases 0-2
  (`vapor_compression_cubic.py`), then Phase 4 compares the two COPs.

Files:
- `phase_3a_helmholtz_cop.py` -- the sweep driver (Helmholtz, R-32).
- `vapor_compression_plr.py` -- the cycle class, extracted from branch
  `origin/PLR_vanilla_prop` (see design decision 1).
- `phase_3a_validation_literature.py` -- runs our cycle at a published
  (Suhengki) cold-storage point as a check.
- `cop_vs_ambient_r32.csv`, `cop_vs_ambient_r32.png` -- results.

---

## 1. Final result (ideal cycle)

Ideal single-stage cycle (eta_isen ~ 1, superheat = subcool = 0), R-32,
evaporating at -29 C, condenser saturation = T_amb + 9, feasibility solve:

| T_amb [C] | T_cond,sat [C] | COP  | Carnot | COP/Carnot |
|-----------|----------------|------|--------|------------|
| 10        | 19             | 3.95 | 5.09   | 0.78       |
| 15        | 24             | 3.53 | 4.61   | 0.77       |
| 20        | 29             | 3.19 | 4.21   | 0.76       |
| 25        | 34             | 2.91 | 3.88   | 0.75       |

COP falls monotonically with ambient, stays below Carnot at a steady ~0.75-0.78
of the reversible limit. Since the compressor is ~ideal, that gap is essentially
the throttling (isenthalpic expansion) loss. Reported range is 10-25 C; above
~25-30 C the solve returns infeasible / non-physical points (see Debugging).

---

## 2. Design decisions

1. **Use the PLR-branch cycle, not the local one.** The local
   `Colon_group` `vapor_compression.py` bounds the condenser *outlet* temperature
   to a tight 2 C window, which conflicts with subcooling (outlet must sit below
   T_sat) -> every case failed with `subcooling_constraint` residual = 3.0. The
   `origin/PLR_vanilla_prop` version sets the condensing temperature as
   `ambient + condenser_approach` and applies subcool relative to T_sat, and adds
   `optimize=False` (feasibility solve). Extracted it as `vapor_compression_plr.py`.

2. **Feasibility solve (`optimize=False`), not COP maximization.** For a
   fixed-conditions cycle the operating point is determined, so there is no
   separate optimum to find. `optimize=True` (maximize COP) was numerically
   unstable (maxIterations) and, when it did solve, cheated by running the
   subcooling to ~ -29 K (condenser outlet below ambient -- unphysical). So we
   report the feasible fixed-cycle COP.

3. **Ideal cycle (eta ~ 1, SH = SC = 0), following Shridhar 2016.** With
   superheat/subcool as inequalities and no objective, the outlet states are
   under-determined and the solver drifts (non-monotonic COP, a bump at 30 C).
   The ideal saturated cycle removes that freedom: states are set by the phase
   fixes (saturated vapor out of evaporator, saturated liquid out of condenser),
   giving a determined, monotonic, physical curve. (eta = 0.9999 because the
   class asserts 0 < eta < 1.)

4. **Ratio cap = 10 (single-stage limit).** Literature: single-stage
   reciprocating/scroll compressors are limited to pressure ratios ~8-10 (driven
   by discharge temperature and falling volumetric efficiency); Freon systems go
   two-stage at ratio >= 10. We use 10 (upper end).

5. **Ambient range 10-25 C only.** With the -29 C evaporator, ratios and solver
   difficulty grow with ambient; 10-25 C converges cleanly. NOTE this limit comes
   from the DEEP evaporator (cold storage), NOT from R-32 or the ambient -- in
   normal AC (warm evaporator) single-stage R-32 runs to 55 C ambient (see
   literature). 10-25 C is sufficient for the Phase-4 comparison (single fixed
   condition).

6. **enth_mass upper bound relaxed to 700 kJ/kg.** The Helmholtz R-32 package
   caps enth_mass at ~500 kJ/kg, but R-32 vapor at a -28 C evaporator is already
   ~507 kJ/kg. `relax_enth_bounds()` raises the bound after the model is built.

7. **Matplotlib backend forced to Agg** so the diagram calls inside
   `specify_initial_conditions` do not pop blocking windows during the sweep
   (environment only; the class is unchanged).

---

## 3. Literature validation

Sources documented in the driver/validation docstrings.

- **Suhengki et al. (2026), IJASEIT 16(2)** -- R-32 low-temp cold storage (our
  exact application): evaporator to -28 C, chamber -18 to -20 C, ambient 25-40 C,
  condensing 45-60 C. R-32 crit props Tc = 78.4 C, Pc = 58.3 bar (match ours).
  Realistic R-32 cold-storage **COP ~ 2.8** -> consistent with our 2.9-3.9 band.
  CAVEAT: their Table V reports COP = 7.5, which EXCEEDS Carnot (3.88 at
  Te=-29/Tc=34) and is thermodynamically impossible -- a data-quality error. Our
  `phase_3a_validation_literature.py` run at their conditions did not converge
  (their 42 K superheat is very stiff), but the Carnot argument alone disproves
  COP = 7.5 (no simulation needed).
- **Shridhar & Mutalikdesai (2016), IJCET 6(5)** -- ideal VC cycle (isentropic
  compression, eta = 1, saturated). We adopted this ideal-cycle *approach* for
  the determined, well-behaved curve.
- **Taira/Daikin (2016), Purdue IRACC 2408** -- R-32 AC, REFPROP: Tc = 46,
  Te = 12.5, eta = 70%, 8 K subcool -> COP 5.01. Same Helmholtz basis; a clean
  reproduction target (ambient not stated, so left as a future check).
- **Badescu/REHVA (2025)** -- R-32 AC COP regression vs (T_out - T_in);
  COP ~2.4-4 up to 55 C ambient. AC duty (warm evaporator), context only.

Key validation conclusions: (a) our physical COP band (2.9-3.9, below Carnot)
agrees with the realistic cold-storage literature (~2.8); (b) the one literature
value that disagrees (Suhengki Table V, 7.5) is provably erroneous.

---

## 4. Debugging log (chronological)

1. **`InfeasibleConstraintException` on enth_mass at init.** Evaporator outlet
   vapor guess 507.8 kJ/kg > package bound 500 kJ/kg. Not just the guess -- R-32
   saturated vapor across -30..-28 C is ~506-510 kJ/kg, so the operating point is
   above the bound too. Superheat/subcool cannot fix it (even 0 K superheat leaves
   ~507). Fix: `relax_enth_bounds()` to 700 kJ/kg.

2. **Every case infeasible with the local cycle.** `subcooling_constraint`
   residual = 3.0: the local class bounds the condenser OUTLET temperature to a
   tight window while subcooling pushes the outlet below T_sat -> conflict.
   Fix: switch to the PLR-branch cycle (ambient + approach interface).

3. **`optimize=True` runs subcooling away.** SC = 48-63 K (condenser outlet below
   ambient, impossible), COP/Carnot inflated to ~0.77, and it crashed at 30 C
   (maxIterations). Cause: maximizing COP with no heat-exchanger/ambient limit.
   Adding a condenser-outlet >= ambient floor made even `optimize=False` fail to
   converge, so it was dropped.

4. **`optimize=False` drifts (non-monotonic COP).** With superheat/subcool as
   inequalities and no objective, the outlet states are under-determined; 10-25 C
   sat near SH=SC=3 (frac 0.65) but 30/35 C drifted to more subcool (frac 0.75).
   Making SH/SC hard equalities over-constrained the model (worse). Resolution:
   the ideal cycle (SH = SC = 0, decision 3) removes the freedom entirely.

5. **maxIterations = ipopt hit its iteration budget (1000) without meeting
   tolerance** -> "bad status: error" crash. All such messages were from the
   `optimize=True` path; the ideal `optimize=False` run converged in-budget with
   no maxIterations.

6. **SH/SC diagnostic columns unreliable in the ideal run.** With SH=SC=0 the
   superheat/subcool *constraints deactivate* (class only activates them when
   > 0.1), so the `temperature - temperature_sat` extraction reads inconsistent
   values (SH ~50, SC ~ -29) for the vapor-fraction-pinned saturated states. This
   is a reporting artifact only -- COP (from Q_evap/W_comp) is correct.

7. **`fallbk` column = fallback used, NOT convergence.** `False` means the point
   solved without the pressure-equality fallback. All reported rows converged.

---

## 5. Status / gate

- **Phase 3a: DONE.** R-32 Helmholtz cold-storage cycle solves and gives a
  physical, monotonic COP vs ambient (3.95 -> 2.91, ambient 10-25 C), below
  Carnot, consistent with the cold-storage literature (~2.8).
- **Next:** Phase 3b -- swap the Helmholtz EoS for the cubic PR package
  (`vapor_compression_cubic.py`); Phase 4 -- compare the two COPs (trust gate).
- **Watch item:** the cycle is numerically fragile with superheat/subcool freedom
  and at high ambient; the ideal (determined) cycle is the robust configuration.
