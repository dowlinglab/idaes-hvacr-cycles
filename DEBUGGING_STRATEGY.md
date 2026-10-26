# Debugging strategy — process models that converge and are still wrong

Written 2026-10-08, after a flowsheet defect survived roughly five months of
active use. Every run reported `converged = True`. Every trend was monotonic.
Every number looked plausible. The loop did not conserve energy.

This file is the procedure that would have caught it in a week instead of
five months. It is written to be read by a person or by Claude at the start
of a session.

---

## The core rule

**`converged = True` is not validation.** It means the solver found a point
satisfying the constraints it was handed. If a constraint is *missing*, a
physically impossible answer converges perfectly and reports success.

A model can be:

- **solvable** — the solver terminates
- **self-consistent** — the equations present are satisfied
- **wrong** — the equations present are not all the equations

All three at once is the dangerous case, because nothing in the output
distinguishes it from a correct model.

---

## Check against physics, not against other models

The failure mode that hid this defect: every check in the project compared one
model to another. NIST vs GCGP vs SPGP, "% difference vs the Helmholtz
reference." Those are *relative* checks. All the methods ran through the same
flowsheet, so they shared the defect and moved together. The comparison looked
healthy because everything was wrong in the same direction.

Relative checks cannot detect a defect in shared machinery. You need at least
one check that does not use the machinery at all.

**Absolute checks, in rough order of value:**

| check | catches |
|---|---|
| Conservation (mass, energy, momentum) around every closed loop | streams that do not connect; sources and sinks that should not exist |
| Thermodynamic limits (Carnot, second law, equilibrium) | answers that are impossible rather than merely wrong |
| Continuity across every connection | the state leaving one unit differs from the state entering the next |
| Read specifications back out of the solved model | a specification that was silently ignored |
| Closed-form or analytic solution of a simplified case | everything, when one is available |

Most of these are three lines and cost nothing per run.

---

## Read the specification back out

Setting a specification and trusting it is the single most common way a
process model lies to you.

```python
vc.set_specifications(superheating=0, subcooling=0, ...)
# ... solve ...
# What actually came out:  superheat = +59.35 K,  subcool = -29.54 K
```

In this project `superheating=0` **deactivated** the constraint rather than
enforcing saturation. The guard read `if superheating > 0.1: activate()
else: deactivate()`, which looks like a sensible shortcut and is not.

There were even `assert superheating >= 0` lines nearby — which validate the
**input** and pass cleanly at zero, making the parameter look checked when
nothing checked what the model did with it.

**Rule: never assert on an input where you can assert on the solved state.**

Compare with the property side of the same project, which did this correctly:
`verify_0924_chain.py` reads `temperature_crit`, `pressure_crit`, `omega` and
the Cp coefficients back out of a *built* `GenericParameterBlock`, rather than
trusting the dict that was passed in. That discipline existed for properties
and not for the cycle. The defect lived in the half that had no read-back.

---

## Watch for mechanisms that hide symptoms

A retry loop that runs until the answer looks right is a filter that removes
evidence.

This project's `phase6` retried up to 10 times until the COP trend came out
monotonic. A run that converged to a visibly wrong answer got retried until it
stopped being visibly wrong. The loop was built to improve robustness and
functioned as a way to launder bad results.

Anything of this shape deserves suspicion:

- retry until the result passes a plausibility test
- fall back to a different starting point on failure, silently
- clip, clamp or bound a variable into a reasonable range
- discard non-converged points and report only the survivors

None are wrong in themselves. All of them hide information. Log what they
suppressed, and look at the log.

---

## Be suspicious of forgiving metrics

COP is a ratio of enthalpy *differences*. Mass flow cancels. Errors partly
cancel between numerator and denominator. A loop that failed to close by
1.3 % produced a COP comfortably inside the Carnot limit and looked fine, while
individual connections were off by 12–43 %.

**Check the raw quantities, not just the headline number.** If the reported
figure is a ratio, a difference, an efficiency or a normalised error, compute
the absolute quantities underneath and check those too.

---

## Pressure-test by making the problem harder

The defect here became undeniable only when the optimiser was turned loose.
With a square solve the model sat near a plausible answer; with an objective
active it found the loose variable and exploited it, landing 34–44 % above
Carnot.

**An optimiser is a defect detector.** If a model has an unconstrained degree
of freedom, an objective will find it. Running a square simulation and running
the same model as an optimisation are two different tests, and the second one
is much harder to pass.

Other ways to make it harder:

- widen the operating range until something breaks, then ask why there
- run the same case on a different platform or solver version
- perturb the initial guess and check the answer does not move
- turn on an objective that should not be improvable

---

## Do not let "pre-existing" close an investigation

The worst reasoning error in this episode was mine, and it is worth naming.

A 20–60 % disagreement with a closed-form reference turned up on a new fluid.
Checking the old fluid showed the same gap. I concluded "pre-existing, so the
new work is clean" and filed it as resolved — then kept building for the rest
of the day.

**"Not new" is not "not a problem."** Finding that a discrepancy is older means
the defect has been contaminating results for longer, not that it is
acceptable. A disagreement that size with an independent reference should stop
work on both the old and the new case.

Related traps that each got a plausible local explanation here:

- "it's a platform difference" (1–6 % spread across macOS and Linux)
- "it's solver noise" (non-monotonic trends)
- "it's environment, not the chain" (one ambient failing on one machine)

Each explanation was reasonable. None was tested against "what if the model is
wrong?" When several independent oddities each need their own excuse, the
excuses are probably all symptoms of one cause.

---

## Know what your connections actually enforce

The specific defect, as an example of a general class:

The flowsheet's arcs enforced `flow_mol`, `mole_frac_comp`, `pressure` and
`temperature`. For a **pure fluid at saturation** (T, P) are not independent —
they lie on the saturation line — so the vapour fraction was left free. Two
blocks could satisfy every arc equality while being in different
thermodynamic states:

```
block              T [C]   P [bar]   phase_frac Vap    h [J/mol]
valve_out        -28.503   2.83315         0.579633      -7051.0
evap_in          -28.503   2.83315         0.000003     -18127.3
```

The Helmholtz path in the same project used `StateVars.TPX`, where vapour
fraction *is* a state variable, and its loop closed to 0.00 W. Same flowsheet
topology, different state-variable set, sound versus broken.

**General rule: enumerate what each connection transmits, and check that it is
enough to determine the state.** A variable set that is sufficient in one
region (single phase) can be degenerate in another (two-phase, pure
component). The place to check is the hardest region the model will enter.

---

## A minimal standing check

Cheap, model-independent, and would have caught this in phase 4:

```python
# 1. physical limit
assert cop <= T_evap_K / (T_cond_K - T_evap_K), "COP exceeds Carnot"

# 2. conservation around the loop
assert abs(Q_evap + W_comp + Q_cond) < tol, "loop does not close"

# 3. continuity across every connection
for upstream, downstream in connections:
    assert abs(h_out[upstream] - h_in[downstream]) < tol, "stream discontinuity"

# 4. the specification you asked for is the one you got
assert abs(measured_superheat - requested_superheat) < tol
```

In this project these live in `R32/cycle_guards.py` as `check_cycle()`.

Put the guards in **before** attempting a fix, so the fix is verified by
something that already runs rather than by eyeballing output.

---

## Order of operations when something looks off

1. **Run the discriminating test first.** Do not try fifteen variations of a
   setting hoping one works. Find the single test that splits the hypothesis
   space in half and run that. (Sixteen solver configurations were tried here
   before the one diagnostic that showed the operating point was the problem.)
2. **Get an independent reference.** Closed form, hand calculation, a
   different tool, a published value. Without one you can only compare the
   model to itself.
3. **Check conservation and limits before anything subtle.** They are cheap
   and they fail loudly.
4. **Read the solved state, not the inputs.**
5. **Work outward from the smallest failing case.** One method, one operating
   point, one unit.
6. **When you find the cause, check whether it explains *all* the symptoms.**
   If three oddities were waved off earlier, go back and confirm the cause
   accounts for each. If it does not, there is a second defect.

---

## What was and was not affected, as an illustration

Worth noting how the blast radius fell out, because it shows where the
verification effort had gone:

| | status | why |
|---|---|---|
| Property predictions, p-h / T-s / saturation figures | **sound** | standalone, no flowsheet, no solver; validated against a datasheet and against a paper's own supporting information |
| Helmholtz cycle results | energy-consistent | vapour fraction was a state variable, so connections held |
| Cubic (PR) cycle results | **unusable** | connections did not transmit state |

The half of the project with three verification scripts and a cross-check
against published reference code was fine. The half with none was not. That
is not a coincidence, and it is the whole argument for this file.
