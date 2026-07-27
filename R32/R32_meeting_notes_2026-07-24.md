# R-32 Modeling — Meeting Notes, 2026-07-24

## Background: what is a property package?

A property package is the piece of a process model that answers "given a
fluid's temperature, pressure, and composition, what are its physical
properties (density, enthalpy, entropy, phase behavior, etc.)?" Every unit
operation in a flowsheet (evaporator, compressor, condenser, valve) calls
into the property package to know what state the fluid is actually in at
each point. Swapping the property package (e.g. from a fluid-specific
Helmholtz-energy formulation to a general cubic equation of state) changes
how those properties are predicted, but not the structure of the cycle
model itself — that's why it's useful as a way to test how sensitive a
result like COP is to the underlying thermodynamic model.

## History (results presented before 07/15)

- Predictions up to that point used the **Ambrose-Walton correlation** for
  Psat(T) — a simple, direct fit for saturation pressure vs. temperature.
  It only applies to pure components, so a different approach would
  eventually be needed for mixtures.
- Next step at the time: interface Shilpa's new equation-of-state (EoS)
  models with the cycle (DVCT) to estimate COP, rather than relying on the
  Ambrose-Walton correlation alone.

## Today's update: NIST vs. Helmholtz validation

Reporting **NIST only** today, compared against the existing Helmholtz
reference model. GCGP and SPGP are not included in today's report (see
Next Steps).

- Built and validated a cubic (Peng-Robinson) equation-of-state property
  package for R-32 in IDAES, using the NIST-fitted critical properties
  (Tc, Pc, acentric factor).
- Assembled the full four-unit vapor-compression cycle (evaporator,
  compressor, condenser, expansion valve) and got it solving end-to-end
  with this new cubic EoS.
- Computed cycle COP across ambient temperatures 10–25°C, using the ideal
  vapor-compression cycle assumption (isentropic compression, no
  superheat/subcooling) from **Raskar & Mutalikdesai (2016), International
  Journal of Current Engineering and Technology (IJCET), Vol. 6, Issue 5**.
- **Validation result:** the NIST-based cubic EoS model tracks the
  Helmholtz reference model within about 2–10% across the ambient range
  (closest agreement at 20°C, within 2%), confirming the new cubic-EoS
  modeling approach reproduces the established reference model's cycle
  performance.

## Next steps

- Extend the same COP comparison to the other two published R-32
  parameter sets, **GCGP and SPGP** — not reported today.
- One open numerical issue to resolve first: at one ambient condition, the
  compressor's internal reference-state calculation is landing on an
  inconsistent solution for GCGP; this is being isolated before GCGP
  results are reported.
- SPGP's fitted parameters produce a thermodynamically invalid equation
  of state (negative acentric factor); expect to report this as a
  structural limitation rather than a numerical bug once GCGP is closed
  out.
