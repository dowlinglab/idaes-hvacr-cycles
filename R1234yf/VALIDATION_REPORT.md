# R1234yf Helmholtz EOS Property Package — Validation Report

**Component:** R1234yf (2,3,3,3-tetrafluoropropene, HFO-1234yf)
**Equation of state:** Lemmon, E.W., Akasaka, R. (2022), "Fundamental equation of state for 2,3,3,3-tetrafluoropropene (HFO-1234yf)," Int. J. Thermophys. 43, 171. https://doi.org/10.1007/s10765-022-03206-9
**Files covered:** `R1234yf.py` (EOS implementation), `R1234yf_validation.py` (saturation/critical-point solver and p-H diagram driver)
**Date:** August 17, 2026

## Bottom line

The property package's critical point and reference state are validated against independent literature sources, and both checks pass within stated experimental uncertainty. The saturation solver converges everywhere it was tested. The generated p-H diagram is internally self-consistent and visually matches the real Danfoss Coolselector®2 reference chart's shape and isentrope values. **This codebase has now also been checked point-by-point against 135 independent experimental measurements** (30 vapor-pressure + 105 p-ρ-T points from Richter, McLinden, Lemmon (2011), read directly from the primary source) — agreement is within a few tenths of a percent almost everywhere, with the one region needing a caveat (near-critical density) explained physically in Section 6, not hand-waved away.

"Validated" below is used narrowly, for exactly the claims each section supports — not as a blanket label for the whole package.

## 1. Critical point

The critical point was solved directly from the EOS by finding the state where `dP/d`delta` = 0` and `d²P/d`delta`² = 0` (the two defining conditions of a true critical point), via a 2D nonlinear least-squares solve. Result:

| Quantity | Solved value |
|---|---|
| delta_c (reduced density) | 1.000000 |
| tau_c (reduced inverse temperature) | 1.000000 |
| T_c | 367.8500 K |
| rho_c | 476.6900 kg/m^3 |
| P_c | 3,384,349.91 Pa (3.38435 MPa) |

**Comparison against the EOS's own published fit** (Lemmon & Akasaka Table 1: Tc=367.85 K, rhoc=476.69 kg/m^3, Pc=3,384,400 Pa) — matches to within 0.0015% on pressure. This check is **circular** (the EOS was fit to reproduce these exact numbers) and only confirms the solver itself is accurate, not that the EOS is physically correct.

**Comparison against Tanaka & Higashi (2010)** — an independently measured critical point (visual meniscus-disappearance method, not from this EOS's own fit): Tc = 367.85 ± 0.01 K, rhoc = 478 ± 3 kg/m^3, Pc = 3.382 ± 0.003 MPa.

| Quantity | Solved | T&H (2010) | T&H uncertainty band | Solved value within band? |
|---|---|---|---|---|
| T_c | 367.8500 K | 367.85 K | 367.84 – 367.86 K | Yes |
| rho_c | 476.69 kg/m^3 | 478 kg/m^3 | 475 – 481 kg/m^3 | Yes |
| P_c | 3.38435 MPa | 3.382 MPa | 3.379 – 3.385 MPa | Yes |

All three solved critical-point values fall inside Tanaka & Higashi's stated experimental uncertainty. This is a genuine, non-circular validation of the EOS's critical point.

## 2. Reference state

The property package's `enthalpy()`/`entropy()` methods carry an additive offset (`h_offset`, `s_offset`), currently both set to `0.0`. Rather than assume this is correct, it was checked directly: the true saturated-liquid state at 273.15 K (0°C) was solved using the saturation solver described in Section 3, and the raw (`offset=0`) enthalpy/entropy were evaluated there.

| Quantity | Computed (offset=0) | IIR reference target | Relative error |
|---|---|---|---|
| h at 273.15 K sat. liquid | 200,000.2257 J/kg | 200,000 J/kg | 0.00011% |
| s at 273.15 K sat. liquid | 1000.0011 J/(kg·K) | 1000 J/(kg·K) (1.00 kJ/(kg·K)) | 0.00011% |

This confirms `h_offset = s_offset = 0.0` is already correct — the paper's own ideal-gas constants happen to reproduce the industry-standard IIR reference convention (h=200 kJ/kg, s=1.00 kJ/(kg·K) at 0°C saturated liquid) without any additional rebasing. `enthalpy()`/`entropy()` outputs are directly comparable to manufacturer charts and other IIR-convention data with no offset correction needed.

## 3. Saturation dome solver

Saturated liquid/vapor densities and pressure are solved at each temperature via a 4-step method: (1) scan for the true stability limits (spinodals) at that temperature, (2) solve each phase's density independently within its own stability region (structurally prevents the two phases from ever converging to the same state), (3) adjust a trial pressure until both phases' Gibbs energies agree (the Maxwell equal-area condition), (4) assemble the converged result.

**Convergence:** 100 of 100 temperature points converged, spanning T = 200 K up to T_c − 0.5 K = 367.35 K.

Sample converged values:

| T (K) | rho_liq (kg/m^3) | rho_vap (kg/m^3) | P_sat (Pa) |
|---|---|---|---|
| 200.0 | 1378.8 | 0.62 | 8,960 |
| 250.0 | 1245.3 | 7.7 | 132,638 |
| 300.0 | 1085.9 | 39.6 | 713,539 |
| 340.0 | 900.2 | 121.6 | 1,922,913 |
| 367.35 | 572.9 | 386.8 | 3,350,212 |

An earlier version of this solver had a defect (a trivial mathematical solution where the "solved" liquid and vapor states silently collapsed to the same density everywhere) that was found, diagnosed, and fixed before this validation; the numbers above are from the corrected solver.

## 4. p-H diagram construction vs. the Danfoss reference chart

A full pressure–enthalpy diagram was built and compared against a real reference chart (Danfoss Coolselector®2, v3.3.1, "Aserep" database v3.5.0). The diagram includes:

- The saturation dome (bubble/dew curves), spliced through the solved critical point.
- Constant-quality (vapor-fraction) lines at x = 0.1 through 0.9.
- Isotherms at T = 230–360 K (10 K steps, subcritical) and T = 370–420 K (10 K steps, supercritical).
- Isentropes at 38 entropy values (775–1575 J/(kg·K) in 50 J/(kg·K) steps, and 1625–2125 J/(kg·K) in 25 J/(kg·K) steps) — these exact values were read directly off the uploaded Danfoss chart's own isentrope labels, not chosen arbitrarily.
- Isochores (constant-density lines) are present on the real Danfoss chart but were deliberately left out of this reproduction, by request.

**On the isentrope values specifically:** because these 38 numbers were transcribed from the Danfoss chart and used as direct inputs (targets to solve for), their presence on both charts is by construction, not an independent cross-check. What *is* an independent check is that the property package's own entropy calculation, using no adjustable parameters beyond what's described in Sections 1–2, produced convergent, physically sensible (monotonic, smoothly varying) curves at every one of those 38 target values, including several very close to the critical entropy where the solve is numerically difficult (see Section 5).

**Visual comparison:** the constructed diagram reproduces the Danfoss chart's qualitative shape — dome shape and peak location, isotherms fanning densely on the vapor/supercritical side while running steep and nearly vertical through the liquid side, and isentropes following the same pattern (near-vertical through the liquid region, fanning through the vapor/supercritical region). The plotted axis range was deliberately extended above the Danfoss chart's own 3 MPa ceiling (to 3.6 MPa) so the full critical point remains visible — this is an intentional deviation from the reference chart's own cropping, not a mismatch.

## 5. Bugs found and fixed during this validation effort

For transparency, three defects were found and corrected while building and checking this diagram (none are still present in the current files):

1. **Saturation solver trivial-solution collapse.** The original 2D joint solve for liquid/vapor density had a mathematical trap where comparing a state to itself trivially satisfied the solver's own equations. Fixed by restructuring into the 4-step method described in Section 3, which makes that trap structurally unreachable.
2. **Isentrope two-phase segment plotted in reverse order**, which produced a visible double-line/zigzag artifact where an isentrope crossed the two-phase dome. Fixed by correcting the point ordering before joining it to the rest of the curve.
3. **Isentrope solver silently returning wrong answers very close to the critical point.** The solver reported "converged" while actually being off by roughly 100 J/(kg·K) from the intended target — `sol.success` from the underlying numerical solver means it terminated normally, not that it found the correct answer, and that distinction was being missed. Fixed by explicitly checking that the solved state's entropy actually matches the target (not just trusting the solver's own success flag), and by trying alternate starting guesses for entropies close to the critical value, since the correct continuation there passes through the supercritical region rather than compressed liquid.

## 6. Point-by-point comparison against independent NIST experimental data

This section was added after obtaining the actual primary-source paper: Richter, M., McLinden, M.O., Lemmon, E.W. (2011), "Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf): Vapor Pressure and p-ρ-T Measurements and an Equation of State," J. Chem. Eng. Data 56, 3254–3264. This paper's own tables (Table 1: vapor pressure; Tables 2–3: p-ρ-T density) were read directly from the PDF, not summarized secondhand.

**Important distinction:** this 2011 paper fits its *own* 15-term equation of state (its Table 6), which is different from the 17-term Lemmon & Akasaka (2022) EOS this codebase implements (it also adopts slightly different critical constants: rho_c ≈ 475.55 kg/m^3, P_c = 3382.2 kPa, vs. this codebase's solved 476.69 kg/m^3 / 3384.35 kPa). The paper's own printed deviation columns (Δp, Δρ) are relative to *its* EOS, not this codebase's. To keep the comparison non-circular and directly relevant, the raw experimental (T, P) and (T, P, ρ) measurements were taken from the tables and compared directly against this codebase's own computed values — not against the paper's reported deviations.

### 6.1 Vapor pressure (Table 1, 30 points, T = 250–366 K)

For each experimental temperature, `solve_saturation_at_tau` computed P_sat independently; compared against the paper's measured value.

| Statistic | This codebase vs. Richter et al. (2011) data | Paper's own EOS vs. its own data (for reference) |
|---|---|---|
| Mean absolute relative error, all 30 points | 0.056% | — |
| Std. dev., all 30 points | 0.120% | 0.11% |
| Mean absolute relative error, T ≥ 270 K (26 points) | 0.027% | — |
| Std. dev., T ≥ 270 K | 0.040% | 0.06% |
| Worst single point | 0.541% at T = 250.002 K | — |

This codebase's agreement is as good as or slightly better than the paper's own EOS fit against its own data, despite comparing against a newer/different EOS entirely. The single outlier (T = 250.002 K, the coldest point measured) is still under 0.6% — well within engineering tolerance, and consistent with the paper's own Figure 2, which shows increasing scatter (across all compared EOS/literature datasets, not just this one) toward the low-T end of the measured range.

### 6.2 p-ρ-T density (Tables 2–3, 105 points, T = 232–400 K, compressed liquid through supercritical)

For each experimental (T, P) pair, density was solved directly from the EOS (`pressure(delta, tau) = P_target`, root-found for delta) and compared against the paper's measured density.

| Region | n points | Mean absolute relative error | Std. dev. | Max |
|---|---|---|---|---|
| Normal region (ρ outside 285–761 kg/m^3) | 93 | 0.045% | 0.072% | 0.277% |
| Near-critical region (285 < ρ < 761 kg/m^3) | 12 | 0.486% (density basis) | 0.627% | 1.456% |
| Near-critical region, re-checked on a **pressure** basis instead | 12 | 0.169% | 0.184% | 0.473% |

The near-critical region needs a caveat, and the paper explains why: within roughly 285–761 kg/m^3 of the critical density (476.69 kg/m^3), compressibility (∂P/∂ρ) approaches zero as the critical point is approached, so a tiny, entirely normal pressure difference between two slightly-different EOS fits corresponds to a much larger *density* difference at the same pressure. The paper's own Figure 5 handles this by comparing *pressure* deviations in this region instead of density deviations — re-running the comparison the same way (evaluating this codebase's pressure at the reported (T, ρ) point and comparing to the reported pressure) brings the near-critical agreement back down to 0.169% mean / 0.473% max, consistent with the rest of the data set. The larger density-basis numbers above are not a defect; they are the expected, physically-explained consequence of comparing density right at the point where density stops being a sensitive measure of agreement.

### 6.3 Conclusion of this section

Across 135 independent experimental data points (30 vapor pressure + 105 density), spanning T = 232–400 K and P up to nearly 10 MPa, this codebase's EOS implementation agrees with the primary experimental measurements to within a few tenths of a percent in the well-conditioned regions, and within a few tenths of a percent in the near-critical region once compared on the physically appropriate (pressure) basis. This is now a genuine, non-circular, point-by-point numeric validation against independent experimental data — not merely a visual or aggregate-statistic comparison.

## 7. What is still NOT validated (remaining open items)

- **The paper's own Table 7 verification points** (from the separate Lemmon & Akasaka 2022 paper — published (T, rho) → (P, cv, cp, w) values) have not been re-run against the current code.
- **No pixel-level or coordinate-level overlay** of the generated diagram against the Danfoss chart has been performed — the comparison in Section 4 is visual/qualitative (shape, density of lines, general pattern), not a numeric overlay.
- **Heat capacity and speed of sound** are implemented in `R1234yf.py` but are not exercised by this diagram, this section's comparison, or this report at all; they are outside its scope. (The Richter et al. 2011 paper does report speed-of-sound and heat-capacity comparisons against its own EOS, but those were not re-run against this codebase.)
- The 105 p-ρ-T points and 30 vapor-pressure points used here are the paper's own tabulated *averages* of 4–8 replicate measurements each (per the paper's own Table 1/2/3 notes) — the full underlying replicate-level dataset (216 vapor-pressure replicates, 557 density replicates) was not obtained or used.

## Verdict

Calling the **critical point and reference state** "validated" is accurate and supported by independent literature comparison (Section 1–2). Calling the **saturation solver, EOS implementation, and diagram construction** "validated" is now supported by a genuine point-by-point comparison against 135 independent experimental measurements (Section 6), in addition to being internally self-consistent and bug-free as currently known (Sections 3–5) and visually consistent with the real Danfoss reference chart. The specific items in Section 7 remain open and are not covered by this claim.
