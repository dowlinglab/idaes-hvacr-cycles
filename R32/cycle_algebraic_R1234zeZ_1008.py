"""
cycle_algebraic_1001.py -- the ideal vapor-compression cycle computed
ALGEBRAICALLY, with no IDAES, no Pyomo and no IPOPT.

Why (2026-09-17): the Tc/Pc sensitivity sweep goes through the full IDAES
flowsheet, which solves every unit balance, property correlation and phase-
equilibrium relation simultaneously as one large nonlinear system. That is
the right tool for the DVCT deliverable, but for a sensitivity study it
injects solver behaviour into the answer: convergence flags that depend on
the starting point, a loosened tolerance (tol=1e-4, see vapor_compression_
cubic_refstate_0903.py line 1633) that sets a floor on resolvable COP
differences, and 251 of 900 case-ambient solves reported non-converged.

But with the evaporator and condenser saturation temperatures fixed and
superheat = subcool = 0, this cycle is NOT a simultaneous system. It is
sequential:

    P_evap  = Psat(T_evap)
    state 1 = saturated vapour at T_evap          -> h1, s1
    P_cond  = Psat(T_cond)
    state 2 = compress isentropically to P_cond   -> h2s, then h2 via eta
    state 3 = saturated liquid at T_cond          -> h3
    state 4 = isenthalpic expansion to P_evap     -> h4 = h3
    COP     = (h1 - h4) / (h2 - h1)

Only step 2 is implicit, and only in one variable: find T2 such that
s_vap(T2, P_cond) = s1. That is a 1-D root find on a monotone function.

Consequences for the sensitivity work:
  - deterministic: one answer per (Tc, Pc), every time. No convergence
    flag, no path dependence, no warmstart, no tolerance floor -- so the
    small deviations (+/-0.01%, +/-0.1%) become meaningful instead of
    measuring solver residual.
  - fast: the whole 15x15 grid runs in well under a second.
  - fails only where the physics genuinely fails (e.g. a condenser
    temperature at or above the perturbed Tc, where no saturation state
    exists), which is a real result rather than a solver shrug.

Physics is the same as the IDAES package and as compare_cp_methods_*.py:
Peng-Robinson cubic EoS, Ambrose-Walton saturation pressure (Poling
sec. 7-4), PR departure functions for h and s (Poling Table 6-3), and
Shomate ideal-gas Cp integrated for the ideal-gas h and s (t = T/1000).

Reference state does not matter here: COP is a ratio of enthalpy
DIFFERENCES, so any constant offset cancels exactly. No IIR anchoring
needed.

VALIDATION (run this file; it checks itself): NIST baseline COPs are
compared against the IDAES pipeline's own converged results, and the
saturated-dome pressures against the Linde datasheet.

Added 2026-09-21: curve-generation and curve-shape-comparison utilities
(ts_curve, ph_curve, curve_area, curve_shape_metrics), for the sensitivity
analysis redesign requested by the Colon group -- COP error alone can look
acceptable while the underlying T-S/P-H behaviour is not, so these give a
second, independent way to score a property-prediction method: how closely
its saturated-liquid/vapour curves track a baseline (NIST) fluid's, in
both a pointwise (mean absolute error) and a shape (area-under-curve) sense.
These do not change run_cycle or Fluid in any way -- T_evap/T_cond are
still always passed in explicitly by the caller, never solved for.

Also fixed 2026-09-21: the Linde datasheet's 10 C saturation-pressure
entry in the validation table below was 10.065 bar, a transposed-digit
typo -- the real Linde datasheet (finer-resolution scan, confirmed
against enthalpy/entropy columns that were already correct) gives
11.069 bar, matching CoolProp (11.0691 bar) to four significant figures.
Corrected here; Validation 1's 10 C error dropped from 10.69% to 0.65%,
in line with every other point. The same wrong value is still present
in ~15 other files across this project (sensitivity_algebraic_0917.py,
compare_cp_methods_*.py, phase_2_saturation_dome.py, pr_eos_lib.py,
PR_EOS.py) -- not yet propagated there.


Added 2026-09-24 (this _0924 copy): the six datasets from the Colon group's
09/19 delivery ("spgp r32.xlsx", the yellow-highlighted blocks under the red
2026-09-19 band), plus a COP sweep and a Linde saturation comparison across
every method. 0921 had only NIST and GCGP defined and no method sweep at all
-- its __main__ was a self-validation harness, so new fluids had nothing to
run through. Both additions are purely additive: Fluid, run_cycle and the
curve utilities are untouched, and Validations 1 and 2 still run first and
still compare NIST/GCGP exactly as before.

The 09/10 "(7).xlsx" fits are NOT carried into this file. That delivery is a
separate branch of the collaborators' work rather than an earlier version of
the 09/19 one -- neither of its pvap values (7853.748, 1670.067 mmHg) appears
anywhere in the 09/19 sheet, while the 09/06 value 7409.643 does.

Two conversions were needed to get the sheet's numbers into Fluid; both are
documented at the fluid definitions below. In short: the sheet's Cp fit is in
the raw T-polynomial convention and has to be scaled into Shomate's t = T/1000
form, and the sheet supplies no acentric factor at all -- only "pvap (mmHg)",
which is Psat at Tr = 0.7, from which omega is derived here.

Author: Shilpa Narasimhan. Support: Claude AI.

Date created: 09/24/2026 (copied from cycle_algebraic_0921.py)

--- 10/08/2026 -- INCOMPLETE -- READ THIS BEFORE USING ---
This file is NOT retargeted to R-1234ze(Z). Only the molar mass
was changed. Every Fluid row below still carries R-32 critical
properties and R-32 Shomate coefficients, and the __main__ block
still sweeps the R-32 duty (T_evap = -29 C, ambients 10/15/20/25).
Running it as-is gives R-32 numbers with an R-1234ze(Z) molar mass,
which is meaningless.

It is shipped only so the copy set is complete. The algebraic
PR-vs-Helmholtz comparison for this fluid already exists and is
verified -- use compare_cp_methods_R1234zeZ_1008.py instead.

To finish this file: replace the Fluid rows with the REFPROP PR row
from phase_1_cubic_eos_validation_refstate_R1234zeZ_1008.py, and set
the __main__ sweep to T_evap = 30 C, T_cond = 70/80/90/100 C.
"""

import numpy as np
import pandas as pd
import os
from scipy.optimize import brentq

R = 8.314          # J/mol/K  -- matches phase_1_cubic_eos_validation_refstate_0903.py
MW = 1.140416e-01     # kg/mol, R-1234ze(Z)
P_REF = 1.0e5      # Pa, only used inside s^ig; cancels out of COP


# =============================================================================
# Fluid definition
# =============================================================================

class Fluid:
    """One property parameter set: critical properties, acentric factor and
    Shomate ideal-gas Cp coefficients (A-E on the t = T/1000 convention, i.e.
    already rescaled the same way phase_1_cubic_eos_validation_refstate_*.py
    stores them)."""

    def __init__(self, name, Tc, Pc, omega, A, B, C, D, E):
        self.name = name
        self.Tc = float(Tc)          # K
        self.Pc = float(Pc)          # Pa
        self.omega = float(omega)
        self.A, self.B, self.C, self.D, self.E = map(float, (A, B, C, D, E))
        self.kappa = 0.37464 + 1.54226*self.omega - 0.26992*self.omega**2

    # ---- PR EoS ------------------------------------------------------------

    def alpha(self, T):
        return (1.0 + self.kappa*(1.0 - np.sqrt(T/self.Tc)))**2

    def a(self, T):
        return 0.45724 * R**2 * self.Tc**2 * self.alpha(T) / self.Pc

    def dadT(self, T):
        """d a/dT. From alpha = [1 + k(1 - sqrt(T/Tc))]^2,
        dalpha/dT = -k*sqrt(alpha)/sqrt(T*Tc)."""
        dalpha = -self.kappa*np.sqrt(self.alpha(T))/np.sqrt(T*self.Tc)
        return 0.45724 * R**2 * self.Tc**2 * dalpha / self.Pc

    @property
    def b(self):
        return 0.07780 * R * self.Tc / self.Pc

    def z_roots(self, T, P):
        """Real roots of the PR cubic in Z, ascending."""
        A_ = self.a(T)*P/(R*T)**2
        B_ = self.b*P/(R*T)
        coeffs = [1.0,
                  -(1.0 - B_),
                  A_ - 2.0*B_ - 3.0*B_**2,
                  -(A_*B_ - B_**2 - B_**3)]
        roots = np.roots(coeffs)
        real = np.sort(roots[np.abs(roots.imag) < 1e-9].real)
        real = real[real > B_]           # physically admissible only
        if real.size == 0:
            raise ValueError(f"{self.name}: no admissible Z root at T={T:.2f} K, P={P:.1f} Pa")
        return real

    def z_vap(self, T, P):
        return self.z_roots(T, P)[-1]

    def z_liq(self, T, P):
        return self.z_roots(T, P)[0]

    # ---- departure functions (Poling Table 6-3, PR: u=2, w=-1) -------------

    def _dep_log(self, Z, B_):
        s2 = np.sqrt(2.0)
        return np.log((Z + (1.0 + s2)*B_) / (Z + (1.0 - s2)*B_))

    def h_dep(self, T, P, Z):
        """h - h_ideal_gas, J/mol."""
        B_ = self.b*P/(R*T)
        s2 = np.sqrt(2.0)
        return (R*T*(Z - 1.0)
                + (T*self.dadT(T) - self.a(T))/(2.0*s2*self.b)*self._dep_log(Z, B_))

    def s_dep(self, T, P, Z):
        """s - s_ideal_gas(T, P), J/mol/K."""
        B_ = self.b*P/(R*T)
        s2 = np.sqrt(2.0)
        return (R*np.log(Z - B_)
                + self.dadT(T)/(2.0*s2*self.b)*self._dep_log(Z, B_))

    # ---- ideal-gas Shomate integrals (t = T/1000) --------------------------

    def h_ig(self, T):
        """Ideal-gas enthalpy, J/mol, to an arbitrary constant (cancels in COP)."""
        t = T/1000.0
        return 1000.0*(self.A*t + self.B*t**2/2.0 + self.C*t**3/3.0
                       + self.D*t**4/4.0 - self.E/t)

    def s_ig(self, T, P):
        """Ideal-gas entropy at (T, P), J/mol/K, to an arbitrary constant."""
        t = T/1000.0
        s_T = (self.A*np.log(t) + self.B*t + self.C*t**2/2.0
               + self.D*t**3/3.0 - self.E/(2.0*t**2))
        return s_T - R*np.log(P/P_REF)

    # ---- saturation --------------------------------------------------------

    def psat(self, T):
        """Ambrose-Walton saturation pressure, Pa (Poling sec. 7-4)."""
        Tr = T/self.Tc
        if Tr >= 1.0:
            raise ValueError(f"{self.name}: T={T:.2f} K is at or above Tc={self.Tc:.2f} K "
                             f"-- no saturation state exists")
        tau = 1.0 - Tr
        f0 = (-5.97616*tau + 1.29874*tau**1.5 - 0.60394*tau**2.5 - 1.06841*tau**5)/Tr
        f1 = (-5.03365*tau + 1.11505*tau**1.5 - 5.41217*tau**2.5 - 7.46628*tau**5)/Tr
        f2 = (-0.64771*tau + 2.41539*tau**1.5 - 4.26979*tau**2.5 + 3.25259*tau**5)/Tr
        return self.Pc*np.exp(f0 + self.omega*f1 + self.omega**2*f2)

    # ---- state properties --------------------------------------------------

    def h_vap(self, T, P):
        return self.h_ig(T) + self.h_dep(T, P, self.z_vap(T, P))

    def s_vap(self, T, P):
        return self.s_ig(T, P) + self.s_dep(T, P, self.z_vap(T, P))

    def h_liq(self, T, P):
        return self.h_ig(T) + self.h_dep(T, P, self.z_liq(T, P))

    def s_liq(self, T, P):
        return self.s_ig(T, P) + self.s_dep(T, P, self.z_liq(T, P))

    def sat_vapour(self, T):
        P = self.psat(T)
        return P, self.h_vap(T, P), self.s_vap(T, P)

    def sat_liquid(self, T):
        P = self.psat(T)
        return P, self.h_liq(T, P), self.s_liq(T, P)


# =============================================================================
# Curve generation -- T-S and P-H diagrams for shape comparison
# =============================================================================

# TS_GRID_C must be defined before the functions below, since it is used
# as a default argument value -- Python evaluates default arguments once,
# when the `def` line itself runs, not when the function is later called.
# Same 20 points as the Linde datasheet table used throughout this project
# (sensitivity_algebraic_0917.py's LINDE_SAT, compare_cp_methods_*.py's
# error tables), reused here so curve comparisons stay on the same footing
# as the existing P_MAPE/hf_MAE/etc. conventions.
TS_GRID_C = [-130, -110, -90, -70, -50, -30, -10, 0, 10, 20, 30, 40, 50, 58, 62, 66, 70, 74, 76, 78]


def ts_curve(fluid, T_grid_C=TS_GRID_C):
    """Saturated-liquid and saturated-vapour entropy of `fluid` over the
    shared temperature grid -- the raw data behind a T-S diagram.

    Temperatures at or above fluid.Tc are skipped (no saturation state
    exists there, same guard as Fluid.psat). Since T_grid_C is ascending,
    once a temperature is skipped every later one is too, so the returned
    arrays are always a temperature-ordered PREFIX of the full grid, never
    a scattered subset -- this is what makes it safe for
    curve_shape_metrics() to line two fluids' curves up by simple
    truncation below.

    Returns (T_K, sf, sg): temperature [K], saturated-liquid entropy
    [J/mol/K], saturated-vapour entropy [J/mol/K].
    """
    T_K, sf, sg = [], [], []
    for T_C in T_grid_C:
        T = T_C + 273.15
        if T >= fluid.Tc:
            continue
        _, _, s_liq = fluid.sat_liquid(T)
        _, _, s_vap = fluid.sat_vapour(T)
        T_K.append(T)
        sf.append(s_liq)
        sg.append(s_vap)
    return np.array(T_K), np.array(sf), np.array(sg)


def ph_curve(fluid, T_grid_C=TS_GRID_C):
    """Saturation pressure and saturated-liquid/vapour enthalpy of `fluid`
    over the shared temperature grid -- the raw data behind a P-H diagram.
    Same Tc-skip behaviour and prefix guarantee as ts_curve() above.

    Returns (T_K, P, hf, hg): temperature [K], saturation pressure [Pa],
    saturated-liquid enthalpy [J/mol], saturated-vapour enthalpy [J/mol].
    """
    T_K, P_arr, hf, hg = [], [], [], []
    for T_C in T_grid_C:
        T = T_C + 273.15
        if T >= fluid.Tc:
            continue
        P, h_liq, _ = fluid.sat_liquid(T)
        _, h_vap, _ = fluid.sat_vapour(T)
        T_K.append(T)
        P_arr.append(P)
        hf.append(h_liq)
        hg.append(h_vap)
    return np.array(T_K), np.array(P_arr), np.array(hf), np.array(hg)


def curve_area(x, y):
    """Trapezoidal area under y(x). Returns NaN for fewer than 2 points
    (e.g. a fluid whose Tc is so low almost every grid point is skipped)
    rather than raising, so a single bad case does not stop a sweep."""
    if len(x) < 2:
        return float("nan")
    return float(np.trapz(y, x))


def curve_shape_metrics(fluid, baseline, T_grid_C=TS_GRID_C):
    """How much `fluid`'s T-S and P-H curves differ from `baseline`'s
    (e.g. a perturbed Tc/Pc case vs. the NIST reference), as a full-range,
    fluid-level metric -- independent of any single cycle/ambient.

    Two kinds of disagreement are reported, since they catch different
    failure modes:
      - *_mae: pointwise mean absolute error (liquid + vapour branches
        pooled). Good at catching a curve that is uniformly shifted off
        the baseline.
      - *_area_diff: absolute difference in the area under each branch.
        Good at catching a curve that wiggles above and below the
        baseline -- pointwise errors can partly cancel in the mean, but
        the enclosed area still shows the shape is wrong. This is the
        metric the Colon group's "COP error alone can hide bad
        thermodynamic behaviour" concern was specifically about.

    fluid and baseline may have different Tc, so their curves (from
    ts_curve/ph_curve) can be different lengths; both are truncated to
    the shorter one. This is a valid match, not just an approximation --
    ts_curve/ph_curve always return a temperature-ordered PREFIX of the
    shared grid (see ts_curve docstring), so truncating to the shorter
    length keeps the same, correctly-matched leading temperatures on both
    sides.

    Returns a dict: TS_mae [J/mol/K], PH_mae [J/mol], TS_area_diff,
    PH_area_diff, and n_points (the number of matched grid points actually
    used -- 0 if the two fluids share no valid temperature at all).
    """
    T_f, sf_f, sg_f = ts_curve(fluid, T_grid_C)
    T_b, sf_b, sg_b = ts_curve(baseline, T_grid_C)
    n = min(len(T_f), len(T_b))
    T_ts = T_f[:n]
    sf_f, sg_f, sf_b, sg_b = sf_f[:n], sg_f[:n], sf_b[:n], sg_b[:n]

    _, _, hf_f, hg_f = ph_curve(fluid, T_grid_C)
    _, _, hf_b, hg_b = ph_curve(baseline, T_grid_C)
    n2 = min(len(hf_f), len(hf_b))
    T_ph = T_f[:n2]
    hf_f, hg_f, hf_b, hg_b = hf_f[:n2], hg_f[:n2], hf_b[:n2], hg_b[:n2]

    if n == 0 or n2 == 0:
        return dict(TS_mae=np.nan, PH_mae=np.nan,
                     TS_area_diff=np.nan, PH_area_diff=np.nan, n_points=0)

    TS_mae = float(np.mean(np.abs(np.concatenate([sf_f - sf_b, sg_f - sg_b]))))
    PH_mae = float(np.mean(np.abs(np.concatenate([hf_f - hf_b, hg_f - hg_b]))))

    TS_area_diff = (abs(curve_area(T_ts, sg_f) - curve_area(T_ts, sg_b))
                     + abs(curve_area(T_ts, sf_f) - curve_area(T_ts, sf_b)))
    PH_area_diff = (abs(curve_area(T_ph, hg_f) - curve_area(T_ph, hg_b))
                     + abs(curve_area(T_ph, hf_f) - curve_area(T_ph, hf_b)))

    return dict(TS_mae=TS_mae, PH_mae=PH_mae,
                 TS_area_diff=TS_area_diff, PH_area_diff=PH_area_diff,
                 n_points=n)








# =============================================================================
# The cycle -- sequential, no simultaneous solve
# =============================================================================

def run_cycle(fluid, T_evap_C=-29.0, T_amb_C=20.0, cond_approach=9.0,
              eta_isentropic=0.9999, superheat=0.0, subcool=0.0):
    """Ideal vapour-compression cycle, computed explicitly.

    Returns a dict with the four state points and the COP, or a dict with
    'error' set if the case is physically impossible (e.g. condenser
    temperature at or above Tc).

    superheat/subcool are in K and default to 0, matching the ideal-cycle
    specification used throughout this project (Raskar & Mutalikdesai,
    IJCET 6(5), 2016). They are exposed so the effect of moving slightly
    off the saturation dome can be tested.
    """
    T_evap = T_evap_C + 273.15
    T_cond = T_amb_C + cond_approach + 273.15

    try:
        P_evap = fluid.psat(T_evap)
        P_cond = fluid.psat(T_cond)

        # 1) evaporator outlet: saturated (or slightly superheated) vapour
        T1 = T_evap + superheat
        h1 = fluid.h_vap(T1, P_evap)
        s1 = fluid.s_vap(T1, P_evap)

        # 2) compressor outlet: isentropic to P_cond, then apply efficiency.
        #    s_vap(T, P_cond) increases monotonically with T, so this is a
        #    clean 1-D root find -- the only implicit step in the whole cycle.
        def ds(T):
            return fluid.s_vap(T, P_cond) - s1

        T_lo, T_hi = T_cond + 0.01, T_cond + 400.0
        if ds(T_lo) > 0:
            # already above the target entropy at the bubble point: walk down
            T_lo = T_cond - 60.0
        if ds(T_lo)*ds(T_hi) > 0:
            return {"error": "isentropic discharge temperature not bracketed"}
        T2s = brentq(ds, T_lo, T_hi, xtol=1e-8, rtol=1e-12)
        h2s = fluid.h_vap(T2s, P_cond)
        h2 = h1 + (h2s - h1)/eta_isentropic

        # actual discharge temperature for the non-ideal case
        if abs(eta_isentropic - 1.0) < 1e-9:
            T2 = T2s
        else:
            def dh(T):
                return fluid.h_vap(T, P_cond) - h2
            T2 = brentq(dh, T_cond + 0.01, T_cond + 400.0, xtol=1e-8) \
                if dh(T_cond + 0.01)*dh(T_cond + 400.0) < 0 else T2s

        # 3) condenser outlet: saturated (or slightly subcooled) liquid
        T3 = T_cond - subcool
        h3 = fluid.h_liq(T3, P_cond)
        s3 = fluid.s_liq(T3, P_cond)

        # 4) expansion valve: isenthalpic to P_evap. Two-phase at the outlet;
        #    quality from the saturated liquid/vapour enthalpies at T_evap.
        h4 = h3
        hf4 = fluid.h_liq(T_evap, P_evap)
        hg4 = fluid.h_vap(T_evap, P_evap)
        x4 = (h4 - hf4)/(hg4 - hf4)
        sf4 = fluid.s_liq(T_evap, P_evap)
        sg4 = fluid.s_vap(T_evap, P_evap)
        s4 = sf4 + x4*(sg4 - sf4)

        q_evap = h1 - h4
        w_comp = h2 - h1
        if w_comp <= 0:
            return {"error": f"non-physical compressor work ({w_comp:.3e} J/mol)"}

        return {
            "error": None,
            "COP": q_evap/w_comp,
            "q_evap": q_evap, "w_comp": w_comp,
            "P_evap": P_evap, "P_cond": P_cond,
            "pressure_ratio": P_cond/P_evap,
            "points": {
                "evap_out":  {"T": T1,     "P": P_evap, "h": h1,  "s": s1,  "x": 1.0},
                "comp_out":  {"T": T2,     "P": P_cond, "h": h2,  "s": s1,  "x": 1.0},
                "cond_out":  {"T": T3,     "P": P_cond, "h": h3,  "s": s3,  "x": 0.0},
                "valve_out": {"T": T_evap, "P": P_evap, "h": h4,  "s": s4,  "x": x4},
            },
        }
    except (ValueError, ZeroDivisionError, FloatingPointError) as e:
        return {"error": str(e)}


# =============================================================================
# Method definitions -- same numbers as phase_1_cubic_eos_validation_refstate_0903
# =============================================================================

NIST = dict(Tc=351.3, Pc=57.82e5, omega=0.2769,
            A=-6.098682, B=179.2200, C=-122.3682, D=32.30207, E=0.491361)
GCGP = dict(Tc=355.354, Pc=50.730e5, omega=0.1711,
            A=14.161, B=0.124, C=-6.340e-05, D=1.190e-8, E=0.0)


# -----------------------------------------------------------------------------
# 09/19/2026 delivery -- "spgp r32.xlsx", the six yellow-highlighted blocks
# under the red 09/19 date band. Added 09/24/2026.
# -----------------------------------------------------------------------------
# Two conversions are applied to every row below, because this file's Fluid
# takes neither of the forms the spreadsheet supplies:
#
#   (a) Shomate. The sheet gives a-e for the RAW convention
#           Cp = a + b*T + c*T^2 + d*T^3 + e/T^2        (T in K)
#       while Fluid.h_ig/s_ig integrate the Shomate form with t = T/1000:
#           A = a     B = b*1e3     C = c*1e6     D = d*1e9     E = e*1e-6
#       Verified to machine precision against the already-correct rows in
#       phase_1_cubic_eos_validation_refstate_GCN.py, and independently by
#       recomputing the sheet's own Cp column from its own coefficients.
#
#   (b) omega. The sheet has NO acentric-factor column -- it gives
#       "pvap (mmHg)", which is Psat at Tr = 0.7. Unlike
#       compare_cp_methods_0924.py, Fluid.__init__ takes omega directly and
#       immediately forms kappa from it, so there is no derivation step here.
#       omega_from_psat_mmHg() below closes that gap, and each row records
#       the resulting value in its comment so it can be eyeballed against the
#       compare_cp run without executing anything.
#
# Note Pc here is in Pa (Fluid's unit), where compare_cp_methods_0924.py uses
# bar. That transposition is the easiest mistake to make between the two files.

MMHG_TO_PA = 133.322


def omega_from_psat_mmHg(Pc_pa, psat_mmHg):
    """Pitzer acentric factor, omega = -1 - log10(Psat(Tr=0.7)/Pc).

    Same definition compare_cp_methods_0924.py's omega_from_psat() uses; the
    only difference is that this one takes the sheet's raw mmHg reading and
    Pc in Pa, matching Fluid's units.
    """
    return -1.0 - np.log10(psat_mmHg*MMHG_TO_PA/Pc_pa)


# --- block 1: "single output after thermoreconciliation" (sheet I13:L23) -----
# omega = 0.352253 -- physical, and the closest of the six to real R32 (~0.277).
# WARNING: Tc = 325.387 K. run_cycle's condenser sits at T_amb + 9 C, so at the
# 25 C ambient T_cond = 307.15 K and Tr = 0.944 -- near-critical, where the PR
# liquid root is poorly conditioned. Expect a degraded or non-computable result
# at that ambient. That is the fit's physics, not a solver failure.
FP_SINGLE_0919 = dict(
    Tc=325.387, Pc=52.168e5,
    omega=omega_from_psat_mmHg(52.168e5, 1738.797),
    A=32.6437, B=26.9797, C=591.354, D=-597.447, E=-0.0474595)

# omega = 0.433440 -- physical. Cp @298.15 K = 172.4 J/mol/K, roughly 3x real
# R32 (~55): the compressor work this produces will be correspondingly wrong.
GCN_SINGLE_0919 = dict(
    Tc=342.704, Pc=63.686e5,
    omega=omega_from_psat_mmHg(63.686e5, 1760.77),
    A=115.569, B=154.661, C=67.2036, D=174.449, E=0.00882475)

# --- block 2: "agentic approach" (sheet N13:Q21), labelled "spg" and "gcn" ---
# omega = 0.139628. Tc = 344.129 K -> Tr = 0.893 at the hottest condenser,
# comfortably sub-critical across all four ambients.
SPG_AGENTIC_0919 = dict(
    Tc=344.12921, Pc=52.183553e5,
    omega=omega_from_psat_mmHg(52.183553e5, 2837.94505436),
    A=32.849609, B=26.23072, C=590.20445, D=-592.38319, E=-0.048854552)

# omega = 0.448579 -- physical. Same ~3x Cp problem as GCN_SINGLE_0919; the two
# fits are near-identical in every coefficient.
GCN_AGENTIC_0919 = dict(
    Tc=342.89307, Pc=63.693509e5,
    omega=omega_from_psat_mmHg(63.693509e5, 1700.64828868),
    A=115.50434, B=155.06984, C=66.413405, D=174.76964, E=0.0092029555)

# --- block 3: "multi output after recon" (sheet I38:L48) --------------------
# omega = 0.051900 -- physical but low. Tc = 475.885 K is 125 K above real
# R32, so the whole cycle runs at unusually low reduced temperature.
FP_MULTI_0919 = dict(
    Tc=475.885, Pc=59.621e5,
    omega=omega_from_psat_mmHg(59.621e5, 3968.231),
    A=-6.96979, B=321.086, C=480.505, D=-1492.61, E=0.220907)

# omega = 0.396637 -- physical, unlike the 09/10 gcn_multi's unusable 1.254451.
# Cp @298.15 K = 239.5 J/mol/K, ~4x real R32 -- the worst Cp fit of the six.
# Exploratory only.
GCN_MULTI_0919 = dict(
    Tc=385.667, Pc=35.705e5,
    omega=omega_from_psat_mmHg(35.705e5, 1074.461),
    A=224.051, B=58.5856, C=-21.4895, D=-5.11394, E=0.000132741)


# ---- 09/25/2026 delivery ("spgp r32 8.xlsx"), red band I52 ----------------
# I53:L63 "single output after thermoreconciliation"        -> *_single_0925
# I65:L73 "single output after AGENTIC thermoreconciliation" -> *_agentic_0925
# All four discharge BELOW the condensing temperature (-1.4 to -7.2 K); see
# the superheat column this file prints. fp_* also have negative omega.

FP_SINGLE_0925 = dict(
    Tc=339.618, Pc=57.007e5,
    omega=omega_from_psat_mmHg(57.007e5, 4931.696),
    A=32.4646, B=11.4676, C=656.149, D=-665.778, E=-0.0544686)

GCN_SINGLE_0925 = dict(
    Tc=341.815, Pc=65.11e5,
    omega=omega_from_psat_mmHg(65.11e5, 3724.852),
    A=60.7245, B=108.396, C=-53.2265, D=356.844, E=0.0288859)

FP_AGENTIC_0925 = dict(
    Tc=345.56529, Pc=57.05028e5,
    omega=omega_from_psat_mmHg(57.05028e5, 5473.73681307),
    A=32.688237, B=10.68016, C=654.7169, D=-659.9117, E=-0.05605179)

GCN_AGENTIC_0925 = dict(
    Tc=336.81219, Pc=64.7458e5,
    omega=omega_from_psat_mmHg(64.7458e5, 3296.24475319),
    A=60.674957, B=108.7004, C=-53.91022, D=357.4103, E=0.02930463)


# Every fluid this file reports on: the two baselines, then the 09/19 six.
# The 09/10 "(7).xlsx" fits are deliberately absent -- that delivery is a
# different branch of the collaborators' work, not an earlier version of this
# one (none of its pvap values appear in the 09/19 sheet).
METHODS_0919 = [
    ("NIST",             NIST),
    ("GCGP",             GCGP),
    ("fp_single_0919",   FP_SINGLE_0919),
    ("fp_single_0925",   FP_SINGLE_0925),
    ("gcn_single_0925",  GCN_SINGLE_0925),
    ("fp_agentic_0925",  FP_AGENTIC_0925),
    ("gcn_agentic_0925", GCN_AGENTIC_0925),
]


def perturbed(tc_dev_pct, pc_dev_pct, base=NIST, name=None):
    """NIST baseline with Tc and/or Pc shifted by a percentage -- exactly the
    perturbation register_case() applies in the sweep. omega and the Shomate
    coefficients are held at baseline, same as the sweep."""
    p = dict(base)
    p["Tc"] = base["Tc"]*(1 + tc_dev_pct/100.0)
    p["Pc"] = base["Pc"]*(1 + pc_dev_pct/100.0)
    return Fluid(name or f"Tc{tc_dev_pct:+.2f}%_Pc{pc_dev_pct:+.2f}%", **p)


def show(res, label):
    print(f"\n  {label}")
    if res.get("error"):
        print(f"    NOT COMPUTABLE: {res['error']}")
        return
    print(f"    COP = {res['COP']:.4f}   "
          f"(q_evap = {res['q_evap']/MW/1e3:.2f} kJ/kg, "
          f"w_comp = {res['w_comp']/MW/1e3:.2f} kJ/kg, "
          f"P_cond/P_evap = {res['pressure_ratio']:.3f})")
    print(f"    {'point':<11}{'T [C]':>10}{'P [bar]':>11}"
          f"{'h [kJ/kg]':>12}{'s [kJ/kg/K]':>13}{'x':>8}")
    for k in ["evap_out", "comp_out", "cond_out", "valve_out"]:
        p = res["points"][k]
        print(f"    {k:<11}{p['T']-273.15:10.2f}{p['P']/1e5:11.4f}"
              f"{p['h']/MW/1e3:12.2f}{p['s']/MW/1e3:13.4f}{p['x']:8.4f}")


if __name__ == "__main__":
    AMBIENTS = [10, 15, 20, 25]

    # ---- validation 1: saturation pressures vs the Linde datasheet --------
    print("=" * 74)
    print("  VALIDATION 1 -- NIST saturation pressure vs Linde datasheet")
    print("=" * 74)
    nist = Fluid("NIST", **NIST)
    print(f"  {'T [C]':>8}{'this model [bar]':>20}{'Linde [bar]':>14}{'error':>10}")
    for T_C, linde in [(-30, 2.7344), (-10, 5.8263), (0, 8.131),
                       (10, 11.069), (20, 14.746), (30, 19.275)]:
        p = nist.psat(T_C + 273.15)/1e5
        print(f"  {T_C:>8}{p:>20.4f}{linde:>14.4f}{100*(p-linde)/linde:>9.2f}%")

    # ---- validation 2: COP vs the IDAES pipeline's own converged results --
    print("\n" + "=" * 74)
    print("  VALIDATION 2 -- COP vs the IDAES pipeline (phase6_final_0910.py run)")
    print("=" * 74)
    IDAES_COP = {"NIST": {10: 3.5677, 15: 3.3974, 20: 3.1316, 25: 2.8005},
                 "GCGP": {10: 3.5094, 15: 3.2581, 20: 3.0549, 25: 2.8910}}
    for label, params in [("NIST", NIST), ("GCGP", GCGP)]:
        f = Fluid(label, **params)
        print(f"\n  {label}:")
        print(f"    {'T_amb':>7}{'algebraic':>12}{'IDAES':>10}{'diff':>10}")
        for Tamb in AMBIENTS:
            r = run_cycle(f, T_amb_C=Tamb)
            if r.get("error"):
                print(f"    {Tamb:>7}{'--':>12}{IDAES_COP[label][Tamb]:>10.4f}   {r['error']}")
                continue
            ref = IDAES_COP[label][Tamb]
            print(f"    {Tamb:>7}{r['COP']:>12.4f}{ref:>10.4f}"
                  f"{100*(r['COP']-ref)/ref:>9.2f}%")

    # ---- Linde saturation check, every method ----------------------------
    # Validation 1 above is NIST-only and deliberately stays that way (it is
    # the file's self-check). This is the comparison across the full method
    # set: the same Linde datasheet points, every fluid, so the 09/19 fits can
    # be read against the same reference compare_cp_methods_0924.py uses.
    print("\n" + "=" * 100)
    print("  SATURATION PRESSURE vs LINDE DATASHEET -- all methods")
    print("=" * 100)
    LINDE_PSAT = [(-30, 2.7344), (-10, 5.8263), (0, 8.131),
                  (10, 11.069), (20, 14.746), (30, 19.275)]
    sat_rows = []
    hdr = "".join(f"{str(T) + ' C':>11}" for T, _ in LINDE_PSAT)
    print(f"  {'Psat [bar] at T =':<18}{hdr}{'MAPE':>10}")
    print(f"  {'Linde datasheet':<18}" + "".join(f"{lp:>11.4f}" for _, lp in LINDE_PSAT) + f"{'ref':>10}")
    print("  " + "-" * 96)
    for label, params in METHODS_0919:
        f = Fluid(label, **params)
        cells, errs = [], []
        for T_C, linde in LINDE_PSAT:
            try:
                pv = f.psat(T_C + 273.15)/1e5
                cells.append(f"{pv:>11.4f}")
                errs.append(abs(pv - linde)/linde)
            except ValueError:
                cells.append(f"{'--':>11}")
        mape = f"{100*np.mean(errs):>9.2f}%" if errs else f"{'--':>10}"
        print(f"  {label:<18}" + "".join(cells) + mape)
        rec = {"method": label, "Tc_K": f.Tc, "Pc_bar": f.Pc/1e5, "omega": f.omega,
               "Psat_MAPE_pct": 100*np.mean(errs) if errs else None}
        for T_C, linde in LINDE_PSAT:
            try:
                rec[f"Psat_{T_C}C_bar"] = f.psat(T_C + 273.15)/1e5
            except ValueError:
                rec[f"Psat_{T_C}C_bar"] = None
            rec[f"Linde_{T_C}C_bar"] = linde
        sat_rows.append(rec)

    # ---- COP sweep, every method, every ambient ---------------------------
    # run_cycle returns {'error': ...} rather than raising when a case is
    # physically impossible -- most often T_cond >= Tc, i.e. no saturation
    # state exists at the condenser. Those are reported as "--" with the
    # reason alongside; they are a real result about the fit, not a failure
    # of this script.
    print("\n" + "=" * 100)
    print("  COP vs AMBIENT -- all methods")
    print(f"  (T_evap = -29 C, condenser approach = +9 K, eta_isentropic = 0.9999,")
    print(f"   superheat = subcool = 0 -- run_cycle defaults)")
    print("=" * 100)
    print(f"  {'method':<18}{'Tc [K]':>9}{'omega':>10}" +
          "".join(f"{'COP@' + str(t) + 'C':>11}" for t in AMBIENTS) + "   note")
    print("  " + "-" * 96)
    cop_table = {}
    cop_rows, point_rows = [], []
    for label, params in METHODS_0919:
        f = Fluid(label, **params)
        cells, notes = [], []
        cop_table[label] = {}
        for Tamb in AMBIENTS:
            r = run_cycle(f, T_amb_C=Tamb)
            if r.get("error"):
                cells.append(f"{'--':>11}")
                notes.append(f"{Tamb}C: {r['error']}")
                cop_rows.append({"method": label, "T_amb_C": Tamb, "COP": None,
                                 "q_evap_kJ_per_kg": None, "w_comp_kJ_per_kg": None,
                                 "P_evap_bar": None, "P_cond_bar": None,
                                 "pressure_ratio": None, "error": r["error"]})
            else:
                cells.append(f"{r['COP']:>11.4f}")
                cop_table[label][Tamb] = r["COP"]
                cop_rows.append({"method": label, "T_amb_C": Tamb, "COP": r["COP"],
                                 "q_evap_kJ_per_kg": r["q_evap"]/MW/1e3,
                                 "w_comp_kJ_per_kg": r["w_comp"]/MW/1e3,
                                 "P_evap_bar": r["P_evap"]/1e5,
                                 "P_cond_bar": r["P_cond"]/1e5,
                                 "pressure_ratio": r["pressure_ratio"], "error": None})
                for pt, v in r["points"].items():
                    point_rows.append({"method": label, "T_amb_C": Tamb, "point": pt,
                                       "T_C": v["T"]-273.15, "P_bar": v["P"]/1e5,
                                       "h_kJ_per_kg": v["h"]/MW/1e3,
                                       "s_kJ_per_kgK": v["s"]/MW/1e3, "quality_x": v["x"]})
        note = "; ".join(notes)
        print(f"  {label:<18}{f.Tc:>9.3f}{f.omega:>10.6f}" + "".join(cells) +
              (f"   {note}" if note else ""))

    # ---- monotonicity check ----------------------------------------------
    # COP must fall as ambient rises (hotter condenser -> higher pressure
    # ratio -> more compressor work for the same duty). A method that
    # violates this is not physical, whatever its COP magnitude. This is the
    # same property phase6_final.py chases with warmstarts; here it is just
    # read off, since the calculation is deterministic.
    print("\n" + "=" * 100)
    print("  MONOTONICITY -- COP must decrease as ambient increases")
    print("=" * 100)
    mono_rows = []
    for label, _ in METHODS_0919:
        vals = [cop_table[label].get(t) for t in AMBIENTS]
        got = [v for v in vals if v is not None]
        if len(got) < 2:
            print(f"  {label:<18} insufficient converged points ({len(got)}/4)")
            continue
        ok = all(b < a for a, b in zip(got, got[1:]))
        span = f"{got[0]:.4f} -> {got[-1]:.4f}"
        print(f"  {label:<18}{'MONOTONIC' if ok else 'NOT MONOTONIC':>16}   "
              f"{span}   ({len(got)}/4 points)")
        mono_rows.append({"method": label, "monotonic": bool(ok),
                          "points_converged": len(got), "points_total": len(AMBIENTS),
                          "COP_first": got[0], "COP_last": got[-1]})

    # ---- Excel export (added 09/24/2026) ---------------------------------
    # Same multi-sheet pd.ExcelWriter convention as optimizer_noise_free_
    # sensitivity_0921.py. Everything printed above, plus the four state
    # points behind every COP so the cycle can be inspected or re-plotted
    # without re-running.
    OUT_XLSX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "cycle_algebraic_1001.xlsx")
    method_rows = [{"method": n, "Tc_K": p["Tc"], "Pc_bar": p["Pc"]/1e5,
                    "omega": p["omega"], "A": p["A"], "B": p["B"], "C": p["C"],
                    "D": p["D"], "E": p["E"]} for n, p in METHODS_0919]
    with pd.ExcelWriter(OUT_XLSX) as w:
        pd.DataFrame(method_rows).to_excel(w, sheet_name="methods", index=False)
        pd.DataFrame(cop_rows).to_excel(w, sheet_name="COP_sweep", index=False)
        pd.DataFrame(point_rows).to_excel(w, sheet_name="state_points", index=False)
        pd.DataFrame(sat_rows).to_excel(w, sheet_name="saturation_vs_Linde", index=False)
        pd.DataFrame(mono_rows).to_excel(w, sheet_name="monotonicity", index=False)
    print(f"\n  wrote {OUT_XLSX}")

    # ---- the case that fails in the IDAES sweep ---------------------------
    print("\n" + "=" * 74)
    print("  THE FAILING SWEEP POINT -- Tc +3.00%, Pc +0.00%, T_amb = 20 C")
    print("=" * 74)
    show(run_cycle(perturbed(3.0, 0.0), T_amb_C=20), "Tc +3.00%, Pc +0.00%, T_amb = 20 C")
    show(run_cycle(perturbed(0.0, 0.0), T_amb_C=20), "baseline  Tc +0.00%, Pc +0.00%, T_amb = 20 C")

    # ---- the whole deviation range, at one ambient ------------------------
    print("\n" + "=" * 74)
    print("  FULL Tc SWEEP AT T_amb = 20 C (Pc at baseline)")
    print("=" * 74)
    print(f"  {'Tc dev %':>10}{'Tc [C]':>10}{'COP':>10}{'P_evap':>10}"
          f"{'P_cond':>10}{'ratio':>9}   note")
    DEV_PCTS = [-10,  10]
    for d in DEV_PCTS:
        f = perturbed(d, 0.0)
        r = run_cycle(f, T_amb_C=20)
        if r.get("error"):
            print(f"  {d:>10.2f}{f.Tc-273.15:>10.2f}{'--':>10}{'--':>10}"
                  f"{'--':>10}{'--':>9}   {r['error']}")
        else:
            print(f"  {d:>10.2f}{f.Tc-273.15:>10.2f}{r['COP']:>10.4f}"
                  f"{r['P_evap']/1e5:>10.4f}{r['P_cond']/1e5:>10.4f}"
                  f"{r['pressure_ratio']:>9.3f}")
    print()
# -5, -3, -1, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, 1, 3, 5,