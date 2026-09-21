"""
optimizer_noise_free_sensitivity_0921.py -- the full Tc/Pc sensitivity
sweep, computed ALGEBRAICALLY with no IDAES, no Pyomo and no IPOPT, and
scored on four targets instead of COP alone: COP error, T-S curve error,
P-H curve error, and the area under both curves.

Why no optimizer (2026-09-17 finding, still the reason this file exists):
the original IDAES Tc/Pc sensitivity sweep solves every unit balance,
property correlation and phase-equilibrium relation simultaneously as one
large nonlinear system via IPOPT. That injects solver behaviour into the
answer: convergence flags that depend on the starting point, a loosened
tolerance (tol=1e-4) that sets a floor on resolvable COP differences, and
251 of 900 case-ambient solves reported non-converged. With the evaporator
and condenser saturation temperatures fixed (never optimized -- always the
NIST-derived operating point, for every case) and superheat = subcool = 0,
the cycle is instead a sequential, deterministic calculation:

    P_evap  = Psat(T_evap)
    state 1 = saturated vapour at T_evap          -> h1, s1
    P_cond  = Psat(T_cond)
    state 2 = compress isentropically to P_cond   -> h2s, then h2 via eta
    state 3 = saturated liquid at T_cond          -> h3
    state 4 = isenthalpic expansion to P_evap     -> h4 = h3
    COP     = (h1 - h4) / (h2 - h1)

Only step 2 is implicit, and only in one variable: find T2 such that
s_vap(T2, P_cond) = s1 -- a 1-D root find on a monotone function. One
answer per (Tc, Pc), every time; no convergence flag, no warmstart, no
solver-noise floor.

Why four targets, not just COP (2026-09-21, per the Colon group's actual
meeting notes): an acceptable COP error can arise from unacceptable T-S
and P-H diagrams -- i.e. a property-prediction method can get the cycle's
energy balance roughly right while getting the underlying thermodynamic
behaviour wrong. So alongside COP error, every case also reports how far
its saturated-liquid/vapour T-S and P-H curves are from the NIST baseline
curve, in two complementary senses: TS_mae/PH_mae (pointwise mean absolute
error, catches a uniformly shifted curve) and TS_area_diff/PH_area_diff
(absolute difference in area under each branch, catches a curve that
wiggles around the baseline without a large mean error). These curve
metrics are computed once per (Tc, Pc) case over the full temperature
range (not windowed per ambient) since they are a property of the fluid
itself, not of any one cycle operating point -- decided explicitly rather
than assumed.

This file does NOT compare anything against IDAES and has no validation
section -- that comparison lives in cycle_algebraic_0921.py, which this
file's Fluid/run_cycle/curve tools are a direct copy of. This file's only
job is to run the sweep and write results.

Physics: Peng-Robinson cubic EoS, Ambrose-Walton saturation pressure
(Poling sec. 7-4), PR departure functions for h and s (Poling Table 6-3),
Shomate ideal-gas Cp integrated for the ideal-gas h and s (t = T/1000).
Reference state does not matter for COP (a ratio of enthalpy DIFFERENCES,
so any constant offset cancels) -- but DOES matter for TS_mae/PH_mae/
area_diff, which is why those compare each fluid against the same NIST
baseline rather than against an externally-anchored datasheet value.

Author: Shilpa Narasimhan. Support: Claude AI.

Date created: 09/21/2026
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import brentq

R = 8.314          # J/mol/K  -- matches phase_1_cubic_eos_validation_refstate_0903.py
MW = 52.024e-3     # kg/mol
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
    shared temperature grid (the 20 Linde datasheet points, fixed in
    absolute temperature -- the same physical points for every fluid) --
    the raw data behind a T-S diagram.

    CHANGED 2026-09-21, per Shilpa's explicit instruction: every one of
    the 20 grid points now always gets a slot in the returned arrays, for
    every fluid -- nothing is omitted/shrunk. Where a temperature sits at
    or above fluid.Tc, the slot holds NaN instead of a real value, since
    above the critical temperature liquid and vapour are the same phase
    and "saturated vapour entropy" is not a defined quantity there for
    ANY fluid -- that part is physics, not something a different
    implementation could avoid. What changed is only how the "no value
    here" case is represented: an explicit NaN at a fixed position,
    instead of the point disappearing from the array entirely.

    Returns (T_K, sf, sg), each length len(T_grid_C): temperature [K]
    (always fully populated), saturated-liquid entropy [J/mol/K],
    saturated-vapour entropy [J/mol/K] (NaN where undefined).
    """
    T_K = np.array([T_C + 273.15 for T_C in T_grid_C])
    sf = np.full(len(T_grid_C), np.nan)
    sg = np.full(len(T_grid_C), np.nan)
    for i, T in enumerate(T_K):
        if T >= fluid.Tc:
            continue
        _, _, s_liq = fluid.sat_liquid(T)
        _, _, s_vap = fluid.sat_vapour(T)
        sf[i] = s_liq
        sg[i] = s_vap
    return T_K, sf, sg


def ph_curve(fluid, T_grid_C=TS_GRID_C):
    """Saturation pressure and saturated-liquid/vapour enthalpy of `fluid`
    over the shared temperature grid -- the raw data behind a P-H diagram.
    Same always-20-slots-with-NaN behaviour as ts_curve() above.

    Returns (T_K, P, hf, hg), each length len(T_grid_C): temperature [K]
    (always fully populated), saturation pressure [Pa], saturated-liquid
    enthalpy [J/mol], saturated-vapour enthalpy [J/mol] (NaN where
    undefined).
    """
    T_K = np.array([T_C + 273.15 for T_C in T_grid_C])
    P_arr = np.full(len(T_grid_C), np.nan)
    hf = np.full(len(T_grid_C), np.nan)
    hg = np.full(len(T_grid_C), np.nan)
    for i, T in enumerate(T_K):
        if T >= fluid.Tc:
            continue
        P, h_liq, _ = fluid.sat_liquid(T)
        _, h_vap, _ = fluid.sat_vapour(T)
        P_arr[i] = P
        hf[i] = h_liq
        hg[i] = h_vap
    return T_K, P_arr, hf, hg


def curve_area(x, y):
    """Trapezoidal area under y(x), using only the points where both x
    and y are finite (NaN entries -- points above a fluid's Tc -- are
    dropped for this calculation only; the arrays themselves stay full
    length everywhere else). Returns NaN if fewer than 2 finite points
    remain, rather than raising, so a single bad case does not stop a
    sweep."""
    x, y = np.asarray(x), np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
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
        pooled), ignoring any grid point where either fluid lacks a value
        (np.nanmean). Good at catching a curve that is uniformly shifted
        off the baseline.
      - *_area_diff: absolute difference in the area under each branch,
        computed over whichever points are finite for that fluid
        (curve_area already drops NaNs). Good at catching a curve that
        wiggles above and below the baseline -- pointwise errors can
        partly cancel in the mean, but the enclosed area still shows the
        shape is wrong. This is the metric the Colon group's "COP error
        alone can hide bad thermodynamic behaviour" concern was
        specifically about.

    ts_curve/ph_curve now always return the full 20-point grid for both
    fluid and baseline (see ts_curve docstring, 2026-09-21) -- no
    truncation is needed or performed; every point is compared at its own
    fixed grid position, NaN-aware.

    Returns a dict: TS_mae [J/mol/K], PH_mae [J/mol], TS_area_diff,
    PH_area_diff, and n_points (how many of the 20 grid points `fluid`
    itself had a real value at -- 0 if none).
    """
    T, sf_f, sg_f = ts_curve(fluid, T_grid_C)
    _, sf_b, sg_b = ts_curve(baseline, T_grid_C)
    _, _, hf_f, hg_f = ph_curve(fluid, T_grid_C)
    _, _, hf_b, hg_b = ph_curve(baseline, T_grid_C)

    n_points = int(np.isfinite(sf_f).sum())
    if n_points == 0:
        return dict(TS_mae=np.nan, PH_mae=np.nan,
                     TS_area_diff=np.nan, PH_area_diff=np.nan, n_points=0)

    with np.errstate(invalid="ignore"):
        TS_mae = float(np.nanmean(np.abs(np.concatenate([sf_f - sf_b, sg_f - sg_b]))))
        PH_mae = float(np.nanmean(np.abs(np.concatenate([hf_f - hf_b, hg_f - hg_b]))))

    TS_area_diff = (abs(curve_area(T, sg_f) - curve_area(T, sg_b))
                     + abs(curve_area(T, sf_f) - curve_area(T, sf_b)))
    PH_area_diff = (abs(curve_area(T, hg_f) - curve_area(T, hg_b))
                     + abs(curve_area(T, hf_f) - curve_area(T, hf_b)))

    return dict(TS_mae=TS_mae, PH_mae=PH_mae,
                 TS_area_diff=TS_area_diff, PH_area_diff=PH_area_diff,
                 n_points=n_points)








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


# =============================================================================
# Sweep
# =============================================================================

DEV_PCTS = [-10, -5, -3, -1, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, 1, 3, 5, 10]
AMBIENTS = [10, 15, 20, 25]

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_XLSX = os.path.join(HERE, "optimizer_noise_free_sensitivity_0921.xlsx")

NIST_BASELINE = Fluid("NIST", **NIST)


def rows_for_case(sweep, tc_dev, pc_dev):
    """One row per ambient. The curve-shape metrics (TS_mae, PH_mae,
    TS_area_diff, PH_area_diff) depend only on the parameter set, not on
    T_amb, so they are computed once per case and repeated on every row --
    tidy long-format output, no separate join needed downstream.
    """
    f = perturbed(tc_dev, pc_dev)
    shape = curve_shape_metrics(f, NIST_BASELINE)

    rows = []
    for Tamb in AMBIENTS:
        r = run_cycle(f, T_amb_C=Tamb)
        cop = None if r.get("error") else r["COP"]
        rows.append({
            "sweep": sweep,
            "Tc_dev_pct": tc_dev, "Pc_dev_pct": pc_dev,
            "Tc_K": f.Tc, "Tc_C": f.Tc - 273.15, "Pc_Pa": f.Pc,
            "omega": f.omega,
            "T_amb_C": Tamb,
            "COP": cop,
            "TS_mae": shape["TS_mae"], "PH_mae": shape["PH_mae"],
            "TS_area_diff": shape["TS_area_diff"],
            "PH_area_diff": shape["PH_area_diff"],
            "shape_n_points": shape["n_points"],
        })
    return rows


if __name__ == "__main__":
    tc_only, pc_only, grid = [], [], []

    print("running Tc_only ...")
    for d in DEV_PCTS:
        tc_only += rows_for_case("Tc_only", d, 0)

    print("running Pc_only ...")
    for d in DEV_PCTS:
        pc_only += rows_for_case("Pc_only", 0, d)

    print(f"running Tc_Pc_grid ({len(DEV_PCTS)}x{len(DEV_PCTS)}) ...")
    for dtc in DEV_PCTS:
        for dpc in DEV_PCTS:
            grid += rows_for_case("Tc_Pc_grid", dtc, dpc)

    dfs = {"Tc_only": pd.DataFrame(tc_only),
           "Pc_only": pd.DataFrame(pc_only),
           "Tc_Pc_grid": pd.DataFrame(grid)}

    with pd.ExcelWriter(OUT_XLSX) as w:
        for name, df in dfs.items():
            df.to_excel(w, sheet_name=name, index=False)

    print(f"\nwrote {OUT_XLSX}")
    for name, df in dfs.items():
        ok = df["COP"].notna().sum()
        print(f"  {name:<12} {len(df):>5} rows   COP computed: {ok}/{len(df)} "
              f"({100.0*ok/len(df):.1f}%)")
    print()