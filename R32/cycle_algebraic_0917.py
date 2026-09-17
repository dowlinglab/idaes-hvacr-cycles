"""
cycle_algebraic_0917.py -- the ideal vapor-compression cycle computed
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

Author: Shilpa Narasimhan. Support: Claude AI.
"""

import numpy as np
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


if __name__ == "__main__":
    AMBIENTS = [10, 15, 20, 25]

    # ---- validation 1: saturation pressures vs the Linde datasheet --------
    print("=" * 74)
    print("  VALIDATION 1 -- NIST saturation pressure vs Linde datasheet")
    print("=" * 74)
    nist = Fluid("NIST", **NIST)
    print(f"  {'T [C]':>8}{'this model [bar]':>20}{'Linde [bar]':>14}{'error':>10}")
    for T_C, linde in [(-30, 2.7344), (-10, 5.8263), (0, 8.131),
                       (10, 10.065), (20, 14.746), (30, 19.275)]:
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

    # ---- the case that fails in the IDAES sweep ---------------------------
    print("\n" + "=" * 74)
    print("  THE FAILING SWEEP POINT -- Tc +3.00%, Pc +0.00%, T_amb = 20 C")
    print("  (IDAES: not converged at any ambient, all 10 warmstart attempts)")
    print("=" * 74)
    show(run_cycle(perturbed(3.0, 0.0), T_amb_C=20), "Tc +3.00%, Pc +0.00%, T_amb = 20 C")
    show(run_cycle(perturbed(0.0, 0.0), T_amb_C=20), "baseline  Tc +0.00%, Pc +0.00%, T_amb = 20 C")

    # ---- the whole deviation range, at one ambient ------------------------
    print("\n" + "=" * 74)
    print("  FULL Tc SWEEP AT T_amb = 20 C (Pc at baseline)")
    print("=" * 74)
    print(f"  {'Tc dev %':>10}{'Tc [C]':>10}{'COP':>10}{'P_evap':>10}"
          f"{'P_cond':>10}{'ratio':>9}   note")
    DEV_PCTS = [-10, -5, -3, -1, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, 1, 3, 5, 10]
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
