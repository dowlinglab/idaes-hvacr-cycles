"""
verify_akasaka_lemmon_SI.py -- 10/08/2026.

QUESTION: is CoolProp's R1234ze(Z) really the published Akasaka & Lemmon
equation of state, and is the ideal-gas Cp used by
compare_cp_methods_R1234zeZ_1008.py the published one?

The Supporting Information of

    R. Akasaka and E.W. Lemmon, "Fundamental Equations of State for
    cis-1,3,3,3-Tetrafluoropropene [R-1234ze(Z)] and 3,3,3-Trifluoropropene
    (R-1243zf)", J. Chem. Eng. Data (2019). DOI 10.1021/acs.jced.9b00007

carries a reference Python implementation written by Ian Bell (NIST), together
with the coefficients and a set of verification state points. That code is
transcribed below EXACTLY as published -- the EOS class, the coefficient
arrays, and the six (T, rho) test points -- and then:

  STAGE 1  the published implementation is evaluated at its own test points
           and compared against CoolProp. If CoolProp is the same equation,
           p, cv, cp and w must agree to round-off.
  STAGE 2  the ideal-gas Cp used for the Shomate refit in
           compare_cp_methods_R1234zeZ_1008.py is compared against the
           ideal-gas part of the published equation. This checks the
           reading of the alpha0 terms, in particular that the two Einstein
           temperatures are 20 K and 1335 K and that cp0/R = 4 + (PE terms).
  STAGE 3  the Shomate A-E actually compiled into the comparison script are
           evaluated against that same published Cp0.

Run from the R32 directory in an environment with CoolProp.
"""

import os
import sys
from math import exp, log

import numpy as np

# so the comparison script next to this file is importable however it is run
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PASS = True


def fail(msg):
    global PASS
    PASS = False
    print(f"   ** FAIL: {msg}")


def banner(t):
    print(f"\n{'='*72}\n  {t}\n{'='*72}")


# =============================================================================
# Transcribed verbatim from the SI (Ian Bell, NIST)
# =============================================================================

class EOS(object):
    def __init__(self, *, n, t, d, l, eta, beta, gamma, epsilon, Tc, rhoc,
                 a1, a2, c0, u, v, Nterms=(5, 5, 6)):
        assert (sum([abs(len(v_) - len(n))
                     for v_ in (n, t, d, l, eta, beta, gamma, epsilon)]) == 0)
        self.n = n
        self.t = t
        self.d = d
        self.l = l
        self.eta = eta
        self.beta = beta
        self.gamma = gamma
        self.epsilon = epsilon
        self.Tc = Tc
        self.rhoc = rhoc
        self.a1 = a1
        self.a2 = a2
        self.c0 = c0
        self.u = u
        self.v = v
        self.Npoly, self.Nexp, self.Ngauss = Nterms

    def resid_derivs(self, T, rho):
        tau = self.Tc / T
        delta = rho / self.rhoc
        A00, A10, A01, A20, A11, A02 = 0, 0, 0, 0, 0, 0

        for i in range(0, self.Npoly):
            _s = self.n[i] * delta**self.d[i] * tau**self.t[i]
            A00 += _s
            A10 += self.t[i] * _s
            A01 += self.d[i] * _s
            A20 += self.t[i] * (self.t[i] - 1) * _s
            A11 += self.d[i] * self.t[i] * _s
            A02 += self.d[i] * (self.d[i] - 1) * _s

        for i in range(self.Npoly, self.Npoly + self.Nexp):
            _s = (self.n[i] * delta**self.d[i] * tau**self.t[i]
                  * exp(-delta**self.l[i]))
            A00 += _s
            A10 += _s * self.t[i]
            A01 += _s * (self.d[i] - self.l[i] * delta**self.l[i])
            A20 += _s * self.t[i] * (self.t[i] - 1)
            A11 += _s * self.t[i] * (self.d[i] - self.l[i] * delta**self.l[i])
            A02 += _s * ((self.d[i] - self.l[i] * delta**self.l[i])
                         * (self.d[i] - 1 - self.l[i] * delta**self.l[i])
                         - self.l[i]**2 * delta**self.l[i])

        for i in range(self.Npoly + self.Nexp,
                       self.Npoly + self.Nexp + self.Ngauss):
            _s = (self.n[i] * delta**self.d[i] * tau**self.t[i]
                  * exp(-self.eta[i] * (delta - self.epsilon[i])**2
                        - self.beta[i] * (tau - self.gamma[i])**2))
            A00 += _s
            A10 += _s * (self.t[i] - 2 * self.beta[i] * tau * (tau - self.gamma[i]))
            A01 += _s * (self.d[i] - 2 * self.eta[i] * delta * (delta - self.epsilon[i]))
            A20 += _s * ((self.t[i] - 2 * self.beta[i] * tau * (tau - self.gamma[i]))**2
                         - self.t[i] - 2 * self.beta[i] * tau**2)
            A11 += _s * ((self.t[i] - 2 * self.beta[i] * tau * (tau - self.gamma[i]))
                         * (self.d[i] - 2 * self.eta[i] * delta * (delta - self.epsilon[i])))
            A02 += _s * ((self.d[i] - 2 * self.eta[i] * delta * (delta - self.epsilon[i]))**2
                         - self.d[i] - 2 * self.eta[i] * delta**2)

        return A00, A10, A01, A20, A11, A02

    def ideal_derivs(self, T, rho):
        tau = self.Tc / T
        delta = rho / self.rhoc

        A01 = 1
        A02 = -1
        A11 = 0
        A00 = ((0 if delta == 0 else log(delta)) + self.a1 + self.a2 * tau
               + (self.c0 - 1) * log(tau)
               + sum([v_i * log(1 - exp(-u_i * tau / self.Tc))
                      for u_i, v_i in zip(self.u, self.v)]))
        A10 = (self.a2 * tau + (self.c0 - 1)
               + tau * sum([v_i * (u_i / self.Tc) / (exp(u_i * tau / self.Tc) - 1)
                            for u_i, v_i in zip(self.u, self.v)]))
        A20 = (-(self.c0 - 1)
               - tau**2 * sum([v_i * (u_i / self.Tc)**2 * exp(u_i * tau / self.Tc)
                               / (exp(u_i * tau / self.Tc) - 1)**2
                               for u_i, v_i in zip(self.u, self.v)]))
        return A00, A10, A01, A20, A11, A02

    def derivs(self, T, rho):
        return [resid + ideal for resid, ideal
                in zip(self.resid_derivs(T, rho), self.ideal_derivs(T, rho))]


# ---- R1234ze(Z) coefficients, verbatim from the SI --------------------------
_n = [0.03194509, 1.394592, -2.300799, -0.2556693, 0.1282934, -1.335381,
      -1.366494, 0.2004912, -0.6489709, -0.02220033, 1.66538, 0.3427048,
      -0.6510217, -0.5067066, -0.1231787, 0.08828106]
_t = [1.0, 0.333, 1.0, 1.0, 0.38, 2.85, 3.16, 0.607, 2.2, 1.0, 1.83, 3.3,
      1.9, 2.6, 2.9, 3.0]
_d = [4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 2, 3, 2]
_l = [0]*5 + [2, 2, 1, 2, 1] + [0]*6
_eta = [0]*10 + [1.108, 1.579, 1.098, 0.672, 3.38, 1.6]
_beta = [0]*10 + [0.563, 1.724, 0.806, 0.505, 26.4, 8.82]
_gamma = [0]*10 + [1.246, 1.05, 1.0, 0.677, 1.302, 1.274]
_epsilon = [0]*10 + [0.933, 0.786, 0.496, 0.327, 0.523, 0.308]

ZeZ = EOS(n=_n, t=_t, d=_d, l=_l, eta=_eta, beta=_beta, gamma=_gamma,
          epsilon=_epsilon, Tc=423.27, rhoc=4000,
          a1=-2.422442259, a2=8.190539844, c0=4,
          v=[4.2365, 13.063], u=[20, 1335], Nterms=(5, 5, 6))

R_SI = 8.3144598          # J/mol/K, the value the SI code uses
M_SI = 114.0416 / 1000.0  # kg/mol

# The SI's own verification state points, (T [K], rho [mol/L])
TEST_POINTS = [(300, 0), (300, 11), (300, 0.05), (400, 8), (400, 1), (424, 4)]


def si_properties(T, rho_molL):
    """p [Pa], cv [J/mol/K], cp [J/mol/K], w [m/s] from the published code."""
    rho = rho_molL * 1000.0
    Ar00, Ar10, Ar01, Ar20, Ar11, Ar02 = ZeZ.resid_derivs(T, rho)
    A00, A10, A01, A20, A11, A02 = ZeZ.derivs(T, rho)
    p = rho * R_SI * T * (Ar01 + 1)
    cv = -R_SI * A20
    cp = cv + R_SI * (1 + Ar01 - Ar11)**2 / (1 + 2*Ar01 + Ar02)
    w = (R_SI*T/M_SI * ((1 + 2*Ar01 + Ar02)
                        - (1 + Ar01 - Ar11)**2 / (A20)))**0.5
    return p, cv, cp, w


# =============================================================================
banner("STAGE 1 -- published SI implementation vs CoolProp")
# =============================================================================
try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    print("   CoolProp not importable -- cannot run STAGE 1.")
    sys.exit(2)

FLUID = "R1234ze(Z)"
print(f"   {'T [K]':>7}{'rho[mol/L]':>12}"
      f"{'p SI [MPa]':>13}{'p CP [MPa]':>13}{'d %':>9}"
      f"{'cp SI':>10}{'cp CP':>10}{'d %':>9}")
for T, rho_molL in TEST_POINTS:
    p_si, cv_si, cp_si, w_si = si_properties(T, rho_molL)
    if rho_molL == 0:
        # ideal-gas limit: CoolProp cannot take rho = 0; compare cp0 instead
        cp_cp = PropsSI("Cp0molar", "T", T, "D", 1e-10, FLUID)
        p_cp = 0.0
        dp = 0.0
    else:
        p_cp = PropsSI("P", "T", T, "Dmolar", rho_molL*1000.0, FLUID)
        cp_cp = PropsSI("Cpmolar", "T", T, "Dmolar", rho_molL*1000.0, FLUID)
        dp = 100.0*(p_cp - p_si)/p_si if p_si != 0 else 0.0
    dcp = 100.0*(cp_cp - cp_si)/cp_si
    print(f"   {T:>7}{rho_molL:>12}"
          f"{p_si/1e6:>13.6f}{p_cp/1e6:>13.6f}{dp:>9.4f}"
          f"{cp_si:>10.4f}{cp_cp:>10.4f}{dcp:>9.4f}")
    if abs(dp) > 1e-3 or abs(dcp) > 1e-3:
        fail(f"SI vs CoolProp mismatch at T={T}, rho={rho_molL}")
print("\n   Agreement to <0.001 % means CoolProp implements this exact equation.")


# =============================================================================
banner("STAGE 2 -- ideal-gas Cp: this project's reading vs the published SI")
# =============================================================================
# The comparison script reads alpha0 out of CoolProp's JSON and evaluates
#     cp0/R = 4 + sum_i n_i x_i^2 e^{x_i} / (e^{x_i} - 1)^2 ,  x_i = theta_i/T
# with theta = [20, 1335] K.  The SI gives c0 = 4, u = [20, 1335],
# v = [4.2365, 13.063] -- so the two should be identical, not merely close.
from compare_cp_methods_R1234zeZ_1008 import (     # noqa: E402
    cp0_eos, CP0_THETA, CP0_N, PR_ROW, _MAXREL)

print(f"   Einstein temperatures used here : {np.round(CP0_THETA, 4)}")
print(f"   Einstein temperatures in the SI : {ZeZ.u}")
print(f"   n coefficients used here        : {CP0_N}")
print(f"   n coefficients in the SI        : {ZeZ.v}")
if not np.allclose(CP0_THETA, ZeZ.u, rtol=1e-6):
    fail("Einstein temperatures do not match the SI")
if not np.allclose(CP0_N, ZeZ.v, rtol=1e-12):
    fail("Planck-Einstein n coefficients do not match the SI")


def cp0_from_si(T):
    """cp0 from the published ideal_derivs: cp0 = cv0 + R, cv0 = -R*A20."""
    _, _, _, A20, _, _ = ZeZ.ideal_derivs(T, 0.0)
    return R_SI * (-A20 + 1.0)


print(f"\n   {'T [K]':>8}{'cp0 SI':>12}{'cp0 here':>12}{'d %':>10}")
for T in [238.13, 273.15, 300.0, 350.0, 400.0, 423.27]:
    a, b = cp0_from_si(T), float(cp0_eos(T)[0])
    d = 100.0*(b - a)/a
    print(f"   {T:>8.2f}{a:>12.5f}{b:>12.5f}{d:>10.6f}")
    if abs(d) > 1e-6:
        fail(f"cp0 mismatch at T={T}")


# =============================================================================
banner("STAGE 3 -- the Shomate A-E actually compiled into the script")
# =============================================================================
A, B, C, D, E = (PR_ROW["A"], PR_ROW["B"], PR_ROW["C"], PR_ROW["D"], PR_ROW["E"])
print(f"   A = {A:.6f}   B = {B:.6f}   C = {C:.6f}")
print(f"   D = {D:.6f}   E = {E:.6f}")
# The fit window is the experimental range, 238.13-372.61 K. Points beyond it
# are extrapolation and are shown but not asserted on -- a Shomate polynomial
# fitted over one interval is not expected to hold outside it, and the cycle
# never leaves the interval anyway.
FIT_LO, FIT_HI = 238.13, 372.61

print(f"\n   {'T [K]':>8}{'cp0 SI':>12}{'Shomate':>12}{'d %':>10}  note")
worst = 0.0
for T in [238.13, 273.15, 300.0, 350.0, 372.61, 400.0]:
    t = T/1000.0
    sh = A + B*t + C*t**2 + D*t**3 + E/t**2
    a = cp0_from_si(T)
    d = 100.0*(sh - a)/a
    inside = FIT_LO <= T <= FIT_HI
    if inside:
        worst = max(worst, abs(d))
    print(f"   {T:>8.2f}{a:>12.5f}{sh:>12.5f}{d:>10.5f}  "
          f"{'' if inside else 'EXTRAPOLATION, not asserted'}")
print(f"\n   worst deviation inside the fit window: {worst:.5f} %")
print(f"   (script reports {_MAXREL:.5f} % over its full fit window)")
if worst > 0.01:
    fail(f"Shomate refit deviates {worst:.4f} % from the published Cp0")


# =============================================================================
banner("VERDICT")
# =============================================================================
if PASS:
    print("   ALL CHECKS PASSED.")
    print("   CoolProp is the published Akasaka & Lemmon equation; the")
    print("   ideal-gas Cp used for the Shomate refit is the published one;")
    print("   and the fitted A-E reproduce it to well under 0.01 %.")
    sys.exit(0)
else:
    print("   ONE OR MORE CHECKS FAILED -- see ** FAIL lines above.")
    sys.exit(1)
