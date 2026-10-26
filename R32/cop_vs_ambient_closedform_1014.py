"""
cop_vs_ambient_closedform_1014.py -- R-32 ambient sweep from the CLOSED-FORM
cycle. No IDAES, no solver.

WHY CLOSED FORM
    On 2026-10-08 the IDAES flowsheet was found to produce converged,
    plausible COPs from a loop that does not conserve energy: the arcs
    transmit flow, composition, pressure and temperature, but NOT phase
    fraction, so the stream leaving one unit is not the stream entering the
    next. Breaks of 12-67 % of the evaporator duty, in single-phase regions
    as well as two-phase. See BREADCRUMB_09-28.md, 10/08 evening section.

    Until the arcs are fixed, the only trustworthy cycle numbers come from
    the closed form: four state points, each a direct property lookup,
    nothing to solve. These values pass every check in cycle_guards by
    construction.

CYCLE
    T_evap = -29 C,  T_cond = T_amb + 9 C
    superheat = 3 K, subcool = 8 K,  compressor isentropic efficiency 0.9999
    saturated states at 0/0 are NOT used: set_specifications deactivates its
    own superheat/subcool constraints when passed 0, so the model never ran
    the cycle that was specified. 3 K / 8 K activates them.

SERIES
    Helmholtz    CoolProp, Tillner-Roth & Yokozeki
    PR / NIST    Peng-Robinson, NIST critical properties and Shomate Cp
    Carnot       T_evap / (T_cond - T_evap), the ceiling

OMITTED, deliberately: PR/GCGP and PR/SPGP. Their closed-form discharge
temperatures are 306-441 C and 36-55 C respectively, against NIST's physical
75-103 C for the same compression. That points at the Cp convention (the
GCGP-family coefficients are fitted at T = t/1000, not Shomate's t) and is
unresolved. Their COP values are not reportable until it is.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/14/2026
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as CP

import compare_cp_methods_0929_Helmholtz_SPGPlabel as M

FLUID = "R32"
T_EVAP = -29.0
APPROACH = 9.0
ETA = 0.9999
SUPERHEAT = 3.0
SUBCOOL = 8.0
AMBIENTS = np.array([10.0, 15.0, 20.0, 25.0])


def _solve_T(f, lo):
    """Expand upward until a sign change, then brentq. Avoids guessing a
    bracket, which fails for methods whose Tc is close to the condensing
    temperature."""
    a, fa = lo, f(lo)
    for hi in np.arange(lo + 2.0, lo + 800.0, 2.0):
        try:
            fb = f(hi)
        except Exception:
            continue
        if np.isfinite(fb) and fa * fb < 0:
            return brentq(f, a, hi)
        a, fa = hi, fb
    raise RuntimeError("no bracket found")


def cop_helmholtz(T_cond_C):
    Te, Tc = T_EVAP + 273.15, T_cond_C + 273.15
    pe = CP("P", "T", Te, "Q", 1, FLUID)
    pc = CP("P", "T", Tc, "Q", 0, FLUID)
    h1 = CP("H", "T", Te + SUPERHEAT, "P", pe, FLUID)       # 1 evap outlet
    s1 = CP("S", "T", Te + SUPERHEAT, "P", pe, FLUID)
    h2 = h1 + (CP("H", "P", pc, "S", s1, FLUID) - h1) / ETA  # 2 comp outlet
    h3 = CP("H", "T", Tc - SUBCOOL, "P", pc, FLUID)          # 3 cond outlet
    T2 = CP("T", "P", pc, "H", h2, FLUID) - 273.15
    return (h1 - h3) / (h2 - h1), T2                         # 4 valve: h4 = h3


def cop_pr(m, T_cond_C):
    Te, Tc = T_EVAP + 273.15, T_cond_C + 273.15
    pe, pc = m.saturation_pressure(Te), m.saturation_pressure(Tc)
    h1 = m._h_mass(Te + SUPERHEAT, pe, "vapor") + m.h_off
    s1 = m._s_mass(Te + SUPERHEAT, pe, "vapor") + m.s_off
    T2s = _solve_T(lambda T: m._s_mass(T, pc, "vapor") + m.s_off - s1, Tc + 0.01)
    h2 = h1 + ((m._h_mass(T2s, pc, "vapor") + m.h_off) - h1) / ETA
    T2 = _solve_T(lambda T: m._h_mass(T, pc, "vapor") + m.h_off - h2, Tc + 0.01)
    h3 = m._h_mass(Tc - SUBCOOL, pc, "liquid") + m.h_off
    return (h1 - h3) / (h2 - h1), T2 - 273.15


nist = {x.name: x for x in M.methods}["NIST"]

Tcond = AMBIENTS + APPROACH
Te_K = T_EVAP + 273.15
carnot = np.array([Te_K / ((t + 273.15) - Te_K) for t in Tcond])
helm = np.array([cop_helmholtz(t)[0] for t in Tcond])
pr = np.array([cop_pr(nist, t)[0] for t in Tcond])
helm_T = np.array([cop_helmholtz(t)[1] for t in Tcond])
pr_T = np.array([cop_pr(nist, t)[1] for t in Tcond])

print(f"R-32 closed form   T_evap = {T_EVAP} C   superheat {SUPERHEAT} K   "
      f"subcool {SUBCOOL} K\n")
print(f"{'T_amb':>6}{'T_cond':>8}{'Carnot':>9}{'Helmholtz':>11}{'PR/NIST':>10}"
      f"{'PR-Helm %':>11}{'T_dis Helm':>12}{'T_dis PR':>10}")
for a, t, c, h, p, ht, pt in zip(AMBIENTS, Tcond, carnot, helm, pr, helm_T, pr_T):
    print(f"{a:>6.0f}{t:>8.0f}{c:>9.4f}{h:>11.4f}{p:>10.4f}"
          f"{100*(p-h)/h:>11.2f}{ht:>12.1f}{pt:>10.1f}")

# sanity: both must sit under Carnot and fall monotonically
assert (helm < carnot).all() and (pr < carnot).all(), "COP exceeds Carnot"
assert (np.diff(helm) < 0).all() and (np.diff(pr) < 0).all(), "not monotonic"
print("\nchecks: both series below Carnot, both monotonically decreasing")

fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.plot(AMBIENTS, carnot, ":", color="0.45", marker="^", ms=4.5, mfc="none",
        mew=1.1, lw=1.3, label="Carnot limit")
ax.plot(AMBIENTS, helm, "-", color="k", marker="o", ms=5, mew=1.1, lw=2.4,
        label="Helmholtz EoS / CoolProp")
ax.plot(AMBIENTS, pr, "--", color="tab:blue", marker="s", ms=4.5, mfc="none",
        mew=1.1, lw=1.3, label="PR EoS / NIST")

ax.set_xlabel("Ambient temperature [deg C]")
ax.set_ylabel("COP")
ax.set_xticks(AMBIENTS)
_lo = min(helm.min(), pr.min())
_hi = carnot.max()
_pad = 0.05 * (_hi - _lo)
ax.set_ylim(_lo - _pad, _hi + _pad)
ax.set_xlim(AMBIENTS.min() - 1.0, AMBIENTS.max() + 1.0)
ax.grid(True, ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
fig.tight_layout()
fig.savefig("cop_vs_ambient_closedform_1014.png", dpi=200)
plt.close(fig)
print("\nwrote cop_vs_ambient_closedform_1014.png")
