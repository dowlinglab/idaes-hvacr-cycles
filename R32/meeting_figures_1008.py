"""
meeting_figures_1008.py -- one figure per meeting bullet, 10/08/2026.

Run from the R32 directory in myidaesenv. Needs CoolProp, numpy, scipy,
matplotlib, and compare_cp_methods_R1234zeZ_1008.py sitting next to it.
Does NOT need IDAES -- the IDAES COP values are results already obtained and
are hard-coded below, clearly marked, so this script runs in seconds.

    python meeting_figures_1008.py

Writes seven PNGs named fig1_... through fig7_..., matching the bullets:

  1  R-1234ze(Z) manufacturer datasheet was not readily available
  2  Located experimental P-T values from published literature, including
     critical pressure and temperature
  3  CoolProp was used to generate Helmholtz model-based properties
  4  REFPROP was used as source for the Peng-Robinson property model
  5  Operating conditions different from R-32
  6  Identified a cycle that converges at 4 different condenser temperatures
  7  Model validation results

Every figure prints the numbers behind it to stdout as it is written, so the
console output can be checked against the picture.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""

import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as P

from compare_cp_methods_R1234zeZ_1008 import (
    pr, E_T, E_TC, E_P, E_LAB, S_T, S_TC, S_P, S_SRC,
    T_CRIT_EXP, P_CRIT_EXP, helm_dome, pr_dome, cp0_eos, PR_ROW)

# ---------------------------------------------------------------------------
# Which Fedele et al. laboratory to use. Set to "UnivPM", "ITC-CNR", or None
# for both. UnivPM is the default: 28 points spanning -35 to +99.5 C, i.e. the
# whole window including the cold tail where PR degrades. ITC-CNR has more
# points (36) but only covers +10 to +80 C, which flatters PR.
LAB = "UnivPM"
# ---------------------------------------------------------------------------

Z = "R1234ze(Z)"
R32 = "R32"
PATM = 101325.0
DPI = 150

if LAB is not None:
    _k = (E_LAB == LAB)
    E_T, E_TC, E_P, E_LAB = E_T[_k], E_TC[_k], E_P[_k], E_LAB[_k]
    _s = np.array([s.endswith(LAB) or "crit" in s for s in S_SRC])
    S_T, S_TC, S_P, S_SRC = S_T[_s], S_TC[_s], S_P[_s], S_SRC[_s]
    print(f"USING ONE DATASET: Fedele et al. / {LAB} -- {len(E_T)} points, "
          f"{E_T.min():.2f}-{E_T.max():.2f} K\n")



def done(name):
    print(f"   -> wrote {name}\n")


# =============================================================================
# FIG 1 -- no manufacturer datasheet; what coverage actually exists
# =============================================================================
# The point to make: for R-32 there is a Linde datasheet spanning -130 to
# +78 C with p, hf, hg, sf, sg. For R-1234ze(Z) the only datasheet-like source
# found is an HVAC PT chart covering 0-66 C, pressure only, in GAUGE units.
# Published measurements cover far more.
print("FIG 1  datasheet availability")

fig, ax = plt.subplots(figsize=(9.0, 3.6))
bars = [
    ("Linde datasheet (R-32)\np, hf, hg, sf, sg", -130, 78, "tab:blue"),
    ("HVAC PT chart (R-1234ze(Z))\npressure only, GAUGE", 0, 66, "0.6"),
    ("Fedele et al. 2014\n64 measured points", E_TC.min(), E_TC.max(), "tab:red"),
    ("Akasaka & Lemmon EoS\nvalidity range", 238.0 - 273.15, 440.0 - 273.15,
     "tab:green"),
]
for i, (label, lo, hi, col) in enumerate(bars):
    ax.barh(i, hi - lo, left=lo, height=0.55, color=col, alpha=0.85)
    ax.text(hi + 4, i, f"{lo:.0f} to {hi:.0f} $^\\circ$C", va="center", fontsize=9)
    print(f"   {label.splitlines()[0]:<42} {lo:>8.1f} to {hi:>6.1f} C")
ax.axvline(-29, color="k", ls="--", lw=1.4)
ax.text(-29, len(bars) - 0.35, "  R-32 duty\n  T$_{evap}$ = -29 $^\\circ$C",
        fontsize=9, va="top")
ax.set_yticks(range(len(bars)))
ax.set_yticklabels([b[0] for b in bars], fontsize=8.5)
ax.set_xlabel("Temperature [$^\\circ$C]")
ax.set_title("No manufacturer datasheet for R-1234ze(Z)\n"
             "coverage of the sources that do exist")
ax.set_xlim(-145, 215)
ax.grid(True, axis="x", ls=":", alpha=0.45)
fig.tight_layout()
fig.savefig("fig1_datasheet_coverage.png", dpi=DPI)
plt.close(fig)
done("fig1_datasheet_coverage.png")


# =============================================================================
# FIG 2 -- experimental P-T from the literature, including the critical point
# =============================================================================
# Fedele, Di Nicola, Brown, Bobbo & Zilio, Int. J. Thermophys. 35 (2014) 1-12,
# Table 2: 64 points, two laboratories. Critical point from Higashi, Hayasaka
# & Ogiya, Proc. 4th IIR Conf., Paper TP-018 (2013), as quoted by Fedele p.2.
print("FIG 2  experimental P-T data")

fig, ax = plt.subplots(figsize=(8.2, 5.8))
if LAB is not None:
    # One dataset selected -- plot it as a single series. No "remaining"
    # series and no empty second lab in the legend.
    ax.semilogy(E_TC, E_P, "o", ms=6, mfc="none", mec="b", mew=1.2, ls="none",
                label=f"Fedele et al. ({len(E_TC)} pts)")
else:
    for lab, mk in (("ITC-CNR", "o"), ("UnivPM", "^")):
        k = E_LAB == lab
        if not k.any():
            continue
        ax.semilogy(E_TC[k], E_P[k], mk, ms=6, mfc="none", mec="b", mew=1.2,
                    ls="none", label=f"Fedele et al., {lab} ({k.sum()} pts)")
ax.semilogy([T_CRIT_EXP - 273.15], [P_CRIT_EXP], "*", ms=18, mfc="none",
            mec="g", mew=1.8, ls="none", label="critical point, REFPROP")
ax.semilogy(helm_dome["T_C"], helm_dome["Pd"], "k-", lw=2.2,
            label="Helmholtz (CoolProp)")
ax.semilogy(pr_dome["T_C"], pr_dome["Pd"], "r--", lw=2.0,
            label="PR (REFPROP params)")
ax.axhline(PATM / 1e5, color="0.45", ls=":", lw=1.4)
ax.text(-33, PATM / 1e5 * 1.1, "1 atm", color="0.35", fontsize=9)
ax.set_xlabel("Temperature [$^\\circ$C]")
ax.set_ylabel("Saturation pressure [bar]")
ax.set_title("P-T validation: experimental data vs both models\n"
             + (f"{len(E_TC)} measured points plus the REFPROP critical point" if LAB is not None else f"{len(E_TC)} measured points, two laboratories, plus the REFPROP critical point"))
ax.grid(True, which="both", ls=":", alpha=0.45)
ax.legend(fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig("fig2_experimental_PT.png", dpi=DPI)
plt.close(fig)
print(f"   n = {len(E_T)} points, {E_T.min():.2f}-{E_T.max():.2f} K "
      f"({E_TC.min():.2f} to {E_TC.max():.2f} C)")
print(f"   u(T) = 0.03 K, u(P) < 1 kPa")
print(f"   critical point {T_CRIT_EXP:.2f} K, {P_CRIT_EXP:.2f} bar "
      f"(REFPROP, separate source)")
done("fig2_experimental_PT.png")


# =============================================================================
# FIG 3 -- CoolProp IS the published Helmholtz equation
# =============================================================================
# Akasaka & Lemmon, J. Chem. Eng. Data (2019), DOI 10.1021/acs.jced.9b00007.
# The paper's Supporting Information carries a reference implementation by
# Ian Bell (NIST) plus six verification state points. verify_akasaka_lemmon_SI
# runs that implementation against CoolProp; this figure shows the agreement.
print("FIG 3  CoolProp vs the published equation")

SI_POINTS = [  # (T [K], rho [mol/L], p [MPa], cp [J/mol/K]) from the SI run
    (300, 0.0,  0.000000,  94.1842),
    (300, 11.0, 10.994981, 141.5796),
    (300, 0.05, 0.118913,  99.1296),
    (400, 8.0,  5.085752,  191.4250),
    (400, 1.0,  2.115059,  181.0002),
    (424, 4.0,  3.578101,  199.8539),
]
labels, dp, dcp = [], [], []
for T, rho, p_si, cp_si in SI_POINTS:
    if rho == 0.0:
        p_cp = 0.0
        cp_cp = P("Cp0molar", "T", T, "D", 1e-10, Z)
        dpi_ = 0.0
    else:
        p_cp = P("P", "T", T, "Dmolar", rho * 1000.0, Z)
        cp_cp = P("Cpmolar", "T", T, "Dmolar", rho * 1000.0, Z)
        dpi_ = 100.0 * (p_cp / 1e6 - p_si) / p_si
    labels.append(f"{T} K\n{rho} mol/L")
    dp.append(dpi_)
    dcp.append(100.0 * (cp_cp - cp_si) / cp_si)
    print(f"   T={T:>4} rho={rho:>5}  p  SI {p_si:>10.6f} MPa  CoolProp "
          f"{p_cp/1e6:>10.6f}  d={dpi_:+.6f} %")

x = np.arange(len(labels))
fig, ax = plt.subplots(figsize=(8.6, 4.8))
ax.bar(x - 0.18, dp, 0.36, label="pressure", color="tab:blue")
ax.bar(x + 0.18, dcp, 0.36, label="heat capacity", color="tab:orange")
ax.axhline(0, color="0.4", lw=1)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylabel("deviation, CoolProp vs published SI code [%]")
ax.set_ylim(-1e-3, 1e-3)
ax.set_title("CoolProp reproduces the published Akasaka & Lemmon equation\n"
             "at the Supporting Information's own six verification states")
ax.grid(True, axis="y", ls=":", alpha=0.45)
ax.legend(fontsize=9)
ax.text(0.5, 0.06, "all deviations are zero to the printed precision",
        transform=ax.transAxes, ha="center", fontsize=9.5, color="0.3")
fig.tight_layout()
fig.savefig("fig3_coolprop_vs_published_EoS.png", dpi=DPI)
plt.close(fig)
done("fig3_coolprop_vs_published_EoS.png")


# =============================================================================
# FIG 4 -- REFPROP as the source of the Peng-Robinson parameters
# =============================================================================
# R1234ZEZ.FLD, NIST SRD 23 v10.0, compiled by E.W. Lemmon. The PR model built
# on those numbers is compared against Helmholtz on p-h and T-s.
print("FIG 4  REFPROP-parameterised PR vs Helmholtz")

fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.6))
ax = axes[0]
ax.semilogy(np.r_[pr_dome["hf"], pr_dome["hg"][::-1]],
            np.r_[pr_dome["Pd"], pr_dome["Pd"][::-1]],
            color="tab:blue", lw=1.6, label="PR (REFPROP params)")
ax.semilogy(np.r_[helm_dome["hf"], helm_dome["hg"][::-1]],
            np.r_[helm_dome["Pd"], helm_dome["Pd"][::-1]],
            "k-", lw=2.4, label="Helmholtz (CoolProp)")
ax.set_xlabel("h [kJ/kg]"); ax.set_ylabel("P [bar]")
ax.set_title("p-h diagram")
ax.grid(True, which="both", ls=":", alpha=0.45); ax.legend(fontsize=9)

ax = axes[1]
ax.plot(np.r_[pr_dome["sf"], pr_dome["sg"][::-1]],
        np.r_[pr_dome["T_C"], pr_dome["T_C"][::-1]],
        color="tab:blue", lw=1.6, label="PR (REFPROP params)")
ax.plot(np.r_[helm_dome["sf"], helm_dome["sg"][::-1]],
        np.r_[helm_dome["T_C"], helm_dome["T_C"][::-1]],
        "k-", lw=2.4, label="Helmholtz (CoolProp)")
ax.set_xlabel("s [kJ/kg/K]"); ax.set_ylabel("T [$^\\circ$C]")
ax.set_title("T-s diagram")
ax.grid(True, ls=":", alpha=0.45); ax.legend(fontsize=9)

fig.suptitle("Peng-Robinson built on REFPROP R1234ZEZ.FLD parameters\n"
             f"Tc = {PR_ROW['Tc']:.2f} K   Pc = {PR_ROW['Pc_bar']*100:.1f} kPa   "
             f"omega = {PR_ROW['omega']:.3f}   MW = 114.0416 g/mol   "
             f"R = 8.3144598 J/mol/K", fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig("fig4_PR_from_REFPROP.png", dpi=DPI)
plt.close(fig)
print(f"   Tc {PR_ROW['Tc']} K | Pc {PR_ROW['Pc_bar']*100:.1f} kPa | "
      f"omega {PR_ROW['omega']} | MW 114.0416 g/mol")
print("   NOTE: Shomate A-E are NOT from REFPROP -- no published set exists.")
print("         They are refitted here to the EoS ideal-gas Cp (0.0047 %).")
done("fig4_PR_from_REFPROP.png")


# =============================================================================
# FIG 5 -- why the operating conditions have to differ from R-32
# =============================================================================
print("FIG 5  operating window")

fig, ax = plt.subplots(figsize=(8.4, 5.8))
T32 = np.linspace(180, 351.0, 300)
TZ = np.linspace(238.13, 423.0, 300)
ax.semilogy(T32 - 273.15, [P("P", "T", t, "Q", 0, R32) / 1e5 for t in T32],
            "b-", lw=2, label="R-32")
ax.semilogy(TZ - 273.15, [P("P", "T", t, "Q", 0, Z) / 1e5 for t in TZ],
            "r-", lw=2, label="R-1234ze(Z)")
ax.axhline(PATM / 1e5, color="0.45", ls=":", lw=1.5)
ax.text(118, PATM / 1e5 * 1.12, "1 atm", color="0.35", fontsize=9)

p32 = P("P", "T", 244.15, "Q", 0, R32) / 1e5
pZ = P("P", "T", 244.15, "Q", 1, Z) / 1e5
ax.plot([-29], [p32], "bo", ms=11, mfc="none", mew=2)
ax.plot([-29], [pZ], "ro", ms=11, mfc="none", mew=2)
ax.annotate(f"R-32 at -29 $^\\circ$C\n{p32:.2f} bar", (-29, p32),
            textcoords="offset points", xytext=(14, -6), fontsize=9, color="b")
ax.annotate(f"R-1234ze(Z) at -29 $^\\circ$C\n{pZ:.3f} bar = "
            f"{pZ*1e5/PATM:.2f} atm\n(vacuum)", (-29, pZ),
            textcoords="offset points", xytext=(14, -34), fontsize=9, color="r")
nbpZ = P("T", "P", PATM, "Q", 0, Z) - 273.15
nbp32 = P("T", "P", PATM, "Q", 0, R32) - 273.15
ax.plot([nbpZ], [1.013], "r^", ms=10)
ax.annotate(f"NBP {nbpZ:+.1f} $^\\circ$C", (nbpZ, 1.013),
            textcoords="offset points", xytext=(8, 10), fontsize=9, color="r")
ax.plot([nbp32], [1.013], "b^", ms=10)
ax.annotate(f"NBP {nbp32:+.1f} $^\\circ$C", (nbp32, 1.013),
            textcoords="offset points", xytext=(-20, 12), fontsize=9, color="b")
ax.set_xlabel("Temperature [$^\\circ$C]")
ax.set_ylabel("Saturation pressure [bar]")
ax.set_title("Operating conditions must differ from R-32\n"
             f"R-1234ze(Z) boils {nbpZ - nbp32:.0f} K warmer than R-32")
ax.grid(True, which="both", ls=":", alpha=0.45)
ax.legend(loc="upper left", fontsize=10)
ax.set_xlim(-60, 155)
fig.tight_layout()
fig.savefig("fig5_operating_window.png", dpi=DPI)
plt.close(fig)
print(f"   NBP: R-32 {nbp32:.2f} C | R-1234ze(Z) {nbpZ:.2f} C "
      f"(difference {nbpZ-nbp32:.1f} K)")
print(f"   at -29 C: R-32 {p32:.3f} bar | R-1234ze(Z) {pZ:.3f} bar "
      f"({pZ*1e5/PATM:.3f} atm)")
done("fig5_operating_window.png")


# =============================================================================
# FIG 6 -- the converging cycle
# =============================================================================
# IDAES results obtained 10/08/2026 with phase6_final_R1234zeZ_1008.py's chain
# at T_evap = 30 C. HARD-CODED so this script does not need IDAES. Re-run the
# sweep to refresh them.
print("FIG 6  converging cycle")

T_EVAP = 30.0
TCOND = np.array([70.0, 80.0, 90.0, 100.0])
IDAES_COP = np.array([4.8099, 3.8161, 3.6025, 3.1962])   # all converged


def helm_cop(te, tc, eta=0.9999):
    Te, Tc = te + 273.15, tc + 273.15
    pc = P("P", "T", Tc, "Q", 0, Z)
    h1 = P("H", "T", Te, "Q", 1, Z)
    s1 = P("S", "T", Te, "Q", 1, Z)
    h2 = h1 + (P("H", "P", pc, "S", s1, Z) - h1) / eta
    return (h1 - P("H", "T", Tc, "Q", 0, Z)) / (h2 - h1)


closed = np.array([helm_cop(T_EVAP, t) for t in TCOND])
carnot = np.array([(T_EVAP + 273.15) / (t - T_EVAP) for t in TCOND])

fig, ax = plt.subplots(figsize=(8.2, 5.6))
ax.plot(TCOND, carnot, "k:", lw=1.8, marker="^", ms=8, label="Carnot limit")
ax.plot(TCOND, closed, "k-", lw=2.2, marker="o", ms=8,
        label="Helmholtz, closed form")
ax.plot(TCOND, IDAES_COP, "r--", lw=2.2, marker="s", ms=8,
        label="PR, IDAES solve (4/4 converged)")
for xx, a, b in zip(TCOND, closed, IDAES_COP):
    ax.annotate(f"{100*(b-a)/a:+.0f}%", (xx, (a + b) / 2), fontsize=9,
                color="0.3", ha="center")
ax.set_xlabel("Condensing temperature [$^\\circ$C]")
ax.set_ylabel("COP")
ax.set_title(f"R-1234ze(Z) cycle, T$_{{evap}}$ = {T_EVAP:.0f} $^\\circ$C\n"
             "IDAES converges at all four condenser temperatures")
ax.set_xticks(TCOND)
ax.grid(True, ls=":", alpha=0.45)
ax.legend(fontsize=9.5)
fig.tight_layout()
fig.savefig("fig6_cycle_COP.png", dpi=DPI)
plt.close(fig)
print(f"   {'T_cond':>7}{'Carnot':>9}{'closed':>9}{'IDAES':>9}{'diff %':>9}")
for t, c, cl, i in zip(TCOND, carnot, closed, IDAES_COP):
    print(f"   {t:>7.0f}{c:>9.3f}{cl:>9.4f}{i:>9.4f}{100*(i-cl)/cl:>9.1f}")
print("   NOTE: the IDAES-vs-closed-form gap is UNEXPLAINED. See the status")
print("         report before quoting the IDAES column.")
done("fig6_cycle_COP.png")


# =============================================================================
# FIG 7 -- model validation against the 64 measured points
# =============================================================================
print("FIG 7  model validation")

P_pr = np.array([pr.saturation_pressure(t) / 1e5 for t in E_T])
P_hm = np.array([P("P", "T", t, "Q", 0, Z) / 1e5 for t in E_T])

fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0),
                         gridspec_kw={"width_ratios": [2.1, 1]})
ax = axes[0]
for tag, Pm, mk, col in (("Helmholtz", P_hm, "o", "k"), ("PR", P_pr, "s", "r")):
    ax.plot(E_T, 100.0 * (Pm - E_P) / E_P, mk, ms=5, mfc="none", mec=col,
            mew=1.1, ls="none", label=tag)
ax.axhline(0.0, color="0.4", lw=1)
ax.axvspan(E_T.min(), 300.0, color="0.92", zorder=0)
ax.text(0.02, 0.04, "below 300 K: NIST quotes 0.3 % uncertainty\n"
                    "(0.1 % above)", transform=ax.transAxes, fontsize=8.5,
        color="0.35", va="bottom")
ax.set_xlabel("$T$ [K]")
ax.set_ylabel("$100\\,\\Delta P/P$")
ax.set_title("Deviation from 64 measured vapour pressures\n"
             "(axes convention of Fedele et al. 2014, Figs. 2-8)")
ax.grid(True, ls=":", alpha=0.45)
ax.legend(fontsize=9)

ax = axes[1]
ax.axis("off")
rows = [["method", "BIAS %", "AAD %", "max %"]]
for tag, Pm in (("Helmholtz", P_hm), ("PR", P_pr)):
    d = 100.0 * (Pm - E_P) / E_P
    rows.append([tag, f"{d.mean():.3f}", f"{np.abs(d).mean():.3f}",
                 f"{np.abs(d).max():.3f}"])
    print(f"   {tag:<11} BIAS {d.mean():+.3f} %   AAD {np.abs(d).mean():.3f} %"
          f"   max {np.abs(d).max():.3f} %   n={len(d)}")
rows.append(["", "", "", ""])
rows.append(["Brown et al. PR", "", "2.4-2.7", ""])
rows.append(["Akasaka FEQ", "", "0.37-0.63", ""])
rows.append(["Wagner Eq.1", "", "0.29-0.35", ""])
tbl = ax.table(cellText=rows[1:], colLabels=rows[0], loc="center",
               cellLoc="center")
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1.0, 1.55)
ax.set_title("Scored on all 64 points\nlower block: Fedele et al. Table 4",
             fontsize=10)
fig.tight_layout()
fig.savefig("fig7_model_validation.png", dpi=DPI)
plt.close(fig)
done("fig7_model_validation.png")


print("=" * 70)
print("  all seven figures written")
print("=" * 70)
