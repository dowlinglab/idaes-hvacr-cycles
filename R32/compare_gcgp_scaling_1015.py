"""
compare_gcgp_scaling_1015.py -- which reading of the GCGP Shomate row is the
right one? Four PR series scored against two independent references.

THE QUESTION
    The GCGP delivery gives Cp coefficients as

        a = 14.161   b = 0.124   c = -6.340e-05   d = 1.190e-8   e = 0.0

    Read them as Shomate (Cp = A + Bt + Ct^2 + Dt^3 + E/t^2, t = T/1000) and
    Cp(300 K) = 14.20 J/mol/K. Read them as a raw polynomial in T and convert
    to Shomate (A=a, B=b*1e3, C=c*1e6, D=d*1e9, E=e*1e-6) and Cp(300 K) =
    45.98 J/mol/K. Both run. Both converge. Only one can be right.

    This is OPEN QUESTION 1 in BREADCRUMB_07-20.md, raised 09/17 and still
    waiting on Fathya. This script settles what the DATA can settle without
    her, and marks clearly what it cannot.

SERIES
    NIST            NIST Shomate + NIST Pc/Tc          (the control)
    GCGP_unscaled   coefficients used as entered        Cp(300) = 14.20
    GCGP_scaled     coefficients converted to Shomate   Cp(300) = 45.98
    SPGP            fp_single_0919 under its slide label

REFERENCES
    Linde datasheet   20 saturation points
    Helmholtz         CoolProp, Tillner-Roth & Yokozeki 1997

THREE TESTS, in increasing order of discriminating power

    1. SATURATION PROPERTIES.  Cannot tell the two GCGP readings apart on
       pressure or on latent heat, and that is not a defect in the test --
       it is the physics. Saturation pressure comes from Ambrose-Walton,
       which sees only (Tc, Pc, omega); Cp never enters. Latent heat is
       hg - hf at one temperature, and the ideal-gas Cp integral is the SAME
       integral on both branches, so it cancels exactly. The two readings
       therefore give IDENTICAL P and IDENTICAL hfg, and differ only in the
       absolute h and s, which the IIR anchoring then partly absorbs.

       This is why the p-h diagram cannot settle the question. Dome height is
       pressure; dome width is latent heat. Both are blind to Cp.

    2. IDEAL-GAS Cp AGAINST A PHYSICAL FLOOR.  A non-linear polyatomic
       molecule has at least 3/2 R of translation and 3/2 R of rotation, plus
       R for Cp - Cv, so Cp0 >= 4R = 33.258 J/mol/K for ANY such molecule at
       ANY temperature, before a single vibrational mode is counted. R-32
       (CH2F2, 5 atoms, 9 vibrational modes) is non-linear. This is a
       theoretical limit, not a fit, and nothing about the GCGP method can
       move it.

    3. COMPRESSOR DISCHARGE TEMPERATURE.  The closed-form cycle at the poster
       operating point. Cp sets how much the vapour heats on compression, so
       here the two readings separate by hundreds of degrees and one of them
       leaves the physical range entirely.

Test 1 is the one that has been run for months. Tests 2 and 3 are new here.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/15/2026
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

import compare_cp_methods_0929_Helmholtz_SPGPlabel as M
from compare_cp_methods_0929_Helmholtz_SPGPlabel import PRMethod

FLUID = "R32"
R = 8.314                       # J/mol/K, as the rest of the R-32 chain uses
CP0_FLOOR = 4.0 * R             # 33.258 J/mol/K, non-linear polyatomic minimum

# Linde reference arrays, re-used rather than re-typed
L_T, L_P = M.L_T, M.L_P
L_HF, L_HG, L_SF, L_SG = M.L_HF, M.L_HG, M.L_SF, M.L_SG


# =============================================================================
# The four series
# =============================================================================
# GCGP row as delivered.  Converting to Shomate is A=a, B=b*1e3, C=c*1e6,
# D=d*1e9, E=e*1e-6 -- the same conversion applied to the SPGP row in
# compare_cp_methods_0929_Helmholtz_SPGPlabel.py, where the sheet says
# explicitly that its a-e are raw. The GCGP sheet does not say.
G_RAW = dict(Pc_bar=50.730, Tc=355.354, omega=0.1711)
G_ABCDE_ASIS = (14.161, 0.124, -6.340e-05, 1.190e-8, 0.0)
G_ABCDE_CONV = (14.161, 124.0, -63.40, 11.90, 0.0)

_by_name = {x.name: x for x in M.methods}

series = [
    ("NIST",          _by_name["NIST"],                                  "tab:blue",  "-"),
    ("GCGP_unscaled", PRMethod("GCGP_unscaled", *(
        (G_RAW["Pc_bar"], G_RAW["Tc"]) + G_ABCDE_ASIS),
        omega=G_RAW["omega"]),                                           "tab:orange", "--"),
    ("GCGP_scaled",   PRMethod("GCGP_scaled", *(
        (G_RAW["Pc_bar"], G_RAW["Tc"]) + G_ABCDE_CONV),
        omega=G_RAW["omega"]),                                           "tab:green",  "-."),
    ("SPGP",          _by_name["SPGP"],                                  "tab:red",   ":"),
]
NAMES = [s[0] for s in series]
MODELS = {n: m for n, m, _, _ in series}
COLOR = {n: c for n, _, c, _ in series}
STYLE = {n: s for n, _, _, s in series}


def cp0(m, T):
    """Ideal-gas Cp from a method's Shomate row, J/mol/K."""
    t = np.asarray(T, float) / 1000.0
    return m.A + m.B*t + m.C*t**2 + m.D*t**3 + np.divide(
        m.E, t**2, out=np.zeros_like(t), where=(m.E != 0))


def cp0_helm(T):
    return np.array([CP.PropsSI("Cp0molar", "T", t, "P", 101325.0, FLUID)
                     for t in np.atleast_1d(T)])


# =============================================================================
# TEST 1 -- saturation properties vs Linde and vs Helmholtz
# =============================================================================

def _sat_pack(m, T):
    return (np.array([m.saturation_pressure(t)/1e5 for t in T]),
            np.array([m.hf(t) for t in T]),
            np.array([m.hg(t) for t in T]),
            np.array([m.sf(t) for t in T]),
            np.array([m.sg(t) for t in T]))


def _sat_pack_helm(T):
    return (np.array([CP.PropsSI("P", "T", t, "Q", 0, FLUID)/1e5 for t in T]),
            np.array([CP.PropsSI("H", "T", t, "Q", 0, FLUID)/1e3 for t in T]),
            np.array([CP.PropsSI("H", "T", t, "Q", 1, FLUID)/1e3 for t in T]),
            np.array([CP.PropsSI("S", "T", t, "Q", 0, FLUID)/1e3 for t in T]),
            np.array([CP.PropsSI("S", "T", t, "Q", 1, FLUID)/1e3 for t in T]))


def test1():
    print("=" * 78)
    print("TEST 1   saturation properties")
    print("=" * 78)
    for ref_name in ("Linde datasheet", "Helmholtz (CoolProp)"):
        print(f"\nvs {ref_name}, at the 20 Linde temperatures")
        print(f"{'series':>14}{'P MAPE%':>10}{'hf MAE':>9}{'hg MAE':>9}"
              f"{'hfg MAE':>9}{'sf MAE':>10}{'sg MAE':>10}{'n':>6}")
        print(f"{'':>14}{'':>10}{'kJ/kg':>9}{'kJ/kg':>9}{'kJ/kg':>9}"
              f"{'kJ/kgK':>10}{'kJ/kgK':>10}")
        for n in NAMES:
            m = MODELS[n]
            mask = L_T < 0.99 * m.Tc
            T = L_T[mask]
            P, hf, hg, sf, sg = _sat_pack(m, T)
            if ref_name.startswith("Linde"):
                rP, rhf, rhg = L_P[mask], L_HF[mask], L_HG[mask]
                rsf, rsg = L_SF[mask], L_SG[mask]
            else:
                rP, rhf, rhg, rsf, rsg = _sat_pack_helm(T)
            print(f"{n:>14}{np.mean(np.abs((P-rP)/rP))*100:>10.2f}"
                  f"{np.mean(np.abs(hf-rhf)):>9.2f}"
                  f"{np.mean(np.abs(hg-rhg)):>9.2f}"
                  f"{np.mean(np.abs((hg-hf)-(rhg-rhf))):>9.2f}"
                  f"{np.mean(np.abs(sf-rsf)):>10.4f}"
                  f"{np.mean(np.abs(sg-rsg)):>10.4f}"
                  f"{int(mask.sum()):>4}/{len(L_T)}")

    # The identity that makes test 1 useless for this question, stated as a
    # check rather than as a claim.
    a, b = MODELS["GCGP_unscaled"], MODELS["GCGP_scaled"]
    T = L_T[L_T < 0.99 * a.Tc]
    dP = np.max(np.abs([a.saturation_pressure(t) - b.saturation_pressure(t)
                        for t in T]))
    dfg = np.max(np.abs([(a.hg(t) - a.hf(t)) - (b.hg(t) - b.hf(t)) for t in T]))
    print(f"\n  max |P_unscaled - P_scaled|   = {dP:.3e} Pa")
    print(f"  max |hfg_unscaled - hfg_scaled| = {dfg:.3e} kJ/kg")
    assert dP < 1e-6, "pressure should not depend on Cp"
    assert dfg < 1e-6, "latent heat should not depend on Cp"
    print("  -> the two readings are INDISTINGUISHABLE on pressure and on")
    print("     latent heat. No p-h or T-s diagram can separate them.")


# =============================================================================
# TEST 2 -- ideal-gas Cp against the 4R floor
# =============================================================================

def test2():
    print("\n" + "=" * 78)
    print("TEST 2   ideal-gas Cp against the 4R floor")
    print("=" * 78)
    print(f"\n  4R = {CP0_FLOOR:.3f} J/mol/K is the minimum Cp0 of ANY non-linear")
    print("  polyatomic molecule: 3/2 R translation + 3/2 R rotation + R for")
    print("  Cp - Cv, with every vibrational mode still frozen out. R-32 is")
    print("  CH2F2, non-linear, 9 vibrational modes.\n")

    T = np.array([250.0, 300.0, 350.0, 400.0])
    print(f"{'series':>14}" + "".join(f"{f'Cp0({t:.0f}K)':>13}" for t in T)
          + f"{'min over 200-500 K':>22}")
    Tw = np.linspace(200.0, 500.0, 301)
    rows = {}
    for n in NAMES:
        v = cp0(MODELS[n], T)
        vw = cp0(MODELS[n], Tw)
        rows[n] = vw
        print(f"{n:>14}" + "".join(f"{x:>13.2f}" for x in v)
              + f"{vw.min():>22.2f}")
    vh = cp0_helm(T)
    rows["Helmholtz"] = cp0_helm(Tw)
    print(f"{'Helmholtz':>14}" + "".join(f"{x:>13.2f}" for x in vh)
          + f"{rows['Helmholtz'].min():>22.2f}")
    print(f"{'4R floor':>14}" + "".join(f"{CP0_FLOOR:>13.2f}" for _ in T))

    print()
    for n in NAMES:
        below = rows[n].min() < CP0_FLOOR
        tag = "VIOLATES the 4R floor" if below else "physical"
        err = np.mean(np.abs(rows[n] - rows["Helmholtz"])
                      / rows["Helmholtz"]) * 100
        print(f"  {n:<14} {tag:<24} MAPE vs Helmholtz Cp0 = {err:6.2f} %")

    assert rows["GCGP_unscaled"].min() < CP0_FLOOR, \
        "unscaled GCGP was expected to fall below 4R"
    assert rows["GCGP_scaled"].min() > CP0_FLOOR, \
        "scaled GCGP was expected to clear 4R"
    return rows, Tw


# =============================================================================
# TEST 3 -- compressor discharge temperature, closed form
# =============================================================================
# Poster cycle: T_evap = -29 C, T_cond = T_amb + 9 C, superheat 3 K,
# subcool 8 K, isentropic efficiency 0.9999. Four state-point lookups, no
# solver, so nothing here can be blamed on the flowsheet defect of 10/08.

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.9999
SUPERHEAT, SUBCOOL = 3.0, 8.0
AMBIENTS = np.array([10.0, 15.0, 20.0, 25.0])


def _solve_T(f, lo):
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


def cycle_pr(m, T_cond_C):
    Te, Tc = T_EVAP + 273.15, T_cond_C + 273.15
    pe, pc = m.saturation_pressure(Te), m.saturation_pressure(Tc)
    h1 = m._h_mass(Te + SUPERHEAT, pe, "vapor") + m.h_off
    s1 = m._s_mass(Te + SUPERHEAT, pe, "vapor") + m.s_off
    T2s = _solve_T(lambda T: m._s_mass(T, pc, "vapor") + m.s_off - s1, Tc + 0.01)
    h2 = h1 + ((m._h_mass(T2s, pc, "vapor") + m.h_off) - h1) / ETA
    T2 = _solve_T(lambda T: m._h_mass(T, pc, "vapor") + m.h_off - h2, Tc + 0.01)
    h3 = m._h_mass(Tc - SUBCOOL, pc, "liquid") + m.h_off
    return (h1 - h3) / (h2 - h1), T2 - 273.15


def cycle_helm(T_cond_C):
    Te, Tc = T_EVAP + 273.15, T_cond_C + 273.15
    pe = CP.PropsSI("P", "T", Te, "Q", 1, FLUID)
    pc = CP.PropsSI("P", "T", Tc, "Q", 0, FLUID)
    h1 = CP.PropsSI("H", "T", Te + SUPERHEAT, "P", pe, FLUID)
    s1 = CP.PropsSI("S", "T", Te + SUPERHEAT, "P", pe, FLUID)
    h2 = h1 + (CP.PropsSI("H", "P", pc, "S", s1, FLUID) - h1) / ETA
    h3 = CP.PropsSI("H", "T", Tc - SUBCOOL, "P", pc, FLUID)
    return (h1 - h3) / (h2 - h1), CP.PropsSI("T", "P", pc, "H", h2, FLUID) - 273.15


def test3():
    print("\n" + "=" * 78)
    print("TEST 3   compressor discharge temperature, closed form")
    print("=" * 78)
    print(f"\n  T_evap = {T_EVAP} C, T_cond = T_amb + {APPROACH} C, "
          f"superheat {SUPERHEAT} K, subcool {SUBCOOL} K, eta = {ETA}\n")
    Tcond = AMBIENTS + APPROACH
    res = {}
    print(f"{'series':>14}" + "".join(f"{f'T_dis @{a:.0f}C':>14}" for a in AMBIENTS)
          + f"{'COP @20C':>11}")
    for n in NAMES + ["Helmholtz"]:
        T2, cops = [], []
        for t in Tcond:
            try:
                c, d = (cycle_helm(t) if n == "Helmholtz"
                        else cycle_pr(MODELS[n], t))
            except Exception:
                c, d = np.nan, np.nan
            cops.append(c); T2.append(d)
        res[n] = (np.array(cops), np.array(T2))
        print(f"{n:>14}" + "".join(f"{x:>14.1f}" if np.isfinite(x) else f"{'--':>14}"
                                   for x in T2)
              + f"{cops[2]:>11.4f}")

    print("\n  Helmholtz is the physical answer. A discharge temperature far")
    print("  above it for the same compression means the Cp is wrong, because")
    print("  the pressure ratio is set by (Tc, Pc, omega) and is shared.")
    for n in NAMES:
        d = res[n][1] - res["Helmholtz"][1]
        print(f"    {n:<14} T_dis - T_dis(Helmholtz) = "
              f"{np.nanmin(d):+7.1f} to {np.nanmax(d):+7.1f} K")
    return res


# =============================================================================
# Figures
# =============================================================================

def figures(cp_rows, Tw, cyc):
    # --- fig A: the discriminating one -----------------------------------
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    for n in NAMES:
        ax.plot(Tw, cp_rows[n], STYLE[n], color=COLOR[n], lw=1.8, label=n)
    ax.plot(Tw, cp_rows["Helmholtz"], "-", color="k", lw=2.4,
            label="Helmholtz (CoolProp)")
    ax.axhline(CP0_FLOOR, color="0.3", lw=1.4, ls=(0, (6, 3)))
    ax.fill_between(Tw, 0, CP0_FLOOR, color="0.85", alpha=0.6, zorder=0)
    ax.annotate(f"4R = {CP0_FLOOR:.1f} J/mol/K\nphysical floor, non-linear polyatomic",
                xy=(Tw[len(Tw)//2], CP0_FLOOR), xytext=(0, -34),
                textcoords="offset points", ha="center", fontsize=7.5,
                color="0.25")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel(r"Ideal-gas $C_p^0$ [J/mol/K]")
    ax.set_xlim(Tw.min(), Tw.max())
    ax.set_ylim(0, max(90.0, max(cp_rows[n].max() for n in NAMES) * 1.08))
    ax.grid(True, ls=":", alpha=0.35)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    fig.tight_layout()
    fig.savefig("gcgp_scaling_cp0_1015.png", dpi=200)
    plt.close(fig)

    # --- fig B: the test that cannot tell them apart ----------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    Td = np.linspace(M.T_MIN, 0.999 * MODELS["NIST"].Tc, 120)
    for ax, (xf, yf, xl, yl, logy) in zip(axes, [
            (lambda m, T: (np.array([m.hf(t) for t in T]),
                           np.array([m.hg(t) for t in T])),
             lambda m, T: np.array([m.saturation_pressure(t)/1e5 for t in T]),
             "h [kJ/kg]", "P [bar]", True),
            (lambda m, T: (np.array([m.sf(t) for t in T]),
                           np.array([m.sg(t) for t in T])),
             lambda m, T: T - 273.15, "s [kJ/kg/K]", "T [deg C]", False)]):
        for n in NAMES:
            m = MODELS[n]
            Tg = np.linspace(M.T_MIN, 0.999 * m.Tc, 120)
            lo, hi = xf(m, Tg)
            y = yf(m, Tg)
            ax.plot(np.concatenate([lo, hi[::-1]]),
                    np.concatenate([y, y[::-1]]), STYLE[n], color=COLOR[n],
                    lw=1.8, label=n)
        hd = M.helm_dome
        ax.plot(np.concatenate([hd["hf" if logy else "sf"],
                                hd["hg" if logy else "sg"][::-1]]),
                np.concatenate([hd["Pd"], hd["Pd"][::-1]]) if logy
                else np.concatenate([hd["T_C"], hd["T_C"][::-1]]),
                "-", color="k", lw=2.4, label="Helmholtz (CoolProp)")
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.grid(True, which="both", ls=":", alpha=0.35)
        ax.legend(fontsize=8, framealpha=0.9)
    axes[0].set_title("Same dome HEIGHT (P identical to 0 Pa) and same WIDTH\n"
                      "(hfg identical to 1e-13 kJ/kg) -- only the h datum moves",
                      fontsize=9)
    axes[1].set_title("T-s: same shape, shifted along s.  Neither figure can\n"
                      "separate the two readings -- see TEST 2 and TEST 3",
                      fontsize=9)
    fig.tight_layout()
    fig.savefig("gcgp_scaling_domes_1015.png", dpi=200)
    plt.close(fig)

    # --- fig C: discharge temperature -------------------------------------
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    for n in NAMES:
        ax.plot(AMBIENTS, cyc[n][1], STYLE[n], color=COLOR[n], marker="o",
                ms=5, mfc="none", mew=1.3, lw=1.8, label=n)
    ax.plot(AMBIENTS, cyc["Helmholtz"][1], "-", color="k", marker="o", ms=5,
            lw=2.4, label="Helmholtz (CoolProp)")
    ax.axhspan(60, 130, color="0.88", zorder=0)
    ax.annotate("typical R-32 discharge range", xy=(AMBIENTS.mean(), 125),
                ha="center", va="top", fontsize=7.5, color="0.25")
    ax.set_xlabel("Ambient temperature [deg C]")
    ax.set_ylabel("Compressor discharge temperature [deg C]")
    ax.set_xticks(AMBIENTS)
    ax.set_yscale("symlog", linthresh=100)
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    fig.tight_layout()
    fig.savefig("gcgp_scaling_discharge_1015.png", dpi=200)
    plt.close(fig)

    print("\nwrote gcgp_scaling_cp0_1015.png")
    print("      gcgp_scaling_domes_1015.png")
    print("      gcgp_scaling_discharge_1015.png")


# =============================================================================

if __name__ == "__main__":
    test1()
    cp_rows, Tw = test2()
    cyc = test3()
    figures(cp_rows, Tw, cyc)

    print("\n" + "=" * 78)
    print("WHAT THIS SETTLES AND WHAT IT DOES NOT")
    print("=" * 78)
    print("""
  SETTLED by the data:
    The unscaled reading is not physically admissible. Cp0 = 14.2 J/mol/K
    at 300 K is below 4R = 33.26, which no non-linear polyatomic molecule
    can be, and the closed-form discharge temperature it implies is far
    outside any range a real R-32 compressor operates in. The scaled
    reading clears the floor and tracks Helmholtz's Cp0.

  NOT SETTLED, and only Fathya can settle it:
    Whether the GCGP sheet INTENDED raw coefficients in T. The physics says
    the scaled reading is the admissible one; it does not say the sheet was
    written that way, and a third possibility -- that the coefficients are
    correct as Shomate but for a different quantity, or carry a unit other
    than J/mol/K -- is not excluded by anything here. Report the scaled
    reading with the convention stated explicitly until she confirms.

  CARRIED FORWARD: OPEN QUESTION 1, BREADCRUMB_07-20.md, raised 09/17.
""")
