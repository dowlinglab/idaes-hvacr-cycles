"""
Compare p-h and T-s diagrams for R-32 across several Cp / critical-property
methods, benchmarked against the Linde datasheet. This code adds Helmholtz 
model predictions as a baseline.

Five data series are shown on each diagram:
    1. Linde datasheet         (reference; from NIST REFPROP 9.0)
    2. NIST                     (NIST Shomate + NIST Pc/Tc)
    3. GCGP                     (GCGP Shomate + GCGP Pc/Tc)
    4. SPGP                     (SPGP Shomate + SPGP predicted Tc)
    5. Helmholtz model          (real R-32 EOS, via CoolProp -- same
                                  Tillner-Roth & Yokozeki 1997 correlation
                                  IDAES's general_helmholtz package uses)

Each method supplies only its ROW DATA: Pc [bar], Tc [K], and Shomate A..E.
Shared/fixed: R-32 molar mass, gas constant, the Linde saturation table
(pressure column used to back out each method's acentric factor at Tr = 0.7;
enthalpy/entropy columns used as the validation reference), and the IIR
reference state (sat. liquid at 0 C: h = 200 kJ/kg, s = 1.0 kJ/kg/K).

Model per method:
    1. Ambrose-Walton saturation pressure (Poling et al., 2001, sec. 7-4).
    2. Peng-Robinson cubic EOS departure h and s (Poling Table 6-3).
    3. Shomate ideal-gas Cp integration (t = T/1000).

Author: Shilpa Narasimhan and Claude AI
Date Created: 08/11/2026
QA/testing: Shilpa Narasimhan
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP


# =============================================================================
# Shared constants
# =============================================================================

R = 8.314          # J/mol/K
MW = 52.023        # kg/kmol (numerically g/mol), R-32

u = 2.0
w = -1.0
Omega_A = 0.45724
Omega_B = 0.0778

T_REF = 273.15     # K       IIR reference temperature (0 C)
H_REF = 200.0      # kJ/kg   IIR reference enthalpy (sat. liquid)
S_REF = 1.0        # kJ/kg/K IIR reference entropy  (sat. liquid)

T_MIN = -130.0 + 273.15   # K, low end of the dome sweep


# =============================================================================
# Linde saturation table (NIST REFPROP 9.0) at the selected points
#   (T[C], Psat[bar], hf[kJ/kg], hg[kJ/kg], sf[kJ/kg/K], sg[kJ/kg/K])
#   hf/hg = saturated liquid/vapor enthalpy ; sf/sg = liquid/vapor entropy
# =============================================================================

LINDE_SAT = [
    (-130, 0.001312,  -8.26, 448.77, -0.028, 3.165),
    (-110, 0.014525,  23.20, 461.86,  0.178, 2.867),
    (-90,  0.07556,   54.42, 474.61,  0.359, 2.653),
    (-70,  0.36067,   85.66, 486.57,  0.520, 2.494),
    (-50,  1.014,    117.22, 497.27,  0.668, 2.371),
    (-30,  2.7344,   149.45, 506.27,  0.806, 2.274),
    (-10,  5.8263,   182.76, 513.02,  0.937, 2.192),
    (0,    8.131,    200.00, 515.30,  1.000, 2.154),
    (10,   10.065,   217.74, 516.66,  1.063, 2.119),
    (20,   14.746,   236.12, 516.90,  1.125, 2.083),
    (30,   19.275,   255.32, 515.72,  1.188, 2.047),
    (40,   24.783,   275.61, 512.71,  1.252, 2.009),
    (50,   31.412,   297.49, 507.10,  1.318, 1.967),
    (58,   37.635,   316.75, 499.82,  1.375, 1.928),
    (62,   41.089,   327.30, 494.76,  1.405, 1.905),
    (66,   44.793,   338.78, 488.26,  1.438, 1.879),
    (70,   48.768,   351.73, 479.52,  1.474, 1.846),
    (74,   53.046,   367.53, 466.41,  1.518, 1.803),
    (76,   55.315,   378.03, 455.86,  1.547, 1.770),
    (78,   57.697,   400.38, 428.90,  1.610, 1.691),
]

# Linde arrays for convenience
L_TC = np.array([r[0] for r in LINDE_SAT])           # deg C
L_T = L_TC + 273.15                                  # K
L_P = np.array([r[1] for r in LINDE_SAT])            # bar
L_HF = np.array([r[2] for r in LINDE_SAT])
L_HG = np.array([r[3] for r in LINDE_SAT])
L_SF = np.array([r[4] for r in LINDE_SAT])
L_SG = np.array([r[5] for r in LINDE_SAT])


# =============================================================================
# Method row data (the only per-method input)
# =============================================================================

METHODS = [
    {"name": "NIST", "Pc_bar": 57.85,   "Tc": 351.3,
     "A": -6.098682, "B": 179.2200, "C": -122.3682, "D": 32.30207, "E": 0.491361},
    {"name": "GCGP", "Pc_bar": 50.730,  "Tc": 355.354,
     "A": 14.161,    "B": 0.124,    "C": -6.340e-05, "D": 1.190e-8, "E": 0.0},
    {"name": "SPGP", "Pc_bar": 50.8106, "Tc": 400.898,
     "A": 129.687,   "B": 171.303,  "C": 146.2,      "D": 61.9837,  "E": -0.0000361638},
]
# SPGP hybrid variants (SPGP_NIST, SPGP_Fathya, SPGP_Shomate_GCGP,
# SPGP_NIST_crit_prps, SPGP_Fathya_crit_prps) removed per the updated header
# above -- this file now tracks exactly the five series it documents:
# Linde, NIST, GCGP, SPGP, Helmholtz.


# =============================================================================
# Property model for one method
# =============================================================================

class PRMethod:
    """Peng-Robinson + Shomate property model for one (Pc, Tc, Shomate) row."""

    def __init__(self, name, Pc_bar, Tc, A, B, C, D, E):
        self.name = name
        self.Pc = Pc_bar * 1e5
        self.Tc = Tc
        self.A, self.B, self.C, self.D, self.E = A, B, C, D, E

        self.omega = self._omega_from_linde()
        self.kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        self.b = Omega_B * R * self.Tc / self.Pc

        # IIR anchoring offsets (sat. liquid at 0 C)
        P0 = self.saturation_pressure(T_REF)
        self.h_off = H_REF - self._h_mass(T_REF, P0, "liquid")
        self.s_off = S_REF - self._s_mass(T_REF, P0, "liquid")

    # ---- acentric factor from Linde at Tr = 0.7 ------------------------
    def _omega_from_linde(self):
        T_target = 0.7 * self.Tc
        lnP = np.interp(1.0 / T_target, (1.0 / L_T)[::-1], np.log(L_P * 1e5)[::-1])
        return -1.0 - np.log10(np.exp(lnP) / self.Pc)

    # ---- saturation pressure (Ambrose-Walton) --------------------------
    def saturation_pressure(self, T):
        Tr = T / self.Tc
        tau = 1.0 - Tr
        f0 = (-5.97616*tau + 1.29874*tau**1.5 - 0.60394*tau**2.5 - 1.06841*tau**5) / Tr
        f1 = (-5.03365*tau + 1.11505*tau**1.5 - 5.41217*tau**2.5 - 7.46628*tau**5) / Tr
        f2 = (-0.64771*tau + 2.41539*tau**1.5 - 4.26979*tau**2.5 + 3.25259*tau**5) / Tr
        return self.Pc * np.exp(f0 + self.omega*f1 + self.omega**2*f2)

    # ---- Peng-Robinson parameters --------------------------------------
    def a_param(self, T):
        return Omega_A * R**2 * self.Tc**2 * (1 + self.kappa*(1 - np.sqrt(T/self.Tc)))**2 / self.Pc

    def dadT(self, T):
        Tr = T / self.Tc
        return Omega_A * R**2 * self.Tc**2 / self.Pc * (
            -self.kappa * (1 + self.kappa*(1 - np.sqrt(Tr))) / (self.Tc*np.sqrt(Tr)))

    def z_roots(self, T, P):
        a = self.a_param(T)
        A = a*P/(R**2*T**2)
        B = self.b*P/(R*T)
        coeffs = [1.0, -(1.0 + B - u*B), A - u*B - (u - w)*B**2, -(A*B + w*B**2 + w*B**3)]
        roots = sorted(r.real for r in np.roots(coeffs) if abs(r.imag) < 1e-9 and r.real > B)
        if not roots:
            raise ValueError(f"{self.name}: no Z root at T={T:.4g} K, P={P:.4g} Pa")
        return roots, A, B

    # ---- Shomate ideal-gas terms ---------------------------------------
    def _int_h(self, T):
        t = T/1000.0
        return 1000.0*(self.A*t + self.B*t**2/2 + self.C*t**3/3 + self.D*t**4/4 - self.E/t)

    def _int_s(self, T):
        t = T/1000.0
        return (self.A*np.log(t) + self.B*t + self.C*t**2/2 + self.D*t**3/3 - self.E/(2*t**2))

    def h_ideal(self, T):
        return self._int_h(T) - self._int_h(T_REF)

    def s_ideal(self, T, P):
        return (self._int_s(T) - self._int_s(T_REF)) - R*np.log(P/1e5)

    # ---- PR departures (Poling Table 6-3) ------------------------------
    def _departure(self, T, P, phase):
        a = self.a_param(T)
        da = self.dadT(T)
        roots, A, B = self.z_roots(T, P)
        Z = roots[-1] if phase == "vapor" else roots[0]
        s = np.sqrt(u**2 - 4.0*w)
        log_arg = (2.0*Z + B*(u + s)) / (2.0*Z + B*(u - s))
        dh = R*T*(Z - 1.0) + (T*da - a)/(self.b*s)*np.log(log_arg)
        ds = R*np.log(Z - B) + da/(self.b*s)*np.log(log_arg)
        return dh, ds

    def _h_mass(self, T, P, phase):
        dh, _ = self._departure(T, P, phase)
        return (self.h_ideal(T) + dh) / MW

    def _s_mass(self, T, P, phase):
        _, ds = self._departure(T, P, phase)
        return (self.s_ideal(T, P) + ds) / MW

    # ---- anchored saturated properties (kJ/kg, kJ/kg/K) ----------------
    def hf(self, T):
        return self._h_mass(T, self.saturation_pressure(T), "liquid") + self.h_off

    def hg(self, T):
        return self._h_mass(T, self.saturation_pressure(T), "vapor") + self.h_off

    def sf(self, T):
        return self._s_mass(T, self.saturation_pressure(T), "liquid") + self.s_off

    def sg(self, T):
        return self._s_mass(T, self.saturation_pressure(T), "vapor") + self.s_off

    # ---- full dome for plotting ----------------------------------------
    def build_dome(self, n=120):
        T_grid = np.linspace(T_MIN, 0.99*self.Tc, n)
        return {
            "name": self.name,
            "T_C": T_grid - 273.15,
            "Pd": np.array([self.saturation_pressure(T)/1e5 for T in T_grid]),
            "hf": np.array([self.hf(T) for T in T_grid]),
            "hg": np.array([self.hg(T) for T in T_grid]),
            "sf": np.array([self.sf(T) for T in T_grid]),
            "sg": np.array([self.sg(T) for T in T_grid]),
        }

    
# =============================================================================
# Build every method
# =============================================================================

methods = [PRMethod(**row) for row in METHODS]
domes = [m.build_dome() for m in methods]

TC_REAL = 351.255  # real R-32 critical temperature (K), from CoolProp
# Same 0.99*Tc stopping rule as PRMethod.build_dome() above, applied to
# Helmholtz's own (real) Tc -- kept consistent so no series is allowed to
# approach its critical point more closely than any other. Used ONLY for
# the smooth plotted dome curve, not for the quantitative error table below
# (that uses the actual Linde points, L_T, same as every other method).
T_grid_helm = np.linspace(T_MIN, 0.99*TC_REAL, 120)
helm_dome = {
    "name": "Helmholtz",
    "T_C": T_grid_helm - 273.15,
    "Pd": np.array([CP.PropsSI("P", "T", T, "Q", 0, "R32")/1e5 for T in T_grid_helm]),
    "hf": np.array([CP.PropsSI("H", "T", T, "Q", 0, "R32")/1000.0 for T in T_grid_helm]),
    "hg": np.array([CP.PropsSI("H", "T", T, "Q", 1, "R32")/1000.0 for T in T_grid_helm]),
    "sf": np.array([CP.PropsSI("S", "T", T, "Q", 0, "R32")/1000.0 for T in T_grid_helm]),
    "sg": np.array([CP.PropsSI("S", "T", T, "Q", 1, "R32")/1000.0 for T in T_grid_helm]),
}

print(f"{'method':>6}{'omega':>10}{'Tc[C]':>10}")
for m in methods:
    print(f"{m.name:>6}{m.omega:>10.4f}{m.Tc - 273.15:>10.2f}")


# =============================================================================
# Error summary vs Linde (5 series: Linde reference + NIST/GCGP/SPGP/Helmholtz)
# =============================================================================
# Pressure error: percent (depends on Tc, Pc, omega only).
# Enthalpy/entropy error: mean absolute error (MAE). h and s use the IIR datum
# and pass near zero at the cold end, which makes percent error unstable there;
# MAE (kJ/kg, kJ/kg/K) is the meaningful metric and is what actually tests the
# Cp fits.
#
# All methods, including Helmholtz below, are evaluated at L_T -- the Linde
# datasheet's own 20 temperature points -- NOT the dense synthetic plotting
# grid used for build_dome()/helm_dome above. The dense grid is only for a
# smooth-looking curve; the actual accuracy comparison has to be done at the
# same points the reference data itself provides.

def report_errors(pr_methods):
    print("\nSaturation-property error vs Linde datasheet (selected points, -130..78 C)")
    print(f"{'series':>7}{'P_MAPE%':>10}{'hf_MAE':>9}{'hg_MAE':>9}{'sf_MAE':>10}{'sg_MAE':>10}")
    print(f"{'':>7}{'':>10}{'[kJ/kg]':>9}{'[kJ/kg]':>9}{'[kJ/kgK]':>10}{'[kJ/kgK]':>10}")
    print(f"{'Linde':>7}{'ref':>10}{'ref':>9}{'ref':>9}{'ref':>10}{'ref':>10}")
    for m in pr_methods:
        P = np.array([m.saturation_pressure(t)/1e5 for t in L_T])
        hf = np.array([m.hf(t) for t in L_T])
        hg = np.array([m.hg(t) for t in L_T])
        sf = np.array([m.sf(t) for t in L_T])
        sg = np.array([m.sg(t) for t in L_T])
        p_mape = np.mean(np.abs((P - L_P) / L_P)) * 100
        print(f"{m.name:>7}{p_mape:>10.2f}"
              f"{np.mean(np.abs(hf - L_HF)):>9.2f}"
              f"{np.mean(np.abs(hg - L_HG)):>9.2f}"
              f"{np.mean(np.abs(sf - L_SF)):>10.4f}"
              f"{np.mean(np.abs(sg - L_SG)):>10.4f}")

    # Helmholtz (CoolProp) -- same L_T points, same metrics, same format.
    # No h_off/s_off needed: CoolProp is already on the IIR datum.
    P_helm = np.array([CP.PropsSI("P", "T", t, "Q", 0, "R32")/1e5 for t in L_T])
    hf_helm = np.array([CP.PropsSI("H", "T", t, "Q", 0, "R32")/1000.0 for t in L_T])
    hg_helm = np.array([CP.PropsSI("H", "T", t, "Q", 1, "R32")/1000.0 for t in L_T])
    sf_helm = np.array([CP.PropsSI("S", "T", t, "Q", 0, "R32")/1000.0 for t in L_T])
    sg_helm = np.array([CP.PropsSI("S", "T", t, "Q", 1, "R32")/1000.0 for t in L_T])
    p_mape_helm = np.mean(np.abs((P_helm - L_P) / L_P)) * 100
    print(f"{'Helmholtz':>7}{p_mape_helm:>10.2f}"
          f"{np.mean(np.abs(hf_helm - L_HF)):>9.2f}"
          f"{np.mean(np.abs(hg_helm - L_HG)):>9.2f}"
          f"{np.mean(np.abs(sf_helm - L_SF)):>10.4f}"
          f"{np.mean(np.abs(sg_helm - L_SG)):>10.4f}")


report_errors(methods)


# =============================================================================
# Plot 1: p-h diagram (5 series: Linde + NIST + GCGP + SPGP + Helmholtz)
# =============================================================================

colors = ["tab:blue", "tab:green", "tab:red"]

plt.figure()
# Linde reference (datasheet) as black dashed with markers
plt.plot(np.concatenate([L_HF, L_HG[::-1]]),
         np.concatenate([L_P, L_P[::-1]]),
         "k--o", ms=3, lw=1, label="Linde (datasheet)")
for d, c in zip(domes, colors):
    h_dome = np.concatenate([d["hf"], d["hg"][::-1]])
    P_dome = np.concatenate([d["Pd"], d["Pd"][::-1]])
    plt.plot(h_dome, P_dome, color=c, label=d["name"])
plt.plot(np.concatenate([helm_dome["hf"], helm_dome["hg"][::-1]]),
         np.concatenate([helm_dome["Pd"], helm_dome["Pd"][::-1]]),
         "k-", lw=2, label="Helmholtz (CoolProp)")
plt.yscale("log")
plt.xlabel("h [kJ/kg]")
plt.ylabel("P [bar]")
plt.title("R-32 p-h diagram: methods vs Linde datasheet")
plt.legend()
plt.grid(True, which="both", ls=":")
plt.tight_layout()


# =============================================================================
# Plot 2: T-s diagram (5 series: Linde + NIST + GCGP + SPGP + Helmholtz)
# =============================================================================

plt.figure()
plt.plot(np.concatenate([L_SF, L_SG[::-1]]),
         np.concatenate([L_TC, L_TC[::-1]]),
         "k--o", ms=3, lw=1, label="Linde (datasheet)")
for d, c in zip(domes, colors):
    s_dome = np.concatenate([d["sf"], d["sg"][::-1]])
    T_dome = np.concatenate([d["T_C"], d["T_C"][::-1]])
    plt.plot(s_dome, T_dome, color=c, label=d["name"])
plt.plot(np.concatenate([helm_dome["sf"], helm_dome["sg"][::-1]]),
         np.concatenate([helm_dome["T_C"], helm_dome["T_C"][::-1]]),
         "k-", lw=2, label="Helmholtz (CoolProp)")
plt.xlabel("s [kJ/kg/K]")
plt.ylabel("T [deg C]")
plt.title("R-32 T-s diagram: methods vs Linde datasheet")
plt.legend()
plt.grid(True, ls=":")
plt.tight_layout()

plt.show()
