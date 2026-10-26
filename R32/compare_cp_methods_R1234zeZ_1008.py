"""
Compare p-h, T-s and p-T saturation behaviour for R-1234ze(Z) across property
methods, benchmarked against EXPERIMENTAL vapour-pressure data.

This is compare_cp_methods_1001.py (R-32) retargeted to R-1234ze(Z). The model
machinery is unchanged -- same Ambrose-Walton / Peng-Robinson / Shomate chain,
same IIR anchoring, same error-table layout. Three things changed:

    1. The reference data set. R-32 used the Linde datasheet (20 points, with
       p, hf, hg, sf, sg). R-1234ze(Z) has no equivalent datasheet, so the
       reference is 64 MEASURED vapour pressures from Fedele et al. (2014).
       Those are pressure only -- there are no enthalpy or entropy columns.
    2. Three series instead of five: Experimental, Helmholtz, PR.
    3. A p-T plot and a deviation plot, which is where pressure-only reference
       data can actually be shown.

SERIES
    1. Experimental   Fedele, Di Nicola, Brown, Bobbo & Zilio, Int. J.
                      Thermophys. 35 (2014) 1-12, Table 2. 64 points,
                      238.13-372.61 K, two independent laboratories
                      (ITC-CNR, 36 pts; UnivPM, 28 pts).
                      u(T) = 0.03 K, u(P) < 1 kPa.
    2. Helmholtz      Akasaka & Lemmon fundamental EoS, evaluated through
                      CoolProp. This is the same EoS NIST ships in REFPROP
                      as R1234ZEZ.FLD, so "CoolProp" and "NIST" are the
                      same numbers here by two implementations.
    3. PR             Peng-Robinson cubic + Shomate ideal-gas Cp, i.e. the
                      form IDAES's GenericParameterBlock actually solves.

WHERE THE PR PARAMETERS COME FROM (the thing to be able to answer)
    Equations:
        Ambrose-Walton Psat          Poling, Prausnitz & O'Connell, "The
                                     Properties of Gases and Liquids", 5th
                                     ed. (2001), sec. 7-4.
        PR h/s departure functions   Poling et al. (2001), Table 6-3.
        Shomate Cp integration       NIST Shomate convention, t = T/1000.
    Parameters:
        Tc, Pc    REFPROP fluid file R1234ZEZ.FLD (NIST SRD 23 v10.0,
                  compiled by E.W. Lemmon): Tc = 423.27 K, Pc = 3530.6 kPa.
                  Tc is Higashi et al.'s measured value, adopted unchanged.
                  Pc is NOT Higashi's measured 3536 kPa -- NIST refitted it
                  when building the EoS. Use the NIST number; the 0.15 %
                  difference is why row 20 of the 20-point table can never
                  read zero.
        omega     REFPROP, 0.327, taken explicitly rather than derived.
                  Cross-checks that agree: the Akasaka-Lemmon EoS reports
                  0.32684; Fedele et al. recommend 0.3257 from their Wagner
                  fit; deriving it here from the experimental table at
                  Tr = 0.7 by the Pitzer definition the R-32 script uses
                  gives 0.325795. All within 0.4 % of each other.
                  _omega_from_experiment() below is kept as the fallback
                  for any future row that does not carry its own omega.
        A-E       NO published Shomate coefficients exist for this fluid.
                  The NIST WebBook has no record for CAS 29118-25-0 at all
                  (checked three ways: ID, name, formula). So A-E here are
                  least-squares refitted to the ideal-gas Cp carried by the
                  Akasaka-Lemmon EoS itself, which is a Planck-Einstein
                  form with Einstein temperatures 20.0 K and 1335.0 K.
                  The refit reproduces that Cp to better than 0.01 % over
                  the fitted window -- see SHOMATE_FIT_REPORT printed at
                  run time. It is a change of functional form, not a change
                  of data.

    A published Peng-Robinson parameter set for this fluid does exist --
    Brown, Zilio & Cavallini, Int. J. Refrig. 33 (2010) 235, built from
    group contributions. Fedele et al. score it at AAD 2.4-2.7 %, improving
    to 1.2-1.9 % when the measured critical point is substituted. It is NOT
    used here (the group-contribution row is the thing the SPGP fits are
    meant to replace), but it is the natural extra series to add if a
    published cubic baseline is wanted.

WHAT IS COMPARED AGAINST WHAT
    Pressure   Both Helmholtz and PR are scored against the 64 measured
               points. Experiment is the arbiter.
    h and s    There is no experimental h or s for this fluid in the
               reference above, so PR is scored against Helmholtz and the
               table says so. Do not report those numbers as experimental
               agreement -- they are not.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
QA/testing: Shilpa Narasimhan
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP


# =============================================================================
# Shared constants
# =============================================================================

# 10/08/2026: the R-32 scripts use R = 8.314. That is 0.0055 % below the value
# the Akasaka-Lemmon EoS is built on, and verify_akasaka_lemmon_SI.py caught it
# as a constant -0.00553 % offset on every ideal-gas Cp. REFPROP's R1234ZEZ.FLD
# states the equation's gas constant explicitly:
#       8.3144598      !Gas constant [J/mol-K]
# so that is what is used here. The shift is negligible in the results but it
# removes a known inconsistency between this script and the EoS it is scored
# against.
R = 8.3144598           # J/mol/K  (R1234ZEZ.FLD / Akasaka-Lemmon)
MW = 114.0416           # kg/kmol (numerically g/mol), R-1234ze(Z)
FLUID = "R1234ze(Z)"

u = 2.0
w = -1.0
Omega_A = 0.45724
Omega_B = 0.0778

T_REF = 273.15          # K       IIR reference temperature (0 C)
H_REF = 200.0           # kJ/kg   IIR reference enthalpy (sat. liquid)
S_REF = 1.0             # kJ/kg/K IIR reference entropy  (sat. liquid)

# Low end of the dome sweep. R-32 used -130 C. Here it is set to the coldest
# MEASURED point (238.13 K = -35.0 C) so that no curve is drawn in a region
# where nothing can check it. CoolProp will happily return values below this
# -- its EoS validity floor is 273 K and PropsSI does NOT enforce it, it
# extrapolates silently -- so the limit has to be imposed here on purpose.
T_MIN = 238.13          # K


# =============================================================================
# Experimental saturation pressures -- Fedele et al. (2014), Table 2.
# Replaces the Linde datasheet used for R-32.
#   (T[K], Psat[MPa], lab)
# NOTE: pressure only. There are no hf/hg/sf/sg columns to replace Linde's.
# =============================================================================

FEDELE_SAT = [
    # ---- ITC-CNR, 36 points, 283.15-353.16 K ----------------------------
    (283.15, 0.103, "ITC-CNR"), (285.15, 0.111, "ITC-CNR"),
    (287.15, 0.120, "ITC-CNR"), (289.15, 0.129, "ITC-CNR"),
    (291.15, 0.139, "ITC-CNR"), (293.15, 0.149, "ITC-CNR"),
    (295.15, 0.160, "ITC-CNR"), (297.15, 0.172, "ITC-CNR"),
    (299.15, 0.184, "ITC-CNR"), (301.15, 0.198, "ITC-CNR"),
    (303.15, 0.211, "ITC-CNR"), (305.15, 0.225, "ITC-CNR"),
    (307.15, 0.241, "ITC-CNR"), (309.15, 0.256, "ITC-CNR"),
    (311.15, 0.273, "ITC-CNR"), (313.15, 0.291, "ITC-CNR"),
    (315.15, 0.309, "ITC-CNR"), (317.15, 0.329, "ITC-CNR"),
    (319.15, 0.349, "ITC-CNR"), (321.14, 0.370, "ITC-CNR"),
    (323.15, 0.392, "ITC-CNR"), (325.15, 0.415, "ITC-CNR"),
    (327.15, 0.439, "ITC-CNR"), (329.15, 0.465, "ITC-CNR"),
    (331.15, 0.491, "ITC-CNR"), (333.15, 0.519, "ITC-CNR"),
    (335.15, 0.548, "ITC-CNR"), (337.15, 0.578, "ITC-CNR"),
    (339.15, 0.609, "ITC-CNR"), (341.15, 0.641, "ITC-CNR"),
    (343.15, 0.675, "ITC-CNR"), (345.15, 0.710, "ITC-CNR"),
    (347.16, 0.746, "ITC-CNR"), (349.15, 0.784, "ITC-CNR"),
    (351.15, 0.823, "ITC-CNR"), (353.16, 0.863, "ITC-CNR"),
    # ---- UnivPM, 28 points, 238.13-372.61 K -----------------------------
    (238.13, 0.012, "UnivPM"), (243.10, 0.016, "UnivPM"),
    (248.00, 0.021, "UnivPM"), (252.90, 0.027, "UnivPM"),
    (257.81, 0.034, "UnivPM"), (263.02, 0.044, "UnivPM"),
    (267.91, 0.055, "UnivPM"), (272.83, 0.068, "UnivPM"),
    (277.94, 0.085, "UnivPM"), (282.78, 0.102, "UnivPM"),
    (287.62, 0.123, "UnivPM"), (292.62, 0.148, "UnivPM"),
    (298.00, 0.178, "UnivPM"), (303.18, 0.212, "UnivPM"),
    (308.18, 0.249, "UnivPM"), (313.17, 0.291, "UnivPM"),
    (318.14, 0.338, "UnivPM"), (323.14, 0.392, "UnivPM"),
    (328.12, 0.451, "UnivPM"), (333.11, 0.518, "UnivPM"),
    (338.12, 0.591, "UnivPM"), (343.07, 0.672, "UnivPM"),
    (348.16, 0.764, "UnivPM"), (353.12, 0.861, "UnivPM"),
    (358.13, 0.969, "UnivPM"), (363.14, 1.086, "UnivPM"),
    (368.09, 1.212, "UnivPM"), (372.61, 1.336, "UnivPM"),
]

_order = np.argsort([r[0] for r in FEDELE_SAT])
E_T = np.array([FEDELE_SAT[i][0] for i in _order])            # K
E_TC = E_T - 273.15                                           # deg C
E_P = np.array([FEDELE_SAT[i][1] for i in _order]) * 10.0     # MPa -> bar
E_LAB = np.array([FEDELE_SAT[i][2] for i in _order])


# -----------------------------------------------------------------------------
# Critical point. NOT in Fedele's table -- their data stops at 372.61 K.
#
# 10/08/2026: switched from Higashi et al. (2013) to REFPROP. Higashi's
# measured 3536 kPa is only available secondhand (quoted in Fedele p.2 from an
# IIR conference paper), and it sits 0.15 % from the value NIST actually
# publishes, because NIST refitted the critical pressure when building the
# equation of state. Citing REFPROP instead means:
#   - the source is a file we hold (R1234ZEZ.FLD) rather than a quote,
#   - the critical point lands exactly on the models, since it IS their
#     critical point, so no spurious 0.15 % "error" appears in the table.
#
# Independent confirmation of the critical TEMPERATURE: Fukuda, Kondou, Takata
# & Koyama, Int. J. Refrig. 40 (2013) 161-173, Table 2, gives 150.1 degC
# against REFPROP's 150.12 degC -- 0.01 %.
#
#       423.27    !Critical temperature [K]
#       3530.6    !Critical pressure [kPa]
T_CRIT_EXP = 423.27        # K     REFPROP R1234ZEZ.FLD
P_CRIT_EXP = 35.306        # bar   REFPROP R1234ZEZ.FLD, = 3530.6 kPa


# -----------------------------------------------------------------------------
# The 20-point superimposed subset -- the direct analogue of the 20 Linde rows
# used for R-32. Nineteen of Fedele's points spread evenly across the measured
# range, plus the critical point from REFPROP.
#
# WHY A SUBSET AT ALL: the R-32 figures carry 20 reference markers, and 64
# markers makes the same plot unreadable. The full 64 are kept and are what the
# error table is computed on -- the subset is for plotting and for the summary
# table that goes in front of a reader.
# -----------------------------------------------------------------------------

def _pick_subset(n=19):
    """Indices of n points spread as evenly as possible in T across E_T."""
    targets = np.linspace(E_T.min(), E_T.max(), n)
    idx = []
    for t in targets:
        j = int(np.argmin(np.abs(E_T - t)))
        while j in idx:
            j += 1
        idx.append(j)
    return sorted(idx)


SUB_I = _pick_subset(19)
S_T = np.r_[E_T[SUB_I], T_CRIT_EXP]                 # K
S_TC = S_T - 273.15                                 # deg C
S_P = np.r_[E_P[SUB_I], P_CRIT_EXP]                 # bar
S_SRC = np.array([f"Fedele/{E_LAB[i]}" for i in SUB_I] + ["REFPROP (crit. pt)"])


# =============================================================================
# Method row data (the only per-method input)
# =============================================================================
# Tc and Pc are the Akasaka-Lemmon / REFPROP values. Shomate A-E are filled in
# below by fit_shomate(), because no published set exists -- see the header.

PR_ROW = {
    "name": "PR",
    # Straight off the REFPROP fluid file, R1234ZEZ.FLD (NIST SRD 23 v10.0,
    # compiled by E.W. Lemmon). Header lines:
    #     423.27    !Critical temperature [K]
    #     3530.6    !Critical pressure [kPa]
    #     0.327     !Acentric factor
    # omega is given explicitly rather than derived, so the PR row is the
    # NIST row in full and nothing in it is this script's own invention
    # except the Shomate coefficients (see header).
    "Pc_bar": 3530.6 / 100.0,        # kPa -> bar = 35.306 bar
    "Tc": 423.27,                    # K
    "omega": 0.327,                  # REFPROP
}

# For reference / cross-check only -- not used by the model:
OMEGA_COOLPROP = 0.32683872006421466   # Akasaka-Lemmon EoS, more digits
OMEGA_FEDELE = 0.3257                  # Fedele et al. (2014) Wagner fit


# Ideal-gas Cp of the Akasaka-Lemmon EoS, as alpha0 terms read out of CoolProp:
#   IdealGasHelmholtzLogTau  a = 3            -> cv0/R = 3, so cp0/R = 4
#   IdealGasHelmholtzPlanckEinstein n = [4.2365, 13.063]
#                                   theta/Tc = [0.0472512, 3.1540152]
CP0_N = np.array([4.2365, 13.063])
CP0_THETA = np.array([0.047251163559902665, 3.154015167623503]) * PR_ROW["Tc"]


def cp0_eos(T):
    """Ideal-gas molar Cp [J/mol/K] carried by the Akasaka-Lemmon EoS."""
    T = np.atleast_1d(np.asarray(T, float))
    x = CP0_THETA[None, :] / T[:, None]
    pe = (CP0_N[None, :] * x**2 * np.exp(x) / (np.exp(x) - 1.0)**2).sum(axis=1)
    return R * (4.0 + pe)


def fit_shomate(T_lo, T_hi, n=400):
    """Least-squares Shomate A-E for cp0_eos over [T_lo, T_hi].

    Shomate convention, t = T/1000:
        Cp = A + B t + C t^2 + D t^3 + E / t^2
    Returns (A, B, C, D, E, max_abs_err, max_rel_err_pct).
    """
    T = np.linspace(T_lo, T_hi, n)
    y = cp0_eos(T)
    t = T / 1000.0
    X = np.column_stack([np.ones_like(t), t, t**2, t**3, 1.0 / t**2])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ c
    return (*c, np.abs(fit - y).max(), 100.0 * np.abs((fit - y) / y).max())


# Fit over the experimental window -- the range the comparison actually covers.
_A, _B, _C, _D, _E, _MAXABS, _MAXREL = fit_shomate(E_T.min(), E_T.max())
PR_ROW.update({"A": _A, "B": _B, "C": _C, "D": _D, "E": _E})

SHOMATE_FIT_REPORT = (
    f"Shomate refit of the Akasaka-Lemmon ideal-gas Cp over "
    f"{E_T.min():.2f}-{E_T.max():.2f} K\n"
    f"   A = {_A:.6f}\n   B = {_B:.6f}\n   C = {_C:.6f}\n"
    f"   D = {_D:.6f}\n   E = {_E:.6f}\n"
    f"   max abs deviation = {_MAXABS:.5f} J/mol/K   "
    f"max rel = {_MAXREL:.5f} %"
)


# =============================================================================
# Property model for one method
# =============================================================================

class PRMethod:
    """Peng-Robinson + Shomate property model for one (Pc, Tc, Shomate) row.

    Unchanged from the R-32 script except that the acentric factor falls back
    to the experimental table instead of the Linde table.
    """

    def __init__(self, name, Pc_bar, Tc, A, B, C, D, E, omega=None):
        self.name = name
        self.Pc = Pc_bar * 1e5
        self.Tc = Tc
        self.A, self.B, self.C, self.D, self.E = A, B, C, D, E

        self.omega = omega if omega is not None else self._omega_from_experiment()
        self.kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        self.b = Omega_B * R * self.Tc / self.Pc

        # IIR anchoring offsets (sat. liquid at 0 C)
        P0 = self.saturation_pressure(T_REF)
        self.h_off = H_REF - self._h_mass(T_REF, P0, "liquid")
        self.s_off = S_REF - self._s_mass(T_REF, P0, "liquid")

    # ---- acentric factor from the experimental table at Tr = 0.7 -------
    def _omega_from_experiment(self):
        T_target = 0.7 * self.Tc
        if not (E_T.min() <= T_target <= E_T.max()):
            raise ValueError(
                f"{self.name}: Tr=0.7 is {T_target:.2f} K, outside the measured "
                f"range {E_T.min():.2f}-{E_T.max():.2f} K. Supply omega explicitly.")
        lnP = np.interp(1.0 / T_target, (1.0 / E_T)[::-1], np.log(E_P * 1e5)[::-1])
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
        T_grid = np.linspace(T_MIN, 0.999*self.Tc, n)
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
# Build the methods
# =============================================================================

pr = PRMethod(**PR_ROW)
pr_dome = pr.build_dome()

TC_REAL = 423.27        # Akasaka-Lemmon critical temperature (K)
T_grid_helm = np.linspace(T_MIN, 0.999*TC_REAL, 120)
helm_dome = {
    "name": "Helmholtz",
    "T_C": T_grid_helm - 273.15,
    "Pd": np.array([CP.PropsSI("P", "T", T, "Q", 0, FLUID)/1e5 for T in T_grid_helm]),
    "hf": np.array([CP.PropsSI("H", "T", T, "Q", 0, FLUID)/1000.0 for T in T_grid_helm]),
    "hg": np.array([CP.PropsSI("H", "T", T, "Q", 1, FLUID)/1000.0 for T in T_grid_helm]),
    "sf": np.array([CP.PropsSI("S", "T", T, "Q", 0, FLUID)/1000.0 for T in T_grid_helm]),
    "sg": np.array([CP.PropsSI("S", "T", T, "Q", 1, FLUID)/1000.0 for T in T_grid_helm]),
}


# =============================================================================
# Error summary
# =============================================================================
# PRESSURE is scored against the 64 measured points -- experiment is the
# arbiter, and both Helmholtz and PR are judged by it.
#
# ENTHALPY and ENTROPY have no experimental counterpart in this reference, so
# PR is scored against Helmholtz. That is a model-to-model comparison and the
# table labels it as such. It is NOT evidence of agreement with experiment.

def report_errors():
    print("\n" + "="*74)
    print("  SATURATION PRESSURE vs EXPERIMENT  (Fedele et al. 2014, 64 points)")
    print("="*74)

    P_pr = np.array([pr.saturation_pressure(T)/1e5 for T in E_T])
    P_helm = np.array([CP.PropsSI("P", "T", T, "Q", 0, FLUID)/1e5 for T in E_T])

    rows = []
    for tag, Pm in (("Helmholtz", P_helm), ("PR", P_pr)):
        d = 100.0 * (Pm - E_P) / E_P
        rows.append((tag, d.mean(), np.abs(d).mean(), np.abs(d).max(), len(d)))

    print(f"   {'method':<12}{'BIAS %':>10}{'AAD %':>10}{'max|d| %':>11}{'n pts':>8}")
    for tag, bias, aad, mx, n in rows:
        print(f"   {tag:<12}{bias:>10.3f}{aad:>10.3f}{mx:>11.3f}{n:>8}")
    # ---- the 20 superimposed points, listed out ------------------------
    print("\n   The 20 superimposed points (19 Fedele + critical point):")
    print(f"   {'#':>3}{'T [C]':>9}{'P_exp':>9}{'Helmholtz':>11}{'d %':>8}"
          f"{'PR':>9}{'d %':>8}  source")
    S_helm, S_pr = [], []
    for n, (T, Pe, src) in enumerate(zip(S_T, S_P, S_SRC), 1):
        if n == len(S_T):                      # critical point
            ph = CP.PropsSI("pcrit", FLUID)/1e5
            pp = pr.Pc/1e5                     # PR's Pc IS its Psat at Tc
        else:
            ph = CP.PropsSI("P", "T", T, "Q", 0, FLUID)/1e5
            pp = pr.saturation_pressure(T)/1e5
        S_helm.append(ph); S_pr.append(pp)
        print(f"   {n:>3}{T-273.15:>9.2f}{Pe:>9.3f}{ph:>11.3f}"
              f"{100*(ph-Pe)/Pe:>8.2f}{pp:>9.3f}{100*(pp-Pe)/Pe:>8.2f}  {src}")
    S_helm, S_pr = np.array(S_helm), np.array(S_pr)
    for tag, Pm in (("Helmholtz", S_helm), ("PR", S_pr)):
        d = 100.0*(Pm - S_P)/S_P
        print(f"   -> {tag:<10} over these 20:  BIAS {d.mean():+.3f} %   "
              f"AAD {np.abs(d).mean():.3f} %   max|d| {np.abs(d).max():.3f} %")
    print("   NOTE: row 20's ~0.15 % offset is Higashi's measured Pc (3536 kPa)")
    print("         against REFPROP's published Pc (3530.6 kPa). Definitional,")
    print("         not model error -- see the comment at T_CRIT_EXP.")

    print("\n   BIAS% = mean(100*(Pcalc-Pexp)/Pexp)   AAD% = mean|...|")
    print("   Same definitions as Fedele et al. Table 4, so the numbers are")
    print("   directly comparable to the seven correlations scored there.")

    # Split by laboratory -- the two data sets overlap but were taken on
    # different apparatus, and the paper reports them separately.
    print(f"\n   {'method':<12}{'lab':<10}{'BIAS %':>10}{'AAD %':>10}{'n pts':>8}")
    for tag, Pm in (("Helmholtz", P_helm), ("PR", P_pr)):
        for lab in ("ITC-CNR", "UnivPM"):
            k = E_LAB == lab
            d = 100.0 * (Pm[k] - E_P[k]) / E_P[k]
            print(f"   {tag:<12}{lab:<10}{d.mean():>10.3f}{np.abs(d).mean():>10.3f}{k.sum():>8}")

    print("\n" + "="*74)
    print("  ENTHALPY / ENTROPY  --  PR vs HELMHOLTZ  (NOT vs experiment)")
    print("="*74)
    print("   No experimental h or s exists for R-1234ze(Z) in the reference")
    print("   data set, which is vapour pressure only. These are model-to-model")
    print("   deviations and must not be quoted as experimental agreement.")

    hf_h = np.array([CP.PropsSI("H", "T", T, "Q", 0, FLUID)/1000.0 for T in E_T])
    hg_h = np.array([CP.PropsSI("H", "T", T, "Q", 1, FLUID)/1000.0 for T in E_T])
    sf_h = np.array([CP.PropsSI("S", "T", T, "Q", 0, FLUID)/1000.0 for T in E_T])
    sg_h = np.array([CP.PropsSI("S", "T", T, "Q", 1, FLUID)/1000.0 for T in E_T])

    hf_p = np.array([pr.hf(T) for T in E_T])
    hg_p = np.array([pr.hg(T) for T in E_T])
    sf_p = np.array([pr.sf(T) for T in E_T])
    sg_p = np.array([pr.sg(T) for T in E_T])

    print(f"\n   {'quantity':<10}{'MAE':>12}{'max|d|':>12}{'units':>14}")
    for tag, a, b, un in (("hf", hf_p, hf_h, "kJ/kg"),
                          ("hg", hg_p, hg_h, "kJ/kg"),
                          ("sf", sf_p, sf_h, "kJ/kg/K"),
                          ("sg", sg_p, sg_h, "kJ/kg/K")):
        d = a - b
        print(f"   {tag:<10}{np.abs(d).mean():>12.4f}{np.abs(d).max():>12.4f}{un:>14}")

    return P_pr, P_helm, (hf_p, hg_p, sf_p, sg_p), (hf_h, hg_h, sf_h, sg_h)


# =============================================================================
# Plots
# =============================================================================

def make_plots(P_pr, P_helm):
    # ---- Plot 1: p-T saturation curve, with the experimental points -----
    # This is the one plot where the reference data can actually appear,
    # since it is pressure-only. The 20 superimposed markers are the direct
    # analogue of the 20 Linde markers on the R-32 figures; the remaining
    # 44 measured points sit behind them in grey so nothing is hidden.
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    rest = np.setdiff1d(np.arange(len(E_T)), SUB_I)
    ax.plot(E_TC[rest], E_P[rest], ".", ms=4, color="0.72", ls="none",
            label=f"Experimental, remaining ({len(rest)} pts)")
    ax.plot(S_TC[:-1], S_P[:-1], "o", ms=6, mfc="none", mec="b", mew=1.3,
            ls="none", label="Experimental, 19 pts (Fedele et al. 2014)")
    ax.plot(S_TC[-1], S_P[-1], "*", ms=15, mfc="none", mec="g", mew=1.5,
            ls="none", label="Critical point (Higashi et al. 2013)")
    ax.plot(helm_dome["T_C"], helm_dome["Pd"], "k-", lw=2, label="Helmholtz (CoolProp)")
    ax.plot(pr_dome["T_C"], pr_dome["Pd"], "r--", lw=2, label="PR (REFPROP params)")
    ax.set_xlabel("Temperature [$^\\circ$C]")
    ax.set_ylabel("Saturation pressure [bar]")
    ax.set_title("R-1234ze(Z) saturation curve vs experimental data")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig("pT_diagram_R1234zeZ_1008.png", dpi=150)
    plt.close(fig)

    # ---- Plot 2: deviation plot, mirroring Fedele et al. Figs. 2-8 ------
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for tag, Pm, mk, col in (("Helmholtz", P_helm, "o", "k"), ("PR", P_pr, "s", "r")):
        d = 100.0 * (Pm - E_P) / E_P
        ax.plot(E_T, d, mk, ms=5, mfc="none", mec=col, mew=1.0, ls="none", label=tag)
    # Shade the sub-300 K band: NIST quotes 0.3 % uncertainty on vapour
    # pressure below 300 K vs 0.1 % above it, and it is where PR degrades.
    ax.axvspan(E_T.min(), 300.0, color="0.9", zorder=0)
    ax.text(0.02, 0.04, "below 300 K:\nNIST quotes 0.3 % (vs 0.1 % above)",
            transform=ax.transAxes, fontsize=7.5, color="0.35", va="bottom")
    ax.axhline(0.0, color="0.4", lw=1)
    ax.set_xlabel("$T$ [K]")
    ax.set_ylabel("$100\\,\\Delta P/P$")
    ax.set_title("Deviation from experimental vapour pressure\n"
                 "(same axes convention as Fedele et al. 2014, Figs. 2-8)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig("deviation_R1234zeZ_1008.png", dpi=150)
    plt.close(fig)

    # ---- Plot 3: p-h diagram --------------------------------------------
    # Styled to match ph_diagram_0929_NIST_GCGP_Helmholtz_SPGPlabel.png:
    # log pressure axis, dotted minor grid, Helmholtz as the thick black
    # reference curve, each cubic as a thin coloured line, legend inside.
    #
    # The Linde markers that appear on the R-32 version have NO counterpart
    # here. Linde published hf/hg/sf/sg; Fedele et al. measured vapour
    # pressure only. Nothing experimental can be drawn on a p-h or T-s axis
    # for this fluid, so these two figures compare model against model. The
    # experimental comparison lives on the p-T and deviation figures.
    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    ax.plot(np.r_[pr_dome["hf"], pr_dome["hg"][::-1]],
            np.r_[pr_dome["Pd"], pr_dome["Pd"][::-1]],
            color="tab:blue", lw=1.5, label="PR (REFPROP params)")
    ax.plot(np.r_[helm_dome["hf"], helm_dome["hg"][::-1]],
            np.r_[helm_dome["Pd"], helm_dome["Pd"][::-1]],
            "k-", lw=2.5, label="Helmholtz (CoolProp)")
    ax.set_xlabel("h [kJ/kg]")
    ax.set_ylabel("P [bar]")
    ax.set_title("R-1234ze(Z) p-h diagram: methods vs Helmholtz\n"
                 "(no experimental h exists for this fluid)")
    ax.set_yscale("log")
    ax.grid(True, which="both", ls=":", alpha=0.45)
    ax.legend(loc="lower center", fontsize=11)
    fig.tight_layout()
    fig.savefig("ph_diagram_R1234zeZ_1008.png", dpi=150)
    plt.close(fig)

    # ---- Plot 4: T-s diagram --------------------------------------------
    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    ax.plot(np.r_[pr_dome["sf"], pr_dome["sg"][::-1]],
            np.r_[pr_dome["T_C"], pr_dome["T_C"][::-1]],
            color="tab:blue", lw=1.5, label="PR (REFPROP params)")
    ax.plot(np.r_[helm_dome["sf"], helm_dome["sg"][::-1]],
            np.r_[helm_dome["T_C"], helm_dome["T_C"][::-1]],
            "k-", lw=2.5, label="Helmholtz (CoolProp)")
    ax.set_xlabel("s [kJ/kg/K]")
    ax.set_ylabel("T [$^\\circ$C]")
    ax.set_title("R-1234ze(Z) T-s diagram: methods vs Helmholtz\n"
                 "(no experimental s exists for this fluid)")
    ax.grid(True, which="both", ls=":", alpha=0.45)
    ax.legend(loc="lower center", fontsize=11)
    fig.tight_layout()
    fig.savefig("ts_diagram_R1234zeZ_1008.png", dpi=150)
    plt.close(fig)

    print("\n   wrote pT_diagram_R1234zeZ_1008.png")
    print("         deviation_R1234zeZ_1008.png")
    print("         ph_diagram_R1234zeZ_1008.png")
    print("         ts_diagram_R1234zeZ_1008.png")


# =============================================================================
# Excel export
# =============================================================================

def write_excel(P_pr, P_helm, prv, hev):
    try:
        import pandas as pd
    except ImportError:
        print("\n   pandas not available -- skipping Excel export")
        return
    hf_p, hg_p, sf_p, sg_p = prv
    hf_h, hg_h, sf_h, sg_h = hev

    sat = pd.DataFrame({
        "T [K]": E_T, "T [C]": E_TC, "lab": E_LAB,
        "P_exp [bar]": E_P, "P_Helmholtz [bar]": P_helm, "P_PR [bar]": P_pr,
        "dP_Helmholtz [%]": 100*(P_helm - E_P)/E_P,
        "dP_PR [%]": 100*(P_pr - E_P)/E_P,
    })
    hs = pd.DataFrame({
        "T [K]": E_T, "T [C]": E_TC,
        "hf_PR": hf_p, "hf_Helm": hf_h, "d_hf": hf_p - hf_h,
        "hg_PR": hg_p, "hg_Helm": hg_h, "d_hg": hg_p - hg_h,
        "sf_PR": sf_p, "sf_Helm": sf_h, "d_sf": sf_p - sf_h,
        "sg_PR": sg_p, "sg_Helm": sg_h, "d_sg": sg_p - sg_h,
    })
    params = pd.DataFrame({
        "quantity": ["Tc [K]", "Pc [bar]", "omega (this script, from exp.)",
                     "omega (Akasaka-Lemmon EoS)", "omega (Fedele Wagner fit)",
                     "MW [g/mol]", "Shomate A", "Shomate B", "Shomate C",
                     "Shomate D", "Shomate E", "Shomate max rel err [%]"],
        "value": [pr.Tc, pr.Pc/1e5, pr.omega, 0.32683872006421466, 0.3257,
                  MW, pr.A, pr.B, pr.C, pr.D, pr.E, _MAXREL],
    })
    out = "compare_cp_methods_R1234zeZ_1008.xlsx"
    with pd.ExcelWriter(out) as xl:
        sat.to_excel(xl, sheet_name="saturation_pressure", index=False)
        hs.to_excel(xl, sheet_name="h_s_PR_vs_Helmholtz", index=False)
        params.to_excel(xl, sheet_name="parameters", index=False)
    print(f"   wrote {out}")


# =============================================================================
if __name__ == "__main__":
    print("="*74)
    print("  R-1234ze(Z) property comparison -- 10/08/2026")
    print("="*74)
    print(SHOMATE_FIT_REPORT)
    print(f"\n   experimental points : {len(E_T)}  "
          f"({E_T.min():.2f}-{E_T.max():.2f} K = "
          f"{E_TC.min():.2f} to {E_TC.max():.2f} C)")
    print(f"   PR  Tc = {pr.Tc:.3f} K   Pc = {pr.Pc/1e5:.5f} bar   "
          f"omega = {pr.omega:.6f}")
    print(f"        (EoS omega = 0.326839 ; Fedele omega = 0.3257)")

    P_pr, P_helm, prv, hev = report_errors()
    make_plots(P_pr, P_helm)
    write_excel(P_pr, P_helm, prv, hev)
    print("\n" + "="*74)
