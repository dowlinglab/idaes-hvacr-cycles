"""
Compare p-h and T-s diagrams for R-32 across several Cp / critical-property
methods, benchmarked against the Linde datasheet.

Each method supplies only its ROW DATA: Pc [bar], Tc [K], and Shomate A..E.
Shared/fixed: R-32 molar mass, gas constant, the Linde saturation table
(pressure column used to back out each method's acentric factor at Tr = 0.7;
enthalpy/entropy columns used as the validation reference), and the IIR
reference state (sat. liquid at 0 C: h = 200 kJ/kg, s = 1.0 kJ/kg/K).

Model per method:
    1. Ambrose-Walton saturation pressure (Poling et al., 2001, sec. 7-4).
    2. Peng-Robinson cubic EOS departure h and s (Poling Table 6-3).
    3. Shomate ideal-gas Cp integration (t = T/1000).

--- Before 08/26/2026 ---
Nine data series were compared: Linde datasheet (reference, 20 points),
NIST (NIST Shomate + NIST Pc/Tc), GCGP (GCGP Shomate + GCGP Pc/Tc), SPGP
(SPGP Shomate + SPGP predicted Tc), and five SPGP mix-and-match variants
(SPGP_NIST, SPGP_Fathya, SPGP_Shomate_GCGP, SPGP_NIST_crit_prps,
SPGP_Fathya_crit_prps) pairing SPGP's Shomate fit against other methods'
Tc/Pc, built to probe whether SPGP's poor fit traced to its Shomate curve
or its critical properties.

--- 08/26/2026 ---
SPGP and its five mix-and-match variants are commented out (not deleted)
in METHODS, replaced by two new methods from collaborator-supplied Shomate
fits: "first_principle" and "gcn" (from spgp_r32.xlsx's a-e coefficients,
rescaled from raw-T form to this file's t = T/1000 Shomate convention).

Fixed, was crashing: "gcn"'s given Tc (329.12 K / 55.97 C) falls below the
top of Linde's own comparison range (78 C), so Ambrose-Walton's tau =
1 - T/Tc went negative there, and its fractional powers returned NaN,
breaking the PR cubic-root solver inside report_errors(). Fixed by
masking each method's Linde comparison points to T < 0.99*Tc and
reporting how many of the 20 points were actually used (see "npts"
column) -- mirrors the same "flag partial coverage" approach already
used for SPGP's limited T-range in
compare_cp_methods_numerical_integration.py.

omega for first_principle/gcn no longer comes from _omega_from_linde().
Both carry a "psat_pa" row entry instead of a hardcoded omega -- the
Colon group's own "pvap" spreadsheet column (mislabeled mmHg, actually
Pa), taken as Psat at Tr = 0.7*Tc, the standard Pitzer input. The new
omega_from_psat(name, Pc_bar, psat_pa) function computes AND PRINTS
omega live from that value (omega = -1 - log10(Psat/Pc)) rather than
hardcoding the result; NIST/GCGP (no "psat_pa" key) are unaffected and
still resolve via _omega_from_linde() as before. Current values:
first_principle = 0.191312, gcn = 0.667570 (both replace the earlier,
now-obsolete Linde-derived placeholders of -0.2385/0.6234).

Even with both fixes in, first_principle/gcn still show large errors vs
Linde (first_principle: 65.99% P_MAPE; gcn: 46.72% P_MAPE, hf_MAE alone
192 kJ/kg) -- this is no longer a units or omega-computation bug, it
points to the underlying Tc/Pc/Shomate values themselves still being
questionable, the same category of problem the old SPGP had. Not yet
resolved; awaiting Colon-group confirmation.

--- 09/03/2026 ---
first_principle updated to the Colon group's post-thermoreconciliation
values from spgp_r32 (4).xlsx: Pc=49.62 bar, Tc=392.094 K, pvap=4579.217
(confirmed today as an mmHg reading -- gives a physical omega=-0.090036,
unlike the Pa reading which gives omega=2.03), Shomate rescaled A-E. This
is a real correctness improvement, independently verified: this same
first_principle data already passed the IDAES-vs-vanilla-PR cubic-EOS
gate in phase_1_cubic_eos_validation_refstate_0903.py (|dZ|=4.10e-07)
before being carried over here. Note this file's omega_from_psat() takes
psat_pa directly (Pa), unlike _0903.py's mmHg-native version -- the mmHg
value is converted to Pa inline in the METHODS row below.

--- 09/03/2026 (later same day) ---
gcn ALSO updated to (4).xlsx thermoreconciled values, per Shilpa's
explicit direction that today's computations use ONLY the
thermoreconciled numbers, not a mix of old/new: Pc=63.572 bar,
Tc=328.614 K, pvap=36566.042. omega uses the SAME mmHg-reading
convention as first_principle (same spreadsheet/column -- the unit
convention belongs to the column, not the method), giving omega=-0.88.
This is well outside any physically normal range (real fluids run
roughly -0.2 to 0.9), but is used anyway for consistency with
first_principle rather than switching conventions per-method. This is
corroborating evidence of gcn's already-flagged data problem (Colon
group's own "4.5% feasible within CI" note), not a reason to read the
units differently. Treat gcn's curve here as provisional/exploratory,
not a validated result.

--- 09/09/2026 ---
New spreadsheet delivery (spgp_r32 (5).xlsx, data section itself marked
"updated on sep 6"). first_principle and gcn's "after thermoreconciliation"
values were revised (superseding 09/03's (4).xlsx numbers), and two new
sections appeared -- "multi-output raw" and "multi after recon" -- added
here as first_principle_multi / gcn_multi. "Multi-output" not yet defined
by the collaborators; question sent, reply pending.

first_principle: Pc=60.921 bar (was 49.62), Tc=355.615 K (was 392.094),
pvap=7409.643 (was 4579.217). Shomate a-e UNCHANGED from 09/03. omega
(mmHg reading): -0.209932 -- physically reasonable, looks better than
09/03's value.

gcn: Pc=62.836 bar (was 63.572), Tc=339.914 K (was 328.614),
pvap=12585.513 (was 36566.042). Shomate a-e ALSO changed this time.
omega (mmHg): -0.426564 -- still negative but a real improvement over
09/03's -0.88.

first_principle_multi (new): Pc=64.965 bar, Tc=454.253 K, pvap=4243.056.
omega (mmHg): 0.060099 -- comfortably physical, the best-behaved omega
of any first_principle/gcn variant so far.

gcn_multi (new): Pc=35.704 bar, Tc=385.664 K, pvap=149.061. omega (mmHg):
1.254451 -- unphysical, same category of problem as gcn's original
09/03 value. Treat as provisional/exploratory only.

Also see phase_1_cubic_eos_validation_refstate_0909.py's own 09/06 note
for the still-unconfirmed Shomate T-vs-t=T/1000 scaling question
(tested empirically there; question sent to the collaborator, reply
pending) -- that assumption applies to every A-E value in this file too.

--- 09/10/2026 ---
New spreadsheet delivery (spgp_r32 (7).xlsx). Confirmed by direct cell
comparison that (7).xlsx restates the 09/06 numbers unchanged in one
block and adds a genuinely revised fit in a second block ("additional
vapor information") -- wired in the revised block here, matching
phase_1_cubic_eos_validation_refstate_0910.py's own 09/10 note:

first_principle: Pc/Tc unchanged, pvap=7853.748 (was 7409.643). Shomate
a-e revised to a larger raw scale. omega (mmHg): -0.235212 (was
-0.209932). Cp sanity check @298.15K: 84.6 J/mol/K -- sane.

gcn: Pc/Tc unchanged, pvap=1670.067 (was 12585.513). Shomate a-e
revised. omega (mmHg): +0.450573 (was -0.426564) -- sign flip, now the
closest to real R32's literature omega (~0.277) of any gcn revision.

first_principle_multi: Pc=64.94 (was 64.965), Tc=455.096 (was 454.253),
pvap=3392.137 (was 4243.056). omega (mmHg): 0.157137 (was 0.060099).

gcn_multi: fit FAILED this round (pvap="failed", no Pc/Tc/Shomate).
Kept the 09/06 fit as a fallback rather than dropping the method.

Also added plt.savefig() calls for both diagrams (ph_diagram_0910.png,
ts_diagram_0910.png) right before plt.show(), so a run produces files
on disk in addition to the interactive windows -- purely additive, does
not change the plots themselves.

--- 09/24/2026 ---
New delivery: "spgp r32.xlsx", red date band 2026-09-19. Six yellow-
highlighted datasets, all registered in METHODS below alongside the NIST
and GCGP baselines. The 09/10 "(7).xlsx" rows (first_principle, gcn,
first_principle_multi, gcn_multi) are NOT carried forward -- see the note
above METHODS for why they are a different branch, not a prior version.

Series in this run: Linde (datasheet reference) + NIST + GCGP + the six
09/19 sets = 9 curves.

UNIT CONVENTION -- the one thing to get right. The spreadsheet reports its
ideal-gas Cp fit in the RAW polynomial convention,

    Cp = a + b*T + c*T^2 + d*T^3 + e/T^2          (T in K)

while this file (and every other file in the project) uses Shomate with
t = T/1000. The conversion, verified to machine precision against the
already-correct rows in phase_1_cubic_eos_validation_refstate_GCN.py:

    A = a          B = b * 1e3          C = c * 1e6
    D = d * 1e9    E = e * 1e-6

Confirmed independently by recomputing the sheet's own Cp column (B2:B141)
from its E5:E9 coefficients under the raw convention -- agreement to ~1e-4,
consistent with the 6-significant-figure coefficients printed in the sheet.

ACENTRIC FACTOR -- not supplied. The sheet has no omega/acentric column;
its "pvap (mmHg)" row is Psat at Tr = 0.7, which is exactly the input
omega_from_psat() needs. Every omega below is derived, not transcribed.
All six come out physical (0.052 to 0.449), a clear improvement on the
09/10 delivery (-0.235, +0.451, +0.157, and an unusable +1.254).

Figures now save as ph_diagram_0924.png / ts_diagram_0924.png so a run
does not overwrite the 09/10 figures.

Author: Shilpa Narasimhan and Claude AI
Date Created: 07/07/2026 (original compare_cp_methods.py)
This _GCN copy created: 08/26/2026
This _0903 copy created: 09/03/2026
This _0909 copy created: 09/09/2026
This _0910 copy created: 09/10/2026
This _0924 copy created: 09/24/2026
QA/testing: Shilpa Narasimhan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


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
# Linde saturation table -- 20 points extracted directly from the Linde
# datasheet.
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
    (10,   11.069,   217.74, 516.66,  1.063, 2.119),
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
    # -------------------------------------------------------------------
    # Baselines (unchanged)
    # -------------------------------------------------------------------
    {"name": "NIST", "Pc_bar": 57.85,   "Tc": 351.3,
     "A": -6.098682, "B": 179.2200, "C": -122.3682, "D": 32.30207, "E": 0.491361},
    {"name": "GCGP", "Pc_bar": 50.730,  "Tc": 355.354,
     "A": 14.161,    "B": 0.124,    "C": -6.340e-05, "D": 1.190e-8, "E": 0.0},

    # -------------------------------------------------------------------
    # 09/19/2026 delivery -- "spgp r32.xlsx", the six yellow-highlighted
    # blocks under the red 09/19 date band. Every row below is transcribed
    # from that sheet; see the docstring for the a-e -> A-E conversion.
    #
    # NOTE: this delivery is NOT a continuation of the 09/10 "(7).xlsx"
    # fits that compare_cp_methods_0910.py carries. Neither of that file's
    # pvap values (7853.748, 1670.067) appears anywhere in this sheet;
    # what does appear, un-highlighted, is the 09/06 "(5).xlsx" value
    # 7409.643. The two are separate branches of the collaborators' work,
    # so the 09/10 rows are dropped here rather than "updated".
    # -------------------------------------------------------------------

    # --- block 1: "single output after thermoreconciliation" (I13:L23) ---
    # omega (from pvap) = 0.352253 -- physical, closest of the six to real
    # R32 (~0.277). Cp @298.15K = 76.9 J/mol/K (real R32 ~55, so high but
    # the same order). CAUTION: Tc = 325.387 K puts Tr = 0.944 at the
    # 25 C-ambient condenser (T_cond = 307.15 K) -- near-critical, where the
    # PR liquid root is ill-conditioned. Expect trouble at that ambient in
    # the cycle run; the dome here is cut at 0.99*Tc so plotting is fine.
    {"name": "fp_single_0919", "Pc_bar": 52.168, "Tc": 325.387,
     "psat_pa": 1738.797 * 133.322,   # mmHg reading, converted to Pa
     "A": 32.6437, "B": 26.9797, "C": 591.354, "D": -597.447, "E": -0.0474595},
    # omega = 0.433440 -- physical. Cp @298.15K = 172.4 J/mol/K, ~3x real
    # R32. This Cp fit is suspect; flag to the Colon group.
    {"name": "gcn_single_0919", "Pc_bar": 63.686, "Tc": 342.704,
     "psat_pa": 1760.77 * 133.322,
     "A": 115.569, "B": 154.661, "C": 67.2036, "D": 174.449, "E": 0.00882475},

    # --- block 2: "agentic approach" (N13:Q21), labelled "spg" and "gcn" ---
    # omega = 0.139628. Cp @298.15K = 76.9 J/mol/K. Tc = 344.129 K gives
    # Tr = 0.893 at the hottest condenser -- comfortable.
    {"name": "spg_agentic_0919", "Pc_bar": 52.183553, "Tc": 344.12921,
     "psat_pa": 2837.94505436 * 133.322,
     "A": 32.849609, "B": 26.23072, "C": 590.20445, "D": -592.38319,
     "E": -0.048854552},
    # omega = 0.448579 -- physical. Cp @298.15K = 172.4 J/mol/K, same ~3x
    # problem as gcn_single_0919 (the two fits are near-identical).
    {"name": "gcn_agentic_0919", "Pc_bar": 63.693509, "Tc": 342.89307,
     "psat_pa": 1700.64828868 * 133.322,
     "A": 115.50434, "B": 155.06984, "C": 66.413405, "D": 174.76964,
     "E": 0.0092029555},

    # --- block 3: "multi output after recon" (I38:L48) -------------------
    # omega = 0.051900 -- physical but low. Cp @298.15K = 94.4 J/mol/K.
    # Tc = 475.885 K is far above real R32's 351.3 K, so its dome extends
    # well past the Linde table; report_errors() masks to L_T < 0.99*Tc,
    # which for this method is the whole table.
    {"name": "fp_multi_0919", "Pc_bar": 59.621, "Tc": 475.885,
     "psat_pa": 3968.231 * 133.322,
     "A": -6.96979, "B": 321.086, "C": 480.505, "D": -1492.61, "E": 0.220907},
    # omega = 0.396637 -- physical (the 09/10 gcn_multi's 1.254451 was not).
    # Cp @298.15K = 239.5 J/mol/K, ~4x real R32 -- the worst Cp fit of the
    # six. Treat results as exploratory.
    {"name": "gcn_multi_0919", "Pc_bar": 35.705, "Tc": 385.667,
     "psat_pa": 1074.461 * 133.322,
     "A": 224.051, "B": 58.5856, "C": -21.4895, "D": -5.11394,
     "E": 0.000132741},
]


# =============================================================================
# Property model for one method
# =============================================================================

class PRMethod:
    """Peng-Robinson + Shomate property model for one (Pc, Tc, Shomate) row."""

    def __init__(self, name, Pc_bar, Tc, A, B, C, D, E, omega = None):
        self.name = name
        self.Pc = Pc_bar * 1e5
        self.Tc = Tc
        self.A, self.B, self.C, self.D, self.E = A, B, C, D, E

        self.omega = omega if omega is not None else self._omega_from_linde()
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
# Acentric factor from a directly-supplied Psat (Colon group's own data)
# =============================================================================
# Same Pitzer definition _omega_from_linde() uses (omega = -1 - log10(Psat/Pc)
# at Tr = 0.7), but here Psat comes directly from the Colon group's own
# "pvap" column instead of interpolating Linde's table. Used for
# first_principle/gcn, which carry a "psat_pa" row instead of a Linde lookup.

def omega_from_psat(name, Pc_bar, psat_pa):
    Pc_pa = Pc_bar * 1e5
    omega = -1.0 - np.log10(psat_pa / Pc_pa)
    print(f"omega_from_psat: {name}: Psat={psat_pa/1e5:.5f} bar, Pc={Pc_bar:.5f} bar, "
          f"Psat/Pc={psat_pa/Pc_pa:.6f}, omega={omega:.6f}")
    return omega


# =============================================================================
# Build every method
# =============================================================================

METHODS_RESOLVED = []
for row in METHODS:
    row = dict(row)  # don't mutate METHODS itself
    if "psat_pa" in row:
        psat_pa = row.pop("psat_pa")
        row["omega"] = omega_from_psat(row["name"], row["Pc_bar"], psat_pa)
    METHODS_RESOLVED.append(row)

methods = [PRMethod(**row) for row in METHODS_RESOLVED]
domes = [m.build_dome() for m in methods]

print(f"{'method':>18}{'omega':>10}{'Tc[C]':>10}{'Tc[K]':>10}")
for m in methods:
    print(f"{m.name:>18}{m.omega:>10.6f}{m.Tc - 273.15:>10.2f}{m.Tc:>10.3f}")


# =============================================================================
# Error summary vs Linde (4 series: Linde reference + 3 methods)
# =============================================================================
# Pressure error: percent (depends on Tc, Pc, omega only).
# Enthalpy/entropy error: mean absolute error (MAE). h and s use the IIR datum
# and pass near zero at the cold end, which makes percent error unstable there;
# MAE (kJ/kg, kJ/kg/K) is the meaningful metric and is what actually tests the
# Cp fits.

def report_errors(pr_methods):
    """Print the error table and return it as a list of dicts for export."""
    rows = []
    print("\nSaturation-property error vs Linde datasheet (selected points, -130..78 C)")
    print(f"{'series':>18}{'P_MAPE%':>10}{'hf_MAE':>9}{'hg_MAE':>9}{'sf_MAE':>10}{'sg_MAE':>10}")
    print(f"{'':>18}{'':>10}{'[kJ/kg]':>9}{'[kJ/kg]':>9}{'[kJ/kgK]':>10}{'[kJ/kgK]':>10}")
    print(f"{'Linde':>18}{'ref':>10}{'ref':>9}{'ref':>9}{'ref':>10}{'ref':>10}")
    for m in pr_methods:
        mask = L_T < 0.99 * m.Tc   # only compare where this method's own dome is defined
        T_use, P_use = L_T[mask], L_P[mask]
        HF_use, HG_use = L_HF[mask], L_HG[mask]
        SF_use, SG_use = L_SF[mask], L_SG[mask]
        P = np.array([m.saturation_pressure(t)/1e5 for t in T_use])
        hf = np.array([m.hf(t) for t in T_use])
        hg = np.array([m.hg(t) for t in T_use])
        sf = np.array([m.sf(t) for t in T_use])
        sg = np.array([m.sg(t) for t in T_use])
        p_mape = np.mean(np.abs((P - P_use) / P_use)) * 100
        rows.append({"method": m.name, "P_MAPE_pct": p_mape,
                     "hf_MAE_kJ_per_kg": float(np.mean(np.abs(hf - HF_use))),
                     "hg_MAE_kJ_per_kg": float(np.mean(np.abs(hg - HG_use))),
                     "sf_MAE_kJ_per_kgK": float(np.mean(np.abs(sf - SF_use))),
                     "sg_MAE_kJ_per_kgK": float(np.mean(np.abs(sg - SG_use))),
                     "n_Linde_points_used": int(mask.sum())})
        print(f"{m.name:>18}{p_mape:>10.2f}"
              f"{np.mean(np.abs(hf - HF_use)):>9.2f}"
              f"{np.mean(np.abs(hg - HG_use)):>9.2f}"
              f"{np.mean(np.abs(sf - SF_use)):>10.4f}"
              f"{np.mean(np.abs(sg - SG_use)):>10.4f}")
    return rows


error_rows = report_errors(methods)


# =============================================================================
# Excel export (added 09/24/2026)
# =============================================================================
# Same multi-sheet pd.ExcelWriter convention this project already uses in
# optimizer_noise_free_sensitivity_0921.py and phase6_final_sensitivity_Tc_
# 0917.py. Everything the script prints, plus the full dome curves behind the
# two diagrams, so the figures can be rebuilt or replotted without re-running.
#
#   methods          -- each method's inputs and its derived omega
#   errors_vs_Linde  -- the error table printed above
#   linde_reference  -- the 20-point Linde datasheet this file benchmarks on
#   domes            -- tidy (method, T, P, hf, hg, sf, sg), 120 points each

OUT_XLSX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "compare_cp_methods_0924.xlsx")

_raw_by_name = {r["name"]: r for r in METHODS}
methods_rows = []
for m in methods:
    raw = _raw_by_name[m.name]
    methods_rows.append({
        "method": m.name,
        "Pc_bar": m.Pc/1e5, "Tc_K": m.Tc, "Tc_C": m.Tc - 273.15,
        "omega": m.omega,
        # NIST/GCGP carry omega implicitly via the Linde interpolation; the
        # 09/19 sets carry psat_pa (Psat at Tr=0.7) and have omega derived.
        "omega_source": "derived from psat (Tr=0.7)" if "psat_pa" in raw
                        else "Linde interpolation at Tr=0.7",
        "psat_Tr0.7_bar": raw["psat_pa"]/1e5 if "psat_pa" in raw else None,
        "A": m.A, "B": m.B, "C": m.C, "D": m.D, "E": m.E,
    })

dome_rows = []
for d in domes:
    for i in range(len(d["T_C"])):
        dome_rows.append({"method": d["name"], "T_C": d["T_C"][i],
                          "P_bar": d["Pd"][i],
                          "hf_kJ_per_kg": d["hf"][i], "hg_kJ_per_kg": d["hg"][i],
                          "sf_kJ_per_kgK": d["sf"][i], "sg_kJ_per_kgK": d["sg"][i]})

linde_rows = [{"T_C": r[0], "Psat_bar": r[1], "hf_kJ_per_kg": r[2],
               "hg_kJ_per_kg": r[3], "sf_kJ_per_kgK": r[4], "sg_kJ_per_kgK": r[5]}
              for r in LINDE_SAT]

with pd.ExcelWriter(OUT_XLSX) as w:
    pd.DataFrame(methods_rows).to_excel(w, sheet_name="methods", index=False)
    pd.DataFrame(error_rows).to_excel(w, sheet_name="errors_vs_Linde", index=False)
    pd.DataFrame(linde_rows).to_excel(w, sheet_name="linde_reference", index=False)
    pd.DataFrame(dome_rows).to_excel(w, sheet_name="domes", index=False)
print(f"\nwrote {OUT_XLSX}")


# =============================================================================
# Plot 1: p-h diagram (4 curves)
# =============================================================================

colors = ["tab:blue", "tab:green", "tab:red", "tab:purple", "tab:orange",
          "tab:brown", "tab:pink", "tab:gray", "tab:olive"]

plt.figure()
# Linde reference (datasheet) as black dashed with markers
plt.plot(np.concatenate([L_HF, L_HG[::-1]]),
         np.concatenate([L_P, L_P[::-1]]),
         "k--o", ms=3, lw=1, label="Linde (datasheet)")
for d, c in zip(domes, colors):
    h_dome = np.concatenate([d["hf"], d["hg"][::-1]])
    P_dome = np.concatenate([d["Pd"], d["Pd"][::-1]])
    plt.plot(h_dome, P_dome, color=c, label=d["name"])
plt.yscale("log")
plt.xlabel("h [kJ/kg]")
plt.ylabel("P [bar]")
plt.title("R-32 p-h diagram: methods vs Linde datasheet")
plt.legend()
plt.grid(True, which="both", ls=":")
plt.tight_layout()
plt.savefig("ph_diagram_0924.png", dpi=150)


# =============================================================================
# Plot 2: T-s diagram (4 curves)
# =============================================================================

plt.figure()
plt.plot(np.concatenate([L_SF, L_SG[::-1]]),
         np.concatenate([L_TC, L_TC[::-1]]),
         "k--o", ms=3, lw=1, label="Linde (datasheet)")
for d, c in zip(domes, colors):
    s_dome = np.concatenate([d["sf"], d["sg"][::-1]])
    T_dome = np.concatenate([d["T_C"], d["T_C"][::-1]])
    plt.plot(s_dome, T_dome, color=c, label=d["name"])
plt.xlabel("s [kJ/kg/K]")
plt.ylabel("T [deg C]")
plt.title("R-32 T-s diagram: methods vs Linde datasheet")
plt.legend()
plt.grid(True, ls=":")
plt.tight_layout()
plt.savefig("ts_diagram_0924.png", dpi=150)

plt.show()
