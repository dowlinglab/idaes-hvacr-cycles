"""
Sensitivity analysis (OAT + Sobol + scrambled Monte Carlo) for the
first_principle/gcn Shomate methods, against the standalone
Peng-Robinson + Shomate saturation-property model.

This file started as a copy of compare_cp_methods_GCN.py and was
repurposed for sensitivity analysis on 08/28/2026 -- it no longer does
p-h/T-s comparison plotting against Linde (that stays in
compare_cp_methods_GCN.py); the old METHODS list, report_errors(), the
Linde datasheet table, and both plot blocks were removed here as part of
that repurposing, along with a duplicate PRMethod class and a duplicate
omega_from_psat() left over from the original copy-paste.

Scope decisions (08/28/2026):
  - Shomate coefficients (A-E) are FROZEN -- no uncertainty data exists
    for them from the Colon group (confirmed: no CI columns in either
    spgp_r32 (2) or (3).xlsx; (3)'s "residual" fields are a
    self-consistency check against the Shomate formula, not real fit
    uncertainty). Only Pc_bar, Tc, and (first_principle only) psat_pa
    are varied.
  - NIST/GCGP are excluded from this file entirely: they have no CI
    data either (they resolve omega via a Linde-table lookup, not a
    supplied Psat), so there is nothing to run a sensitivity analysis
    against for them. (This is also why the Linde datasheet table and
    _omega_from_linde() were removed -- nothing left in this file needs
    either one.)
  - Property-prediction analysis only, at two representative fixed
    temperatures -- NOT COP, and NOT a full T_amb sweep. COP needs a
    third state (compressor discharge, an isentropic solve) this file
    does not build.
  - Both spgp_r32 (2).xlsx and (3).xlsx are carried side by side
    (DATASETS["v2"], DATASETS["v3"]) because the comparison between
    them is itself a finding: gcn's psat_pa confidence interval is
    degenerate (~70 orders of magnitude wide) in BOTH files, and (3)'s
    pvap point estimate gives an unphysical omega for first_principle
    under either unit reading (Pa or mmHg) -- see the DATASETS comment
    block below for the numbers. Used anyway per the decision to report
    exactly what the Colon group supplied rather than silently
    substituting a value that looks more reasonable.

Author: Shilpa Narasimhan and Claude AI
Date Created: 07/07/2026 (original compare_cp_methods.py)
This _GCN copy created: 08/28/2026 as compare_cp_methods_GCN.py;
repurposed for sensitivity analysis: 08/28/2026
QA/testing: Shilpa Narasimhan
"""

import csv
import numpy as np
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from scipy.stats import qmc, norm

# Output CSV: one long-format row per (dataset, method, point, analysis,
# parameter, qoi, metric) -- "parameter" is "ALL" for run_mc()'s joint
# rows, since Monte Carlo doesn't attribute its result to one parameter
# the way OAT/Sobol do. Written once, at the very end of the driver loop.
CSV_PATH = "sensitivity_results.csv"
CSV_FIELDS = ["dataset", "method", "point", "T_K", "analysis",
              "parameter", "qoi", "metric", "value", "n_ok", "n_total"]


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

# Two representative fixed points, not a full dome sweep -- see EDIT 5 in
# the 08/28/2026 breadcrumb entry for the full reasoning. Not T_REF
# (hf/hg are degenerate there by construction -- h_off/s_off are
# calibrated per-instance to force the IIR values exactly at T_REF, for
# ANY parameter values, so nothing shows sensitivity at that point).
#
# T_EVAP corresponds to cycle state 1 (saturated vapor leaving
# evaporator); T_COND to state 3 (saturated liquid leaving condenser) --
# see QOIS_BY_POINT below for why only one of {hf, hg} is reported at
# each.
#
# T_EVAP was moved from -10 C to -20 C (253.15 K): this lowers gcn's
# approximate evaporator-point failure rate from ~3.6% to ~1.9% (gcn's
# Tc has to drop further below its own nominal value to cross below a
# colder T_EVAP) -- a deliberate improvement, not an arbitrary change.
#
# IMPORTANT finding from testing this file: at T_COND, gcn fails ~38-42%
# of sampled draws (both datasets) -- far worse than T_EVAP's ~2-11%.
# Cause: gcn's fitted Tc (329.12 K / 55.97 C) is only ~11 C above a
# realistic R32 condenser temperature, and gcn's OWN Tc confidence
# interval reaches down to 257.05 K -- so a large fraction of gcn's
# plausible parameter space puts its critical point AT OR BELOW a normal
# condenser operating temperature, which is physically impossible for a
# working refrigerant. This is a genuine data-quality finding about gcn,
# not a bug in this script.
T_EVAP = 253.15   # -20 C
T_COND = 318.15   # 45 C

QOIS_BY_POINT = {"T_EVAP": ["Psat", "hg"], "T_COND": ["Psat", "hf"]}


# =============================================================================
# Property model -- identical physics to compare_cp_methods_GCN.py's
# PRMethod, reimplemented standalone here so this file has no dependency
# on that file's module-level report_errors()/plt.show() side effects.
# =============================================================================

class PRMethod:
    """Peng-Robinson + Shomate property model for one (Pc, Tc, Shomate) row.

    Physics: Ambrose-Walton saturation pressure correlation (Poling et al.,
    2001 sec. 7-4), Peng-Robinson cubic EOS enthalpy/entropy departure
    functions (Poling Table 6-3), and Shomate ideal-gas Cp integration
    (t = T/1000 convention). Saturated liquid/vapor h and s are anchored to
    the IIR reference state (sat. liquid at 0 C: h=200 kJ/kg, s=1.0 kJ/kg/K)
    via h_off/s_off, computed once at construction time.
    """

    def __init__(self, name, Pc_bar, Tc, A, B, C, D, E, omega):
        """Build one method's property model.

        Parameters
        ----------
        name : str
            Method label (e.g. "first_principle", "gcn") -- used only in
            error messages and printed reports.
        Pc_bar : float
            Critical pressure, bar. Converted to Pa internally (self.Pc).
        Tc : float
            Critical temperature, K.
        A, B, C, D, E : float
            Shomate ideal-gas Cp coefficients, t = T/1000 convention.
            FROZEN in the sensitivity analysis -- never sampled.
        omega : float
            Pitzer acentric factor. Must be supplied by the caller (from
            omega_from_psat() below) -- this class does not compute it
            itself; there is no Linde-table fallback in this file since
            NIST/GCGP are out of scope here.
        """
        self.name = name
        self.Pc = Pc_bar * 1e5
        self.Tc = Tc
        self.A, self.B, self.C, self.D, self.E = A, B, C, D, E
        self.omega = omega
        self.kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        self.b = Omega_B * R * self.Tc / self.Pc
        P0 = self.saturation_pressure(T_REF)
        self.h_off = H_REF - self._h_mass(T_REF, P0, "liquid")
        self.s_off = S_REF - self._s_mass(T_REF, P0, "liquid")

    # ---- saturation pressure (Ambrose-Walton) --------------------------
    def saturation_pressure(self, T):
        """Saturation pressure (Pa) at temperature T (K), Ambrose-Walton."""
        Tr = T / self.Tc
        tau = 1.0 - Tr
        f0 = (-5.97616*tau + 1.29874*tau**1.5 - 0.60394*tau**2.5 - 1.06841*tau**5) / Tr
        f1 = (-5.03365*tau + 1.11505*tau**1.5 - 5.41217*tau**2.5 - 7.46628*tau**5) / Tr
        f2 = (-0.64771*tau + 2.41539*tau**1.5 - 4.26979*tau**2.5 + 3.25259*tau**5) / Tr
        return self.Pc * np.exp(f0 + self.omega*f1 + self.omega**2*f2)

    # ---- Peng-Robinson parameters --------------------------------------
    def a_param(self, T):
        """PR attraction parameter a(T)."""
        return Omega_A * R**2 * self.Tc**2 * (1 + self.kappa*(1 - np.sqrt(T/self.Tc)))**2 / self.Pc

    def dadT(self, T):
        """da/dT, needed by the departure-function enthalpy term."""
        Tr = T / self.Tc
        return Omega_A * R**2 * self.Tc**2 / self.Pc * (
            -self.kappa * (1 + self.kappa*(1 - np.sqrt(Tr))) / (self.Tc*np.sqrt(Tr)))

    def z_roots(self, T, P):
        """Real compressibility-factor roots of the PR cubic at (T, P).

        Returns (roots, A, B) with roots sorted ascending; roots[0] is the
        liquid-like root, roots[-1] the vapor-like root. Raises ValueError
        if no physically valid root (Z > B) exists -- this is how an
        unphysical sampled parameter draw is surfaced to the caller.
        """
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
        """Indefinite Shomate enthalpy integral (kJ/kmol) at T (K)."""
        t = T/1000.0
        return 1000.0*(self.A*t + self.B*t**2/2 + self.C*t**3/3 + self.D*t**4/4 - self.E/t)

    def _int_s(self, T):
        """Indefinite Shomate entropy integral (kJ/kmol/K) at T (K)."""
        t = T/1000.0
        return (self.A*np.log(t) + self.B*t + self.C*t**2/2 + self.D*t**3/3 - self.E/(2*t**2))

    def h_ideal(self, T):
        """Ideal-gas enthalpy change from T_REF to T (kJ/kmol)."""
        return self._int_h(T) - self._int_h(T_REF)

    def s_ideal(self, T, P):
        """Ideal-gas entropy change from (T_REF, 1 bar) to (T, P) (kJ/kmol/K)."""
        return (self._int_s(T) - self._int_s(T_REF)) - R*np.log(P/1e5)

    # ---- PR departures (Poling Table 6-3) ------------------------------
    def _departure(self, T, P, phase):
        """Enthalpy/entropy departure (dh, ds) from ideal-gas at (T, P)."""
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
        """Specific enthalpy (kJ/kg), ideal-gas + departure, unanchored."""
        dh, _ = self._departure(T, P, phase)
        return (self.h_ideal(T) + dh) / MW

    def _s_mass(self, T, P, phase):
        """Specific entropy (kJ/kg/K), ideal-gas + departure, unanchored."""
        _, ds = self._departure(T, P, phase)
        return (self.s_ideal(T, P) + ds) / MW

    # ---- anchored saturated properties (kJ/kg, kJ/kg/K) ----------------
    def hf(self, T):
        """Saturated liquid specific enthalpy at T (K), IIR-anchored."""
        return self._h_mass(T, self.saturation_pressure(T), "liquid") + self.h_off

    def hg(self, T):
        """Saturated vapor specific enthalpy at T (K), IIR-anchored."""
        return self._h_mass(T, self.saturation_pressure(T), "vapor") + self.h_off

    def sf(self, T):
        """Saturated liquid specific entropy at T (K), IIR-anchored."""
        return self._s_mass(T, self.saturation_pressure(T), "liquid") + self.s_off

    def sg(self, T):
        """Saturated vapor specific entropy at T (K), IIR-anchored."""
        return self._s_mass(T, self.saturation_pressure(T), "vapor") + self.s_off


def omega_from_psat(name, Pc_bar, psat_pa):
    """Pitzer acentric factor from a directly-supplied Psat (Colon group's
    own 'pvap' data): omega = -1 - log10(Psat/Pc) at Tr = 0.7.

    Pc_bar is part of this formula's definition, not an optional extra --
    the acentric factor is defined via the REDUCED pressure Psat/Pc, so
    Pc necessarily belongs here. This also means Pc_bar plays two roles
    in the overall model (this ratio, and the direct EOS terms in
    PRMethod), and evaluate() below always recomputes omega from
    whichever Pc_bar is active for a given call so both roles stay
    consistent with each other.
    """
    Pc_pa = Pc_bar * 1e5
    omega = -1.0 - np.log10(psat_pa / Pc_pa)
    print(f"omega_from_psat: {name}: Psat={psat_pa/1e5:.5f} bar, Pc={Pc_bar:.5f} bar, "
          f"Psat/Pc={psat_pa/Pc_pa:.6f}, omega={omega:.6f}")
    return omega


# =============================================================================
# Sensitivity-analysis input data: nominal values + 95% confidence intervals
# from the Colon group's spreadsheets, spgp_r32 (2).xlsx and (3).xlsx.
#
# Only first_principle and gcn appear here -- NIST/GCGP have no supplied
# uncertainty (their omega comes from a Linde-table lookup, not a Psat
# data point), so there is nothing to vary for them.
#
# Shomate A-E are FROZEN: they sit in "nominal" only and never appear in
# "ci", so no OAT/Sobol/Monte-Carlo routine built on this dict will ever
# sample them. Decision made 08/28/2026: no uncertainty data exists for
# these coefficients in either spreadsheet (no CI columns for a-e; (3)'s
# added "residual mean/std" fields turned out to be a self-consistency
# check against the Shomate formula, not real fit uncertainty).
#
# gcn's psat_pa has NO usable CI in either file -- (2): (1.1e-30, 1.6e+40);
# (3): (1.7e-31, 2.4e+39) -- both spans ~70 orders of magnitude, almost
# certainly a bug in the Colon group's CI-computation method for gcn
# specifically (first_principle's CI is fine in both). Left out of "ci"
# for gcn in both datasets below.
#
# v3's pvap point estimate itself is also flagged as likely wrong: reading
# it as Pa (this file's convention) gives omega=1.5485 for first_principle,
# well outside the physically normal range for a refrigerant (roughly
# -0.2 to 0.4; R32's literature value is ~0.277). Kept in DATASETS anyway,
# per the decision to report exactly what the Colon group supplied rather
# than silently substituting a value that looks more reasonable.
# =============================================================================

DATASETS = {
    "v2": {  # spgp_r32 (2).xlsx
        "first_principle": dict(
            nominal=dict(Pc_bar=52.5421618866304, Tc=397.613602773852, psat_pa=338217.0,
                         A=29.1173, B=139.601, C=62.8355, D=25.3128, E=-0.000791509),
            ci=dict(Pc_bar=(42.93621073000689, 63.89123549858172),
                    Tc=(300.26931099694696, 494.9578945507577),
                    psat_pa=(117500.74741569083, 973531.9662558418)),
        ),
        "gcn": dict(
            nominal=dict(Pc_bar=62.2266425311274, Tc=329.124473597174, psat_pa=133784.622853502,
                         A=165.655, B=104.894, C=26.4086, D=-4.04051, E=-0.000122942),
            ci=dict(Pc_bar=(56.404459223309054, 68.04882583894577),
                    Tc=(257.05124641847226, 401.1977007758769)),
                    # psat_pa CI excluded: (1.1e-30, 1.6e+40), degenerate.
        ),
    },
    "v3": {  # spgp_r32 (3).xlsx -- see caveat above re: pvap value
        "first_principle": dict(
            nominal=dict(Pc_bar=52.5421618866304, Tc=397.613602773852, psat_pa=14860.742,
                         A=29.1173, B=139.601, C=62.8355, D=25.3128, E=-0.000791509),
            ci=dict(Pc_bar=(42.93621073000689, 63.89123549858172),
                    Tc=(300.26931099694696, 494.9578945507577),
                    psat_pa=(5259.81362824021, 41986.594650004125)),
        ),
        "gcn": dict(
            nominal=dict(Pc_bar=62.2266425311274, Tc=329.124473597174, psat_pa=20226.848,
                         A=165.655, B=104.894, C=26.4086, D=-4.04051, E=-0.000122942),
            ci=dict(Pc_bar=(56.404459223309054, 68.04882583894577),
                    Tc=(257.05124641847226, 401.1977007758769)),
                    # psat_pa CI excluded: (1.7e-31, 2.4e+39), still degenerate.
        ),
    },
}


# =============================================================================
# EDIT 6: single shared evaluation function. Every one of OAT/Sobol/MC
# below calls this and only this -- the omega-consistency rule (omega is
# DERIVED from Pc_bar + psat_pa, never independent) lives in exactly one
# place.
# =============================================================================

def evaluate(dataset, name, T, **overrides):
    """Build one PRMethod instance and evaluate it at temperature T.

    The single shared evaluation point for OAT, Sobol, and Monte Carlo --
    every call to any of them goes through this function, not a direct
    PRMethod construction, so the omega-consistency rule below only has
    to live in one place.

    Parameters
    ----------
    dataset : str
        "v2" or "v3" -- which spreadsheet's nominal/CI data to use.
    name : str
        "first_principle" or "gcn".
    T : float
        Temperature (K) to evaluate at -- T_EVAP or T_COND.
    **overrides
        Zero or more of {Pc_bar, Tc, psat_pa} set to a perturbed value;
        anything not passed here stays at its nominal value (Shomate
        A-E are never in overrides -- they are frozen).

    omega consistency: omega is DERIVED from Pc_bar and psat_pa
    (omega = -1 - log10(psat_pa / Pc)), not independent of them. This
    function always recomputes omega from whichever Pc_bar/psat_pa are
    active for THIS call (perturbed or nominal), so a perturbed Pc_bar is
    never paired with a stale, nominal-derived omega -- Pc_bar plays two
    roles (the EOS term and the omega ratio) and both must see the same
    value for the resulting PRMethod instance to be physically
    self-consistent.

    Returns
    -------
    dict with keys "Psat" (bar), "hf", "hg" (kJ/kg), or None if the
    sampled parameters were unphysical (e.g. sampled Tc < T, or no valid
    PR cubic root) -- callers count these as failures rather than
    letting the exception propagate.
    """
    row = dict(DATASETS[dataset][name]["nominal"])
    row.update(overrides)
    psat_pa = row.pop("psat_pa")
    omega = omega_from_psat(name, row["Pc_bar"], psat_pa)
    try:
        m = PRMethod(name=name, omega=omega, **row)
        return dict(Psat=m.saturation_pressure(T)/1e5, hf=m.hf(T), hg=m.hg(T))
    except Exception:
        return None


# =============================================================================
# EDIT 7: OAT -- vary one parameter, freeze the rest at nominal. Cheapest,
# but blind to interactions between parameters (see run_sobol() for that).
# =============================================================================

def run_oat(dataset, name, point, T, qois, N=50):
    """One-at-a-time (OAT) sensitivity: vary one parameter, freeze the rest
    at nominal, repeat per parameter in DATASETS[dataset][name]["ci"].

    sigma conversion: each parameter's "ci" entry is a 95% confidence
    interval, (lo, hi). Converting that to a standard deviation for
    rng.normal() uses

        sigma = (hi - lo) / (2 * 1.959963985)

    1.959963985 is the z-score for a 95% CI under a NORMAL distribution
    (norm.ppf(0.975) -- 95% of a standard normal's mass lies within
    +/-1.96 sigma of the mean, so the full interval width (hi - lo) spans
    2*1.96 sigma). This conversion is only correct if the Colon group's
    interval actually IS a symmetric 95% CI computed under an
    asymptotic-normality assumption -- unconfirmed, and known to be wrong
    for at least one case: gcn's psat_pa "CI" spans ~70 orders of
    magnitude and is wildly asymmetric around its point estimate, which a
    genuinely normal quantity would never produce (see DATASETS comments
    above -- gcn's psat_pa is excluded from "ci" for exactly this reason).

    Cheap, but blind to interactions between parameters -- see run_sobol()
    for the interaction-aware version.

    Returns a list of CSV row dicts (see CSV_FIELDS), one per
    (parameter, qoi), in addition to printing the same info to console.
    """
    rng = np.random.default_rng(0)
    ci = DATASETS[dataset][name]["ci"]
    nominal = DATASETS[dataset][name]["nominal"]
    print(f"  OAT (N={N}):")
    rows = []
    for pname, (lo, hi) in ci.items():
        sigma = (hi - lo) / (2 * 1.959963985)
        draws = rng.normal(nominal[pname], sigma, size=N)
        vals = {q: [] for q in qois}
        n_ok = 0
        for v in draws:
            r = evaluate(dataset, name, T, **{pname: v})
            if r is not None:
                n_ok += 1
                for q in qois:
                    vals[q].append(r[q])
        stds = {q: (np.std(vals[q]) if len(vals[q]) > 1 else float("nan")) for q in qois}
        print(f"    {pname:9s}: {n_ok}/{N} ok, " + ", ".join(f"{q}_std={stds[q]:.4f}" for q in qois))
        for q in qois:
            rows.append(dict(dataset=dataset, method=name, point=point, T_K=T, analysis="OAT",
                              parameter=pname, qoi=q, metric="std", value=stds[q],
                              n_ok=n_ok, n_total=N))
    return rows


# =============================================================================
# EDIT 8: Sobol -- all parameters move together (Saltelli sampling, scrambled
# by default in this SALib version). Decomposes output variance into S1
# (parameter acting alone) and ST (parameter alone + all its interactions).
# Skips analyze() outright if any evaluation failed -- an index computed on
# an incomplete/NaN-containing array isn't trustworthy, so it isn't reported.
# =============================================================================

def run_sobol(dataset, name, point, T, qois, N=512):
    """Sobol global sensitivity: all parameters sampled together via
    Saltelli sampling (scrambled by default in this SALib version),
    decomposing output variance into S1 (parameter acting alone) and ST
    (parameter alone + all its interactions with the others).

    Unlike OAT, this can show a parameter mattering mostly through
    interaction (ST >> S1) rather than on its own -- something OAT
    cannot detect since it never lets two parameters move together.

    Skips analyze() entirely (prints how many evaluations failed instead)
    if any evaluation in the sample returned None -- an index computed on
    an incomplete/NaN-containing array is not trustworthy and is not
    reported. This is expected to happen often for gcn, especially at
    T_COND (see the T_EVAP/T_COND comment above).

    Returns a list of CSV row dicts (see CSV_FIELDS): two per
    (parameter, qoi) (metric="S1" and metric="ST") when analyze() ran,
    or one row per qoi with metric="skipped" (value=NaN) when it didn't.
    """
    ci = DATASETS[dataset][name]["ci"]
    pnames = list(ci.keys())
    problem = {"num_vars": len(pnames), "names": pnames, "bounds": [list(ci[p]) for p in pnames]}
    X = sobol_sample.sample(problem, N, calc_second_order=False)   # scramble=True by default
    Y = {q: [] for q in qois}
    n_fail = 0
    for row in X:
        overrides = dict(zip(pnames, row))
        r = evaluate(dataset, name, T, **overrides)
        if r is None:
            n_fail += 1
            for q in qois:
                Y[q].append(np.nan)
        else:
            for q in qois:
                Y[q].append(r[q])
    n_ok = len(X) - n_fail
    print(f"  Sobol: {len(X)} evals, {n_fail} failed ({100*n_fail/len(X):.0f}%)")
    rows = []
    for q in qois:
        y = np.array(Y[q])
        if np.isnan(y).any():
            print(f"    {q}: SKIPPED -- {int(np.isnan(y).sum())} NaNs present")
            rows.append(dict(dataset=dataset, method=name, point=point, T_K=T, analysis="Sobol",
                              parameter="ALL", qoi=q, metric="skipped", value=float("nan"),
                              n_ok=n_ok, n_total=len(X)))
            continue
        Si = sobol_analyze.analyze(problem, y, calc_second_order=False, print_to_console=False)
        s1 = ", ".join(f"{n}:{s:.3f}" for n, s in zip(pnames, Si["S1"]))
        st = ", ".join(f"{n}:{s:.3f}" for n, s in zip(pnames, Si["ST"]))
        print(f"    {q}: S1=[{s1}]  ST=[{st}]")
        for pname, s1_val, st_val in zip(pnames, Si["S1"], Si["ST"]):
            rows.append(dict(dataset=dataset, method=name, point=point, T_K=T, analysis="Sobol",
                              parameter=pname, qoi=q, metric="S1", value=s1_val,
                              n_ok=n_ok, n_total=len(X)))
            rows.append(dict(dataset=dataset, method=name, point=point, T_K=T, analysis="Sobol",
                              parameter=pname, qoi=q, metric="ST", value=st_val,
                              n_ok=n_ok, n_total=len(X)))
    return rows


# =============================================================================
# EDIT 9: plain joint Monte Carlo, using a SCRAMBLED quasi-random (Sobol
# sequence) sampler rather than naive rng.normal -- fills the parameter
# space more evenly for the same N, which matters at N as low as 64.
# Unlike Sobol above, this doesn't decompose anything -- it answers "how
# uncertain is the prediction" (a percentile band), not "why".
# =============================================================================

def run_mc(dataset, name, point, T, qois, N=64, seed=2):
    """Plain joint Monte Carlo, using a scrambled quasi-random (Sobol
    sequence) sampler via scipy.stats.qmc rather than naive rng.normal --
    fills the parameter space more evenly for the same N, which matters
    at N as low as 64 (vs. plain pseudorandom draws, which converge more
    slowly and clump unevenly at small sample sizes).

    Unlike run_sobol(), this does not decompose variance by parameter --
    it answers "how uncertain is the prediction" (an actual percentile
    band on Psat/hf/hg), complementary to Sobol's "why is it uncertain".
    Rows are written with parameter="ALL" for this reason -- the result
    isn't attributable to one parameter.

    norm.ppf(unit_cube[i, j], loc=..., scale=...) converts each dimension
    of the scrambled uniform-cube sample into a normal draw with the
    given mean/sigma -- the same sigma-from-CI conversion used in
    run_oat() (see that docstring for the 1.959963985 caveat).

    Returns a list of CSV row dicts (see CSV_FIELDS), one per (qoi, metric)
    with metric in {"mean", "std", "p2.5", "p97.5"}.
    """
    ci = DATASETS[dataset][name]["ci"]
    nominal = DATASETS[dataset][name]["nominal"]
    pnames = list(ci.keys())
    sampler = qmc.Sobol(d=len(pnames), scramble=True, seed=seed)
    unit_cube = sampler.random(N)   # N should be a power of 2 for QMC's guarantees
    vals = {q: [] for q in qois}
    n_fail = 0
    for i in range(N):
        overrides = {}
        for j, pname in enumerate(pnames):
            lo, hi = ci[pname]
            sigma = (hi - lo) / (2 * 1.959963985)
            overrides[pname] = norm.ppf(unit_cube[i, j], loc=nominal[pname], scale=sigma)
        r = evaluate(dataset, name, T, **overrides)
        if r is None:
            n_fail += 1
        else:
            for q in qois:
                vals[q].append(r[q])
    n_ok = N - n_fail
    print(f"  MC (scrambled QMC, N={N}): {n_ok}/{N} ok")
    rows = []
    for q in qois:
        y = np.array(vals[q])
        if len(y) < 2:
            print(f"    {q}: too few successful draws to report")
            rows.append(dict(dataset=dataset, method=name, point=point, T_K=T, analysis="MC",
                              parameter="ALL", qoi=q, metric="insufficient_data", value=float("nan"),
                              n_ok=n_ok, n_total=N))
            continue
        p2_5, p97_5 = np.percentile(y, [2.5, 97.5])
        print(f"    {q}: mean={y.mean():.4f}, std={y.std():.4f}, 95% band=[{p2_5:.4f}, {p97_5:.4f}]")
        for metric, value in [("mean", y.mean()), ("std", y.std()), ("p2.5", p2_5), ("p97.5", p97_5)]:
            rows.append(dict(dataset=dataset, method=name, point=point, T_K=T, analysis="MC",
                              parameter="ALL", qoi=q, metric=metric, value=value,
                              n_ok=n_ok, n_total=N))
    return rows


# =============================================================================
# Run everything: both datasets x both methods x both temperature points.
# =============================================================================

if __name__ == "__main__":
    all_rows = []
    for dataset in DATASETS:
        for name in DATASETS[dataset]:
            for point, T in [("T_EVAP", T_EVAP), ("T_COND", T_COND)]:
                qois = QOIS_BY_POINT[point]
                print(f"\n=== {dataset}/{name} @ {point}={T} K ===")
                all_rows += run_oat(dataset, name, point, T, qois)
                all_rows += run_sobol(dataset, name, point, T, qois)
                all_rows += run_mc(dataset, name, point, T, qois)

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nWrote {len(all_rows)} rows to {CSV_PATH}")

    print("\n" + "="*70)
    print("FINDING: gcn's pvap confidence interval is degenerate in BOTH")
    print("(2).xlsx (1.1e-30, 1.6e+40) and (3).xlsx (1.7e-31, 2.4e+39) --")
    print("same ~70-order-of-magnitude span, only the point estimate moved.")
    print("Two independent occurrences of the same failure mode points to a")
    print("bug in the Colon group's CI-computation method for gcn, not a")
    print("one-off data-entry error. Report this with the numbers above.")
    print("="*70)
