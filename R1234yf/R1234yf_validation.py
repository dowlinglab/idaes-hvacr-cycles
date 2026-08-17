"""
# R1234yf IDAES Helmholtz Property Package -- Validation / Diagram Driver Basis
#
# Author: Shilpa Narasimhan
# Support: Claude AI
# QA / testing: Shilpa Narasimhan
# Date Created: 2026-08-17
# Organization: Dowling Lab, University of Notre Dame
#
# Description:
#   Thin driver file for validating and extending the R1234yf (HFO-1234yf)
#   Helmholtz property package defined in R1234yf.py. Imports the already
#   validated R1234yfPropertyParameterBlock rather than duplicating its EOS
#   math, so any future bug fix only ever needs to happen in one place.
#
# References:
#   [1] Lemmon, E.W., Akasaka, R. (2022)
#       "Fundamental equation of state for 2,3,3,3-tetrafluoropropene
#       (HFO-1234yf)"
#       Int. J. Thermophys. 43, 171.
#   [2] CoolProp fluid definition for R1234yf.
#   [3] IDAES-PSE general Helmholtz expression implementations.
#   [4] Danfoss technical datasheet p-H diagram (R1234yf SI Units.pdf,
#       Coolselector(R)2 v3.3.1, "Aserep" database v3.5.0) -- diagram
#       reproduction target, not an EOS data source.
#   [5] Tanaka, K., Higashi, Y. (2010)
#       "Thermodynamic properties of HFO-1234yf
#       (2,3,3,3-tetrafluoropropene)"
#       International Journal of Refrigeration, Volume 33, Issue 3,
#       Pages 474-479. https://doi.org/10.1016/j.ijrefrig.2009.10.003
#       -- independent critical-point measurement (Tc=367.85+-0.01 K,
#       rhoc=478+-3 kg/m^3, Pc=3382+-3 kPa, via visual meniscus-
#       disappearance observation), used as a non-circular cross-check
#       for the critical-point solver's result, separate from [1]'s own
#       EOS-fit values.
#   [6] Richter, M., McLinden, M.O., Lemmon, E.W. (2011)
#       "Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene
#       (R1234yf): Vapor Pressure and p-rho-T Measurements and an
#       Equation of State"
#       J. Chem. Eng. Data 56, 3254-3264. https://doi.org/10.1021/je200369m
#       -- primary source of independent experimental data used to
#       validate this codebase's EOS point-by-point (not a source of any
#       coefficients used in the EOS itself). This paper's own 15-term
#       EOS (its Table 6) is DIFFERENT from [1]'s 17-term EOS used
#       throughout this codebase, and uses slightly different critical
#       constants (rhoc=475.55 kg/m^3, Pc=3382.2 kPa, both taken directly
#       from [5] rather than fit, vs. this codebase's solved 476.69
#       kg/m^3 / 3384.35 kPa) -- so the paper's own printed deviation
#       columns (relative to ITS EOS) were not reused. Instead, the raw
#       experimental (T, P) and (T, P, rho) triples (its Tables 1-3) were
#       compared directly against this codebase's own
#       solve_saturation_at_tau (vapor pressure) and a direct density
#       root-find via find_stable_bounds + brentq (p-rho-T).
#
#       Results (135 independent points total, T=232-400 K, P up to
#       ~10 MPa):
#         - Vapor pressure (30 points, T=250-366 K): mean abs error
#           0.056%, std dev 0.120% (T>=270 K subset, 26 points: mean
#           0.027%, std dev 0.040%). Worst point 0.541% at T=250.002 K,
#           the coldest point measured -- consistent with the paper's
#           own Figure 2, which shows increased scatter for ALL compared
#           EOS/literature datasets (not just this one) at the low-T end.
#         - p-rho-T density (105 points, T=232-400 K): normal region (93
#           points, rho outside 285-761 kg/m^3) mean abs error 0.045%,
#           std dev 0.072%, max 0.277%. Near-critical region (12 points,
#           285<rho<761 kg/m^3) checked on a PRESSURE basis instead of
#           density (evaluate pressure(delta_exp, tau) and compare to the
#           reported P) since compressibility (dP/drho) approaches zero
#           near Tc, making density hypersensitive to tiny pressure
#           differences between any two EOS fits -- this is the same
#           approach the paper's own Figure 5 uses for this region, for
#           the same reason: mean abs error 0.169%, std dev 0.184%, max
#           0.473%.
#       All figures are this codebase's own computed values vs. the raw
#       measured values -- not this paper's own accuracy claims. See
#       VALIDATION_REPORT.md for the full per-point tables.
#
# Reference state:
#   h_offset/s_offset in R1234yfPropertyParameterBlock are both 0.0, and
#   this has been verified computationally, not merely assumed: solving
#   for the true saturated-liquid state at 273.15 K and evaluating raw
#   enthalpy()/entropy() there gives h=200000.23 J/kg and
#   s=1000.0011 J/(kg K), matching the IIR reference convention
#   (h=200 kJ/kg, s=1.00 kJ/(kg K) at saturated liquid, 273.15 K) to
#   within solver tolerance. The paper's own ideal-gas constants already
#   bake in this reference state, so no additive rebasing is needed --
#   enthalpy()/entropy() values, and the isentrope target values used in
#   this file, are directly comparable to NIST/REFPROP/CoolProp/
#   IIR-convention charts (including the Danfoss diagram, reference [4])
#   with no offset correction required.
#
"""

import numpy as np
from R1234yf import R1234yfPropertyParameterBlock
from scipy.optimize import least_squares, brentq
import matplotlib.pyplot as plt

params = R1234yfPropertyParameterBlock()
params.build()

def critical_point_residuals(x,params):
    """Return the 2 criticality-condition residuals at (delta, tau).

    The critical point needs 2 equations for 2 unknowns.

    r1: dP/ddelta = 0
        Where the pressure-density isotherm goes flat.
        Below Tc, this also happens at a fake local max/min
        (the van der Waals loop) -- so this alone isn't enough.

    r2: d2P/ddelta2 = 0
        Where that flat point ALSO loses its curvature.
        Only the true critical point satisfies both at once --
        that's where the loop's max and min have merged.

    Why only alphar (residual), never alpha0 (ideal gas):
        P = delta*rhoc*R*Tc/tau * (1 + delta*alphar_delta)
        The "1 +" is the entire ideal-gas contribution -- a
        constant, always exactly 1, no delta-dependence at all.
        Its derivative is always zero. All curvature (and so all
        of criticality) comes from the residual term. Physically:
        an ideal gas has no critical point (P=rho*R*T is a
        straight line) -- criticality is a non-ideal phenomenon.

    Why r2 needs a finite difference:
        r1 only needs alphar_delta and alphar_delta_delta (2
        delta-derivatives) -- both already exist as real methods.
        r2 differentiates once more, needing a 3rd delta-
        derivative (alphar_delta_delta_delta), which isn't
        implemented analytically. Approximated here via a
        central finite difference of alphar_delta_delta.

    Why CENTRAL difference, not forward/backward:
        Central differencing cancels the leading error term,
        giving O(h^2) error instead of O(h) -- much more accurate
        for the same step size h. Less numerical noise in r2
        means the solved (delta_c, tau_c) lands closer to the
        true physical critical point.

    Parameters
    ----------
    x : sequence of float
        (delta, tau) -- reduced density, inverse reduced
        temperature.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.

    Returns
    -------
    list of float
        [r1, r2] -- both approximately 0 at the true critical
        point.
    """
    delta, tau = x
    ar_d = params.alphar_delta(delta,tau)
    ar_dd = params.alphar_delta_delta(delta,tau)

    r1 = 1.0 + 2.0*delta*ar_d + delta**2* ar_dd

    h = 1e-6*delta
    ar_dd_plus = params.alphar_delta_delta(delta+h,tau)
    ar_dd_minus = params.alphar_delta_delta(delta-h,tau)
    ar_ddd = (ar_dd_plus-ar_dd_minus)/(2.0*h)

    r2 = 2.0*ar_d + 4.0*delta*ar_dd + delta**2*ar_ddd

    return[r1,r2]

def gibbs_over_RT(delta,tau,params):
    """Return g/(R*T) at (delta, tau).

    g/(R*T) = alpha_total + Z

    Used for the saturation solver's 2nd Maxwell condition:
    equal Gibbs energy between liquid and vapor branches.

    Why G, not U:
        G is the natural potential at fixed T, P -- exactly
        what's controlled here (we sweep T, enforce equal P).
        U is natural at fixed S, V instead, so it doesn't apply.
        G = U - T*S + P*V is the Legendre transform that swaps
        to T,P as the natural variables.

    Parameters
    ----------
    delta, tau : float
        Reduced density, inverse reduced temperature.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.

    Returns
    -------
    float
        g/(R*T) -- equal between two phases means equal fugacity.
    """
    alpha = params.alpha_total(delta,tau)
    Z = params.compressibility_factor(delta,tau)
    return alpha + Z

def saturation_residuals(x,tau,params):
    """Return the 2 Maxwell-condition residuals at fixed tau.

    x = (delta_liq, delta_vap) -- 2 unknowns, 2 residuals:

    r1: equal pressure between liquid and vapor branches
    r2: equal Gibbs energy (fugacity) between branches

    Both branches share the same T (same tau) by construction.

    NOTE: this 2-unknown joint solve is KNOWN BROKEN -- it has a
    spurious trivial solution at delta_liq == delta_vap (comparing
    a state to itself trivially zeroes both residuals: plugging
    x=[d, d] into this function returns exactly [0.0, 0.0] for
    ANY d and ANY tau, verified directly). Running the old
    (un-normalized) version of this solve across T=200-367 K
    landed on that trivial solution almost everywhere -- delta_liq
    and delta_vap stayed within ~1e-8 of 1.0 for nearly the whole
    sweep, instead of spreading apart as T drops (e.g. at T=250 K
    they should reach roughly delta_liq~2.61, delta_vap~0.016, per
    a hand-picked-guess solve that DID converge correctly:
    rho_liq=1245 kg/m^3, rho_vap=7.7 kg/m^3, P~132,638 Pa).
    Normalizing the residuals (as done below) does NOT fix this --
    at the trivial point both raw residuals are exactly 0, so
    dividing by anything nonzero still gives exactly 0. This
    function has since been REPLACED by the nested branch-specific
    solve (dPddelta / find_stable_bounds / solve_branch_density /
    saturation_pressure_residual / solve_saturation_at_tau, verified
    end-to-end at T=200/250/300/340/360/367 K, 100/100 points across
    the full sweep) -- sweep_saturation_dome no longer calls this
    function. Kept here only as a record of the broken approach; do
    not call this function directly.

    Parameters
    ----------
    x : sequence of float
        (delta_liq, delta_vap).
    tau : float
        Fixed inverse reduced temperature for this T.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.

    Returns
    -------
    list of float
        [r1, r2] -- both approximately 0 at the true saturation
        state for this T.
    """
    delta_liq, delta_vap = x
    r1 = (params.pressure(delta_liq,tau) - params.pressure(delta_vap,tau))/params.pressure(delta_liq,tau)
    r2 = (gibbs_over_RT(delta_liq,tau,params) - gibbs_over_RT(delta_vap,tau,params))/gibbs_over_RT(delta_liq,tau,params)
    return [r1,r2]

def dPddelta(delta,tau,params):
    """Return dP/ddelta at fixed tau (unscaled by rhoc*R*Tc/tau).

    Same expression as r1 in critical_point_residuals, but here
    tau is FIXED and we scan over delta, instead of solving for
    the one tau where two delta-roots merge into the critical
    point.

    Sign tells us whether a density is physically real at this
    tau:
        positive -> stable, real state (compress it, P goes up,
                    like any normal fluid)
        negative -> unstable, unphysical state (compress it, P
                    would drop -- never occurs in nature; an
                    artifact of this EOS's math in the two-phase
                    interior, not a real physical branch)

    Parameters
    ----------
    delta, tau : float
        Reduced density, inverse reduced temperature.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.

    Returns
    -------
    float
        dP/ddelta (dimensionless prefactor form -- same units
        convention as r1 in critical_point_residuals) at
        (delta, tau).
    """
    ar_d = params.alphar_delta(delta,tau)
    ar_dd = params.alphar_delta_delta(delta,tau)
    return 1.0 + 2.0*delta*ar_d + delta**2*ar_dd

def find_stable_bounds(tau,params,delta_max = 5.0, n_scan =3000):
    """Find the true outer vapor/liquid spinodal points at fixed tau.

    Dense-scans dPddelta over delta and finds every sign change.
    Below Tc, there can be MORE than 2 (this EOS is only fit to
    be well-behaved in the stable, single-phase region -- nothing
    constrains the unstable/metastable interior to be simple, so
    it can wiggle). Confirmed directly at n_scan=3000, delta_max=5:
        T=367 K: 2 crossings (vapor_spinodal=0.8461, liquid_spinodal=1.1696)
        T=360 K: 4 crossings (vapor_spinodal=0.6232, liquid_spinodal=1.3943)
        T=340 K: 4 crossings (vapor_spinodal=0.4512, liquid_spinodal=1.6270)
        T=300 K: 2 crossings (vapor_spinodal=0.2717, liquid_spinodal=1.9393)
        T=250 K: 4 crossings (vapor_spinodal=0.1503, liquid_spinodal=2.2338)
        T=200 K: 4 crossings (vapor_spinodal=0.0794, liquid_spinodal=2.4896)
    i.e. the number of crossings varies (2 or 4) but the outermost
    pair always comes out physically sensible regardless -- vapor
    spinodal shrinking and liquid spinodal growing as T drops away
    from Tc, exactly as expected. The TRUE stable-branch limits are
    always the FIRST crossing (vapor spinodal, scanning up from
    low delta) and the LAST crossing
    (liquid spinodal, scanning down from high delta) -- whatever
    happens in between is interior noise to avoid searching in,
    not a boundary to solve at.

    This is piece 1 of the Bug 7 fix: these two bounds are what
    let piece 2 (the branch-specific density solve) search each
    branch separately, without ever being able to land on the
    liquid branch while looking for the vapor branch or vice
    versa.

    Returns (None, None) if no sign change is found anywhere in
    the scan. That means no two-phase loop exists at this tau --
    i.e. supercritical (tau <= tau_c). This should never actually
    happen for any tau used by the saturation sweep, since that
    sweep is built to stay below tau_c by construction (starts at
    Tc-0.5 K and marches down). If (None, None) ever comes back
    for a tau that ISN'T right at the very top of the sweep,
    that's a sign n_scan needs to increase (the loop got too
    narrow to resolve), not a temperature to silently skip.

    Parameters
    ----------
    tau : float
        Fixed inverse reduced temperature.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.
    delta_max : float
        Upper end of the scan range.
    n_scan : int
        Number of scan points -- finer resolution catches thinner
        loops (e.g. close to tau_c, where the loop is very
        narrow).

    Returns
    -------
    tuple
        (vapor_spinodal, liquid_spinodal), or (None, None) if
        supercritical / no loop found in this scan.
    """
    deltas = np.linspace(1e-6,delta_max,n_scan)
    vals = np.array([dPddelta(d,tau,params) for d in deltas])
    signs = np.sign(vals)
    crossings = np.where(np.diff(signs)!=0)[0]

    if len(crossings) == 0:
        return None, None
    i = crossings[0]
    vapor_spinodal = brentq(lambda d: dPddelta(d,tau,params),deltas[i],deltas[i+1])
    j = crossings[-1]
    liquid_spinodal = brentq(lambda d: dPddelta(d,tau,params), deltas[j], deltas[j+1])
    return vapor_spinodal, liquid_spinodal

def solve_branch_density(P_target, tau, delta_lo, delta_hi, params):
    """Solve ONE branch's density for a target pressure, bounded by a wall.

    This is piece 2 of the Bug 7 fix. Given a trial pressure and a
    (delta_lo, delta_hi) bracket that stays entirely on one side of
    the unstable region (e.g. (1e-6, vapor_spinodal) for the vapor
    branch, or (liquid_spinodal, 5.0) for the liquid branch, both
    from find_stable_bounds), finds the single delta in that bracket
    where pressure(delta, tau) equals P_target.

    Because the bracket never crosses into the other branch's
    territory, this can never return "liquid" when asked for
    "vapor" or vice versa -- unlike the old saturation_residuals
    2D solve, which had no such separation and collapsed onto
    delta_liq == delta_vap almost everywhere (see that function's
    docstring). Verified directly: at T=300 K, P=900,000 Pa,
    solve_branch_density(900000, tau, 1e-6, vap_spin, params) gives
    delta=0.1131 while solve_branch_density(900000, tau, liq_spin,
    5.0, params) gives delta=2.2790 -- two clearly distinct branches,
    both independently confirmed to satisfy pressure(delta,tau) ==
    900,000 Pa to solver tolerance.

    Parameters
    ----------
    P_target : float
        Target pressure [Pa] to match on this branch.
    tau : float
        Fixed inverse reduced temperature.
    delta_lo, delta_hi : float
        Search bracket, staying entirely within one physically
        stable branch (see find_stable_bounds).
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.

    Returns
    -------
    float
        delta in [delta_lo, delta_hi] where pressure(delta,tau) ==
        P_target.
    """
    return brentq(lambda d: params.pressure(d, tau) - P_target,
                  delta_lo, delta_hi)

def saturation_pressure_residual(P_target, tau, vap_spin, liq_spin, params):
    """Return the Gibbs-energy mismatch between branches at P_target.

    This is piece 3 of the Bug 7 fix -- the residual the OUTER
    solve drives to zero by adjusting the trial pressure. For a
    given P_target:
      1. Solve the vapor branch's density at that pressure (piece
         2), bounded between near-zero delta and vap_spin.
      2. Solve the liquid branch's density at that SAME pressure
         (piece 2 again), bounded between liq_spin and a generous
         upper delta.
      3. Compare their Gibbs energies (the equal-fugacity
         condition).

    Both branches already agree on pressure by construction (piece
    2 solves each one TO that exact P_target) -- so this residual
    only needs to check the second Maxwell condition, equal Gibbs
    energy. Zero only at the true saturation pressure for this tau.

    Verified end-to-end at T=200/250/300/340/360/367 K -- e.g. at
    T=250 K the outer solve converges to Psat=132,638 Pa,
    rho_liq=1245.3 kg/m^3, rho_vap=7.7 kg/m^3, matching the
    hand-verified point referenced in saturation_residuals'
    docstring.

    Parameters
    ----------
    P_target : float
        Trial pressure [Pa] -- what the outer solve adjusts.
    tau : float
        Fixed inverse reduced temperature for this T.
    vap_spin, liq_spin : float
        Vapor/liquid spinodal bounds, from
        find_stable_bounds(tau, params).
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.

    Returns
    -------
    float
        gibbs_over_RT(liquid) - gibbs_over_RT(vapor) at P_target --
        zero at the true saturation pressure.
    """
    d_vap = solve_branch_density(P_target, tau, 1e-6, vap_spin, params)
    d_liq = solve_branch_density(P_target, tau, liq_spin, 5.0, params)
    return gibbs_over_RT(d_liq, tau, params) - gibbs_over_RT(d_vap, tau, params)

def solve_saturation_at_tau(tau, params, delta_vap_floor=1e-6, delta_liq_max=5.0):
    """Solve the true saturation state at one fixed tau -- piece 4.

    Assembles pieces 1-3 into a single call:
      1. find_stable_bounds -- get this tau's true vapor/liquid
         spinodal walls.
      2. Work out the safe trial-pressure range: the vapor branch
         can only match a pressure between P at delta_vap_floor
         (its lowest achievable, near-zero-density value) and P at
         vap_spin (its highest, the local max of the loop); the
         liquid branch can only match a pressure at or above P at
         liq_spin (its lowest achievable value there, the local
         min of the loop -- which can come out strongly NEGATIVE
         far below Tc, e.g. -7.69 MPa at T=300 K, since that's deep
         in the unphysical interior). Taking the max of the two
         branches' floors and the vapor branch's ceiling gives a
         bracket where BOTH branches are guaranteed solvable.
      3. brentq on saturation_pressure_residual within that
         bracket, to find the one pressure where both branches'
         Gibbs energies agree -- the true saturation pressure.
      4. Re-solve both branch densities at that pressure (piece 2)
         to report back.

    Verified end-to-end across the full sweep range (T from Tc-0.5K
    down to 200 K, 100 points): 100/100 converged, no failures.
    Sample points: T=200K -> rho_liq=1378.8, rho_vap=0.62,
    Psat=8,960 Pa; T=300K -> rho_liq=1085.9, rho_vap=39.6,
    Psat=713,539 Pa; T=340K -> rho_liq=900.2, rho_vap=121.6,
    Psat=1,922,913 Pa; T=367K (near Tc) -> rho_liq=572.9,
    rho_vap=386.8, Psat=3,350,212 Pa -- branches properly separate
    at low T and pinch together approaching Tc, as a real dome
    should, instead of staying stuck together everywhere like the
    old saturation_residuals solve did.

    Parameters
    ----------
    tau : float
        Fixed inverse reduced temperature for this T.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.
    delta_vap_floor : float
        Lower edge of the vapor-branch search bracket (near-zero
        density, not exactly 0 to avoid a singularity).
    delta_liq_max : float
        Upper edge of the liquid-branch search bracket.

    Returns
    -------
    tuple
        (delta_liq, delta_vap, P_sat) at this tau.
    """
    vap_spin, liq_spin = find_stable_bounds(tau, params)

    P_vap_max = params.pressure(vap_spin, tau)
    P_vap_floor = params.pressure(delta_vap_floor, tau)
    P_liq_min = params.pressure(liq_spin, tau)

    P_lo = max(P_vap_floor, P_liq_min) * (1.0 + 1e-6) + 1e-9
    P_hi = P_vap_max * (1.0 - 1e-6)

    P_sat = brentq(saturation_pressure_residual, P_lo, P_hi,
                   args=(tau, vap_spin, liq_spin, params))

    delta_vap = solve_branch_density(P_sat, tau, delta_vap_floor, vap_spin, params)
    delta_liq = solve_branch_density(P_sat, tau, liq_spin, delta_liq_max, params)

    return delta_liq, delta_vap, P_sat

def quality_line(x, dome, params):
    """Return (H, P) arrays for one constant-quality line through the dome.

    A "quality" x is the vapor mass fraction (x=0 -> saturated
    liquid, x=1 -> saturated vapor, 0<x<1 -> a mix of both). This
    computes where a fixed quality line sits inside the dome, at
    every T already solved by sweep_saturation_dome.

    Why no new solving is needed here (unlike Steps 1-2): inside
    the two-phase dome, a PURE fluid only has 1 degree of freedom
    (Gibbs phase rule: F = C - P + 2 = 1 - 2 + 2 = 1, C=1 component,
    P=2 phases). Fixing T alone therefore already fixes the
    pressure -- every quality at that T sits at the SAME P as the
    dome's own bubble/dew points there. So P for a quality line is
    just the dome's own P at that T, reused directly.

    Enthalpy is different -- it's an EXTENSIVE property (depends on
    how much mass you have), so it mixes linearly by quality:
        H(x) = (1-x)*H_liquid + x*H_vapor
    at that T. x=0 collapses to the liquid branch, x=1 collapses to
    the vapor branch, 0.5 sits exactly halfway between them in H
    (at that T's P) -- purely arithmetic, no root-finding.

    Parameters
    ----------
    x : float
        Vapor quality, between 0 and 1.
    dome : list of tuple
        (T, delta_liq, delta_vap) triples, from
        sweep_saturation_dome -- reuses that solve, doesn't repeat
        it.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.

    Returns
    -------
    tuple of np.ndarray
        (H, P) arrays, one point per T in dome, tracing this
        quality line from low-T up to near the critical point.
    """
    H_list, P_list = [], []
    for T, delta_liq, delta_vap in dome:
        tau = params.Tc / T
        H_liq_T = params.enthalpy(delta_liq, tau)
        H_vap_T = params.enthalpy(delta_vap, tau)
        H_list.append((1.0 - x) * H_liq_T + x * H_vap_T)
        P_list.append(params.pressure(delta_liq, tau))
    return np.array(H_list), np.array(P_list)

def compute_isotherm(T, params, delta_max = 5.0, delta_floor = 1e-06, n_points = 160):
    """Return (H, P) arrays tracing one full isotherm at fixed T.

    An isotherm below Tc has 3 physically distinct pieces, but they
    can be built as just 2 delta-sweeps plus one implicit straight
    line -- no extra special-casing needed:

    1. Compressed/subcooled liquid: delta swept from delta_max
       (highest density -- most compressed) DOWN to delta_liq_sat
       (the exact saturated-liquid density at this T, from
       solve_saturation_at_tau -- re-solved exactly here, not
       looked up from the coarser dome sweep grid, so the endpoint
       lines up exactly with the true bubble point).

    2. The 2-phase segment isn't swept in delta at all -- a 2-phase
       state is a MIXTURE of two densities, not a single one. It
       doesn't need its own code: since the liquid branch ends
       exactly at (H_liq_sat, Psat) and the vapor branch (below)
       starts exactly at (H_vap_sat, Psat) -- same pressure both
       ends, because it's a fixed-T saturation state -- simply
       plotting these two branches back-to-back as ONE line means
       matplotlib draws that connecting segment automatically, and
       it comes out perfectly flat because both endpoints share the
       same Psat.

    3. Superheated vapor: delta swept from delta_vap_sat DOWN to
       delta_floor (near-zero, most expanded/dilute).

    Both single-phase sweeps stay entirely within their own stable
    region (delta_max down to delta_liq_sat is always above the
    liquid spinodal; delta_vap_sat down to delta_floor is always
    below the vapor spinodal -- see find_stable_bounds), so
    pressure is guaranteed monotonic across each sweep -- no risk
    of clipping into the unstable interior.

    Verified visually: rendered together with the dome and quality
    lines (dome_isotherms.png) -- liquid branch comes out steeply
    near-vertical (liquid barely compresses), the 2-phase segment
    comes out flat as expected, and the vapor branch curves the way
    a superheated-vapor isotherm should.

    Parameters
    ----------
    T : float
        Isotherm temperature [K]. Must be below the solved Tc (this
        function does not handle supercritical isotherms).
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.
    delta_max : float
        Upper (most compressed) end of the liquid sweep.
    delta_floor : float
        Lower (most dilute) end of the vapor sweep.
    n_points : int
        Number of points per single-phase branch (2 branches, so
        2*n_points total points returned).

    Returns
    -------
    tuple of np.ndarray
        (H, P) arrays tracing the full isotherm: compressed liquid
        -> saturated liquid -> (implicit flat 2-phase segment) ->
        saturated vapor -> superheated vapor. For T at or above Tc,
        instead traces a single continuous supercritical sweep (see
        below) -- no dome, no 2-phase segment, since none exists.

    Supercritical case (T >= Tc): find_stable_bounds returns
    (None, None) here -- confirmed directly (T=368 K and above give
    None, None; T=367 K and below still give real spinodal pairs,
    consistent with the solved Tc=367.85 K). No two-phase dome
    exists above Tc, so there's nothing to re-solve via
    solve_saturation_at_tau -- the whole isotherm is just ONE
    continuous single-phase sweep from delta_max down to
    delta_floor. Added so more isotherms could be drawn crowding
    the vapor side above the critical point (matching the Danfoss
    chart's denser vapor-side isotherm spacing), without needing a
    saturation solve that wouldn't have anything to converge to.
    """
    tau = params.Tc/T
    vap_spin, liq_spin = find_stable_bounds(tau, params)
    if vap_spin is None:
        deltas = np.linspace(delta_max, delta_floor, 2 * n_points)
        H = np.array([params.enthalpy(d, tau) for d in deltas])
        P = np.array([params.pressure(d, tau) for d in deltas])
        return H, P
    delta_liq_sat,delta_vap_sat, Psat = solve_saturation_at_tau(tau, params)
    deltas_liq = np.linspace(delta_max,delta_liq_sat,n_points)
    H_liq = np.array([params.enthalpy(d,tau) for d in deltas_liq])
    P_liq = np.array([params.pressure(d,tau) for d in deltas_liq])

    deltas_vap = np.linspace(delta_vap_sat,delta_floor,n_points)
    H_vap = np.array([params.enthalpy(d,tau) for d in deltas_vap])
    P_vap = np.array([params.pressure(d,tau) for d in deltas_vap])
    return np.concatenate([H_liq,H_vap]), np.concatenate([P_liq,P_vap])

def isentrope_point_residuals(x,P_target,s_target,params):
    """Return the 2-equation residual for a single-phase isentrope point.

    ``x = [delta, tau]``. Zero when BOTH the pressure and the entropy at
    this (delta, tau) simultaneously match the targets -- unlike an
    isotherm (where tau is fixed and only pressure needs solving for),
    an isentrope needs a genuine 2-unknown solve because entropy can't
    be inverted directly the way temperature can.

    Only valid for a homogeneous single-phase state. Inside the
    two-phase dome there is no (delta, tau) that satisfies both targets,
    since a real two-phase mixture isn't a single density -- see
    ``isentrope_two_phase_segment`` for that region instead.
    """
    delta,tau = x
    P = params.pressure(delta,tau)
    s = params.entropy(delta,tau)
    return [P - P_target, s - s_target]

def solve_isentrope_point(P_target, s_target, params, x0):
    """Solve isentrope_point_residuals for one (P_target, s_target) pair.

    ``x0 = [delta_guess, tau_guess]`` is the warm-start seed (from the
    previous P step in a sweep, or a hand-picked starting guess for the
    first step). Returns ``(delta, tau, sol.success)`` -- the solved
    state and whether least_squares actually converged.
    """
    sol = least_squares(isentrope_point_residuals, x0, args =(P_target,s_target,params))
    delta,tau = sol.x
    return delta, tau, sol.success

def compute_isentrope(s_target, params, P_min, P_max, n_points,x0):
    """Sweep P (log-spaced) at fixed entropy, solving (delta, tau) at each step.

    Each step is warm-started from the previous step's converged
    ``(delta, tau)`` (``x_current``), which keeps the solve well-behaved
    across the sweep. ``x0`` seeds only the very first P value.

    Only valid for a single continuous single-phase branch (liquid OR
    vapor/supercritical, not both) -- a sweep that tries to cross the
    two-phase dome will fail or converge to a nonphysical state, since
    no homogeneous-phase solution exists in there. Splitting a full
    isentrope into its single-phase piece(s) plus the two-phase piece is
    handled by ``compute_full_isentrope``, not by this function.

    Returns
    -------
    tuple of np.ndarray
        (deltas, taus, P_values), one triple per P step.
    """
    P_values = np.geomspace(P_min,P_max, n_points)
    deltas = []
    taus = []
    x_current = x0
    for P in P_values:
        delta, tau, success = solve_isentrope_point(P, s_target, params, x_current)
        deltas.append(delta)
        taus.append(tau)
        x_current = [delta, tau]
    return np.array(deltas), np.array(taus), P_values

def isentrope_two_phase_segment(s_target, dome, params):
    """Return (H, P) arrays for the part of an isentrope inside the dome.

    Reuses dome's already-solved per-T saturation triples (T, delta_liq,
    delta_vap) from sweep_saturation_dome -- no new root-finding needed,
    same idea as quality_line. At each T, checks whether s_target falls
    between that T's saturated-liquid and saturated-vapor entropy
    (S_liq, S_vap); if so, solves directly by arithmetic which quality
    x hits s_target at that T (x = (s_target-S_liq)/(S_vap-S_liq)), then
    mixes H the same way (H = (1-x)*H_liq + x*H_vap) and reuses that T's
    saturation pressure (same P for both phases). T's where s_target
    doesn't fall in range are skipped entirely -- since the saturated
    liquid/vapor entropy gap narrows to zero at the critical point, only
    a sub-range of T's near the middle of the dome brackets any given
    s_target.

    Returns an empty pair of arrays if s_target never falls inside the
    dome across the whole T range sweep_saturation_dome covered (e.g.
    very high entropy targets that are always superheated vapor).
    """
    H_list, P_list = [],[]
    for T, delta_liq, delta_vap in dome:
        tau = params.Tc/T
        S_liq = params.entropy(delta_liq,tau)
        S_vap = params.entropy(delta_vap,tau)
        if S_liq<= s_target <= S_vap:
            x = (s_target-S_liq)/(S_vap-S_liq)
            H_liq = params.enthalpy(delta_liq,tau)
            H_vap = params.enthalpy(delta_vap,tau)
            H_list.append((1.0-x)*H_liq + x*H_vap)
            P_list.append(params.pressure(delta_liq,tau))
    return np.array(H_list), np.array(P_list)

def compute_full_isentrope(s_target, dome, params, P_max=3600000.0, P_low=50000.0, n_points=60):
    """Return one full (H, P) isentrope curve, stitching single-phase + two-phase pieces.

    This is the assembly step compute_isentrope's own docstring flags as
    needed: a real isentrope can pass through the two-phase dome, where
    compute_isentrope's homogeneous-phase solve has nothing to converge
    to. Three cases, distinguished using isentrope_two_phase_segment's
    own output:

    1. s_target never falls inside the dome anywhere in the solved T
       range (empty two-phase segment) -- the whole visible curve is
       single-phase vapor/supercritical. One continuous compute_isentrope
       sweep from P_low to P_max, seeded near a low-density vapor guess.
    2. The two-phase segment exists but its highest pressure is still
       below P_low (the chart's own floor) -- the crossing itself is
       invisible on this chart, and the entire visible range is already
       single-phase compressed liquid. One continuous compute_isentrope
       sweep from P_low to P_max, trying the dome-edge liquid state as a
       seed first and falling back to the vapor-style seed if that
       fails (this matters right at the boundary between the two entropy
       groups, e.g. s=1625 J/(kg K), where the dome-edge state is right
       next to the critical point and behaves more like the vapor seed).
    3. The two-phase segment is visible (its highest pressure exceeds
       P_low) -- return the two-phase points as-is, plus (if there's
       room below P_max) a compressed-liquid continuation above it,
       seeded from the dome's own liquid state at the crossing's
       highest-pressure edge.

    Parameters
    ----------
    s_target : float
        Target entropy for this isentrope [J/(kg K)].
    dome : list of tuple
        (T, delta_liq, delta_vap) triples from sweep_saturation_dome.
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.
    P_max, P_low : float
        Visible pressure range to cover [Pa] -- matches the plot's own
        axis limits.
    n_points : int
        Number of points for each single-phase sweep piece.

    Returns
    -------
    tuple of np.ndarray
        (H, P) arrays tracing this isentrope across whichever
        pieces actually apply.
    """
    H2, P2 = isentrope_two_phase_segment(s_target, dome, params)
    # isentrope_two_phase_segment walks `dome` from high P (near Tc) down
    # to low P (T_min) -- reverse to low-to-high P so it concatenates
    # cleanly with the liquid continuation below, which sweeps low-to-
    # high P. Without this, the plotted line draws down through the
    # dome then jumps back up to resume the liquid branch, showing as a
    # spurious second line crossing back over the first.
    H2, P2 = H2[::-1], P2[::-1]

    edge_delta_liq, edge_delta_vap, edge_tau = None, None, None
    for T, delta_liq, delta_vap in dome:
        tau = params.Tc / T
        S_liq = params.entropy(delta_liq, tau)
        S_vap = params.entropy(delta_vap, tau)
        if S_liq <= s_target <= S_vap:
            edge_delta_liq, edge_delta_vap, edge_tau = delta_liq, delta_vap, tau
            break

    if len(H2) == 0:
        deltas, taus, Pvals = compute_isentrope(s_target, params, P_low, P_max, n_points, [0.005, 1.0])
        H = np.array([params.enthalpy(d, t) for d, t in zip(deltas, taus)])
        return H, Pvals

    if P2.max() < P_low:
        for x0 in ([edge_delta_liq, edge_tau], [0.005, 1.0]):
            try:
                deltas, taus, Pvals = compute_isentrope(s_target, params, P_low, P_max, n_points, x0)
                H = np.array([params.enthalpy(d, t) for d, t in zip(deltas, taus)])
                return H, Pvals
            except Exception:
                continue
        return H2, P2

    H_all, P_all = H2, P2
    if edge_delta_liq is not None and P2.max() < P_max:
        # Stepping in too close to the saturation edge (e.g. *1.001) is
        # numerically ill-conditioned -- least_squares reports success
        # but silently lands on a garbage state (verified directly: it
        # returned entropy off by ~100 J/(kg K) from the target while
        # still reporting converged). A slightly bigger first step
        # (*1.01) avoids that, but the seed itself also matters: the
        # edge's saturated-LIQUID state is only a good starting guess
        # when s_target sits close to the liquid side of the bracket.
        # For s_target close to the VAPOR side (common near the critical
        # entropy, where liquid and vapor branches nearly coincide), the
        # true continuation above the dome passes through the
        # supercritical region at much LOWER density, not compressed
        # liquid -- confirmed directly for s=1575 J/(kg K), where only
        # the vapor-state/critical-point seeds converged to the correct
        # target entropy; the liquid seed converged to a nearby but
        # wrong state every time. Try liquid seed, then vapor seed, then
        # the critical point itself, keeping the first one that both
        # reports success AND actually reproduces s_target.
        P_start = P2.max() * 1.01
        for x0 in ([edge_delta_liq, edge_tau], [edge_delta_vap, edge_tau], [1.0, 1.0]):
            try:
                delta0, tau0, ok = solve_isentrope_point(P_start, s_target, params, x0)
                if not ok or abs(params.entropy(delta0, tau0) - s_target) > 1.0:
                    continue
                deltas, taus, Pvals = compute_isentrope(s_target, params, P_start, P_max, 40, [delta0, tau0])
                H_liq = np.array([params.enthalpy(d, t) for d, t in zip(deltas, taus)])
                H_all = np.concatenate([H2, H_liq])
                P_all = np.concatenate([P2, Pvals])
                break
            except Exception:
                continue
    return H_all, P_all

def sweep_saturation_dome(tau_c,params,T_min = 200.0, n_points = 100):
    """Sweep T downward from just below Tc, solving saturation at each step.

    WHAT USED TO HAPPEN (Bug 7): this function used to call
    saturation_residuals through a blind 2D least_squares solve on
    (delta_liq, delta_vap) together, warm-started from the previous
    T step. That joint solve had a spurious trivial solution at
    delta_liq == delta_vap (comparing a state to itself trivially
    zeroes both Maxwell residuals -- verified directly:
    saturation_residuals([d,d], tau, params) returns exactly
    [0.0, 0.0] for any d, any tau). Running it across T=200-367 K
    landed on that trivial solution almost everywhere: delta_liq
    and delta_vap stayed within ~1e-8 of 1.0 for nearly the whole
    sweep instead of spreading apart as T dropped. Normalizing the
    residuals was tried and did NOT fix it, for the same reason --
    0 divided by anything is still 0.

    WHAT FIXED IT: replaced the joint 2D solve with
    solve_saturation_at_tau (piece 4), which never lets liquid and
    vapor be parameters of the same solve. At each T it (1) scans
    for this tau's true spinodal walls (find_stable_bounds), (2)
    solves each branch's density independently, bounded so it can't
    cross into the other branch (solve_branch_density), and (3)
    adjusts a trial pressure in its own 1D solve until the two
    branches' Gibbs energies agree (saturation_pressure_residual).
    Since neither branch's solve can ever see the other branch's
    density, the delta_liq==delta_vap trap is structurally
    unreachable now. No warm-starting is needed anymore either --
    each T is solved fresh and independently, since
    find_stable_bounds finds fresh spinodal walls at every T rather
    than relying on continuity from the previous step.

    Also avoids needing any ancillary saturated-density correlation
    for an initial guess -- an earlier version of this project's
    r1234yf.json had one (aux.delta_l_sat_approx/delta_v_sat_approx)
    that was removed after it was found to predict a saturated-
    liquid density of 1780 kg/m^3 at 273.15 K, denser than this
    fluid's own triple-point liquid density of 1550 kg/m^3
    (basic.rhot_l) -- a physical impossibility, since liquid density
    can only rise as T drops toward the triple point. (CoolProp's
    real ancillary, fetched directly from its GitHub source for
    comparison, gives 1176 kg/m^3 at the same T -- consistent with
    R1234yf's known real-world density -- confirming the removed
    data was fabricated/placeholder, not a real fit.) The fix above
    doesn't need an ancillary OR an initial guess of any kind --
    each T's bounds come directly from find_stable_bounds's own
    spinodal scan.

    Verified: 100/100 points converged across the full T range
    (Tc-0.5 K down to 200 K). See solve_saturation_at_tau's
    docstring for sample (T, rho_liq, rho_vap, Psat) values.

    Parameters
    ----------
    tau_c : float
        Solved critical tau from Step 1 (use the SOLVED value,
        not an assumed 1.0 -- it may differ slightly).
    params : R1234yfPropertyParameterBlock
        Already-built EOS parameter block.
    T_min : float
        Lowest temperature to sweep down to [K].
    n_points : int
        Number of T points in the sweep.

    Returns
    -------
    list of tuple
        (T, delta_liq, delta_vap) for each converged point.
    """
    T_c = params.Tc/tau_c
    T_start = T_c -0.5

    T_values = np.linspace(T_start,T_min,n_points)
    results = []
    for T in T_values:
        tau = params.Tc/T
        delta_liq, delta_vap, P_sat = solve_saturation_at_tau(tau, params)
        results.append((T,delta_liq,delta_vap))

    return results


result = least_squares(
    critical_point_residuals, x0=[1.0, 1.0], args=(params,),
    bounds=([0.5, 0.5], [2.0, 2.0]),
)
delta_c, tau_c = result.x

T_c = params.Tc / tau_c
rho_c = delta_c * params.rhoc
P_c = params.pressure(delta_c, tau_c)

print(f"Solver converged: {result.success}, cost: {result.cost:.3e}, residual: {result.fun}")
print(f"delta_c = {delta_c:.6f}, tau_c = {tau_c:.6f}")
print(f"T_c   = {T_c:.4f} K")
print(f"rho_c = {rho_c:.4f} kg/m^3")
print(f"P_c   = {P_c:.2f} Pa  ({P_c/1e6:.4f} MPa)")
print()
print(f"JSON values:              Tc={params.Tc} K, rhoc={params.rhoc} kg/m^3, Pc={params.Pc} Pa")
print(f"Tanaka & Higashi (2010):  Tc=367.85+-0.01 K, rhoc=478+-3 kg/m^3, Pc=3.382+-0.003 MPa (independent)")


dome = sweep_saturation_dome(tau_c, params, T_min=200.0, n_points=100)

H_liq, P_liq, H_vap, P_vap = [], [], [], []
for T, d_liq, d_vap in dome:
    tau = params.Tc / T
    H_liq.append(params.enthalpy(d_liq, tau))
    P_liq.append(params.pressure(d_liq, tau))
    H_vap.append(params.enthalpy(d_vap, tau))
    P_vap.append(params.pressure(d_vap, tau))

H_liq = np.array(H_liq)
P_liq = np.array(P_liq)
H_vap = np.array(H_vap)
P_vap = np.array(P_vap)

H_c = params.enthalpy(delta_c, tau_c)

H_dome = np.concatenate([H_liq[::-1], [H_c], H_vap])
P_dome = np.concatenate([P_liq[::-1], [P_c], P_vap])

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(H_dome, P_dome, 'k-', linewidth=2)

for x in np.arange(0.1, 0.95, 0.1):
    H_q, P_q = quality_line(x, dome, params)
    ax.plot(H_q, P_q, 'k-', linewidth=0.5)

for T in list(range(230, 361, 10)) + list(range(370, 421, 10)):
    H_iso, P_iso = compute_isotherm(float(T), params)
    ax.plot(H_iso, P_iso, 'r-', linewidth=0.6)

s_targets = list(range(775, 1576, 50)) + list(range(1625, 2126, 25))
for s in s_targets:
    H_ise, P_ise = compute_full_isentrope(float(s), dome, params)
    ax.plot(H_ise, P_ise, 'b-', linewidth=0.5)

ax.set_yscale('log')
ax.set_xlim(125000, 500000)
ax.set_ylim(50000, 3600000)
ax.set_xlabel('H (J/kg)')
ax.set_ylabel('P (Pa)')
plt.savefig('dome.png', dpi=100)
