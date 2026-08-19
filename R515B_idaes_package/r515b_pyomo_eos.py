"""
r515b_pyomo_eos.py -- Stage L (part 1): native Pyomo expressions for the
R-515B Helmholtz EOS kernel, for use INSIDE the active IDAES/Pyomo NLP.

Why this file exists, separately from `r515b_helmholtz_core.py`
-------------------------------------------------------------------
MASTER TASK spec rule 44 ("no hidden SciPy solves inside the active IDAES
NLP") means the property package's Pyomo Constraints/Expressions must be
built from native Pyomo operators (pyomo.environ.exp/log, Pyomo Var
objects) so Pyomo's own automatic differentiation and IPOPT can solve them
as part of the flowsheet's equation-oriented NLP -- NOT by calling a
plain-Python/NumPy function (like `r515b_helmholtz_core.mix_state`) from
inside a Constraint body, which would hide a black-box evaluation from the
solver's Jacobian.

`r515b_helmholtz_core.py` (SciPy/NumPy-based, Stage D-K, already validated
exact vs. the oracle) is NOT deleted or superseded by this file -- per spec
rule 44's own carve-out, it remains the correct tool for INITIALIZATION
(computing good starting T/rho/composition guesses via its validated
bubble/dew/critical-point solvers before handing off to IPOPT) and for
ongoing regression/validation. This file is the separate, from-scratch
translation of the SAME validated equations into Pyomo-native form for use
inside StateBlockData's actual Constraints (Stage L, still in progress --
see this file's own STATUS note at the bottom).

Scope narrowing confirmed by direct inspection (2026-08-17): both
R-1234ze(E) (r1234ze.json) and R-227ea (r227ea.json) use
`phi_ideal_type=1` and `phi_residual_type=2` -- CONFIRMED by reading each
JSON's `eos` section directly, not assumed. `alpha0_idaes_with_derivs`/
`alphar_idaes_with_derivs` (in `linear_model_codex.py`) implement 3-4
branches each for generality across IDAES's whole fluid library; THIS
file implements ONLY the type-1 (ideal) / type-2 (residual) branches,
since those are the only ones R-515B's two pure components actually need.
If a future change ever swapped in a pure-fluid JSON using a different
phi_*_type, this file would need a corresponding new branch -- it does
NOT generalize to arbitrary phi types by design (consistent with the
spec's explicit "no generalized...framework" scope restriction, rule 1).

Traceability: every expression here is a direct line-by-line Pyomo
transcription of `linear_model_codex.py`'s `alpha0_idaes_with_derivs`
(phi==1 branch), `alphar_idaes_with_derivs` (phi==2 branch),
`bell2023_Tred_vred`, `bell2023_departure_base`/`_alphar` -- NOT an
independent re-derivation. Validated in `validate_pyomo_eos_vs_core.py`
by evaluating these Pyomo expressions at fixed (non-Var, i.e. Pyomo
Param-like float-valued) inputs via `pyomo.environ.value()` and comparing
against `r515b_helmholtz_core.py`'s own NumPy evaluation at the same
inputs.
"""

from typing import Dict

from pyomo.environ import exp, log, Expression
from pyomo.core.expr.calculus.derivatives import differentiate, Modes

from linear_model_codex import BELL_2023_R1234ZE_R227EA, BELL_2023_DEP_COEFFS, R_u  # noqa: F401

# Copied verbatim from the oracle (mixture_fully_validated.py) / from
# r515b_helmholtz_core.py's own copy, NOT imported from that file -- this
# module is deliberately kept independent of r515b_helmholtz_core.py (see
# module docstring: the two files are separate, from-scratch translations
# of the same validated equations, one SciPy/NumPy, one native-Pyomo; they
# do not import each other).
ENTROPY_REFERENCE_OFFSET_JMOLK = -1.4595


def pure_alpha0_expr(eos: Dict, tau, delta):
    """
    Native-Pyomo transcription of `alpha0_idaes_with_derivs`'s phi==1
    branch (the ONLY branch R-1234ze(E)/R-227ea need -- confirmed by
    direct JSON inspection). `tau`/`delta` may be Pyomo expressions
    (e.g. built from Vars T, rho, x1) or plain floats; this function
    itself does not care which, since Pyomo expression trees compose.

    Returns alpha0 (a Pyomo expression) ONLY -- unlike the NumPy version,
    this does NOT also return alpha0_tau, because Pyomo/IPOPT computes
    derivatives of expressions automatically via its own AD; a
    hand-coded companion derivative expression would be redundant and a
    second place for the two to silently drift apart. Any caller needing
    d(alpha0)/d(tau) should let Pyomo differentiate this expression, not
    call a separate hand-written derivative expression.
    """
    n0 = eos["n0"]
    g0 = eos["g0"]
    phi = int(eos["phi_ideal_type"])
    if phi != 1:
        raise NotImplementedError(
            f"pure_alpha0_expr only implements phi_ideal_type=1 (R-515B's two "
            f"pure components both use type 1, confirmed by direct JSON "
            f"inspection); got phi_ideal_type={phi}. This module deliberately "
            f"does not generalize to other types (spec rule 1 scope limit)."
        )
    last = int(eos["last_term_ideal"])

    alpha0 = log(delta) + float(n0["1"]) + float(n0["2"]) * tau + float(n0["3"]) * log(tau)
    for k in range(4, last + 1):
        nk = float(n0[str(k)])
        ak = float(g0[str(k)])
        alpha0 = alpha0 + nk * log(1.0 - exp(-ak * tau))
    return alpha0


def pure_alphar_expr(eos: Dict, tau, delta):
    """
    Native-Pyomo transcription of `alphar_idaes_with_derivs`'s phi==2
    branch (the ONLY branch R-1234ze(E)/R-227ea need). Returns alphar
    ONLY (a Pyomo expression) -- see `pure_alpha0_expr`'s docstring for
    why no hand-coded derivative companion is returned; Pyomo/IPOPT's own
    AD differentiates this expression directly.
    """
    phi = int(eos["phi_residual_type"])
    if phi != 2:
        raise NotImplementedError(
            f"pure_alphar_expr only implements phi_residual_type=2 (R-515B's "
            f"two pure components both use type 2, confirmed by direct JSON "
            f"inspection); got phi_residual_type={phi}."
        )
    hlist = eos["last_term_residual"]
    h1, h2, h3 = int(hlist[0]), int(hlist[1]), int(hlist[2])
    n, d, t = eos["n"], eos["d"], eos["t"]
    c, a, b, e, g = eos["c"], eos["a"], eos["b"], eos["e"], eos["g"]

    alphar = 0.0
    for i in range(1, h1 + 1):
        ni, di, ti = float(n[str(i)]), float(d[str(i)]), float(t[str(i)])
        alphar = alphar + ni * (delta ** di) * (tau ** ti)
    for i in range(h1 + 1, h2 + 1):
        ni, di, ti, ci = float(n[str(i)]), float(d[str(i)]), float(t[str(i)]), float(c[str(i)])
        alphar = alphar + ni * (delta ** di) * (tau ** ti) * exp(-(delta ** ci))
    for i in range(h2 + 1, h3 + 1):
        ni, di, ti = float(n[str(i)]), float(d[str(i)]), float(t[str(i)])
        ai, bi, ei, gi = float(a[str(i)]), float(b[str(i)]), float(e[str(i)]), float(g[str(i)])
        alphar = alphar + ni * (delta ** di) * (tau ** ti) * exp(-ai * ((delta - ei) ** 2) - bi * ((tau - gi) ** 2))
    return alphar


def bell2023_Tred_vred_expr(x1, x2, Tc1: float, Tc2: float, vc1: float, vc2: float, params):
    """Native-Pyomo transcription of `bell2023_Tred_vred` (Bell 2023 Eq.
    3-5) -- purely algebraic already, no branching, direct port. x1/x2 may
    be Pyomo expressions (e.g. a composition Var)."""
    beta_T, beta_v, gamma_T, gamma_v = params.beta_T, params.beta_v, params.gamma_T, params.gamma_v
    theta_T = (x1 + x2) / ((beta_T ** 2) * x1 + x2)
    theta_v = (x1 + x2) / ((beta_v ** 2) * x1 + x2)
    Tc12 = beta_T * gamma_T * (Tc1 * Tc2) ** 0.5
    vc12 = beta_v * gamma_v * ((vc1 ** (1.0 / 3.0) + vc2 ** (1.0 / 3.0)) ** 3) / 8.0
    Tred = (x1 ** 2) * Tc1 + (x2 ** 2) * Tc2 + 2.0 * x1 * x2 * Tc12 * theta_T
    vred = (x1 ** 2) * vc1 + (x2 ** 2) * vc2 + 2.0 * x1 * x2 * vc12 * theta_v
    return Tred, vred


def bell2023_departure_expr(x1, x2, tau, delta, dep_coeffs=BELL_2023_DEP_COEFFS):
    """Native-Pyomo transcription of the Bell (2023) Eq. 6-7 departure
    function (3-term Table 7 form, R1234ze(E)/R227ea-specific
    coefficients): alphar_dep = x1*x2*sum_k[ n_k * tau^t_k * delta^d_k *
    exp(-delta^l_k) ].

    IMPORTANT tuple order: `BELL_2023_DEP_COEFFS` entries unpack as
    (n, t, d, l) -- CONFIRMED by reading the oracle's own
    `bell2023_departure_alphar`/`bell2023_departure_base` source in
    `linear_model_codex.py` (both loop `for nk, tk, dk, lk in
    BELL_2023_DEP_COEFFS`), NOT (n, d, t, unused). An earlier draft of this
    function used the wrong order (n, d, t, unused) and also omitted the
    exp(-delta^l) damping factor entirely; that bug was caught by
    `validate_pyomo_eos_vs_core.py` (dep term off by ~13x at the
    liquid_direct state) and is recorded in helmholtz_prop_validation.md's
    failure log. This corrected version matches the oracle term-for-term:
    term_k = n_k * tau^t_k * delta^d_k * exp(-delta^l_k).

    Matches `bell2023_departure_base`/`bell2023_departure_alphar`'s
    algebraic content exactly (those functions split out tau/delta
    derivatives by hand for the SciPy path; here Pyomo/IPOPT differentiates
    the single expression automatically)."""
    base = 0.0
    for (nk, tk, dk, lk) in dep_coeffs:
        base = base + nk * (tau ** tk) * (delta ** dk) * exp(-(delta ** lk))
    return x1 * x2 * base


def _mixture_ideal_residual_parts_expr(d1: Dict, d2: Dict, x1, T, rho, Tc1: float, Tc2: float,
                                        vc1: float, vc2: float, params=BELL_2023_R1234ZE_R227EA):
    """
    Native-Pyomo assembly of the mixture ideal (a0_mix, including the
    entropy-of-mixing term) and residual (ar_mix) Helmholtz contributions
    SEPARATELY, mirroring `r515b_helmholtz_core._mix_alpha_and_derivs` /
    oracle's `_mix_alpha_and_derivs` (which sums them before returning).
    Kept as two separate expressions here -- rather than pre-summed into
    one `alpha_mix` -- because `mixture_state_expr` and
    `mixture_chemical_potentials_expr` (below) each need ar_mix alone (the
    residual-only quantity that P/Z/fugacity are built from), and
    presenting a single already-summed `alpha_mix` would force callers to
    re-derive a0_mix by subtraction. `mixture_alpha_and_derivs_expr` (kept
    below, unchanged signature, for backward compatibility with
    `validate_pyomo_eos_vs_core.py`) simply sums these two.

    x1/T/rho may be Pyomo Vars (or expressions of Vars); Tc1/Tc2/vc1/vc2
    are plain floats (fixed pure-component critical parameters, not model
    unknowns).

    Returns (tau, delta, Tred, vred, a0_mix, ar_mix).
    """
    x2 = 1.0 - x1
    Tred, vred = bell2023_Tred_vred_expr(x1, x2, Tc1, Tc2, vc1, vc2, params)
    tau = Tred / T
    delta = rho * vred
    rho_red = 1.0 / vred

    c1 = Tc1 / Tred
    c2 = Tc2 / Tred
    k1 = rho_red / (1.0 / vc1)
    k2 = rho_red / (1.0 / vc2)
    tau1 = c1 * tau
    tau2 = c2 * tau
    delta1 = k1 * delta
    delta2 = k2 * delta

    a01 = pure_alpha0_expr(d1["eos"], tau1, delta1)
    a02 = pure_alpha0_expr(d2["eos"], tau2, delta2)
    ar1 = pure_alphar_expr(d1["eos"], tau1, delta1)
    ar2 = pure_alphar_expr(d2["eos"], tau2, delta2)
    ar_dep = bell2023_departure_expr(x1, x2, tau, delta)

    a0_mix = x1 * a01 + x2 * a02 + x1 * log(x1) + x2 * log(x2)
    ar_mix = x1 * ar1 + x2 * ar2 + ar_dep
    return tau, delta, Tred, vred, a0_mix, ar_mix


def mixture_alpha_and_derivs_expr(d1: Dict, d2: Dict, x1, T, rho, Tc1: float, Tc2: float,
                                   vc1: float, vc2: float, params=BELL_2023_R1234ZE_R227EA):
    """
    Native-Pyomo assembly of the FULL mixture alpha_mix (a0_mix + ar_mix,
    including the entropy-of-mixing term), mirroring
    `r515b_helmholtz_core._mix_alpha_and_derivs` / oracle's
    `_mix_alpha_and_derivs`, but as ONE composed Pyomo expression tree
    (built from x1, T, rho -- which may be Pyomo Vars) instead of returning
    separately hand-differentiated tau/delta derivative expressions. Pyomo/
    IPOPT differentiates the returned `alpha_mix` expression automatically
    with respect to whatever Vars it was built from.

    Returns (tau, delta, alpha_mix, Tred, vred) -- Tred/vred returned too
    since callers need them again for delta1/delta2/tau1/tau2 chain-rule
    mapping (fugacity) and for P=rho*R*T*Z.
    """
    tau, delta, Tred, vred, a0_mix, ar_mix = _mixture_ideal_residual_parts_expr(
        d1, d2, x1, T, rho, Tc1, Tc2, vc1, vc2, params
    )
    alpha_mix = a0_mix + ar_mix
    return tau, delta, alpha_mix, Tred, vred


def mixture_state_expr(d1: Dict, d2: Dict, x1, T, rho, Tc1: float, Tc2: float,
                        vc1: float, vc2: float, params=BELL_2023_R1234ZE_R227EA,
                        apply_entropy_offset: bool = False):
    """
    Native-Pyomo P/h/s/g/Z expressions (Stage L part 2), built from
    `_mixture_ideal_residual_parts_expr`'s a0_mix/ar_mix using Pyomo's own
    exact SYMBOLIC differentiation (`pyomo.core.expr.calculus.derivatives.
    differentiate`, mode=reverse_symbolic) instead of hand-transcribing the
    oracle's analytic tau/delta-derivative formulas term-by-term. This is
    deliberate: a Helmholtz EOS's alpha(tau,delta,x) is, by construction, a
    genuine closed-form algebraic function of tau and delta at fixed
    composition, so its exact partial derivatives w.r.t. tau and delta ARE
    well-defined mathematical objects that symbolic differentiation
    computes exactly (not an approximation) -- using it here removes an
    entire class of hand-transcription bugs (like the tuple-order bug
    already caught and fixed in `bell2023_departure_expr`, see
    helmholtz_prop_validation.md Section 21) at the cost of one extra
    differentiation call per partial derivative needed.

    IMPORTANT precondition: T and rho MUST be actual Pyomo Var objects (not
    plain floats or arbitrary expressions) -- `differentiate(..., wrt=...)`
    requires a Var to differentiate against. x1 should also be a Var for
    consistency (used elsewhere in the module for composition
    derivatives), though this specific function does not differentiate
    w.r.t. x1.

    Derivation of the tau/rho-derivative identities used below (recorded
    here, not just in the head, per the master task's traceability
    requirement -- every non-obvious formula gets a citation or a derivation):
      tau = Tred(x1)/T  => d(tau)/dT|_{rho,x1} = -Tred/T^2 = -tau/T
        (delta = rho*vred(x1) has NO T-dependence at fixed rho,x1)
        => d(alpha_mix)/dT|_{rho,x1} = alpha_tau_mix * (-tau/T)
        => alpha_tau_mix = -(T/tau) * d(alpha_mix)/dT|_{rho,x1}
      delta = rho*vred(x1) => d(delta)/drho|_{T,x1} = vred
        (tau = Tred(x1)/T has NO rho-dependence at fixed T,x1)
        => d(ar_mix)/drho|_{T,x1} = ar_del_mix * vred
        => ar_del_mix = (1/vred) * d(ar_mix)/drho|_{T,x1}
    These are the SAME chain-rule identities implicitly used by the
    oracle's own `mix_state` (Z=1+delta*ar_del_mix, h/(RT)=1+tau*alpha_tau_
    mix+delta*ar_del_mix, g/(RT)=1+alpha+delta*ar_del_mix) and
    `_mix_entropy_direct` (s=R*(tau*alpha_tau_mix-alpha)); only the route to
    alpha_tau_mix/ar_del_mix differs (symbolic-AD here vs. hand-derived
    term sums in `mixture_alpha0_alphar_derivs`).

    Returns a dict: tau, delta, Tred, vred, alpha_mix, alpha_tau_mix,
    ar_del_mix, Z, P_pa, h_jmol, g_jmol, s_jmolk.
    """
    tau, delta, Tred, vred, a0_mix, ar_mix = _mixture_ideal_residual_parts_expr(
        d1, d2, x1, T, rho, Tc1, Tc2, vc1, vc2, params
    )
    alpha_mix = a0_mix + ar_mix

    dalpha_dT = differentiate(alpha_mix, wrt=T, mode=Modes.reverse_symbolic)
    dar_drho = differentiate(ar_mix, wrt=rho, mode=Modes.reverse_symbolic)

    alpha_tau_mix = -(T / tau) * dalpha_dT
    ar_del_mix = dar_drho / vred

    z = 1.0 + delta * ar_del_mix
    p_pa = rho * R_u * T * z
    h_jmol = R_u * T * (1.0 + tau * alpha_tau_mix + delta * ar_del_mix)
    g_jmol = R_u * T * (1.0 + alpha_mix + delta * ar_del_mix)
    s_jmolk = R_u * (tau * alpha_tau_mix - alpha_mix)
    if apply_entropy_offset:
        s_jmolk = s_jmolk + ENTROPY_REFERENCE_OFFSET_JMOLK

    return {
        "tau": tau, "delta": delta, "Tred": Tred, "vred": vred,
        "alpha_mix": alpha_mix, "alpha_tau_mix": alpha_tau_mix, "ar_del_mix": ar_del_mix,
        "Z": z, "P_pa": p_pa, "h_jmol": h_jmol, "g_jmol": g_jmol, "s_jmolk": s_jmolk,
    }


def mixture_chemical_potentials_expr(d1: Dict, d2: Dict, x1, T, rho, Tc1: float, Tc2: float,
                                      vc1: float, vc2: float, params=BELL_2023_R1234ZE_R227EA):
    """
    Native-Pyomo mu1/mu2 (chemical potential) expressions (Stage L part 2),
    mirroring `r515b_helmholtz_core.chemical_potentials_analytic` / the
    oracle's `chemical_potentials_analytic`, but derived via Pyomo's exact
    symbolic differentiation of `ar_mix` w.r.t. x1 and rho instead of
    hand-transcribing `_bell2023_reducing_derivs_binary` +
    `bell2023_departure_base`/`bell2023_departure_alphar`'s separate
    tau/delta-derivative formulas.

    Derivation (recorded here for traceability -- standard binary-mixture
    fugacity-coefficient identity for a molar-Helmholtz-explicit EOS
    parameterized by (T, rho, x1), n = n1+n2 total moles, V = n/rho volume):
      x1 = n1/n  =>  d(x1)/d(n1)|_{V,n2} = x2/n ;  d(x1)/d(n2)|_{V,n1} = -x1/n
      rho = n/V  =>  d(rho)/d(n1)|_{V,n2} = 1/V = rho/n ; same for n2
      d(n*ar_molar)/d(n1)|_{V,T,n2}
        = ar_molar + n*(d(ar_molar)/d(n1))|_{V,T,n2}
        = ar_molar + n*[ (d ar/d x1)|_{T,rho} * (x2/n) + (d ar/d rho)|_{T,x1} * (rho/n) ]
        = ar_molar + x2*(d ar/d x1)|_{T,rho} + rho*(d ar/d rho)|_{T,x1}
      and symmetrically for n2 (x2/n1-derivative gives a MINUS x1 term instead of PLUS x2,
      since d(x1)/d(n2)|_{n1} = -x1/n):
        d(n*ar_molar)/d(n2)|_{V,T,n1} = ar_molar + rho*(d ar/d rho)|_{T,x1} - x1*(d ar/d x1)|_{T,rho}
    This matches the oracle's own `d_na_dn1 = ar_mix + rho*dar_drho + x2*dar_dx1` /
    `d_na_dn2 = ar_mix + rho*dar_drho - x1*dar_dx1` term-for-term (their dar_dx1,
    dar_drho are hand-derived equivalents of the symbolic d(ar)/dx1, d(ar)/drho
    computed here). Fugacity f_i = x_i*rho*R*T*exp(d(n*ar)/dn_i); mu_i=R*T*ln(f_i).

    IMPORTANT precondition: x1 and rho MUST be actual Pyomo Var objects (T
    need not be, since this function does not differentiate w.r.t. T).

    Returns (mu1_jmol, mu2_jmol) as Pyomo expressions.
    """
    x2 = 1.0 - x1
    tau, delta, Tred, vred, a0_mix, ar_mix = _mixture_ideal_residual_parts_expr(
        d1, d2, x1, T, rho, Tc1, Tc2, vc1, vc2, params
    )

    dar_dx1 = differentiate(ar_mix, wrt=x1, mode=Modes.reverse_symbolic)
    dar_drho = differentiate(ar_mix, wrt=rho, mode=Modes.reverse_symbolic)

    d_na_dn1 = ar_mix + rho * dar_drho + x2 * dar_dx1
    d_na_dn2 = ar_mix + rho * dar_drho - x1 * dar_dx1

    f1_pa = x1 * rho * R_u * T * exp(d_na_dn1)
    f2_pa = x2 * rho * R_u * T * exp(d_na_dn2)
    mu1 = R_u * T * log(f1_pa)
    mu2 = R_u * T * log(f2_pa)
    return mu1, mu2


def saturation_residuals_expr(d1: Dict, d2: Dict, x1, T_sat, rho_l_sat, rho_v_sat, P,
                               Tc1: float, Tc2: float, vc1: float, vc2: float,
                               params=BELL_2023_R1234ZE_R227EA):
    """
    Native-Pyomo saturation-curve residuals (Stage L part 3) for the
    pseudo-pure (x1=y1=z1 FIXED) two-phase treatment adopted per the
    2026-08-18 rule-22 resolution (see helmholtz_prop_validation.md
    Section 22 and `r515b_helmholtz_core.solve_pseudopure_saturation_at_t`,
    the SciPy/initialize()-time counterpart to this function). Returns
    THREE residual expressions that should each be driven to zero by a
    Constraint in the calling StateBlockData/validation model:
        res_P_liq  = P(T_sat, rho_l_sat, x1) - P      [mechanical equilibrium, liquid side]
        res_P_vap  = P(T_sat, rho_v_sat, x1) - P      [mechanical equilibrium, vapor side]
        res_gibbs  = g(T_sat, rho_l_sat, x1) - g(T_sat, rho_v_sat, x1)   [Maxwell/equal-Gibbs]

    Given a fixed pressure P (typically the StateBlockData's own pressure
    state Var), solving these three residuals to zero for the three
    unknowns (T_sat, rho_l_sat, rho_v_sat) gives the saturation temperature
    and saturated liquid/vapor densities AT THAT PRESSURE -- the
    native-Pyomo, always-solvable-inside-the-active-NLP counterpart to
    `general_helmholtz`'s externally-compiled `t_sat_func`/saturated-
    density correlations (see this project's own research note in
    PROJECT_CONTEXT.md's 2026-08-18 entries: `general_helmholtz` uses
    precompiled smooth correlations for T_sat(P)/rho_sat(P); we have no
    such correlation library, so we solve the EXACT same defining
    equations implicitly, as three ordinary Pyomo Constraints, which is
    fully rule-44-compliant -- no SciPy call happens here).

    T_sat, rho_l_sat, rho_v_sat, P, x1 may be Pyomo Vars or expressions of
    Vars (P is typically the StateBlockData's pressure Var; T_sat/rho_l_sat/
    rho_v_sat are typically new auxiliary Vars introduced specifically to
    carry this system). Tc1/Tc2/vc1/vc2 are plain floats.

    Derivation of why the third condition is g_l=g_v (not separate mu1/mu2
    equalities): see the module-level comment above
    `r515b_helmholtz_core.solve_pseudopure_saturation_at_t` for the full
    argument (Euler relation g=x1*mu1+x2*mu2; with x FIXED identical on
    both phases, the compositional degree of freedom that made two
    independent mu-equalities necessary no longer exists, and requiring
    both would over-determine this 3-unknown system).
    """
    state_l = mixture_state_expr(d1, d2, x1, T_sat, rho_l_sat, Tc1, Tc2, vc1, vc2, params)
    state_v = mixture_state_expr(d1, d2, x1, T_sat, rho_v_sat, Tc1, Tc2, vc1, vc2, params)
    res_p_liq = state_l["P_pa"] - P
    res_p_vap = state_v["P_pa"] - P
    res_gibbs = state_l["g_jmol"] - state_v["g_jmol"]
    return res_p_liq, res_p_vap, res_gibbs


def mixture_ph_flash_residuals_expr(d1: Dict, d2: Dict, x1, P, H,
                                     T_sat, rho_l_sat, rho_v_sat, vapor_frac,
                                     T_liq, rho_liq, T_vap, rho_vap,
                                     Tc1: float, Tc2: float, vc1: float, vc2: float,
                                     params=BELL_2023_R1234ZE_R227EA, eps_jmol: float = 1.0):
    """
    Native-Pyomo smooth single-phase/two-phase PH-flash Constraint system
    (Stage L part 3) for the pseudo-pure (x1=y1=z1 FIXED) treatment. Given
    fixed pressure P and molar enthalpy H (the StateBlockData's own state
    Vars), determines whether the state is subcooled liquid, two-phase, or
    superheated vapor -- WITHOUT explicit if/then branching, so the whole
    system stays smooth and differentiable for IPOPT.

    Design (recorded here for traceability; mirrors the SPIRIT of
    `general_helmholtz`'s pure-fluid smooth complementarity, adapted to our
    situation -- we have no compiled saturation-correlation library, so
    everything below is built from our own already-validated
    `mixture_state_expr`/`saturation_residuals_expr`, using
    `idaes.core.util.math.smooth_max` for the smoothing):

    8 unknowns: T_sat, rho_l_sat, rho_v_sat (saturation curve), vapor_frac
    (quality), T_liq/rho_liq (liquid single-phase branch), T_vap/rho_vap
    (vapor single-phase branch). 8 equations (all returned by this
    function, to be wrapped in Constraints by the caller):
      1-3. the saturation-curve system (`saturation_residuals_expr`)
      4. smooth complementarity pinning vapor_frac:
         h_over_sat  = smooth_max(0, H - h_v_sat, eps)   (>0 iff superheated)
         h_under_sat = smooth_max(0, h_l_sat - H, eps)   (>0 iff subcooled)
         0 == vapor_frac*h_over_sat - (1-vapor_frac)*h_under_sat
         (when truly subcooled, h_over_sat~0 so this forces vapor_frac~0;
         when truly superheated, h_under_sat~0 so this forces vapor_frac~1;
         when truly two-phase, BOTH slacks are ~0 so this equation is
         trivially satisfied and vapor_frac is instead pinned by equation 6
         below via the lever-rule enthalpy balance)
      5-6. LIQUID branch, fed a CLIPPED target enthalpy so it never has to
         represent an enthalpy outside the real liquid domain:
         H_liq_target = H - smooth_max(0, H - h_l_sat, eps)   (clips at h_l_sat)
         P(T_liq,rho_liq) == P ; h(T_liq,rho_liq) == H_liq_target
         (when H<=h_l_sat, H_liq_target=H exactly -- real subcooled solve;
         otherwise H_liq_target=h_l_sat exactly, so T_liq=T_sat, rho_liq=
         rho_l_sat identically -- sits AT the boundary, never extrapolates)
      7-8. VAPOR branch, symmetric clipping at h_v_sat (floors instead of
         caps): H_vap_target = H + smooth_max(0, h_v_sat - H, eps)
         P(T_vap,rho_vap) == P ; h(T_vap,rho_vap) == H_vap_target

    Given a solution, ANY extensive-like mixture property Y is then
    recovered via the SAME always-valid blend Y_actual = (1-vapor_frac)*
    Y(T_liq,rho_liq) + vapor_frac*Y(T_vap,rho_vap) -- correct in all three
    regions, since outside the two-phase dome one of the two branches
    carries zero weight and the OTHER branch holds the real single-phase
    answer (its own target was never clipped there). In particular
    T_actual = (1-vapor_frac)*T_liq + vapor_frac*T_vap collapses to T_sat
    identically throughout the two-phase region (both branches sit at the
    boundary there) and to the real single-phase T outside it.

    Returns a dict of 8 residual expressions (keys: "sat_p_liq",
    "sat_p_vap", "sat_gibbs", "complementarity", "liq_p", "liq_h", "vap_p",
    "vap_h") plus the derived "T_actual" expression, for the caller to wrap
    in Constraints (and to define its own density/entropy/etc. blends via
    the same (1-vapor_frac)/vapor_frac pattern).
    """
    from idaes.core.util.math import smooth_max

    res_sat_p_liq, res_sat_p_vap, res_sat_gibbs = saturation_residuals_expr(
        d1, d2, x1, T_sat, rho_l_sat, rho_v_sat, P, Tc1, Tc2, vc1, vc2, params
    )
    state_l_sat = mixture_state_expr(d1, d2, x1, T_sat, rho_l_sat, Tc1, Tc2, vc1, vc2, params)
    state_v_sat = mixture_state_expr(d1, d2, x1, T_sat, rho_v_sat, Tc1, Tc2, vc1, vc2, params)
    h_l_sat = state_l_sat["h_jmol"]
    h_v_sat = state_v_sat["h_jmol"]

    h_over_sat = smooth_max(0.0, H - h_v_sat, eps_jmol)
    h_under_sat = smooth_max(0.0, h_l_sat - H, eps_jmol)
    # IMPORTANT pairing (self-caught bug, see helmholtz_prop_validation.md
    # Section 22's matching entry): h_under_sat is the SUBCOOLED indicator
    # (active when H<h_l_sat) and must multiply vapor_frac (forcing it to
    # 0 when active); h_over_sat is the SUPERHEATED indicator (active when
    # H>h_v_sat) and must multiply (1-vapor_frac) (forcing it to 0, i.e.
    # vapor_frac->1, when active) -- mirroring general_helmholtz's own
    # `eq_complementarity` pairing (their pressure_over_sat, the SUBCOOLED
    # indicator, multiplies vf; pressure_under_sat, the SUPERHEATED
    # indicator, multiplies (1-vf)). An earlier draft paired these
    # backwards (h_over_sat*vf - h_under_sat*(1-vf)), which forced
    # vapor_frac to exactly 1-(the correct value) at every test point --
    # caught by `validate_pyomo_flash_vs_core.py`'s exact 1-vf inversion
    # pattern across all 4 region tests.
    res_complementarity = vapor_frac * h_under_sat - (1.0 - vapor_frac) * h_over_sat

    h_liq_target = H - smooth_max(0.0, H - h_l_sat, eps_jmol)
    h_vap_target = H + smooth_max(0.0, h_v_sat - H, eps_jmol)

    state_liq = mixture_state_expr(d1, d2, x1, T_liq, rho_liq, Tc1, Tc2, vc1, vc2, params)
    state_vap = mixture_state_expr(d1, d2, x1, T_vap, rho_vap, Tc1, Tc2, vc1, vc2, params)
    res_liq_p = state_liq["P_pa"] - P
    res_liq_h = state_liq["h_jmol"] - h_liq_target
    res_vap_p = state_vap["P_pa"] - P
    res_vap_h = state_vap["h_jmol"] - h_vap_target

    t_actual = (1.0 - vapor_frac) * T_liq + vapor_frac * T_vap

    # Entropy, additive extension (2026-08-18, Stage N prep): Honeywell-
    # rebased molar entropy on each branch, via the SAME `mixture_state_
    # expr`("s_jmolk") already computed above for state_liq/state_vap --
    # no new EOS machinery, just exposing an already-validated field
    # (Section 24b) plus the same offset convention `r515b_helmholtz_core.
    # mix_entropy_direct` applies unconditionally. Backward compatible:
    # only NEW dict keys added, nothing above this line changed, so no
    # re-validation risk to the existing 8-residual system
    # (`validate_pyomo_flash_vs_core.py`/`validate_state_block_
    # construction.py` both still pass unchanged after this addition).
    s_liq_jmolk = state_liq["s_jmolk"] + ENTROPY_REFERENCE_OFFSET_JMOLK
    s_vap_jmolk = state_vap["s_jmolk"] + ENTROPY_REFERENCE_OFFSET_JMOLK
    s_actual_jmolk = (1.0 - vapor_frac) * s_liq_jmolk + vapor_frac * s_vap_jmolk

    return {
        "sat_p_liq": res_sat_p_liq, "sat_p_vap": res_sat_p_vap, "sat_gibbs": res_sat_gibbs,
        "complementarity": res_complementarity,
        "liq_p": res_liq_p, "liq_h": res_liq_h,
        "vap_p": res_vap_p, "vap_h": res_vap_h,
        "T_actual": t_actual, "h_l_sat": h_l_sat, "h_v_sat": h_v_sat,
        "s_liq_jmolk": s_liq_jmolk, "s_vap_jmolk": s_vap_jmolk, "S_actual_jmolk": s_actual_jmolk,
    }


# =============================================================================
# STATUS (Stage L, part 2 -- 2026-08-18): this file now implements, as
# native Pyomo expressions:
#   - the EOS KERNEL (alpha0, alphar, Bell reducing functions, departure
#     function, assembled mixture alpha_mix / a0_mix / ar_mix) -- Stage L
#     part 1, validated in `validate_pyomo_eos_vs_core.py`
#     (rel_err 0 to 2.2e-16 vs. r515b_helmholtz_core.py at 3 states, after
#     fixing the departure-function tuple-order/missing-exp-term bug
#     recorded in helmholtz_prop_validation.md Section 21)
#   - P/h/s/g/Z (`mixture_state_expr`) -- Stage L part 2, built via Pyomo's
#     own exact symbolic differentiation of alpha_mix/ar_mix w.r.t. the T
#     and rho Vars (NOT hand-transcribed tau/delta-derivative formulas),
#     using the chain-rule identities derived in that function's docstring
#   - mu1/mu2 chemical potentials (`mixture_chemical_potentials_expr`) --
#     Stage L part 2, built via Pyomo's symbolic differentiation of ar_mix
#     w.r.t. the x1 and rho Vars, using the standard binary-mixture
#     fugacity-coefficient identity derived in that function's docstring
# Both new functions REQUIRE their differentiated-against arguments (T,rho
# for mixture_state_expr; x1,rho for mixture_chemical_potentials_expr) to
# be actual Pyomo Var objects, not plain floats or composed expressions --
# `pyomo.core.expr.calculus.derivatives.differentiate` needs a Var to
# differentiate against. Validated in `validate_pyomo_state_vs_core.py`
# against `r515b_helmholtz_core.py`'s `mix_state`/`mix_entropy_direct`/
# `chemical_potentials_analytic` at the same 3 representative states used
# throughout this task.
#
# NOT YET DONE, deliberately deferred to a follow-up continuation:
#   - the actual PhysicalParameterBlock/StateBlockData classes wiring
#     these expressions into Constraints/Expressions on a real Pyomo
#     Block, with PH state variables (matching vapor_compression.py's
#     PH mode), phase-equilibrium constraints (P_L=P_V, mu_i_L=mu_i_V),
#     and an initialize() method that calls r515b_helmholtz_core.py's
#     validated SciPy solvers to seed the NLP (not yet written)
#   - metadata, units, scaling (not yet written)
#   - a numerical-safety guard on log(f_i) analogous to the oracle's
#     max(1e-300, ...) (not yet added -- an initialize()/scaling-time
#     concern per spec rule 44's own carve-out, not a kernel-correctness
#     issue; noted here so it isn't forgotten)
# See helmholtz_prop_validation.md Section 22 for the full architecture
# decision record and PROJECT_CONTEXT.md for the next-action note.
# =============================================================================
