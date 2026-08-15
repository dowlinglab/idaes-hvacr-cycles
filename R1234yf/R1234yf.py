################################################################################
# R1234yf IDAES Helmholtz Property Package
#
# Author: Shilpa Narasimhan
# Support: Claude AI, ChatGPT, Codex
# QA / testing: Shilpa Narasimhan
# Date Created: 2026-08-14
# Organization: Dowling Lab, University of Notre Dame
#
# Description:
#   Custom R1234yf (HFO-1234yf) Helmholtz-property implementation based on the
#   Lemmon & Akasaka (2022) fundamental equation of state. The ideal and
#   residual Helmholtz structures are written to match the corresponding
#   IDAES general-Helmholtz Type-01 ideal and Type-02 residual forms.
#
# Context breadcrumb:
#   This revision corrects the ideal Helmholtz expression, restores the
#   exponentially damped residual terms, replaces phi notation with alpha,
#   and corrects the standard Helmholtz formulas for pressure, enthalpy,
#   entropy, heat capacities, and speed of sound.
#
# References:
#   [1] Lemmon, E.W., Akasaka, R. (2022)
#       "Fundamental equation of state for 2,3,3,3-tetrafluoropropene
#       (HFO-1234yf)"
#       Int. J. Thermophys. 43, 171.
#   [2] CoolProp fluid definition for R1234yf.
#   [3] IDAES-PSE general Helmholtz expression implementations.
#
################################################################################

import json
import os

import pyomo.environ as pyo

import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)


################################################################################
# R1234yf Property Parameter Block Class
################################################################################


class R1234yfPropertyParameterBlock(pyo.Block):
    r"""Store R1234yf Helmholtz EOS parameters and property expressions.

    The dimensionless Helmholtz free energy is written as

    .. math::

        \alpha(\delta, \tau) = \alpha^0(\delta, \tau)
        + \alpha^r(\delta, \tau),

    with reduced density and inverse reduced temperature

    .. math::

        \delta = \rho / \rho_c, \qquad \tau = T_c / T.

    The ideal contribution follows the IDAES ideal Type-01 algebra, while the
    residual contribution follows the IDAES residual Type-02 algebra. The
    method names in this file use ``alpha`` because it is the conventional symbol
    for the dimensionless Helmholtz free energy in the CoolProp formulation.

    The JSON file is expected to contain the coefficient arrays used by the
    original development file, including ``n0``, ``g0``, ``n``, ``d``, ``t``,
    ``c``, ``a``, ``b``, ``e``, and ``g`` under ``params['eos']``. For R1234yf,
    the residual structure contains five undamped power terms, five
    exponentially damped power terms, and seven Gaussian terms.

    Notes
    -----
    The reference targets used here are the common IIR-style refrigerant
    reference values for saturated liquid at 0 degC:

    * ``h_ref = 200 kJ/kg``
    * ``s_ref = 1 kJ/(kg K)``

    The offsets are recomputed from the corrected EOS during ``build`` rather
    than being hard-coded from an earlier, inconsistent implementation.
    """

    def build(self):
        """Load EOS data, construct parameters, and compute reference offsets.

        The build sequence is deliberately ordered so that all Helmholtz
        coefficients exist before the reference-state enthalpy and entropy are
        evaluated. This avoids the circular/hard-coded reference-state logic in
        the earlier implementation.
        """
        # ========== Load JSON parameter file ==========
        json_file = os.path.join(os.path.dirname(__file__), "r1234yf.json")
        with open(json_file, "r", encoding="utf-8") as handle:
            params = json.load(handle)

        # ========== Basic constants ==========
        # NOTE (2026-08-14): stored as plain Python floats, not mutable
        # pyo.Param objects. This class subclasses bare pyo.Block and calls
        # build() manually (not via Pyomo's own construct() lifecycle), so
        # pyo.Param/pyo.Set sub-components never actually get constructed --
        # any attempt to read them (even later in this same build() method)
        # raises "Cannot iterate/evaluate before it has been constructed."
        # Plain Python objects have no such lifecycle requirement and were
        # already the intended design per BREADCRUMB.md's Session 2 design
        # decision; this restores that decision for the constants that had
        # drifted back to pyo.Param during the residual-restructuring pass.
        self.Tc = float(params["basic"]["Tc"])  # Critical temperature [K]
        self.Pc = float(params["basic"]["Pc"])  # Critical pressure [Pa]
        self.rhoc = float(params["basic"]["rhoc"])  # Critical mass density [kg/m^3]
        self.R_gas = float(params["basic"]["R"])  # Specific gas constant [J/(kg K)]
        self.MW = float(params["basic"]["MW"])  # Molecular weight [kg/kmol]

        # ========== Ideal Helmholtz coefficients: IDAES Type 01 ==========
        n0_dict = {int(k): v for k, v in params["eos"]["n0"].items()}
        g0_dict = {int(k): v for k, v in params["eos"]["g0"].items()}

        self.ideal_index = sorted(n0_dict.keys())
        self.planck_einstein_index = sorted(i for i in n0_dict if i >= 4)
        self.n0 = dict(n0_dict)  # Ideal Helmholtz amplitude coefficients
        self.g0 = {
            i: g0_dict[i] for i in sorted(g0_dict) if i >= 4
        }  # Planck-Einstein characteristic-temperature coefficients

        # ========== Residual Helmholtz coefficients: IDAES Type 02 ==========
        eos = params["eos"]
        n_dict = {int(k): v for k, v in eos["n"].items()}
        d_dict = {int(k): v for k, v in eos["d"].items()}
        t_dict = {int(k): v for k, v in eos["t"].items()}
        c_dict = {int(k): v for k, v in eos.get("c", {}).items()}
        a_dict = {int(k): v for k, v in eos.get("a", {}).items()}
        b_dict = {int(k): v for k, v in eos.get("b", {}).items()}
        e_dict = {int(k): v for k, v in eos.get("e", {}).items()}
        g_dict = {int(k): v for k, v in eos.get("g", {}).items()}

        # CoolProp R1234yf contains 5 polynomial + 5 damped + 7 Gaussian
        # residual terms. If an IDAES-style boundary list is provided in the
        # JSON, use it; otherwise use the verified R1234yf boundaries.
        last_terms = list(eos.get("last_term_residual", [5, 10, 17]))
        if len(last_terms) != 3:
            raise ValueError(
                "R1234yf residual Type-02 requires three term boundaries: "
                "[last_polynomial, last_damped, last_gaussian]."
            )

        poly_last, damped_last, gaussian_last = [int(v) for v in last_terms]
        if not (1 <= poly_last < damped_last < gaussian_last):
            raise ValueError(
                "Invalid residual term boundaries in r1234yf.json: "
                f"{last_terms}."
            )

        self.residual_polynomial_index = list(range(1, poly_last + 1))
        self.residual_damped_index = list(range(poly_last + 1, damped_last + 1))
        self.residual_gaussian_index = list(range(damped_last + 1, gaussian_last + 1))

        self.residual_all_index = list(range(1, gaussian_last + 1))
        self.n_res = {i: n_dict[i] for i in self.residual_all_index}  # amplitudes
        self.d_res = {i: d_dict[i] for i in self.residual_all_index}  # density exponents
        self.t_res = {i: t_dict[i] for i in self.residual_all_index}  # inverse-temp exponents

        missing_c = [i for i in self.residual_damped_index if i not in c_dict]
        if missing_c:
            raise KeyError(
                "Missing residual damping exponents c[i] for terms "
                f"{missing_c}."
            )

        self.c_res = {
            i: c_dict[i] for i in self.residual_damped_index
        }  # density-exponential damping exponents

        gaussian_dicts = {
            "a": a_dict,
            "b": b_dict,
            "e": e_dict,
            "g": g_dict,
        }
        for name, coefficient_dict in gaussian_dicts.items():
            missing = [
                i for i in self.residual_gaussian_index if i not in coefficient_dict
            ]
            if missing:
                raise KeyError(
                    f"Missing Gaussian residual coefficients {name}[i] for "
                    f"terms {missing}."
                )

        self.a_res = {
            i: a_dict[i] for i in self.residual_gaussian_index
        }  # Gaussian width coefficient in reduced density
        self.b_res = {
            i: b_dict[i] for i in self.residual_gaussian_index
        }  # Gaussian width coefficient in inverse reduced temperature
        self.e_res = {
            i: e_dict[i] for i in self.residual_gaussian_index
        }  # Gaussian center in reduced density
        self.g_res = {
            i: g_dict[i] for i in self.residual_gaussian_index
        }  # Gaussian center in inverse reduced temperature

        # ========== Saturated-liquid density approximation ==========
        sat_liq = params["aux"]["delta_l_sat_approx"]
        delta_l_sat_n = {int(k): v for k, v in sat_liq["n"].items()}
        delta_l_sat_t = {int(k): v for k, v in sat_liq["t"].items()}
        delta_l_sat_c = sat_liq["c"]

        sat_index = sorted(set(delta_l_sat_n).intersection(delta_l_sat_t))
        if not sat_index:
            raise ValueError(
                "No saturated-liquid density approximation coefficients were found."
            )

        self.delta_l_sat_index = list(sat_index)
        self.delta_l_sat_c = float(delta_l_sat_c)  # sat-liquid reduced-density approx constant
        self.delta_l_sat_n = {i: delta_l_sat_n[i] for i in sat_index}  # amplitudes
        self.delta_l_sat_t = {i: delta_l_sat_t[i] for i in sat_index}  # exponents

        # ========== Reference-state offsets ==========
        T_ref = 273.15
        theta_ref = 1.0 - T_ref / self.Tc
        delta_sat_liq_ref = self.delta_l_sat_c + sum(
            self.delta_l_sat_n[i] * theta_ref ** self.delta_l_sat_t[i]
            for i in self.delta_l_sat_index
        )
        tau_ref = self.Tc / T_ref

        h_ref_target = 200_000.0
        s_ref_target = 1_000.0
        h_eos_ref = pyo.value(self.enthalpy_eos(delta_sat_liq_ref, tau_ref))
        s_eos_ref = pyo.value(self.entropy_eos(delta_sat_liq_ref, tau_ref))

        self.h_offset = h_eos_ref - h_ref_target  # additive enthalpy reference offset [J/kg]
        self.s_offset = s_eos_ref - s_ref_target  # additive entropy reference offset [J/(kg K)]

        # ========== Available property names ==========
        self._property_list = [
            "pressure",
            "enthalpy",
            "entropy",
            "heat_capacity_v",
            "heat_capacity_p",
            "density",
            "speed_of_sound",
            "compressibility_factor",
        ]

    @property
    def property_list(self):
        """Return the thermodynamic properties exposed by this parameter block."""
        return self._property_list

    # ======================================================================
    # Ideal dimensionless Helmholtz energy: IDAES Type 01
    # ======================================================================

    def alpha_ideal(self, delta, tau):
        """Return the ideal-gas contribution ``alpha0(delta, tau)``.

        The implemented expression is

        ``ln(delta) + n0[1] + n0[2]*tau + n0[3]*ln(tau)``
        ``+ sum(n0[i]*ln(1 - exp(-g0[i]*tau)))``.

        This is the same algebraic structure as IDAES ideal Type 01.
        """
        return (
            pyo.log(delta)
            + self.n0[1]
            + self.n0[2] * tau
            + self.n0[3] * pyo.log(tau)
            + sum(
                self.n0[i] * pyo.log(1.0 - pyo.exp(-self.g0[i] * tau))
                for i in self.planck_einstein_index
            )
        )

    def alpha0_delta(self, delta, tau):
        """Return ``d(alpha0)/d(delta)`` for the Type-01 ideal term."""
        del tau
        return 1.0 / delta

    def alpha0_delta_delta(self, delta, tau):
        """Return ``d2(alpha0)/d(delta)2`` for the Type-01 ideal term."""
        del tau
        return -1.0 / delta**2

    def alpha0_tau(self, delta, tau):
        """Return ``d(alpha0)/d(tau)`` for the Type-01 ideal term."""
        del delta
        return (
            self.n0[2]
            + self.n0[3] / tau
            + sum(
                self.n0[i]
                * self.g0[i]
                / (pyo.exp(self.g0[i] * tau) - 1.0)
                for i in self.planck_einstein_index
            )
        )

    def alpha0_tau_tau(self, delta, tau):
        """Return ``d2(alpha0)/d(tau)2`` for the Type-01 ideal term."""
        del delta
        return (
            -self.n0[3] / tau**2
            - sum(
                self.n0[i]
                * self.g0[i] ** 2
                * pyo.exp(-self.g0[i] * tau)
                / (1.0 - pyo.exp(-self.g0[i] * tau)) ** 2
                for i in self.planck_einstein_index
            )
        )

    def alpha0_delta_tau(self, delta, tau):
        """Return the mixed ideal derivative ``d2(alpha0)/(d(delta)d(tau))``."""
        del delta, tau
        return 0.0

    # ======================================================================
    # Residual dimensionless Helmholtz energy: IDAES Type 02
    # ======================================================================

    def alpha_residual(self, delta, tau):
        """Return the residual contribution ``alphar(delta, tau)``.

        The residual R1234yf EOS contains three groups:

        1. undamped power terms,
        2. density-exponentially damped power terms, and
        3. two-dimensional Gaussian terms.
        """
        polynomial = sum(
            self.n_res[i] * delta ** self.d_res[i] * tau ** self.t_res[i]
            for i in self.residual_polynomial_index
        )
        damped = sum(
            self.n_res[i]
            * delta ** self.d_res[i]
            * tau ** self.t_res[i]
            * pyo.exp(-(delta ** self.c_res[i]))
            for i in self.residual_damped_index
        )
        gaussian = sum(
            self.n_res[i]
            * delta ** self.d_res[i]
            * tau ** self.t_res[i]
            * pyo.exp(
                -self.a_res[i] * (delta - self.e_res[i]) ** 2
                - self.b_res[i] * (tau - self.g_res[i]) ** 2
            )
            for i in self.residual_gaussian_index
        )
        return polynomial + damped + gaussian

    def alphar_delta(self, delta, tau):
        """Return ``d(alphar)/d(delta)`` for residual Type 02."""
        polynomial = sum(
            self.n_res[i]
            * self.d_res[i]
            * delta ** (self.d_res[i] - 1.0)
            * tau ** self.t_res[i]
            for i in self.residual_polynomial_index
        )
        damped = sum(
            self.n_res[i]
            * pyo.exp(-(delta ** self.c_res[i]))
            * delta ** (self.d_res[i] - 1.0)
            * tau ** self.t_res[i]
            * (self.d_res[i] - self.c_res[i] * delta ** self.c_res[i])
            for i in self.residual_damped_index
        )
        gaussian = sum(
            self.n_res[i]
            * delta ** self.d_res[i]
            * tau ** self.t_res[i]
            * pyo.exp(
                -self.a_res[i] * (delta - self.e_res[i]) ** 2
                - self.b_res[i] * (tau - self.g_res[i]) ** 2
            )
            * (
                self.d_res[i] / delta
                - 2.0 * self.a_res[i] * (delta - self.e_res[i])
            )
            for i in self.residual_gaussian_index
        )
        return polynomial + damped + gaussian

    def alphar_delta_delta(self, delta, tau):
        """Return ``d2(alphar)/d(delta)2`` for residual Type 02."""
        polynomial = sum(
            self.n_res[i]
            * self.d_res[i]
            * (self.d_res[i] - 1.0)
            * delta ** (self.d_res[i] - 2.0)
            * tau ** self.t_res[i]
            for i in self.residual_polynomial_index
        )
        damped = sum(
            self.n_res[i]
            * pyo.exp(-(delta ** self.c_res[i]))
            * delta ** (self.d_res[i] - 2.0)
            * tau ** self.t_res[i]
            * (
                (
                    self.d_res[i]
                    - self.c_res[i] * delta ** self.c_res[i]
                )
                * (
                    self.d_res[i]
                    - 1.0
                    - self.c_res[i] * delta ** self.c_res[i]
                )
                - self.c_res[i] ** 2 * delta ** self.c_res[i]
            )
            for i in self.residual_damped_index
        )
        gaussian = sum(
            self.n_res[i]
            * tau ** self.t_res[i]
            * pyo.exp(
                -self.a_res[i] * (delta - self.e_res[i]) ** 2
                - self.b_res[i] * (tau - self.g_res[i]) ** 2
            )
            * (
                -2.0 * self.a_res[i] * delta ** self.d_res[i]
                + 4.0
                * self.a_res[i] ** 2
                * delta ** self.d_res[i]
                * (delta - self.e_res[i]) ** 2
                - 4.0
                * self.d_res[i]
                * self.a_res[i]
                * delta ** (self.d_res[i] - 1.0)
                * (delta - self.e_res[i])
                + self.d_res[i]
                * (self.d_res[i] - 1.0)
                * delta ** (self.d_res[i] - 2.0)
            )
            for i in self.residual_gaussian_index
        )
        return polynomial + damped + gaussian

    def alphar_tau(self, delta, tau):
        """Return ``d(alphar)/d(tau)`` for residual Type 02."""
        polynomial = sum(
            self.n_res[i]
            * self.t_res[i]
            * delta ** self.d_res[i]
            * tau ** (self.t_res[i] - 1.0)
            for i in self.residual_polynomial_index
        )
        damped = sum(
            self.n_res[i]
            * self.t_res[i]
            * delta ** self.d_res[i]
            * tau ** (self.t_res[i] - 1.0)
            * pyo.exp(-(delta ** self.c_res[i]))
            for i in self.residual_damped_index
        )
        gaussian = sum(
            self.n_res[i]
            * delta ** self.d_res[i]
            * tau ** self.t_res[i]
            * pyo.exp(
                -self.a_res[i] * (delta - self.e_res[i]) ** 2
                - self.b_res[i] * (tau - self.g_res[i]) ** 2
            )
            * (
                self.t_res[i] / tau
                - 2.0 * self.b_res[i] * (tau - self.g_res[i])
            )
            for i in self.residual_gaussian_index
        )
        return polynomial + damped + gaussian

    def alphar_tau_tau(self, delta, tau):
        """Return ``d2(alphar)/d(tau)2`` for residual Type 02."""
        polynomial = sum(
            self.n_res[i]
            * self.t_res[i]
            * (self.t_res[i] - 1.0)
            * delta ** self.d_res[i]
            * tau ** (self.t_res[i] - 2.0)
            for i in self.residual_polynomial_index
        )
        damped = sum(
            self.n_res[i]
            * self.t_res[i]
            * (self.t_res[i] - 1.0)
            * delta ** self.d_res[i]
            * tau ** (self.t_res[i] - 2.0)
            * pyo.exp(-(delta ** self.c_res[i]))
            for i in self.residual_damped_index
        )
        gaussian = sum(
            self.n_res[i]
            * delta ** self.d_res[i]
            * tau ** self.t_res[i]
            * pyo.exp(
                -self.a_res[i] * (delta - self.e_res[i]) ** 2
                - self.b_res[i] * (tau - self.g_res[i]) ** 2
            )
            * (
                (
                    self.t_res[i] / tau
                    - 2.0 * self.b_res[i] * (tau - self.g_res[i])
                )
                ** 2
                - self.t_res[i] / tau**2
                - 2.0 * self.b_res[i]
            )
            for i in self.residual_gaussian_index
        )
        return polynomial + damped + gaussian

    def alphar_delta_tau(self, delta, tau):
        """Return the mixed residual derivative ``d2(alphar)/(d(delta)d(tau))``."""
        polynomial = sum(
            self.n_res[i]
            * self.d_res[i]
            * self.t_res[i]
            * delta ** (self.d_res[i] - 1.0)
            * tau ** (self.t_res[i] - 1.0)
            for i in self.residual_polynomial_index
        )
        damped = sum(
            self.n_res[i]
            * self.t_res[i]
            * delta ** (self.d_res[i] - 1.0)
            * tau ** (self.t_res[i] - 1.0)
            * (
                self.d_res[i]
                - self.c_res[i] * delta ** self.c_res[i]
            )
            * pyo.exp(-(delta ** self.c_res[i]))
            for i in self.residual_damped_index
        )
        gaussian = sum(
            self.n_res[i]
            * delta ** self.d_res[i]
            * tau ** self.t_res[i]
            * pyo.exp(
                -self.a_res[i] * (delta - self.e_res[i]) ** 2
                - self.b_res[i] * (tau - self.g_res[i]) ** 2
            )
            * (
                self.d_res[i] / delta
                - 2.0 * self.a_res[i] * (delta - self.e_res[i])
            )
            * (
                self.t_res[i] / tau
                - 2.0 * self.b_res[i] * (tau - self.g_res[i])
            )
            for i in self.residual_gaussian_index
        )
        return polynomial + damped + gaussian

    def alpha_total(self, delta, tau):
        """Return total dimensionless Helmholtz free energy ``alpha0 + alphar``."""
        return self.alpha_ideal(delta, tau) + self.alpha_residual(delta, tau)

    # ======================================================================
    # Thermodynamic properties from Helmholtz derivatives
    # ======================================================================

    def pressure(self, delta, tau):
        """Return absolute pressure in Pa.

        The standard Helmholtz relation is

        ``p = rho*R*T*(1 + delta*alphar_delta)``.
        """
        return (
            delta
            * self.rhoc
            * self.R_gas
            * self.Tc
            / tau
            * (1.0 + delta * self.alphar_delta(delta, tau))
        )

    def enthalpy_eos(self, delta, tau):
        """Return raw EOS specific enthalpy before the reference-state offset.

        The standard reduced-Helmholtz relation is

        ``h/(R*T) = 1 + tau*alpha_tau + delta*alphar_delta``.
        """
        alpha_tau = self.alpha0_tau(delta, tau) + self.alphar_tau(delta, tau)
        return (
            self.R_gas
            * self.Tc
            / tau
            * (
                1.0
                + tau * alpha_tau
                + delta * self.alphar_delta(delta, tau)
            )
        )

    def enthalpy(self, delta, tau):
        """Return reference-adjusted specific enthalpy in J/kg."""
        return self.enthalpy_eos(delta, tau) - self.h_offset

    def entropy_eos(self, delta, tau):
        """Return raw EOS specific entropy before the reference-state offset.

        The standard reduced-Helmholtz relation is

        ``s/R = tau*alpha_tau - alpha``.
        """
        alpha_tau = self.alpha0_tau(delta, tau) + self.alphar_tau(delta, tau)
        return self.R_gas * (
            tau * alpha_tau - self.alpha_total(delta, tau)
        )

    def entropy(self, delta, tau):
        """Return reference-adjusted specific entropy in J/(kg K)."""
        return self.entropy_eos(delta, tau) - self.s_offset

    def internal_energy_eos(self, delta, tau):
        """Return raw EOS specific internal energy in J/kg.

        The Helmholtz relation is ``u/(R*T) = tau*alpha_tau``. Because
        ``T = Tc/tau``, this is also ``u = R*Tc*alpha_tau``.
        """
        alpha_tau = self.alpha0_tau(delta, tau) + self.alphar_tau(delta, tau)
        return self.R_gas * self.Tc * alpha_tau

    def heat_capacity_v(self, delta, tau):
        """Return isochoric heat capacity ``cv`` in J/(kg K)."""
        alpha_tt = self.alpha0_tau_tau(delta, tau) + self.alphar_tau_tau(
            delta, tau
        )
        return -self.R_gas * tau**2 * alpha_tt

    def heat_capacity_p(self, delta, tau):
        """Return isobaric heat capacity ``cp`` in J/(kg K).

        The expression is the standard single-phase Helmholtz EOS identity

        ``cp = cv + R*A^2/B``, where

        ``A = 1 + delta*alphar_delta - delta*tau*alphar_delta_tau``

        and

        ``B = 1 + 2*delta*alphar_delta + delta^2*alphar_delta_delta``.

        Near the critical point or a spinodal, ``B`` approaches zero and the
        property becomes numerically ill-conditioned, as expected physically.
        """
        ar_d = self.alphar_delta(delta, tau)
        ar_dd = self.alphar_delta_delta(delta, tau)
        ar_dt = self.alphar_delta_tau(delta, tau)
        cv = self.heat_capacity_v(delta, tau)

        numerator = 1.0 + delta * ar_d - delta * tau * ar_dt
        denominator = 1.0 + 2.0 * delta * ar_d + delta**2 * ar_dd
        return cv + self.R_gas * numerator**2 / denominator

    def speed_of_sound(self, delta, tau):
        """Return the thermodynamic speed of sound in m/s.

        The implemented single-phase Helmholtz relation is

        ``w^2 = R*T*[B - A^2/(tau^2*alpha_tt)]``,

        with ``A`` and ``B`` defined as in :meth:`heat_capacity_p` and
        ``alpha_tt = alpha0_tt + alphar_tt``.

        The expression is singular or ill-conditioned on stability limits and
        in the immediate critical region. It should not be evaluated inside
        the two-phase region as a homogeneous single-phase property.
        """
        ar_d = self.alphar_delta(delta, tau)
        ar_dd = self.alphar_delta_delta(delta, tau)
        ar_dt = self.alphar_delta_tau(delta, tau)
        alpha_tt = self.alpha0_tau_tau(delta, tau) + self.alphar_tau_tau(
            delta, tau
        )

        numerator = 1.0 + delta * ar_d - delta * tau * ar_dt
        mechanical = 1.0 + 2.0 * delta * ar_d + delta**2 * ar_dd
        thermal_correction = numerator**2 / (tau**2 * alpha_tt)
        w_squared = (
            self.R_gas
            * self.Tc
            / tau
            * (mechanical - thermal_correction)
        )
        return pyo.sqrt(w_squared)

    def density(self, delta):
        """Return mass density ``rho = delta*rhoc`` in kg/m^3."""
        return delta * self.rhoc

    def compressibility_factor(self, delta, tau):
        """Return compressibility factor ``Z = p/(rho*R*T)``."""
        return 1.0 + delta * self.alphar_delta(delta, tau)


################################################################################
# R1234yf Property State Block Class
################################################################################


class R1234yfPropertyStateBlock(pyo.Block):
    """Helper state block for solving an R1234yf state from pressure and enthalpy.

    This class preserves the structure of the original development file: it is
    intended to be nested directly inside an ``R1234yfPropertyParameterBlock``
    so that ``self.parent_block()`` is the parameter block. It is a Pyomo helper
    block, not a replacement for IDAES ``StateBlockData`` infrastructure.

    Primary state variables are pressure ``P`` and specific enthalpy ``H``.
    Temperature, density, reduced density, and inverse reduced temperature are
    auxiliary unknowns constrained by the Helmholtz EOS.
    """

    def build(self):
        """Construct state variables, EOS constraints, and property expressions."""
        params = self.parent_block()
        if not hasattr(params, "Tc") or not hasattr(params, "pressure"):
            raise RuntimeError(
                "R1234yfPropertyStateBlock must be nested directly inside an "
                "R1234yfPropertyParameterBlock."
            )

        # ========== Primary state variables ==========
        self.P = pyo.Var(
            initialize=1.0e5,
            bounds=(1.0e3, 1.0e8),
            doc="Absolute pressure [Pa]",
        )
        self.H = pyo.Var(
            initialize=3.0e5,
            bounds=(-1.0e5, 1.5e6),
            doc="Reference-adjusted specific enthalpy [J/kg]",
        )

        # ========== Auxiliary state variables ==========
        self.T = pyo.Var(
            initialize=300.0,
            bounds=(params.Tc * 0.30, params.Tc * 1.50),
            doc="Temperature [K]",
        )
        self.rho = pyo.Var(
            initialize=100.0,
            bounds=(1.0e-6, params.rhoc * 5.0),
            doc="Mass density [kg/m^3]",
        )
        self.delta = pyo.Var(
            initialize=0.5,
            bounds=(1.0e-10, 5.0),
            doc="Reduced density rho/rhoc [-]",
        )
        self.tau = pyo.Var(
            initialize=1.2,
            bounds=(0.3, 10.0),
            doc="Inverse reduced temperature Tc/T [-]",
        )

        def _delta_constraint(block):
            """Enforce the definition ``delta = rho/rhoc``."""
            return block.delta == block.rho / params.rhoc

        self.delta_constraint = pyo.Constraint(
            rule=_delta_constraint,
            doc="Reduced-density definition",
        )

        def _tau_constraint(block):
            """Enforce the definition ``tau = Tc/T``."""
            return block.tau == params.Tc / block.T

        self.tau_constraint = pyo.Constraint(
            rule=_tau_constraint,
            doc="Inverse reduced-temperature definition",
        )

        def _pressure_constraint(block):
            """Enforce pressure consistency with the Helmholtz EOS."""
            return block.P == params.pressure(block.delta, block.tau)

        self.pressure_constraint = pyo.Constraint(
            rule=_pressure_constraint,
            doc="Pressure from Helmholtz EOS",
        )

        def _enthalpy_constraint(block):
            """Enforce enthalpy consistency with the Helmholtz EOS."""
            return block.H == params.enthalpy(block.delta, block.tau)

        self.enthalpy_constraint = pyo.Constraint(
            rule=_enthalpy_constraint,
            doc="Specific enthalpy from Helmholtz EOS",
        )

        def _entropy(block):
            """Return specific entropy from the converged Helmholtz state."""
            return params.entropy(block.delta, block.tau)

        self.S = pyo.Expression(
            rule=_entropy,
            doc="Specific entropy [J/(kg K)]",
        )

        def _heat_capacity_v(block):
            """Return isochoric heat capacity from the Helmholtz EOS."""
            return params.heat_capacity_v(block.delta, block.tau)

        self.Cv = pyo.Expression(
            rule=_heat_capacity_v,
            doc="Isochoric heat capacity [J/(kg K)]",
        )

        def _heat_capacity_p(block):
            """Return isobaric heat capacity from the Helmholtz EOS."""
            return params.heat_capacity_p(block.delta, block.tau)

        self.Cp = pyo.Expression(
            rule=_heat_capacity_p,
            doc="Isobaric heat capacity [J/(kg K)]",
        )

        def _speed_of_sound(block):
            """Return the thermodynamic speed of sound from the Helmholtz EOS."""
            return params.speed_of_sound(block.delta, block.tau)

        self.W = pyo.Expression(
            rule=_speed_of_sound,
            doc="Speed of sound [m/s]",
        )

        def _compressibility_factor(block):
            """Return compressibility factor ``Z`` from the Helmholtz EOS."""
            return params.compressibility_factor(block.delta, block.tau)

        self.Z = pyo.Expression(
            rule=_compressibility_factor,
            doc="Compressibility factor [-]",
        )

        def _ideal_gas_temperature_diagnostic(block):
            """Return ``P/(rho*R) = Z*T`` as an ideal-gas-limit diagnostic.

            This quantity is not expected to equal the actual temperature for a
            dense real fluid. It approaches ``T`` only when ``Z`` approaches 1.
            """
            return block.P / (block.rho * params.R_gas)

        self.T_from_ideal_gas = pyo.Expression(
            rule=_ideal_gas_temperature_diagnostic,
            doc="Ideal-gas-limit diagnostic P/(rho R) = Z*T [K]",
        )


################################################################################
# Convenience factory function
################################################################################


def create_r1234yf_property_package():
    """Create and return an R1234yf parameter block component.

    The returned object must still be attached to a Pyomo model before normal
    component construction. To preserve the parent lookup used by the helper
    state block, add any ``R1234yfPropertyStateBlock`` as a child of the
    parameter block.

    Examples
    --------
    ``model.params = create_r1234yf_property_package()``

    ``model.params.state = R1234yfPropertyStateBlock()``
    """
    return R1234yfPropertyParameterBlock()
