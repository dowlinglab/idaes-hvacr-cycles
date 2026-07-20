"""
This code generates the p-H diagram based on the Peng-Robinson equation
of state per the IDAES property package cubic EOS convention.

The model of the PR equation of state used is based on:
https://idaes-pse.readthedocs.io/en/1.8.0/user_guide/components/property_package/general/eos/cubic.html

The Shomate expressions are based on the NIST WebBook.

Author: Shilpa Narasimhan
Date Created: 07/20/2026
Support from: Claude AI, ChatGPT, and Codex (Modifications for readability)
QA/testing: Shilpa Narasimhan

    This script estimates a saturated p-H diagram for R-32 using:
        1. Ambrose-Walton saturation pressure correlation.
        2. Peng-Robinson cubic EOS departure enthalpy.
        3. Shomate ideal-gas heat-capacity integration.
        4. IIR enthalpy reference: saturated liquid at 0 deg C has h = 200 kJ/kg.
"""

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Pure-component constants
# =============================================================================

R = 8.314  # J/mol/K, universal gas constant

Tc = 351.3  # K
Pc = 57.82e5  # Pa
MW = 52.023  # kg/kmol, numerically equivalent to g/mol

omega_ref = 0.2769  # REFPROP acentric factor reference


# =============================================================================
# VLE data from Linde datasheet
# =============================================================================

linde = [
    (-130, 0.001312),
    (-110, 0.014525),
    (-90, 0.07556),
    (-70, 0.36067),
    (-50, 1.014),
    (-30, 2.7344),
    (-10, 5.8263),
    (0, 8.131),
    (10, 10.065),
    (20, 14.746),
    (30, 19.275),
    (40, 24.783),
    (50, 31.412),
    (58, 37.635),
    (62, 41.089),
    (66, 44.793),
    (70, 48.768),
    (74, 53.046),
    (76, 55.315),
    (78, 57.697),
]


# =============================================================================
# Acentric factor estimate from Linde data near Tr = 0.7
# =============================================================================

T_r_target = 0.7
T_target = T_r_target * Tc  # K

# Tr = 0.7 corresponds to:
# T_target = 245.91 K = -27.24 deg C
#
# The Linde datasheet does not have a point exactly at -27.24 deg C.
# This pressure uses the nearby -28 deg C value noted by Shilpa.
P_target = 2.9675e5  # Pa

P_r_target = P_target / Pc
omega_act = -np.log10(P_r_target) - 1.0

# Choose which acentric factor to use.
# Use omega_act if intentionally matching the Linde-derived estimate.
# Use omega_ref if intentionally matching REFPROP.
omega = omega_act


# =============================================================================
# Saturation pressure from Ambrose-Walton
# =============================================================================

def saturation_pressure(T):
    """
    Return saturation pressure from the Ambrose-Walton correlation.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Saturation pressure in Pa.

    Notes
    -----
    The correlation uses critical properties and the acentric factor.
    """
    Tr = T / Tc
    tau = 1.0 - Tr

    f0 = (
        -5.97616 * tau
        + 1.29874 * tau**1.5
        - 0.60394 * tau**2.5
        - 1.06841 * tau**5
    ) / Tr

    f1 = (
        -5.03365 * tau
        + 1.11505 * tau**1.5
        - 5.41217 * tau**2.5
        - 7.46628 * tau**5
    ) / Tr

    f2 = (
        -0.64771 * tau
        + 2.41539 * tau**1.5
        - 4.26979 * tau**2.5
        + 3.25259 * tau**5
    ) / Tr

    return Pc * np.exp(f0 + omega * f1 + omega**2 * f2)


# =============================================================================
# Temperature grid and saturation pressure check
# =============================================================================

T_max = 78 + 273.15  # K
Tmin = -130.0 + 273.15  # K

T_grid = np.linspace(Tmin, T_max, 100)
P_sat = np.array([saturation_pressure(T) for T in T_grid])

# Saturation-pressure vs Linde check -- only when run directly, not on import.
if __name__ == "__main__":
    print(f"{'T[C]':>6}{'Linde[bar]':>12}{'Model[bar]':>12}{'error[%]':>10}")

    errors = []

    for T_C, P_linde in linde:
        P_model = saturation_pressure(T_C + 273.15) / 1e5
        err = 100.0 * (P_model - P_linde) / P_linde
        errors.append(err)
        print(f"{T_C:6}{P_linde:12.4f}{P_model:12.4f}{err:10.2f}")

    errors = np.array(errors)

    print(
        f"\nMAPE = {np.mean(np.abs(errors)):.2f}%"
        f"   max|err| = {np.max(np.abs(errors)):.2f}%"
    )


# =============================================================================
# Peng-Robinson EOS parameters
# =============================================================================

u = 2.0
w = -1.0

Omega_A = 0.45724
Omega_B = 0.0778

kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2

b = Omega_B * R * Tc / Pc


def alpha_pr(T):
    """
    Return Peng-Robinson alpha(T).

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Dimensionless alpha function.
    """
    Tr = T / Tc
    return (1.0 + kappa * (1.0 - np.sqrt(Tr)))**2


def a_param(T):
    """
    Return Peng-Robinson attraction parameter a(T).

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Peng-Robinson attraction parameter in Pa m^6/mol^2.
    """
    return Omega_A * R**2 * Tc**2 * alpha_pr(T) / Pc


def dadT(T):
    """
    Return derivative da/dT for Peng-Robinson a(T).

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Temperature derivative of a(T).
    """
    Tr = T / Tc

    return (
        Omega_A
        * R**2
        * Tc**2
        / Pc
        * (
            -kappa
            * (1.0 + kappa * (1.0 - np.sqrt(Tr)))
            / (Tc * np.sqrt(Tr))
        )
    )


# =============================================================================
# Shomate ideal-gas heat capacity
# =============================================================================

A_shomate = -6.098682
B_shomate = 179.22
C_shomate = -122.3682
D_shomate = 32.30207
E_shomate = 0.491361


def Cp_ideal(T):
    """
    Return ideal-gas heat capacity from the Shomate expression.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Ideal-gas heat capacity in J/mol/K.

    Notes
    -----
    The Shomate variable is t = T/1000.
    """
    t = T / 1000.0

    return (
        A_shomate
        + B_shomate * t
        + C_shomate * t**2
        + D_shomate * t**3
        + E_shomate / t**2
    )


def int_shomate(T):
    """
    Return the indefinite integral of the Shomate Cp expression.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Integrated ideal-gas enthalpy expression in J/mol.

    Notes
    -----
    The integration constant is omitted because only enthalpy differences
    are used.
    """
    return (
        A_shomate * T
        + 0.5e-3 * B_shomate * T**2
        + (1.0 / 3.0) * 1e-6 * C_shomate * T**3
        + 0.25e-9 * D_shomate * T**4
        - 1e6 * E_shomate / T
    )


def h_ideal(T):
    """
    Return ideal-gas enthalpy change from 273.15 K to T.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Ideal-gas enthalpy change in J/mol.

    Notes
    -----
    The 273.15 K reference is used because the final enthalpy offset is set
    to the IIR convention: saturated liquid at 0 deg C has h = 200 kJ/kg.
    """
    return int_shomate(T) - int_shomate(273.15)


# =============================================================================
# PR compressibility roots and departure enthalpy
# =============================================================================

def z_roots(T, P):
    """
    Return real vapor/liquid compressibility roots for the cubic EOS.

    Parameters
    ----------
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.

    Returns
    -------
    tuple
        roots, A, B where:
            roots : list[float]
                Real compressibility roots greater than B.
            A : float
                Dimensionless PR attraction parameter.
            B : float
                Dimensionless PR co-volume parameter.
    """
    a = a_param(T)

    A = a * P / (R**2 * T**2)
    B = b * P / (R * T)

    coefficients = [
        1.0,
        -(1.0 + B - u * B),
        A - u * B - (u - w) * B**2,
        -(A * B + w * B**2 + w * B**3),
    ]

    roots_raw = np.roots(coefficients)

    roots = sorted(
        root.real
        for root in roots_raw
        if abs(root.imag) < 1e-9 and root.real > B
    )

    if not roots:
        raise ValueError(
            f"No valid real Z roots found at T = {T:.6g} K, P = {P:.6g} Pa."
        )

    return roots, A, B


def dh_dep(T, P, phase):
    """
    Return Peng-Robinson departure enthalpy.

    Parameters
    ----------
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    phase : str
        Either "liquid" or "vapor".

    Returns
    -------
    float
        Departure enthalpy in J/mol.
    """
    a = a_param(T)
    da = dadT(T)

    roots, A, B = z_roots(T, P)

    if phase == "vapor":
        Z = roots[-1]
    elif phase == "liquid":
        Z = roots[0]
    else:
        raise ValueError("phase must be either 'liquid' or 'vapor'.")

    s = np.sqrt(u**2 - 4.0 * w)

    log_argument = (
        (2.0 * Z + B * (u + s))
        / (2.0 * Z + B * (u - s))
    )

    return (
        R * T * (Z - 1.0)
        + (T * da - a)
        / (b * s)
        * np.log(log_argument)
    )


def h_mass(T, P, phase):
    """
    Return total specific enthalpy from ideal and departure contributions.

    Parameters
    ----------
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    phase : str
        Either "liquid" or "vapor".

    Returns
    -------
    float
        Specific enthalpy in kJ/kg.

    Notes
    -----
    The numerator is J/mol. Dividing by MW in kg/kmol gives kJ/kg because
    MW is numerically equal to g/mol.
    """
    return (h_ideal(T) + dh_dep(T, P, phase)) / MW


# =============================================================================
# Shomate ideal-gas entropy integral
# =============================================================================

def int_shomate_S(T):
    """
    Return the integral of Cp_ideal/T dT for the Shomate expression.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Integrated ideal-gas entropy expression in J/mol/K.

    Notes
    -----
    The Shomate variable is t = T/1000. Because dT/T = dt/t, the integral is
    A*ln(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2). The integration constant
    is omitted because only entropy differences are used.
    """
    t = T / 1000.0

    return (
        A_shomate * np.log(t)
        + B_shomate * t
        + C_shomate * t**2 / 2.0
        + D_shomate * t**3 / 3.0
        - E_shomate / (2.0 * t**2)
    )


def s_ideal(T, P):
    """
    Return ideal-gas entropy change relative to 273.15 K and 1 bar.

    Parameters
    ----------
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.

    Returns
    -------
    float
        Ideal-gas entropy change in J/mol/K.

    Notes
    -----
    Unlike enthalpy, ideal-gas entropy depends on pressure through the
    -R*ln(P/P_ref) term. P_ref = 1 bar; that choice only shifts a constant,
    which is absorbed by the IIR entropy offset.
    """
    return (int_shomate_S(T) - int_shomate_S(273.15)) - R * np.log(P / 1e5)


# =============================================================================
# PR departure entropy
# =============================================================================

def ds_dep(T, P, phase):
    """
    Return Peng-Robinson departure entropy.

    Parameters
    ----------
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    phase : str
        Either "liquid" or "vapor".

    Returns
    -------
    float
        Departure entropy in J/mol/K.

    Notes
    -----
    Alignment with Poling Table 6-3 (entropy row). Table 6-3 lists the departure
    (using the sign convention F_d = F_ig - F) as:

        (S_ig - S)/R = [(dTheta/dT) / (R*sqrt(delta**2 - 4*eps))] * ln(RATIO_V)
                       - ln[Z*(1 - b/V)]

    with the PR mapping  Theta = a,  dTheta/dT = da/dT,  delta = u*b,
    eps = w*b**2, so sqrt(delta**2 - 4*eps) = b*sqrt(u**2 - 4*w) = b*s, and
    RATIO_V = (2V + delta - b*s) / (2V + delta + b*s).

    Steps to reach the expression coded below:
        1. Multiply by R and negate to convert (S_ig - S) -> (S - S_ig).
        2. Negating flips the log ratio to (2V + delta + b*s)/(2V + delta - b*s).
        3. Convert V -> Z using V = Z*R*T/P and b = B*R*T/P; the R*T/P factor
           cancels in the ratio, giving (2Z + B*(u + s)) / (2Z + B*(u - s)).
        4. Z*(1 - b/V) = Z - B, so -ln[Z*(1 - b/V)] -> +R*ln(Z - B).

    Result:
        S - S_ig = R*ln(Z - B) + (da/dT)/(b*s) * ln[(2Z + B(u+s))/(2Z + B(u-s))]

    which is coded below. The -R*ln(P/P_ref) pressure term is carried in
    s_ideal, consistent with Table 6-1's F_d = F_ig(T, P) - F(T, P) definition
    (ideal gas evaluated at the same T and P).
    """
    da = dadT(T)

    roots, A, B = z_roots(T, P)

    if phase == "vapor":
        Z = roots[-1]
    elif phase == "liquid":
        Z = roots[0]
    else:
        raise ValueError("phase must be either 'liquid' or 'vapor'.")

    s = np.sqrt(u**2 - 4.0 * w)

    log_argument = (
        (2.0 * Z + B * (u + s))
        / (2.0 * Z + B * (u - s))
    )

    return (
        R * np.log(Z - B)
        + da
        / (b * s)
        * np.log(log_argument)
    )


def s_mass(T, P, phase):
    """
    Return total specific entropy from ideal and departure contributions.

    Parameters
    ----------
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    phase : str
        Either "liquid" or "vapor".

    Returns
    -------
    float
        Specific entropy in kJ/kg/K.

    Notes
    -----
    The numerator is J/mol/K. Dividing by MW in kg/kmol gives kJ/kg/K because
    MW is numerically equal to g/mol.
    """
    return (s_ideal(T, P) + ds_dep(T, P, phase)) / MW


# =============================================================================
# Apply IIR enthalpy reference and generate saturated p-H data
# =============================================================================

T0 = 273.15
P0 = saturation_pressure(T0)

h_offset = 200.0 - h_mass(T0, P0, "liquid")

hf = np.array(
    [
        h_mass(T, saturation_pressure(T), "liquid") + h_offset
        for T in T_grid
    ]
)

hg = np.array(
    [
        h_mass(T, saturation_pressure(T), "vapor") + h_offset
        for T in T_grid
    ]
)

Pd = P_sat / 1e5


# =============================================================================
# Apply IIR entropy reference and generate saturated T-s data
# =============================================================================

s_offset = 1.0 - s_mass(T0, P0, "liquid")

sf = np.array(
    [
        s_mass(T, saturation_pressure(T), "liquid") + s_offset
        for T in T_grid
    ]
)

sg = np.array(
    [
        s_mass(T, saturation_pressure(T), "vapor") + s_offset
        for T in T_grid
    ]
)

Ts = T_grid - 273.15  # deg C

# =============================================================================
# Plots (only when run directly, not on import)
# =============================================================================
# Guarded under __main__ so this module can be imported as a library
# (e.g. by phase_1_cubic_eos_validation.py) without popping up the diagrams
# or blocking on plt.show(). Run `python pr_eos_lib.py` to see the plots.

if __name__ == "__main__":

    # ---- p-H diagram ----
    plt.figure()
    plt.plot(hf, Pd, "tab:blue", label="saturated liquid")
    plt.plot(hg, Pd, "tab:red", label="saturated vapor")
    plt.yscale("log")
    plt.xlabel("h [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title("R-32 p-h diagram (Peng-Robinson)")
    plt.legend()
    plt.grid(True, which="both", ls=":")
    plt.tight_layout()
    plt.show()

    # ---- T-s diagram ----
    plt.figure()
    plt.plot(sf, Ts, "tab:blue", label="saturated liquid")
    plt.plot(sg, Ts, "tab:red", label="saturated vapor")
    plt.xlabel("s [kJ/kg/K]")
    plt.ylabel("T [deg C]")
    plt.title("R-32 T-s diagram (Peng-Robinson)")
    plt.legend()
    plt.grid(True, ls=":")
    plt.tight_layout()
    plt.show()