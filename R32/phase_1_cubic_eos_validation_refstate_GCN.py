"""
COPY for reference-state isolation testing -- see BREADCRUMB_07-20.md,
2026-08-10 update. Identical to phase_1_cubic_eos_validation.py except:
METHODS now carries per-method "F" (kJ/mol) and "G" (J/mol/K) Shomate
offsets, calibrated so each method's own predicted saturated liquid at
0 C lands on the IIR reference state (h=200 kJ/kg, s=1.00 kJ/kg-K) --
the same convention Helmholtz/CoolProp already use. make_config() reads
F/G from p["F"]/p["G"] instead of the hardcoded 0.0. Original file is
untouched.

This code serves as a validation to check if IDAES-generated cubic equations of
state properties align with the vanilla python code PR_EOS
 for all methods within the Colon group collaboration.
This file calls pr_eos_lib.py which is a version of the vanilla python code created
for this project
Date created: 07/20/2026

Author: Shilpa Narasimhan

Support: Claude AI

--- 08/26/2026 ---
This is the working copy for registering the Colon group's new
collaborator-supplied Shomate fits. Changes made to METHODS below:
  - SPGP commented out (not deleted), per "comment out old SPGP."
  - "first_principle" and "gcn" added -- Shomate-style Cp coefficients
    from spgp_r32.xlsx's a-e columns, rescaled from raw-T form to this
    file's t = T/1000 Shomate convention (A=a, B=1000b, C=1e6*c,
    D=1e9*d, E=e/1e6).
  - Both new methods' omega is computed from the Colon group's own
    "pvap" spreadsheet column, read as Pa (not mmHg as labeled -- that
    literal reading gives Psat > Pc for both methods, physically
    impossible) and taken as Psat at Tr=0.7*Tc, the standard Pitzer
    input: omega = -1 - log10(Psat/Pc). Values: first_principle =
    0.191312, gcn = 0.667570 (replacing earlier Linde-derived
    placeholders of -0.2385/0.6234). Tr=0.7 is an assumption, not yet
    confirmed with the collaborator.
  - F/G set to 0.0 for every method (not just the new ones) -- confirmed
    algebraically COP-invariant, see phase_1_cubic_eos_validation_
    refstate.py's 08/26 note for the full reasoning. NIST's F/G updated
    from its old calibrated 25.9022/84.9031 to 0.0/0.0 for consistency
    with GCGP/SPGP(commented)/first_principle/gcn.
  - Even with omega/F-G resolved, first_principle/gcn still show large
    errors vs Linde in compare_cp_methods_GCN.py (65.99%/46.72% pressure
    MAPE) -- likely a genuine Tc/Pc data-quality issue, not a units or
    omega-computation bug. Not yet resolved.
"""

## Importing packages

import numpy as np
from pyomo.environ import units as pyunits ## To enable IDAES to attach physical units

from idaes.core import Component, LiquidPhase, VaporPhase # To declare the component phase
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
## Generic parameter block reads the configuration file provided as dict and builds all the
## thermodynamic property knowledge that will be used within IDAES
from idaes.models.properties.modular_properties.state_definitions import FTPx
## Tell IDAES which variables describe the thermodynamic state of each stream       /point
# F: Molar flow rate, T: Temperature, P = Pressure, x = mole fraction
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
## Type of equations of state, here we use cubic equations of state: Cubic
## CubicType: which cubic equation, e.g., SRK/PR
from idaes.models.properties.modular_properties.pure import NIST
# For IDAES built-in expressions for Shomate 
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
#the method that lets the solver cross those edges smoothly —
#  going from all-liquid → two-phase → all-vapor — 
# without a hard jump when a phase appears or disappears.
#  "Smooth" so ipopt doesn't choke at the transition.
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew
#finds the edges of the two-phase region: the bubble point (where the first vapor bubble
#  appears) and the dew point (where the first liquid drop appears). For a pure fluid
#  like R-32, both are just the saturation point. It's where boiling/condensing happens.
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
# the rule for equilibrium: liquid and vapor coexist when their fugacities are equal
#  (ln f_liquid = ln f_vapor). 


##############################################################################################
#           Defining all the constants
##############################################################################################
R = 8.314 # Real gas constant, J/mol/K
MW = 52.024e-03 # kg/mol Molar mass of CH₂F₂, from standard atomic weights
#(NIST WebBook: 52.024 g/mol).
# Antoine: log10(P_bar) = A - B/(T_K + C)  -- saturation-pressure init guess, fit to Linde
Antoine_A = 4.60123
Antoine_B = 959.89766   # K
Antoine_C = -13.71589   # K
## Defining all the properties
# "omega" = Pitzer acentric factor: omega = -1 - log10(Psat(Tr=0.7) / Pc).
# NIST/GCGP/SPGP: Psat(Tr=0.7) from Linde's own vapor-pressure table, using
# each method's own Tc and Pc. Values: NIST 0.2769, GCGP 0.1711, SPGP -0.2741.
# (Negative for SPGP because its Tc is far too high -> a symptom of bad
#  critical properties, not a real acentric factor.)
# first_principle/gcn: Psat(Tr=0.7) instead comes directly from the Colon
# group's own "pvap" spreadsheet column (mislabeled mmHg, actually Pa) --
# same formula, different Psat source. Values: first_principle 0.191312,
# gcn 0.667570.

METHODS = {
    "NIST": {"Pc": 57.82e5,  "Tc": 351.3,   "omega": 0.2769,
             "A": -6.098682, "B": 179.2200, "C": -122.3682, "D": 32.30207, "E": 0.491361,
             "F": 0.0, "G": 0.0},
    "GCGP": {"Pc": 50.730e5, "Tc": 355.354, "omega": 0.1711,
             "A": 14.161,    "B": 0.124,    "C": -6.340e-05, "D": 1.190e-8, "E": 0.0,
             "F": 0.0, "G": 0.0},
    # "SPGP": {"Pc": 50.8106e5,"Tc": 400.898, "omega": -0.2741,
    #          "A": 129.687,   "B": 171.303,  "C": 146.2,      "D": 61.9837,
    #          "E": -0.0000361638,
    #          "F": 0.0, "G": 0.0},
    "first_principle": {"Pc": 52.54216e5, "Tc": 397.6136, "omega": 0.191312,
             "A": 29.1173, "B": 139.601, "C": 62.8355, "D": 25.3128, "E": -0.000791509,
             "F": 0.0, "G": 0.0},
    "gcn": {"Pc": 62.22664e5, "Tc": 329.1245, "omega": 0.667570,
             "A": 165.655, "B": 104.894, "C": 26.4086, "D": -4.04051, "E": -0.000122942,
             "F": 0.0, "G": 0.0},
}
# F/G: set to 0.0 for every method (08/26/2026) -- NIST/GCGP/SPGP's F/G were
# never supplied externally either; they were back-calculated offline via
# calibrate_FG.py to anchor each method's saturated liquid at 0 C to the IIR
# reference state (h=200 kJ/kg, s=1.00 kJ/kg-K). Confirmed algebraically that
# zeroing F/G changes no COP number -- F/G are per-method constants added
# uniformly to every h/s a method computes, so they cancel out of every h/s
# DIFFERENCE (compressor work, evaporator/condenser heat) and out of the
# compressor's isentropic equality (s_out=s_in), which is all COP is built
# from. Only effect: absolute h/s values are no longer IIR-anchored.

##############################################################################################
#           Build the IDAES generic-property configuration for parameters p from methods
##############################################################################################


def make_config(p):
    ## Return which property methods + parameter values
    return{
        "components":{
            "R32":{
                "type": Component,
                "cp_mol_ig_comp": NIST, # Shomate ideal-gas Cp
                "enth_mol_ig_comp": NIST, # Shomate ideal-gas enthalpy
                "entr_mol_ig_comp": NIST, # Shomate ideal-gas entropy
                "pressure_sat_comp": NIST,  # Antoine saturation pressure (init guess)
                "phase_equilibrium_form": {("Vap","Liq"): log_fugacity}, ## VLE form
                "parameter_data": {
                    "mw": (MW, pyunits.kg/pyunits.mol),
                    "pressure_crit": (p["Pc"], pyunits.Pa),
                    "temperature_crit": (p["Tc"], pyunits.K),
                    "omega": p["omega"],
                    "cp_mol_ig_comp_coeff":{
                        "A": (p["A"], pyunits.J/pyunits.mol/pyunits.K),
                        "B": (p["B"], pyunits.J/pyunits.mol/pyunits.K/pyunits.kiloK),
                        "C": (p["C"], pyunits.J/pyunits.mol/pyunits.K/pyunits.kiloK**2),
                        "D": (p["D"], pyunits.J/pyunits.mol/pyunits.K/pyunits.kiloK**3),
                        "E": (p["E"], pyunits.J*pyunits.kiloK**2/pyunits.mol/pyunits.K),
                        "F": (p["F"], pyunits.kJ/pyunits.mol),
                        "G": (p["G"], pyunits.J/pyunits.mol/pyunits.K),
                        "H": (0.0, pyunits.kJ/pyunits.mol),
                    },
                    "pressure_sat_comp_coeff":{
                        "A": (Antoine_A, None),
                        "B": (Antoine_B, pyunits.K),
                        "C": (Antoine_C, pyunits.K),
                    },
                },
            }
        },
        "phases": {
            "Liq": {"type": LiquidPhase, "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR}},
            "Vap": {"type": VaporPhase, "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR}},
        },
        "base_units": {"time": pyunits.s, "length": pyunits.m, "mass": pyunits.kg,
                       "amount": pyunits.mol, "temperature": pyunits.K},
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0.0, 1.0, 1000.0, pyunits.mol/pyunits.s),
            "temperature": (200.0, 300.0, 450.0, pyunits.K),
            "pressure": (1.0e4, 1.0e5, 1.0e7, pyunits.Pa),
        },
        "pressure_ref": (1.0e5, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "include_enthalpy_of_formation": False,
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
        "bubble_dew_method": LogBubbleDew,
        "parameter_data": {"PR_kappa": {("R32", "R32"): 0.000}},
    }

##############################################################################################
#   Vanilla (hand-written) Peng-Robinson references, parameterized per method
##############################################################################################
# Same physics as pr_eos_lib.py, but taking Tc, Pc, omega, and the Shomate
# coefficients as arguments so ONE function works for NIST, GCGP, and SPGP.
# Verified against pr_eos_lib for NIST: Z_vap 0.876, Cp_ideal 42.455.

def vanilla_z_vap(T, P, Tc, Pc, omega):
    """Vapor-root compressibility Z from the PR cubic (largest real root)."""
    kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
    alpha = (1.0 + kappa*(1.0 - np.sqrt(T/Tc)))**2
    a = 0.45724 * R**2 * Tc**2 * alpha / Pc
    b = 0.07780 * R * Tc / Pc
    A = a*P/(R*T)**2
    B = b*P/(R*T)
    coeffs = [1.0, -(1.0 - B), A - 2.0*B - 3.0*B**2, -(A*B - B**2 - B**3)]
    roots = np.roots(coeffs)
    real = roots[np.abs(roots.imag) < 1e-8].real
    return max(real)

def vanilla_cp_ideal(T, A, B, C, D, E):
    """Shomate ideal-gas Cp [J/mol/K], t = T/1000."""
    t = T/1000.0
    return A + B*t + C*t**2 + D*t**3 + E/t**2


##############################################################################################
#   Phase-1 validation: IDAES cubic-EoS Z (and Cp) vs the vanilla PR code
##############################################################################################

if __name__ == "__main__":
    from pyomo.environ import ConcreteModel, value
    from idaes.core import FlowsheetBlock
    from idaes.core.util.model_statistics import degrees_of_freedom
    from idaes.core.solvers import get_solver
    import pr_eos_lib as pr   # NIST-only, hand-written reference (the "vanilla" code)
    import idaes.logger as idaeslog

    T_test, P_test = 293.15, 10e5   # 20 deg C, 10 bar

    print(f"\n{'method':>6}{'Z_idaes':>10}{'Z_vanilla':>11}{'dZ':>9}"
          f"{'Cp_idaes':>11}{'Cp_ideal':>11}")
    print("-" * 58)

    z_idaes = {} # dictionary that stashes the IDAES-computed vapor Z for each method so you can use it after the loop.
    z_van = {} # disctionary that stasges the vanilla python-based vapor Z for each method so you can use it after the loop
    for name, p in METHODS.items():
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = GenericParameterBlock(**make_config(p))
        m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
        state_r32 = m.fs.state[0]
        state_r32.flow_mol.fix(1.0)
        state_r32.mole_frac_comp["R32"].fix(1.0)
        state_r32.temperature.fix(T_test)
        state_r32.pressure.fix(P_test)
        assert degrees_of_freedom(m) == 0, f"{name}: DOF != 0"
        m.fs.state.initialize(outlvl=idaeslog.WARNING)
        get_solver().solve(m)

        zi = value(state_r32.compress_fact_phase["Vap"])
        cpi = value(state_r32.cp_mol_phase["Vap"])            # total (ideal + departure)
        zv = vanilla_z_vap(T_test, P_test, p["Tc"], p["Pc"], p["omega"])
        cpv = vanilla_cp_ideal(T_test, p["A"], p["B"], p["C"], p["D"], p["E"])
        z_idaes[name] = zi
        z_van[name] = zv
        print(f"{name:>6}{zi:10.4f}{zv:11.4f}{zi-zv:9.4f}{cpi:11.3f}{cpv:11.3f}")

    # Per-method gate: IDAES cubic solver vs the vanilla property prediction
    TOL = 1e-3
    print()
    all_pass = True
    for name in METHODS:
        dZ = abs(z_idaes[name] - z_van[name])
        ok = dZ< TOL
        all_pass &= ok ## Boolean AND on the truth table: &=
        print(f"{name} gate: |dZ| = {dZ:.2e}  -> {'PASS' if ok else 'FAIL'}")
    # Anchor: the vanilla helper reproduces the actual standalone code (NIST only)
    z_pr = pr.z_roots(T_test, P_test)[0][-1]
    print(f"\nHelper anchor (NIST): vanilla {z_van['NIST']:.4f} vs "
          f"pr_eos_lib {z_pr:.4f}  |dZ| = {abs(z_van['NIST'] - z_pr):.2e}")
    print("\nALL GATES PASSED" if all_pass else "\nSOME GATES FAILED")