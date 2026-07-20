"""
This is the first accuracy gate: does the IDAES package reproduce the real
 R-32 saturation dome (Linde)? 
This code reuses make_config/METHODS from Phase 1

Trick to getting getting saturated liquid and vapor 
properties out of IDAES: At saturation a pure fluid 
sits at one pressure (Psat) where both phases coexist. 
The clean IDAES trick: fix T, then free pressure 
and force two phases by fixing the vapor fraction
 — the solver then finds P = Psat 
 and hands you both phase enthalpies in one solve.

 For a pure fluid at saturation: C = 1, P = 2 (liquid + vapor), 
 so F = 1 − 2 + 2 = 1. On the saturation curve there is exactly 
 one free variable — fix the temperature and everything else follows:
the pressure must be Psat(T), and the saturated-liquid and saturated-vapor 
states are both fully determined. This is why you can't independently 
choose both T and P and still get two-phase coexistence; you only get to pick one.

 IIR conditions are assumed (sat. liquid at 0 °C → h = 200 kJ/kg, s = 1.0 kJ/kg·K)
Date created: 07/20/2026

Author: Shilpa Narasimhan

Support: Claude AI

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
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
import pr_eos_lib as pr   # NIST-only, hand-written reference (the "vanilla" code)
import idaes.logger as idaeslog


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
# "omega" = Pitzer acentric factor, computed per method from the Linde vapor
# pressure at Tr = 0.7:  omega = -1 - log10(Psat(Tr=0.7) / Pc), using each
# method's own Tc and Pc. Values: NIST 0.2769, GCGP 0.1711, SPGP -0.2741.
# (Negative for SPGP because its Tc is far too high -> a symptom of bad
#  critical properties, not a real acentric factor.)

METHODS = {
    "NIST": {"Pc": 57.82e5,  "Tc": 351.3,   "omega": 0.2769,
             "A": -6.098682, "B": 179.2200, "C": -122.3682, "D": 32.30207, "E": 0.491361},
    "GCGP": {"Pc": 50.730e5, "Tc": 355.354, "omega": 0.1711,
             "A": 14.161,    "B": 0.124,    "C": -6.340e-05, "D": 1.190e-8, "E": 0.0},
    "SPGP": {"Pc": 50.8106e5,"Tc": 400.898, "omega": -0.2741,
             "A": 129.687,   "B": 171.303,  "C": 146.2,      "D": 61.9837,  
             "E": -0.0000361638},
}
## To provide initial guess pressure to the solver
def antoine_psat(T):
    """Saturation pressure [Pa] from the Antoine fit: log10(P_bar) = A - B/(T+C)."""
    return 10.0**(Antoine_A - Antoine_B/(T + Antoine_C)) * 1e5
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
                        "F": (0.0, pyunits.kJ/pyunits.mol),
                        "G": (0.0, pyunits.J/pyunits.mol/pyunits.K),
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
#   Extract saturated liquid/vapor properties from IDAES
##############################################################################################


def sat_point(p,T):
    """
    Saturated molar h, s for liquid and vapor at temperature T [K].
    Evaluate each phase's EoS root at (T, Psat): enth_mol_phase["Liq"] is the
    saturated liquid, enth_mol_phase["Vap"] the saturated vapor. Psat from the
    Antoine fit -- sitting on the saturation curve gives both roots with no flash,
    avoiding the trivial-solution collapse of a free-pressure solve.
    """
    Psat = antoine_psat(T)
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic = False)
    m.fs.properties = GenericParameterBlock(**make_config(p))
    m.fs.state = m.fs.properties.build_state_block([0], defined_state = True)
    state_r32 = m.fs.state[0]

    state_r32.flow_mol.fix(1.0)
    state_r32.mole_frac_comp["R32"].fix(1.0)
    state_r32.temperature.fix(T)
    state_r32.pressure.fix(Psat)
    m.fs.state.initialize(outlvl = idaeslog.WARNING)

    get_solver().solve(m)

    Psat = value(state_r32.pressure)
    hl = value(state_r32.enth_mol_phase["Liq"])
    hg = value(state_r32.enth_mol_phase["Vap"])
    sl = value(state_r32.entr_mol_phase["Liq"])
    sg = value(state_r32.entr_mol_phase["Vap"])
    return Psat, hl, hg, sl,sg

##############################################################################################
#   Single-point test: confirm the saturation recipe at 0 deg C (NIST)
##############################################################################################

if __name__ == "__main__":
    T0 = 273.15 # 0 deg C
    Psat,hl,hg,sl,sg = sat_point(METHODS["NIST"],T0)
    print("NIST saturated at 0 C:")
    print(f"  Psat = {Psat/1e5:.3f} bar        (Linde ~ 8.13)")
    print(f"  hl = {hl:.1f} J/mol   hg = {hg:.1f} J/mol   (raw, zero datum)")
    print(f"  latent heat hg-hl = {(hg-hl)/MW/1000:.1f} kJ/kg   (Linde ~ 316)")