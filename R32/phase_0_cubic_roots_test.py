"""
This code serves as a test to see that IDAES is able to solve the cubic equations
of state for all methods within the Colon group collaboration with 0 degrees of freedom.

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
## Tell IDAES which variables describe the thermodynamic state of each stream/point
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



