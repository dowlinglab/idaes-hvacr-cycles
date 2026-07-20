"""
This code serves as a test to see that IDAES is able to solve the cubic equations
of state with 0 degrees of freedom.

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
#  (ln f_liquid = ln f_vapor). This is the same fugacity condition we talked about — 
# it's when the two phases balance.

