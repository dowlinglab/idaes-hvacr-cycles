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
from idaes.models.properties.modular_properties.pure import NIST
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity

