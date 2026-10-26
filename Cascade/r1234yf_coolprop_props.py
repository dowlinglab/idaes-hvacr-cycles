from pyomo.environ import units as pyunits
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.general_helmholtz import registered_components
from CoolProp.CoolProp as CP
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    AmountBasis,
)
from pyomo.environ import ConcreteModel, value

m = ConcreteModel()
m.r1234yf = HelmholtzParameterBlock(
    pure_component = "r1234yf",
    amount_basis = AmountBasis.MASS,
)

h_custom = value(m.r1234yf.htpx(T= T * pyunits.K, x = 1))

print("Custom   enthalpy (J/kg):", h_custom)
print("CoolProp enthalpy (J/kg):", h_coolprop)
print("Difference (%):", 100 * (h_custom - h_coolprop) / h_coolprop)
