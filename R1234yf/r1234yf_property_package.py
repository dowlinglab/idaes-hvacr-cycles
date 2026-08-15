################################################################################
# R1234yf IDAES Helmholtz Property Package
#
# Author: Shilpa Narasimhan (snarasi2@nd.edu)
# Support: Claude AI (claude@anthropic.com)
# Date Created: 2026-08-14
# 
# Description:
#   Custom IDAES property package for R1234yf (HFO-1234yf) refrigerant
#   implementing the Helmholtz EOS from Lemmon & Akasaka (2022).
#   Coefficients extracted from CoolProp validation.
#
# References:
#   - Lemmon, E.W., Akasaka, R. (2022). Fundamental equation of state for 
#     2,3,3,3-tetrafluoropropene (HFO-1234yf). Int. J. Thermophys. 43, 171
#   - CoolProp: https://coolprop.org/fluid_properties/fluids/R1234yf.html
#   - NIST Experimental Data: Richter et al., J. Phys. Chem. Ref. Data
#
# BLOCK 1: Imports & Setup
# ============================================================================
# Purpose: Load parameters from r1234yf.json and prepare for property calculation
#
################################################################################

import json
import os

import pyomo.environ as pyo
from pyomo.environ import units as pu
from pyomo.environ import Constraint, Param, Var, Expression

from idaes.core import (
    PropertyParameterBlock,
    PropertyStateBlock,
)

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

# ============================================================================
# BLOCK 1: Load Parameters
# ============================================================================

"""
Load R1234yf EOS parameters from r1234yf.json file located in same directory.
The JSON file contains:
  - basic: Critical constants (Tc, Pc, rhoc, MW, R, Tt, Pt, etc.)
  - eos: EOS coefficients for ideal gas (n0, g0) and residual (n, d, t, a, b, e, g)
  - aux: Saturation approximations
  - transport: Thermal conductivity, viscosity, surface tension (if available)
"""

json_file = os.path.join(os.path.dirname(__file__), 'r1234yf.json')
with open(json_file, 'r') as f:
    params = json.load(f)

# ============================================================================
# Extract Critical Constants (from params['basic'])
# ============================================================================
"""
Critical constants define the reduced variables:
  δ = ρ / ρc  (reduced density)
  τ = Tc / T  (reduced inverse temperature)
  
These are fundamental to the Helmholtz EOS formulation.
"""

Tc = params['basic']['Tc']        # Critical temperature [K]
Pc = params['basic']['Pc']        # Critical pressure [Pa]
rhoc = params['basic']['rhoc']    # Critical density [kg/m³]
R_gas = params['basic']['R']      # Specific gas constant [J/(kg·K)]
MW = params['basic']['MW']        # Molar weight [g/mol]

# ============================================================================
# Extract Ideal Gas Coefficients (phi_ideal_type = 1)
# ============================================================================
"""
Ideal gas Helmholtz energy:
  φ⁰(δ, τ) = n₀[1]·ln(τ) + Σ(n₀[i]) + Σ(n₀[j]·g₀[j]·τ^g₀[j])
  
where:
  n0[1] = log(τ) coefficient
  n0[2:3] = linear and constant terms
  n0[4:6] = Planck-Einstein oscillator numerators
  g0[4:6] = Planck-Einstein oscillator tau exponents
"""

n0 = {int(k): v for k, v in params['eos']['n0'].items()}
g0 = {int(k): v for k, v in params['eos']['g0'].items()}
last_term_ideal = params['eos']['last_term_ideal']

# ============================================================================
# Extract Residual Coefficients (phi_residual_type = 2)
# ============================================================================
"""
Residual Helmholtz energy (power law + Gaussian bell):
  φʳ(δ, τ) = Σ(n[i]·δ^d[i]·τ^t[i])                                    [Terms 1-10]
           + Σ(n[i]·δ^d[i]·τ^t[i]·exp(-a[i]·(δ-e[i])² - b[i]·(τ-g[i])²))  [Terms 11-17]

Power law terms (1-10):
  n, d, t = amplitude, delta exponent, tau exponent

Gaussian bell terms (11-17):
  n = amplitude
  d = delta exponent (typically 1 for bell terms)
  t = tau exponent
  a = Gaussian width in δ-direction (from CoolProp eta)
  e = Gaussian center in δ-direction (from CoolProp beta)
  b = Gaussian width in τ-direction (from CoolProp gamma)
  g = Gaussian center in τ-direction (from CoolProp epsilon)
"""

n = {int(k): v for k, v in params['eos']['n'].items()}
d = {int(k): v for k, v in params['eos']['d'].items()}
t = {int(k): v for k, v in params['eos']['t'].items()}
a = {int(k): v for k, v in params['eos']['a'].items()}
b = {int(k): v for k, v in params['eos']['b'].items()}
e = {int(k): v for k, v in params['eos']['e'].items()}
g = {int(k): v for k, v in params['eos']['g'].items()}
last_term_residual = params['eos']['last_term_residual']

# ============================================================================
# Print Summary (Validation)
# ============================================================================
"""
Print parameter summary to verify correct loading.
Useful for debugging and documentation.
"""

print("\n" + "="*70)
print("R1234yf IDAES Property Package - Block 1: Parameter Loading")
print("="*70)

print(f"\n✓ Loaded R1234yf parameters from: {json_file}")

print(f"\n--- Critical Constants ---")
print(f"  Tc (Critical Temperature):    {Tc:>10.2f} K")
print(f"  Pc (Critical Pressure):       {Pc:>10.0f} Pa = {Pc/1e6:.4f} MPa")
print(f"  ρc (Critical Density):        {rhoc:>10.2f} kg/m³")
print(f"  R  (Specific Gas Constant):   {R_gas:>10.6f} J/(kg·K)")
print(f"  MW (Molar Weight):            {MW:>10.3f} g/mol")

print(f"\n--- Ideal Gas Terms (n0, g0) ---")
print(f"  Number of terms: {len(n0)}")
print(f"  last_term_ideal: {last_term_ideal}")
print(f"  n0 keys: {sorted(n0.keys())}")
print(f"  g0 keys: {sorted(g0.keys())}")
print(f"\n  n0 values (ideal gas amplitudes):")
for key in sorted(n0.keys()):
    print(f"    n0[{key}] = {n0[key]:>15.10f}")
print(f"\n  g0 values (Planck-Einstein exponents):")
for key in sorted(g0.keys()):
    print(f"    g0[{key}] = {g0[key]:>15.10f}")

print(f"\n--- Residual Terms (n, d, t, a, b, e, g) ---")
print(f"  Total number of terms: {len(n)}")
print(f"  last_term_residual: {last_term_residual}")
print(f"    → Power law terms:    1-{last_term_residual[0]}")
print(f"    → Gaussian bell terms: {last_term_residual[0]+1}-{last_term_residual[1]}")

print(f"\n  Power law terms (1-10) - d, t exponents:")
for i in range(1, last_term_residual[0]+1):
    print(f"    Term {i:>2d}: d={d[i]:>4.1f}, t={t[i]:>6.3f}, n={n[i]:>12.8f}")

print(f"\n  Gaussian bell terms (11-17) - Gaussian parameters:")
for i in range(last_term_residual[0]+1, last_term_residual[1]+1):
    print(f"    Term {i:>2d}: a={a[i]:>7.3f}, e={e[i]:>7.3f}, b={b[i]:>7.3f}, g={g[i]:>7.3f}")

print(f"\n{'='*70}")
print(f"✓ Block 1 Complete: All parameters loaded successfully")
print(f"{'='*70}\n")

