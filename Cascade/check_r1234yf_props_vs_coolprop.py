"""
check_r1234yf_props_vs_coolprop.py

Author:  Shilpa Narasimhan
Support: Claude (Anthropic)

Purpose
-------
Compare our custom R-1234yf property package against CoolProp, point by
point, to find numerical irregularities in how our package is built.

Background: IPOPT struggles with the R-1234yf/CO2 cascade at 45 C. The
R-134a/CO2 cascade has the same equations and converges fine. The main
difference is the property package:
    R-134a   -> IDAES's built-in general_helmholtz file
    R-1234yf -> our own file, R1234yf/r1234yf_general_helmholtz.json
Whatever we find here should also guide fixes for R-515B.

Status
------
Step 1  Register the custom R-1234yf package          done
Step 2  Get a reference value from CoolProp           done
Step 3  Get the same value from the custom package    done
Step 4  Values at the cascade hot-loop states, 45 C   done
Step 5  Slopes (cp) at the same states                done
Step 6  T_sat(P) values and slopes, hot-loop P       done
Step 7  Curvature (2nd derivatives) -- IPOPT's Hessian  done

Result (2026-10-07): at the 45 C cascade hot-loop states, values,
cp slopes and T_sat(P) slopes all match CoolProp to ~1e-4 % and stay
smooth down to tiny steps. No irregularity found at these states.

Requirements
------------
- Run inside myidaesenv.
- CoolProp 8.0 (upgraded 2026-10-07). It uses the same EOS as our
  package, Lemmon (2022), so any difference is from how our package is
  built. CoolProp 7.x uses Richter (2011), which differs by 0.1-0.3%
  and would show up as a false mismatch.

Run
---
    cd ~/idaes-hvacr-cycles/Cascade
    python check_r1234yf_props_vs_coolprop.py
"""

import os                    # build file paths
import sys                   # control where Python looks for imports
from idaes.models.properties.general_helmholtz import (
    registered_components,
    HelmholtzParameterBlock,     # the IDAES property package
    AmountBasis,                 # per-kg or per-mol output
)
from pyomo.environ import ConcreteModel, value, units as pyunits
                                 # model container, expression -> number, units
import CoolProp.CoolProp as CP   # CoolProp property calculator

from idaes.models.properties.general_helmholtz.helmholtz_functions import HelmholtzThermoExpressions


# --- Step 1: make the custom R-1234yf package available ---------------
# The package file lives in ../R1234yf, which Python does not search by
# default. Add that folder to the search path, then import the module.
# Importing it registers our fluid with IDAES under the name "r1234yf".
# The "Check approx sat delta curves" table printed on run comes from
# that registration, not from this script.

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
R1234YF_DIR = os.path.join(THIS_DIR, "..", "R1234yf")

sys.path.insert(0, R1234YF_DIR)

import vapor_compression_plr_r1234yf

print("Custom R-1234yf registered:", "r1234yf" in registered_components())


# --- Step 2: CoolProp reference value ---------------------------------
# Saturated vapor (Q = 1) at -29 C, the evaporator outlet.
# PropsSI returns SI units: P in Pa, H in J/kg.
# Fluid name: CoolProp uses "R1234yf" (case-sensitive, no hyphen);
# IDAES uses "r1234yf".

T = -29 + 273.15

P_sat = CP.PropsSI("P", "T", T, "Q", 1, "R1234yf")
h_coolprop = CP.PropsSI("H","T", T, "Q", 1, "R1234yf")

print("Saturation pressure (kPa):", P_sat / 1000)
print("CoolProp enthalpy (J/kg):", h_coolprop)


# --- Check: which CoolProp version and EOS are in use ------------------
# Should print 8.0.0 and Lemmon-IJT-2022 (see Requirements above).

import CoolProp
print("CoolProp version:", CoolProp.__version__)
print("R1234yf EOS source:", CP.get_fluid_param_string("R1234yf", "BibTeX-EOS"))


# --- Step 3: same value from our custom package -----------------------
# Build a tiny model that holds only the R-1234yf property package.
# MASS basis -> J/kg, same units as CoolProp.
# htpx takes exactly two of (T, p, x); saturated vapor -> T and x = 1.
# IDAES needs units attached (T * pyunits.K).

m = ConcreteModel()
m.r1234yf = HelmholtzParameterBlock(
    pure_component="r1234yf",
    amount_basis=AmountBasis.MASS,
)

h_custom = value(m.r1234yf.htpx(T=T * pyunits.K, x=1))

print("Custom   enthalpy (J/kg):", h_custom)
print("CoolProp enthalpy (J/kg):", h_coolprop)
print("Difference (%):", 100 * (h_custom - h_coolprop) / h_coolprop)


# --- Step 4: values at the cascade hot-loop states, 45 C ---------------
# R-1234yf side of the cascade, from cop_vs_ambient_plr_benchmark_r1234yf.csv
# (45 C row). (T in C, P in kPa); single-phase, so T and P set the state.
states_45C = {
    "evaporator outlet": (11.037, 373.369),
    "compressor outlet": (61.163, 1431.063),
    "condenser outlet": (48.086, 1431.063),
}

for name, (T_C, P_kPa) in states_45C.items():
    T_K = T_C + 273.15
    P_Pa = P_kPa* 1000
    h_custom = value(m.r1234yf.htpx(T = T_K * pyunits.K, p = P_Pa * pyunits.Pa))
    h_cp = CP.PropsSI("H", "T", T_K, "P", P_Pa, "R1234yf")
    diff = 100*(h_custom - h_cp)/h_cp
    print(f"{name}: diff = {diff:.6f}%")

# --- Step 5: slopes, cp = dh/dT at constant P ---------------------------
# Central differences with shrinking steps. Smooth = values settle and
# stay; irregular = they jump at small steps.
def h_custom(T,P):
    return value(m.r1234yf.htpx(T = T*pyunits.K, p = P* pyunits.Pa))

for name, (T_C, P_KPa) in states_45C.items():
    T_K = T_C + 273.15
    P_Pa = P_KPa*1000
    cp_cp = CP.PropsSI("d(Hmass)/d(T)|P", "T", T_K, "P", P_Pa, "R1234yf")
    print(f"{name}: CoolProp cp = {cp_cp:.4f}")
    for dT in [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
        cp_custom = (h_custom(T_K + dT, P_Pa) - h_custom(T_K -dT, P_Pa))/ (2*dT)
        print(f" step {dT}: diff = {100*(cp_custom - cp_cp)/cp_cp:.5f}%")




te = HelmholtzThermoExpressions(m, m.r1234yf)

def Tsat_custom(P):
    return value(te.T_sat(p = P*pyunits.Pa))

pressures_kPa = {"hot evaporator": 373.369, "hot condenser": 1431.063}

for name, P_kPa in pressures_kPa.items():
    P = P_kPa * 1000
    T_cp = CP.PropsSI("T", "P", P, "Q", 1, "R1234yf")
    print(f"{name}: T_sat diff = {Tsat_custom(P) - T_cp:.6f} K")
    slope_cp = CP.PropsSI("d(T)/d(P)|sigma", "P", P, "Q", 1, "R1234yf")
    for dP in [1000.0, 100.0, 10.0, 1.0, 0.1, 0.01]:
        s = (Tsat_custom(P+dP) - Tsat_custom(P-dP))/(2*dP)
        print(f" step {dP} Pa: slope diff = {100*(s-slope_cp)/slope_cp:.5f}%")


# --- Step 7: curvature (second derivatives) ----------------------------
# IPOPT also uses second derivatives (its Hessian). Compare curvature
# from our package against CoolProp's exact value. As a fair reference,
# do the same finite difference on CoolProp too: if ours turns noisy at
# a step where CoolProp's stays clean, that is an irregularity.
# AbstractState is CoolProp's lower-level interface; it gives exact
# second derivatives that PropsSI does not.

AS = CP.AbstractState("HEOS", "R1234yf")

def h_coolprop(T, P):
    return CP.PropsSI("H", "T", T, "P", P, "R1234yf")

def Tsat_coolprop(P):
    return CP.PropsSI("T", "P", P, "Q", 1, "R1234yf")

print("d2h/dT2 at constant P")
for name, (T_C, P_kPa) in states_45C.items():
    T_K = T_C + 273.15
    P_Pa = P_kPa * 1000
    AS.update(CP.PT_INPUTS, P_Pa, T_K)
    exact = AS.second_partial_deriv(CP.iHmass, CP.iT, CP.iP, CP.iT, CP.iP)
    print(f"{name}: CoolProp exact = {exact:.6f}")
    for d in [1.0, 0.1, 0.01, 0.001]:
        ours = (h_custom(T_K + d, P_Pa) - 2 * h_custom(T_K, P_Pa) + h_custom(T_K - d, P_Pa)) / d**2
        cpfd = (h_coolprop(T_K + d, P_Pa) - 2 * h_coolprop(T_K, P_Pa) + h_coolprop(T_K - d, P_Pa)) / d**2
        print(f"   step {d}: ours diff = {100 * (ours - exact) / exact:.4f} %, "
              f"CoolProp-FD diff = {100 * (cpfd - exact) / exact:.4f} %")

print("d2Tsat/dP2")
for name, P_kPa in pressures_kPa.items():
    P = P_kPa * 1000
    AS.update(CP.PQ_INPUTS, P, 1)
    exact = AS.second_saturation_deriv(CP.iT, CP.iP, CP.iP)
    print(f"{name}: CoolProp exact = {exact:.6e}")
    for d in [1000.0, 100.0, 10.0, 1.0]:
        ours = (Tsat_custom(P + d) - 2 * Tsat_custom(P) + Tsat_custom(P - d)) / d**2
        cpfd = (Tsat_coolprop(P + d) - 2 * Tsat_coolprop(P) + Tsat_coolprop(P - d)) / d**2
        print(f"   step {d} Pa: ours diff = {100 * (ours - exact) / exact:.4f} %, "
              f"CoolProp-FD diff = {100 * (cpfd - exact) / exact:.4f} %")
