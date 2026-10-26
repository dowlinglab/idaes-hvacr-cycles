"""
sweep_helmholtz_maximized_1015.py -- R-32 COP vs ambient, MAXIMIZED.

Uses vapor_compression_plr_fixed.py (Helmholtz, StateVars.TPX), which has
2 genuine degrees of freedom: the evaporator and condenser outlet
temperatures, i.e. superheat and subcooling. IPOPT picks both inside their
3-8 K windows.

WHY THIS MODEL AND NOT THE CUBIC
    Verified 10/15: with this spec the Helmholtz model reports DOF = 2 and
    fixes only five variables --

        evaporator.properties_in.flow_mass        the 1 kg/s basis
        evaporator.properties_out.vapor_frac      1.0
        compressor.properties_out.vapor_frac      1.0
        compressor.efficiency_isentropic          0.80
        condenser.properties_out.vapor_frac       0.0

    The cubic (generic PR, FTPx) fixes those five PLUS two composition
    anchors, because a generic package carries mole fractions that a pure
    Helmholtz component does not have. Those two are not spare -- unfixing
    them lets composition drift (measured mole_frac = 1.00082, above 1) and
    inflates COP above the closed form. So the cubic has 0 degrees of freedom
    and cannot maximize. See the 10/15 breadcrumb.

    Consequence: this file gives maximized COP but ONLY on Helmholtz
    properties. It cannot run NIST / GCGP / SPGP -- there is nowhere to put
    their Tc, Pc and Shomate coefficients. The PR method comparison is a
    separate, square sweep.

SPEC (per Shilpa, 10/15)
    T_evap = -29 C, condenser approach 9 K, eta = 0.80, mass flow 1 kg/s,
    superheat 3-8 K, subcooling 3-8 K, pressure ratio <= 20.

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/15/2026
"""

import sys
sys.path.insert(0, "/mnt/user-data/outputs")
import logging
import warnings
import io
import contextlib
import csv

logging.disable(logging.WARNING)
warnings.filterwarnings("ignore")

import CoolProp.CoolProp as _CP

_orig_propssi = _CP.PropsSI


def _propssi(*args):
    """The model hands CoolProp the IDAES lowercase component name."""
    args = tuple("R32" if isinstance(a, str) and a.lower() == "r32" else a
                 for a in args)
    return _orig_propssi(*args)


_CP.PropsSI = _propssi
import vapor_compression_plr_fixed as _M
_M.CP.PropsSI = _propssi

from pyomo.environ import value
from idaes.core.util.model_statistics import degrees_of_freedom
from vapor_compression_plr_fixed import SimpleVaporCompressionCyclePLRFixed as C, Mode

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.80
SH_LO, SH_HI, SC_LO, SC_HI = 3.0, 8.0, 3.0, 8.0
AMBIENTS = [10.0, 15.0, 20.0, 25.0]   # 30 C does not converge -- see 10/15 breadcrumb


def closed(Tc_C, sh, sc):
    """Closed-form Helmholtz COP at the SOLVED superheat / subcool."""
    Te, Tc = T_EVAP + 273.15, Tc_C + 273.15
    pe = _orig_propssi("P", "T", Te, "Q", 1, "R32")
    pc = _orig_propssi("P", "T", Tc, "Q", 0, "R32")
    h1 = _orig_propssi("H", "T", Te + sh, "P", pe, "R32")
    s1 = _orig_propssi("S", "T", Te + sh, "P", pe, "R32")
    h2 = h1 + (_orig_propssi("H", "P", pc, "S", s1, "R32") - h1) / ETA
    h3 = _orig_propssi("H", "T", Tc - sc, "P", pc, "R32")
    return (h1 - h3) / (h2 - h1)


rows = []
hdr = (f"{'Tamb':>5}{'Tcond':>6}{'DOF':>4}{'COP':>9}{'closed':>9}{'diff%':>8}"
       f"{'Carnot':>8}{'/Carn':>7}{'conv':>7}{'SH':>7}{'SC':>7}")
print(f"R-32, MAXIMIZED, Helmholtz.  T_evap={T_EVAP} C, approach={APPROACH} K, "
      f"eta={ETA}, SH {SH_LO:.0f}-{SH_HI:.0f} K, SC {SC_LO:.0f}-{SC_HI:.0f} K\n",
      flush=True)
print(hdr, flush=True)
print("-" * len(hdr), flush=True)

for Ta in AMBIENTS:
    Tc = Ta + APPROACH
    rec = dict(ambient_C=Ta, T_cond_C=Tc)
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc = C("r32", compressor_efficiency=ETA, mode=Mode.IMPROVED_TPX)
            vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                          high_side_temperature=Tc)
            vc.initialize(verbose=False)
            vc.set_specifications(low_side_pressure=(60.0, 200.0),
                                  high_side_pressure=(500.0, 4000.0),
                                  evaporator_temperature=(-55.0, 0.0),
                                  condenser_temperature=(-40.0, 70.0),
                                  evap_sat_temperature=T_EVAP,
                                  condenser_approach=APPROACH,
                                  superheating=SH_LO, subcooling=SC_LO,
                                  superheating_max=SH_HI, subcooling_max=SC_HI,
                                  max_pressure_ratio=20,
                                  ambient_temperature=Ta)
            dof = degrees_of_freedom(vc.model)
            cop, conv = vc.optimize_COP(verbose=False, initialize=True,
                                        optimize=True)
            m = vc.model
            eo = m.fs.evaporator.control_volume.properties_out[0]
            ko = m.fs.condenser.control_volume.properties_out[0]
            sh = value(eo.temperature) - _orig_propssi("T", "P", value(eo.pressure),
                                                       "Q", 1, "R32")
            sc = _orig_propssi("T", "P", value(ko.pressure), "Q", 0, "R32") \
                - value(ko.temperature)
            Q = value(m.fs.evaporator.heat_duty[0])
            W = value(m.fs.compressor.work_mechanical[0])
            K = value(m.fs.condenser.heat_duty[0])
        carn = (T_EVAP + 273.15) / (Tc - T_EVAP)
        cf = closed(Tc, sh, sc)
        rec.update(dof=dof, cop=cop, converged=conv, carnot=carn,
                   cop_over_carnot=cop / carn, closed_form_cop=cf,
                   diff_pct=100 * (cop - cf) / cf,
                   evap_superheat_actual_K=sh, cond_subcool_actual_K=sc,
                   Q_evap_W=Q, W_comp_W=W, Q_cond_W=K, loop_W=Q + W + K)
        print(f"{Ta:>5.0f}{Tc:>6.0f}{dof:>4}{cop:>9.4f}{cf:>9.4f}"
              f"{rec['diff_pct']:>8.2f}{carn:>8.4f}{cop / carn:>7.3f}"
              f"{str(conv):>7}{sh:>7.2f}{sc:>7.2f}", flush=True)
    except Exception as exc:
        rec.update(cop=None, converged=f"EXC {type(exc).__name__}")
        print(f"{Ta:>5.0f}{Tc:>6.0f}   FAILED {type(exc).__name__}: "
              f"{str(exc)[:45]}", flush=True)
    rows.append(rec)

keys = sorted({k for r in rows for k in r})
with open("/mnt/user-data/outputs/sweep_helmholtz_maximized_1015.csv",
          "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(rows)
print("\nwrote sweep_helmholtz_maximized_1015.csv", flush=True)
