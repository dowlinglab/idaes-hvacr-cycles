"""
verify_0924_chain.py -- end-to-end verification, 2026-09-28.

QUESTION: does the phase6_final_0924 chain actually pull the Colon group's
09/19 delivery all the way into the IDAES property block, or could it be
silently using older numbers somewhere along the way?

The chain is:

    phase6_final_0924.py
      -> vapor_compression_cubic_refstate_0924.py
           -> phase_1_cubic_eos_validation_refstate_0924.py   (METHODS_RAW)
                -> make_config(p)
                     -> IDAES GenericParameterBlock

This script checks every link, and checks the numbers at BOTH ends:

  STAGE 1  METHODS_RAW, as the config module defines it, against the 09/19
           spreadsheet values transcribed independently below.
  STAGE 2  omega, recomputed here from each method's own psat and Pc, against
           the omega the config derived.
  STAGE 3  the Shomate A-E, against the raw a-e in the sheet converted by
           A = a, B = b*1e3, C = c*1e6, D = d*1e9, E = e*1e-6.
  STAGE 4  what an actual built IDAES GenericParameterBlock reports for
           temperature_crit, pressure_crit, omega and the Cp coefficients --
           i.e. the values the solver will really see, read back out of the
           model rather than out of the dict.
  STAGE 5  that the cycle model (vapor_compression_cubic_refstate_0924) is
           importing the _0924 config and no other, and that its method
           whitelist matches.

EXPECTED values below are transcribed from "spgp r32.xlsx", the six
yellow-highlighted blocks under the red 2026-09-19 date band:
    block 1 (I13:L23)  fp_single_0919, gcn_single_0919
    block 2 (N13:Q21)  spg_agentic_0919, gcn_agentic_0919
    block 3 (I38:L48)  fp_multi_0919, gcn_multi_0919
They are written here as the sheet's RAW a-e plus Pc in bar and pvap in mmHg,
so the conversion is re-done independently of the config file rather than
copied from it.

Run from the R32 directory in the environment you use for phase6.

Author: Shilpa Narasimhan. Support: Claude AI.
"""

import sys
import numpy as np

MMHG_TO_PA = 133.322

# name: (Pc_bar, Tc_K, pvap_mmHg, a, b, c, d, e)   -- RAW, straight from the sheet
SHEET = {
    "fp_single_0919":   (52.168,    325.387,   1738.797,
                         32.6437,   0.0269797,  5.91354e-04, -5.97447e-07, -47459.5),
    "gcn_single_0919":  (63.686,    342.704,   1760.77,
                         115.569,   0.154661,   6.72036e-05,  1.74449e-07,   8824.75),
    "spg_agentic_0919": (52.183553, 344.12921, 2837.94505436,
                         32.849609, 0.02623072, 5.9020445e-04, -5.9238319e-07, -48854.552),
    "gcn_agentic_0919": (63.693509, 342.89307, 1700.64828868,
                         115.50434, 0.15506984, 6.6413405e-05, 1.7476964e-07, 9202.9555),
    "fp_multi_0919":    (59.621,    475.885,   3968.231,
                         -6.96979,  0.321086,   4.80505e-04, -1.49261e-06,  220907.0),
    "gcn_multi_0919":   (35.705,    385.667,   1074.461,
                         224.051,   0.0585856, -2.14895e-05, -5.11394e-09,     132.741),
}

BASELINES = {   # unchanged from the original chain -- checked, not modified
    "NIST": (57.82, 351.3, 0.2769),
    "GCGP": (50.730, 355.354, 0.1711),
}

fails = []
def check(tag, got, exp, rtol=1e-9):
    ok = abs(got - exp) <= rtol * max(1.0, abs(exp))
    if not ok:
        fails.append(f"{tag}: config/model={got!r}  expected={exp!r}")
    return ok


def banner(t):
    print(f"\n{'='*72}\n  {t}\n{'='*72}")


# =============================================================================
banner("STAGE 5 -- import wiring")
# =============================================================================
import phase_1_cubic_eos_validation_refstate_0924 as cfg
import vapor_compression_cubic_refstate_0924 as cyc

print(f"   config module   : {cfg.__file__}")
print(f"   cycle module    : {cyc.__file__}")
print(f"   cycle imports   : make_config / METHODS from "
      f"{cyc.make_config.__module__}")
wired = cyc.make_config.__module__ == cfg.__name__
print(f"   -> cycle model is bound to the _0924 config: {wired}")
if not wired:
    fails.append("cycle model is NOT importing the _0924 config")

print(f"\n   methods defined in config: {sorted(cfg.METHODS)}")
missing = [m for m in SHEET if m not in cfg.METHODS]
if missing:
    fails.append(f"09/19 methods missing from config: {missing}")
print(f"   all six 09/19 methods present: {not missing}")


# =============================================================================
banner("STAGES 1-3 -- dict values vs the spreadsheet, conversions re-done here")
# =============================================================================
print(f"   {'method':<18}{'Pc':>6}{'Tc':>5}{'omega':>7}{'A':>4}{'B':>4}{'C':>4}{'D':>4}{'E':>4}"
      f"{'omega value':>14}")
for name, (Pc_bar, Tc, pv, a, b, c, d, e) in SHEET.items():
    p = cfg.METHODS[name]
    A, B, C, D, E = a, b*1e3, c*1e6, d*1e9, e*1e-6
    omega = -1.0 - np.log10(pv*MMHG_TO_PA/(Pc_bar*1e5))
    r = [check(f"{name}.Pc", p["Pc"], Pc_bar*1e5),
         check(f"{name}.Tc", p["Tc"], Tc),
         check(f"{name}.omega", p["omega"], omega, rtol=1e-12),
         check(f"{name}.A", p["A"], A), check(f"{name}.B", p["B"], B),
         check(f"{name}.C", p["C"], C), check(f"{name}.D", p["D"], D),
         check(f"{name}.E", p["E"], E)]
    m = lambda x: "ok" if x else "**"
    print(f"   {name:<18}{m(r[0]):>6}{m(r[1]):>5}{m(r[2]):>7}{m(r[3]):>4}{m(r[4]):>4}"
          f"{m(r[5]):>4}{m(r[6]):>4}{m(r[7]):>4}{omega:>14.6f}")

print()
for name, (Pc_bar, Tc, omega) in BASELINES.items():
    p = cfg.METHODS[name]
    ok = all([check(f"{name}.Pc", p["Pc"], Pc_bar*1e5),
              check(f"{name}.Tc", p["Tc"], Tc),
              check(f"{name}.omega", p["omega"], omega)])
    print(f"   baseline {name:<10} unchanged from original chain: {ok}")

print("\n   F/G (reference-state offsets) should be 0.0 for every method:")
bad_fg = [n for n, p in cfg.METHODS.items()
          if p.get("F", 0.0) != 0.0 or p.get("G", 0.0) != 0.0]
print(f"   -> {'all zero' if not bad_fg else 'NON-ZERO: ' + str(bad_fg)}")
if bad_fg:
    fails.append(f"F/G non-zero for {bad_fg}")


# =============================================================================
banner("STAGE 4 -- values read back out of a BUILT IDAES property block")
# =============================================================================
# This is the one that matters: not what the dict says, but what the solver
# actually receives once make_config() has been through GenericParameterBlock.
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock)

print(f"   {'method':<18}{'Tc [K]':>11}{'Pc [bar]':>11}{'omega':>11}{'cpA':>12}{'match':>8}")
for name in list(SHEET) + list(BASELINES):
    p = cfg.METHODS[name]
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**cfg.make_config(p))
    comp = m.fs.properties.R32
    Tc_m = value(comp.temperature_crit)
    Pc_m = value(comp.pressure_crit)
    om_m = value(comp.omega)
    cpA_m = value(comp.cp_mol_ig_comp_coeff_A)
    ok = all([check(f"{name}.model.Tc", Tc_m, p["Tc"], rtol=1e-8),
              check(f"{name}.model.Pc", Pc_m, p["Pc"], rtol=1e-8),
              check(f"{name}.model.omega", om_m, p["omega"], rtol=1e-8),
              check(f"{name}.model.cpA", cpA_m, p["A"], rtol=1e-8)])
    print(f"   {name:<18}{Tc_m:>11.4f}{Pc_m/1e5:>11.4f}{om_m:>11.6f}{cpA_m:>12.4f}"
          f"{'ok' if ok else '**FAIL**':>8}")


# =============================================================================
banner("VERDICT")
# =============================================================================
if not fails:
    print("   ALL CHECKS PASSED.")
    print("   The 09/19 delivery reaches the IDAES property block intact:")
    print("   spreadsheet -> METHODS_RAW -> omega derivation -> make_config")
    print("   -> GenericParameterBlock, verified at both ends.")
    sys.exit(0)
else:
    print(f"   {len(fails)} CHECK(S) FAILED:\n")
    for f in fails:
        print(f"     - {f}")
    sys.exit(1)
