"""
Phase 3a validation: run OUR Helmholtz cycle at published R-32 conditions and
compare COP to the literature, with a Carnot sanity check.

Source (this file):
  Suhengki, Manik H.M., Hestirianoto T., Susilohadi, Yulianto M. (2026).
  "Performance Parameter Study of R32 for Low Temperature Cold Storage."
  International Journal on Advanced Science, Engineering and Information
  Technology (IJASEIT), Vol. 16 No. 2, ISSN 2088-5334. Steady-state Excel-VBA
  simulation ("The Poci") + experimental cold-storage rig; ambient 25-40 C,
  condensing 45-60 C, evaporator to -28 C, chamber -18 to -20 C. R-32 crit.
  props Tc = 78.4 C, Pc = 58.3 bar (match ours). Reported R-32 COP ~ 2.8
  (abstract); MAPE 2.11%.
R-32, cold storage -- our exact application. Validation point (their Table V):
    evaporating ~ -29 C (2.81 bar), condensing ~ 34 C (21 bar),
    superheat ~ 42 K, subcool ~ 0 K, discharge 84.8 C.
    Paper COP: Table V = 7.5  (EXCEEDS Carnot 3.9 -> erroneous);
               abstract  = 2.8 (realistic system COP).
Carnot at Te = -29 / Tc = 34 is 244.15/(307.15-244.15) = 3.88, so any physical
cycle MUST give COP < 3.88. If our model returns ~3 (below Carnot) we have both
a physical result and a demonstration that their Table V COP = 7.5 is impossible.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 07/20/2026
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
from pyomo.environ import Var, value

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vapor_compression_plr import SimpleVaporCompressionCycle, Mode


def relax_enth_bounds(vc, hmax=700e3):
    for v in vc.model.component_data_objects(Var, active=True, descend_into=True):
        if v.parent_component().local_name == "enth_mass" and v.ub is not None and v.ub < hmax:
            v.setub(hmax)


def carnot_cop(te_c, tc_c):
    return (te_c + 273.15) / ((tc_c + 273.15) - (te_c + 273.15))


# name, Te_sat, Tc_sat, ambient, approach, superheat, subcool, eta, paper_COP
CASES = [
    ("Suhengki 2026 cold storage", -29, 34, 32, 2, 42, 0, 0.70,
     "7.5 (Table V, > Carnot!) / 2.8 (abstract)"),
]

print(f"\n{'case':<30}{'COP':>7}{'Carnot':>8}{'frac':>6}   paper COP")
print("-" * 78)
for name, Te, Tc, Tamb, appr, SH, SC, eta, paper in CASES:
    vc = SimpleVaporCompressionCycle("R32", compressor_efficiency=eta,
                                     mode=Mode.IMPROVED_TPX)
    relax_enth_bounds(vc)
    vc.specify_initial_conditions(low_side_temperature=Te, high_side_temperature=Tc)
    vc.initialize(verbose=False)
    vc.set_specifications(
        ambient_temperature=Tamb,
        condenser_approach=appr,        # Tc_sat = Tamb + appr = 34
        evap_sat_temperature=Te,
        superheating=SH,
        subcooling=SC,
        max_pressure_ratio=12,
    )
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    carnot = carnot_cop(Te, Tc)
    if converged:
        print(f"{name:<30}{cop:>7.3f}{carnot:>8.2f}{cop/carnot:>6.2f}   {paper}")
    else:
        print(f"{name:<30}{'FAIL':>7}{carnot:>8.2f}{'--':>6}   {paper}")

print("\nNote: physical COP must be < Carnot. Their Table V COP 7.5 > 3.88 = Carnot,")
print("so it violates the second law -- our model should land below Carnot (~abstract 2.8).")
