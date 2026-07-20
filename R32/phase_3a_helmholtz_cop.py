"""
Phase 3a: COP vs ambient for R-32 using the working PLR-branch Helmholtz cycle
(vapor_compression_plr.py, extracted from origin/PLR_vanilla_prop, unchanged).
Baseline before Phase 3b swaps in the cubic PR package.

Uses the improved spec interface:
  - evap_sat_temperature = -29 C   (T_evap in the -30..-28 band, midpoint)
  - ambient_temperature = T_amb, condenser_approach = 9  -> T_cond_sat = T_amb+9
  - superheating / subcooling = 3 / 3 K   (subcool -> outlet below T_sat)
  - max pressure ratio = 10 (single-stage limit, upper end of the 8-10 band)
  - optimize=False, Mode.IMPROVED_TPX (as the PLR harness uses)
  - Reported range T_amb 10-25 C only; above ~25 C the fixed-cycle solve returns
    infeasible / non-physical points (see AMBIENTS note and PHASE3_NOTES).

Validation reference (to reproduce as a cycle check, separate from the -29 C
cold-storage sweep above):
  Taira, S., Minamida, T., Haikawa, T., Ohta, F. (2016). "Performance Evaluation
  of Heat Pump System using R32 and HFO-mixed Refrigerant in High Ambient
  Temperature." 16th Int. Refrigeration and Air Conditioning Conf. at Purdue,
  Paper 2408/1736.  R-32 (pure), cooling, computed with REFPROP 9.1:
    Tc = 46 C, Te = 12.5 C, suction 15 C (~2.5 K superheat),
    condenser outlet 38 C (8 K subcool), eta_comp = 70%  ->  COP = 5.01
    (refrigerating effect 249.3 kJ/kg, compressor work 49.8 kJ/kg, Td = 83.4 C).
  Same Helmholtz/REFPROP basis as IDAES, so our cycle at these conditions should
  reproduce COP ~ 5.0.  (Tc = 46 C corresponds to ~35 C outdoor + ~11 K condenser
  approach; the paper's experiments ran to 52 C outdoor.)

Ambient (outdoor) temperature ranges in the literature (AC duty, warm indoor
evaporator -- context for "how hot single-stage R-32 runs", NOT our -29 C case):
  - Daikin/Taira 2016 experiments (ISO 5151): outdoor DB 27 / 35 / 52 C
    (T2 / T1 / T3(H)); indoor DB 21 / 27 / 32 C.
  - REHVA/Badescu 2025 regression (from Shen et al. 2016): outdoor 27.8-55 C
    at indoor 26.7-29 C. R-32 fit COP = 1/(0.2526 + 0.0003*dT^2),
    dT = T_outdoor - T_indoor (COP ~2.4 at dT=23 up to ~4.0 at dT=0).
  - Payne & Domanski 2002 pushed R410A tests to 68.3 C ambient (single stage).
This Phase-3a cold-storage sweep instead uses T_amb 10-25 C with a -29 C
evaporator -- a different (deep-refrigeration) regime, so these AC ambient
ranges are context, not a direct comparison.

Most-relevant reference (SAME application, R-32 cold storage):
  Suhengki, Manik H.M., Hestirianoto T., Susilohadi, Yulianto M. (2026).
  "Performance Parameter Study of R32 for Low Temperature Cold Storage,"
  IJASEIT 16(2), ISSN 2088-5334. Ambient 25-40 C, condensing 45-60 C,
  evaporator to -28 C, chamber -18 to -20 C; R-32 crit props Tc = 78.4 C,
  Pc = 58.3 bar (match ours). Realistic R-32 cold-storage COP ~ 2.8 (abstract),
  consistent with our physical 2.5-3.3. NOTE their Table V COP = 7.5 exceeds
  Carnot (3.9) and is erroneous -- see phase_3a_validation_literature.py.

Author: Shilpa Narasimhan   Support: Claude AI
Date created: 07/20/2026
"""
import os, sys
import matplotlib
matplotlib.use("Agg")   # no blocking plot windows during the sweep (env only, class unchanged)

from pyomo.environ import Var, value

# working PLR-branch cycle, extracted into this folder
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vapor_compression_plr import SimpleVaporCompressionCycle, Mode

FLUID = "R32"
AMBIENTS = [10, 15, 20, 25]   # deg C, 5 deg intervals.
# Only 10-25 C is reported: above ~25 C the fixed-cycle feasibility solve returns
# infeasible / non-physical operating points (the superheat/subcool inequalities
# leave the state under-determined and the solver drifts off the physical branch;
# optimize=True instead runs subcooling below ambient). 10-25 C converges cleanly
# and monotonically -- sufficient for the Phase 4 property-package comparison,
# which is done at a single fixed condition.
H_MAX = 700e3   # J/kg -- relaxed enth_mass ceiling (default ~500 kJ/kg is below
                # R-32's saturated-vapor enthalpy at a -28 C evaporator)


def relax_enth_bounds(vc, hmax=H_MAX):
    """Raise the upper bound on every enth_mass state variable in the model."""
    for v in vc.model.component_data_objects(Var, active=True, descend_into=True):
        if v.parent_component().local_name == "enth_mass":
            if v.ub is not None and v.ub < hmax:
                v.setub(hmax)


print(f"\n{'T_amb':>6}{'T_cond_sat':>12}{'COP':>8}{'converged':>11}")
T_EVAP_SAT = -29.0   # deg C, evaporating temperature (fixed across the sweep)


def carnot_cop(t_evap_c, t_cond_c):
    """Reversible (Carnot) cooling COP = T_L / (T_H - T_L), temperatures in K."""
    T_L = t_evap_c + 273.15
    T_H = t_cond_c + 273.15
    return T_L / (T_H - T_L)


print(f"\n{'T_amb':>6}{'Tcond':>7}{'COP':>7}{'Carnot':>8}{'frac':>6}"
      f"{'SH':>6}{'SC':>6}{'Qevap':>8}{'Wcomp':>8}{'ratio':>7}{'fallbk':>7}")
print("-" * 78)

rows = []
for Tamb in AMBIENTS:
    Tcond_sat = Tamb + 9
    vc = SimpleVaporCompressionCycle(FLUID, compressor_efficiency=0.9999,
                                     mode=Mode.IMPROVED_TPX)   # ~ideal (Shridhar)
    relax_enth_bounds(vc)
    vc.specify_initial_conditions(low_side_temperature=-29, high_side_temperature=Tcond_sat)
    vc.initialize(verbose=False)

    def apply_specs(disable_arc_p):
        vc.set_specifications(
            ambient_temperature=Tamb,
            condenser_approach=9,          # T_cond_sat = T_amb + 9
            evap_sat_temperature=-29,      # T_evap in the -30..-28 band
            superheating=0,                # ideal cycle: saturated vapor out of evap
            subcooling=0,                  # ideal cycle: saturated liquid out of cond
            max_pressure_ratio=10,
            debug_disable_arc_pressure_eq=disable_arc_p,
        )

    apply_specs(False)
    cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
    fallback = False
    if not converged:                      # PLR fallback: relax the loop pressure eq
        fallback = True
        apply_specs(True)
        cop, converged = vc.optimize_COP(verbose=False, initialize=True, optimize=False)

    m = vc.model
    Qe = value(m.fs.evaporator.heat_duty[0]) / 1e3          # kW (per kg/s)
    Wc = value(m.fs.compressor.work_mechanical[0]) / 1e3    # kW
    Plow = value(m.fs.compressor.inlet.pressure[0]) / 1e5   # bar
    Phigh = value(m.fs.compressor.outlet.pressure[0]) / 1e5 # bar
    ratio = Phigh / Plow
    # actual superheat / subcool (K) at the solved point
    e_out = m.fs.evaporator.control_volume.properties_out[0]
    c_out = m.fs.condenser.control_volume.properties_out[0]
    SH = value(e_out.temperature) - value(e_out.temperature_sat)
    SC = value(c_out.temperature_sat) - value(c_out.temperature)
    carnot = carnot_cop(T_EVAP_SAT, Tcond_sat)
    frac = (cop / carnot) if converged else None

    rows.append({
        "ambient_C": Tamb, "T_cond_sat_C": Tcond_sat,
        "cop": cop if converged else "",
        "carnot_cop": round(carnot, 4),
        "cop_over_carnot": round(frac, 4) if frac is not None else "",
        "superheat_K": round(SH, 2), "subcool_K": round(SC, 2),
        "Q_evap_kW": round(Qe, 3), "W_comp_kW": round(Wc, 3),
        "P_low_bar": round(Plow, 3), "P_high_bar": round(Phigh, 3),
        "pressure_ratio": round(ratio, 3),
        "fallback_arc_p_disabled": fallback, "converged": converged,
    })

    cop_str = f"{cop:.3f}" if converged else "--"
    frac_str = f"{frac:.2f}" if converged else "--"
    print(f"{Tamb:>6}{Tcond_sat:>7}{cop_str:>7}{carnot:>8.2f}{frac_str:>6}"
          f"{SH:>6.1f}{SC:>6.1f}{Qe:>8.2f}{Wc:>8.2f}{ratio:>7.2f}{str(fallback):>7}")

# --- write CSV next to this script ---
import csv
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "cop_vs_ambient_r32.csv")
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
print(f"\nwrote {out_path}")

# --- plot COP and Carnot COP vs ambient ---
import matplotlib.pyplot as plt

conv = [r for r in rows if r["converged"]]
amb = [r["ambient_C"] for r in conv]
cops = [r["cop"] for r in conv]
carn = [r["carnot_cop"] for r in conv]

plt.figure(figsize=(7, 5))
plt.plot(amb, carn, "k--o", label="Carnot COP (reversible limit)")
plt.plot(amb, cops, "tab:blue", marker="o", label="Cycle COP (R-32, Helmholtz)")
plt.xlabel("Ambient temperature [°C]")
plt.ylabel("COP (cooling)")
plt.title("R-32 ideal single-stage VC cycle: COP vs ambient\n"
          "T_evap = -29 °C, condenser approach 9 °C, eta_isen ~ 1, SH=SC=0 (Shridhar)")
plt.legend()
plt.grid(True, ls=":")
plt.tight_layout()
png_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "cop_vs_ambient_r32.png")
plt.savefig(png_path, dpi=150)
print(f"wrote {png_path}")
