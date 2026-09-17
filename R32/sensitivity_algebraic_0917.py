"""
sensitivity_algebraic_0917.py -- the full Tc/Pc sensitivity sweep, run through
the ALGEBRAIC cycle (cycle_algebraic_0917.py) instead of the IDAES flowsheet,
with the deviation against the Linde datasheet recorded on every row.

Two things this produces that the IDAES sweep could not:

  1. A COP at EVERY grid point. The IDAES sweep
     (phase6_sensitivity_NIST.xlsx) reports 251 of 900 case-ambient solves
     as non-converged, including every single case at Tc >= +3%. The
     algebraic cycle is sequential -- one 1-D root find, no simultaneous
     solve -- so it is deterministic and returns an answer unless the
     physics genuinely forbids one (e.g. a condenser temperature at or
     above the perturbed Tc, where no saturation state exists).

  2. Paired columns: for each perturbed (Tc, Pc), both the resulting COP
     and that same parameter set's deviation from the Linde datasheet.
     This is what turns the sensitivity from "COP per % of Tc" -- which
     the collaborators cannot measure for themselves -- into "COP per unit
     of property-prediction error", which they can, since Linde is public.

The Linde deviation is computed twice per case:
    full   : all 20 datasheet points, -130 to 78 C
    window : only points inside the cycle's own operating range,
             T_WINDOW below. The full-range metric is dominated by regions
             the cycle never visits, including the near-critical points at
             58-78 C where the cubic EoS is known to be weakest (see the
             08/26 breadcrumb note). n_points is written out for both so
             the sample size behind each number is visible.

NOTE on the Linde table: taken as-entered from the datasheet, per Shilpa's
instruction on 2026-09-17. The 10 C entry (10.065 bar) disagrees with
CoolProp (11.0691 bar) by 9.07% while its neighbours match to 0.00%, and
it sits inside the operating window -- so it contributes roughly the same
fixed error to every method's windowed P_MAPE. Not corrected here; the
n_points columns let it be excluded downstream if wanted.

Enthalpy/entropy deviations use the IIR anchoring the datasheet is on:
each method is offset so its own saturated liquid at 0 C reads
h = 200 kJ/kg and s = 1.0 kJ/kg/K, matching compare_cp_methods_*.py.
COP itself needs no anchoring -- it is a ratio of enthalpy differences.

Output: sensitivity_algebraic_0917.xlsx, three sheets (Tc_only, Pc_only,
Tc_Pc_grid), long format -- one row per (case, ambient), case-level fields
repeated. No header rows, no "not feasible" strings; a case that cannot be
computed gets COP = None and a populated 'error' column.

Author: Shilpa Narasimhan. Support: Claude AI.
"""

import os
import numpy as np
import pandas as pd

from cycle_algebraic_0917 import Fluid, run_cycle, NIST, GCGP, MW

# =============================================================================
# Config
# =============================================================================

DEV_PCTS = [-10, -5, -3, -1, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, 1, 3, 5, 10]
AMBIENTS = [10, 15, 20, 25]
T_EVAP_C = -29.0
COND_APPROACH = 9.0

# Operating window for the "window" deviation metric. -29 C is the
# evaporator setpoint; the condenser sits at T_amb + 9, so 19-34 C over
# the ambient range. Widened slightly to -30/+40 to land on datasheet
# temperatures rather than between them.
T_WINDOW = (-30.0, 40.0)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_XLSX = os.path.join(HERE, "sensitivity_algebraic_0917.xlsx")

# Linde datasheet, entered from the datasheet and kept as-is.
# (T [C], P [bar], hf [kJ/kg], hg [kJ/kg], sf [kJ/kg/K], sg [kJ/kg/K])
LINDE_SAT = [
    (-130, 0.001312,  -8.26, 448.77, -0.028, 3.165),
    (-110, 0.014525,  23.20, 461.86,  0.178, 2.867),
    (-90,  0.07556,   54.42, 474.61,  0.359, 2.653),
    (-70,  0.36067,   85.66, 486.57,  0.520, 2.494),
    (-50,  1.014,    117.22, 497.27,  0.668, 2.371),
    (-30,  2.7344,   149.45, 506.27,  0.806, 2.274),
    (-10,  5.8263,   182.76, 513.02,  0.937, 2.192),
    (0,    8.131,    200.00, 515.30,  1.000, 2.154),
    (10,   10.065,   217.74, 516.66,  1.063, 2.119),
    (20,   14.746,   236.12, 516.90,  1.125, 2.083),
    (30,   19.275,   255.32, 515.72,  1.188, 2.047),
    (40,   24.783,   275.61, 512.71,  1.252, 2.009),
    (50,   31.412,   297.49, 507.10,  1.318, 1.967),
    (58,   37.635,   316.75, 499.82,  1.375, 1.928),
    (62,   41.089,   327.30, 494.76,  1.405, 1.905),
    (66,   44.793,   338.78, 488.26,  1.438, 1.879),
    (70,   48.768,   351.73, 479.52,  1.474, 1.846),
    (74,   53.046,   367.53, 466.41,  1.518, 1.803),
    (76,   55.315,   378.03, 455.86,  1.547, 1.770),
    (78,   57.697,   400.38, 428.90,  1.610, 1.691),
]

IIR_H0 = 200.0    # kJ/kg, saturated liquid at 0 C
IIR_S0 = 1.0      # kJ/kg/K


# =============================================================================
# Linde deviation
# =============================================================================

def linde_deviation(fluid, T_lo=None, T_hi=None):
    """Saturation-property deviation of `fluid` from the Linde datasheet.

    Returns P_MAPE [%], hf/hg MAE [kJ/kg], sf/sg MAE [kJ/kg/K] and the
    number of datasheet points actually used. Points above the fluid's own
    Tc are skipped (no saturation state exists there) -- which is why
    n_points matters: a perturbed fluid with a low Tc is scored on fewer
    points than one with a high Tc, and the metrics are not comparable
    across very different n.

    h and s are anchored to the IIR convention the datasheet uses: the
    fluid's own saturated liquid at 0 C is offset to read exactly
    200 kJ/kg and 1.0 kJ/kg/K.
    """
    blank = dict(P_MAPE=np.nan, hf_MAE=np.nan, hg_MAE=np.nan,
                 sf_MAE=np.nan, sg_MAE=np.nan, n_points=0)

    # --- IIR anchor at 0 C ------------------------------------------------
    T0 = 273.15
    if T0 >= fluid.Tc:
        return blank
    try:
        P0, h0, s0 = fluid.sat_liquid(T0)
    except Exception:
        return blank
    h_off = IIR_H0*1e3*MW - h0        # J/mol
    s_off = IIR_S0*1e3*MW - s0        # J/mol/K

    dP, dhf, dhg, dsf, dsg = [], [], [], [], []
    for T_C, P_ref, hf_ref, hg_ref, sf_ref, sg_ref in LINDE_SAT:
        if T_lo is not None and T_C < T_lo:
            continue
        if T_hi is not None and T_C > T_hi:
            continue
        T = T_C + 273.15
        if T >= fluid.Tc:
            continue
        try:
            P = fluid.psat(T)
            hf = (fluid.h_liq(T, P) + h_off)/MW/1e3      # kJ/kg
            hg = (fluid.h_vap(T, P) + h_off)/MW/1e3
            sf = (fluid.s_liq(T, P) + s_off)/MW/1e3      # kJ/kg/K
            sg = (fluid.s_vap(T, P) + s_off)/MW/1e3
        except Exception:
            continue
        if not np.isfinite([P, hf, hg, sf, sg]).all():
            continue
        dP.append(abs(P/1e5 - P_ref)/P_ref)
        dhf.append(abs(hf - hf_ref))
        dhg.append(abs(hg - hg_ref))
        dsf.append(abs(sf - sf_ref))
        dsg.append(abs(sg - sg_ref))

    if not dP:
        return blank
    return dict(P_MAPE=100.0*float(np.mean(dP)),
                hf_MAE=float(np.mean(dhf)), hg_MAE=float(np.mean(dhg)),
                sf_MAE=float(np.mean(dsf)), sg_MAE=float(np.mean(dsg)),
                n_points=len(dP))


# =============================================================================
# Sweep
# =============================================================================

def make_case(tc_dev, pc_dev, base=NIST):
    p = dict(base)
    p["Tc"] = base["Tc"]*(1 + tc_dev/100.0)
    p["Pc"] = base["Pc"]*(1 + pc_dev/100.0)
    return Fluid(f"Tc{tc_dev:+.2f}_Pc{pc_dev:+.2f}", **p)


def rows_for_case(sweep, tc_dev, pc_dev):
    """One row per ambient. Case-level fields (the deviation metrics, which
    depend only on the parameter set) are repeated on each row so the sheet
    is tidy long-format and needs no joining."""
    f = make_case(tc_dev, pc_dev)
    full = linde_deviation(f)
    wind = linde_deviation(f, *T_WINDOW)

    rows = []
    for Tamb in AMBIENTS:
        r = run_cycle(f, T_evap_C=T_EVAP_C, T_amb_C=Tamb,
                      cond_approach=COND_APPROACH)
        row = {
            "sweep": sweep,
            "Tc_dev_pct": tc_dev, "Pc_dev_pct": pc_dev,
            "Tc_K": f.Tc, "Tc_C": f.Tc - 273.15, "Pc_Pa": f.Pc,
            "omega": f.omega,
            "T_amb_C": Tamb,
            "COP": None if r.get("error") else r["COP"],
            "error": r.get("error"),
            "P_evap_bar": None if r.get("error") else r["P_evap"]/1e5,
            "P_cond_bar": None if r.get("error") else r["P_cond"]/1e5,
            "pressure_ratio": None if r.get("error") else r["pressure_ratio"],
            "q_evap_kJkg": None if r.get("error") else r["q_evap"]/MW/1e3,
            "w_comp_kJkg": None if r.get("error") else r["w_comp"]/MW/1e3,
        }
        row.update({f"full_{k}": v for k, v in full.items()})
        row.update({f"wind_{k}": v for k, v in wind.items()})
        rows.append(row)
    return rows


if __name__ == "__main__":
    tc_only, pc_only, grid = [], [], []

    print("running Tc_only ...")
    for d in DEV_PCTS:
        tc_only += rows_for_case("Tc_only", d, 0)

    print("running Pc_only ...")
    for d in DEV_PCTS:
        pc_only += rows_for_case("Pc_only", 0, d)

    print(f"running Tc_Pc_grid ({len(DEV_PCTS)}x{len(DEV_PCTS)}) ...")
    for dtc in DEV_PCTS:
        for dpc in DEV_PCTS:
            grid += rows_for_case("Tc_Pc_grid", dtc, dpc)

    dfs = {"Tc_only": pd.DataFrame(tc_only),
           "Pc_only": pd.DataFrame(pc_only),
           "Tc_Pc_grid": pd.DataFrame(grid)}

    with pd.ExcelWriter(OUT_XLSX) as w:
        for name, df in dfs.items():
            df.to_excel(w, sheet_name=name, index=False)

    print(f"\nwrote {OUT_XLSX}")
    for name, df in dfs.items():
        ok = df["COP"].notna().sum()
        print(f"  {name:<12} {len(df):>5} rows   COP computed: {ok}/{len(df)} "
              f"({100.0*ok/len(df):.1f}%)")

    # --- quick look: does COP error track Linde deviation? ----------------
    g = dfs["Tc_Pc_grid"]
    base = g[(g.Tc_dev_pct == 0) & (g.Pc_dev_pct == 0)].set_index("T_amb_C")["COP"]
    sub = g[g["COP"].notna()].copy()
    sub["COP_err_pct"] = sub.apply(
        lambda r: 100.0*(r["COP"] - base[r["T_amb_C"]])/base[r["T_amb_C"]], axis=1)
    print("\n  correlation of |COP error vs baseline| with each deviation metric")
    print("  (Pearson r over the full grid, all ambients pooled)")
    for col in ["full_P_MAPE", "wind_P_MAPE", "full_hg_MAE", "wind_hg_MAE",
                "full_sg_MAE", "wind_sg_MAE"]:
        d = sub[[col, "COP_err_pct"]].dropna()
        if len(d) > 2:
            r = np.corrcoef(d[col], d["COP_err_pct"].abs())[0, 1]
            print(f"    {col:<16} r = {r:+.3f}   (n = {len(d)})")
    print()
