#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-03
Purpose of file: Validate current R515A mixture model against published
saturation reference data at 277.6 K from NIST TN 2063 Table 7.
Dependencies: csv, pathlib, linear_model_codex, mixture_true_vle_copy
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-03
# Assumptions:
# - Uses R515A composition as w1(R1234ze)=0.88 mass fraction.
# - Uses NIST TN 2063 Table 7 saturation point at T=277.6 K.
# - Evaluates model at reference liquid/vapor densities and at attempted
#   bubble/dew equilibrium solve for diagnostics.
# TODO: Extend to multi-point saturation table validation once additional
#       R515A saturation rows are available in-repo.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from linear_model_codex import load_idaes_helmholtz_json, mw_from_json
from mixture_true_vle_copy import mix_state, solve_bubble_at_t, solve_dew_at_t, w1_to_x1


PSIA_PER_KPA = 0.145037737730
BTU_PER_KJ = 0.429922614


def _to_psia(kpa: float) -> float:
    """
    Purpose
    -------
    Convert pressure from kPa to psia.

    Inputs
    ------
    kpa : float [kPa]

    Outputs
    -------
    float [psia]

    Assumptions
    -----------
    Assumes absolute pressure values.

    Failure modes
    -------------
    No explicit failure handling; propagates NaN/Inf.

    References
    ----------
    Unit conversion constants.

    Numerical stability notes
    -------------------------
    Linear scale conversion; numerically stable.
    """
    return kpa * PSIA_PER_KPA


def _to_btulbm(kjkg: float) -> float:
    """
    Purpose
    -------
    Convert specific enthalpy from kJ/kg to Btu/lbm.

    Inputs
    ------
    kjkg : float [kJ/kg]

    Outputs
    -------
    float [Btu/lbm]

    Assumptions
    -----------
    Uses thermochemical Btu/lbm conversion.

    Failure modes
    -------------
    No explicit failure handling; propagates NaN/Inf.

    References
    ----------
    Unit conversion constants.

    Numerical stability notes
    -------------------------
    Linear scale conversion; numerically stable.
    """
    return kjkg * BTU_PER_KJ


def main() -> None:
    """
    Purpose
    -------
    Run R515A saturation validation at a published reference point and write
    machine-readable output.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    Writes CSV to verification/r515a_saturation_validation_nist_tn2063.csv.

    Assumptions
    -----------
    Composition is fixed at w1=0.88 (R1234ze mass fraction).

    Failure modes
    -------------
    Raises exceptions if model functions fail or files cannot be written.

    References
    ----------
    NIST TN 2063 Table 7 (saturation at 277.6 K).

    Numerical stability notes
    -------------------------
    No iterative method here except solver diagnostics from imported module.
    """
    # === SECTION: Reference Inputs ===
    # Rationale: Centralize published data and blend definition for traceability.
    t_k = 277.6
    w1 = 0.88
    ref_p_kpa = 252.31
    ref_rho_l_kgm3 = 1250.9
    ref_rho_v_kgm3 = 14.28
    ref_hfg_kjkg = 175.43

    # === SECTION: Load Fluids and Basis Conversion ===
    # Rationale: Convert published mass-basis densities to molar basis expected by EOS.
    d1 = load_idaes_helmholtz_json("r1234ze")
    d2 = load_idaes_helmholtz_json("r227ea")
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    x1 = w1_to_x1(w1, mw1, mw2)
    x2 = 1.0 - x1
    mw_mix = x1 * mw1 + x2 * mw2

    rho_l_molm3 = ref_rho_l_kgm3 / mw_mix
    rho_v_molm3 = ref_rho_v_kgm3 / mw_mix

    # === SECTION: Fixed-Density Saturation Check ===
    # Rationale: Compare model pressure and latent heat at published saturation densities.
    st_l = mix_state(d1, d2, t_k, rho_l_molm3, x1)
    st_v = mix_state(d1, d2, t_k, rho_v_molm3, x1)
    p_l_kpa = st_l.p_pa / 1e3
    p_v_kpa = st_v.p_pa / 1e3
    hfg_kjkg = ((st_v.h_jmol - st_l.h_jmol) / mw_mix) / 1e3

    # === SECTION: Bubble/Dew Solve Diagnostic ===
    # Rationale: Keep visibility into current equilibrium solver behavior.
    bubble = solve_bubble_at_t(d1, d2, t_k, x1, rho_l_molm3, rho_v_molm3, x1)
    dew = solve_dew_at_t(d1, d2, t_k, x1, rho_l_molm3, rho_v_molm3, x1)

    # === SECTION: Assemble Output Rows ===
    # Rationale: Single CSV artifact for quick audit and reproducibility.
    rows = [
        {
            "case": "reference_nist_tn2063_table7",
            "status": "REFERENCE",
            "T_K": t_k,
            "w1_mass": w1,
            "x1_mol": x1,
            "p_kPa": ref_p_kpa,
            "p_psia": _to_psia(ref_p_kpa),
            "rho_l_kgm3": ref_rho_l_kgm3,
            "rho_v_kgm3": ref_rho_v_kgm3,
            "hfg_kJkg": ref_hfg_kjkg,
            "hfg_Btu_lbm": _to_btulbm(ref_hfg_kjkg),
            "r_P": "",
            "r_mu": "",
            "notes": "NIST TN 2063 Table 7 reference value",
        },
        {
            "case": "model_fixed_comp_at_ref_rho_liq",
            "status": "MODEL",
            "T_K": t_k,
            "w1_mass": w1,
            "x1_mol": x1,
            "p_kPa": p_l_kpa,
            "p_psia": _to_psia(p_l_kpa),
            "rho_l_kgm3": st_l.rho_mass,
            "rho_v_kgm3": "",
            "hfg_kJkg": "",
            "hfg_Btu_lbm": "",
            "r_P": abs(p_l_kpa - ref_p_kpa) / max(1.0, ref_p_kpa),
            "r_mu": "",
            "notes": "Model liquid-side pressure at reference saturated liquid density",
        },
        {
            "case": "model_fixed_comp_at_ref_rho_vap",
            "status": "MODEL",
            "T_K": t_k,
            "w1_mass": w1,
            "x1_mol": x1,
            "p_kPa": p_v_kpa,
            "p_psia": _to_psia(p_v_kpa),
            "rho_l_kgm3": "",
            "rho_v_kgm3": st_v.rho_mass,
            "hfg_kJkg": "",
            "hfg_Btu_lbm": "",
            "r_P": abs(p_v_kpa - ref_p_kpa) / max(1.0, ref_p_kpa),
            "r_mu": "",
            "notes": "Model vapor-side pressure at reference saturated vapor density",
        },
        {
            "case": "model_fixed_comp_hfg_from_ref_rho",
            "status": "MODEL",
            "T_K": t_k,
            "w1_mass": w1,
            "x1_mol": x1,
            "p_kPa": "",
            "p_psia": "",
            "rho_l_kgm3": st_l.rho_mass,
            "rho_v_kgm3": st_v.rho_mass,
            "hfg_kJkg": hfg_kjkg,
            "hfg_Btu_lbm": _to_btulbm(hfg_kjkg),
            "r_P": "",
            "r_mu": "",
            "notes": "Model latent heat using reference densities and fixed composition",
        },
        {
            "case": "model_bubble_solver",
            "status": bubble["status"],
            "T_K": t_k,
            "w1_mass": w1,
            "x1_mol": x1,
            "p_kPa": bubble["P_Pa"] / 1e3,
            "p_psia": _to_psia(bubble["P_Pa"] / 1e3),
            "rho_l_kgm3": bubble["rho_l_molm3"] * mw_mix,
            "rho_v_kgm3": bubble["rho_v_molm3"] * mw_mix,
            "hfg_kJkg": (bubble["h_v_Jmol"] - bubble["h_l_Jmol"]) / mw_mix / 1e3,
            "hfg_Btu_lbm": _to_btulbm((bubble["h_v_Jmol"] - bubble["h_l_Jmol"]) / mw_mix / 1e3),
            "r_P": bubble["r_P"],
            "r_mu": bubble["r_mu"],
            "notes": bubble["notes"],
        },
        {
            "case": "model_dew_solver",
            "status": dew["status"],
            "T_K": t_k,
            "w1_mass": w1,
            "x1_mol": x1,
            "p_kPa": dew["P_Pa"] / 1e3,
            "p_psia": _to_psia(dew["P_Pa"] / 1e3),
            "rho_l_kgm3": dew["rho_l_molm3"] * mw_mix,
            "rho_v_kgm3": dew["rho_v_molm3"] * mw_mix,
            "hfg_kJkg": (dew["h_v_Jmol"] - dew["h_l_Jmol"]) / mw_mix / 1e3,
            "hfg_Btu_lbm": _to_btulbm((dew["h_v_Jmol"] - dew["h_l_Jmol"]) / mw_mix / 1e3),
            "r_P": dew["r_P"],
            "r_mu": dew["r_mu"],
            "notes": dew["notes"],
        },
    ]

    out = Path("verification/r515a_saturation_validation_nist_tn2063.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Saved validation CSV: {out}")


if __name__ == "__main__":
    main()
