#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Shilpa Narasimhan
Technical support: Codex (version GPT-5)
QA/Testing Responsibility: Shilpa
Creation date: 2026-03-03
Purpose of file: Run a consolidated R515A-oriented validation suite using
available published reference points in local sources.
Dependencies: csv, pathlib, linear_model_codex, mixture_true_vle_copy
Context reference: PROJECT_CONTEXT.md

Version: v0.1.0

# BREADCRUMB:
# Date: 2026-03-03
# Assumptions:
# - Uses Bell 2023 Table 13 check values available for R1234ze(E), R227ea,
#   and R1234ze(E)/227ea pair.
# - Uses NIST TN 2063 Table 7 saturation point for R515A at 277.6 K.
# - Composition held constant where required; no composition-derivative VLE terms
#   are inferred beyond current model implementation.
# TODO: Add multi-temperature R515A saturation table validation when additional
#       digitized reference points are provided.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from linear_model_codex import (  # noqa: E402
    BELL_2023_R1234ZE_R227EA,
    alphar_idaes_with_derivs,
    bell2023_Tred_vred,
    load_idaes_helmholtz_json,
    mixture_alpha0_alphar_derivs,
    mw_from_json,
)
from mixture_true_vle_copy import mix_state, solve_bubble_at_t, solve_dew_at_t, w1_to_x1  # noqa: E402


def _safe_rel_err(model: float, ref: float) -> float:
    """
    Purpose
    -------
    Compute relative error with safe denominator.

    Inputs
    ------
    model : float [property units]
    ref : float [same property units]

    Outputs
    -------
    float [fraction]

    Assumptions
    -----------
    Reference is finite.

    Failure modes
    -------------
    Propagates NaN for non-finite inputs.

    References
    ----------
    Standard relative error definition.

    Numerical stability notes
    -------------------------
    Denominator is clamped at 1e-30 to avoid divide-by-zero.
    """
    denom = max(abs(ref), 1e-30)
    return (model - ref) / denom


def bell_table13_checks() -> List[Dict[str, float | str]]:
    """
    Purpose
    -------
    Validate residual Helmholtz alpha_r against Bell 2023 Table 13 check values.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    list[dict]
      Rows containing reference and model alpha_r for selected check states.

    Assumptions
    -----------
    Table-13 values are copied verbatim from Bell 2023 PDF text.

    Failure modes
    -------------
    Raises if component parameter files are unavailable.

    References
    ----------
    Bell, J. Phys. Chem. Ref. Data 52, 013101 (2023), Table 13.

    Numerical stability notes
    -------------------------
    Uses direct analytic expressions without iterative solving.
    """
    rows: List[Dict[str, float | str]] = []

    # === SECTION: Pure-fluid check states ===
    # Rationale: Validate pure EOS residual mapping at published check points.
    pure_cases = [
        {"fluid": "r1234ze", "name": "R1234ZEE", "T_K": 478.0, "rho_molm3": 3432.0, "Tred_K": 382.513, "rhored_molm3": 4290.0, "ar_ref": -0.46340978447230},
        {"fluid": "r227ea", "name": "R227EA", "T_K": 469.0, "rho_molm3": 2796.0, "Tred_K": 374.900, "rhored_molm3": 3495.0, "ar_ref": -0.44238576197982},
    ]
    for c in pure_cases:
        d = load_idaes_helmholtz_json(c["fluid"])
        tau = c["Tred_K"] / c["T_K"]
        delta = c["rho_molm3"] / c["rhored_molm3"]
        ar_model, _, _ = alphar_idaes_with_derivs(d["eos"], tau, delta)
        rows.append(
            {
                "suite": "Bell2023_Table13",
                "case": c["name"],
                "property": "alpha_r",
                "ref_value": c["ar_ref"],
                "model_value": float(ar_model),
                "rel_err": _safe_rel_err(float(ar_model), c["ar_ref"]),
                "units": "dimensionless",
                "notes": "Pure-fluid check value",
            }
        )

    # === SECTION: Binary check state for R1234ze(E)/227ea ===
    # Rationale: Validate mixture residual contribution at published pair check point.
    d1 = load_idaes_helmholtz_json("r1234ze")
    d2 = load_idaes_helmholtz_json("r227ea")
    x1 = 0.4
    x2 = 0.6
    t_k = 470.0
    rho_mol = 3023.0
    ar_ref = -0.45378834770736

    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    tc1 = float(d1["basic"]["Tc"])
    tc2 = float(d2["basic"]["Tc"])
    rhoc1_mol = float(d1["basic"]["rhoc"]) / mw1
    rhoc2_mol = float(d2["basic"]["rhoc"]) / mw2
    vc1 = 1.0 / rhoc1_mol
    vc2 = 1.0 / rhoc2_mol

    tred, vred = bell2023_Tred_vred(x1, x2, tc1, tc2, vc1, vc2, BELL_2023_R1234ZE_R227EA)
    tau = tred / t_k
    delta = rho_mol * vred
    rho_red_mol = 1.0 / vred

    _, _, ar_model, _, _ = mixture_alpha0_alphar_derivs(
        d1=d1,
        d2=d2,
        x1=x1,
        x2=x2,
        tau=tau,
        delta=delta,
        Tred=tred,
        rho_red_mol=rho_red_mol,
        pair_key="r1234ze|r227ea",
    )

    rows.append(
        {
            "suite": "Bell2023_Table13",
            "case": "R1234ZEE/R227EA_z1_0.4",
            "property": "alpha_r",
            "ref_value": ar_ref,
            "model_value": float(ar_model),
            "rel_err": _safe_rel_err(float(ar_model), ar_ref),
            "units": "dimensionless",
            "notes": "Binary pair check value",
        }
    )

    return rows


def r515a_nist_saturation_checks() -> List[Dict[str, float | str]]:
    """
    Purpose
    -------
    Validate R515A saturation predictions against NIST TN 2063 reference point.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    list[dict]
      Rows of pressure and latent heat mismatches for bubble/dew solve outputs.

    Assumptions
    -----------
    R515A composition is approximated as w1=0.88 for R1234ze(E) mass fraction.

    Failure modes
    -------------
    Raises if solver/model calls fail.

    References
    ----------
    NIST TN 2063 Table 7 point at T=277.6 K.

    Numerical stability notes
    -------------------------
    Uses existing nonlinear solver; convergence status is preserved in notes.
    """
    rows: List[Dict[str, float | str]] = []

    # Reference
    t_k = 277.6
    w1 = 0.88
    p_ref_kpa = 252.31
    hfg_ref_kjkg = 175.43

    d1 = load_idaes_helmholtz_json("r1234ze")
    d2 = load_idaes_helmholtz_json("r227ea")
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    x1 = w1_to_x1(w1, mw1, mw2)
    x2 = 1.0 - x1
    mw_mix = x1 * mw1 + x2 * mw2

    # Seed from reported saturated densities
    rho_l_seed = 1250.9 / mw_mix
    rho_v_seed = 14.28 / mw_mix

    bubble = solve_bubble_at_t(d1, d2, t_k, x1, rho_l_seed, rho_v_seed, x1)
    dew = solve_dew_at_t(d1, d2, t_k, x1, rho_l_seed, rho_v_seed, x1)

    for label, r in [("bubble", bubble), ("dew", dew)]:
        p_model_kpa = float(r["P_Pa"]) / 1e3
        hfg_model_kjkg = (float(r["h_v_Jmol"]) - float(r["h_l_Jmol"])) / mw_mix / 1e3
        rows.append(
            {
                "suite": "NIST_TN2063_Table7",
                "case": f"R515A_{label}_T277.6K",
                "property": "pressure_kPa",
                "ref_value": p_ref_kpa,
                "model_value": p_model_kpa,
                "rel_err": _safe_rel_err(p_model_kpa, p_ref_kpa),
                "units": "kPa",
                "notes": f"status={r['status']}; r_P={r['r_P']:.3e}; r_mu={r['r_mu']:.3e}",
            }
        )
        rows.append(
            {
                "suite": "NIST_TN2063_Table7",
                "case": f"R515A_{label}_T277.6K",
                "property": "hfg_kJkg",
                "ref_value": hfg_ref_kjkg,
                "model_value": hfg_model_kjkg,
                "rel_err": _safe_rel_err(hfg_model_kjkg, hfg_ref_kjkg),
                "units": "kJ/kg",
                "notes": f"status={r['status']}; r_P={r['r_P']:.3e}; r_mu={r['r_mu']:.3e}",
            }
        )

    return rows


def main() -> None:
    """
    Purpose
    -------
    Execute consolidated reference-suite validation and save CSV summary.

    Inputs
    ------
    None [unitless]

    Outputs
    -------
    Writes verification/r515a_reference_suite_validation.csv.

    Assumptions
    -----------
    Required fluid parameter files are present in IDAES installation.

    Failure modes
    -------------
    Raises exceptions on missing dependencies/data.

    References
    ----------
    Bell 2023 Table 13 and NIST TN 2063 Table 7.

    Numerical stability notes
    -------------------------
    No additional iterative method beyond imported VLE solver.
    """
    rows = []
    rows.extend(bell_table13_checks())
    rows.extend(r515a_nist_saturation_checks())

    out = Path("verification/r515a_reference_suite_validation.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["suite", "case", "property", "ref_value", "model_value", "rel_err", "units", "notes"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    # Console summary
    abs_rel = [abs(float(r["rel_err"])) for r in rows if isinstance(r["rel_err"], (float, int)) and math.isfinite(float(r["rel_err"]))]
    print(f"Saved: {out}")
    print(f"Points: {len(rows)}")
    if abs_rel:
        print(f"Mean |rel_err|: {sum(abs_rel)/len(abs_rel):.6e}")
        print(f"Max  |rel_err|: {max(abs_rel):.6e}")


if __name__ == "__main__":
    main()
