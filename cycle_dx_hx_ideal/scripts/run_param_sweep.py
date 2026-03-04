# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Parameter sweep entrypoint for the HX-ideal variant.

Context
-------
This script provides a uniform CSV output contract across variants.
It runs three sweep families:
1) COP vs evaporator saturation temperature
2) COP vs condenser saturation temperature
3) COP/capacity vs global UA scaling factor
"""

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from cycle_dx_hx_ideal.config import CycleConfig
from cycle_dx_hx_ideal.cycle_model import solve_cycle_point


def _row_from_result(sweep_type, sweep_value, cfg, result):
    row = {
        "sweep_type": sweep_type,
        "sweep_value": sweep_value,
        "fluid": result.fluid,
        "T_evap_sat": result.t_evap_sat_c,
        "T_cond_sat": result.t_cond_sat_c,
        "UA_evap_total": cfg.UA_evap_total,
        "UA_cond_total": cfg.UA_cond_total,
        "UA_evap_tp": cfg.UA_evap_tp,
        "UA_evap_sh": cfg.UA_evap_sh,
        "UA_cond_ds": cfg.UA_cond_ds,
        "UA_cond_tp": cfg.UA_cond_tp,
        "UA_cond_sc": cfg.UA_cond_sc,
        "Q_evap": result.q_evap_w,
        "Q_cond": result.q_cond_w,
        "W_comp": result.w_comp_w,
        "COP": result.cop,
        "m_dot_ref": result.m_dot_ref,
        "p_evap_pa": result.p_evap_pa,
        "p_cond_pa": result.p_cond_pa,
    }
    for k, v in result.diagnostics.items():
        row[f"diag_{k}"] = v
    return row


def main():
    parser = argparse.ArgumentParser(description="Run standardized cycle sweeps.")
    parser.add_argument("--fluid", default="R134a", help='e.g. "R134a" or "R1234ze(E)"')
    parser.add_argument("--out", default="results.csv", help="Output CSV path")
    args = parser.parse_args()

    cfg0 = CycleConfig()
    rows = []

    t_cond_fixed = 35.0
    for t_evap in [-40, -35, -30, -25, -20]:
        res = solve_cycle_point(args.fluid, t_evap, t_cond_fixed, cfg0, t_evap + 10.0, t_cond_fixed - 10.0)
        rows.append(_row_from_result("cop_vs_T_evap_sat", t_evap, cfg0, res))

    t_evap_fixed = -30.0
    for ambient_c in [15, 20, 25, 30, 35, 40, 45]:
        t_cond = ambient_c + 10.0
        res = solve_cycle_point(args.fluid, t_evap_fixed, t_cond, cfg0, t_evap_fixed + 10.0, ambient_c)
        rows.append(_row_from_result("cop_vs_ambient_15_45", ambient_c, cfg0, res))

    for ua_scale in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]:
        cfg = cfg0.scaled_ua(ua_scale)
        res = solve_cycle_point(args.fluid, t_evap_fixed, 35.0, cfg, t_evap_fixed + 10.0, 25.0)
        rows.append(_row_from_result("capacity_vs_UA_scale", ua_scale, cfg, res))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({k for row in rows for k in row.keys()})
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()
