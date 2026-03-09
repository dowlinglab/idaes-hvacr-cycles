# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Single-state zoned HX playground using IDAES/Pyomo.

Context
-------
This script mirrors the lightweight UA playground but formulates the
problem on an IDAES FlowsheetBlock and solves it with an NLP/MIP solver.
It is intentionally pseudo-physical (no refrigerant property package);
you provide pseudo-property inputs to understand UA->NTU->epsilon->Q
interactions for a single operating state.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    NonNegativeReals,
    Objective,
    Param,
    Var,
    maximize,
    value,
)


def eps_counterflow(ntu: float, cr: float) -> float:
    """Counterflow effectiveness for 0 <= cr <= 1 and ntu >= 0."""
    ntu = max(0.0, ntu)
    cr = min(max(cr, 0.0), 1.0)
    if abs(1.0 - cr) < 1e-9:
        return ntu / (1.0 + ntu)
    den = 1.0 - cr * math.exp(-ntu * (1.0 - cr))
    if abs(den) < 1e-12:
        return 0.0
    return (1.0 - math.exp(-ntu * (1.0 - cr))) / den


def add_q_bounds(fs, name: str, q_req: float, q_max: float):
    """Create zone duty var and cap constraints q <= q_req, q <= q_max."""
    q = Var(domain=NonNegativeReals, initialize=max(0.0, min(q_req, q_max)))
    setattr(fs, f"q_{name}", q)
    setattr(fs, f"q_{name}_req", Param(initialize=max(0.0, q_req), mutable=True))
    setattr(fs, f"q_{name}_max", Param(initialize=max(0.0, q_max), mutable=True))

    def _c_req(_):
        return q <= getattr(fs, f"q_{name}_req")

    def _c_max(_):
        return q <= getattr(fs, f"q_{name}_max")

    setattr(fs, f"q_{name}_req_con", Constraint(rule=_c_req))
    setattr(fs, f"q_{name}_max_con", Constraint(rule=_c_max))


def save_outputs(out_prefix: str, zone_rows: list[dict], air_points: list[tuple[str, float]]) -> None:
    """Save CSV + two PNG plots for easier interpretation."""
    out_base = Path(out_prefix)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_base.with_suffix(".csv")
    duty_png = out_base.with_name(out_base.name + "_duty.png")
    tair_png = out_base.with_name(out_base.name + "_air_t_profile.png")

    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "zone",
                "UA_W_per_K",
                "NTU",
                "epsilon",
                "Q_req_W",
                "Q_max_W",
                "Q_used_W",
                "limited",
            ]
        )
        for r in zone_rows:
            w.writerow(
                [
                    r["zone"],
                    f"{r['ua']:.6g}",
                    f"{r['ntu']:.6g}",
                    f"{r['eps']:.6g}",
                    f"{r['q_req']:.6g}",
                    f"{r['q_max']:.6g}",
                    f"{r['q_used']:.6g}",
                    int(r["q_used"] + 1e-9 < max(0.0, r["q_req"])),
                ]
            )

    zones = [r["zone"] for r in zone_rows]
    q_req = [r["q_req"] for r in zone_rows]
    q_max = [r["q_max"] for r in zone_rows]
    q_used = [r["q_used"] for r in zone_rows]
    x = list(range(len(zones)))
    width = 0.26

    fig, ax = plt.subplots(figsize=(9.0, 4.8), dpi=150)
    ax.bar([i - width for i in x], q_req, width=width, label="Q_req")
    ax.bar(x, q_max, width=width, label="Q_max")
    ax.bar([i + width for i in x], q_used, width=width, label="Q_used")
    ax.set_xticks(x)
    ax.set_xticklabels(zones, rotation=20)
    ax.set_ylabel("Heat duty (W)")
    ax.set_title("IDAES Single-State Zoned HX Duties")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(duty_png, bbox_inches="tight")
    plt.close(fig)

    labels = [k for k, _ in air_points]
    temps = [v for _, v in air_points]
    fig, ax = plt.subplots(figsize=(9.0, 4.4), dpi=150)
    ax.plot(labels, temps, marker="o", linewidth=2.0)
    ax.set_ylabel("Air temperature (C)")
    ax.set_title("IDAES Air Temperature Progression Through Zones")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(tair_png, bbox_inches="tight")
    plt.close(fig)

    print(f"\nSaved computation table: {csv_path}")
    print(f"Saved duty plot:        {duty_png}")
    print(f"Saved air-T plot:       {tair_png}")


def main() -> None:
    p = argparse.ArgumentParser(description="IDAES single-state UA playground.")
    p.add_argument("--ambient-c", type=float, default=20.0)
    p.add_argument("--cold-storage-c", type=float, default=-20.0)
    p.add_argument("--tevap-sat-c", type=float, default=-30.0)
    p.add_argument("--tcond-sat-c", type=float, default=30.0)
    p.add_argument("--m-ref", type=float, default=0.75, help="kg/s")
    p.add_argument("--m-air-evap", type=float, default=1.2, help="kg/s")
    p.add_argument("--m-air-cond", type=float, default=1.5, help="kg/s")
    p.add_argument("--cp-air", type=float, default=1006.0, help="J/kg-K")
    p.add_argument("--cp-vap", type=float, default=1200.0, help="J/kg-K")
    p.add_argument("--cp-liq", type=float, default=1400.0, help="J/kg-K")
    p.add_argument("--h4-in", type=float, default=240e3, help="J/kg at evap inlet")
    p.add_argument("--h-g-evap", type=float, default=390e3, help="J/kg sat vapor at evap P")
    p.add_argument("--h2-in", type=float, default=470e3, help="J/kg at cond inlet")
    p.add_argument("--h-g-cond", type=float, default=430e3, help="J/kg sat vapor at cond P")
    p.add_argument("--h-f-cond", type=float, default=260e3, help="J/kg sat liquid at cond P")
    p.add_argument("--sh-target-k", type=float, default=5.0)
    p.add_argument("--sc-target-k", type=float, default=3.0)
    p.add_argument("--ua-e-tp", type=float, default=220.0, help="W/K")
    p.add_argument("--ua-e-sh", type=float, default=110.0, help="W/K")
    p.add_argument("--ua-c-ds", type=float, default=110.0, help="W/K")
    p.add_argument("--ua-c-tp", type=float, default=220.0, help="W/K")
    p.add_argument("--ua-c-sc", type=float, default=66.0, help="W/K")
    p.add_argument("--out-prefix", type=str, default="diagnostics/hx_single_state_idaes_playground")
    a = p.parse_args()

    # Capacity rates
    c_air_e = a.m_air_evap * a.cp_air
    c_air_c = a.m_air_cond * a.cp_air
    c_ref_v = a.m_ref * a.cp_vap
    c_ref_l = a.m_ref * a.cp_liq

    # Zone q_req values
    q_e_tp_req = max(0.0, a.m_ref * (a.h_g_evap - a.h4_in))
    q_e_sh_req = max(0.0, a.m_ref * a.cp_vap * a.sh_target_k)
    q_c_ds_req = max(0.0, a.m_ref * (a.h2_in - a.h_g_cond))
    q_c_tp_req = max(0.0, a.m_ref * (a.h_g_cond - a.h_f_cond))
    q_c_sc_req = max(0.0, a.m_ref * a.cp_liq * a.sc_target_k)

    # Zone epsilon/NTU and q_max
    ntu_e_tp = a.ua_e_tp / max(c_air_e, 1e-9)
    eps_e_tp = 1.0 - math.exp(-ntu_e_tp)  # Cr ~ 0
    q_e_tp_max = max(0.0, eps_e_tp * c_air_e * (a.cold_storage_c - a.tevap_sat_c))

    cmin_e_sh = min(c_air_e, c_ref_v)
    cr_e_sh = cmin_e_sh / max(c_air_e, c_ref_v)
    ntu_e_sh = a.ua_e_sh / max(cmin_e_sh, 1e-9)
    eps_e_sh = eps_counterflow(ntu_e_sh, cr_e_sh)

    # ds uses a simple superheated ref inlet proxy at Tcond+15 C
    t_ref_ds_in = a.tcond_sat_c + 15.0
    cmin_c_ds = min(c_air_c, c_ref_v)
    cr_c_ds = cmin_c_ds / max(c_air_c, c_ref_v)
    ntu_c_ds = a.ua_c_ds / max(cmin_c_ds, 1e-9)
    eps_c_ds = eps_counterflow(ntu_c_ds, cr_c_ds)
    q_c_ds_max = max(0.0, eps_c_ds * cmin_c_ds * (t_ref_ds_in - a.ambient_c))

    ntu_c_tp = a.ua_c_tp / max(c_air_c, 1e-9)
    eps_c_tp = 1.0 - math.exp(-ntu_c_tp)  # Cr ~ 0
    # conservative max with inlet air at ambient (upper zone coupling appears in constraints)
    q_c_tp_max = max(0.0, eps_c_tp * c_air_c * (a.tcond_sat_c - a.ambient_c))

    cmin_c_sc = min(c_air_c, c_ref_l)
    cr_c_sc = cmin_c_sc / max(c_air_c, c_ref_l)
    ntu_c_sc = a.ua_c_sc / max(cmin_c_sc, 1e-9)
    eps_c_sc = eps_counterflow(ntu_c_sc, cr_c_sc)
    q_c_sc_max = max(0.0, eps_c_sc * cmin_c_sc * (a.tcond_sat_c - a.ambient_c))

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    fs = m.fs

    add_q_bounds(fs, "e_tp", q_e_tp_req, q_e_tp_max)
    add_q_bounds(fs, "e_sh", q_e_sh_req, 1e9)  # capped later with temp-updated q_max
    add_q_bounds(fs, "c_ds", q_c_ds_req, q_c_ds_max)
    add_q_bounds(fs, "c_tp", q_c_tp_req, q_c_tp_max)
    add_q_bounds(fs, "c_sc", q_c_sc_req, q_c_sc_max)

    # Air outlet coupling for more realistic zoned sequence
    fs.t_air_e1 = Var(initialize=a.cold_storage_c)
    fs.t_air_e2 = Var(initialize=a.cold_storage_c)
    fs.t_air_c1 = Var(initialize=a.ambient_c)
    fs.t_air_c2 = Var(initialize=a.ambient_c)
    fs.t_air_c3 = Var(initialize=a.ambient_c)

    fs.ev_air1 = Constraint(expr=fs.t_air_e1 == a.cold_storage_c - fs.q_e_tp / max(c_air_e, 1e-9))
    fs.ev_air2 = Constraint(expr=fs.t_air_e2 == fs.t_air_e1 - fs.q_e_sh / max(c_air_e, 1e-9))
    fs.cd_air1 = Constraint(expr=fs.t_air_c1 == a.ambient_c + fs.q_c_ds / max(c_air_c, 1e-9))
    fs.cd_air2 = Constraint(expr=fs.t_air_c2 == fs.t_air_c1 + fs.q_c_tp / max(c_air_c, 1e-9))
    fs.cd_air3 = Constraint(expr=fs.t_air_c3 == fs.t_air_c2 + fs.q_c_sc / max(c_air_c, 1e-9))

    # Dynamic q_max constraints for zones whose driving temp depends on upstream air state.
    fs.e_sh_qmax = Constraint(
        expr=fs.q_e_sh <= max(0.0, eps_e_sh * cmin_e_sh) * (fs.t_air_e1 - a.tevap_sat_c)
    )
    fs.c_tp_qmax = Constraint(
        expr=fs.q_c_tp <= max(0.0, eps_c_tp * c_air_c) * (a.tcond_sat_c - fs.t_air_c1)
    )
    fs.c_sc_qmax = Constraint(
        expr=fs.q_c_sc <= max(0.0, eps_c_sc * cmin_c_sc) * (a.tcond_sat_c - fs.t_air_c2)
    )

    fs.obj = Objective(
        expr=fs.q_e_tp + fs.q_e_sh + fs.q_c_ds + fs.q_c_tp + fs.q_c_sc,
        sense=maximize,
    )

    solver = get_solver()
    res = solver.solve(m, tee=False)

    q_e_tp = value(fs.q_e_tp)
    q_e_sh = value(fs.q_e_sh)
    q_c_ds = value(fs.q_c_ds)
    q_c_tp = value(fs.q_c_tp)
    q_c_sc = value(fs.q_c_sc)

    h1_out = a.h4_in + (q_e_tp + q_e_sh) / max(a.m_ref, 1e-9)
    h3_out = a.h2_in - (q_c_ds + q_c_tp + q_c_sc) / max(a.m_ref, 1e-9)

    print("=== IDAES UA PLAYGROUND (single state) ===")
    print(f"Solver status: {res.solver.status}, termination: {res.solver.termination_condition}")
    print(
        f"Ambient={a.ambient_c:.2f} C ColdStorage={a.cold_storage_c:.2f} C "
        f"T_evap_sat={a.tevap_sat_c:.2f} C T_cond_sat={a.tcond_sat_c:.2f} C"
    )
    print()
    print("Evaporator zones")
    print(
        f"  evap_tp: NTU={ntu_e_tp:.3f} eps={eps_e_tp:.3f} "
        f"Q_req={q_e_tp_req:.1f} Q_max={q_e_tp_max:.1f} Q_used={q_e_tp:.1f}"
    )
    print(
        f"  evap_sh: NTU={ntu_e_sh:.3f} eps={eps_e_sh:.3f} "
        f"Q_req={q_e_sh_req:.1f} Q_used={q_e_sh:.1f} (q_max from coupled constraint)"
    )
    print(f"  Air: {a.cold_storage_c:.2f} -> {value(fs.t_air_e1):.2f} -> {value(fs.t_air_e2):.2f} C")
    print(f"  h4_in={a.h4_in:.1f} -> h1_out={h1_out:.1f} J/kg, Q_evap={q_e_tp + q_e_sh:.1f} W")
    print()
    print("Condenser zones")
    print(
        f"  cond_ds: NTU={ntu_c_ds:.3f} eps={eps_c_ds:.3f} "
        f"Q_req={q_c_ds_req:.1f} Q_max={q_c_ds_max:.1f} Q_used={q_c_ds:.1f}"
    )
    print(
        f"  cond_tp: NTU={ntu_c_tp:.3f} eps={eps_c_tp:.3f} "
        f"Q_req={q_c_tp_req:.1f} Q_used={q_c_tp:.1f} (q_max from coupled constraint)"
    )
    print(
        f"  cond_sc: NTU={ntu_c_sc:.3f} eps={eps_c_sc:.3f} "
        f"Q_req={q_c_sc_req:.1f} Q_used={q_c_sc:.1f} (q_max from coupled constraint)"
    )
    print(f"  Air: {a.ambient_c:.2f} -> {value(fs.t_air_c1):.2f} -> {value(fs.t_air_c2):.2f} -> {value(fs.t_air_c3):.2f} C")
    print(f"  h2_in={a.h2_in:.1f} -> h3_out={h3_out:.1f} J/kg, Q_cond={q_c_ds + q_c_tp + q_c_sc:.1f} W")
    print()
    print("Approach checks")
    print(f"  DeltaT_evap = T_cold - T_evap_sat = {a.cold_storage_c - a.tevap_sat_c:.2f} C")
    print(f"  DeltaT_cond = T_cond_sat - T_ambient = {a.tcond_sat_c - a.ambient_c:.2f} C")

    q_e_sh_max = max(0.0, eps_e_sh * cmin_e_sh * (value(fs.t_air_e1) - a.tevap_sat_c))
    q_c_tp_max_dyn = max(0.0, eps_c_tp * c_air_c * (a.tcond_sat_c - value(fs.t_air_c1)))
    q_c_sc_max_dyn = max(0.0, eps_c_sc * cmin_c_sc * (a.tcond_sat_c - value(fs.t_air_c2)))
    save_outputs(
        out_prefix=a.out_prefix,
        zone_rows=[
            {"zone": "evap_tp", "ua": a.ua_e_tp, "ntu": ntu_e_tp, "eps": eps_e_tp, "q_req": q_e_tp_req, "q_max": q_e_tp_max, "q_used": q_e_tp},
            {"zone": "evap_sh", "ua": a.ua_e_sh, "ntu": ntu_e_sh, "eps": eps_e_sh, "q_req": q_e_sh_req, "q_max": q_e_sh_max, "q_used": q_e_sh},
            {"zone": "cond_ds", "ua": a.ua_c_ds, "ntu": ntu_c_ds, "eps": eps_c_ds, "q_req": q_c_ds_req, "q_max": q_c_ds_max, "q_used": q_c_ds},
            {"zone": "cond_tp", "ua": a.ua_c_tp, "ntu": ntu_c_tp, "eps": eps_c_tp, "q_req": q_c_tp_req, "q_max": q_c_tp_max_dyn, "q_used": q_c_tp},
            {"zone": "cond_sc", "ua": a.ua_c_sc, "ntu": ntu_c_sc, "eps": eps_c_sc, "q_req": q_c_sc_req, "q_max": q_c_sc_max_dyn, "q_used": q_c_sc},
        ],
        air_points=[
            ("evap_in", a.cold_storage_c),
            ("evap_after_tp", value(fs.t_air_e1)),
            ("evap_out", value(fs.t_air_e2)),
            ("cond_in", a.ambient_c),
            ("cond_after_ds", value(fs.t_air_c1)),
            ("cond_after_tp", value(fs.t_air_c2)),
            ("cond_out", value(fs.t_air_c3)),
        ],
    )


if __name__ == "__main__":
    main()
