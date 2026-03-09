# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Single-state zoned HX playground for intuition and manual checks.

Context
-------
This script is intentionally simple and educational. It models one evaporator
(2-zone: two-phase + superheat) and one condenser (3-zone: desuperheat +
two-phase + subcool) using epsilon-NTU equations and prints all intermediate
quantities: NTU, epsilon, Q per zone, and enthalpy progression.

It does not use a property package. You provide pseudo-property numbers
(enthalpies and cp values) so you can quickly explore UA sensitivity.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


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


@dataclass
class ZoneResult:
    name: str
    ua_w_per_k: float
    ntu: float
    epsilon: float
    q_req_w: float
    q_max_w: float
    q_used_w: float
    limited: bool
    t_air_in_c: float
    t_air_out_c: float


def two_phase_zone(
    name: str,
    ua_w_per_k: float,
    m_air_kg_s: float,
    cp_air_j_kgk: float,
    t_air_in_c: float,
    t_sat_c: float,
    q_req_w: float,
    air_is_hot: bool,
) -> ZoneResult:
    """Solve a two-phase zone with Cr ~ 0 and epsilon = 1 - exp(-NTU)."""
    c_air = m_air_kg_s * cp_air_j_kgk
    ntu = ua_w_per_k / max(c_air, 1e-9)
    epsilon = 1.0 - math.exp(-ntu)
    driving = (t_air_in_c - t_sat_c) if air_is_hot else (t_sat_c - t_air_in_c)
    q_max = max(0.0, epsilon * c_air * driving)
    q_used = max(0.0, min(max(0.0, q_req_w), q_max))
    sign = -1.0 if air_is_hot else 1.0
    t_air_out = t_air_in_c + sign * q_used / max(c_air, 1e-9)
    return ZoneResult(
        name=name,
        ua_w_per_k=ua_w_per_k,
        ntu=ntu,
        epsilon=epsilon,
        q_req_w=q_req_w,
        q_max_w=q_max,
        q_used_w=q_used,
        limited=(q_used + 1e-9 < max(0.0, q_req_w)),
        t_air_in_c=t_air_in_c,
        t_air_out_c=t_air_out,
    )


def single_phase_zone(
    name: str,
    ua_w_per_k: float,
    m_air_kg_s: float,
    cp_air_j_kgk: float,
    m_ref_kg_s: float,
    cp_ref_j_kgk: float,
    t_air_in_c: float,
    t_ref_in_c: float,
    q_req_w: float,
    air_is_hot: bool,
) -> ZoneResult:
    """Solve a single-phase zone with full counterflow epsilon-NTU."""
    c_air = m_air_kg_s * cp_air_j_kgk
    c_ref = m_ref_kg_s * cp_ref_j_kgk
    c_hot = c_air if air_is_hot else c_ref
    c_cold = c_ref if air_is_hot else c_air
    c_min = max(min(c_hot, c_cold), 1e-9)
    c_max = max(c_hot, c_cold)
    cr = c_min / max(c_max, 1e-9)
    ntu = ua_w_per_k / c_min
    epsilon = eps_counterflow(ntu, cr)
    t_hot_in = t_air_in_c if air_is_hot else t_ref_in_c
    t_cold_in = t_ref_in_c if air_is_hot else t_air_in_c
    q_max = max(0.0, epsilon * c_min * (t_hot_in - t_cold_in))
    q_used = max(0.0, min(max(0.0, q_req_w), q_max))
    sign = -1.0 if air_is_hot else 1.0
    t_air_out = t_air_in_c + sign * q_used / max(c_air, 1e-9)
    return ZoneResult(
        name=name,
        ua_w_per_k=ua_w_per_k,
        ntu=ntu,
        epsilon=epsilon,
        q_req_w=q_req_w,
        q_max_w=q_max,
        q_used_w=q_used,
        limited=(q_used + 1e-9 < max(0.0, q_req_w)),
        t_air_in_c=t_air_in_c,
        t_air_out_c=t_air_out,
    )


def print_zone(z: ZoneResult) -> None:
    print(
        f"  {z.name:14s} UA={z.ua_w_per_k:8.1f} NTU={z.ntu:6.3f} eps={z.epsilon:6.3f} "
        f"Q_req={z.q_req_w:9.1f}W Q_max={z.q_max_w:9.1f}W Q_used={z.q_used_w:9.1f}W "
        f"limited={str(z.limited):5s} T_air:{z.t_air_in_c:6.2f}->{z.t_air_out_c:6.2f}C"
    )


def save_outputs(out_prefix: str, zone_rows: list[ZoneResult], air_points: list[tuple[str, float]]) -> None:
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
                "T_air_in_C",
                "T_air_out_C",
            ]
        )
        for z in zone_rows:
            w.writerow(
                [
                    z.name,
                    f"{z.ua_w_per_k:.6g}",
                    f"{z.ntu:.6g}",
                    f"{z.epsilon:.6g}",
                    f"{z.q_req_w:.6g}",
                    f"{z.q_max_w:.6g}",
                    f"{z.q_used_w:.6g}",
                    int(z.limited),
                    f"{z.t_air_in_c:.6g}",
                    f"{z.t_air_out_c:.6g}",
                ]
            )

    zones = [z.name for z in zone_rows]
    q_req = [z.q_req_w for z in zone_rows]
    q_max = [z.q_max_w for z in zone_rows]
    q_used = [z.q_used_w for z in zone_rows]
    x = list(range(len(zones)))
    width = 0.26

    fig, ax = plt.subplots(figsize=(9.0, 4.8), dpi=150)
    ax.bar([i - width for i in x], q_req, width=width, label="Q_req")
    ax.bar(x, q_max, width=width, label="Q_max")
    ax.bar([i + width for i in x], q_used, width=width, label="Q_used")
    ax.set_xticks(x)
    ax.set_xticklabels(zones, rotation=20)
    ax.set_ylabel("Heat duty (W)")
    ax.set_title("Single-State Zoned HX Duties")
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
    ax.set_title("Air Temperature Progression Through Zones")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(tair_png, bbox_inches="tight")
    plt.close(fig)

    print(f"\nSaved computation table: {csv_path}")
    print(f"Saved duty plot:        {duty_png}")
    print(f"Saved air-T plot:       {tair_png}")


def main() -> None:
    p = argparse.ArgumentParser(description="Single-state zoned HX playground.")
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
    p.add_argument("--out-prefix", type=str, default="diagnostics/hx_single_state_playground")
    a = p.parse_args()

    print("=== INPUTS ===")
    print(
        f"Ambient={a.ambient_c:.2f} C  ColdStorage={a.cold_storage_c:.2f} C  "
        f"T_evap_sat={a.tevap_sat_c:.2f} C  T_cond_sat={a.tcond_sat_c:.2f} C"
    )
    print(
        f"m_ref={a.m_ref:.4f} kg/s  m_air_evap={a.m_air_evap:.3f} kg/s  m_air_cond={a.m_air_cond:.3f} kg/s"
    )
    print(
        f"UA evap(tp/sh)=({a.ua_e_tp:.1f}, {a.ua_e_sh:.1f}) W/K  "
        f"UA cond(ds/tp/sc)=({a.ua_c_ds:.1f}, {a.ua_c_tp:.1f}, {a.ua_c_sc:.1f}) W/K"
    )

    # Evaporator: 2-zone
    print("\n=== EVAPORATOR (2-zone) ===")
    q_tp_req_e = max(0.0, a.m_ref * (a.h_g_evap - a.h4_in))
    z_e_tp = two_phase_zone(
        name="evap_tp",
        ua_w_per_k=a.ua_e_tp,
        m_air_kg_s=a.m_air_evap,
        cp_air_j_kgk=a.cp_air,
        t_air_in_c=a.cold_storage_c,
        t_sat_c=a.tevap_sat_c,
        q_req_w=q_tp_req_e,
        air_is_hot=True,
    )
    h_after_tp = a.h4_in + z_e_tp.q_used_w / max(a.m_ref, 1e-9)

    q_sh_req_e = max(0.0, a.m_ref * a.cp_vap * a.sh_target_k)
    z_e_sh = single_phase_zone(
        name="evap_sh",
        ua_w_per_k=a.ua_e_sh,
        m_air_kg_s=a.m_air_evap,
        cp_air_j_kgk=a.cp_air,
        m_ref_kg_s=a.m_ref,
        cp_ref_j_kgk=a.cp_vap,
        t_air_in_c=z_e_tp.t_air_out_c,
        t_ref_in_c=a.tevap_sat_c,
        q_req_w=q_sh_req_e,
        air_is_hot=True,
    )
    h1_out = h_after_tp + z_e_sh.q_used_w / max(a.m_ref, 1e-9)
    q_evap = z_e_tp.q_used_w + z_e_sh.q_used_w
    print_zone(z_e_tp)
    print_zone(z_e_sh)
    print(f"  h4_in={a.h4_in:9.1f} J/kg -> h1_out={h1_out:9.1f} J/kg  Q_evap={q_evap:9.1f} W")

    # Condenser: 3-zone
    print("\n=== CONDENSER (3-zone) ===")
    q_ds_req_c = max(0.0, a.m_ref * (a.h2_in - a.h_g_cond))
    z_c_ds = single_phase_zone(
        name="cond_ds",
        ua_w_per_k=a.ua_c_ds,
        m_air_kg_s=a.m_air_cond,
        cp_air_j_kgk=a.cp_air,
        m_ref_kg_s=a.m_ref,
        cp_ref_j_kgk=a.cp_vap,
        t_air_in_c=a.ambient_c,
        t_ref_in_c=a.tcond_sat_c + 15.0,  # simple proxy for hot superheated inlet
        q_req_w=q_ds_req_c,
        air_is_hot=False,
    )
    h_after_ds = a.h2_in - z_c_ds.q_used_w / max(a.m_ref, 1e-9)

    q_tp_req_c = max(0.0, a.m_ref * (a.h_g_cond - a.h_f_cond))
    z_c_tp = two_phase_zone(
        name="cond_tp",
        ua_w_per_k=a.ua_c_tp,
        m_air_kg_s=a.m_air_cond,
        cp_air_j_kgk=a.cp_air,
        t_air_in_c=z_c_ds.t_air_out_c,
        t_sat_c=a.tcond_sat_c,
        q_req_w=q_tp_req_c,
        air_is_hot=False,
    )
    h_after_tp_c = h_after_ds - z_c_tp.q_used_w / max(a.m_ref, 1e-9)

    q_sc_req_c = max(0.0, a.m_ref * a.cp_liq * a.sc_target_k)
    z_c_sc = single_phase_zone(
        name="cond_sc",
        ua_w_per_k=a.ua_c_sc,
        m_air_kg_s=a.m_air_cond,
        cp_air_j_kgk=a.cp_air,
        m_ref_kg_s=a.m_ref,
        cp_ref_j_kgk=a.cp_liq,
        t_air_in_c=z_c_tp.t_air_out_c,
        t_ref_in_c=a.tcond_sat_c,
        q_req_w=q_sc_req_c,
        air_is_hot=False,
    )
    h3_out = h_after_tp_c - z_c_sc.q_used_w / max(a.m_ref, 1e-9)
    q_cond = z_c_ds.q_used_w + z_c_tp.q_used_w + z_c_sc.q_used_w
    print_zone(z_c_ds)
    print_zone(z_c_tp)
    print_zone(z_c_sc)
    print(f"  h2_in={a.h2_in:9.1f} J/kg -> h3_out={h3_out:9.1f} J/kg  Q_cond={q_cond:9.1f} W")

    print("\n=== QUICK CHECKS ===")
    delta_t_evap = a.cold_storage_c - a.tevap_sat_c
    delta_t_cond = a.tcond_sat_c - a.ambient_c
    print(f"Approach evap (T_cold - T_sat_evap): {delta_t_evap:.2f} C")
    print(f"Approach cond (T_sat_cond - T_ambient): {delta_t_cond:.2f} C")
    print("Tip: vary UA values and watch NTU/epsilon/Q_used changes by zone.")

    save_outputs(
        out_prefix=a.out_prefix,
        zone_rows=[z_e_tp, z_e_sh, z_c_ds, z_c_tp, z_c_sc],
        air_points=[
            ("evap_in", a.cold_storage_c),
            ("evap_after_tp", z_e_tp.t_air_out_c),
            ("evap_out", z_e_sh.t_air_out_c),
            ("cond_in", a.ambient_c),
            ("cond_after_ds", z_c_ds.t_air_out_c),
            ("cond_after_tp", z_c_tp.t_air_out_c),
            ("cond_out", z_c_sc.t_air_out_c),
        ],
    )


if __name__ == "__main__":
    main()
