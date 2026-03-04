# Author: Shilpa Narasimhan (project owner) with Codex implementation support.
# QA/testing and production validation are intentionally left to Shilpa Narasimhan.
"""Overlay previous Carnot/COP curves with new HX-variant COP results.

Context
-------
This utility overlays historical cold-storage Carnot/COP values produced by
`run_plr_cold_storage_*_copy.py` with COP results from the lumped-UA variant.
The comparison axis is ambient temperature over 15..45 C at fixed evaporator
saturation (-30 C).
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def _read_csv(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _to_float(v: str) -> float:
    return float(v)


def _previous_series(fluid_key: str):
    """Return ambient temperature grid, prior COP, and prior Carnot."""
    file_map = {
        "R134a": ROOT / "cop_vs_ambient_plr_cold_storage_r134a_copy.csv",
        "R1234ze(E)": ROOT / "cop_vs_ambient_plr_cold_storage_r1234zee_copy.csv",
    }
    rows = _read_csv(file_map[fluid_key])
    x = [_to_float(r["ambient_C"]) for r in rows]
    y_cop = [_to_float(r["cop_part_load"]) for r in rows]
    y_carnot = [_to_float(r["cop_carnot_ref"]) for r in rows]
    return x, y_cop, y_carnot


def _plr_hx_series(fluid_key: str):
    """Return ambient-sweep x/y for IDAES PLR+HX run outputs."""
    file_map = {
        "R134a": ROOT / "cop_vs_ambient_plr_hx_cold_storage_r134a_copy.csv",
        "R1234ze(E)": ROOT / "cop_vs_ambient_plr_hx_cold_storage_r1234zee_copy.csv",
    }
    rows = _read_csv(file_map[fluid_key])
    x = [_to_float(r["ambient_C"]) for r in rows]
    y = [_to_float(r["cop_part_load"]) for r in rows]
    return x, y


def _plot_for_fluid(fluid_key: str, out_name: str):
    prev_x, prev_cop, prev_carnot = _previous_series(fluid_key)

    x_hx, y_hx = _plr_hx_series(fluid_key)

    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=170)
    ax.plot(prev_x, prev_carnot, "k--", lw=2.2, label="carnot cop")
    ax.plot(prev_x, prev_cop, color="#4c78a8", marker="o", lw=2.2, label="PLR cop")
    ax.plot(x_hx, y_hx, color="#f58518", marker="^", lw=2.0, label="PLR + HX COP")

    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("COP")
    ax.set_title(f"{fluid_key}: COP overlay vs previous Carnot/COP")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()

    out_png = ROOT / f"{out_name}.png"
    out_pdf = ROOT / f"{out_name}.pdf"
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


def main():
    _plot_for_fluid("R134a", "overlay_carnot_vs_cop_r134a")
    _plot_for_fluid("R1234ze(E)", "overlay_carnot_vs_cop_r1234zee")


if __name__ == "__main__":
    main()
