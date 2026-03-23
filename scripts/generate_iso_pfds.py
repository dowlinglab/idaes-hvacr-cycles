"""
ISO-Style PFD Figure Generator (Lumped vs Zoned Refrigeration Models)

Author: Shilpa Narasimhan
Codex Support: OpenAI Codex
QA/Testing: Shilpa Narasimhan

Description:
    Generate ISO-style process flow diagrams (PFD-like) for two model
    structures:
    1) Lumped surrogate with a single condenser block
    2) Zoned IDAES condenser train (DS -> COND -> SC)
    Also generates a side-by-side comparison PNG.

Context Breadcrumb:
    This utility is for communication/diagnostics only and does not alter
    the simulation model equations. It provides publication-ready schematic
    views that align with the current brainstorming thread.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path("/Users/snarasi2/idaes-hvacr-cycles")
OUTDIR = ROOT / "diagnostics" / "pfd"


def _unit(ax, x, y, w, h, label):
    ax.add_patch(Rectangle((x, y), w, h, fill=False, linewidth=1.8, edgecolor="black"))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=9)


def _arrow(ax, x1, y1, x2, y2, label=None):
    arr = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="->", mutation_scale=12, linewidth=1.5, color="black")
    ax.add_patch(arr)
    if label:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.08, label, ha="center", va="bottom", fontsize=8)


def _format_ax(ax, title):
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.text(0.1, 5.7, title, fontsize=11, fontweight="bold", ha="left")
    ax.text(0.1, 5.45, "ISO-style PFD schematic (conceptual)", fontsize=8, ha="left")


def draw_lumped(ax):
    _format_ax(ax, "Model A: Lumped Surrogate (Single COND)")

    _unit(ax, 1.0, 3.9, 1.5, 0.9, "EVAP")
    _unit(ax, 3.2, 3.9, 1.5, 0.9, "COMP")
    _unit(ax, 5.4, 3.9, 1.7, 0.9, "COND")
    _unit(ax, 7.8, 3.9, 1.5, 0.9, "VALVE")

    _arrow(ax, 2.5, 4.35, 3.2, 4.35, "S1")
    _arrow(ax, 4.7, 4.35, 5.4, 4.35, "S2")
    _arrow(ax, 7.1, 4.35, 7.8, 4.35, "S3")
    _arrow(ax, 8.55, 3.9, 2.0, 3.9, "S4")
    _arrow(ax, 2.0, 3.9, 2.0, 3.9)  # closure visual

    # Air side annotations
    _arrow(ax, 5.0, 2.8, 5.4, 3.9, "Ambient air in")
    _arrow(ax, 6.3, 3.9, 6.7, 2.8, "Air out")
    _arrow(ax, 0.7, 4.8, 1.0, 4.35, "Cold-space air in")
    _arrow(ax, 2.5, 4.35, 2.8, 4.8, "Cold-space air out")


def draw_zoned(ax):
    _format_ax(ax, "Model B: IDAES Zoned Condenser (DS-COND-SC)")

    _unit(ax, 0.8, 3.9, 1.4, 0.9, "EVAP")
    _unit(ax, 2.7, 3.9, 1.4, 0.9, "COMP")
    _unit(ax, 4.6, 3.9, 1.2, 0.9, "DS")
    _unit(ax, 6.2, 3.9, 1.4, 0.9, "COND")
    _unit(ax, 8.0, 3.9, 1.2, 0.9, "SC")
    _unit(ax, 9.6, 3.9, 1.2, 0.9, "VALVE")

    _arrow(ax, 2.2, 4.35, 2.7, 4.35, "S1")
    _arrow(ax, 4.1, 4.35, 4.6, 4.35, "S2")
    _arrow(ax, 5.8, 4.35, 6.2, 4.35, "S2a")
    _arrow(ax, 7.6, 4.35, 8.0, 4.35, "S3")
    _arrow(ax, 9.2, 4.35, 9.6, 4.35, "S3a")
    _arrow(ax, 10.2, 3.9, 1.5, 3.9, "S4")

    # Condenser air chain SC -> COND -> DS
    _arrow(ax, 7.6, 2.7, 8.0, 3.9, "Ambient air in")
    _arrow(ax, 8.6, 3.9, 6.9, 3.3, "to COND")
    _arrow(ax, 6.9, 3.3, 5.2, 3.9, "to DS")
    _arrow(ax, 5.0, 3.9, 4.6, 2.7, "Exhaust")

    # Evaporator air loop
    _arrow(ax, 0.3, 4.8, 0.8, 4.35, "Cold-space air in")
    _arrow(ax, 2.2, 4.35, 2.5, 4.8, "Cold-space air out")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Individual figures
    fig1, ax1 = plt.subplots(figsize=(12, 5), dpi=220)
    draw_lumped(ax1)
    p1 = OUTDIR / "pfd_iso_lumped_surrogate.png"
    fig1.tight_layout()
    fig1.savefig(p1, bbox_inches="tight")
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(12, 5), dpi=220)
    draw_zoned(ax2)
    p2 = OUTDIR / "pfd_iso_idaes_zoned.png"
    fig2.tight_layout()
    fig2.savefig(p2, bbox_inches="tight")
    plt.close(fig2)

    # Side-by-side
    fig3, (ax3, ax4) = plt.subplots(2, 1, figsize=(14, 10), dpi=220)
    draw_lumped(ax3)
    draw_zoned(ax4)
    p3 = OUTDIR / "pfd_iso_side_by_side.png"
    fig3.tight_layout()
    fig3.savefig(p3, bbox_inches="tight")
    plt.close(fig3)

    print(str(p1))
    print(str(p2))
    print(str(p3))


if __name__ == "__main__":
    main()

