"""
Simple PFD Figure Generator (Exact Topology)

Author: Shilpa Narasimhan
Codex Support: OpenAI Codex
QA/Testing: Shilpa Narasimhan
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


OUTDIR = Path("/Users/snarasi2/idaes-hvacr-cycles/diagnostics/pfd")


def box(ax, x, y, w, h, label):
    ax.add_patch(Rectangle((x, y), w, h, fill=False, linewidth=1.6, edgecolor="black"))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=10)


def arr(ax, x1, y1, x2, y2, label=None):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="->", mutation_scale=12, linewidth=1.5))
    if label:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.06, label, ha="center", va="bottom", fontsize=8)


def style(ax, title):
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.axis("off")
    ax.text(0.2, 6.6, title, fontsize=12, fontweight="bold")


def draw_lumped(ax):
    style(ax, "Lumped Surrogate (Pre-ISO Style)")

    # Refrigerant train
    box(ax, 1.0, 4.5, 1.5, 0.9, "EVAP")
    box(ax, 3.3, 4.5, 1.5, 0.9, "COMP")
    box(ax, 5.6, 4.5, 1.8, 0.9, "COND")
    box(ax, 8.2, 4.5, 1.5, 0.9, "VALVE")

    arr(ax, 2.5, 4.95, 3.3, 4.95, "S1")
    arr(ax, 4.8, 4.95, 5.6, 4.95, "S2")
    arr(ax, 7.4, 4.95, 8.2, 4.95, "S3")
    arr(ax, 9.0, 4.5, 9.0, 3.5)
    arr(ax, 9.0, 3.5, 1.8, 3.5, "S4")
    arr(ax, 1.8, 3.5, 1.8, 4.5)

    # Air paths
    arr(ax, 4.9, 3.9, 5.6, 4.5, "Ambient air in")
    arr(ax, 7.0, 4.5, 7.7, 3.9, "Air out")
    arr(ax, 0.2, 5.5, 1.0, 4.95, "Cold-space air in")
    arr(ax, 2.5, 4.95, 3.1, 5.5, "Cold-space air out")


def draw_zoned(ax):
    style(ax, "IDAES Zoned Condenser (Pre-ISO Style)")

    # Refrigerant train
    box(ax, 0.8, 4.5, 1.3, 0.9, "EVAP")
    box(ax, 2.5, 4.5, 1.3, 0.9, "COMP")
    box(ax, 4.2, 4.5, 1.1, 0.9, "DS")
    box(ax, 5.8, 4.5, 1.3, 0.9, "COND")
    box(ax, 7.6, 4.5, 1.1, 0.9, "SC")
    box(ax, 9.2, 4.5, 1.2, 0.9, "VALVE")

    arr(ax, 2.1, 4.95, 2.5, 4.95, "S1")
    arr(ax, 3.8, 4.95, 4.2, 4.95, "S2")
    arr(ax, 5.3, 4.95, 5.8, 4.95, "S2a")
    arr(ax, 7.1, 4.95, 7.6, 4.95, "S3")
    arr(ax, 8.7, 4.95, 9.2, 4.95, "S3a")
    arr(ax, 9.8, 4.5, 9.8, 3.5)
    arr(ax, 9.8, 3.5, 1.5, 3.5, "S4")
    arr(ax, 1.5, 3.5, 1.5, 4.5)

    # Condenser air chain: SC -> COND -> DS
    arr(ax, 6.9, 3.9, 7.6, 4.5, "Ambient air in")
    arr(ax, 8.2, 4.5, 6.5, 3.9, "to COND")
    arr(ax, 6.5, 3.9, 4.8, 4.5, "to DS")
    arr(ax, 4.7, 4.5, 4.1, 3.9, "Exhaust")

    # Evaporator air path
    arr(ax, 0.1, 5.5, 0.8, 4.95, "Cold-space air in")
    arr(ax, 2.1, 4.95, 2.7, 5.5, "Cold-space air out")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    p1 = OUTDIR / "pfd_old_lumped.png"
    p2 = OUTDIR / "pfd_old_idaes_zoned.png"
    p3 = OUTDIR / "pfd_old_side_by_side.png"

    f1, a1 = plt.subplots(figsize=(12, 5), dpi=220)
    draw_lumped(a1)
    f1.tight_layout()
    f1.savefig(p1, bbox_inches="tight")
    plt.close(f1)

    f2, a2 = plt.subplots(figsize=(12, 5), dpi=220)
    draw_zoned(a2)
    f2.tight_layout()
    f2.savefig(p2, bbox_inches="tight")
    plt.close(f2)

    f3, (a3, a4) = plt.subplots(2, 1, figsize=(14, 10), dpi=220)
    draw_lumped(a3)
    draw_zoned(a4)
    f3.tight_layout()
    f3.savefig(p3, bbox_inches="tight")
    plt.close(f3)

    print(p1)
    print(p2)
    print(p3)


if __name__ == "__main__":
    main()

