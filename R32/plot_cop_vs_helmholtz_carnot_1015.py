"""
plot_cop_vs_helmholtz_carnot_1015.py -- R-32 COP vs ambient: Carnot, Helmholtz,
PR/NIST, PR/GCGP, PR/SPGP; full load and part load.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    Reads sweep_eq_PR_1015.csv (PR, equation-oriented cycle) and
    sweep_helmholtz_maximized_1015.csv (Helmholtz, IDAES, maximised).
    For each PR method and ambient takes the BEST COP over the square run
    and both optimiser starts (SPGP's optimum is not always SH 3).
    Panels: full load; part load, COP x PLF, PLF = 1 - CD(1 - PLR) = 0.9675
    (CD 0.13, PLR 0.75), Carnot unchanged; % deviation from Helmholtz (PLF
    cancels, so one panel serves both).
Why
    Shilpa, 10/15: plot all methods against Helmholtz and Carnot.
Basis
    T_evap -29 C, T_cond = T_amb + 9, SH/SC 3..8 K, eta 0.80, linear axes.
Outputs
    cop_vs_ambient_eq_PR_1015.png
Run
    python plot_cop_vs_helmholtz_carnot_1015.py
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLF = 1 - 0.13 * (1 - 0.75)
AMB = [10.0, 15.0, 20.0, 25.0]

pr = {}
for r in csv.DictReader(open("sweep_eq_PR_1015.csv")):
    if r["term"] != "optimal": continue
    k = (r["method"], float(r["T_amb"]))
    pr[k] = max(pr.get(k, 0.0), float(r["cop"]))
helm, carnot = {}, {}
for r in csv.DictReader(open("sweep_helmholtz_maximized_1015.csv")):
    a = float(r["ambient_C"]); helm[a] = float(r["cop"]); carnot[a] = float(r["carnot"])

INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6a64", "#e6e5e0", "#fcfcfb"
series = [  # label, values, color, marker, linestyle, width
    ("Carnot", carnot, MUTED, None, (0, (5, 3)), 2.0),
    ("Helmholtz", helm, INK, "o", "-", 2.2),
    ("PR / NIST", {a: pr[("NIST", a)] for a in AMB}, "#2a78d6", "s", "-", 2.0),
    ("PR / GCGP", {a: pr[("GCGP", a)] for a in AMB}, "#eb6834", "^", "-", 2.0),
    ("PR / SPGP", {a: pr[("SPGP", a)] for a in AMB}, "#1baf7a", "D", "-", 2.0),
]

plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), facecolor=SURF,
                         gridspec_kw={"width_ratios": [1, 1, 1]})
handles = []
for ax, (title, f) in zip(axes[:2], (("Full load", 1.0), (f"Part load (PLF = {PLF:.4f})", PLF))):
    ax.set_facecolor(SURF)
    for lab, vals, col, mk, ls, lw in series:
        y = [vals[a] * (1.0 if lab == "Carnot" else f) for a in AMB]
        z = 5 if lab == "Helmholtz" else 3
        h, = ax.plot(AMB, y, color=col, lw=lw, ls=ls, marker=mk, ms=8 if lab != "Helmholtz" else 10,
                     markerfacecolor=(SURF if lab == "Helmholtz" else col),
                     markeredgecolor=(col if lab == "Helmholtz" else SURF), markeredgewidth=1.8,
                     label=lab, zorder=z)
        if ax is axes[0]: handles.append(h)
    ax.set_title(title, color=INK, loc="left", fontsize=12)
    ax.set_ylim(1.5, 5.5); ax.set_ylabel("COP [–]")
# deviation from Helmholtz (PLF cancels, so one panel serves both loads)
ax = axes[2]; ax.set_facecolor(SURF)
ax.axhline(0, color=INK, lw=2.2, zorder=4)
for lab, vals, col, mk, ls, lw in series[2:]:
    ax.plot(AMB, [100 * (vals[a] - helm[a]) / helm[a] for a in AMB], color=col, lw=lw,
            marker=mk, ms=8, markeredgecolor=SURF, markeredgewidth=1.5, zorder=3)
ax.set_title("Deviation from Helmholtz (full and part load)", color=INK, loc="left", fontsize=12)
ax.set_ylabel("(COP − COP_Helmholtz) / COP_Helmholtz  [%]")
for ax in axes:
    ax.set_xlabel("Ambient temperature [°C]")
    ax.set_xticks(AMB); ax.set_xlim(8.5, 26.5)
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=10.5,
           bbox_to_anchor=(0.5, 0.06))
fig.suptitle("R-32 COP vs ambient  ·  T_evap −29 °C, T_cond = T_amb + 9 K, η = 0.80, "
             "SH/SC optimised in 3–8 K", x=0.01, ha="left", color=INK, fontsize=12)
fig.text(0.01, 0.01, "PR: equation-oriented cycle on IDAES PR expressions (cycle_eq_PR_1015.py), "
         "best of square + two optimiser starts. Helmholtz: IDAES, maximised "
         "(sweep_helmholtz_maximized_1015.csv). SPGP = fp_single_0919; GCGP = scaled Shomate.",
         fontsize=8.5, color=MUTED)
fig.tight_layout(rect=(0, 0.12, 1, 0.95))
fig.savefig("cop_vs_ambient_eq_PR_1015.png", dpi=200, facecolor=SURF)

print(f"{'T_amb':>6}{'Carnot':>9}{'Helm':>9}{'NIST':>9}{'GCGP':>9}{'SPGP':>9}   (full load)")
for a in AMB:
    print(f"{a:6.0f}{carnot[a]:9.4f}{helm[a]:9.4f}" + "".join(f"{pr[(m, a)]:9.4f}" for m in ("NIST", "GCGP", "SPGP")))
