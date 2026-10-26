"""
plot_cop_fullload_only_1015.py -- R-32 COP vs ambient, FULL LOAD ONLY.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/22/2026

What this does
    The left-hand panel of plot_cop_vs_helmholtz_carnot_1015.py on its own:
    Carnot, Helmholtz, PR/NIST, PR/GCGP, PR/SPGP at full load. Same data,
    same styling, same sources -- the part-load and deviation panels dropped.
Why
    Shilpa, 10/22: only the left-hand figure is wanted.
Basis
    T_evap -29 C, T_cond = T_amb + 9, SH/SC optimised in 3..8 K, eta 0.80.
Inputs
    sweep_eq_PR_1015.csv              PR, equation-oriented cycle
    sweep_helmholtz_maximized_1015.csv  Helmholtz, IDAES, maximised
Outputs
    cop_vs_ambient_fullload_1015.png  (and .pdf)
Run
    python plot_cop_fullload_only_1015.py
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AMB = [10.0, 15.0, 20.0, 25.0]

# Best COP over the square run and both optimiser starts -- SPGP's optimum is
# not always at SH 3, so taking the max is what makes the comparison fair.
pr = {}
for r in csv.DictReader(open("sweep_eq_PR_1015.csv")):
    if r["term"] != "optimal":
        continue
    k = (r["method"], float(r["T_amb"]))
    pr[k] = max(pr.get(k, 0.0), float(r["cop"]))

helm, carnot = {}, {}
for r in csv.DictReader(open("sweep_helmholtz_maximized_1015.csv")):
    a = float(r["ambient_C"])
    helm[a] = float(r["cop"])
    carnot[a] = float(r["carnot"])

INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6a64", "#e6e5e0", "#fcfcfb"
series = [  # label, values, colour, marker, linestyle, width
    ("Carnot", carnot, MUTED, None, (0, (5, 3)), 2.0),
    ("Helmholtz", helm, INK, "o", "-", 2.2),
    ("PR / NIST", {a: pr[("NIST", a)] for a in AMB}, "#2a78d6", "s", "-", 2.0),
    ("PR / GCGP", {a: pr[("GCGP", a)] for a in AMB}, "#eb6834", "^", "-", 2.0),
    ("PR / SPGP", {a: pr[("SPGP", a)] for a in AMB}, "#1baf7a", "D", "-", 2.0),
]

plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED})
fig, ax = plt.subplots(figsize=(7.2, 5.6), facecolor=SURF)
ax.set_facecolor(SURF)

handles = []
for lab, vals, col, mk, ls, lw in series:
    y = [vals[a] for a in AMB]
    z = 5 if lab == "Helmholtz" else 3
    h, = ax.plot(AMB, y, color=col, lw=lw, ls=ls, marker=mk,
                 ms=10 if lab == "Helmholtz" else 8,
                 markerfacecolor=(SURF if lab == "Helmholtz" else col),
                 markeredgecolor=(col if lab == "Helmholtz" else SURF),
                 markeredgewidth=1.8, label=lab, zorder=z)
    handles.append(h)

ax.set_title("Full load", color=INK, loc="left", fontsize=12)
ax.set_ylim(1.5, 5.5)
ax.set_ylabel("COP [–]")
ax.set_xlabel("Ambient temperature [°C]")
ax.set_xticks(AMB)
ax.set_xlim(8.5, 26.5)
ax.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
           fontsize=10.5, bbox_to_anchor=(0.5, 0.025))
fig.suptitle("R-32 COP vs ambient  ·  T_evap −29 °C, T_cond = T_amb + 9 K,\n"
             "η = 0.80, SH/SC optimised in 3–8 K",
             x=0.01, y=0.985, va="top", ha="left", color=INK, fontsize=12)
fig.text(0.01, 0.005,
         "PR: equation-oriented cycle on IDAES PR expressions "
         "(cycle_eq_PR_1015.py), best of square + two optimiser starts.\n"
         "Helmholtz: IDAES, maximised. SPGP = fp_single_0919; "
         "GCGP = scaled Shomate.",
         fontsize=8.0, color=MUTED)
fig.tight_layout(rect=(0, 0.13, 1, 0.945))
fig.savefig("cop_vs_ambient_fullload_1015.png", dpi=200, facecolor=SURF)
fig.savefig("cop_vs_ambient_fullload_1015.pdf", facecolor=SURF)

print(f"{'T_amb':>6}{'Carnot':>9}{'Helm':>9}{'NIST':>9}{'GCGP':>9}{'SPGP':>9}")
for a in AMB:
    print(f"{a:6.0f}{carnot[a]:9.4f}{helm[a]:9.4f}"
          + "".join(f"{pr[(m, a)]:9.4f}" for m in ("NIST", "GCGP", "SPGP")))
print("\nwrote cop_vs_ambient_fullload_1015.png and .pdf")
