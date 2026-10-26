"""
plot_cop_fullload_paper_1015.py -- R-32 COP vs ambient, full load, paper style.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/22/2026

What this does
    The full-load comparison -- Carnot, Helmholtz, PR/NIST, PR/GCGP, PR/SPGP
    -- in a publication style: no title, legend top right, full box around
    the axes, LaTeX axis labels, ticks inward on all four sides.
Why
    Shilpa, 10/22.
Basis
    T_evap -29 C, T_cond = T_amb + 9 K, SH/SC optimised in 3..8 K, eta 0.80.
Inputs
    sweep_eq_PR_1015.csv                PR, equation-oriented cycle
    sweep_helmholtz_maximized_1015.csv  Helmholtz, IDAES, maximised
Outputs
    cop_vs_ambient_fullload_paper_1015.png / .pdf
Run
    python plot_cop_fullload_paper_1015.py
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AMB = [10.0, 15.0, 20.0, 25.0]

# Best COP over the square run and both optimiser starts. SPGP's optimum is
# not always at SH = 3 K, so taking the max is what makes the comparison fair.
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

# mathtext rather than usetex: renders the same labels without needing a
# LaTeX installation. Set text.usetex=True if a full TeX toolchain is present.
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
    "font.size": 11,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
})

series = [  # label, values, colour, marker, linestyle
    (r"Carnot",          carnot,                                 "0.35",     None, (0, (5, 3))),
    (r"Helmholtz",       helm,                                   "k",        "o",  "-"),
    (r"PR / NIST",       {a: pr[("NIST", a)] for a in AMB},      "#2a78d6",  "s",  "-"),
    (r"PR / GCGP",       {a: pr[("GCGP", a)] for a in AMB},      "#eb6834",  "^",  "-"),
    (r"PR / SPGP",       {a: pr[("SPGP", a)] for a in AMB},      "#1baf7a",  "D",  "-"),
]

fig, ax = plt.subplots(figsize=(6.0, 4.6))

for lab, vals, col, mk, ls in series:
    ax.plot(AMB, [vals[a] for a in AMB], color=col, lw=1.6, ls=ls, marker=mk,
            ms=6, markerfacecolor=("none" if lab == "Helmholtz" else col),
            markeredgecolor=col, markeredgewidth=1.2, label=lab,
            zorder=5 if lab == "Helmholtz" else 3)

ax.set_xlabel(r"$T_{\mathrm{amb}}$ [$^\circ$C]")
ax.set_ylabel(r"$\mathrm{COP}$ [--]")
ax.set_xticks(AMB)
ax.set_xlim(8.5, 26.5)
ax.set_ylim(1.5, 5.5)
ax.minorticks_on()
ax.tick_params(which="minor", top=True, right=True)

# full box
for sp in ax.spines.values():
    sp.set_visible(True)

ax.legend(loc="upper right", frameon=True, framealpha=1.0, edgecolor="0.3",
          fancybox=False, fontsize=9.5, handlelength=2.4, borderpad=0.6)

fig.tight_layout()
fig.savefig("cop_vs_ambient_fullload_paper_1015.png", dpi=300)
fig.savefig("cop_vs_ambient_fullload_paper_1015.pdf")

print(f"{'T_amb':>6}{'Carnot':>9}{'Helm':>9}{'NIST':>9}{'GCGP':>9}{'SPGP':>9}")
for a in AMB:
    print(f"{a:6.0f}{carnot[a]:9.4f}{helm[a]:9.4f}"
          + "".join(f"{pr[(m, a)]:9.4f}" for m in ("NIST", "GCGP", "SPGP")))
print("\nwrote cop_vs_ambient_fullload_paper_1015.png and .pdf")
