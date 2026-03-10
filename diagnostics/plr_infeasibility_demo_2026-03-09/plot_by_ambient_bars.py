import csv
import numpy as np
import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt

csv_path = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_vs_repo_relaxed_scsh_pairs_r134a.csv'
rows = list(csv.DictReader(open(csv_path)))

pairs = [
    (1.0,1.0),(1.0,3.0),(1.0,5.0),(1.0,7.0),
    (3.0,3.0),(3.0,5.0),(3.0,7.0),(5.0,5.0),(5.0,7.0),
]
ambs = sorted({float(r['ambient_C']) for r in rows})

fig, axes = plt.subplots(3, 3, figsize=(14, 10), dpi=170, sharex=True, sharey=True)
axes = axes.flatten()

x = np.arange(len(ambs))
w = 0.38

for ax, (sc, sh) in zip(axes, pairs):
    strict = [r for r in rows if r['band']=='Strict' and float(r['subcooling_C'])==sc and float(r['superheating_C'])==sh]
    relax = [r for r in rows if r['band']=='Relaxed(repo)' and float(r['subcooling_C'])==sc and float(r['superheating_C'])==sh]

    y_strict = np.array([float(r['feasible_fraction']) for r in strict], dtype=float)
    y_relax = np.array([float(r['feasible_fraction']) for r in relax], dtype=float)

    ax.bar(x - w/2, y_strict, width=w, color='#d62728', label='Strict')
    ax.bar(x + w/2, y_relax, width=w, color='#1f77b4', label='Relaxed(repo)')
    ax.set_title(f'(sc,sh)=({int(sc)},{int(sh)})')
    ax.set_ylim(0.0, 1.05)
    ax.grid(axis='y', alpha=0.25)

for ax in axes:
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(a)) for a in ambs], rotation=0)

for i in [0,3,6]:
    axes[i].set_ylabel('Feasible fraction')
for i in [6,7,8]:
    axes[i].set_xlabel('Ambient temperature (C)')

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
fig.suptitle('Feasibility by Ambient: Strict vs Relaxed(repo), R134a, Tsp=-20C', y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])

out_png = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_feasibility_by_ambient_bars_r134a.png'
out_pdf = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_feasibility_by_ambient_bars_r134a.pdf'
fig.savefig(out_png, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')
print(out_png)
print(out_pdf)
