import csv
import math
from collections import defaultdict

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
import numpy as np

csv_path = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_vs_repo_relaxed_scsh_pairs_r134a.csv'
rows = list(csv.DictReader(open(csv_path)))

pairs = [
    (1.0, 1.0), (1.0, 3.0), (1.0, 5.0), (1.0, 7.0),
    (3.0, 3.0), (3.0, 5.0), (3.0, 7.0), (5.0, 5.0), (5.0, 7.0),
]


def summarize(band_name):
    out = []
    for sc, sh in pairs:
        cur = [
            r for r in rows
            if r['band'] == band_name
            and float(r['subcooling_C']) == sc
            and float(r['superheating_C']) == sh
        ]
        cvals = [float(r['best_cop']) for r in cur if math.isfinite(float(r['best_cop']))]
        feasible_ambient = sum(1 for r in cur if int(r['feasible']) > 0)
        out.append({
            'sc': sc,
            'sh': sh,
            'mean_cop': float(np.mean(cvals)) if cvals else float('nan'),
            'n_feasible_ambient': feasible_ambient,
            'n_ambient_total': len(cur),
        })
    return out


strict = summarize('Strict')
relaxed = summarize('Relaxed(repo)')

all_cops = [r['mean_cop'] for r in strict + relaxed if math.isfinite(r['mean_cop'])]
vmin = min(all_cops)
vmax = max(all_cops)

fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), dpi=180, sharex=True, sharey=True)

for ax, data, title in [
    (axes[0], strict, 'Strict band'),
    (axes[1], relaxed, 'Relaxed(repo) band'),
]:
    x = np.array([d['sh'] for d in data], dtype=float)
    y = np.array([d['sc'] for d in data], dtype=float)
    c = np.array([d['mean_cop'] for d in data], dtype=float)
    n = np.array([d['n_feasible_ambient'] for d in data], dtype=float)

    # Marker size also carries robustness across ambient (bigger = feasible at more ambient points).
    sizes = 80 + 50 * n
    sca = ax.scatter(x, y, c=c, s=sizes, cmap='viridis', vmin=vmin, vmax=vmax, edgecolor='black', linewidth=0.4)

    for d in data:
        txt = f"{d['mean_cop']:.2f}\n({d['n_feasible_ambient']}/{d['n_ambient_total']})"
        ax.text(d['sh'], d['sc'], txt, ha='center', va='center', fontsize=7, color='white')

    ax.set_title(title)
    ax.set_xticks([1, 3, 5, 7])
    ax.set_yticks([1, 3, 5, 7])
    ax.grid(True, alpha=0.2)

axes[0].set_ylabel('Subcooling (C)')
for ax in axes:
    ax.set_xlabel('Superheating (C)')

cbar = fig.colorbar(sca, ax=axes, fraction=0.03, pad=0.04)
cbar.set_label('Mean COP over feasible ambient points')

fig.suptitle('COP sensitivity to superheating/subcooling (text: mean COP and feasible ambient count)', y=0.99)
fig.tight_layout(rect=[0, 0, 0.96, 0.95])

out_png = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_cop_sensitivity_scsh_r134a.png'
out_pdf = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_cop_sensitivity_scsh_r134a.pdf'
fig.savefig(out_png, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')

print(out_png)
print(out_pdf)
