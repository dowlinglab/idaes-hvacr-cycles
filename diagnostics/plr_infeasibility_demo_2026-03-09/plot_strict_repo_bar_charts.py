import csv
import math
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
bands = ['Strict', 'Relaxed(repo)']
labels = [f'({int(a)},{int(b)})' for a,b in pairs]

feas_counts = {b: [] for b in bands}
mean_cop = {b: [] for b in bands}

for band in bands:
    for sc, sh in pairs:
        cur = [r for r in rows if r['band'] == band and float(r['subcooling_C']) == sc and float(r['superheating_C']) == sh]
        feasible_ambient_count = sum(1 for r in cur if int(r['feasible']) > 0)
        cvals = [float(r['best_cop']) for r in cur if math.isfinite(float(r['best_cop']))]
        feas_counts[band].append(feasible_ambient_count)
        mean_cop[band].append(float(np.mean(cvals)) if cvals else np.nan)

x = np.arange(len(pairs))
w = 0.38

fig, axes = plt.subplots(2, 1, figsize=(11, 8), dpi=170, sharex=True)

ax0 = axes[0]
ax0.bar(x - w/2, feas_counts['Strict'], width=w, label='Strict', color='#d62728')
ax0.bar(x + w/2, feas_counts['Relaxed(repo)'], width=w, label='Relaxed(repo)', color='#1f77b4')
ax0.set_ylabel('Feasible ambient points (out of 8)')
ax0.set_ylim(0, 8.5)
ax0.grid(axis='y', alpha=0.25)
ax0.legend(loc='best')
ax0.set_title('R134a: Strict vs Repo-Relaxed by (subcooling, superheating)')

ax1 = axes[1]
strict_vals = np.array(mean_cop['Strict'], dtype=float)
relax_vals = np.array(mean_cop['Relaxed(repo)'], dtype=float)
ax1.bar(x - w/2, np.nan_to_num(strict_vals, nan=0.0), width=w, label='Strict', color='#d62728')
ax1.bar(x + w/2, np.nan_to_num(relax_vals, nan=0.0), width=w, label='Relaxed(repo)', color='#1f77b4')
for i, v in enumerate(strict_vals):
    if np.isnan(v):
        ax1.text(i - w/2, 0.05, 'NA', ha='center', va='bottom', fontsize=8, rotation=90)
for i, v in enumerate(relax_vals):
    if np.isnan(v):
        ax1.text(i + w/2, 0.05, 'NA', ha='center', va='bottom', fontsize=8, rotation=90)
ax1.set_ylabel('Mean COP over feasible points')
ax1.set_xlabel('(subcooling, superheating) C')
ax1.grid(axis='y', alpha=0.25)

ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=0)

fig.tight_layout()
out_png = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_vs_repo_relaxed_bars_r134a.png'
out_pdf = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_vs_repo_relaxed_bars_r134a.pdf'
fig.savefig(out_png, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')

print(out_png)
print(out_pdf)
