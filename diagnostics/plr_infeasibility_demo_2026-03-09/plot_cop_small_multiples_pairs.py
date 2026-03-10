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

fig, axes = plt.subplots(3, 3, figsize=(13, 9), dpi=180, sharex=True, sharey=True)
axes = axes.flatten()

for ax, (sc, sh) in zip(axes, pairs):
    cur_strict = [r for r in rows if r['band']=='Strict' and float(r['subcooling_C'])==sc and float(r['superheating_C'])==sh]
    cur_relax = [r for r in rows if r['band']=='Relaxed(repo)' and float(r['subcooling_C'])==sc and float(r['superheating_C'])==sh]

    x_s = np.array([float(r['ambient_C']) for r in cur_strict], dtype=float)
    y_s = np.array([float(r['best_cop']) if math.isfinite(float(r['best_cop'])) else np.nan for r in cur_strict], dtype=float)
    x_r = np.array([float(r['ambient_C']) for r in cur_relax], dtype=float)
    y_r = np.array([float(r['best_cop']) if math.isfinite(float(r['best_cop'])) else np.nan for r in cur_relax], dtype=float)

    ax.plot(x_s, y_s, '-o', color='#d62728', linewidth=1.8, markersize=4, label='Strict')
    ax.plot(x_r, y_r, '-s', color='#1f77b4', linewidth=1.8, markersize=4, label='Relaxed(repo)')

    # Mark infeasible points at baseline for visibility
    bad_s = np.isnan(y_s)
    bad_r = np.isnan(y_r)
    if np.any(bad_s):
        ax.plot(x_s[bad_s], np.full(np.sum(bad_s), 1.0), 'x', color='#d62728', alpha=0.8)
    if np.any(bad_r):
        ax.plot(x_r[bad_r], np.full(np.sum(bad_r), 1.05), 'x', color='#1f77b4', alpha=0.8)

    ax.set_title(f'(sc,sh)=({int(sc)},{int(sh)})')
    ax.grid(True, alpha=0.25)

for i in [0,3,6]:
    axes[i].set_ylabel('COP')
for i in [6,7,8]:
    axes[i].set_xlabel('Ambient (C)')

axes[0].legend(loc='best', fontsize=8)
fig.suptitle('COP vs Ambient by (subcooling, superheating): Strict vs Relaxed(repo)', y=0.995)
fig.tight_layout(rect=[0,0,1,0.97])

out_png = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_cop_small_multiples_pairs_r134a.png'
out_pdf = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_cop_small_multiples_pairs_r134a.pdf'
fig.savefig(out_png, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')
print(out_png)
print(out_pdf)
