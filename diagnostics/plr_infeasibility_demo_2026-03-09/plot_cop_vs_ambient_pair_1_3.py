import csv
import math
import numpy as np
import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt

csv_path = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_vs_repo_relaxed_scsh_pairs_r134a.csv'
rows = list(csv.DictReader(open(csv_path)))

pair_rows = [
    r for r in rows
    if float(r['subcooling_C']) == 1.0 and float(r['superheating_C']) == 3.0
]

strict = sorted([r for r in pair_rows if r['band'] == 'Strict'], key=lambda r: float(r['ambient_C']))
relaxed = sorted([r for r in pair_rows if r['band'] == 'Relaxed(repo)'], key=lambda r: float(r['ambient_C']))

x_s = np.array([float(r['ambient_C']) for r in strict], dtype=float)
y_s = np.array([float(r['best_cop']) if math.isfinite(float(r['best_cop'])) else np.nan for r in strict], dtype=float)

x_r = np.array([float(r['ambient_C']) for r in relaxed], dtype=float)
y_r = np.array([float(r['best_cop']) if math.isfinite(float(r['best_cop'])) else np.nan for r in relaxed], dtype=float)

fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=180)
ax.plot(x_s, y_s, '-o', color='#d62728', linewidth=2.2, markersize=5, label='Strict')
ax.plot(x_r, y_r, '-s', color='#1f77b4', linewidth=2.2, markersize=5, label='Relaxed(repo)')

ax.set_xlabel('Ambient temperature (C)')
ax.set_ylabel('Best feasible COP')
ax.set_title('COP vs Ambient for (subcooling, superheating) = (1,3)')
ax.grid(True, alpha=0.25)
ax.legend(loc='best')

fig.tight_layout()
out_png = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_cop_vs_ambient_pair_1_3_r134a.png'
out_pdf = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_cop_vs_ambient_pair_1_3_r134a.pdf'
fig.savefig(out_png, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')
print(out_png)
print(out_pdf)
