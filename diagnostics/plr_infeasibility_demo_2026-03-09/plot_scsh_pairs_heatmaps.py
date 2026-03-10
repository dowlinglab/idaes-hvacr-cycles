import csv
import math
import numpy as np
import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt

csv_path = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_relaxed_scsh_pairs_r134a.csv'
rows = list(csv.DictReader(open(csv_path)))

pairs = [
    (1.0,1.0),(1.0,3.0),(1.0,5.0),(1.0,7.0),
    (3.0,3.0),(3.0,5.0),(3.0,7.0),
    (5.0,5.0),(5.0,7.0),
]
ambs = sorted({float(r['ambient_C']) for r in rows})


def matrix(band, field):
    m = np.full((len(pairs), len(ambs)), np.nan)
    for i,p in enumerate(pairs):
        for j,a in enumerate(ambs):
            r = next(rr for rr in rows if rr['band']==band and float(rr['subcooling_C'])==p[0] and float(rr['superheating_C'])==p[1] and float(rr['ambient_C'])==a)
            v = float(r[field])
            if field == 'best_cop' and (not math.isfinite(v)):
                v = np.nan
            m[i,j] = v
    return m

strict_f = matrix('Strict', 'feasible_fraction')
relax_f = matrix('Relaxed', 'feasible_fraction')
strict_c = matrix('Strict', 'best_cop')
relax_c = matrix('Relaxed', 'best_cop')

fig, axes = plt.subplots(2,2, figsize=(12,8), dpi=170, sharex=True, sharey=True)

im00 = axes[0,0].imshow(strict_f, aspect='auto', vmin=0, vmax=1, cmap='RdYlGn')
axes[0,0].set_title('Strict: Feasible fraction')
im01 = axes[0,1].imshow(relax_f, aspect='auto', vmin=0, vmax=1, cmap='RdYlGn')
axes[0,1].set_title('Relaxed: Feasible fraction')

cop_min = np.nanmin(np.concatenate([strict_c.flatten(), relax_c.flatten()]))
cop_max = np.nanmax(np.concatenate([strict_c.flatten(), relax_c.flatten()]))
im10 = axes[1,0].imshow(strict_c, aspect='auto', vmin=cop_min, vmax=cop_max, cmap='viridis')
axes[1,0].set_title('Strict: Best COP (NaN masked)')
im11 = axes[1,1].imshow(relax_c, aspect='auto', vmin=cop_min, vmax=cop_max, cmap='viridis')
axes[1,1].set_title('Relaxed: Best COP (NaN masked)')

for ax in axes.flat:
    ax.set_xticks(range(len(ambs)))
    ax.set_xticklabels([str(int(a)) for a in ambs], rotation=0)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([f'({int(p[0])},{int(p[1])})' for p in pairs])

axes[1,0].set_xlabel('Ambient C')
axes[1,1].set_xlabel('Ambient C')
axes[0,0].set_ylabel('(subcooling, superheating) C')
axes[1,0].set_ylabel('(subcooling, superheating) C')

cbar1 = fig.colorbar(im00, ax=[axes[0,0], axes[0,1]], fraction=0.03, pad=0.02)
cbar1.set_label('Feasible fraction')
cbar2 = fig.colorbar(im10, ax=[axes[1,0], axes[1,1]], fraction=0.03, pad=0.02)
cbar2.set_label('Best COP')

fig.suptitle('R134a, Tsp=-20C: Requested sc/sh pairs show widespread infeasibility above ~20-25C ambient', y=0.995)
fig.tight_layout()
out_png = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_scsh_pairs_heatmaps_r134a.png'
out_pdf = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_scsh_pairs_heatmaps_r134a.pdf'
fig.savefig(out_png, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')
print(out_png)
print(out_pdf)
