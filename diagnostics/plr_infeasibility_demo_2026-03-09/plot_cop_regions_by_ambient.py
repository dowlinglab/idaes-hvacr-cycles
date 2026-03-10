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
    (3.0,3.0),(3.0,5.0),(3.0,7.0),
    (5.0,5.0),(5.0,7.0),
]
ambs = sorted({float(r['ambient_C']) for r in rows})

all_cops = [float(r['best_cop']) for r in rows if math.isfinite(float(r['best_cop']))]
vmin, vmax = min(all_cops), max(all_cops)


def make_figure(band_name, out_png, out_pdf):
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), dpi=170, sharex=True, sharey=True)
    axes = axes.flatten()
    last = None

    for ax, amb in zip(axes, ambs):
        sub = [r for r in rows if r['band'] == band_name and float(r['ambient_C']) == amb]
        x = np.array([float(r['superheating_C']) for r in sub], dtype=float)
        y = np.array([float(r['subcooling_C']) for r in sub], dtype=float)
        cop = np.array([float(r['best_cop']) if math.isfinite(float(r['best_cop'])) else np.nan for r in sub], dtype=float)
        feasible = np.array([int(r['feasible']) for r in sub], dtype=int)

        ok = np.isfinite(cop)
        if np.any(ok):
            last = ax.scatter(x[ok], y[ok], c=cop[ok], cmap='viridis', vmin=vmin, vmax=vmax,
                              s=320, marker='s', edgecolor='black', linewidth=0.5)
        if np.any(~ok):
            ax.scatter(x[~ok], y[~ok], c='lightgray', s=320, marker='s', edgecolor='black', linewidth=0.5)
            ax.scatter(x[~ok], y[~ok], c='black', s=45, marker='x')

        for xs, ys, c, fe in zip(x, y, cop, feasible):
            txt = f"{c:.2f}" if math.isfinite(c) else f"NA"
            ax.text(xs, ys, txt, ha='center', va='center', fontsize=7, color='white' if math.isfinite(c) else 'black')

        ax.set_title(f"Tamb={int(amb)}C")
        ax.set_xticks([1,3,5,7])
        ax.set_yticks([1,3,5,7])
        ax.grid(True, alpha=0.2)

    for i in [0,4]:
        axes[i].set_ylabel('Subcooling (C)')
    for i in [4,5,6,7]:
        axes[i].set_xlabel('Superheating (C)')

    if last is not None:
        cbar = fig.colorbar(last, ax=axes, fraction=0.02, pad=0.01)
        cbar.set_label('Best feasible COP')

    fig.suptitle(f'{band_name}: feasible COP regions by ambient (gray+X = infeasible)', y=0.995)
    fig.tight_layout(rect=[0,0,0.98,0.96])
    fig.savefig(out_png, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')


make_figure(
    'Strict',
    '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_cop_regions_by_ambient_r134a.png',
    '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_strict_cop_regions_by_ambient_r134a.pdf',
)
make_figure(
    'Relaxed(repo)',
    '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_relaxed_repo_cop_regions_by_ambient_r134a.png',
    '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/plr_relaxed_repo_cop_regions_by_ambient_r134a.pdf',
)

print('done')
