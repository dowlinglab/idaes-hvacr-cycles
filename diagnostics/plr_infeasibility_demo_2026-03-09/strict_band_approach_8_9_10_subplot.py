import csv
import math

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr import Mode, SimpleVaporCompressionCyclePLR

# Assumed pair from latest runs
SUBCOOLING_C = 5.0
SUPERHEATING_C = 7.0
AMBIENT_GRID = np.arange(10.0, 46.0, 5.0)
TEVAP_GRID = [-30.0, -29.0, -28.0]
APPROACH_LIST = [8.0, 9.0, 10.0]
OUT_DIR = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09'


def solve_point(ambient_c, tevap_sat_c, approach_c):
    cycle = SimpleVaporCompressionCyclePLR(
        fluid_name='R134a', compressor_efficiency=0.75, mode=Mode.PH
    )
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    kwargs = dict(
        low_side_pressure=(60.0, 120.0),
        high_side_pressure=(500.0, 1000.0),
        evaporator_temperature=(-55.0, -20.0),
        condenser_temperature=(15.0, 60.0),
        subcooling=SUBCOOLING_C,
        superheating=SUPERHEATING_C,
        max_pressure_ratio=20.0,
        plr=0.75,
        cd=0.13,
        ambient_temperature=float(ambient_c),
        evap_sat_temperature=float(tevap_sat_c),
        condenser_approach=float(approach_c),
    )

    try:
        cycle.set_specifications(**kwargs)
        _, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        cop = cycle.get_part_load_cop() if ok else float('nan')
        if ok and math.isfinite(cop):
            return True, float(cop)
    except Exception:
        pass

    try:
        retry = dict(kwargs)
        retry['debug_disable_arc_pressure_eq'] = True
        cycle.set_specifications(**retry)
        _, ok = cycle.optimize_COP(verbose=False, initialize=True, optimize=False)
        cop = cycle.get_part_load_cop() if ok else float('nan')
        if ok and math.isfinite(cop):
            return True, float(cop)
    except Exception:
        pass

    return False, float('nan')


def main():
    rows = []
    for approach in APPROACH_LIST:
        for ambient in AMBIENT_GRID:
            for tevap in TEVAP_GRID:
                ok, cop = solve_point(ambient, tevap, approach)
                rows.append(
                    dict(
                        approach_C=float(approach),
                        ambient_C=float(ambient),
                        tevap_sat_C=float(tevap),
                        converged=int(ok),
                        cop=float(cop) if math.isfinite(cop) else float('nan'),
                    )
                )

    out_csv = f'{OUT_DIR}/strict_band_approach_8_9_10_pair57_r134a.csv'
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['approach_C', 'ambient_C', 'tevap_sat_C', 'converged', 'cop'])
        w.writeheader()
        w.writerows(rows)

    colors = {-30.0: '#1f77b4', -29.0: '#ff7f0e', -28.0: '#2ca02c'}
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), dpi=180, sharex=True, sharey=True)

    for ax, approach in zip(axes, APPROACH_LIST):
        subset = [r for r in rows if abs(r['approach_C'] - approach) < 1e-9]
        for tevap in TEVAP_GRID:
            cur = [r for r in subset if abs(r['tevap_sat_C'] - tevap) < 1e-9]
            x = np.array([r['ambient_C'] for r in cur], dtype=float)
            y = np.array([r['cop'] if math.isfinite(r['cop']) else np.nan for r in cur], dtype=float)
            ax.plot(x, y, '-o', color=colors[tevap], linewidth=2.0, markersize=4.5, label=f'Tevap={tevap:.0f} C')
            bad = np.isnan(y)
            if np.any(bad):
                ax.plot(x[bad], np.full(np.sum(bad), 1.0), 'x', color=colors[tevap], alpha=0.8)

        conv = sum(r['converged'] for r in subset)
        total = len(subset)
        ax.set_title(f'Approach={approach:.0f} C  (conv {conv}/{total})')
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel('COP')
    for ax in axes:
        ax.set_xlabel('Ambient (C)')
    axes[0].legend(loc='best', fontsize=8)

    fig.suptitle('Strict Tevap band [-30,-28] C, R134a, (subcool,superheat)=(5,7)', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_png = f'{OUT_DIR}/strict_band_approach_8_9_10_pair57_r134a.png'
    out_pdf = f'{OUT_DIR}/strict_band_approach_8_9_10_pair57_r134a.pdf'
    fig.savefig(out_png, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')

    print(f'Saved: {out_csv}')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_pdf}')


if __name__ == '__main__':
    main()
