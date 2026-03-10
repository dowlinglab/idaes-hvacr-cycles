import csv
import math

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr import Mode, SimpleVaporCompressionCyclePLR

OUT_DIR = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09'
AMBIENT_GRID = np.arange(10.0, 46.0, 5.0)
EVAP_TEMP_GRID = np.arange(-55.0, -19.0, 5.0)

SUBCOOLING_C = 5.0
SUPERHEATING_C = 7.0
APPROACH_C = 15.0


def solve_point(ambient_c, evap_temp_c):
    cycle = SimpleVaporCompressionCyclePLR(fluid_name='R134a', compressor_efficiency=0.75, mode=Mode.PH)
    cycle.specify_initial_conditions(low_side_temperature=-20, high_side_temperature=30)
    cycle.initialize(verbose=False)

    kwargs = dict(
        low_side_pressure=(60.0, 120.0),
        high_side_pressure=(500.0, 1000.0),
        evaporator_temperature=(float(evap_temp_c), float(evap_temp_c)),
        condenser_temperature=(15.0, 60.0),
        subcooling=SUBCOOLING_C,
        superheating=SUPERHEATING_C,
        max_pressure_ratio=20.0,
        plr=0.75,
        cd=0.13,
        ambient_temperature=float(ambient_c),
        condenser_approach=APPROACH_C,
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
    for a in AMBIENT_GRID:
        for te in EVAP_TEMP_GRID:
            ok, cop = solve_point(a, te)
            rows.append(dict(ambient_C=float(a), evap_temp_C=float(te), converged=int(ok), cop=cop))

    out_csv = f'{OUT_DIR}/approach15_evaptemp_m55_m20_pair57_r134a_fresh.csv'
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['ambient_C', 'evap_temp_C', 'converged', 'cop'])
        w.writeheader(); w.writerows(rows)

    # COP vs ambient: one curve per evap temperature
    fig, ax = plt.subplots(figsize=(9.0, 5.4), dpi=180)
    cmap = plt.get_cmap('tab10')
    for i, te in enumerate(EVAP_TEMP_GRID):
        cur = [r for r in rows if abs(r['evap_temp_C'] - float(te)) < 1e-9]
        cur = sorted(cur, key=lambda r: float(r['ambient_C']))
        x = np.array([float(r['ambient_C']) for r in cur], dtype=float)
        y = np.array([float(r['cop']) if int(r['converged']) == 1 and math.isfinite(float(r['cop'])) else np.nan for r in cur], dtype=float)
        ax.plot(x, y, '-o', color=cmap(i % 10), linewidth=1.8, markersize=4, label=f'Evap temp {int(te)} C')

    ax.set_xlabel('Ambient temperature (C)')
    ax.set_ylabel('COP')
    ax.set_title('COP vs Ambient from fresh sweep (approach=15 C, evap temp -55..-20 C)')
    ax.grid(True, alpha=0.25)
    ax.legend(loc='best', ncol=2, fontsize=8)
    fig.tight_layout()

    out_png = f'{OUT_DIR}/cop_vs_ambient_fresh_sweep_approach15_evaptemp_m55_m20_pair57_r134a.png'
    out_pdf = f'{OUT_DIR}/cop_vs_ambient_fresh_sweep_approach15_evaptemp_m55_m20_pair57_r134a.pdf'
    fig.savefig(out_png, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')

    conv = sum(r['converged'] for r in rows)
    print(f'converged={conv}/{len(rows)}')
    print(f'Saved: {out_csv}')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_pdf}')


if __name__ == '__main__':
    main()
