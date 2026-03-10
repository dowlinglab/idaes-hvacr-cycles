import csv
import math

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
import numpy as np

from vapor_compression_plr import Mode, SimpleVaporCompressionCyclePLR


def solve_point(ambient_c, tevap_sat_c):
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
        subcooling=3.0,
        superheating=3.0,
        max_pressure_ratio=20.0,
        plr=0.75,
        cd=0.13,
        ambient_temperature=float(ambient_c),
        evap_sat_temperature=float(tevap_sat_c),
        condenser_approach=10.0,
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
    ambient_grid = np.arange(10.0, 46.0, 5.0)
    tevap_grid = [-30.0, -29.0, -28.0]

    rows = []
    for a in ambient_grid:
        for te in tevap_grid:
            ok, cop = solve_point(a, te)
            rows.append(dict(ambient_C=float(a), tevap_sat_C=float(te), converged=int(ok), cop=cop))

    out_csv = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/strict_fixed10_pair33_r134a.csv'
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['ambient_C', 'tevap_sat_C', 'converged', 'cop'])
        w.writeheader()
        w.writerows(rows)

    fig, ax = plt.subplots(figsize=(8.3, 5.1), dpi=180)
    colors = {-30.0: '#1f77b4', -29.0: '#ff7f0e', -28.0: '#2ca02c'}
    for te in tevap_grid:
        cur = [r for r in rows if abs(r['tevap_sat_C'] - te) < 1e-9]
        x = np.array([r['ambient_C'] for r in cur], dtype=float)
        y = np.array([r['cop'] if math.isfinite(r['cop']) else np.nan for r in cur], dtype=float)
        ax.plot(x, y, '-o', color=colors[te], linewidth=2.0, markersize=4.8, label=f'Tevap_sat={te:.0f} C')

    ax.set_title('Strict band (3,3): COP vs Ambient\nTevap_sat in [-30,-28] C, approach=10 C')
    ax.set_xlabel('Ambient temperature (C)')
    ax.set_ylabel('COP')
    ax.grid(True, alpha=0.25)
    ax.legend(loc='best')
    fig.tight_layout()

    out_png = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/strict_fixed10_pair33_r134a.png'
    out_pdf = '/Users/snarasi2/idaes-hvacr-cycles/diagnostics/plr_infeasibility_demo_2026-03-09/strict_fixed10_pair33_r134a.pdf'
    fig.savefig(out_png, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')

    total = len(rows)
    conv = sum(r['converged'] for r in rows)
    print(f'converged={conv}/{total}')
    print(f'Saved: {out_csv}')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_pdf}')


if __name__ == '__main__':
    main()
