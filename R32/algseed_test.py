"""Does seeding IDAES from the algebraic state points reduce the slope noise?"""
import logging, warnings, io, contextlib, sys, json
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use('Agg')

import compare_cp_methods_0929_Helmholtz_SPGPlabel as M   # PRMethod objects
from vapor_compression_cubic_refstate_1001 import SimpleVaporCompressionCycle as C, Mode
# Inlined rather than imported: phase6_final_1001 has no __main__ guard, so
# importing it would run the entire sweep. This is a verbatim copy.
def warmstart_from_state(model, points):
    unit_map = {"evap_out": model.fs.evaporator,
                "comp_out": model.fs.compressor,
                "cond_out": model.fs.condenser,
                "valve_out": model.fs.expansion_valve}
    for stream, unit in unit_map.items():
        state = unit.control_volume.properties_out[0]
        state.temperature.set_value(points[stream]["T"])
        state.pressure.set_value(points[stream]["P"])

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.9999
AMB = [10, 15, 20, 25]
BY_NAME = {m.name: m for m in M.methods}

def solve_T(f, lo):
    a, fa = lo, f(lo)
    for hi in np.arange(lo+2.0, lo+600.0, 2.0):
        try: fb = f(hi)
        except Exception: continue
        if np.isfinite(fb) and fa*fb < 0: return brentq(f, a, hi)
        a, fa = hi, fb
    raise RuntimeError('no bracket')

def alg_points(m, Tev_C, Tcd_C):
    """Closed-form state points in the shape warmstart_from_state wants."""
    Te, Tc = Tev_C+273.15, Tcd_C+273.15
    pe, pc = m.saturation_pressure(Te), m.saturation_pressure(Tc)
    s1 = m.sg(Te); h1 = m.hg(Te)
    T2s = solve_T(lambda T: m._s_mass(T, pc, 'vapor')+m.s_off-s1, Tc+0.01)
    h2 = h1 + ((m._h_mass(T2s, pc, 'vapor')+m.h_off)-h1)/ETA
    T2 = solve_T(lambda T: m._h_mass(T, pc, 'vapor')+m.h_off-h2, Tc+0.01)
    return {"evap_out":  {"T": Te, "P": pe},
            "comp_out":  {"T": T2, "P": pc},
            "cond_out":  {"T": Tc, "P": pc},
            "valve_out": {"T": Te, "P": pe}}

def run(method, Tamb, seed):
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc = C('R32', compressor_efficiency=ETA, mode=Mode.IMPROVED_TPX, method=method)
            vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                          high_side_temperature=Tamb+APPROACH)
            vc.initialize(verbose=False)
            if seed is not None:
                warmstart_from_state(vc.model, seed)
            vc.set_specifications(ambient_temperature=Tamb, condenser_approach=APPROACH,
                                  evap_sat_temperature=T_EVAP, superheating=0,
                                  subcooling=0, max_pressure_ratio=10)
            cop, conv = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
        return cop, conv
    except Exception as e:
        return None, False

CFG = {"NIST": "NIST", "SPGP": "fp_single_0919"}
UNSEEDED = {"NIST":[3.567653,3.397422,3.131608,2.800516],
            "SPGP":[3.739092,3.549206,3.184143,3.041566]}

for disp, cfgname in CFG.items():
    m = BY_NAME[disp]
    seeded = []
    for Tamb in AMB:
        pts = alg_points(m, T_EVAP, Tamb+APPROACH)
        cop, conv = run(cfgname, Tamb, pts)
        seeded.append(cop if conv else np.nan)
        print(f'  {disp:<6} Tamb={Tamb:>3}  seeded COP = '
              f'{cop:.4f} conv={conv}' if cop is not None else
              f'  {disp:<6} Tamb={Tamb:>3}  FAILED', flush=True)
    u = np.array(UNSEEDED[disp]); s = np.array(seeded)
    print(f'\n  {disp}  unseeded : {np.round(u,4)}   curvature {np.round(np.diff(np.diff(u)),4)}')
    print(f'  {disp}  seeded   : {np.round(s,4)}   curvature {np.round(np.diff(np.diff(s)),4)}\n')
