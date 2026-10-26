"""NIST at Tamb=10 failed with the algebraic seed. Three restarts tried."""
import logging, warnings, io, contextlib
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use('Agg')
from pyomo.environ import value
import compare_cp_methods_0929_Helmholtz_SPGPlabel as M
from vapor_compression_cubic_refstate_1001 import SimpleVaporCompressionCycle as C, Mode

T_EVAP, APPROACH, ETA = -29.0, 9.0, 0.9999
m = {x.name: x for x in M.methods}['NIST']

def warmstart(model, pts):
    for s, u in (("evap_out", model.fs.evaporator), ("comp_out", model.fs.compressor),
                 ("cond_out", model.fs.condenser), ("valve_out", model.fs.expansion_valve)):
        st = u.control_volume.properties_out[0]
        st.temperature.set_value(pts[s]["T"]); st.pressure.set_value(pts[s]["P"])

def grab(model):
    out = {}
    for s, u in (("evap_out", model.fs.evaporator), ("comp_out", model.fs.compressor),
                 ("cond_out", model.fs.condenser), ("valve_out", model.fs.expansion_valve)):
        st = u.control_volume.properties_out[0]
        out[s] = {"T": value(st.temperature), "P": value(st.pressure)}
    return out

def solve_T(f, lo):
    a, fa = lo, f(lo)
    for hi in np.arange(lo+2.0, lo+600.0, 2.0):
        try: fb = f(hi)
        except Exception: continue
        if np.isfinite(fb) and fa*fb < 0: return brentq(f, a, hi)
        a, fa = hi, fb
    raise RuntimeError('no bracket')

def alg_points(Tev_C, Tcd_C):
    Te, Tc = Tev_C+273.15, Tcd_C+273.15
    pe, pc = m.saturation_pressure(Te), m.saturation_pressure(Tc)
    s1, h1 = m.sg(Te), m.hg(Te)
    T2s = solve_T(lambda T: m._s_mass(T, pc, 'vapor')+m.s_off-s1, Tc+0.01)
    h2 = h1 + ((m._h_mass(T2s, pc, 'vapor')+m.h_off)-h1)/ETA
    T2 = solve_T(lambda T: m._h_mass(T, pc, 'vapor')+m.h_off-h2, Tc+0.01)
    return {"evap_out":{"T":Te,"P":pe}, "comp_out":{"T":T2,"P":pc},
            "cond_out":{"T":Tc,"P":pc}, "valve_out":{"T":Te,"P":pe}}

def run(Tamb, seed):
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc = C('R32', compressor_efficiency=ETA, mode=Mode.IMPROVED_TPX, method='NIST')
            vc.specify_initial_conditions(low_side_temperature=T_EVAP,
                                          high_side_temperature=Tamb+APPROACH)
            vc.initialize(verbose=False)
            if seed is not None: warmstart(vc.model, seed)
            vc.set_specifications(ambient_temperature=Tamb, condenser_approach=APPROACH,
                                  evap_sat_temperature=T_EVAP, superheating=0,
                                  subcooling=0, max_pressure_ratio=10)
            cop, conv = vc.optimize_COP(verbose=False, initialize=True, optimize=False)
            pts = grab(vc.model)
        return cop, conv, pts
    except Exception as e:
        return None, False, None

print("NIST at T_amb = 10 C -- restart strategies\n")

a = alg_points(T_EVAP, 19.0)
c, v, _ = run(10, a)
print(f"  1. algebraic seed (repeat)        COP={c if c is None else round(c,4)}  converged={v}")

c15, v15, p15 = run(15, alg_points(T_EVAP, 24.0))
print(f"     [15 C seeded solve            COP={round(c15,4)}  converged={v15}]")
if v15:
    c, v, _ = run(10, p15)
    print(f"  2. seed from converged 15 C       COP={c if c is None else round(c,4)}  converged={v}")

c, v, p0 = run(10, None)
print(f"  3. no seed (cold start)           COP={c if c is None else round(c,4)}  converged={v}")
if v:
    c2, v2, _ = run(10, p0)
    print(f"  4. re-solve seeded from its own   COP={c2 if c2 is None else round(c2,4)}  converged={v2}")
