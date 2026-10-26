"""nist_tolcheck_1008.py -- NIST: does a tighter IPOPT tolerance (1e-8) make cold,
T,P-seeded and full-state-seeded solves agree? Usage: python nist_tolcheck_1008.py 15
Prints COP and evap/cond outlet T minus saturation T (K) for each start.
Author: Shilpa Narasimhan and Claude AI, 10/08/2026"""
import logging, warnings, io, contextlib, sys
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
from pyomo.environ import value, Var
from vapor_compression_cubic_refstate_1001 import SimpleVaporCompressionCycle as C, Mode
UN=("evap_out","comp_out","cond_out","valve_out")
OPTS={'tol':1e-8,'constr_viol_tol':1e-8,'acceptable_tol':1e-8,'max_iter':3000}
def units(m): return dict(zip(UN,(m.fs.evaporator,m.fs.compressor,m.fs.condenser,m.fs.expansion_valve)))
def run(Tamb, mode=None, src=None):
    buf=io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc=C('R32',compressor_efficiency=0.9999,mode=Mode.IMPROVED_TPX,method='NIST')
            vc.specify_initial_conditions(low_side_temperature=-29,high_side_temperature=Tamb+9)
            vc.initialize(verbose=False)
            if mode=="tp":
                for s,u in units(vc.model).items():
                    st=u.control_volume.properties_out[0]; st.temperature.set_value(src[s][0]); st.pressure.set_value(src[s][1])
            vc.set_specifications(ambient_temperature=Tamb,condenser_approach=9,evap_sat_temperature=-29,superheating=0,subcooling=0,max_pressure_ratio=10)
            if mode=="full":
                for v in vc.model.component_data_objects(Var,descend_into=True):
                    if v.fixed: continue
                    s=src.find_component(v.name)
                    if s is not None and s.value is not None: v.set_value(s.value,skip_validation=True)
            cop,conv=vc.optimize_COP(verbose=False,initialize=True,optimize=False,solver_options=OPTS)
    except Exception as e:
        return None,f"{type(e).__name__}",None
    m=vc.model; ev=m.fs.evaporator.control_volume.properties_out[0]; cd=m.fs.condenser.control_volume.properties_out[0]
    sh=value(ev.temperature)-value(ev.temperature_bubble["Vap","Liq"]); sc=value(cd.temperature)-value(cd.temperature_bubble["Vap","Liq"])
    st={s:(value(u.control_volume.properties_out[0].temperature),value(u.control_volume.properties_out[0].pressure)) for s,u in units(m).items()}
    return cop,conv,(sh,sc,st,m)
for Tamb in [int(a) for a in sys.argv[1:]]:
    c0,v0,i0=run(Tamb); line=f"Tamb={Tamb} cold {c0} {v0}"
    if i0: line+=f" [evap {i0[0]:+.3f} cond {i0[1]:+.3f}]"
    if v0 is True:
        c1,v1,i1=run(Tamb,"tp",i0[2]); c2,v2,i2=run(Tamb,"full",i0[3])
        line+=f" | tp {c1} {v1}" + (f" [evap {i1[0]:+.3f} cond {i1[1]:+.3f}]" if i1 else "")
        line+=f" | full {c2} {v2}" + (f" [evap {i2[0]:+.3f} cond {i2[1]:+.3f}]" if i2 else "")
    print(line,flush=True)
