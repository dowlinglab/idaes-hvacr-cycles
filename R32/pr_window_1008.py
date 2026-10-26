"""pr_window_1008.py -- NIST, GCGP, SPGP (fp_single_0919): optimize=True with superheat and
subcooling each free in a 3-8 K window (lower bound via the model's own constraints with
superheating=3, subcooling=3; 8 K ceilings added here). Do cold,
T,P-seeded and full-state-seeded solves agree? Usage: ONLY=NIST python pr_window_1008.py 15
Carnot limit at T_evap=-29 C, T_cond=T_amb+9 is printed for comparison.   (several ambients: 10 15 20 25)
Prints COP and evap/cond outlet T minus saturation T (K) for each start.
Author: Shilpa Narasimhan and Claude AI, 10/08/2026"""
import logging, warnings, io, contextlib, sys
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
from pyomo.environ import value, Var, Constraint
from vapor_compression_cubic_refstate_1001 import SimpleVaporCompressionCycle as C, Mode
UN=("evap_out","comp_out","cond_out","valve_out")
OPTS={'tol':1e-8,'constr_viol_tol':1e-8,'acceptable_tol':1e-8,'max_iter':3000}
def units(m): return dict(zip(UN,(m.fs.evaporator,m.fs.compressor,m.fs.condenser,m.fs.expansion_valve)))
def run(method, Tamb, mode=None, src=None):
    buf=io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            vc=C('R32',compressor_efficiency=0.9999,mode=Mode.IMPROVED_TPX,method=method)
            vc.specify_initial_conditions(low_side_temperature=-29,high_side_temperature=Tamb+9)
            vc.initialize(verbose=False)
            if mode=="tp":
                for s,u in units(vc.model).items():
                    st=u.control_volume.properties_out[0]; st.temperature.set_value(src[s][0]); st.pressure.set_value(src[s][1])
            vc.set_specifications(ambient_temperature=Tamb,condenser_approach=9,evap_sat_temperature=-29,superheating=3,subcooling=3,max_pressure_ratio=10)
            _ev=vc.model.fs.evaporator.control_volume.properties_out[0]; _cd=vc.model.fs.condenser.control_volume.properties_out[0]
            vc.model.fs.sh_max=Constraint(expr=_ev.temperature <= _ev.temperature_bubble["Vap","Liq"] + 8)
            vc.model.fs.sc_max=Constraint(expr=_cd.temperature >= _cd.temperature_bubble["Vap","Liq"] - 8)
            if mode=="full":
                for v in vc.model.component_data_objects(Var,descend_into=True):
                    if v.fixed: continue
                    s=src.find_component(v.name)
                    if s is not None and s.value is not None: v.set_value(s.value,skip_validation=True)
            cop,conv=vc.optimize_COP(verbose=False,initialize=True,optimize=True)
    except Exception as e:
        return None,f"{type(e).__name__}",None
    m=vc.model; ev=m.fs.evaporator.control_volume.properties_out[0]; cd=m.fs.condenser.control_volume.properties_out[0]
    sh=value(ev.temperature)-value(ev.temperature_bubble["Vap","Liq"]); sc=value(cd.temperature)-value(cd.temperature_bubble["Vap","Liq"])
    st={s:(value(u.control_volume.properties_out[0].temperature),value(u.control_volume.properties_out[0].pressure)) for s,u in units(m).items()}
    return cop,conv,(sh,sc,st,m)
ALL={"NIST":"NIST","GCGP":"GCGP","SPGP":"fp_single_0919"}
import os
METHODS={k:v for k,v in ALL.items() if os.environ.get("ONLY","") in ("",k)}
for lab,cfg in METHODS.items():
  for Tamb in [int(a) for a in sys.argv[1:]]:
    carnot=(273.15-29)/((Tamb+9)-(-29)); c0,v0,i0=run(cfg,Tamb); line=f"{lab:<5} Tamb={Tamb} Carnot {carnot:.3f} | cold {c0} {v0}"
    if i0: line+=f" [SH {i0[0]:+.3f} cond-Tsat {i0[1]:+.3f}]"
    if v0 is True:
        c1,v1,i1=run(cfg,Tamb,"tp",i0[2]); c2,v2,i2=run(cfg,Tamb,"full",i0[3])
        line+=f" | tp {c1} {v1}" + (f" [SH {i1[0]:+.3f} cond-Tsat {i1[1]:+.3f}]" if i1 else "")
        line+=f" | full {c2} {v2}" + (f" [SH {i2[0]:+.3f} cond-Tsat {i2[1]:+.3f}]" if i2 else "")
    print(line,flush=True)
