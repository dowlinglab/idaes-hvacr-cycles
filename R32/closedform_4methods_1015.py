"""
closedform_4methods_1015.py -- R-32 COP vs ambient, closed form, four property methods.

Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/15/2026

What this does
    COP from four state points, no solver, for Helmholtz (CoolProp), PR/NIST,
    PR/GCGP_scaled, PR/SPGP (= fp_single_0919).
Why
    The IDAES PR blocks can land on impossible phase states (10/15 night).
    Closed form picks the phase by fugacity sign, so it cannot.
Basis
    PR states: PRProbe (IDAES's own PR expressions, METHODS from
    phase_1_cubic_eos_validation_refstate_1001.py, PR's own Psat).
    T_evap -29 C, T_cond = T_amb + 9, SH 3 K, SC 8 K, eta 0.80.
Outputs
    closedform_4methods_1015.csv, table on screen.
Run
    python closedform_4methods_1015.py 2>&1 | tee closedform_4methods_1015.log
"""
import csv, logging, warnings
logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
from CoolProp.CoolProp import PropsSI as CP
import vapor_compression_cubic_FPhx_1015 as V
from phase_1_cubic_eos_validation_refstate_1001 import METHODS

T_EVAP, APPROACH, SH, SC, ETA = -29.0, 9.0, 3.0, 8.0, 0.80
AMB = [10.0, 15.0, 20.0, 25.0]

def helm(Tc_C):
    Te, Tc = T_EVAP + 273.15, Tc_C + 273.15
    pe = CP("P","T",Te,"Q",1,"R32"); pc = CP("P","T",Tc,"Q",0,"R32")
    h1 = CP("H","T",Te+SH,"P",pe,"R32"); s1 = CP("S","T",Te+SH,"P",pe,"R32")
    h2 = h1 + (CP("H","P",pc,"S",s1,"R32") - h1)/ETA
    h3 = CP("H","T",Tc-SC,"P",pc,"R32")
    return (h1-h3)/(h2-h1), CP("T","P",pc,"H",h2,"R32")-273.15

rows = []
for a in AMB:
    Tcond = a + APPROACH
    carnot = (T_EVAP+273.15)/(Tcond - T_EVAP)
    cop, Td = helm(Tcond); rows.append(dict(method="Helmholtz", T_amb=a, cop=cop, carnot=carnot, T_dis=Td))
for label, key in [("NIST","NIST"), ("GCGP","GCGP_scaled"), ("SPGP","fp_single_0919")]:
    probe = V.PRProbe(key)
    for a in AMB:
        Tcond = a + APPROACH
        r = V.closed_form(probe, T_EVAP, Tcond, SH, SC, ETA)
        rows.append(dict(method=label, T_amb=a, cop=r["cop"],
                         carnot=(T_EVAP+273.15)/(Tcond-T_EVAP), T_dis=r["T2"]-273.15,
                         Pe_kPa=r["Pe"]/1e3, Pc_kPa=r["Pc"]/1e3))

# checks
for m in ("Helmholtz","NIST","GCGP","SPGP"):
    s = [r for r in rows if r["method"]==m]
    assert all(r["cop"] < r["carnot"] for r in s), f"{m}: above Carnot"
    assert all(s[i]["cop"] > s[i+1]["cop"] for i in range(len(s)-1)), f"{m}: not monotonic"

print(f"{'T_amb':>6}{'Carnot':>9}{'Helmholtz':>11}{'NIST':>9}{'GCGP':>9}{'SPGP':>9}")
for a in AMB:
    g = {r["method"]: r for r in rows if r["T_amb"]==a}
    print(f"{a:>6.0f}{g['NIST']['carnot']:>9.4f}" + "".join(f"{g[m]['cop']:>{11 if m=='Helmholtz' else 9}.4f}" for m in ("Helmholtz","NIST","GCGP","SPGP")))
print("\n% vs Helmholtz")
for a in AMB:
    g = {r["method"]: r for r in rows if r["T_amb"]==a}; h = g["Helmholtz"]["cop"]
    print(f"{a:>6.0f}" + "".join(f"{100*(g[m]['cop']-h)/h:>9.2f}" for m in ("NIST","GCGP","SPGP")))
print("\nPR Psat at T_cond [kPa]  (CoolProp in brackets)")
for a in AMB:
    g = {r["method"]: r for r in rows if r["T_amb"]==a}
    print(f"{a:>6.0f}  [{CP('P','T',a+APPROACH+273.15,'Q',0,'R32')/1e3:7.1f}]" + "".join(f"{g[m]['Pc_kPa']:>9.1f}" for m in ("NIST","GCGP","SPGP")))
print("\ndischarge T [C]")
for a in AMB:
    g = {r["method"]: r for r in rows if r["T_amb"]==a}
    print(f"{a:>6.0f}" + "".join(f"{g[m]['T_dis']:>9.1f}" for m in ("Helmholtz","NIST","GCGP","SPGP")))
print("\nchecks passed: all below Carnot, all monotonic")
with open("closedform_4methods_1015.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["method","T_amb","cop","carnot","T_dis","Pe_kPa","Pc_kPa"])
    w.writeheader(); w.writerows(rows)
