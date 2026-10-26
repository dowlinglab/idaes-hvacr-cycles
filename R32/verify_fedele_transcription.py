"""Independent re-transcription of Fedele et al. (2014) Table 2, diffed against
the FEDELE_SAT table compiled into compare_cp_methods_R1234zeZ_1008.py.
Typed fresh from the PDF page image, not copied from the module."""
import sys; sys.path.insert(0,'/home/claude')
from compare_cp_methods_R1234zeZ_1008 import FEDELE_SAT

ITC = [(283.15,0.103),(285.15,0.111),(287.15,0.120),(289.15,0.129),(291.15,0.139),
       (293.15,0.149),(295.15,0.160),(297.15,0.172),(299.15,0.184),(301.15,0.198),
       (303.15,0.211),(305.15,0.225),(307.15,0.241),(309.15,0.256),(311.15,0.273),
       (313.15,0.291),(315.15,0.309),(317.15,0.329),(319.15,0.349),(321.14,0.370),
       (323.15,0.392),(325.15,0.415),(327.15,0.439),(329.15,0.465),(331.15,0.491),
       (333.15,0.519),(335.15,0.548),(337.15,0.578),(339.15,0.609),(341.15,0.641),
       (343.15,0.675),(345.15,0.710),(347.16,0.746),(349.15,0.784),(351.15,0.823),
       (353.16,0.863)]
UNI = [(238.13,0.012),(243.10,0.016),(248.00,0.021),(252.90,0.027),(257.81,0.034),
       (263.02,0.044),(267.91,0.055),(272.83,0.068),(277.94,0.085),(282.78,0.102),
       (287.62,0.123),(292.62,0.148),(298.00,0.178),(303.18,0.212),(308.18,0.249),
       (313.17,0.291),(318.14,0.338),(323.14,0.392),(328.12,0.451),(333.11,0.518),
       (338.12,0.591),(343.07,0.672),(348.16,0.764),(353.12,0.861),(358.13,0.969),
       (363.14,1.086),(368.09,1.212),(372.61,1.336)]

ref = {("ITC-CNR",T):P for T,P in ITC} | {("UnivPM",T):P for T,P in UNI}
mod = {(lab,T):P for T,P,lab in FEDELE_SAT}

print(f"PDF re-transcription : {len(ITC)} ITC-CNR + {len(UNI)} UnivPM = {len(ref)}")
print(f"module FEDELE_SAT    : {len(mod)}")
bad = []
for k in sorted(ref.keys() | mod.keys()):
    a, b = ref.get(k), mod.get(k)
    if a is None or b is None or abs(a-b) > 1e-12:
        bad.append((k, a, b))
if bad:
    print(f"\n{len(bad)} MISMATCH(ES):")
    for (lab,T),a,b in bad:
        print(f"   {lab:<9}{T:>8.2f}  pdf={a}  module={b}")
    sys.exit(1)
print("\nALL 64 POINTS MATCH (temperature and pressure, exactly).")
print("Ranges: ITC-CNR %.2f-%.2f K | UnivPM %.2f-%.2f K"
      % (ITC[0][0], ITC[-1][0], UNI[0][0], UNI[-1][0]))
