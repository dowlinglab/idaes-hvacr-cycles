"""
cop_vs_ambient_seeded_1008.py -- R-32 COP vs ambient, IDAES solved with and
without algebraic warmstart seeding, for SPGP (fp_single_0919).

SOURCE OF THE NUMBERS
    Unseeded: phase6_combined_report.csv (phase6_final_0924.py, 2026-09-28),
              i.e. the poster values.
    Seeded:   algseed_test.py, DOWLINGLAB-7 / myidaesenv, 2026-10-08, as
              recorded in BREADCRUMB_09-28.md (10/08 "later" section).
              algseed_test.py prints but does not save; values transcribed.
    Helmholtz reference: phase6_combined_report.csv.

    NIST is NOT drawn: seeded NIST failed at 10 C and only the 15 C value
    (3.3428) was recorded. Re-run algseed_test.py to get 20/25 C.

    Closed-form SPGP is NOT drawn: its discharge temperatures are unphysical
    (33-52 C), see breadcrumb caveat. Do not quote its absolute values.

CYCLE: T_evap = -29 C, T_cond = T_amb + 9, eta = 0.9999, no SH/SC.
Formatting follows cop_vs_ambient_poster_1008.py (6.4 x 5.0, no title).

Author: Shilpa Narasimhan and Claude AI
Date Created: 10/08/2026
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AMB = np.array([10, 15, 20, 25], float)
HELM      = [3.948947, 3.526194, 3.190784, 2.914391]
SPGP_UNS  = [3.739092, 3.549206, 3.184143, 3.041566]
SPGP_SEED = [3.4878, 3.2201, 2.9714, 2.8036]

def curv(y):
    return np.round(np.diff(np.diff(np.asarray(y))), 4)

print(f"{'T_amb':>6}{'Helmholtz':>12}{'SPGP unseeded':>16}{'SPGP seeded':>14}")
for i, a in enumerate(AMB):
    print(f"{a:>6.0f}{HELM[i]:>12.4f}{SPGP_UNS[i]:>16.4f}{SPGP_SEED[i]:>14.4f}")
print("curvature  Helmholtz", curv(HELM), " unseeded", curv(SPGP_UNS),
      " seeded", curv(SPGP_SEED))

fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.plot(AMB, HELM, "-", color="k", marker="o", ms=5, lw=2.4,
        label="Helmholtz EoS / CoolProp")
ax.plot(AMB, SPGP_UNS, "--", color="tab:red", marker="D", ms=4.5, lw=1.3,
        mfc="none", mew=1.1, label="PR EoS / SPGP, cold start + chase")
ax.plot(AMB, SPGP_SEED, "-", color="tab:red", marker="D", ms=4.5, lw=1.6,
        label="PR EoS / SPGP, algebraic seed")
ax.set_xlabel("Ambient temperature [deg C]")
ax.set_ylabel("COP")
ax.set_xticks(AMB)
allv = np.concatenate([HELM, SPGP_UNS, SPGP_SEED])
pad = 0.05 * (allv.max() - allv.min())
ax.set_ylim(allv.min() - pad, allv.max() + pad)
ax.set_xlim(AMB.min() - 1.0, AMB.max() + 1.0)
ax.grid(True, ls=":", alpha=0.35)
ax.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
fig.tight_layout()
fig.savefig("cop_vs_ambient_seeded_1008.png", dpi=200)
print("wrote cop_vs_ambient_seeded_1008.png")
