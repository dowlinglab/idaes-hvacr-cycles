"""
One focused test for 0.39 Btu/lb-R's vapor-side continuation: instead of
walking DOWN in pressure from its near-critical anchor (the direction every
previous attempt tried -- adaptive P-stepping, T-stepping, ideal-gas-biased
seeding -- all of which re-enter the two-phase dome because that direction
cools back toward the dew branch's entropy-hump peak), walk UP in pressure
from the anchor instead, toward and past the critical point.

Rationale
---------
0.39's anchor sits at T=378.4K, only ~3-3.5K below Tc (~381.5-381.9K). Above
Tc there is no two-phase region at all, for ANY entropy -- so if the walk can
get through that narrow near-critical gap without jumping onto the wrong
density root, everything above it is structurally guaranteed single-phase.
This is different from the earlier (gated-off) "upward walk" logic in
compute_isentrope_vapor_side, which was blocked for 0.39 based on comparing
its entropy target against the dome's global ceiling (s_v_max_dew) -- a test
that doesn't account for the anchor ALREADY being right at the critical
shoulder, as opposed to a mid-dome anchor being pushed upward from far away.

What this does
---------------
1. Reproduces 0.39's real anchor exactly as compute_isentrope_vapor_side
   finds it (same converged_dew search, same crit_point fallback logic).
2. Walks pressure UPWARD from that anchor toward ISENTROPE_P_MAX_PA (1350
   psia), using the same adaptive step-halving/doubling continuation
   already in the real code, reusing the same dome-reentry check (real
   bubble/dew enthalpy interpolation) at every step.
3. Reports every point's (T, P, h) and whether it lands inside the dome,
   plus how far the walk gets before either succeeding fully or getting
   stuck.

This is read-only / standalone -- it does NOT modify mixture_isentrope_validation.py.
Run from R515B_props_validated/:
    python3 diagnose_039_upward_from_anchor.py
"""

import numpy as np
from scipy.optimize import root

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    run_true_vle_envelope,
    mix_state,
    _mix_entropy_direct,
    _isentrope_2eq_residual,
    BTU_LBMR_TO_JKGK,
    ISENTROPE_P_MAX_PA,
)

FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
S_TARGET_BTU = 0.39
TMIN, TMAX, N = 255.4, 380.0, 80

d1 = load_idaes_helmholtz_json(FLUID1)
d2 = load_idaes_helmholtz_json(FLUID2)
mw1 = mw_from_json(d1)
mw2 = mw_from_json(d2)
z1 = w1_to_x1(W1, mw1, mw2)
mw_mix = z1 * mw1 + (1.0 - z1) * mw2
s_target = S_TARGET_BTU * BTU_LBMR_TO_JKGK * mw_mix
print(f"Target entropy for {S_TARGET_BTU} Btu/lb-R: {s_target:.3f} J/(mol*K)")

t_vals = np.linspace(TMIN, TMAX, N)
print(f"Running dome sweep (T={TMIN}-{TMAX}K, n={N})...")
bubble_rows, dew_rows, _, crit_point = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals)

converged_dew = [r for r in dew_rows if r["status"] == "CONVERGED"]
conv_bubble = sorted([r for r in bubble_rows if r["status"] == "CONVERGED"], key=lambda r: r["P_Pa"])
conv_dew_sorted = sorted(converged_dew, key=lambda r: r["P_Pa"])
bubble_P_arr = np.array([r["P_Pa"] for r in conv_bubble])
bubble_h_arr = np.array([mix_state(d1, d2, r["T_K"], r["rho_l_molm3"], z1).h_jmol for r in conv_bubble])
dew_P_arr = np.array([r["P_Pa"] for r in conv_dew_sorted])
dew_h_arr = np.array([mix_state(d1, d2, r["T_K"], r["rho_v_molm3"], z1).h_jmol for r in conv_dew_sorted])

def inside_dome(p_pa, h_jmol):
    if not (bubble_P_arr.min() <= p_pa <= bubble_P_arr.max()):
        return False
    if not (dew_P_arr.min() <= p_pa <= dew_P_arr.max()):
        return False
    h_bub = float(np.interp(p_pa, bubble_P_arr, bubble_h_arr))
    h_dew = float(np.interp(p_pa, dew_P_arr, dew_h_arr))
    return h_bub < h_jmol < h_dew

# --- reproduce the real anchor exactly as compute_isentrope_vapor_side finds it ---
best_row, best_gap = None, None
for r in converged_dew:
    t_k = r["T_K"]
    s_v = _mix_entropy_direct(d1, d2, t_k, r["rho_v_molm3"], z1)
    if s_v >= s_target:
        continue
    gap = s_target - s_v
    if best_gap is None or gap < best_gap:
        best_gap = gap
        best_row = r

if best_row is not None:
    t_k = best_row["T_K"]
    rho_seed = best_row["rho_v_molm3"]
    s_anchor = _mix_entropy_direct(d1, d2, t_k, rho_seed, z1)
    print(f"Anchor from dew-branch search: T={t_k:.2f}K, rho={rho_seed:.2f}, s_anchor={s_anchor:.3f} (gap={best_gap:.3f})")
else:
    t_crit, rho_crit = crit_point["T_K"], crit_point["rho_molm3"]
    t_k, rho_seed = t_crit, rho_crit
    print("Anchor from crit_point fallback")

# ramp entropy onto exact target at fixed pressure (same as real code)
p_seed = best_row["P_Pa"] if best_row is not None else crit_point["P_Pa"]
s_anchor = _mix_entropy_direct(d1, d2, t_k, rho_seed, z1)
n_substeps = max(1, int(np.ceil(abs(s_target - s_anchor) / 2.0)))
s_ramp = np.linspace(s_anchor, s_target, n_substeps + 1)[1:]
for s_step in s_ramp:
    sol = root(_isentrope_2eq_residual, x0=[t_k, rho_seed], args=(d1, d2, z1, p_seed, s_step), method="hybr", tol=1e-10)
    if not sol.success:
        print(f"Ramp failed at s_step={s_step:.3f}")
        break
    t_k, rho_seed = float(sol.x[0]), float(sol.x[1])

st0 = mix_state(d1, d2, t_k, rho_seed, z1)
p_start = float(st0.p_pa)
print(f"Seed after ramp: T={t_k:.4f}K, P={p_start/6894.76:.2f}psia, h={st0.h_jmol:.2f}, inside_dome={inside_dome(p_start, st0.h_jmol)}")
print(f"Critical point: T_c={crit_point['T_K']:.4f}K, P_c={crit_point['P_Pa']/6894.76:.2f}psia")

# --- walk UPWARD from the seed toward ISENTROPE_P_MAX_PA, adaptive step-halving ---
print("\n=== Walking UPWARD in pressure from the anchor (toward/past critical) ===")
log_from = np.log(p_start)
log_to = np.log(ISENTROPE_P_MAX_PA)
nominal_step = (log_to - log_from) / 78.0
min_step = nominal_step / 1024.0
step = nominal_step
t_cur, rho_cur = t_k, rho_seed
log_p_cur = log_from
n_ok, n_inside, n_fail = 0, 0, 0
while True:
    log_p_next = log_p_cur + step
    if log_p_next >= log_to:
        log_p_next = log_to
    if log_p_next == log_p_cur:
        break
    p_next = float(np.exp(log_p_next))
    sol = root(_isentrope_2eq_residual, x0=[t_cur, rho_cur], args=(d1, d2, z1, p_next, s_target), method="hybr", tol=1e-10)
    ok = bool(sol.success)
    if ok:
        t_new, rho_new = float(sol.x[0]), float(sol.x[1])
        st = mix_state(d1, d2, t_new, rho_new, z1)
        dome_hit = inside_dome(float(st.p_pa), float(st.h_jmol))
        print(f"  P={p_next/6894.76:7.2f}psia -> T={t_new:.3f}K rho={rho_new:.2f} h={st.h_jmol:.1f} inside_dome={dome_hit}")
        if dome_hit:
            n_inside += 1
            ok = False  # treat as failed step, same as real code's guard
        else:
            n_ok += 1
    if ok:
        t_cur, rho_cur = t_new, rho_new
        log_p_cur = log_p_next
        if log_p_cur == log_to:
            break
        grown = step * 2.0
        step = grown if abs(grown) <= abs(nominal_step) else abs(nominal_step)
    else:
        n_fail += 1
        step = step / 2.0
        if abs(step) < abs(min_step):
            print(f"  [stuck: step below min at log_p_cur -> P={np.exp(log_p_cur)/6894.76:.2f}psia, T={t_cur:.3f}K]")
            break

print(f"\nSummary: {n_ok} steps succeeded outside dome, {n_inside} rejected as inside dome, "
      f"{n_fail} total failed/rejected steps, final P reached = {np.exp(log_p_cur)/6894.76:.2f}psia "
      f"(target ceiling = {ISENTROPE_P_MAX_PA/6894.76:.2f}psia)")
