"""
Diagnostic script: trace the VAPOR-side isentrope walk for a single
target entropy (default 0.37 Btu/lb-R) step by step.

Purpose
-------
Mirror of diagnose_isentrope_026.py, but targeting compute_isentrope_vapor_side()
instead of compute_isentrope_liquid_side(). The liquid side (0.22-0.35 Btu/lb-R)
is confirmed working after the 2026-08-14 residual-normalization fix to
_isentrope_2eq_residual(); the vapor side (0.37 Btu/lb-R onward) is visually
broken in the regenerated p-H chart even with the same fixed residual function,
so this script isolates WHERE in the vapor-side anchor-selection / pressure-walk
it goes wrong (same class of check: does rho move monotonically -- vapor density
should DECREASE as pressure drops toward ISENTROPE_P_MIN_PA -- and does the
solver ever fail to converge).

Usage
-----
Run from R515B_props_validated/ (same directory as mixture_isentrope_validation.py
and its own copy of linear_model_codex.py):
    python3 diagnose_isentrope_vapor_037.py
Change S_TARGET_BTU below to inspect a different vapor-side isentrope.

Does not modify mixture_isentrope_validation.py or any of its outputs
-- read-only diagnostic, safe to run repeatedly.
"""

import numpy as np
from scipy.optimize import root

from mixture_isentrope_validation import (
    run_true_vle_envelope,
    load_idaes_helmholtz_json,
    mw_from_json,
    mix_state,
    _mix_entropy_direct,
    _isentrope_2eq_residual,
    BTU_LBMR_TO_JKGK,
    ISENTROPE_P_MIN_PA,
    ISENTROPE_N_POINTS,
)

# ---- Match this to whatever dome run you've been comparing against ----
FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
TMIN, TMAX, N = 255.4, 380.0, 80  # Tmax raised to ~380 K per the 2026-08-14 dome-kink fix
S_TARGET_BTU = 0.37  # <-- first confirmed-bad vapor-side isentrope; change to inspect others
# -------------------------------------------------------------------------

t_vals = np.linspace(TMIN, TMAX, N)
print(f"Running dome sweep ({FLUID1}/{FLUID2}, w1={W1}, T={TMIN}-{TMAX} K, n={N})...")
bubble_rows, dew_rows, z1, crit_point = run_true_vle_envelope(FLUID1, FLUID2, W1, t_vals)

d1 = load_idaes_helmholtz_json(FLUID1)
d2 = load_idaes_helmholtz_json(FLUID2)
mw1 = mw_from_json(d1)
mw2 = mw_from_json(d2)
mw_mix = z1 * mw1 + (1.0 - z1) * mw2

s_target = S_TARGET_BTU * BTU_LBMR_TO_JKGK * mw_mix
print(f"\ns_target = {S_TARGET_BTU} Btu/lb-R = {s_target:.4f} J/(mol*K)  (mw_mix={mw_mix*1000:.3f} g/mol)\n")

converged_dew = [r for r in dew_rows if r["status"] == "CONVERGED"]
print(f"{len(converged_dew)}/{len(dew_rows)} dew rows converged.\n")

# ---- Step 1: anchor selection (mirrors compute_isentrope_vapor_side) ----
print("=== Step 1: anchor candidates (rows where s_v(T) < s_target) ===")
print(f"{'T_K':>10} {'s_v':>12} {'gap':>12}  best?")
best_row = None
best_gap = None
for r in converged_dew:
    t_k = r["T_K"]
    s_v = _mix_entropy_direct(d1, d2, t_k, r["rho_v_molm3"], z1)
    if s_v >= s_target:
        continue
    gap = s_target - s_v
    is_best = best_gap is None or gap < best_gap
    if is_best:
        best_gap = gap
        best_row = r
    print(f"{t_k:10.3f} {s_v:12.4f} {gap:12.4f}  {'<-- BEST' if is_best else ''}")

if best_row is None:
    print("\nNo anchor row found -- this isentrope never reaches the superheated-vapor branch in this T range.")
    raise SystemExit(0)

t_k = best_row["T_K"]
rho_seed = best_row["rho_v_molm3"]
p_seed = best_row["P_Pa"]
print(f"\nAnchor selected: T={t_k:.4f} K   rho_seed={rho_seed:.4f} mol/m3   p_seed={p_seed:.2f} Pa")

# ---- Step 2: initial solve at P_target = p_seed ----
print("\n=== Step 2: initial solve (should be near-trivial since seed IS the answer) ===")
sol0 = root(_isentrope_2eq_residual, x0=[t_k, rho_seed], args=(d1, d2, z1, p_seed, s_target), method="hybr", tol=1.0e-10)
print(f"success={sol0.success}  T_solved={sol0.x[0]:.4f}  rho_solved={sol0.x[1]:.4f}")
print(f"  moved from seed by: dT={sol0.x[0]-t_k:+.4f} K   drho={sol0.x[1]-rho_seed:+.4f} mol/m3")
if not sol0.success:
    print(f"  sol0.message: {sol0.message}")
    print(f"  sol0.fun (residual at returned point): {sol0.fun}")
    print("INITIAL SOLVE FAILED -- stopping here.")
    raise SystemExit(0)

t_walk, rho_walk = float(sol0.x[0]), float(sol0.x[1])
st0 = mix_state(d1, d2, t_walk, rho_walk, z1)
s0 = _mix_entropy_direct(d1, d2, t_walk, rho_walk, z1)
print(f"  Achieved P={st0.p_pa:.2f} Pa (target {p_seed:.2f})   Achieved s={s0:.4f} (target {s_target:.4f})")

# ---- Step 3: the pressure walk itself (DOWN from saturation to ISENTROPE_P_MIN_PA) ----
p_start = float(st0.p_pa)
print(f"\n=== Step 3: pressure walk DOWN, {p_start:.1f} Pa -> {ISENTROPE_P_MIN_PA:.1f} Pa, {ISENTROPE_N_POINTS} log-spaced steps ===")
print(f"{'step':>4} {'P_target':>12} {'T_K':>10} {'rho':>12} {'dRho':>10} {'P_achieved':>12} {'s_achieved':>10} {'ok':>5}")

if p_start <= ISENTROPE_P_MIN_PA:
    print(f"p_start ({p_start:.1f}) already <= ISENTROPE_P_MIN_PA ({ISENTROPE_P_MIN_PA:.1f}) -- no walk needed/possible.")
else:
    prev_rho = rho_walk
    for i, p_target in enumerate(np.geomspace(p_start, ISENTROPE_P_MIN_PA, ISENTROPE_N_POINTS)[1:], start=1):
        sol = root(_isentrope_2eq_residual, x0=[t_walk, rho_walk], args=(d1, d2, z1, p_target, s_target), method="hybr", tol=1.0e-10)
        if not sol.success:
            print(f"{i:4d} {p_target:12.1f}   SOLVER FAILED -- walk stops here")
            print(f"       sol.message: {sol.message}")
            print(f"       last good (T,rho) before failure: T={t_walk:.4f}  rho={rho_walk:.4f}")
            break
        t_walk, rho_walk = float(sol.x[0]), float(sol.x[1])
        st = mix_state(d1, d2, t_walk, rho_walk, z1)
        s_ach = _mix_entropy_direct(d1, d2, t_walk, rho_walk, z1)
        drho = rho_walk - prev_rho
        # For vapor, rho should DECREASE monotonically as pressure drops -- flag any increase or big jump.
        flag = "  <<<< JUMP/REVERSAL" if drho >= 0 or abs(drho) > 0.5 * max(prev_rho, 1.0) else ""
        print(f"{i:4d} {p_target:12.1f} {t_walk:10.3f} {rho_walk:12.3f} {drho:+10.3f} {st.p_pa:12.1f} {s_ach:10.4f}  {str(sol.success)[:1]:>5}{flag}")
        prev_rho = rho_walk

print("\nDone. Density (rho) should be monotonically DECREASING as P_target decreases")
print("(superheated vapor gets less dense as pressure drops at constant entropy).")
print("Any 'JUMP/REVERSAL' flag marks a step where the solver likely left the true vapor branch.")
