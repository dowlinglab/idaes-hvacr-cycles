"""
Diagnostic script: trace the liquid-side isentrope walk for a single
target entropy (default 0.26 Btu/lb-R) step by step.

Purpose
-------
Isolate exactly where compute_isentrope_liquid_side()'s continuation
walk goes wrong: which anchor row gets picked, whether the initial
solve moves the seed far from where it started, and -- most
importantly -- whether the solved density sequence along the pressure
walk is monotonically increasing (physically correct for a subcooled
liquid at constant entropy) or jumps/reverses (signature of the solver
landing on a spurious/degenerate root, same failure class as the
dome solver's pre-fix behavior).

Usage
-----
Run from the same directory as mixture_isentrope_validation.py (i.e.
R515B_props_validated/, which now also carries its own copy of
linear_model_codex.py alongside it -- same plain-import style as
mixture_isentrope_validation.py itself, no path tricks needed):
    python3 diagnose_isentrope_026.py
Change S_TARGET_BTU below to inspect a different isentrope (e.g. 0.22,
which is the one known-good baseline, for comparison).

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
    ISENTROPE_P_MAX_PA,
    ISENTROPE_N_POINTS,
)

# ---- Match this to whatever dome run you've been comparing against ----
FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
TMIN, TMAX, N = 255.4, 349.8, 80
S_TARGET_BTU = 0.26  # <-- change this to inspect a different isentrope
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

converged_bubble = [r for r in bubble_rows if r["status"] == "CONVERGED"]
print(f"{len(converged_bubble)}/{len(bubble_rows)} bubble rows converged.\n")

# ---- Step 1: anchor selection (mirrors compute_isentrope_liquid_side) ----
print("=== Step 1: anchor candidates (rows where s_l(T) > s_target) ===")
print(f"{'T_K':>10} {'s_l':>12} {'gap':>12}  best?")
best_row = None
best_gap = None
for r in converged_bubble:
    t_k = r["T_K"]
    s_l = _mix_entropy_direct(d1, d2, t_k, r["rho_l_molm3"], z1)
    if s_l <= s_target:
        continue
    gap = s_l - s_target
    is_best = best_gap is None or gap < best_gap
    if is_best:
        best_gap = gap
        best_row = r
    print(f"{t_k:10.3f} {s_l:12.4f} {gap:12.4f}  {'<-- BEST' if is_best else ''}")

if best_row is None:
    print("\nNo anchor row found -- this isentrope never reaches the subcooled branch in this T range.")
    raise SystemExit(0)

t_k = best_row["T_K"]
rho_seed = best_row["rho_l_molm3"]
p_seed = best_row["P_Pa"]
print(f"\nAnchor selected: T={t_k:.4f} K   rho_seed={rho_seed:.4f} mol/m3   p_seed={p_seed:.2f} Pa")

# ---- Step 2: initial solve at P_target = p_seed ----
print("\n=== Step 2: initial solve (should be near-trivial since seed IS the answer) ===")
sol0 = root(_isentrope_2eq_residual, x0=[t_k, rho_seed], args=(d1, d2, z1, p_seed, s_target), method="hybr", tol=1.0e-10)
print(f"success={sol0.success}  T_solved={sol0.x[0]:.4f}  rho_solved={sol0.x[1]:.4f}")
print(f"  moved from seed by: dT={sol0.x[0]-t_k:+.4f} K   drho={sol0.x[1]-rho_seed:+.4f} mol/m3")
if not sol0.success:
    print("INITIAL SOLVE FAILED -- stopping here.")
    raise SystemExit(0)

t_walk, rho_walk = float(sol0.x[0]), float(sol0.x[1])
st0 = mix_state(d1, d2, t_walk, rho_walk, z1)
s0 = _mix_entropy_direct(d1, d2, t_walk, rho_walk, z1)
print(f"  Achieved P={st0.p_pa:.2f} Pa (target {p_seed:.2f})   Achieved s={s0:.4f} (target {s_target:.4f})")

# ---- Step 3: the pressure walk itself ----
p_start = float(st0.p_pa)
print(f"\n=== Step 3: pressure walk, {p_start:.1f} Pa -> {ISENTROPE_P_MAX_PA:.1f} Pa, {ISENTROPE_N_POINTS} log-spaced steps ===")
print(f"{'step':>4} {'P_target':>12} {'T_K':>10} {'rho':>12} {'dRho':>10} {'P_achieved':>12} {'s_achieved':>10} {'ok':>5}")

prev_rho = rho_walk
for i, p_target in enumerate(np.geomspace(p_start, ISENTROPE_P_MAX_PA, ISENTROPE_N_POINTS)[1:], start=1):
    sol = root(_isentrope_2eq_residual, x0=[t_walk, rho_walk], args=(d1, d2, z1, p_target, s_target), method="hybr", tol=1.0e-10)
    if not sol.success:
        print(f"{i:4d} {p_target:12.1f}   SOLVER FAILED -- walk stops here")
        break
    t_walk, rho_walk = float(sol.x[0]), float(sol.x[1])
    st = mix_state(d1, d2, t_walk, rho_walk, z1)
    s_ach = _mix_entropy_direct(d1, d2, t_walk, rho_walk, z1)
    drho = rho_walk - prev_rho
    flag = "  <<<< JUMP/REVERSAL" if drho <= 0 or abs(drho) > 0.5 * max(prev_rho, 1.0) else ""
    print(f"{i:4d} {p_target:12.1f} {t_walk:10.3f} {rho_walk:12.3f} {drho:+10.3f} {st.p_pa:12.1f} {s_ach:10.4f}  {str(sol.success)[:1]:>5}{flag}")
    prev_rho = rho_walk

print("\nDone. Density (rho) should be monotonically INCREASING as P_target increases")
print("(subcooled liquid gets denser under higher pressure at constant entropy).")
print("Any 'JUMP/REVERSAL' flag marks a step where the solver likely left the true liquid branch.")
