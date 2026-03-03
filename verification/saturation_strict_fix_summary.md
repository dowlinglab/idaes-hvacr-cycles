# Saturation Solver Strict-Gate Fix Summary

Author: Shilpa Narasimhan  
Technical support: Codex (version GPT-5)  
QA/Testing Responsibility: Shilpa  
Creation date: 2026-03-02  
Purpose of file: Summarize strict saturation acceptance fix, outputs, and diagnostics.  
Dependencies: helmholtz_saturation.py outputs in verification/  
Context reference: PROJECT_CONTEXT.md

## Fix Description
Implemented strict acceptance gates in the pure-fluid Helmholtz saturation solver:
- `r_P = |P_l-P_v| / max(1.0, 0.5*(P_l+P_v)) <= 1e-6`
- `r_mu = |g_l-g_v| / (R_u*T) <= 1e-6`
- `rho_l > rho_v*(1+1e-8)`

Auxiliary fallback values are now used only for retry seeding/diagnostics and are never accepted as converged dome points.

## Run Summary
- Fluid: `r1234ze`
- Grid: `T=240..372 K`, `n=200`
- Attempted: `200`
- Converged: `199`
- Failed: `1`
- CSV: [`verification/saturation_dome_run.csv`](/Users/snarasi2/idaes-hvacr-cycles/verification/saturation_dome_run.csv)
- Log: [`verification/log_saturation_solver.txt`](/Users/snarasi2/idaes-hvacr-cycles/verification/log_saturation_solver.txt)

## Plots
### Clean dome (converged only)
![clean dome](/Users/snarasi2/idaes-hvacr-cycles/verification/ph_dome_clean.png)

### Dome with failures marked (red X)
![dome with failures](/Users/snarasi2/idaes-hvacr-cycles/verification/ph_dome_with_failures.png)

## Three Diagnostics
1. `T=240.0 K` (previously quoted failure case, now fixed):
   - Status: `CONVERGED`
   - `P_sat=6.1495 bar`, `r_P=3.98e-15`, `r_mu=4.84e-16`
   - Jacobian condition number: `1.0839e+03`
   - Final Newton step norm: `1.13e-12`

2. `T=370.010050 K` (full-dome failed point):
   - Status: `DIVERGED`
   - `r_P=3.1735e-01`, `r_mu=9.4479e-02`
   - Jacobian condition number: `1.2813e+16`
   - Final Newton step norm: `9.5649e-02`
   - Notes: `attempt=1; line-search rejected step; fallback_seed_used`

3. Pre-fix baseline reference at `T=240.0 K` (before strict-gate fix):
   - Reported non-physical mismatch: `P_l≈1336.29 bar`, `P_v≈0.53 bar`
   - Normalized mismatches: `r_P≈9.996e-01`, `r_mu≈8.766`
   - This point was previously passing through fallback acceptance and was incorrectly plotted.
