# Saturation Solver Review (Pure Fluid)

## Flowsheet

```text
[1] Inputs: fluid, T-grid
        |
        v
[2] JSON load + critical conversion (MW, Tc, rhoc_mol)
        |
        v
[3] Initial guesses (aux delta_l_sat_approx / delta_v_sat_approx)
        |
        v
[4] Nonlinear solve at each T:
    r1 = P(T,rho_l)-P(T,rho_v)
    r2 = g(T,rho_l)-g(T,rho_v)
    Newton + damping, scipy fallback
        |
        v
[5] Compute outputs:
    p_sat, h_l, h_v
        |
        v
[6] Continuation to next T and save dome arrays
```

## Representative States (r1234ze)

| T [K] | rho_l [mol/m^3] | rho_v [mol/m^3] | p_sat [kPa] | h_l [kJ/kg] | h_v [kJ/kg] |
|---:|---:|---:|---:|---:|---:|
| 240.000 | 11675.943 | 26.928 | 66840.839 | 408.070 | 361.628 |
| 305.333 | 4032.557 | 454.229 | 952.232 | 563.289 | 407.189 |
| 370.667 | 3060.214 | 3060.214 | 3344.801 | 438.018 | 438.018 |

## Notes for Signoff

- Solver equations and derivatives are implemented in `helmholtz_saturation.py`.
- Robust fallback path currently allows auxiliary-density fallback when strict root solve fails.
- Near-critical branch merging is detected when `|rho_l-rho_v|` collapses.
