# PLR Cold-Storage Brainstorming Notes (2026-03-04)

## Scope
This note captures brainstorming decisions for the PLR refrigeration runs so context is preserved for manual mass/energy checks.

## Agreed Modeling Direction
1. Do not hard-fix evaporator saturation temperature to one value.
2. Use an evaporator saturation band tied to cold-storage setpoint:
   - `T_evap_sat = T_cold_storage - (8 to 10 C)`
   - Example for `T_cold_storage = -20 C`: `T_evap_sat in [-30, -28] C`.
3. Use a condenser saturation band tied to ambient:
   - `T_cond_sat = T_ambient + (8 to 10 C)`.
4. Ambient sweep window for this cold-storage analysis:
   - `T_ambient in [10, 25] C`.
5. Keep compressor vapor-only outlet guard (physical and numerical guardrail).
6. Keep compressor isentropic efficiency aligned at `eta_isen = 0.75` for both base and PLR classes.

## Clarified Physical Interpretation
1. Base `vapor_compression.py` couples ambient explicitly at condenser side only.
2. Evaporator side in base model is represented by refrigerant-side constraints (sat/superheat), not an explicit air-side HX equation.
3. Saturation temperature and pressure are thermodynamically linked:
   - `T_sat = Tsat(P)` and `P_sat = Psat(T)`.

## Manual Balance Checks (per solved point)
Use absolute values with consistent sign conventions if needed.

1. Refrigeration COP:
   - `COP = Q_evap / W_comp`

2. First-law closure around cycle:
   - `Q_cond ~= Q_evap + W_comp`
   - Residual:
     - `r_E = Q_cond - (Q_evap + W_comp)`
   - Relative residual:
     - `r_E_rel = r_E / max(1e-6, |Q_cond|)`

3. Evaporator approach check (to cold storage):
   - `DeltaT_evap = T_cold_storage - T_evap_sat`
   - Target band: `8 to 10 C`.

4. Condenser approach check (to ambient):
   - `DeltaT_cond = T_cond_sat - T_ambient`
   - Target band: `8 to 10 C`.

5. Compressor pressure ratio check:
   - `PR = P_high / P_low`
   - Track against practical expected range per refrigerant and operating condition.

## Implementation Options (not yet applied)
1. Hard band constraints on `T_evap_sat` and `T_cond_sat`.
2. Soft-penalty objective around band midpoints.
3. Explicit HX UA/NTU equations so approach values emerge from heat transfer.

## Important Run-Logic Note
If using fallback solve paths, pass the full point specification again.
Do not call `set_specifications()` with only a debug flag, since that can silently revert to defaults.
