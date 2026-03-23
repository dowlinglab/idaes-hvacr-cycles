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

## Addendum (2026-03-10): IDAES Condenser-Train Initialization Debug

### Current Debug Scope
- Active copy model: `vapor_compression_plr_hx_0d_cond3.py`.
- Working mode: one-ambient debug first, then sweeps only after one-point convergence.

### Key Observation
- Nonphysical condenser-air seeds (~500 K) appear during warm-start mapping when air inlet temperatures are reconstructed directly from zoned duties.

### Current Stabilization Direction
1. Fix only the first condenser-train air inlet temperature and composition.
2. Leave downstream condenser-air inlets unfixed so arc equalities determine continuity.
3. Seed downstream air temperatures to ambient-like values (not from raw `Q/C_air` reconstruction).
4. Keep refrigerant enthalpy seeds continuous across DS->COND->SC interfaces.
5. Avoid conflicting hardcoded inlet-enthalpy anchors during zoned warm-start usage.

### Reporting Rule
- For this stage, output PFD + stream table for one ambient point before generating overlay figures.

## Addendum (2026-03-10 Late): Why Air Temperatures Look Nonphysical

### Confirmed Status
- Current cond3 IDAES copy remains non-converged at one-point test (`R134a`, `Tamb=20 C`).
- Structural closure is currently square (`DOF=0`) after pressure/area cleanup.

### What Is Happening
1. The run fails before constraints are satisfied globally.
2. Reported values like `226.85 C` air outlet are failed-iterate values, not physical solutions.
3. Subcooler can appear to "pump" pressure at failed iterate because pressure constraints are not yet driven to zero residual.

### Key Residual Evidence
- `subcooler_hot_dp0` residual is large (order `1e6 Pa`).
- `P_high_comp_out` residual is large (order `1e5 Pa`).
- `subcooler.heat_transfer_equation` and delta-T residuals remain nonzero.
- SC duty at failed state is effectively zero (`Q_hot = Q_cold = 0`).

### Practical Meaning
- The model is not predicting hot air physically; it is numerically stalled.
- Next work should prioritize numerical stabilization of HX delta-T/LMTD path and initialization continuation.

## Addendum (2026-03-10 Late): SC Mass/Energy Balance Snapshot (Failed Iterate)

### Reported SC Balance Values
- Hot side: `m_dot=1.0 kg/s`, `h_in=241722.392 J/kg`, `h_out=236722.392 J/kg`, `Q_hot=0 W`.
- Cold side: `N_dot=356 mol/s`, `h_in=-250 J/mol`, `h_out=-250 J/mol`, `Q_cold=0 W`.
- Cold-side temperatures: `T_in=20.0 C`, `T_out=226.85 C` (nonphysical in this failed state).

### Manual Check Results
- Hot-side energy residual: `m_dot*(h_out-h_in)+Q_hot = -5000 W`.
- Cold-side energy residual (algebraic): `N_dot*(h_out-h_in)+Q_cold = 0`.
- Cross-side duty closure: `Q_hot + Q_cold = 0` (trivial at failed iterate).

### Interpretation
- SC is not physically solved yet.
- Current values are for diagnostic debugging only and should not be used as cycle performance outputs.
