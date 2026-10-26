# IDAES vapor-compression cycle notes (single stage + CO2 cascade)

Author: Shilpa Narasimhan / Support: Claude AI
Repo: ~/idaes-hvacr-cycles. Detailed derivations: Cascade/tutorial/cycle_tutorial.tex.

## Setpoints (collaborator-sourced, Cascade/diagnose_45C.py)
- Box/CO2 evaporating T_sat = -29 C. Condenser T_sat = ambient + 9 C.
- 21 C condensing limit case (head pressure control, Copeland): effective_ambient = max(ambient, 12 C).
- SH, SC windows 3..8 K (CO2 subcool 0..8 K). Compressor isentropic efficiency 0.80.
- Cascade HX: 3 K approach (T_CAS,E vs CO2 condensing), 1 K cold end. Pressure-ratio cap 20.
- Part load: PLF = 1 - CD(1 - PLR) = 0.9675 with CD = 0.13, PLR = 0.75.
- Carnot: 244.15 / (T_cond - 244.15), T_cond = ambient + 9 (or max(ambient + 9, 21)) in K.

## Model structure
- Single stage: evaporator, compressor, condenser, expansion_valve. Flow fixed 1 kg/s.
  DOF = 2 (SH, SC). Pinned: flow, evap T_sat, cond T_sat, efficiency, isenthalpic valve.
- Cascade: 8 units = hot_ and cold_ versions of the four. The cascade HX is
  hot_evaporator + cold_condenser, joined by cascade_energy_balance,
  cascade_approach_constraint (3 K), cascade_cold_end_constraint (1 K).
  flow_fixed_role="cold": CO2 flow = 1 kg/s, hot flow set by the HX balance.
  DOF = 5 (SH, SC each loop + T_CAS,E).
- Isentropic outlet comes from the ghost state properties_isentropic; h2 = h1 + (h2s - h1)/eta.
- Helmholtz state vars: P-h or TPX. R515B: custom R515BParameterBlock, per-loop modes
  (R515B P-h, CO2 TPX). R1234yf: custom JSON (R1234yf/r1234yf_general_helmholtz.json).
  R134a, CO2, R32: IDAES built-ins.
- Pressure bounds are derived from each fluid's envelope (CoolProp) over cascade_T_window_C;
  use (-26, 20) C so the bound never binds at cold ambients.

## Key files
- Cascade/vapor_compression_cascade.py (engine), run_cascade_sweep_hpc.py (R134a, R1234yf),
  run_cascade_sweep_r515b.py, retry_r515b_points.py, run_cascade_sweep_r32.py (separate project).
- R1234yf/run_single_stage_sweep.py (single stage, all fluids).
- Cascade/physics_checks.py: cascade_checks / single_checks.
- Cascade/heat_balance.py: cascade_balance / single_balance (m, T, P, h per pipe, duties).
- Cascade/make_cop_figures.py: COP vs ambient figures (full and part load + Carnot).

## Basis of extensive numbers
- Single stage: per 1 kg/s refrigerant. Cascade: per 1 kg/s CO2 (cold loop); the hot-loop
  flow differs per fluid. COP and ratios (W/Q_box, Q_air/Q_box) don't depend on the basis.

## Results status (10/2026)
- 204 points (3 fluids x cascade/single x 2 cases x 17 ambients) converged, passed physics
  checks, and matched an independent CoolProp brute force (worst 0.0004%).
- Cascade vs single-stage crossover: R134a 26.4 C, R1234yf 19.4 C, R515B 23.0 C.
- Next: CAPEX/OPEX, TEWI (waiting on Zhuorui's weather data, annual kWh, cooling load).
