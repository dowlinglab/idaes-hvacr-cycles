# Variant Notes: cycle_dx_hx_ideal

- Baseline copy with no intended behavior change in HX thermodynamics.
- Introduces common HX API:
  - `solve_evaporator(state_in, air_in, cfg, refrigerant)`
  - `solve_condenser(state_in, air_in, cfg, refrigerant)`
- Cycle solver calls only these functions so HX models are swappable.
