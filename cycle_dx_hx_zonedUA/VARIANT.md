# Variant Notes: cycle_dx_hx_zonedUA

- Implements sequential DX zoned HX models.
- Evaporator zones: `tp`, `sh`.
- Condenser zones: `ds`, `tp`, `sc`.
- Later zones are naturally disabled when earlier zones are UA-limited.
- Keeps compressor and expansion modeling unchanged from other variants.
