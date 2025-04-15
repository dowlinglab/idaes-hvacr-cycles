# idaes-hvacr-cycles

Optimization models for HVAC and R applications in IDAES/Pyomo

This repository for created by Prof. Alexander Dowling at the University of Notre Dame. This work is part of the NSF ERC Environmentally Applied Refrigerant Technology Hub (EARTH).

This repository is an "alpha" version; please use the code with caution.

## Version

`simple_refridgeration.ipynb`
* Initial implementation for R134a
* PH state variables
* Initializes but does not solve... small infeasibility issue

`simple_refridgeration2.ipynb`
* Switching from PH to PTx state variables
* Got it to work!

`simple_refridgeration3.ipynb`
* Refactored to use `vapor_compression.py`
* Determined the complementary-like constraints was causing convergence issues (multipliers cycling)
* Improved initialization
* Implemented PH state variables
* Determined sensitivity analysis did not work because evaporator temperature was too low
* Determined subcooling/supercooling is important to make the sensitivity analysis reliable

Enhancement Ideas:
* Try CoolProp too
* Try other refrigerants