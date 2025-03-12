# idaes-hvacr-cycles
Optimization models for HVAC and R applications in IDAES/Pyomo

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
* TODO: Switch from Tsat to Psat constraint (not sure if this matters)
* TODO: Remove obsolete initialization options, streamline code
* TODO: Add more comments
* TODO: Determine why the sensitivity analysis does not always converge or gets stuck in local solutions. Ideally, every solve in the sensitivity analysis should be very clean.


Enhancement Ideas:
* Clean up/streamline initialization routine such that the user just speficies two temperatures
* Create another version using PH state variables
* Try CoolProp too
* Try other refrigerants