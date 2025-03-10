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
* Sensitivity analysis shows the model is brittle
* TODO: subcool and superheat by 3 deg C
* TODO: clean up/streamline initialization routine such that the user just speficies two temperatures
* TODO: Try CoolProp too
* TODO: try other refrigerants