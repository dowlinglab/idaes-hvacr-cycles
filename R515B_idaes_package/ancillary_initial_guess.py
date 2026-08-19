"""
Fast, non-iterative initial-guess machinery for R-515B saturated
liquid/vapor molar density, built by composing the PURE-COMPONENT
ancillary saturation-density correlations already shipped in IDAES's own
Helmholtz JSON parameter files (r1234ze.json, r227ea.json) through this
project's existing Bell (2023) reducing-function chain-rule mapping.

Why this exists
----------------
Neither `mixture_fully_validated.py` (the oracle) nor `linear_model_codex.py`
defines or uses any ancillary/shortcut correlation for R-515B's own
saturation properties -- confirmed by reading both files in full this
session. Every bubble/dew point in the oracle is found by a real iterative
solve (3-equation VLE via `scipy.optimize.least_squares`), seeded either by
continuation from a neighboring temperature or a rough 80%/1%-of-critical-
density heuristic for the very first point. That works fine inside a single
Python process sweeping a whole temperature grid, but the new IDAES/Pyomo
property package's own `initialize()` routine (MASTER TASK spec rule 47)
will need a good starting guess for a cold solve with no "previous point" to
continue from -- this module builds that guess from real, already-published
physics instead of an arbitrary number.

What "ancillary equation" means here
-------------------------------------
An ancillary equation is a fast, closed-form (non-iterative) approximation
to a saturation property, standard in REFPROP-style multiparameter EOS
packages, used ONLY to seed a real iterative solve -- never to replace it.
IDAES's own Helmholtz JSON files for BOTH pure fluids used here (confirmed
by direct inspection, `/root/.idaes/bin/helm_data/{r1234ze,r227ea}.json`,
`aux` section) ship exactly this: `delta_l_sat_approx` (type 1) and
`delta_v_sat_approx` (type 2), giving each PURE fluid's own saturated
reduced liquid/vapor density as a function of its own reduced temperature.
No R-515B-mixture-level ancillary exists anywhere (published or in this
project) -- there is no Wagner-type fit for this specific blend. This
module builds a mixture-level estimate BY COMPOSING the two pure-fluid
ancillaries, not by inventing a new mixture correlation from scratch.

Exact formulas used (confirmed directly from IDAES's own implementation,
`idaes/models/properties/general_helmholtz/expressions/sat_delta_approx.py`,
NOT guessed or reverse-engineered from the coefficients alone):
    type 1 (delta_l_sat_approx): delta_sat = c + sum_i[ n_i * (1 - 1/tau)^t_i ]
    type 2 (delta_v_sat_approx): delta_sat = c * exp(sum_i[ n_i * (1 - 1/tau)^t_i ])
where `tau` is IDAES's own convention tau = Tc/T (confirmed consistent with
this project's own tau_i = (Tc_i/Tred_mix)*tau_mix chain-rule mapping,
already used throughout `mixture_fully_validated.py`/`linear_model_codex.py`
-- same tau definition, so no unit-convention mismatch). Both fluids here
use c=1 for both types (confirmed by direct inspection, not assumed).

Combination rule (mixture-level guess from the two pure-fluid estimates)
--------------------------------------------------------------------------
This is NEW engineering, not present in the oracle, and is deliberately the
simplest defensible choice: ideal volume-additivity mixing. At a given
(T, z1), each pure fluid's own tau_i=Tc_i/T feeds its own ancillary to get
its own saturated molar volume (1/(delta_sat_i * rhoc_i_mol)), and the
mixture's saturated molar volume guess is the mole-fraction-weighted sum of
those pure molar volumes (v_mix = z1*v1 + z2*v2), inverted back to a molar
density. This is only ever used as a SEED -- the real accepted saturation
state still requires the full VLE solve (this project's own
`solve_bubble_at_t`/`solve_dew_at_t`, or the new package's own IPOPT-based
equivalent). Validated below (see `validate_ancillary_guess.py`) against
the oracle's own real converged bubble/dew densities to characterize how
close a seed this actually produces -- not to claim it as an independently
accurate saturation model.

References
----------
- IDAES `general_helmholtz` ancillary saturated-density expressions,
  `expressions/sat_delta_approx.py` (types 1-3; types 1/2 used here).
- Bell, I. H. (2023), J. Phys. Chem. Ref. Data 52(1), 013101 -- the
  reducing-function tau_i=(Tc_i/Tred_mix)*tau_mix chain rule this module
  reuses (imported from `linear_model_codex.py`, an existing, unmodified
  project dependency -- NOT the oracle file).

Read-only against the oracle in this module: none at import time (only
`validate_ancillary_guess.py`, a separate file, imports the oracle for
comparison purposes).
"""

import json
from pathlib import Path
from typing import Dict, Tuple

from linear_model_codex import (
    BELL_2023_R1234ZE_R227EA,
    bell2023_Tred_vred,
    load_idaes_helmholtz_json,
    mw_from_json,
)


def _eval_delta_sat_type01(aux: Dict, tau: float) -> float:
    """delta_sat = c + sum_i[ n_i * (1 - 1/tau)^t_i ] -- IDAES ancillary type 1
    (used here for delta_l_sat_approx), formula confirmed directly from
    idaes/models/properties/general_helmholtz/expressions/sat_delta_approx.py."""
    c = float(aux["c"])
    n = aux["n"]
    t = aux["t"]
    theta = 1.0 - 1.0 / tau
    return c + sum(float(n[k]) * (theta ** float(t[k])) for k in n)


def _eval_delta_sat_type02(aux: Dict, tau: float) -> float:
    """delta_sat = c * exp(sum_i[ n_i * (1 - 1/tau)^t_i ]) -- IDAES ancillary
    type 2 (used here for delta_v_sat_approx), same source as type01 above."""
    import math
    c = float(aux["c"])
    n = aux["n"]
    t = aux["t"]
    theta = 1.0 - 1.0 / tau
    return c * math.exp(sum(float(n[k]) * (theta ** float(t[k])) for k in n))


_EVAL_BY_TYPE = {1: _eval_delta_sat_type01, 2: _eval_delta_sat_type02}


def pure_fluid_sat_densities_molm3(d: Dict, t_k: float) -> Tuple[float, float]:
    """
    Evaluate ONE pure fluid's own ancillary saturated liquid/vapor molar
    density at temperature t_k, using its own tau_i = Tc_i/T (IDAES's own
    tau convention) and its own aux.delta_{l,v}_sat_approx correlation type.

    Inputs
    ------
    d : dict -- parsed IDAES Helmholtz JSON for one pure fluid.
    t_k : float [K]

    Outputs
    -------
    (rho_l_molm3, rho_v_molm3) : tuple[float,float] [mol/m^3]
      This pure fluid's OWN ancillary-estimated saturated liquid/vapor molar
      density at t_k -- not yet combined into a mixture estimate.

    Failure modes
    -------------
    - KeyError if the JSON lacks an aux/delta_{l,v}_sat_approx section
      (not expected for r1234ze.json/r227ea.json, confirmed present).
    - Physically meaningless (but not raising) if t_k > Tc (tau<1, theta<0,
      fractional powers of a negative theta can go complex/NaN) -- ancillary
      equations are only valid below Tc by construction; callers must not
      use this above the pure fluid's own Tc.
    """
    tc = float(d["basic"]["Tc"])
    mw = mw_from_json(d)
    rhoc_mol = float(d["basic"]["rhoc"]) / mw
    tau = tc / t_k

    aux_l = d["aux"]["delta_l_sat_approx"]
    aux_v = d["aux"]["delta_v_sat_approx"]
    fn_l = _EVAL_BY_TYPE[int(aux_l["type"])]
    fn_v = _EVAL_BY_TYPE[int(aux_v["type"])]
    delta_l = fn_l(aux_l, tau)
    delta_v = fn_v(aux_v, tau)
    return float(delta_l * rhoc_mol), float(delta_v * rhoc_mol)


def mixture_ancillary_saturation_guess_molm3(
    d1: Dict, d2: Dict, z1: float, t_k: float
) -> Tuple[float, float]:
    """
    Fast, non-iterative INITIAL-GUESS ONLY for R-515B's (mixture) saturated
    liquid/vapor molar density at (T, z1), built by combining each pure
    fluid's own ancillary saturation-density estimate (see
    `pure_fluid_sat_densities_molm3`) via ideal volume-additivity mixing:
        v_mix = z1*v1_sat + z2*v2_sat  (molar volume, mole-fraction-weighted)
        rho_mix = 1/v_mix
    applied separately to the liquid-branch and vapor-branch pure estimates.

    THIS IS NOT A VALIDATED SATURATION MODEL -- it is a seed-quality
    estimate for warm-starting a real VLE solve (this project's own
    `solve_bubble_at_t`/`solve_dew_at_t`, or the new IDAES package's
    equivalent). See `validate_ancillary_guess.py` for an empirical
    characterization of how close this gets to the oracle's own real
    converged densities across the practical temperature range.

    Inputs
    ------
    d1, d2 : dict -- parsed IDAES Helmholtz JSON, fluid 1 (r1234ze) / 2 (r227ea).
    z1 : float [mol/mol] -- overall (feed) mole fraction of fluid 1.
    t_k : float [K] -- must be below both pure fluids' critical temperatures
      (and, in practice, comfortably below the true mixture Tc) for the
      underlying pure-fluid ancillaries to be physically meaningful.

    Outputs
    -------
    (rho_l_guess_molm3, rho_v_guess_molm3) : tuple[float,float] [mol/m^3]
    """
    z2 = 1.0 - z1
    rho1_l, rho1_v = pure_fluid_sat_densities_molm3(d1, t_k)
    rho2_l, rho2_v = pure_fluid_sat_densities_molm3(d2, t_k)
    v_l_mix = z1 / rho1_l + z2 / rho2_l
    v_v_mix = z1 / rho1_v + z2 / rho2_v
    return float(1.0 / v_l_mix), float(1.0 / v_v_mix)
