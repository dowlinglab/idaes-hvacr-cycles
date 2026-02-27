"""Unit tests for PLR compressor flow scaling (no solves).

Run:
  pytest -q
  pytest -q -m slow
"""

import pytest
from pyomo.environ import value

idaes = pytest.importorskip("idaes")
pyomo = pytest.importorskip("pyomo")

from vapor_compression_plr import SimpleVaporCompressionCyclePLR, Mode


def test_default_plr_value():
    cycle = SimpleVaporCompressionCyclePLR("R134a", mode=Mode.IMPROVED_TPX)
    assert abs(cycle.plr - 0.75) < 1e-12


def test_cd_validation():
    with pytest.raises(ValueError, match="CD must be in"):
        SimpleVaporCompressionCyclePLR._compute_plf(0.5, 0.04)
    with pytest.raises(ValueError, match="CD must be in"):
        SimpleVaporCompressionCyclePLR._compute_plf(0.5, 0.51)

    assert SimpleVaporCompressionCyclePLR._compute_plf(0.5, 0.05) > 0.0
    assert SimpleVaporCompressionCyclePLR._compute_plf(0.5, 0.5) > 0.0


def test_plf_formula():
    assert SimpleVaporCompressionCyclePLR._compute_plf(1.0, 0.25) == 1.0
    assert abs(SimpleVaporCompressionCyclePLR._compute_plf(0.5, 0.1) - 0.95) < 1e-12
    assert abs(SimpleVaporCompressionCyclePLR._compute_plf(0.0, 0.5) - 0.5) < 1e-12


def test_cop_part_scaling():
    cycle = SimpleVaporCompressionCyclePLR("R134a", mode=Mode.IMPROVED_TPX)
    cycle.model.fs.COP_full.set_value(4.0)

    cycle.set_specifications(plr=0.3, Q_cool_rated=1.0, cd=0.1)
    expected_plf = 1.0 - 0.1 * (1.0 - 0.3)
    assert abs(value(cycle.model.fs.PLF) - expected_plf) < 1e-12
    assert abs(value(cycle.model.fs.COP_part) - (expected_plf * 4.0)) < 1e-12
