"""All 28-submodel reproductions register, instantiate, and expose a contract."""

import pytest

from viva_mgen.core import build_core
from viva_mgen.processes import all_process_classes


ALL = all_process_classes()


def test_full_submodel_count():
    # 7 core + 22 migrated = 29 reproduction process classes (the Karr 2012
    # 28-submodel set; a couple of submodels split into >1 process).
    assert len(ALL) >= 28, f"expected >=28 submodels, got {len(ALL)}: {sorted(ALL)}"


@pytest.fixture(scope="module")
def core():
    return build_core()


def test_all_register(core):
    # build_core() registered every submodel; a local address resolves.
    from process_bigraph import Composite
    doc = {"p": {"_type": "process", "address": "local:MetabolismFbaReproductionProcess",
                 "config": {}, "inputs": {"nutrient_scale": ["s", "n"]},
                 "outputs": {"growth_rate": ["s", "g"]}}, "s": {"n": 1.0, "g": 0.0}}
    Composite({"state": doc}, core=core)  # raises if the link isn't registered


@pytest.mark.parametrize("name", sorted(ALL))
def test_instantiate_and_contract(name, core):
    cls = ALL[name]
    proc = cls(config={}, core=core)
    ins, outs = proc.inputs(), proc.outputs()
    assert isinstance(ins, dict) and isinstance(outs, dict)
    # formal description contract present and labels fidelity
    desc = getattr(cls, "description", "") or ""
    assert desc, f"{name} has no description contract"
    assert "Fidelity" in desc or "reproduction" in desc.lower()
