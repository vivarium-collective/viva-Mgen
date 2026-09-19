from viva_mgen.core import build_core
from viva_mgen.structural.pack_step import MgenStructuralStep, register_structural


def test_step_registers_and_declares_io():
    core = build_core()
    register_structural(core)
    step = MgenStructuralStep({"counts_source": "maritan", "top_n": 8}, core=core)
    outs = step.outputs()
    assert "pack_path" in outs


def test_step_declares_inputs():
    core = build_core()
    register_structural(core)
    step = MgenStructuralStep({"counts_source": "maritan", "top_n": 8}, core=core)
    ins = step.inputs()
    assert isinstance(ins, dict)


def test_register_structural_adds_link():
    core = build_core()
    register_structural(core)
    assert "MgenStructuralStep" in core.link_registry
    assert "viva_mgen.structural.pack_step.MgenStructuralStep" in core.link_registry
