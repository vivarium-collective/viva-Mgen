"""Tests that the composite generators are registered and buildable/runnable."""

import pytest
from process_bigraph import Composite, gather_emitter_results

import viva_mgen  # noqa: F401  (fires generator registration)
from viva_mgen.core import build_core
from viva_mgen import composites as C


GENERATOR_NAMES = [
    "fig1_architecture", "fig2_growth", "fig3_expression", "fig5_energy",
    "fig4_cell_cycle", "fig6_gene_essentiality", "fig7_kinetic_parameters",
]


def test_all_generators_registered():
    from process_bigraph.composite_generator import _REGISTRY
    have = {eid.split(".")[-1] for eid in _REGISTRY}
    for name in GENERATOR_NAMES:
        assert name in have, f"{name} missing from registry; have {sorted(have)[:8]}"


@pytest.fixture(scope="module")
def core():
    return build_core()


@pytest.mark.parametrize("name", GENERATOR_NAMES)
def test_composite_builds_and_runs(name, core):
    gen = getattr(C, name)
    doc = gen(core)
    sim = Composite({"state": doc}, core=core)
    sim.run(3.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    assert len(rows) >= 1


def test_fig6_knockout_collapses_growth(core):
    doc = C.fig6_gene_essentiality(core, disrupted_gene="MG_023")
    sim = Composite({"state": doc}, core=core)
    sim.run(2.0)
    last = gather_emitter_results(sim)[("emitter",)][-1]
    assert last["growth_fraction"] < 0.05
