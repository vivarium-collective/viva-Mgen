"""Tests for the single reusable Mycoplasma genitalium composite generator."""

import pytest
from process_bigraph import Composite, gather_emitter_results

import viva_mgen  # noqa: F401  (fires generator registration)
from viva_mgen.core import build_core
from viva_mgen.composites import mycoplasma_genitalium, build_mgen


def test_generator_registered():
    from process_bigraph.composite_generator import _REGISTRY
    have = {eid.split(".")[-1] for eid in _REGISTRY}
    assert "mycoplasma_genitalium" in have, f"missing; have {sorted(have)[:8]}"


@pytest.fixture(scope="module")
def core():
    return build_core()


def test_composite_builds_and_runs(core):
    doc = mycoplasma_genitalium(core)
    sim = Composite({"state": doc}, core=core)
    sim.run(3.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    assert len(rows) >= 1
    last = rows[-1]
    # the integrated cell exposes metabolism, mass, and cell-cycle observables
    for k in ("mass", "growth_fraction", "atp_production", "replicated_fraction"):
        assert k in last


def test_gene_knockout_collapses_growth(core):
    doc = build_mgen(core, disrupted_genes=["MG_023"])
    sim = Composite({"state": doc}, core=core)
    sim.run(2.0)
    last = gather_emitter_results(sim)[("emitter",)][-1]
    assert last["growth_fraction"] < 0.05


def test_reaction_bound_scale_param(core):
    # a limiting reaction throttled hard reduces growth relative to wild-type
    doc = build_mgen(core, reaction_bound_scale={"EX_leu_DASH_L_e": 0.01})
    sim = Composite({"state": doc}, core=core)
    sim.run(1.0)
    last = gather_emitter_results(sim)[("emitter",)][-1]
    assert last["growth_fraction"] < 1.0
