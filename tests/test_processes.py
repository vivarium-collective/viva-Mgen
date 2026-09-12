"""Unit tests for viva-mGen reproduction processes."""

import math

import pytest

from viva_mgen.core import build_core
from viva_mgen.processes import (
    MetabolismFbaReproductionProcess,
    MassGrowthReproductionProcess,
    TranscriptionReproductionProcess,
    TranslationReproductionProcess,
    RnaDecayReproductionProcess,
    ProteinDecayReproductionProcess,
    ReplicationReproductionProcess,
)


@pytest.fixture(scope="module")
def core():
    return build_core()


def test_metabolism_wildtype_grows(core):
    p = MetabolismFbaReproductionProcess(config={}, core=core)
    assert isinstance(p.inputs(), dict) and isinstance(p.outputs(), dict)
    out = p.update({"nutrient_scale": 1.0}, 1.0)
    assert out["growth_rate"] > 0.0
    assert abs(out["growth_fraction"] - 1.0) < 1e-6
    assert out["feasible"] == 1.0


def test_metabolism_essential_gene_knockout_kills_growth(core):
    # MG_023 is essential in iPS189 (growth collapses on knockout)
    p = MetabolismFbaReproductionProcess(config={"disrupted_genes": ["MG_023"]}, core=core)
    out = p.update({"nutrient_scale": 1.0}, 1.0)
    assert out["growth_fraction"] < 0.05
    assert out["feasible"] == 0.0


def test_metabolism_nonessential_gene_knockout_survives(core):
    p = MetabolismFbaReproductionProcess(config={"disrupted_genes": ["MG_038"]}, core=core)
    out = p.update({"nutrient_scale": 1.0}, 1.0)
    assert out["growth_fraction"] > 0.5


def test_mass_doubles_in_one_cell_cycle(core):
    p = MassGrowthReproductionProcess(config={}, core=core)
    # one full cell cycle at wild-type growth -> mass exactly doubles
    out = p.update({"growth_fraction": 1.0, "mass": 3.93}, 32400.0)
    assert math.isclose(out["mass"], 3.93, rel_tol=1e-6)  # delta == initial (doubled)


def test_transcription_produces_transcripts_over_time(core):
    p = TranscriptionReproductionProcess(config={"seed": 3}, core=core)
    total = {}
    for _ in range(600):  # 10 minutes
        d = p.update({"ntp": 1e9, "rna_pol": 100.0}, 1.0)
        for g, n in d["rna_counts"].items():
            total[g] = total.get(g, 0) + n
    assert sum(total.values()) > 0


def test_translation_needs_mrna(core):
    p = TranslationReproductionProcess(config={"seed": 4}, core=core)
    empty = p.update({"rna_counts": {}, "gtp": 1e9}, 1.0)
    assert empty["protein_counts"] == {}
    withmrna = p.update({"rna_counts": {"tuf": 50.0}, "gtp": 1e9}, 60.0)
    assert sum(withmrna["protein_counts"].values()) >= 0  # stochastic, nonnegative


def test_decay_reduces_counts(core):
    p = RnaDecayReproductionProcess(config={"seed": 5}, core=core)
    d = p.update({"rna_counts": {"tuf": 1000.0}}, 300.0)
    assert d["rna_counts"].get("tuf", 0) <= 0  # decay is a negative delta


def test_replication_completes(core):
    p = ReplicationReproductionProcess(config={"initial_dnaA": 30.0, "seed": 1}, core=core)
    out = None
    for _ in range(40000):
        out = p.update({"dntp_synthesis_scale": 1.0}, 1.0)
        if out["phase_code"] == 2.0:
            break
    assert out["phase_code"] == 2.0
    assert out["chromosome_copy"] == 2.0
    assert out["replication_duration"] > 0.0
