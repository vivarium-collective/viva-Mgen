import math
import pytest
from viva_mgen.processes.allocation import allocate, select_budget, demand_entry, AllocatorProcess
from viva_mgen.core import build_core


@pytest.fixture(scope="module")
def core():
    return build_core()


def test_allocate_surplus_grants_full():
    g = allocate(100.0, {"a": 10.0, "b": 20.0})
    assert g == {"a": 10.0, "b": 20.0}


def test_allocate_scarcity_is_proportional():
    g = allocate(30.0, {"a": 30.0, "b": 60.0})
    assert math.isclose(g["a"], 10.0) and math.isclose(g["b"], 20.0)
    assert sum(g.values()) <= 30.0 + 1e-9


def test_allocate_priority_weights():
    g = allocate(30.0, {"a": 30.0, "b": 30.0}, priorities={"a": 2.0, "b": 1.0})
    assert math.isclose(g["a"], 20.0) and math.isclose(g["b"], 10.0)


def test_allocate_zero_supply_grants_zero():
    g = allocate(0.0, {"a": 5.0})
    assert g == {"a": 0.0}


def test_allocate_ignores_nonpositive_and_conserves():
    g = allocate(10.0, {"a": -3.0, "b": 0.0, "c": 40.0})
    assert g["a"] == 0.0 and g["b"] == 0.0
    assert math.isclose(g["c"], 10.0)


def test_select_budget_absent_is_inf():
    assert select_budget({}, "x") == float("inf")
    assert select_budget({"x": 5.0}, "x") == 5.0


def test_demand_entry_clamps_negative():
    assert demand_entry("x", -2.0) == {"x": 0.0}


def test_allocator_replenishes_and_partitions(core):
    proc = AllocatorProcess(config={"pools": ["atp"]}, core=core)
    # Verify port contract
    ins = proc.inputs()
    outs = proc.outputs()
    assert ins["atp"] == "float"
    assert ins["atp_production"] == "float"
    assert ins["demand__atp"] == "map[float]"
    assert outs["atp"] == "float"
    assert outs["alloc__atp"] == "overwrite[map[float]]"
    # Test arithmetic: supply available = level(5) + production(25) = 30; demands 20+40=60 -> proportional
    state = {"atp": 5.0, "atp_production": 25.0,
             "demand__atp": {"a": 20.0, "b": 40.0}}
    out = proc.update(state, 1.0)
    assert out["alloc__atp"]["a"] == 10.0 and out["alloc__atp"]["b"] == 20.0
    # pool replenished by production this tick (consumers draw it down elsewhere)
    assert out["atp"] == 25.0


def test_metabolism_emits_precursor_supply():
    from viva_mgen.processes.metabolism import MetabolismFbaReproductionProcess
    out = MetabolismFbaReproductionProcess.outputs(
        MetabolismFbaReproductionProcess.__new__(MetabolismFbaReproductionProcess))
    assert "ntp_production" in out and "amino_acid_production" in out


def test_allocator_demand_output_negates_input_no_accumulation(core):
    """Controller ruling: demand__<pool> stores are additive map[float]; without a
    zeroing delta a constant per-tick want would accumulate without bound. The
    allocator must zero exactly what it read, so a re-posted constant demand
    yields the SAME allocation next tick, not a doubled one."""
    proc = AllocatorProcess(config={"pools": ["atp"]}, core=core)
    assert proc.outputs()["demand__atp"] == "map[float]"

    state = {"atp": 100.0, "atp_production": 0.0, "demand__atp": {"a": 10.0}}
    out1 = proc.update(state, 1.0)
    # the allocator's demand__atp output negates exactly what it read
    assert out1["demand__atp"] == {"a": -10.0}
    alloc1 = out1["alloc__atp"]

    # simulate the additive store: zeroing delta + consumer re-posting the same
    # want this tick nets to the SAME value, not a sum of two ticks' wants
    stored_demand = 10.0 + out1["demand__atp"]["a"] + 10.0
    assert stored_demand == 10.0  # no accumulation

    out2 = proc.update({"atp": 100.0, "atp_production": 0.0,
                        "demand__atp": {"a": stored_demand}}, 1.0)
    assert out2["alloc__atp"] == alloc1  # identical allocation, not doubled
    assert out2["demand__atp"] == {"a": -10.0}


def test_composite_builds_and_runs_with_allocator():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    from process_bigraph import Composite
    core = build_core()
    doc = build_mgen(core=core)
    assert "allocator" in doc
    comp = Composite({"state": doc}, core=core)
    comp.run(3.0)
    budget = comp.state["cell"]["budget"]
    assert "alloc__atp" in budget and "alloc__gtp" in budget


# NOTE: the brief's worked example instantiates via
# ``Cls.__new__(Cls); Process.__init__(proc, {}, core=None)`` — the installed
# process_bigraph/bigraph-schema version in this venv now requires a real
# core (raises "must provide a core" for core=None), unrelated to this task's
# changes. These tests build the real process via its normal constructor with
# the module ``core`` fixture instead; the assertions are unchanged from the brief.

def test_atp_consumer_respects_budget(core):
    from viva_mgen.processes.dna import DNASupercoilingReproductionProcess as S
    proc = S(config={"consumer_id": "supercoiling", "initial_sigma": 0.0}, core=core)
    st = {"gyrase": 1000.0, "atp": 1e9,
          "alloc__atp": {"supercoiling": 4.0}}   # budget = 4 ATP -> 2 acts
    out = proc.update(st, 1.0)
    assert out["atp"] >= -4.0            # never spend more than the 4-ATP budget
    assert "demand__atp" in out and out["demand__atp"]["supercoiling"] > 0


def test_dna_repair_respects_atp_budget(core):
    from viva_mgen.processes.dna import DNARepairReproductionProcess as R
    proc = R(config={"consumer_id": "dna_repair"}, core=core)
    atp_per_repair = 4.0  # process default
    budget = 4.0
    # plenty of lesions on the shared per-site map (the sole source of repair
    # demand); enzyme capacity (repair_rate*enzyme*interval = 10) also exceeds
    # the budget, so the ATP budget is the binding constraint.
    st = {"lesion_map": {"3": 1000.0}, "repair_enzyme": 1000.0, "atp": 1e9,
          "alloc__atp": {"dna_repair": budget}}   # budget = 4 ATP, 4 ATP/repair -> 1 repair
    out = proc.update(st, 1.0)
    assert out["atp"] >= -budget            # never spend more than the ATP budget
    assert "demand__atp" in out and out["demand__atp"]["dna_repair"] > 0
    # repair actually happens (map-driven), and is bounded EXACTLY by the budget
    assert out["lesions"] < 0
    assert -out["lesions"] == budget / atp_per_repair
    assert sum(out["lesion_map"].values()) == out["lesions"]


def test_protein_folding_respects_atp_budget(core):
    from viva_mgen.processes.protein import ProteinFoldingReproductionProcess as F
    proc = F(config={"consumer_id": "protein_folding", "seed": 0}, core=core)
    st = {"unfolded": {"g1": 1000.0}, "chaperone_count": 1000.0, "atp": 1e9,
          "alloc__atp": {"protein_folding": 7.0}}   # budget = 7 ATP, 7 ATP/fold -> 1 chap fold
    out = proc.update(st, 1.0)
    assert out["atp"] >= -7.0
    assert "demand__atp" in out and out["demand__atp"]["protein_folding"] > 0


def test_protein_modification_respects_atp_budget(core):
    from viva_mgen.processes.protein import ProteinModificationReproductionProcess as M
    proc = M(config={"consumer_id": "protein_modification", "seed": 0}, core=core)
    st = {"unmodified": {"g1": 1000.0}, "modification_enzyme": 1000.0, "atp": 1e9,
          "alloc__atp": {"protein_modification": 4.0}}   # budget = 4 ATP, 1 ATP/mod -> 4 mods
    out = proc.update(st, 1.0)
    assert out["atp"] >= -4.0
    assert "demand__atp" in out and out["demand__atp"]["protein_modification"] > 0


def test_trna_aminoacylation_respects_atp_budget(core):
    from viva_mgen.processes.rna import TRNAAminoacylationReproductionProcess as T
    proc = T(config={"consumer_id": "trna_aminoacylation", "seed": 0}, core=core)
    st = {"free_trna": {"t1": 1000.0}, "amino_acid": 1e9, "atp": 1e9, "synthetase": 1000.0,
          "alloc__atp": {"trna_aminoacylation": 4.0}}   # budget = 4 ATP, 1 ATP/charge -> 4 charges
    out = proc.update(st, 1.0)
    assert out["atp"] >= -4.0
    assert "demand__atp" in out and out["demand__atp"]["trna_aminoacylation"] > 0


# --- GTP consumers (Task 5) ---
# NOTE: as above, these build via the real constructor + module core fixture
# rather than the brief's Cls.__new__/core=None idiom, which this venv's
# process_bigraph/bigraph-schema now rejects.

def test_translation_respects_gtp_budget(core):
    from viva_mgen.processes.translation import TranslationReproductionProcess as T
    proc = T(config={"consumer_id": "translation", "translation_rates": {"g1": 10.0},
                      "gtp_per_protein": 600.0, "seed": 0}, core=core)
    # unconstrained expectation = rate*copies*interval = 10*1000*1 = 10000 proteins
    # (>> the 1-protein budget allows) so this fails without the budget clamp.
    st = {"rna_counts": {"g1": 1000.0}, "gtp": 1e9,
          "alloc__gtp": {"translation": 600.0}}   # budget = 600 GTP, 600 GTP/protein -> 1 protein
    out = proc.update(st, 1.0)
    assert out["gtp"] >= -600.0
    assert "demand__gtp" in out and out["demand__gtp"]["translation"] > 0


def test_protein_translocation_respects_gtp_budget(core):
    from viva_mgen.processes.protein import ProteinTranslocationReproductionProcess as Tl
    # g1 is a synthetic gene name (not a real KB panel key), so it must be listed
    # explicitly as a translocation substrate for GTP gating to apply to it.
    proc = Tl(config={"consumer_id": "translocation", "seed": 0,
                       "translocated_genes": ["g1"]}, core=core)
    st = {"process_i_done": {"g1": 1000.0}, "translocase": 1000.0, "gtp": 1e9,
          "alloc__gtp": {"translocation": 4.0}}   # budget = 4 GTP, 2 GTP/monomer -> 2 monomers
    out = proc.update(st, 1.0)
    assert out["gtp"] >= -4.0
    assert "demand__gtp" in out and out["demand__gtp"]["translocation"] > 0


def test_ribosome_assembly_respects_gtp_budget(core):
    from viva_mgen.processes.protein import (
        RibosomeAssemblyReproductionProcess as R, _RIBOSOME_30S, _RIBOSOME_50S,
    )
    proc = R(config={"consumer_id": "ribosome_assembly", "seed": 0}, core=core)
    rprot = {p: 1000.0 for p in _RIBOSOME_30S["rproteins"] + _RIBOSOME_50S["rproteins"]}
    rrna_keys = (
        (_RIBOSOME_30S["rrna"] if isinstance(_RIBOSOME_30S["rrna"], list) else [_RIBOSOME_30S["rrna"]])
        + (_RIBOSOME_50S["rrna"] if isinstance(_RIBOSOME_50S["rrna"], list) else [_RIBOSOME_50S["rrna"]])
    )
    rrna = {r: 1000.0 for r in rrna_keys}
    st = {"rprotein_counts": rprot, "rrna_counts": rrna, "assembly_factor": 1000.0,
          "gtp": 1e9, "alloc__gtp": {"ribosome_assembly": 4.0}}   # budget=4 GTP, 2/complex -> 2 complexes
    out = proc.update(st, 1.0)
    assert out["gtp"] >= -4.0
    assert "demand__gtp" in out and out["demand__gtp"]["ribosome_assembly"] > 0


def test_gtp_consumer_respects_budget(core):
    from viva_mgen.processes.cytokinesis import FtsZPolymerizationReproductionProcess as F
    proc = F(config={"consumer_id": "ftsz"}, core=core)
    st = {"ftsz_monomer": 500.0, "gtp": 1e9, "alloc__gtp": {"ftsz": 3.0}}
    out = proc.update(st, 1.0)
    assert out["gtp"] >= -3.0
    assert out["demand__gtp"]["ftsz"] >= 0.0


# --- NTP + amino-acid consumers (Task 6) ---

def test_transcription_respects_ntp_budget(core):
    from viva_mgen.processes.transcription import TranscriptionReproductionProcess as T
    # a strong synthesis rate over a long gene so the unconstrained Poisson
    # expectation (want) is far above the tiny budget, and fails without the clamp.
    proc = T(config={"consumer_id": "transcription", "synthesis_rates": {"g1": 50.0},
                      "gene_lengths": {"g1": 1200.0}, "seed": 0}, core=core)
    st = {"ntp": 1e9, "rna_pol": 100.0, "alloc__ntp": {"transcription": 500.0}}
    out = proc.update(st, 1.0)
    assert -out["ntp"] <= 500.0 + 1e-6
    assert "demand__ntp" in out and out["demand__ntp"]["transcription"] > 500.0


def test_trna_aminoacylation_respects_amino_acid_budget(core):
    from viva_mgen.processes.rna import TRNAAminoacylationReproductionProcess as T
    proc = T(config={"consumer_id": "trna_aminoacylation", "seed": 0}, core=core)
    # amino_acid pool itself is abundant; the allocator's amino_acid budget is
    # the binding constraint (1 AA per charge -> budget caps charges at 4).
    st = {"free_trna": {"t1": 1000.0}, "amino_acid": 1e9, "atp": 1e9, "synthetase": 1000.0,
          "alloc__atp": {"trna_aminoacylation": 1e9},
          "alloc__amino_acid": {"trna_aminoacylation": 4.0}}
    out = proc.update(st, 1.0)
    charged = sum(v for v in out["aminoacylated_trna"].values())
    assert charged <= 4.0 + 1e-6
    assert "demand__amino_acid" in out and out["demand__amino_acid"]["trna_aminoacylation"] > 4.0


def test_trna_aminoacylation_decrements_amino_acid_pool(core):
    """Mass balance: each charged tRNA consumes exactly one free amino acid, so
    the process must emit an "amino_acid" pool delta of -charged (not just cap
    against it), or the pool never drains and the budget stops binding."""
    from viva_mgen.processes.rna import TRNAAminoacylationReproductionProcess as T
    proc = T(config={"consumer_id": "trna_aminoacylation", "seed": 0}, core=core)
    # no alloc keys -> budgets are inf, so nothing but the co-substrate pools caps charging
    st = {"free_trna": {"t1": 1000.0}, "amino_acid": 1e9, "atp": 1e9, "synthetase": 1000.0}
    out = proc.update(st, 1.0)
    charged = sum(v for v in out["aminoacylated_trna"].values())
    assert charged > 0
    assert out["amino_acid"] == -charged


def test_scarcity_partitions_across_consumers():
    """Under a tight shared pool, two consumers scale back PROPORTIONALLY
    rather than one starving the other."""
    from viva_mgen.processes.allocation import allocate
    g = allocate(40.0, {"protein_folding": 100.0, "dna_repair": 300.0})
    assert abs(g["protein_folding"] - 10.0) < 1e-6
    assert abs(g["dna_repair"] - 30.0) < 1e-6
    assert sum(g.values()) <= 40.0 + 1e-9
