"""Whole-cell atom balance (gap #2b): the spent-carrier byproducts of
transcription (PPi) and translation (GDP + Pi) are conserved in pools that
metabolism recycles each tick, rather than being silently dropped."""
from viva_mgen.processes.transcription import TranscriptionReproductionProcess
from viva_mgen.processes.translation import TranslationReproductionProcess
from viva_mgen.core import build_core


def test_transcription_releases_ppi_equal_to_ntp_consumed():
    p = TranscriptionReproductionProcess({"seed": 0}, core=build_core())
    out = p.update({"ntp": 1e9, "rna_pol": 100.0, "alloc__ntp": {}}, 60.0)
    assert out["ppi"] == -out["ntp"]  # one PPi per NTP incorporated (NTP Δ is negative)
    assert out["ppi"] >= 0.0


def test_translation_releases_gdp_and_pi_equal_to_gtp_consumed():
    p = TranslationReproductionProcess({"seed": 1}, core=build_core())
    out = p.update({"rna_counts": {"tuf": 500.0}, "gtp": 1e12, "alloc__gtp": {}}, 60.0)
    assert out["gdp"] == -out["gtp"] and out["pi"] == -out["gtp"]  # GTP -> GDP + Pi
    assert out["gdp"] >= 0.0


def test_byproduct_pools_stay_bounded_not_accumulating():
    # metabolism recycles the byproducts each tick, so the pools hold ~one tick's
    # production over a full 9 h cycle rather than growing ~108x (conservation, not
    # a silent drop or an unbounded pile-up).
    from viva_mgen.composites.mgen import build_mgen
    from process_bigraph import Composite, gather_emitter_results
    core = build_core()
    doc = build_mgen(core, interval=300.0, seed=0)
    doc["metabolism"]["interval"] = 300.0
    for pk in ("transcription", "translation", "rna_decay", "protein_decay", "replication"):
        doc[pk]["interval"] = 300.0
    sim = Composite({"state": doc}, core=core)
    sim.run(3600.0)
    one_tick = None
    sim.run(3600.0 * 8)  # to ~9 h
    m = sim.state["cell"]["metabolism"]
    # after a full cycle the pools are still ~one tick's output, not ~108 ticks
    # (translation GDP per tick is O(1e5-1e6)); assert < 5x a single tick's scale.
    assert 0 <= m["gdp"] < 5e6 and 0 <= m["ppi"] < 5e6 and 0 <= m["pi"] < 5e6


def test_rna_decay_salvages_nmps_to_ntp_pool():
    from viva_mgen.processes.decay import RnaDecayReproductionProcess
    from viva_mgen.expression_defaults import gene_lengths
    p = RnaDecayReproductionProcess({"seed": 3}, core=build_core())
    out = p.update({"rna_counts": {"tuf": 1000.0}}, 300.0)
    decayed = -out["rna_counts"]["tuf"]
    assert out["ntp"] == decayed * gene_lengths().get("tuf", 1000.0)  # NMPs recycled
    assert out["ntp"] >= 0.0


def test_protein_decay_salvages_residues_to_amino_acid_pool():
    from viva_mgen.processes.decay import ProteinDecayReproductionProcess
    from viva_mgen.expression_defaults import gene_lengths
    p = ProteinDecayReproductionProcess({"seed": 4}, core=build_core())
    out = p.update({"protein_counts": {"tuf": 1e5}}, 300.0)
    decayed = -out["protein_counts"]["tuf"]
    assert out["amino_acid"] == decayed * (gene_lengths().get("tuf", 1000.0) / 3.0)
    assert out["amino_acid"] >= 0.0


def test_aminoacylation_releases_amp_and_ppi():
    from viva_mgen.processes.rna import TRNAAminoacylationReproductionProcess
    p = TRNAAminoacylationReproductionProcess({"seed": 5}, core=build_core())
    out = p.update({"free_trna": {f"t{i}": 1000.0 for i in range(5)},
                    "amino_acid": 1e6, "atp": 1e6, "synthetase": 1e4,
                    "alloc__atp": {}, "alloc__amino_acid": {}}, 1.0)
    assert out["amp"] == -out["atp"] and out["ppi"] == -out["atp"]  # ATP -> AMP + PPi
    assert out["amp"] > 0.0
