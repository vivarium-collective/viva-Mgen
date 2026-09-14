from viva_mgen.kb import load_protein_maturation, maturation_by_panel_key

def test_counts():
    m = load_protein_maturation()
    assert sum(1 for v in m.values() if v["signal_type"] == "lipoprotein") == 14
    assert sum(1 for v in m.values() if v["signal_type"] == "secretory") == 20
    assert sum(1 for v in m.values() if v["met_cleavage"]) == 35
    assert m["MG_067"]["signal_type"] == "lipoprotein"

def test_panel_key_map():
    pk = maturation_by_panel_key()
    assert len(pk) > 100 and any(v["signal_type"] == "lipoprotein" for v in pk.values())


def test_processing_ii_lipoprotein_throttled():
    from viva_mgen.processes.protein import ProteinProcessingIIReproductionProcess as P2
    from viva_mgen.core import build_core
    p = P2(config={"lipoprotein_genes": ["lipoA"]}, core=build_core())
    # a lipoprotein and a non-lipoprotein both present in large amounts
    out = p.update({"translocated": {"lipoA": 1000.0, "plainB": 1000.0},
                    "signal_peptidase": 100.0, "diacylglyceryl_transferase": 1.0}, 1.0)
    moved = out["processed_ii"]
    # lipoA limited by transferase (1 * 0.0165) -> ~0; plainB limited only by peptidase -> more
    assert moved.get("plainB", 0) > moved.get("lipoA", 0)


def test_processing_ii_default_lipoprotein_genes_from_kb():
    from viva_mgen.processes.protein import ProteinProcessingIIReproductionProcess as P2
    from viva_mgen.core import build_core
    # empty config -> real KB lipoprotein classification (14 genes)
    p = P2(config={}, core=build_core())
    assert len(p._lipoprotein_genes) == 14


def test_processing_i_met_cleavage_throttled():
    from viva_mgen.processes.protein import ProteinProcessingIReproductionProcess as P1
    from viva_mgen.core import build_core
    p = P1(config={"met_cleavage_genes": ["metA"]}, core=build_core())
    out = p.update({"nascent": {"metA": 1000.0, "otherB": 1000.0},
                    "deformylase": 100.0, "aminopeptidase": 1.0}, 1.0)
    moved = out["process_i_done"]
    # otherB deformylated fully by the shared deformylase limit; metA additionally
    # gated by the near-zero aminopeptidase limit (1 * 6.0) -> far less
    assert moved.get("otherB", 0) > moved.get("metA", 0)
    # every protein still advances -- nothing is dropped from the pipeline
    assert moved.get("otherB", 0) > 0


def test_translocation_only_substrates():
    from viva_mgen.processes.protein import ProteinTranslocationReproductionProcess as T
    from viva_mgen.core import build_core
    p = T(config={"translocated_genes": ["secA"]}, core=build_core())
    out = p.update({"process_i_done": {"secA": 100.0, "cytoC": 100.0},
                    "translocase": 1e12, "gtp": 1e9}, 1.0)
    # both reach 'translocated' (pipeline doesn't stall) but only secA consumed GTP
    assert out["gtp"] < 0     # GTP consumed for secA
    assert out["translocated"].get("secA", 0) > 0
    assert out["translocated"].get("cytoC", 0) == 100.0  # cytoC passes through in full, no gating
