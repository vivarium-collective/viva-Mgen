from viva_mgen.provenance import audit_constants, SOURCE_TIERS

def test_audit_covers_processes_and_shape():
    a = audit_constants()
    assert len(a) >= 10
    for proc, consts in a.items():
        for name, rec in consts.items():
            assert set(rec) >= {"value", "source_tier", "kb_match", "note"}
            assert rec["source_tier"] in SOURCE_TIERS

def test_real_kb_classification():
    a = audit_constants()
    # DNASupercoiling gyrase_rate 1.2 is in karr_parameters.json -> real_kb
    assert a["DNASupercoilingReproductionProcess"]["gyrase_rate"]["source_tier"] == "real_kb"
    # ProteinFolding spontaneous_rate has no KB entry -> not real_kb
    assert a["ProteinFoldingReproductionProcess"]["spontaneous_rate"]["source_tier"] != "real_kb"

def test_no_contradiction():
    a = audit_constants()
    for consts in a.values():
        for rec in consts.values():
            if rec["source_tier"] == "real_kb":
                assert rec["kb_match"] is True
            if rec["source_tier"] == "irreducible":
                assert rec["kb_match"] is False
