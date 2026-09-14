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
