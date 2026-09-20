from viva_mgen.structural.maritan_tables import load_proteins, load_genes


def test_proteins_roster_size_and_known_row():
    rows = load_proteins()
    assert len(rows) > 400                      # ~996 species in S1
    gyr = next(r for r in rows if r.prot_id == "DNA_GYRASE")
    assert gyr.kind == "complex"
    assert gyr.pdb_id == "6RKW"                 # S1 "Structural Model"
    assert gyr.compartment == "c"


def test_genes_known_row_and_uniprot_join():
    genes = load_genes()
    trna = next(g for g in genes if g.gene_id == "MG471")
    assert trna.gtype.lower() == "trna"
    assert trna.essential is True
    # at least some protein-coding genes carry a UniProt accession
    assert any(g.uniprot for g in genes)
