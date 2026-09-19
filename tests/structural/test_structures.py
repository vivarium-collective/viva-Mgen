from pbg_parsimony import StructureRef
from viva_mgen.structural.maritan_tables import load_proteins, load_genes
from viva_mgen.structural.structures import structure_ref_for, uniprot_index


def test_pdb_takes_priority():
    proteins, genes = load_proteins(), load_genes()
    idx = uniprot_index(genes, proteins)
    gyr = next(r for r in proteins if r.prot_id == "DNA_GYRASE")
    ref = structure_ref_for(gyr, idx)
    assert ref.kind in ("pdb", "cif") and ref.ref == "6RKW"


def test_alphafold_fallback_when_no_pdb_but_uniprot():
    from viva_mgen.structural.maritan_tables import ProteinRow
    p = ProteinRow(prot_id="MG_XXX_MONOMER", name="x", function="f",
                   compartment="c", kind="monomer", pdb_id=None,
                   biosynthesis="", seq_length=100)
    ref = structure_ref_for(p, {"MG_XXX_MONOMER": "P47000"})
    assert ref.kind == "alphafold" and ref.ref == "P47000"


def test_none_when_no_structure_source():
    from viva_mgen.structural.maritan_tables import ProteinRow
    p = ProteinRow("NO_STRUCT", "x", "f", "c", "monomer", None, "", None)
    assert structure_ref_for(p, {}) is None
