import re

import pytest

from viva_mgen.structural.maritan_tables import load_proteins, load_genes
from viva_mgen.structural.counts import maritan_counts, mgen_sim_counts, complex_count


def test_maritan_counts_cover_monomers_positive_ints():
    proteins, genes = load_proteins(), load_genes()
    counts = maritan_counts(proteins, genes)
    mono = [p for p in proteins if p.kind == "monomer"]
    covered = [p for p in mono if counts.get(p.prot_id, 0) > 0]
    assert len(covered) > 0.5 * len(mono)          # majority of monomers get a count
    assert all(isinstance(v, int) and v >= 0 for v in counts.values())


def test_complex_count_is_limiting_subunit():
    # "(2.0)MG_003_MONOMER'+(2.0)MG_004_MONOMER'" with subunit stoich 2 each
    mono = {"MG_003_MONOMER": 100, "MG_004_MONOMER": 50}
    bio = "(2.0)MG_003_MONOMER'+(2.0)MG_004_MONOMER'"
    assert complex_count(bio, mono) == 25           # min(100//2, 50//2)


def test_complex_count_zero_when_subunit_missing():
    bio = "(2.0)MG_003_MONOMER'+(2.0)MG_999_MONOMER'"
    assert complex_count(bio, {"MG_003_MONOMER": 100}) == 0


def test_maritan_counts_covers_complexes_including_nested_ones():
    proteins, genes = load_proteins(), load_genes()
    counts = maritan_counts(proteins, genes)
    complexes = [p for p in proteins if p.kind == "complex"]
    assert len(complexes) > 0
    # every complex gets a defined (possibly zero) int count, incl. those
    # whose biosynthesis references other complexes as subunits
    assert all(c.prot_id in counts for c in complexes)
    nested = next(
        c for c in complexes
        if any(not sub.endswith("_MONOMER") for _, sub in
               re.findall(r"\(([\d.]+)\)\s*'?([A-Za-z0-9_]+)'?", c.biosynthesis))
    )
    assert isinstance(counts[nested.prot_id], int)


def test_maritan_counts_scale_sums_near_total_protein_target():
    from viva_mgen.structural.counts import TOTAL_PROTEIN_MOLECULES
    proteins, genes = load_proteins(), load_genes()
    counts = maritan_counts(proteins, genes)
    mono = [p for p in proteins if p.kind == "monomer"]
    monomer_total = sum(counts[p.prot_id] for p in mono)
    # rounding per-gene keeps the sum close to (not exactly) the target
    assert monomer_total == pytest.approx(TOTAL_PROTEIN_MOLECULES, rel=0.05)


def test_mgen_sim_counts_from_synthetic_mapping():
    proteins, genes = load_proteins(), load_genes()
    gene_totals = {"MG_003": 100.0, "MG_004": 50.0}
    counts = mgen_sim_counts(gene_totals, proteins, genes)
    assert counts["MG_003_MONOMER"] == 100
    assert counts["MG_004_MONOMER"] == 50
    gyr = next(p for p in proteins if p.prot_id == "DNA_GYRASE")
    assert counts[gyr.prot_id] == complex_count(gyr.biosynthesis, counts)


def test_mgen_sim_counts_rounds_and_ignores_unmapped_genes():
    proteins, genes = load_proteins(), load_genes()
    gene_totals = {"MG_003": 10.6, "NOT_A_GENE": 999.0}
    counts = mgen_sim_counts(gene_totals, proteins, genes)
    assert counts["MG_003_MONOMER"] == 11
    assert all(k != "NOT_A_GENE" for k in counts)
