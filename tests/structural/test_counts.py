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


def test_maritan_counts_total_protein_is_maritan_plausible():
    # Subunit-conserving assembly: free monomers + assembled complexes should
    # total in the Maritan Table-1 ballpark (~26k), not the ~10x-inflated figure
    # the old independent-complex counting produced. (The volume is conserved, so
    # test_protein_occupancy_matches_maritan pins the 0.144 fraction separately.)
    proteins, genes = load_proteins(), load_genes()
    counts = maritan_counts(proteins, genes)
    total = sum(v for v in counts.values() if v > 0)
    assert 18_000 <= total <= 36_000        # Maritan ~26k
    # ribosomes assemble (rRNA subunits non-limiting)
    assert counts.get("RIBOSOME_70S", 0) > 0


def test_mgen_sim_counts_from_synthetic_mapping():
    proteins, genes = load_proteins(), load_genes()
    gene_totals = {"MG_003": 100.0, "MG_004": 50.0}
    counts = mgen_sim_counts(gene_totals, proteins, genes)
    # DNA_GYRASE = (2)MG_003 + (2)MG_004 -> min(100//2, 50//2) = 25, consuming
    # 50 of each. Subunit-conserving assembly leaves the FREE remainder:
    assert counts["DNA_GYRASE"] == 25
    assert counts["MG_003_MONOMER"] == 50   # 100 synthesized - 50 consumed
    assert counts["MG_004_MONOMER"] == 0    # 50 synthesized - 50 consumed


def test_mgen_sim_counts_rounds_and_ignores_unmapped_genes():
    proteins, genes = load_proteins(), load_genes()
    gene_totals = {"MG_003": 10.6, "NOT_A_GENE": 999.0}
    counts = mgen_sim_counts(gene_totals, proteins, genes)
    assert counts["MG_003_MONOMER"] == 11
    assert all(k != "NOT_A_GENE" for k in counts)
