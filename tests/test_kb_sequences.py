"""Sequence + composition decoding from the knowledge base (genome + coordinates).

Skipped when knowledgeBase.mat is not present (it is external, not vendored).
"""
import os
from pathlib import Path

import pytest

_KB = next((p for p in (
    Path.home() / "code/WholeCell-reference/data/knowledgeBase.mat",
    Path.home() / "Downloads/WholeCell-master/data/knowledgeBase.mat",
) if p.is_file()), None)

pytestmark = pytest.mark.skipif(_KB is None, reason="knowledgeBase.mat not available")


def test_genome_and_gene_sequences():
    from viva_mgen.kb_decode import (aa_composition, decode_genes, gene_dna,
                                     load_genome, rna_base_composition)
    genome = load_genome(_KB)
    assert len(genome) > 500_000 and set(genome[:500].upper()) <= set("ACGTN")

    genes = {g["gene_id"]: g for g in decode_genes(_KB)}
    dnaN = genes["MG_001"]
    dna = gene_dna(dnaN, genome)
    # the sliced CDS length matches the decoded gene length, starts at a start codon
    assert len(dna) == dnaN["length"]
    assert dna[:3] in ("ATG", "GTG", "TTG")

    # RNA base composition accounts for every nucleotide
    rna = rna_base_composition(dna)
    assert sum(rna.values()) == len(dna)

    # translation yields a protein of ~length/3 residues with a single terminal stop
    aa = aa_composition(dna)
    assert abs(sum(aa.values()) - (len(dna) // 3 - 1)) <= 1
