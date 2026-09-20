"""Resolve Maritan proteins to a :class:`pbg_parsimony.StructureRef`.

Priority: an S1-curated PDB id (``ProteinRow.pdb_id``) wins; otherwise fall
back to AlphaFold via the protein's UniProt accession, looked up through
``uniprot_index`` (built from S2's per-gene ``uniprot`` / ``protein_monomer``
columns); otherwise ``None``.
"""

from __future__ import annotations

from pbg_parsimony import StructureRef

from viva_mgen.structural.maritan_tables import GeneRow, ProteinRow

# Legacy RCSB ids are 4 characters (e.g. "6RKW"); the newer "extended" PDB
# identifiers introduced in 2023 are longer and are only served as mmCIF.
_LEGACY_PDB_ID_LENGTH = 4


def uniprot_index(
    genes: list[GeneRow], proteins: list[ProteinRow]
) -> dict[str, str]:
    """Map a protein's ``ProtID`` (S1) to its UniProt accession (S2).

    The join key is ``GeneRow.protein_monomer``, which is the exact
    ``ProtID`` used for that gene's monomer in ``mgen_proteins.csv`` (e.g.
    gene ``MG_003`` has ``protein_monomer="MG_003_MONOMER"``, which is
    itself a ``ProteinRow.prot_id``). ``proteins`` is accepted to match the
    brief's interface but is not needed for the join itself.
    """
    del proteins  # not needed for the join; kept for interface parity
    return {
        gene.protein_monomer: gene.uniprot
        for gene in genes
        if gene.protein_monomer and gene.uniprot
    }


# Curated structures for the large assembled machines Maritan hand-curated
# (their §"molecular ingredient modeling" names the 70S ribosome, RNA polymerase
# and GroEL/ES). S1 leaves no single bare PDB id for the assembled particle, and
# a complex has no UniProt, so without these they'd be skipped. Real bacterial
# homolog structures from RCSB (mmCIF for the big assemblies).
_CURATED_STRUCTURES: dict[str, tuple[str, str]] = {
    "RIBOSOME_70S": ("cif", "4YBB"),               # E. coli 70S ribosome
    "RNA_POLYMERASE": ("cif", "4YG2"),             # bacterial RNA polymerase core
    "RNA_POLYMERASE_HOLOENZYME": ("cif", "4YG2"),
    "GROEL_GROES": ("pdb", "1AON"),                # GroEL/ES chaperonin
    "GROEL": ("pdb", "1GRL"),
}


def structure_ref_for(
    protein: ProteinRow, uniprot_by_prot: dict[str, str]
) -> StructureRef | None:
    """Resolve ``protein`` to a structure source, or ``None`` if it has
    neither a curated PDB id nor a resolvable UniProt accession."""
    if protein.pdb_id:
        kind = "pdb" if len(protein.pdb_id) <= _LEGACY_PDB_ID_LENGTH else "cif"
        return StructureRef(kind=kind, ref=protein.pdb_id)
    uniprot = uniprot_by_prot.get(protein.prot_id)
    if uniprot:
        return StructureRef(kind="alphafold", ref=uniprot)
    curated = _CURATED_STRUCTURES.get(protein.prot_id)
    if curated:
        return StructureRef(kind=curated[0], ref=curated[1])
    return None
