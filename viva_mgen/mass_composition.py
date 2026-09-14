"""Emergent dry-mass composition: Σ(species count × molecular weight).

Computes the cell's dry mass and its RNA/protein/DNA/metabolite breakdown from the
real molecular inventory + real molecular weights, for validation against Karr's
fitted dry-weight fractions. See
docs/superpowers/specs/2026-09-14-emergent-mass-design.md.
"""
from __future__ import annotations

from .constants import GENOME_LENGTH_BP, AVOGADRO

_N_A = AVOGADRO        # Avogadro
_RNA_NT_MW = 340.0     # avg ribonucleotide-monophosphate MW in a chain (g/mol)
_AA_MW = 110.0         # avg amino-acid residue MW in a chain (g/mol)
_BP_MW = 660.0         # avg base-pair MW, both strands (g/mol)
_MET_MW = {"atp": 507.0, "gtp": 523.0, "ntp": 500.0, "amino_acid": 110.0}


def _lengths(lengths):
    if lengths is not None:
        return lengths
    from .expression_defaults import gene_lengths
    return gene_lengths()


def rna_mass_g(rna_counts, lengths=None):
    L = _lengths(lengths)
    return sum(float(c) * float(L.get(g, 1000.0)) * _RNA_NT_MW
               for g, c in (rna_counts or {}).items()) / _N_A


def protein_mass_g(protein_counts, lengths=None):
    L = _lengths(lengths)
    return sum(float(c) * (float(L.get(g, 1000.0)) / 3.0) * _AA_MW
               for g, c in (protein_counts or {}).items()) / _N_A


def dna_mass_g(chromosome_copy):
    return float(chromosome_copy) * float(GENOME_LENGTH_BP) * _BP_MW / _N_A


def metabolite_mass_g(pools):
    return sum(float(c) * _MET_MW.get(k, 300.0)
               for k, c in (pools or {}).items()) / _N_A


def emergent_composition(rna_counts, protein_counts, chromosome_copy, pools, lengths=None):
    rna = rna_mass_g(rna_counts, lengths)
    prot = protein_mass_g(protein_counts, lengths)
    dna = dna_mass_g(chromosome_copy)
    met = metabolite_mass_g(pools)
    # "total" is the INCLUSIVE sum (RNA+protein+DNA+metabolite). Callers validating
    # against Karr's macromolecule fractions should sum RNA+protein+DNA themselves
    # (as the mass process does), NOT use "total".
    return {"RNA": rna, "protein": prot, "DNA": dna, "metabolite": met,
            "total": rna + prot + dna + met}
