"""Readers for the Maritan 2022 supplement, parsed once by
``parse_supplement.py`` into committed CSV data artifacts.

These readers only touch ``data/mgen_proteins.csv`` and ``data/mgen_genes.csv``
at runtime (the source .xlsx files are not re-parsed here).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent / "data"
_PROTEINS_CSV = _DATA_DIR / "mgen_proteins.csv"
_GENES_CSV = _DATA_DIR / "mgen_genes.csv"


@dataclass(frozen=True)
class ProteinRow:
    prot_id: str
    name: str
    function: str
    compartment: str
    kind: str  # "monomer" | "complex"
    pdb_id: str | None
    biosynthesis: str
    seq_length: int | None
    dna_binding: str | None = None    # "dsDNA" | "ssDNA" — a nucleoid-bound protein
    dna_footprint: int | None = None  # bp footprint on the DNA (positional)


@dataclass(frozen=True)
class GeneRow:
    gene_id: str
    gtype: str
    coord: int
    length: int
    direction: str
    tu: str
    essential: bool
    name: str
    uniprot: str | None
    protein_monomer: str | None


def _int_or_none(raw: str) -> int | None:
    return int(raw) if raw not in (None, "") else None


def _str_or_none(raw: str) -> str | None:
    return raw if raw not in (None, "") else None


def load_proteins() -> list[ProteinRow]:
    with _PROTEINS_CSV.open(newline="") as f:
        reader = csv.DictReader(f)
        return [
            ProteinRow(
                prot_id=row["prot_id"],
                name=row["name"],
                function=row["function"],
                compartment=row["compartment"],
                kind=row["kind"],
                pdb_id=_str_or_none(row["pdb_id"]),
                biosynthesis=row["biosynthesis"],
                seq_length=_int_or_none(row["seq_length"]),
                dna_binding=_str_or_none(row.get("dna_binding", "")),
                dna_footprint=_int_or_none(row.get("dna_footprint", "")),
            )
            for row in reader
        ]


def load_genes() -> list[GeneRow]:
    with _GENES_CSV.open(newline="") as f:
        reader = csv.DictReader(f)
        return [
            GeneRow(
                gene_id=row["gene_id"],
                gtype=row["gtype"],
                coord=_int_or_none(row["coord"]),
                length=_int_or_none(row["length"]),
                direction=row["direction"],
                tu=row["tu"],
                essential=row["essential"] == "True",
                name=row["name"],
                uniprot=_str_or_none(row["uniprot"]),
                protein_monomer=_str_or_none(row["protein_monomer"]),
            )
            for row in reader
        ]
