"""One-shot generator: parse the Maritan 2022 supplementary S1/S2 spreadsheets
into committed CSV data artifacts under ``viva_mgen/structural/data/``.

openpyxl is not installed in this environment, so this reads the .xlsx files
directly as zip archives of OOXML: the shared-string table
(``xl/sharedStrings.xml``) plus the first worksheet (``xl/worksheets/sheet1.xml``).

Run once via ``python -m viva_mgen.structural.parse_supplement`` (or
``python viva_mgen/structural/parse_supplement.py``) from the repo root; the
resulting CSVs are committed so that runtime code (``maritan_tables.py``)
never touches the source .xlsx files.
"""

from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path
from typing import Any

import xml.etree.ElementTree as ET

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

_HERE = Path(__file__).resolve().parent
_SUPP_DIR = (
    _HERE.parent.parent
    / "workspace"
    / "references"
    / "papers"
    / "Maritan2022_supplementary"
)
_DATA_DIR = _HERE / "data"

_S1_PATH = _SUPP_DIR / "S1_MG_proteins_ingredients.xlsx"
_S2_PATH = _SUPP_DIR / "S2_MG_genes.xlsx"
_S3_PATH = _SUPP_DIR / "S3_membrane_protein_assignment.xlsx"
_S3_DATA_START = 4          # header at row 3 (Monomer Id, Consensus, …)
_S3_COL_MONOMER = 0
_S3_COL_CONSENSUS = 1

# S1 column indices (0-based), per the curated header row (row index 3).
_S1_HEADER_ROW = 3
_S1_DATA_START = 4
_S1_COL_PROT_ID = 0
_S1_COL_NAME = 1
_S1_COL_FUNCTION = 2
_S1_COL_COMPARTMENT = 3
_S1_COL_TYPE = 4
_S1_COL_SEQ_LENGTH = 6
_S1_COL_BIOSYNTHESIS = 8
_S1_COL_DNA_FOOTPRINT = 9   # bp footprint of a DNA-binding protein (positional)
_S1_COL_DNA_BINDING = 10    # "dsDNA" | "ssDNA" — marks nucleoid-bound proteins
_S1_COL_STRUCTURAL_MODEL = 11  # curated "Structural Model" (first occurrence)

# S2 column indices (0-based), per the header row (row index 2).
_S2_HEADER_ROW = 2
_S2_DATA_START = 3
_S2_COL_GENE_ID = 0
_S2_COL_TYPE = 1
_S2_COL_COORD = 2
_S2_COL_LENGTH = 3
_S2_COL_DIRECTION = 4
_S2_COL_TU = 5
_S2_COL_ESSENTIAL = 6
_S2_COL_NAME = 7
_S2_COL_UNIPROT = 8
_S2_COL_PROTEIN_MONOMER = 9

_PDB_ID_RE = re.compile(r"^[0-9A-Za-z]{4}$")


def _col_to_idx(ref: str) -> int:
    """Convert an Excel cell reference's column letters (e.g. 'AB12') to a
    0-based column index."""
    letters = re.match(r"[A-Z]+", ref).group()
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - ord("A") + 1)
    return idx - 1


def _load_shared_strings(z: zipfile.ZipFile) -> list[str]:
    try:
        data = z.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(data)
    strings = []
    for si in root.findall(f"{_NS}si"):
        # Concatenate all <t> runs (handles rich-text cells split across runs).
        texts = si.findall(f".//{_NS}t")
        strings.append("".join(t.text or "" for t in texts))
    return strings


def _load_sheet_rows(
    z: zipfile.ZipFile, sheet: str, shared: list[str]
) -> list[list[Any]]:
    """Return a list of rows (each a list of cell values, padded with None up
    to the row's last populated column). Row order matches sheet order."""
    root = ET.fromstring(z.read(sheet))
    sheet_data = root.find(f"{_NS}sheetData")
    rows: list[list[Any]] = []
    for row in sheet_data.findall(f"{_NS}row"):
        cells: dict[int, Any] = {}
        for c in row.findall(f"{_NS}c"):
            ref = c.get("r")
            t = c.get("t")
            v = c.find(f"{_NS}v")
            val = v.text if v is not None else None
            if t == "s" and val is not None:
                val = shared[int(val)]
            cells[_col_to_idx(ref)] = val
        if cells:
            max_c = max(cells.keys())
            rows.append([cells.get(i) for i in range(max_c + 1)])
        else:
            rows.append([])
    return rows


def _cell(row: list[Any], idx: int) -> Any:
    return row[idx] if idx < len(row) else None


def _read_workbook_rows(path: Path) -> list[list[Any]]:
    with zipfile.ZipFile(path) as z:
        shared = _load_shared_strings(z)
        return _load_sheet_rows(z, "xl/worksheets/sheet1.xml", shared)


def _normalize_kind(raw: str | None) -> str:
    if raw and "complex" in raw.lower():
        return "complex"
    return "monomer"


def _normalize_pdb_id(raw: str | None) -> str | None:
    if raw and _PDB_ID_RE.match(raw):
        return raw
    return None


def _int_footprint(raw) -> int | None:
    """DNA footprint (bp) as an int, or None when absent/unparseable."""
    if raw in (None, ""):
        return None
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return None


def _s3_consensus() -> dict[str, str]:
    """S3 consensus membrane-compartment assignments (Monomer Id -> c|m|e).

    S3 resolves the ~31 proteins where the cytoplasmic-MG and WC-MG models (and
    the localization predictors: BUSCA, SignalP, SOSUI, Phobius, PSORTb)
    disagreed. We use its consensus to refine S1's compartment (it reclassifies
    a few cytoplasmic proteins as membrane)."""
    if not _S3_PATH.is_file():
        return {}
    out = {}
    for row in _read_workbook_rows(_S3_PATH)[_S3_DATA_START:]:
        mid = _cell(row, _S3_COL_MONOMER)
        cons = (_cell(row, _S3_COL_CONSENSUS) or "").strip().lower()
        if mid and cons in ("c", "m", "e"):
            out[mid] = cons
    return out


def parse_proteins() -> list[dict[str, Any]]:
    rows = _read_workbook_rows(_S1_PATH)
    s3 = _s3_consensus()
    out = []
    for row in rows[_S1_DATA_START:]:
        prot_id = _cell(row, _S1_COL_PROT_ID)
        if not prot_id:
            continue
        seq_length_raw = _cell(row, _S1_COL_SEQ_LENGTH)
        seq_length = None
        if seq_length_raw not in (None, ""):
            try:
                seq_length = int(float(seq_length_raw))
            except ValueError:
                seq_length = None
        compartment = (_cell(row, _S1_COL_COMPARTMENT) or "").lower()
        compartment = s3.get(prot_id, compartment)   # S3 consensus refines S1
        out.append(
            {
                "prot_id": prot_id,
                "name": _cell(row, _S1_COL_NAME) or "",
                "function": _cell(row, _S1_COL_FUNCTION) or "",
                "compartment": compartment,
                "kind": _normalize_kind(_cell(row, _S1_COL_TYPE)),
                "pdb_id": _normalize_pdb_id(_cell(row, _S1_COL_STRUCTURAL_MODEL)),
                "biosynthesis": _cell(row, _S1_COL_BIOSYNTHESIS) or "",
                "seq_length": seq_length,
                "dna_binding": (_cell(row, _S1_COL_DNA_BINDING) or "").strip(),
                "dna_footprint": _int_footprint(_cell(row, _S1_COL_DNA_FOOTPRINT)),
            }
        )
    return out


def parse_genes() -> list[dict[str, Any]]:
    rows = _read_workbook_rows(_S2_PATH)
    out = []
    for row in rows[_S2_DATA_START:]:
        gene_id = _cell(row, _S2_COL_GENE_ID)
        if not gene_id:
            continue
        coord_raw = _cell(row, _S2_COL_COORD)
        length_raw = _cell(row, _S2_COL_LENGTH)
        essential_raw = _cell(row, _S2_COL_ESSENTIAL)
        out.append(
            {
                "gene_id": gene_id,
                "gtype": _cell(row, _S2_COL_TYPE) or "",
                "coord": int(float(coord_raw)) if coord_raw not in (None, "") else None,
                "length": int(float(length_raw)) if length_raw not in (None, "") else None,
                "direction": _cell(row, _S2_COL_DIRECTION) or "",
                "tu": _cell(row, _S2_COL_TU) or "",
                "essential": (essential_raw or "").strip().lower() == "yes",
                "name": _cell(row, _S2_COL_NAME) or "",
                "uniprot": _cell(row, _S2_COL_UNIPROT) or None,
                "protein_monomer": _cell(row, _S2_COL_PROTEIN_MONOMER) or None,
            }
        )
    return out


_PROTEIN_FIELDS = [
    "prot_id",
    "name",
    "function",
    "compartment",
    "kind",
    "pdb_id",
    "biosynthesis",
    "seq_length",
    "dna_binding",
    "dna_footprint",
]

_GENE_FIELDS = [
    "gene_id",
    "gtype",
    "coord",
    "length",
    "direction",
    "tu",
    "essential",
    "name",
    "uniprot",
    "protein_monomer",
]


def _write_csv(path: Path, fieldnames: list[str], records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


def main() -> None:
    proteins = parse_proteins()
    genes = parse_genes()
    _write_csv(_DATA_DIR / "mgen_proteins.csv", _PROTEIN_FIELDS, proteins)
    _write_csv(_DATA_DIR / "mgen_genes.csv", _GENE_FIELDS, genes)
    print(f"wrote {len(proteins)} protein rows -> {_DATA_DIR / 'mgen_proteins.csv'}")
    print(f"wrote {len(genes)} gene rows -> {_DATA_DIR / 'mgen_genes.csv'}")


if __name__ == "__main__":
    main()
