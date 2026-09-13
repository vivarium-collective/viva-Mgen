"""Native decoder for the Karr 2012 knowledge base (``knowledgeBase.mat``).

The knowledge base is a MATLAB MCOS object graph. scipy returns it as an opaque
blob and Octave segfaults reconstructing it, so we decode it directly in Python —
no MATLAB/Octave. Two layers:

1. **Subsystem extraction** — the MAT5 header points to a compressed subsystem
   element (the MCOS ``FileWrapper__``); scipy inflates it, and its
   ``_ObjectMetadata`` is a cell array whose first cell is the MCOS metadata block
   and whose remaining cells hold every property value.

2. **MCOS metadata parse** — the metadata block encodes the class table, the
   per-object class ids, and the per-object property lists (field-name → value).
   Format per Matt Bauman's reverse engineering (as implemented in HebiRobotics/MFL):
   ``version, numStrings, int32[8] segmentIndices``; then ``numStrings``
   null-terminated names; segment 1 = class table; segment 3 = object info
   (classId + segment2/segment4 property-list indices); segments 2 & 4 = the
   property lists (``name, flag, heapIndex`` — flag 1 ⇒ value is cell[heapIndex+2]).

This recovers the genuine per-gene knowledge-base values (expression, half-life,
coordinates, RNA type, Karr's fitted synthesis rate) that the parameter calculator
(:mod:`viva_mgen.parca`) consumes. Verified: MG_001 startCoordinate 686, half-life
2.4 min; 16S rRNA length 1519 nt; the most-expressed genes are the rRNAs.
"""
from __future__ import annotations

import io
import struct
from pathlib import Path

import numpy as np


def _read_subsystem_cells(mat_path: Path):
    """Return the MCOS ``_ObjectMetadata`` value-cell array from a MAT5 file."""
    import scipy.io as sio
    raw = Path(mat_path).read_bytes()
    endian = "<" if raw[126:128] == b"IM" else ">"
    subsys_off = struct.unpack(endian + "q", raw[116:124])[0]
    if subsys_off <= 0:
        raise ValueError("no MCOS subsystem in this .mat")
    src = io.BytesIO(raw)
    rdr = sio.matlab._mio5.MatFile5Reader(src)
    rdr.initialize_read()
    src.seek(subsys_off)
    subsys = np.asarray(rdr.read_var_array(rdr.read_var_header()[0])).tobytes()

    hdr = (b"MATLAB 5.0 subsystem".ljust(116, b" ")
           + (0).to_bytes(8, "little") + (0x0100).to_bytes(2, "little") + b"IM")
    fbuf = io.BytesIO(hdr + subsys[8:])
    frdr = sio.matlab._mio5.MatFile5Reader(fbuf, struct_as_record=True, squeeze_me=False)
    frdr.initialize_read()
    fbuf.seek(128)
    fw = frdr.read_var_array(frdr.read_var_header()[0])
    return np.asarray(fw["MCOS"][0, 0]["_ObjectMetadata"]).ravel()[0].ravel()


def _parse_objects(cells):
    """Parse the MCOS metadata → ``(class_names, objects)`` where each object is
    ``(class_name, {property_name: value})``."""
    meta = np.asarray(cells[0]).ravel().tobytes()
    b = io.BytesIO(meta)
    gi = lambda: struct.unpack("<i", b.read(4))[0]
    _version = gi()
    num_strings = gi()
    seg = [gi() for _ in range(8)]

    strings = []
    for _ in range(num_strings):
        s = bytearray()
        while True:
            c = b.read(1)
            if c in (b"\x00", b""):
                break
            s += c
        strings.append(s.decode("latin1"))
    gs = lambda i: strings[i - 1] if i > 0 else ""

    def value(flag, heap):
        if flag == 0:
            return gs(heap)
        if flag == 1:
            return cells[heap + 2]      # content.get(heapIndex + 2)
        if flag == 2:
            return bool(heap)
        raise ValueError(f"unknown property flag {flag}")

    def parse_props(start, end):
        if start == end:
            return []
        b.seek(start)
        b.read(8)
        out = []
        while b.tell() < end:
            n = gi()
            d = {}
            for _ in range(n):
                name = gs(gi())
                flag = gi()
                heap = gi()
                d[name] = value(flag, heap)
            out.append(d)
            if (n * 3 + 1) % 2 != 0:
                gi()
        return out

    # segment 1 — class table
    b.seek(seg[0])
    b.read(16)
    class_names = []
    while b.tell() < seg[1]:
        _pkg = gs(gi())
        cls = gs(gi())
        b.read(8)
        class_names.append(cls)

    seg2 = parse_props(seg[1], seg[2])
    seg4 = parse_props(seg[3], seg[4])

    # segment 3 — object info
    b.seek(seg[2])
    b.read(24)
    objects = []
    while b.tell() < seg[3]:
        class_id = gi()
        gi(); gi()
        s2 = gi()
        s4 = gi()
        gi()
        props = {}
        if s2 > 0:
            props.update(seg2[s2 - 1])
        if s4 > 0:
            props.update(seg4[s4 - 1])
        cls = class_names[class_id - 1] if 1 <= class_id <= len(class_names) else None
        objects.append((cls, props))
    return class_names, objects


def _scalar(v):
    a = np.asarray(v).ravel()
    return a[0] if a.size else None


def _floats(v):
    return [float(x) for x in np.asarray(v).ravel()]


def load_genome(mat_path: Path) -> str:
    """The full M. genitalium genome sequence (a single ~580 kb ACGT string)."""
    cells = _read_subsystem_cells(mat_path)
    best = ""
    for c in cells:
        a = np.asarray(c)
        if a.dtype.kind in ("U", "S") and a.size == 1:
            s = str(a.ravel()[0])
            if len(s) > len(best) and set(s[:200].upper()) <= set("ACGTUN\n "):
                best = s
    return best


_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# standard genetic code (DNA codons → 1-letter amino acid; * = stop)
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",  # M. genitalium: TGA = Trp
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def gene_dna(gene: dict, genome: str) -> str:
    """The sense-strand DNA sequence of a gene, from the genome + coordinates
    (reverse-complemented for the reverse strand)."""
    s, e, d = gene.get("start"), gene.get("end"), gene.get("direction")
    if not s or not e or e < s:
        return ""
    seq = genome[s - 1:e]                       # 1-based inclusive
    if d is not None and d < 0:
        seq = seq.translate(_COMPLEMENT)[::-1]
    return seq


def rna_base_composition(dna: str) -> dict:
    """RNA nucleotide counts (A/C/G/U) for a transcript from its sense DNA."""
    u = dna.upper()
    return {"A": u.count("A"), "C": u.count("C"), "G": u.count("G"), "U": u.count("T")}


def aa_composition(dna: str) -> dict:
    """Amino-acid counts for a protein by translating its CDS (M. genitalium
    code, TGA=Trp); stops excluded."""
    comp: dict = {}
    for i in range(0, len(dna) - 2, 3):
        aa = _CODON_TABLE.get(dna[i:i + 3].upper())
        if aa and aa != "*":
            comp[aa] = comp.get(aa, 0) + 1
    return comp


def decode_genes(mat_path: Path) -> list:
    """Decode the 525 ``Gene`` objects → a list of dicts with the genuine KB values:
    ``gene_id, symbol, name, rna_type, start, end, length, direction, half_life_min,
    synthesis_rate, expression`` (a 3-list: 32/37/43 C)."""
    cells = _read_subsystem_cells(mat_path)
    _classes, objects = _parse_objects(cells)
    genes = []
    for cls, p in objects:
        if cls != "Gene":
            continue
        start = _scalar(p.get("startCoordinate"))
        end = _scalar(p.get("endCoordinate"))
        hl = _scalar(p.get("halfLife"))
        expr = _floats(p.get("expression")) if p.get("expression") is not None else []
        genes.append({
            "gene_id": _scalar(p.get("wholeCellModelID")),
            "symbol": _scalar(p.get("symbol")) or "",
            "name": _scalar(p.get("name")) or "",
            "rna_type": _scalar(p.get("type")) or "",
            "start": int(start) if start is not None else None,
            "end": int(end) if end is not None else None,
            "length": (int(end) - int(start) + 1) if (start is not None and end is not None) else None,
            "direction": int(_scalar(p.get("direction"))) if _scalar(p.get("direction")) is not None else None,
            "half_life_min": float(hl) if hl is not None else None,
            "synthesis_rate": float(_scalar(p.get("synthesisRate"))) if _scalar(p.get("synthesisRate")) is not None else None,
            "expression": expr,
        })
    return genes
