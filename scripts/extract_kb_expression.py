#!/usr/bin/env python3
"""Decode the Karr 2012 whole-cell knowledge base (``knowledgeBase.mat``) and
export the per-gene observed expression matrix to a plain CSV.

The knowledge base is a MATLAB MCOS object graph that scipy returns as an opaque
blob and that Octave segfaults trying to reconstruct. Rather than run MATLAB, we
decode it *directly*: the MAT5 file header points to a compressed "subsystem"
data element (the MCOS FileWrapper) whose ``_ObjectMetadata`` cell array holds
every property value as an ordinary array. The gene expression matrix is the one
``(nGenes, 3)`` float array in that store — the three columns are the observed
transcription profiles of Weiner et al. 2003 (M. pneumoniae, grown at 32 / 37 /
43 C), mapped onto the M. genitalium genes. Row order matches ``datasets/genes.csv``
(verified: the top rows are 5S/23S/16S rRNA then the tRNAs — the expression peaks).

This is the observed-expression input the Karr parameter fit (FitConstants) uses.

Usage:
    python scripts/extract_kb_expression.py [path/to/knowledgeBase.mat]

Writes datasets/karr_gene_expression.csv. The .mat itself is NOT vendored (4 MB,
external); point this at a WholeCell checkout's data/knowledgeBase.mat.
"""
from __future__ import annotations

import csv
import io
import struct
import sys
import warnings
from pathlib import Path

import numpy as np
import scipy.io as sio

WS = Path(__file__).resolve().parents[1]
_DEFAULT_KB = [
    Path.home() / "code/WholeCell-reference/data/knowledgeBase.mat",
    Path.home() / "Downloads/WholeCell-master/data/knowledgeBase.mat",
]


def _read_subsystem(mat_path: Path) -> bytes:
    """Return the decompressed MCOS subsystem stream from a MAT5 file."""
    raw = mat_path.read_bytes()
    endian = "<" if raw[126:128] == b"IM" else ">"
    subsys_off = struct.unpack(endian + "q", raw[116:124])[0]
    if subsys_off <= 0:
        raise ValueError("no MCOS subsystem in this .mat")
    buf = io.BytesIO(raw)
    rdr = sio.matlab._mio5.MatFile5Reader(buf)
    rdr.initialize_read()
    buf.seek(subsys_off)
    hdr, _ = rdr.read_var_header()
    arr = rdr.read_var_array(hdr)  # scipy transparently inflates the compressed element
    return np.asarray(arr).tobytes()


def _object_metadata_cells(subsys: bytes):
    """Parse the subsystem stream and return the MCOS _ObjectMetadata cell list."""
    # the subsystem is a header-less MAT5 stream (8-byte marker, then elements);
    # wrap it in a synthetic 128-byte header so the standard reader can walk it.
    hdr = (b"MATLAB 5.0 subsystem".ljust(116, b" ")
           + (0).to_bytes(8, "little") + (0x0100).to_bytes(2, "little") + b"IM")
    buf = io.BytesIO(hdr + subsys[8:])
    rdr = sio.matlab._mio5.MatFile5Reader(buf, struct_as_record=True, squeeze_me=False)
    rdr.initialize_read()
    buf.seek(128)
    fw = rdr.read_var_array(rdr.read_var_header()[0])
    mcos = fw["MCOS"][0, 0]
    return np.asarray(mcos["_ObjectMetadata"]).ravel()[0].ravel()


def _find_expression(cells) -> np.ndarray:
    """The gene expression matrix is the sole (nGenes, 3) float64 cell."""
    cands = []
    for c in cells:
        xa = np.asarray(c)
        if xa.dtype == np.float64 and xa.ndim == 2 and xa.shape[1] == 3 and xa.shape[0] > 400:
            cands.append(xa)
    if len(cands) != 1:
        raise ValueError(f"expected exactly one (nGenes,3) float matrix, found {len(cands)}")
    return cands[0]


def main() -> int:
    warnings.filterwarnings("ignore")
    if len(sys.argv) > 1:
        kb = Path(sys.argv[1])
    else:
        kb = next((p for p in _DEFAULT_KB if p.is_file()), None)
    if not kb or not kb.is_file():
        print("knowledgeBase.mat not found; pass its path as an argument.\n"
              f"looked in: {', '.join(str(p) for p in _DEFAULT_KB)}", file=sys.stderr)
        return 2

    expr = _find_expression(_object_metadata_cells(_read_subsystem(kb)))
    genes = list(csv.DictReader((WS / "datasets/genes.csv").open()))
    if len(genes) != expr.shape[0]:
        raise ValueError(f"gene count mismatch: genes.csv={len(genes)} expr={expr.shape[0]}")

    out = WS / "datasets/karr_gene_expression.csv"
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene_id", "symbol", "name",
                    "expression_32C", "expression_37C", "expression_43C", "expression_mean"])
        for g, row in zip(genes, expr):
            w.writerow([g["gene_id"], g.get("symbol", ""), g.get("name", ""),
                        f"{row[0]:.6g}", f"{row[1]:.6g}", f"{row[2]:.6g}", f"{row.mean():.6g}"])
    print(f"wrote {out} ({expr.shape[0]} genes) from {kb}")
    # provenance sanity: the top-expressed genes should be rRNA
    top = np.argsort(expr.mean(1))[::-1][:3]
    print("top-3 expressed:", [genes[int(i)]["gene_id"] for i in top])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
