#!/usr/bin/env python
"""Extract real macromolecular-complex subunit stoichiometry from the Karr 2012
knowledge base (knowledgeBase.mat) into datasets/karr_complexes.json.

The KB stores each ProteinComplex's genuine subunit composition — protein
monomers, sub-complexes, and RNAs each with an integer coefficient. This is the
real structure the whole-cell model assembles (e.g. DNA gyrase = 2 GyrB + 2 GyrA
= A2B2; the 30S ribosomal subunit = 20 r-proteins + 16S rRNA).

Usage: python scripts/extract_kb_complexes.py [path/to/knowledgeBase.mat]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from viva_mgen import kb_decode as k

DEFAULT_MAT = Path.home() / "code/WholeCell-reference/data/knowledgeBase.mat"
OUT = Path(__file__).resolve().parents[1] / "datasets" / "karr_complexes.json"


def main() -> None:
    mat = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_MAT
    raw = k.decode_protein_complexes(mat)

    # Complexes assembled purely from protein monomers → the general
    # MacromolecularComplexation network (Karr forms these from the monomer pool).
    monomer_only = {}
    # Complexes involving RNA or sub-complexes (ribosome, RNA polymerase holoenzyme,
    # etc.) — kept separately with full composition for the assembly processes.
    with_rna_or_subcomplex = {}

    for cid, comp in sorted(raw.items()):
        mono = {str(m): float(n) for m, n in comp["monomers"].items()}
        cxs = {str(c): float(n) for c, n in comp["complexs"].items()}
        rnas = {str(r): float(n) for r, n in comp["rnas"].items()}
        if cxs or rnas:
            entry = {"monomers": mono}
            if cxs:
                entry["complexs"] = cxs
            if rnas:
                entry["rnas"] = rnas
            with_rna_or_subcomplex[cid] = entry
        elif mono:
            monomer_only[cid] = mono

    # Structured ribosome breakdown for RibosomeAssembly (subunit composition of
    # 30S / 50S / 70S as stored in the KB).
    ribosome = {}
    for cid in ("RIBOSOME_30S", "RIBOSOME_50S", "RIBOSOME_70S"):
        if cid in raw:
            comp = raw[cid]
            ribosome[cid] = {
                "monomers": {str(m): float(n) for m, n in comp["monomers"].items()},
                "complexs": {str(c): float(n) for c, n in comp["complexs"].items()},
                "rnas": {str(r): float(n) for r, n in comp["rnas"].items()},
            }

    payload = {
        "source": "Karr et al. 2012 knowledgeBase.mat (ProteinComplex subunit composition)",
        "n_complexes_total": len(raw),
        "monomer_only_complexes": monomer_only,
        "rna_or_subcomplex_complexes": with_rna_or_subcomplex,
        "ribosome": ribosome,
    }
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"wrote {OUT}")
    print(f"  {len(raw)} complexes total")
    print(f"  {len(monomer_only)} monomer-only (general complexation)")
    print(f"  {len(with_rna_or_subcomplex)} with RNA/sub-complex")
    print(f"  ribosome subunits: {list(ribosome)}")


if __name__ == "__main__":
    main()
