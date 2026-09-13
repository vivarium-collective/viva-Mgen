#!/usr/bin/env python3
"""Decode the Karr 2012 knowledge base (``knowledgeBase.mat``) and export the
genuine per-gene parameters to ``datasets/karr_gene_expression.csv``.

No MATLAB/Octave — the decode is done natively by :mod:`viva_mgen.kb_decode`,
which parses the MCOS object graph directly (see that module for the format). For
each of the 525 genes it recovers the real knowledge-base values: RNA type, the
genome coordinates (→ length), transcription direction, the mRNA half-life, Karr's
fitted synthesis rate, and the Weiner et al. 2003 observed expression profile
(32 / 37 / 43 C). Row order follows ``datasets/genes.csv`` (joined by gene id).

Validation printed on run: mRNA half-lives average a few minutes while tRNA/rRNA
are stable; 16S rRNA length ≈ 1519 nt; the most-expressed genes are the rRNAs.

Usage:
    python scripts/extract_kb_genes.py [path/to/knowledgeBase.mat]

The .mat itself is NOT vendored (4 MB, external); point this at a WholeCell
checkout's data/knowledgeBase.mat.
"""
from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[1]
_DEFAULT_KB = [
    Path.home() / "code/WholeCell-reference/data/knowledgeBase.mat",
    Path.home() / "Downloads/WholeCell-master/data/knowledgeBase.mat",
]


def main() -> int:
    sys.path.insert(0, str(WS))
    from viva_mgen.kb_decode import decode_genes, metabolic_demand

    kb = Path(sys.argv[1]) if len(sys.argv) > 1 else next(
        (p for p in _DEFAULT_KB if p.is_file()), None)
    if not kb or not kb.is_file():
        print("knowledgeBase.mat not found; pass its path as an argument.\n"
              f"looked in: {', '.join(map(str, _DEFAULT_KB))}", file=sys.stderr)
        return 2

    decoded = {g["gene_id"]: g for g in decode_genes(kb)}
    genes = list(csv.DictReader((WS / "datasets/genes.csv").open()))

    out = WS / "datasets/karr_gene_expression.csv"
    n = 0
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene_id", "symbol", "name", "rna_type", "start_coordinate",
                    "end_coordinate", "length_nt", "direction", "half_life_min",
                    "synthesis_rate", "expression_32C", "expression_37C",
                    "expression_43C", "expression_mean"])
        for g in genes:                                  # genes.csv order
            d = decoded.get(g["gene_id"])
            if not d:
                continue
            e = d["expression"] or [None, None, None]
            emean = (sum(e) / len(e)) if e and all(x is not None for x in e) else ""
            w.writerow([
                d["gene_id"], d["symbol"] or g.get("symbol", ""),
                d["name"] or g.get("name", ""), d["rna_type"],
                d["start"] if d["start"] is not None else "",
                d["end"] if d["end"] is not None else "",
                d["length"] if d["length"] is not None else "",
                d["direction"] if d["direction"] is not None else "",
                f"{d['half_life_min']:.6g}" if d["half_life_min"] is not None else "",
                f"{d['synthesis_rate']:.6g}" if d["synthesis_rate"] is not None else "",
                f"{e[0]:.6g}" if e[0] is not None else "",
                f"{e[1]:.6g}" if e[1] is not None else "",
                f"{e[2]:.6g}" if e[2] is not None else "",
                f"{emean:.6g}" if emean != "" else "",
            ])
            n += 1
    print(f"wrote {out} ({n} genes) from {kb}")

    # aggregate metabolic demand (NMP + AA composition the expression implies) —
    # the ParCa's forward coupling to metabolism.
    import json
    demand = metabolic_demand(kb)
    dpath = WS / "datasets/karr_metabolic_demand.json"
    json.dump(demand, dpath.open("w"), indent=1)
    print(f"wrote {dpath} (NMP A+U={demand['nmp']['A'] + demand['nmp']['U']:.2f})")

    # provenance sanity
    vals = list(decoded.values())
    top = sorted(vals, key=lambda d: -(sum(d["expression"]) if d["expression"] else 0))[:3]
    print("top-3 expressed:", [d["gene_id"] for d in top])
    for rt in ("mRNA", "tRNA", "rRNA"):
        hls = [d["half_life_min"] for d in vals if d["rna_type"] == rt and d["half_life_min"]]
        if hls:
            print(f"{rt}: n={len(hls)} half-life mean {statistics.mean(hls):.1f} min")
    rrna = next((d for d in vals if "rrn" in d["gene_id"].lower() and d["rna_type"] == "rRNA"), None)
    if rrna:
        print(f"{rrna['gene_id']} length {rrna['length']} nt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
