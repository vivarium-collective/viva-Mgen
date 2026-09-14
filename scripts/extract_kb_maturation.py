#!/usr/bin/env python
"""Extract per-protein maturation classification from knowledgeBase.mat into
datasets/karr_protein_maturation.json. Usage: python scripts/extract_kb_maturation.py [mat]"""
from __future__ import annotations
import json, sys
from pathlib import Path
from viva_mgen import kb_decode as k

DEFAULT_MAT = Path.home() / "code/WholeCell-reference/data/knowledgeBase.mat"
OUT = Path(__file__).resolve().parents[1] / "datasets" / "karr_protein_maturation.json"

def main():
    mat = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_MAT
    m = k.decode_protein_maturation(mat)
    from collections import Counter
    st = Counter(v["signal_type"] for v in m.values())
    payload = {"source": "Karr 2012 knowledgeBase.mat (ProteinMonomer signal/Met classification)",
               "n_proteins": len(m),
               "n_lipoprotein": sum(1 for v in m.values() if v["signal_type"] == "lipoprotein"),
               "n_secretory": sum(1 for v in m.values() if v["signal_type"] == "secretory"),
               "n_met_cleavage": sum(1 for v in m.values() if v["met_cleavage"]),
               "maturation": {g: v for g, v in sorted(m.items())}}
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"wrote {OUT}: {len(m)} proteins, signal_type={dict(st)}")

if __name__ == "__main__":
    main()
