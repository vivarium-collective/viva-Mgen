#!/usr/bin/env python
"""Pack the reproduction-driven 3D M. genitalium cell (s02): same Maritan
roster + structures + geometry as s01, but abundances from a viva_mgen
whole-cell SIMULATION rather than the WC-MG expression proxy.

Runs the ``mycoplasma_genitalium`` composite for a short sim, reads its evolved
per-gene ``protein_counts`` store, maps those to Maritan S1 monomer ProtIDs
(the store is keyed by gene SYMBOL / mixed id, so it is normalized to MG_00N
gene ids first), packs the cell with the real parsimony engine, and compares
occupancy + roster to the s01 baseline.

Usage:
    PARSIMONY_HOME=/path/to/parsimony python sims/run.py [--sim-seconds S] [--top-n N]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

from process_bigraph import Composite

from viva_mgen.core import build_core
from viva_mgen.composites import build_mgen
from viva_mgen.structural.maritan_tables import load_proteins, load_genes
from viva_mgen.structural.counts import maritan_counts, mgen_sim_counts
from viva_mgen.structural.build import build_mgen_pack, mgen_ingredients, estimate_occupancy
from viva_mgen.structural.viewer import relativize_pack_urls

STUDY_DIR = Path(__file__).resolve().parents[1]
# The workbench's built-in Parsimony Viewer discovers packs under <study>/viz/3d/.
PACK_DIR = STUDY_DIR / "viz" / "3d"
_DATASETS = Path(__file__).resolve().parents[4] / "datasets"


def _symbol_to_gene_id(genes) -> dict[str, str]:
    """Map gene SYMBOL -> MG_00N gene id, from karr_gene_expression.csv (gene_id,
    symbol) plus the S2 gene table (gene_id, name). Used to normalize the
    simulation's symbol-keyed ``protein_counts`` to the gene ids mgen_sim_counts
    joins on."""
    m: dict[str, str] = {}
    csv_path = _DATASETS / "karr_gene_expression.csv"
    if csv_path.exists():
        for row in csv.DictReader(open(csv_path)):
            sym, gid = (row.get("symbol") or "").strip(), (row.get("gene_id") or "").strip()
            if sym and gid:
                m[sym] = gid
    for g in genes:
        if g.name and g.gene_id:
            m.setdefault(g.name.strip(), g.gene_id)
    return m


def _normalize_counts(protein_counts: dict, genes) -> dict[str, float]:
    """Turn the sim's ``{symbol|id|"a, b": count}`` store into ``{MG_00N: count}``."""
    valid_ids = {g.gene_id for g in genes}
    sym2id = _symbol_to_gene_id(genes)
    out: dict[str, float] = {}
    for key, val in protein_counts.items():
        if not isinstance(val, (int, float)) or val <= 0:
            continue
        gid = None
        if key in valid_ids:
            gid = key
        else:
            for tok in str(key).replace(";", ",").split(","):
                tok = tok.strip()
                if tok in valid_ids:
                    gid = tok
                    break
                if tok in sym2id:
                    gid = sym2id[tok]
                    break
        if gid:
            out[gid] = out.get(gid, 0.0) + float(val)
    return out


def _find_protein_counts(state) -> dict:
    """cell/proteome/protein_counts, located defensively."""
    stack = [state]
    while stack:
        node = stack.pop()
        if isinstance(node, dict):
            pc = node.get("protein_counts")
            if isinstance(pc, dict) and pc:
                return pc
            stack.extend(v for v in node.values() if isinstance(v, dict))
    return {}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim-seconds", type=float, default=600.0,
                    help="viva_mgen sim time to evolve the proteome (default 600 s)")
    ap.add_argument("--top-n", type=int, default=None)
    args = ap.parse_args()

    if not (os.environ.get("PARSIMONY_HOME") or os.environ.get("PARSIMONY_BIN")):
        print("ERROR: set PARSIMONY_HOME (or PARSIMONY_BIN).", file=sys.stderr)
        return 2

    proteins, genes = load_proteins(), load_genes()

    # 1. Evolve the reproduction's proteome with a short whole-cell simulation.
    print(f"Running mycoplasma_genitalium composite for {args.sim_seconds:.0f}s sim ...")
    core = build_core()
    sim = Composite({"state": build_mgen(core, interval=1.0)}, core=core)
    sim.run(args.sim_seconds)
    raw = _find_protein_counts(sim.state)
    gene_counts = _normalize_counts(raw, genes)
    print(f"  protein_counts: {len(raw)} store keys -> {len(gene_counts)} MG gene ids "
          f"(sum {sum(gene_counts.values()):.0f})")

    # 2. Sim-driven copy numbers over the Maritan roster.
    counts = mgen_sim_counts(gene_counts, proteins, genes)

    # 3. Pack + read out.
    PACK_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Packing reproduction-driven cell (top_n={args.top_n}) ...")
    res = build_mgen_pack(counts, out_dir=PACK_DIR, top_n=args.top_n,
                          name="mgen_reproduction_driven")
    relativize_pack_urls(res["pack_path"])  # for the workbench built-in viewer

    ings = mgen_ingredients(counts, top_n=args.top_n)
    occ = estimate_occupancy(counts, top_n=args.top_n)
    base = maritan_counts(proteins, genes)
    base_species = sum(1 for p in proteins if base.get(p.prot_id, 0) > 0)
    sim_species = sum(1 for p in proteins if counts.get(p.prot_id, 0) > 0)

    checks = {
        "structural-sim-cell-comparable-to-baseline":
            int(res.get("n_placed", 0)) > 0 and sim_species > 0.5 * base_species,
    }
    summary = {
        "n_placed": int(res.get("n_placed", 0)),
        "n_species": len(ings),
        "sim_species_with_counts": sim_species,
        "baseline_species_with_counts": base_species,
        "protein_occupancy": round(occ["occupancy"], 3),
        "pack_path": str(res.get("pack_path", "")),
        "checks": checks,
    }
    (PACK_DIR / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    ok = all(checks.values())
    for name, passed in checks.items():
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
