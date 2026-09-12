"""Knowledge-base loaders for viva-Mgen.

Two genuine data sources back the reproduction:

* **Metabolic network** — the published Suthers et al. 2009 *M. genitalium*
  genome-scale reconstruction **iPS189** (BioModels ``MODEL1507180052``,
  ``datasets/ips189.sbml.xml``; 351 reactions, 346 metabolites, 126 genes).
  This is the model the Karr 2012 whole-cell model's metabolism submodel was
  built on. Loaded via :mod:`cobra`.
* **Gene list + reference essentiality** — extracted from the Karr 2012 repo's
  ``data/allDeletionSimulations.json`` into ``datasets/genes.csv`` (525 genes,
  each with symbol, name, an essentiality classification, and its associated
  metabolic reactions).

The heavy MATLAB knowledge base (``knowledgeBase.mat``) is an undecodable MCOS
object graph, so per-species copy numbers / kcats are not available from the
repo; the FBA network above stands in for the metabolic parts.
"""

from __future__ import annotations

import csv
import functools
import os
from pathlib import Path
from typing import Optional


def dataset_dir() -> Path:
    """Locate the ``datasets/`` directory (env override, else repo-relative)."""
    env = os.environ.get("VIVA_MGEN_DATASETS")
    if env:
        return Path(env)
    here = Path(__file__).resolve().parent
    for cand in (here.parent / "datasets", here / "datasets"):
        if cand.is_dir():
            return cand
    # last resort: alongside the package
    return here.parent / "datasets"


def dataset_path(name: str) -> Path:
    return dataset_dir() / name


def normalize_gene_id(gid: str) -> str:
    """Normalize a locus tag so the JSON form (``MG_001``) and the SBML form
    (``MG005``) compare equal: strip underscores, uppercase."""
    return gid.replace("_", "").upper()


@functools.lru_cache(maxsize=2)
def load_metabolic_model(sbml_path: Optional[str] = None):
    """Load and return the iPS189 cobra model (cached).

    Raises a clear error if cobra or the SBML file is missing.
    """
    try:
        import cobra
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise RuntimeError(
            "viva-Mgen metabolism requires cobra: `uv pip install cobra`"
        ) from exc

    path = Path(sbml_path) if sbml_path else dataset_path("ips189.sbml.xml")
    if not path.is_file():
        raise FileNotFoundError(
            f"iPS189 SBML not found at {path}. Re-download from BioModels "
            "MODEL1507180052 into datasets/ips189.sbml.xml."
        )
    model = cobra.io.read_sbml_model(str(path))
    return model


@functools.lru_cache(maxsize=1)
def load_genes() -> tuple:
    """Return a tuple of gene dicts:
    ``{gene_id, symbol, name, essential_ref (bool|None), reactions (list[str])}``.
    """
    path = dataset_path("genes.csv")
    genes = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            ess = (row.get("essential_ref") or "").strip().upper()
            genes.append({
                "gene_id": row["gene_id"],
                "symbol": row.get("symbol", ""),
                "name": row.get("name", ""),
                "essential_ref": True if ess == "Y" else (False if ess == "N" else None),
                "reactions": [r.strip() for r in (row.get("associated_reactions") or "").split(",") if r.strip()],
            })
    return tuple(genes)


def gene_essentiality_reference() -> dict:
    """Map normalized gene id -> reference essentiality (bool), where known."""
    out = {}
    for g in load_genes():
        if g["essential_ref"] is not None:
            out[normalize_gene_id(g["gene_id"])] = g["essential_ref"]
    return out


def metabolic_gene_ids() -> list:
    """The gene ids present in the iPS189 metabolic model (SBML form)."""
    return [g.id for g in load_metabolic_model().genes]


@functools.lru_cache(maxsize=1)
def biomass_reaction_id() -> str:
    m = load_metabolic_model()
    for r in m.reactions:
        if r.id.lower() == "biomass" or "biomass" in (r.name or "").lower():
            return r.id
    # fall back to the objective
    return str(next(iter(m.objective.variables)).name)
