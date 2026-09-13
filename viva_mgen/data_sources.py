"""Workspace data-source provider for the dashboard Resources tab.

Declared in ``workspace.yaml`` as
``dashboard.data_sources.provider: viva_mgen.data_sources:list_sources``. The
workbench calls ``list_sources()`` (no args) in the env worker and renders the
returned rows in the Resources tab, so every model parameter / knowledge-base
file is browsable there.
"""

from __future__ import annotations

from pathlib import Path

# (key, path relative to repo root, category)
_FILES = [
    ("mgen-parameters", "datasets/parameters.csv", "parameters"),
    ("fitted-constants", "datasets/fitted_constants.json", "parameters"),
    ("ips189-network", "datasets/ips189.sbml.xml", "knowledge-base"),
    ("gene-list", "datasets/genes.csv", "knowledge-base"),
]


def list_sources() -> list:
    root = Path(__file__).resolve().parent.parent
    rows = []
    for key, rel, category in _FILES:
        p = root / rel
        rows.append({
            "key": key,
            "path": rel,
            "category": category,
            "kind": "file",
            "size_bytes": (p.stat().st_size if p.exists() else 0),
            "url": "",
        })
    return rows
