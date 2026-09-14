#!/usr/bin/env python
"""Regenerate the committed composite-state snapshots the read-only dashboard
serves (reports/composite-state/<id>.json).

The read-only bundle and the loom embed render a composite's wiring + per-process
contracts from this committed artifact, NOT from a live build. So whenever the
process descriptions/contracts or the composite wiring change, this snapshot must
be regenerated or the published loom shows stale text.

Run from a checkout that can build the composites (viva_mgen deps installed):
    PYTHONPATH=<repo> python scripts/regen_composite_state.py
"""
from __future__ import annotations

import json
from pathlib import Path

from vivarium_workbench.lib.composite_resolve import resolve_composite

WS_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = WS_ROOT / "reports" / "composite-state"

COMPOSITES = [
    "viva_mgen.composites.mgen.mycoplasma_genitalium",
    "viva_mgen.composites.mgen.mycoplasma_parca",
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for spec_id in COMPOSITES:
        # remove the stale committed artifact first so resolve_composite does a
        # fresh live build (with current process docstrings) instead of reading
        # the old snapshot back.
        art = OUT_DIR / f"{spec_id}.json"
        if art.exists():
            art.unlink()
        payload = resolve_composite(WS_ROOT, spec_id, allow_build=True)
        if not payload or payload.get("state") is None:
            raise SystemExit(
                f"FAILED to build {spec_id}: {payload.get('notice') if payload else 'no payload'}"
            )
        art.write_text(json.dumps(payload, indent=2))
        n_nodes = sum(1 for _ in _walk_processes(payload["state"]))
        print(f"wrote {art.relative_to(WS_ROOT)}  ({n_nodes} process nodes, wiring={payload.get('wiring_status')})")


def _walk_processes(node):
    if isinstance(node, dict):
        if node.get("_type") in ("process", "step"):
            yield node
        for v in node.values():
            yield from _walk_processes(v)
    elif isinstance(node, list):
        for v in node:
            yield from _walk_processes(v)


if __name__ == "__main__":
    main()
