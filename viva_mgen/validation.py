"""Karr-2012 comparison helpers: load reference targets and extract/score the
emergent observables the study report cards gate on. See
docs/superpowers/specs/2026-09-14-karr-validation-design.md.
"""
from __future__ import annotations

import functools
import json

from .kb import dataset_path

_MACRO = ("protein", "DNA", "RNA")


@functools.lru_cache(maxsize=1)
def _reference_values() -> dict:
    path = dataset_path("karr_reference_values.json")
    with open(path) as f:
        return json.load(f).get("values", {})


def reference(vid: str) -> dict:
    vals = _reference_values()
    if vid not in vals:
        raise KeyError(f"no Karr reference value {vid!r}")
    return vals[vid]


def emergent_macro_fractions(row) -> dict:
    """Renormalize a row's emergent_mass_fractions over the macromolecules
    {protein, DNA, RNA} (drop metabolite/other); all-zero → zeros."""
    emf = (row or {}).get("emergent_mass_fractions", {}) or {}
    vals = {k: float(emf.get(k, 0.0)) for k in _MACRO}
    tot = sum(vals.values())
    if tot <= 0:
        return {k: 0.0 for k in _MACRO}
    return {k: vals[k] / tot for k in _MACRO}


def mrna_cv(final_totals) -> float:
    """Population coefficient of variation of per-cell final total mRNA."""
    xs = [float(x) for x in (final_totals or []) if x is not None]
    if len(xs) < 2:
        return 0.0
    mean = sum(xs) / len(xs)
    if mean == 0:
        return 0.0
    var = sum((x - mean) ** 2 for x in xs) / len(xs)
    return (var ** 0.5) / mean


def score(observed, ref) -> dict:
    lo, hi = ref["band"]
    target = float(ref["target"])
    passed = (observed is not None) and (lo <= float(observed) <= hi)
    if observed is None:
        closeness = 0.0
    else:
        denom = max(abs(target), 1e-9)
        closeness = max(0.0, 1.0 - abs(float(observed) - target) / denom)
    return {"observed": (None if observed is None else float(observed)),
            "target": target, "band": [lo, hi], "pass": bool(passed),
            "closeness": float(closeness)}
