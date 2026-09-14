"""Single-cell ensembles: run N independent seeded whole-cell composites and
aggregate observables to mean +/- spread + end-of-cycle distributions — the
infrastructure for reproducing Karr 2012's cell-to-cell variation. The submodels
are already stochastic; this exercises them as a population (no model change).
See docs/superpowers/specs/2026-09-14-ensembles-design.md.
"""
from __future__ import annotations

import numbers
import warnings

import numpy as np
from process_bigraph import Composite, gather_emitter_results

from .composites.mgen import build_mgen
from .core import build_core


def run_ensemble(n_cells, duration, *, build_kwargs=None, seeds=None, core=None):
    """Run ``n_cells`` independent composites (one per seed) for ``duration`` s,
    returning ``{"seeds": [...], "cells": [rows_per_cell]}`` where each entry is
    the cell's emitter rows (list of per-step observable dicts)."""
    if core is None:
        core = build_core()
    if seeds is None:
        seeds = list(range(int(n_cells)))
    seeds = [int(s) for s in seeds]
    base = dict(build_kwargs or {})
    cells = []
    for s in seeds:
        kw = dict(base)
        kw["seed"] = s
        doc = build_mgen(core=core, **kw)
        sim = Composite({"state": doc}, core=core)
        sim.run(float(duration))
        rows = gather_emitter_results(sim).get(("emitter",), [])
        cells.append(rows)
    return {"seeds": seeds, "cells": cells}


def _scalar(row, observable):
    v = row.get(observable)
    if isinstance(v, bool):
        return None
    if isinstance(v, numbers.Real):
        return float(v)
    return None


def aggregate(cells, observable):
    """Mean/std across cells at each timepoint for a SCALAR ``observable``.
    Aligns cells to the common (minimum) row count. A step missing the
    observable in a given cell is left as NaN and skipped (nan-aware
    mean/std), not zero-filled, so it doesn't bias the aggregate toward 0.
    Raises if no cell has it at all."""
    series = []
    for rows in cells:
        vals = [_scalar(r, observable) for r in rows]
        if any(v is not None for v in vals):
            series.append([v if v is not None else float("nan") for v in vals])
    if not series:
        raise ValueError(f"observable {observable!r} not found (or non-scalar) in any cell")
    n_steps = min(len(s) for s in series)
    if n_steps == 0:
        raise ValueError(f"no emitter rows for observable {observable!r}")
    mat = np.array([s[:n_steps] for s in series], dtype=float)  # (cells, steps), NaN = missing
    # nan-aware reduction over cells; a timestep with all-NaN yields NaN (documented)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN slice
        mean = np.nanmean(mat, axis=0)
        std = np.nanstd(mat, axis=0)
    return {
        "observable": observable,
        "n": len(series),
        "t_index": list(range(n_steps)),
        "mean": mean.tolist(),
        "std": std.tolist(),
        "per_cell": mat.tolist(),
    }


def final_values(cells, observable):
    """Last value of a scalar ``observable`` per cell (end-of-cycle distribution)."""
    out = []
    for rows in cells:
        val = None
        for r in reversed(rows):
            v = _scalar(r, observable)
            if v is not None:
                val = v
                break
        out.append(val)
    return out
