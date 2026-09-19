"""Workspace derived-scalar computers for viva-mGen study behavior tests.

Study run scripts (sims/run.py) compute figure observables in-process and
persist them to <workspace>/.pbg/runs.jsonl via append_run_event(..., observables=...).
These computers read those already-computed values back so the study's behavior
tests grade through the framework's derived-scalar seam. Scoped to fig7 for now;
add more fields here as other studies' computers are authored.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_RUN_LOG_RELPATH = ".pbg/runs.jsonl"

# fig7-kinetic-parameters observables (unique field names across all studies)
FIG7_FIELDS = (
    "growth_is_monotonic_in_kcat",
    "growth_saturates_at_wt",
    "growth_is_sigmoidal",
    "growth_dynamic_range",
    "growth_at_max_kcat",
    "growth_at_min_kcat",
    "kcat_half_max",
)


def _latest_observable(ws_root: Any, field: str) -> float:
    """Return `field` from the most-recently-completed run whose observables carry it.

    Reads <ws_root>/.pbg/runs.jsonl (append-only JSONL). Raises ValueError if no
    completed run recorded this observable — the evaluator surfaces that honestly
    rather than fabricating a value.
    """
    log = Path(ws_root) / _RUN_LOG_RELPATH
    if not log.is_file():
        raise ValueError(f"no run log at {log}")
    best_t = None
    best_v = None
    with log.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                ev = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not isinstance(ev, dict):
                continue
            obs = ev.get("observables")
            if not isinstance(obs, dict) or field not in obs:
                continue
            t = ev.get("completed_at") or 0.0
            if best_t is None or t >= best_t:
                best_t = t
                best_v = obs[field]
    if best_v is None:
        raise ValueError(f"observable {field!r} not found in any completed run")
    try:
        return float(best_v)
    except (TypeError, ValueError):
        raise ValueError(f"observable {field!r} is not numeric: {best_v!r}")


def register_derived_scalars(reg: dict) -> None:
    """Register fig7 derived-scalar computers (field -> fn(reader, test, ws_root))."""
    for field in FIG7_FIELDS:
        reg[field] = lambda reader, test, ws_root, _f=field: _latest_observable(ws_root, _f)
