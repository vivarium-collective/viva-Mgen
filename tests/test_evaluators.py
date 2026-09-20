import json
from pathlib import Path

from viva_mgen.evaluators import _latest_observable, register_derived_scalars, ALL_FIELDS, FIG7_FIELDS


def _write_log(ws: Path, events):
    (ws / ".pbg").mkdir(parents=True, exist_ok=True)
    with (ws / ".pbg" / "runs.jsonl").open("w", encoding="utf-8") as fh:
        for ev in events:
            fh.write(json.dumps(ev) + "\n")


def test_reads_latest_observable(tmp_path):
    _write_log(tmp_path, [
        {"run_id": "old", "completed_at": 1.0, "observables": {"growth_dynamic_range": 0.5}},
        {"run_id": "new", "completed_at": 2.0, "observables": {"growth_dynamic_range": 0.9}},
    ])
    assert _latest_observable(tmp_path, "growth_dynamic_range") == 0.9


def test_missing_observable_raises(tmp_path):
    _write_log(tmp_path, [{"run_id": "r", "completed_at": 1.0, "observables": {"x": 1.0}}])
    try:
        _latest_observable(tmp_path, "growth_is_sigmoidal")
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_registers_all_fig7_fields():
    reg = {}
    register_derived_scalars(reg)
    for f in FIG7_FIELDS:
        assert f in reg and callable(reg[f])


def test_registered_fn_reads_value(tmp_path):
    _write_log(tmp_path, [{"run_id": "r", "completed_at": 1.0,
                           "observables": {"growth_is_monotonic_in_kcat": 1.0}}])
    reg = {}
    register_derived_scalars(reg)
    assert reg["growth_is_monotonic_in_kcat"](None, {}, tmp_path) == 1.0


def test_registers_all_studies_fields():
    reg = {}
    register_derived_scalars(reg)
    # every declared field across all 8 studies is registered
    for f in ALL_FIELDS:
        assert f in reg and callable(reg[f])
    # field names are unique across studies (no collision in the union)
    assert len(ALL_FIELDS) == len(set(ALL_FIELDS))


def test_cross_study_field_reads_own_value(tmp_path):
    # a non-fig7 field (fig4) resolves from its own study's run, bound correctly
    _write_log(tmp_path, [
        {"run_id": "fig4-cell-cycle-baseline", "completed_at": 5.0,
         "observables": {"cell_cycle_h": 12.5, "essentiality_accuracy": 0.9}},
    ])
    reg = {}
    register_derived_scalars(reg)
    assert reg["cell_cycle_h"](None, {}, tmp_path) == 12.5
    assert reg["essentiality_accuracy"](None, {}, tmp_path) == 0.9
