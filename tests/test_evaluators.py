import json
from pathlib import Path

from viva_mgen.evaluators import _latest_observable, register_derived_scalars, FIG7_FIELDS


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
