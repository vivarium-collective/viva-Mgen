import numpy as np
import pytest
from viva_mgen.ensemble import run_ensemble, aggregate, final_values


def test_run_ensemble_independent_cells():
    ens = run_ensemble(n_cells=3, duration=5.0)
    assert len(ens["cells"]) == 3 and len(ens["seeds"]) == 3
    assert len(set(ens["seeds"])) == 3
    assert all(len(rows) > 0 for rows in ens["cells"])


def test_cells_diverge():
    ens = run_ensemble(n_cells=3, duration=10.0)
    # total mRNA at the last common step differs across at least two cells (stochastic)
    def total_rna(rows):
        rc = rows[-1].get("rna_counts", {})
        return sum(rc.values()) if isinstance(rc, dict) else 0.0
    totals = [total_rna(c) for c in ens["cells"]]
    assert len(set(totals)) > 1   # not all identical -> seeds genuinely diverge


def test_aggregate_shapes():
    ens = run_ensemble(n_cells=3, duration=10.0)
    agg = aggregate(ens["cells"], "mass")
    n = len(agg["mean"])
    assert n > 0 and len(agg["std"]) == n and len(agg["per_cell"]) == 3
    assert all(s >= 0.0 for s in agg["std"])


def test_final_values():
    ens = run_ensemble(n_cells=3, duration=5.0)
    fv = final_values(ens["cells"], "mass")
    assert len(fv) == 3


def test_missing_observable_raises():
    ens = run_ensemble(n_cells=2, duration=5.0)
    with pytest.raises((KeyError, ValueError)):
        aggregate(ens["cells"], "no_such_observable")


def test_aggregate_skips_missing_step_not_zero_fill():
    # synthetic cells: 3 steps, 3 cells; cell 1 is missing "x" at step 1.
    cells = [
        [{"x": 10.0}, {"x": 10.0}, {"x": 10.0}],
        [{"x": 20.0}, {}, {"x": 20.0}],
        [{"x": 30.0}, {"x": 30.0}, {"x": 30.0}],
    ]
    agg = aggregate(cells, "x")
    # step 0 and 2: mean over all 3 cells -> 20.0
    assert agg["mean"][0] == pytest.approx(20.0)
    assert agg["mean"][2] == pytest.approx(20.0)
    # step 1: cell 1 missing -> mean over remaining cells (10, 30) -> 20.0,
    # NOT dragged toward 0 by zero-filling the missing entry.
    assert agg["mean"][1] == pytest.approx(20.0)
    # the missing entry is preserved as NaN in per_cell, not fabricated as 0.0
    assert np.isnan(agg["per_cell"][1][1])


def test_aggregate_accepts_numpy_scalar():
    cells = [
        [{"x": np.float64(3.0)}, {"x": np.float64(5.0)}],
        [{"x": np.int64(1)}, {"x": np.int64(7)}],
    ]
    agg = aggregate(cells, "x")
    assert agg["n"] == 2
    assert agg["mean"][0] == pytest.approx(2.0)
    assert agg["mean"][1] == pytest.approx(6.0)
