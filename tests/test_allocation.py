import math
import numpy as np
from viva_mgen.processes.allocation import allocate, select_budget, demand_entry
from process_bigraph import Process
from viva_mgen.processes.allocation import AllocatorProcess


def test_allocate_surplus_grants_full():
    g = allocate(100.0, {"a": 10.0, "b": 20.0})
    assert g == {"a": 10.0, "b": 20.0}


def test_allocate_scarcity_is_proportional():
    g = allocate(30.0, {"a": 30.0, "b": 60.0})
    assert math.isclose(g["a"], 10.0) and math.isclose(g["b"], 20.0)
    assert sum(g.values()) <= 30.0 + 1e-9


def test_allocate_priority_weights():
    g = allocate(30.0, {"a": 30.0, "b": 30.0}, priorities={"a": 2.0, "b": 1.0})
    assert math.isclose(g["a"], 20.0) and math.isclose(g["b"], 10.0)


def test_allocate_zero_supply_grants_zero():
    g = allocate(0.0, {"a": 5.0})
    assert g == {"a": 0.0}


def test_allocate_ignores_nonpositive_and_conserves():
    g = allocate(10.0, {"a": -3.0, "b": 0.0, "c": 40.0})
    assert g["a"] == 0.0 and g["b"] == 0.0
    assert math.isclose(g["c"], 10.0)


def test_select_budget_absent_is_inf():
    assert select_budget({}, "x") == float("inf")
    assert select_budget({"x": 5.0}, "x") == 5.0


def test_demand_entry_clamps_negative():
    assert demand_entry("x", -2.0) == {"x": 0.0}


def test_allocator_replenishes_and_partitions():
    proc = AllocatorProcess.__new__(AllocatorProcess)   # skip core-requiring __init__
    # Manually set config and pools without calling __init__
    proc.config = {"pools": ["atp"], "priorities": {}, "pool_cap": 1.0e9}
    proc._pools = ["atp"]
    state = {"atp": 5.0, "atp_supply": 25.0,
             "demand__atp": {"a": 20.0, "b": 40.0}}
    out = proc.update(state, 1.0)
    # supply available = level(5) + production(25) = 30; demands 20+40=60 -> proportional
    assert out["alloc__atp"]["a"] == 10.0 and out["alloc__atp"]["b"] == 20.0
    # pool replenished by production this tick (consumers draw it down elsewhere)
    assert out["atp"] == 25.0
