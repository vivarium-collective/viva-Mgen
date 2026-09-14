import numpy as np
from viva_mgen.chromosome_state import (
    empty_lesion_map, add_lesions, repair_sites, n_lesions, lesion_positions,
)

def test_empty_lesion_map():
    m = empty_lesion_map(10)
    assert len(m) == 10 and set(m.values()) == {0.0} and all(isinstance(k, str) for k in m)

def test_add_lesions_total_weight():
    rng = np.random.default_rng(0)
    d = add_lesions(rng, n_bins=100, count=5)
    assert sum(d.values()) == 5.0
    assert all(0 <= int(b) < 100 for b in d)

def test_add_lesions_zero():
    assert add_lesions(np.random.default_rng(0), 100, 0) == {}

def test_repair_never_exceeds_present():
    m = {"3": 2.0, "7": 1.0}   # 3 lesions present
    delta, repaired = repair_sites(m, capacity=10, rng=np.random.default_rng(1))
    assert repaired == 3.0
    # applying delta zeroes the map, never negative
    for b, dv in delta.items():
        assert m[b] + dv >= -1e-9
    assert sum(delta.values()) == -3.0

def test_repair_partial():
    m = {"3": 5.0}
    delta, repaired = repair_sites(m, capacity=2, rng=np.random.default_rng(2))
    assert repaired == 2.0 and delta == {"3": -2.0}

def test_n_lesions():
    assert n_lesions({"1": 2.0, "2": 0.0, "5": 3.0}) == 5.0

def test_lesion_positions():
    pos = lesion_positions({"0": 1.0, "50": 2.0, "3": 0.0}, genome_length_bp=1000.0, n_bins=100)
    assert 0.0 in pos and 500.0 in pos and len(pos) == 2   # bin 3 has 0 count -> excluded


def test_damage_populates_lesion_map():
    from viva_mgen.processes.dna import DNADamageReproductionProcess
    from viva_mgen.core import build_core
    p = DNADamageReproductionProcess(config={"agent_rate": 1.0, "n_bins": 100}, core=build_core())
    out = p.update({"damaging_agent": 50.0}, 1.0)
    from viva_mgen.chromosome_state import n_lesions
    assert n_lesions(out["lesion_map"]) == out["lesions"] and out["lesions"] > 0

def test_damage_zero_agent_no_lesions():
    from viva_mgen.processes.dna import DNADamageReproductionProcess
    from viva_mgen.core import build_core
    p = DNADamageReproductionProcess(config={"base_rate": 0.0, "n_bins": 100}, core=build_core())
    out = p.update({"damaging_agent": 0.0}, 1.0)
    assert out["lesions"] == 0.0 and out["lesion_map"] == {}

def test_repair_clears_sites():
    from viva_mgen.processes.dna import DNARepairReproductionProcess
    from viva_mgen.core import build_core
    p = DNARepairReproductionProcess(config={"repair_rate": 1.0}, core=build_core())
    out = p.update({"lesion_map": {"3": 2.0, "9": 1.0}, "repair_enzyme": 100.0, "atp": 1e9}, 1.0)
    assert out["lesions"] < 0 and sum(out["lesion_map"].values()) == out["lesions"]

def test_composite_baseline_no_lesions():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    from process_bigraph import Composite
    core = build_core()
    comp = Composite({"state": build_mgen(core=core)}, core=core)
    comp.run(3.0)
    lm = comp.state["cell"]["genome"]["lesion_map"]
    assert sum(lm.values()) == 0.0   # damaging_agent defaults 0 -> no lesions (non-regression)
