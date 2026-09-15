import numpy as np
from viva_mgen.chromosome_state import (
    empty_lesion_map, add_lesions, repair_sites, n_lesions, lesion_positions,
    empty_linking_map, mean_sigma, relax_regions, N_SUPERCOIL_REGIONS,
)


def test_empty_linking_map():
    lm = empty_linking_map(N_SUPERCOIL_REGIONS, -0.06)
    assert len(lm) == N_SUPERCOIL_REGIONS
    assert all(v == -0.06 for v in lm.values())
    assert abs(mean_sigma(lm) - (-0.06)) < 1e-12   # float-sum rounding, not exact
    assert mean_sigma({}) == 0.0


def test_relax_regions_moves_mean_like_single_pool():
    # KEY invariant: distributing `acts` supercoiling acts across R regions moves
    # the MEAN σ by acts/turns_total — exactly as the former single-pool σ did
    # (turns_per_region = turns_total / R). So the superhelical_density observable
    # is unchanged by the per-region migration.
    R, turns_total = 20, 1000.0
    lm = empty_linking_map(R, 0.0)
    setpoint, acts = -0.06, 50.0
    changed = relax_regions(lm, setpoint, acts, turns_total / R, rng=np.random.default_rng(0))
    updated = dict(lm); updated.update(changed)
    # single-pool: dσ = sign(gap)*acts/turns_total
    expected = -1.0 * acts / turns_total
    assert abs(mean_sigma(updated) - expected) < 1e-9
    # never overshoot the setpoint; stays on the correct side
    assert all(setpoint <= v <= 0.0 for v in updated.values())


def test_relax_regions_never_overshoots_setpoint():
    lm = empty_linking_map(4, -0.059)
    # far more acts than needed to reach setpoint -> clamps at setpoint, no overshoot
    changed = relax_regions(lm, -0.06, 1e6, 10.0, rng=np.random.default_rng(1))
    updated = dict(lm); updated.update(changed)
    assert all(abs(v - (-0.06)) < 1e-9 for v in updated.values())

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


def test_chromodynamics_lesions_stall_polymerases():
    # gap #3 fold: ChromosomeDynamics reads the SHARED lesion_map; a lesion on a
    # polymerase's path stalls it and records a *||lesion collision. With NO
    # lesions there is no such collision (baseline non-regression); with lesions
    # blanketing the chromosome, the elongating polymerases hit them.
    from viva_mgen.processes.chromosome import ChromosomeDynamicsReproductionProcess
    from viva_mgen.core import build_core
    core = build_core()
    p = ChromosomeDynamicsReproductionProcess(config={"seed": 0}, core=core)
    clean = p.update({"rna_polymerase": 120.0, "replication_active": 1.0, "lesion_map": {}}, 60.0)
    assert not any("lesion" in k for k in clean["collisions"])  # no lesion collisions at baseline

    p2 = ChromosomeDynamicsReproductionProcess(config={"seed": 0}, core=core)
    lesion_everywhere = {str(b): 1.0 for b in range(0, 580, 3)}  # dense lesions
    out = p2.update({"rna_polymerase": 120.0, "replication_active": 1.0,
                     "lesion_map": lesion_everywhere}, 60.0)
    assert any("lesion" in k for k in out["collisions"])  # polymerases stalled at damaged sites


def test_fork_polymerized_matches_fraction():
    from viva_mgen.chromosome_state import fork_polymerized, empty_polymerized_map
    assert empty_polymerized_map(580) == {str(b): 0.0 for b in range(580)}
    assert fork_polymerized(0.0, 580) == {}
    for frac in (0.1, 0.5, 0.83, 1.0):
        m = fork_polymerized(frac, 580)
        # Σ mask / n_bins reproduces the replicated_fraction (to bin quantization)
        assert abs(sum(m.values()) / 580 - frac) <= 1.0 / 580 + 1e-9
    # bidirectional from oriC(0): both low and high (wrapped) bins fill first
    m = fork_polymerized(0.2, 580)
    assert "0" in m and any(int(b) > 500 for b in m)  # wraps to near-terC on both arms


def test_fork_regions_at_boundaries_only():
    from viva_mgen.chromosome_state import fork_regions, fork_polymerized
    # no replication / complete -> no fork regions (uniform mask)
    assert fork_regions({}, 580, 20) == set()
    assert fork_regions(fork_polymerized(1.0, 580), 580, 20) == set()
    # mid-replication -> a small number of boundary regions (not all 20)
    fr = fork_regions(fork_polymerized(0.4, 580), 580, 20)
    assert 0 < len(fr) < 20
