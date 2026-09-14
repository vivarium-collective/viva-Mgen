"""Whole-cell resource-allocation layer (Karr 2012 hybrid partitioning).

Each finite metabolite pool (ATP, GTP, NTP, amino acids) is partitioned each
tick among the submodels that demand it: demand -> allocate -> run. See
docs/superpowers/specs/2026-09-14-resource-allocation-design.md.
"""
from __future__ import annotations

from process_bigraph import Process

DEFAULT_POOLS = ["atp", "gtp", "ntp", "amino_acid"]


def allocate(supply, demands, priorities=None):
    """Partition ``supply`` (>=0) of one pool among ``demands`` (consumer_id ->
    wanted amount). Non-positive demands get 0. If total weighted demand fits in
    supply every consumer is granted its full demand; otherwise grants are
    proportional to weight*demand. Guarantees ``sum(grants) <= supply``."""
    supply = max(0.0, float(supply))
    pos = {k: float(v) for k, v in (demands or {}).items() if float(v) > 0.0}
    grants = {k: 0.0 for k in (demands or {})}
    if not pos or supply <= 0.0:
        return grants
    total = sum(pos.values())
    if total <= supply:
        grants.update(pos)
        return grants
    w = {k: float((priorities or {}).get(k, 1.0)) for k in pos}
    wsum = sum(w[k] * pos[k] for k in pos)
    if wsum <= 0.0:
        return grants
    for k in pos:
        grants[k] = w[k] * pos[k] / wsum * supply
    return grants


def select_budget(alloc, consumer_id):
    """A consumer's granted budget from an ``alloc__<pool>`` map; ``inf`` when the
    consumer has no entry (tick 0 / unwired) so it behaves as pool-unlimited."""
    if not isinstance(alloc, dict) or consumer_id not in alloc:
        return float("inf")
    return float(alloc[consumer_id])


def demand_entry(consumer_id, want):
    """A ``demand__<pool>`` output entry (non-negative)."""
    return {consumer_id: max(0.0, float(want))}


class AllocatorProcess(Process):
    """Central resource allocator (reproduction of Karr 2012 hybrid partitioning).

    Runs first each tick. For every finite pool it computes the available supply
    (current level + this tick's metabolic production), partitions it among the
    consumers' demands (demand__<pool>) by allocate(), and publishes the grants
    (alloc__<pool>) that each consumer caps its consumption at. Also replenishes
    each pool by this tick's production (consumers draw it down via their own
    negative deltas). The PER-TICK allocation invariant sum(grants) <= supply
    holds exactly; because each consumer's budget is sized one tick ahead of the
    tick it draws against, sustained multi-consumer scarcity can transiently
    drive the pool level slightly below zero (by at most ~one tick's
    over-allocation). This is self-healing within 1-2 ticks — every consumer
    clamps its pool read at max(pool, 0.0) — and the raw pool level is not an
    emitted observable.

    Contract — per pool P: in <P> (level, float), <P>_production (production, float),
    demand__<P> (consumer->want, map). out <P> (replenish delta, float),
    alloc__<P> (consumer->grant, overwrite map).
    Fidelity: FAITHFUL to Karr's demand->allocate->run arbitration (proportional
    partition with a priority-weight hook); one-tick pipeline lag at Δt=1 s.
    """

    description = (
        "Central resource allocator — reproduction of Karr 2012 hybrid partitioning.\n"
        "Runs first each tick: for each finite pool, available supply = pool level + this "
        "tick's metabolic production; partitions it among consumers' demands (demand__<pool>) "
        "by proportional allocation with a priority-weight hook, publishes grants "
        "(alloc__<pool>), and replenishes the pool by production.\n"
        "Contract — per pool P: in <P> (level, float), <P>_production (production, float), "
        "demand__<P> (consumer->want, map). out <P> (replenish delta, float), "
        "alloc__<P> (consumer->grant, overwrite map).\n"
        "Fidelity: FAITHFUL to Karr's demand->allocate->run arbitration (proportional "
        "partition + priority hook); one-tick pipeline lag at dt=1 s."
    )

    config_schema = {
        "pools": {"_type": "list[string]", "_default": DEFAULT_POOLS},
        "priorities": {"_type": "map[map[float]]", "_default": {}},
        "pool_cap": {"_type": "float", "_default": 1.0e9},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._pools = list(self.config["pools"])

    def inputs(self):
        s = {}
        for p in self._pools:
            s[p] = "float"
            s[p + "_production"] = "float"
            s["demand__" + p] = "map[float]"
        return s

    def outputs(self):
        s = {}
        for p in self._pools:
            s[p] = "float"
            s["alloc__" + p] = "overwrite[map[float]]"
            s["demand__" + p] = "map[float]"
        return s

    def initial_state(self):
        return {p + "_production": 0.0 for p in self._pools}

    def update(self, state, interval):
        prio = self.config["priorities"] or {}
        cap = float(self.config["pool_cap"])
        out = {}
        for p in self._pools:
            level = max(0.0, float(state.get(p, 0.0) or 0.0))
            production = max(0.0, float(state.get(p + "_production", 0.0) or 0.0))
            demands = state.get("demand__" + p, {}) or {}
            supply = min(level + production, cap)
            out["alloc__" + p] = allocate(supply, demands, prio.get(p))
            # replenish the pool by this tick's production (bounded by cap)
            out[p] = min(production, max(0.0, cap - level))
            # Zero out exactly what was read this tick, so the additive
            # demand__<pool> store tracks only the latest tick's wants rather
            # than accumulating forever. A consumer's same-tick demand_entry()
            # delta lands alongside this zeroing delta, netting to its latest want.
            out["demand__" + p] = {k: -float(v) for k, v in demands.items()}
        return out
