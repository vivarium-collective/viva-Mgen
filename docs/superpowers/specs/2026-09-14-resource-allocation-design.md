# Resource-allocation layer for viva-Mgen (closing gap #1)

Date: 2026-09-14
Status: approved design, pre-implementation

## Problem

The Karr 2012 *M. genitalium* whole-cell model integrates its 28 submodels with a
**hybrid allocation algorithm**: every 1 s tick each submodel declares its
*demand* for the shared molecule pools (ATP, GTP, NTPs, dNTPs, amino acids, …),
each pool is *partitioned* among the demanding submodels, and each submodel then
runs on its *allocated* share. This demand → allocate → run arbitration is what
makes the model "whole-cell": scarcity in one pathway propagates to every process
that shares a resource.

viva-Mgen currently has **no allocation layer**. Two concrete deficiencies:

1. **The shared pools are effectively infinite.** `cell/metabolism/{atp,gtp,ntp,
   amino_acid}` are initialized at 1e8–1e9 in `_FLOAT_INIT` and never replenished
   or drawn to a limiting level. Metabolism's `atp_production` / `gtp_production`
   are separate readout floats, *not* wired back into the pools. Nothing is ever
   scarce.
2. **Consumers draw opportunistically.** Each pool-consuming process caps its own
   consumption at the raw pool level, in whatever order the composite runs them.
   There is no arbitration: the first process to run can drain a pool the others
   also needed.

Closing gap #1 therefore requires two coupled changes: make the pools a **finite
per-tick budget replenished by metabolism**, and **partition that budget by
demand** each tick.

## Chosen approach

**Central allocator + demand/budget protocol** (selected over a full
Requester/Evolver process split, and over a light post-hoc throttle).

A single `AllocatorProcess` mediates the finite pools. Each consuming process
gains a `demand` output and a `budget` input; its internal biology is unchanged
except that its consumption cap becomes `min(existing kinetic limit, budget)`.

### Timing: one-tick pipeline

process-bigraph calls each `Process.update()` once per tick, so a process cannot
both *declare* its demand and *receive* its allocation within a single call
without being split into two processes (the Requester/Evolver architecture we
ruled out as too heavy). The central-allocator form pipelines the two phases
across consecutive ticks:

- **Tick t, consumer:** consume ≤ its current `budget` (set at t−1), and emit its
  `demand` for tick t (a cheap recomputation of the same unconstrained draw its
  kinetics already imply).
- **Tick t, allocator (scheduled first):** read every consumer's `demand` (from
  t) and each pool's available supply, compute allocations, write each consumer's
  `budget` for tick t+1.

At Δt = 1 s against the ~9 h (32400 s) cell cycle the one-tick lag is
biologically negligible, and the pipeline is a genuine demand → allocate → run
arbitration. Zero-lag within a single tick is only achievable by splitting every
consumer into request/execute processes; that is explicitly out of scope for v1.

### Finite pools replenished by metabolism

Metabolism's FBA output becomes real per-tick supply:

- `atp_production` (ATPS4r flux) → ATP molecules/tick into the ATP pool.
- `gtp_production` (summed GTP-linked flux) → GTP molecules/tick into the GTP pool.
- NTP / dNTP / amino_acid supplied from `growth_fraction ×` the ParCa's fitted
  precursor demand (`datasets/karr_metabolic_demand.json`, already computed).

Pools initialize at a realistic finite level (order-of-magnitude the steady-state
copy number, not 1e9) and are drawn down by consumers and replenished by
metabolism each tick, so scarcity emerges from the dynamics rather than being
hard-coded.

### Allocation math

For each partitioned pool: let `S` = available supply this tick (current pool
level + this tick's metabolic production) and `d_i` = consumer *i*'s demand.

- If `Σ d_i ≤ S`: grant each consumer its full `d_i` (surplus remains in the pool).
- Else: grant proportionally, `g_i = d_i · S / Σ d_i`.

Proportional allocation with an optional per-consumer **priority weight**
(`g_i = w_i d_i · S / Σ w_j d_j`), defaulting to equal weights for v1. The
priority hook is present so Karr's process-priority ordering can be added later
without a redesign. Invariant: `Σ g_i ≤ S` always (conservation).

### Pools and consumers

| Pool | Consumers wired to demand/budget |
|------|----------------------------------|
| ATP  | tRNA-aminoacylation, protein folding, protein modification, DNA repair, DNA supercoiling |
| GTP  | translation, translocation, ribosome assembly, FtsZ polymerization |
| NTP  | transcription |
| dNTP | replication |
| amino_acid | translation, tRNA-aminoacylation |

ATP and GTP are the multi-consumer pools where arbitration truly bites;
single-consumer pools (NTP, dNTP) still receive a finite budget so their pathway
becomes metabolism-limited rather than infinite.

## Components

### `AllocatorProcess` (new — `viva_mgen/processes/allocation.py`)

- **inputs:** the finite pool levels (`atp`, `gtp`, `ntp`, `dntp`, `amino_acid`
  as floats), the per-pool production floats from metabolism (`atp_production`,
  `gtp_production`, …), and the per-pool demand maps
  (`demand__atp`, `demand__gtp`, … : `map[float]` keyed by consumer name).
- **outputs:** the per-pool allocation/budget maps (`alloc__atp`, … :
  `map[float]` keyed by consumer name), and the replenished pool deltas (adds this
  tick's production into each pool).
- **config:** `pools` (list), `priorities` (`map[map[float]]`: pool → consumer →
  weight, default equal), `initial_pool` (finite starting levels).
- **update:** for each pool, `S = level + production`; read `demand__<pool>`;
  compute `alloc__<pool>` by the proportional (+priority) rule; emit the pool
  replenishment delta. Pure and unit-testable.

### Consumer process changes (existing files)

Each consumer in the table above gains:

- a `consumer_id` config value (its allocation key; defaults to a stable name for
  that process).
- an input port for each pool it uses: the full `alloc__<pool>` map (`map[float]`,
  consumer → grant). The consumer selects its own grant as
  `budget = alloc.get(consumer_id, ∞)`.
- an output port for each pool it uses: it writes `{consumer_id: unconstrained_draw}`
  into `demand__<pool>` for the next tick's allocation.
- a one-line change: the consumption cap becomes `min(existing_kinetic_limit,
  budget)`.

**Absent-allocation default:** when `alloc__<pool>` carries no entry for a
consumer (tick 0 before the allocator has run, or a consumer left unwired), its
budget defaults to the full pool (`∞` cap) — so the very first tick is never
falsely starved and an unwired consumer behaves exactly as it does today. This
keeps the change back-compatible and the allocator strictly *additional*
constraint, never a new failure mode.

Reading/writing the whole per-pool map (rather than wiring each consumer's scalar
entry) keeps the store shape uniform and lets the allocator see every demand in
one place — the natural fit for process-bigraph's `map[float]` stores.

No change to any process's kinetic constants or core mechanism.

### Composite (`viva_mgen/composites/mgen.py`)

- New `budget` store group in `_STORE_GROUP` for the `demand__*` / `alloc__*`
  maps; `_MAP_STORES` entries so they initialize to `{}`.
- `AllocatorProcess` added to `_NODE_NAMES`, scheduled to run first each tick.
- `_FLOAT_INIT` pool sizes reduced from 1e8–1e9 to realistic finite values.
- Metabolism outputs wired so `atp_production` / `gtp_production` (and the
  precursor supply) reach the allocator.

## Data flow

```
metabolism ──atp_production/gtp_production/growth_fraction──▶ AllocatorProcess
                                                                │
   cell/metabolism/{atp,gtp,ntp,dntp,amino_acid} ◀─replenish───┤
                                                                │
consumers ──demand__<pool> (map: consumer→want)──▶  AllocatorProcess
consumers ◀──alloc__<pool>  (map: consumer→grant)──  AllocatorProcess
   │
   └─ consume min(kinetic_limit, budget); emit next-tick demand
```

## Error handling / edge cases

- **Zero supply:** all grants 0; consumers make no progress that tick (correct —
  starvation propagates).
- **No demand for a pool:** allocator grants 0, pool accumulates production up to
  an optional cap (avoid unbounded growth of an unused pool).
- **Unknown consumer key in a demand map:** ignored by the allocator (defensive).
- **Missing metabolic production (FBA infeasible, `feasible = 0`):** production 0
  for that tick; pools deplete — this is the correct coupling of a metabolic
  failure to the rest of the cell.
- **Conservation guard:** `Σ grants ≤ S` (the per-tick allocation invariant) asserted
  in the allocator, and holds exactly. It does *not* mean the pool level itself
  stays ≥ 0 across ticks: because a consumer's budget for tick t is sized at
  t−1 against the t−1 supply, sustained multi-consumer scarcity can transiently
  drive the pool slightly below zero (by at most ~one tick's over-allocation).
  This is self-healing within 1-2 ticks since every consumer clamps its pool
  read at `max(pool, 0.0)`, and the raw pool level is not an emitted observable.

## Testing

Unit (`tests/test_allocation.py`):
- surplus → full grant; scarcity → proportional; priority weights honored;
  zero-supply → zero; conservation (`Σ g_i ≤ S`); unknown-key ignored.

Integration (extend `tests/`):
- Build the full composite; drive one pool scarce (high demand / low supply) and
  assert consumers scale back **proportionally** (no single consumer starves the
  rest), which the pre-allocation opportunistic draw would not do.
- Assert mass/growth still tracks under allocation (no regression in the growth
  phenotype).
- `feasible = 0` propagates to pool depletion.

Post-merge:
- Regenerate `reports/composite-state/*.json` (`scripts/regen_composite_state.py`)
  so the read-only loom shows the new `AllocatorProcess` node + `budget` stores,
  and republish.

## Out of scope (v1)

- Within-tick zero-lag allocation (needs the Requester/Evolver split).
- Karr's exact process-priority ordering (the weight hook is present; default
  equal).
- Partitioning non-energy/precursor molecules (water, Pi, PPi) — those remain
  delegated to the pools; this is tracked separately under gap #2 (mass
  conservation).

## Fidelity impact

Moves the whole-cell **integration algorithm** from "opportunistic draw on
infinite pools" to a genuine demand → allocate → run arbitration on finite,
metabolism-replenished pools — the core of Karr's hybrid algorithm. Each affected
process's `description` fidelity line will be updated to state that its
consumption is now allocation-arbitrated.
