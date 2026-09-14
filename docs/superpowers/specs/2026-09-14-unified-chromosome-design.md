# Unified per-site chromosome for viva-Mgen (gap #3, phase 1)

Date: 2026-09-14
Status: design (autonomous ruling — phased), pre-implementation

## Problem

Karr 2012 represents the chromosome as one `CircularSparseMat` shared by every
DNA submodel: per-site bound-protein footprints, per-region linking numbers,
damaged sites, and polymerized (replicated) regions. viva-Mgen instead splits
this: a coordinate-resolved `ChromosomeDynamics` process holds a private per-bin
occupancy map (for the Fig-3 spatial layer), while the six DNA submodels
(ReplicationInitiation, DNASupercoiling, ChromosomeCondensation/Segregation,
DNADamage, DNARepair) evolve **lumped scalars** (σ, condensed/segregated fraction,
a lesion COUNT, an oriC pool). No single shared per-site structure exists, and
e.g. DNADamage produces a lesion count with no location while DNARepair clears an
abstract count — they don't operate on the same sites.

## Key ruling (autonomous, phased)

Fully unifying all six submodels onto one shared per-site structure is XL and
figure-critical: replication/supercoiling drive Fig 3/Fig 4, whose behavior tests
are tuned to the current dynamics, so a big-bang rewrite risks regressing the
reproduction's headline figures. This is therefore **phased**:

**Phase 1 (this spec):** build the shared per-site chromosome structure as a real,
tested foundation, and migrate the **DNADamage → DNARepair** subsystem onto it —
lesions become genuine per-site entries on a shared chromosome store that damage
adds to and repair clears from the SAME sites. Chosen first because it is a real
fidelity gain (located lesions vs. a lumped count) with **low figure risk**: the
baseline `damaging_agent = 0`, so no lesions arise and behavior is unchanged; the
shared structure only becomes populated under an explicit damage stimulus.

**Later phases (staged, in FIDELITY_GAPS.md):** migrate replication (per-site
polymerized regions), DNASupercoiling (per-region linking number), condensation/
segregation (per-site), and finally fold `ChromosomeDynamics`' private occupancy
into the same structure — each with its own figure re-tuning.

## What phase 1 builds

### 1. Shared chromosome-state library (`viva_mgen/chromosome_state.py`, new)

A pure operations module over a bin-indexed chromosome, pbg-native (operates on
plain `map[float]`/`list` representations, no custom bigraph type needed). The
canonical shared representation phase 1 uses is a **lesion map**: `{bin_index(str)
-> lesion_count(float)}` over `n_bins` bins.

Pure functions:
- `empty_lesion_map(n_bins) -> dict` — `{str(b): 0.0 for b in range(n_bins)}` (the pre-seed so additive deltas land; a `map[float]` store drops deltas to absent keys).
- `add_lesions(rng, n_bins, count) -> dict` — pick `count` random bins (with replacement, like independent damage events) and return the additive delta `{bin: +k}`.
- `repair_sites(lesion_map, capacity, rng) -> (delta, repaired)` — choose up to `capacity` lesion instances from bins with count>0 (weighted by count) and return the negative delta `{bin: -k}` plus the number repaired.
- `n_lesions(lesion_map) -> float` — total lesion count (Σ values), the scalar the existing `lesions` observable reports.
- `lesion_positions(lesion_map, genome_length_bp, n_bins) -> list[float]` — bp coordinates of currently-damaged bins (for a per-site observable / Fig-3-style view).

The module also documents (docstring + typed stubs, NOT implemented in phase 1)
the fuller `Chromosome` layers later phases will add on the same bin index:
bound-protein footprints, per-region linking number, polymerized regions.

### 2. DNADamage → per-site lesions (`viva_mgen/processes/dna.py`)

`DNADamageReproductionProcess`:
- Keep the Poisson lesion-count law (`new_lesions ~ Poisson((base + agent·rate)·Δt)`).
- Add config `n_bins` (default = `chromosome`'s 580) and input/output `lesion_map` (`map[float]`).
- Emit `lesion_map` = `add_lesions(rng, n_bins, new_lesions)` (per-site additive delta) AND keep emitting `lesions` (float, +new_lesions) for back-compat.

`DNARepairReproductionProcess`:
- Add input `lesion_map` (`map[float]`); read it. Repair capacity as today
  (`min(enzyme·rate·Δt, atp/atp_per_repair)`), but repair specific SITES via
  `repair_sites(...)`; emit the negative `lesion_map` delta and `lesions` = −repaired and the ATP delta as today. Repair count is bounded by the actual lesions present in the map (Σ), so it can't over-repair.

Net: damage and repair now act on the SAME per-site lesion structure; the `lesions`
scalar remains the Σ of that structure (consistent by construction).

### 3. Composite wiring (`viva_mgen/composites/mgen.py`)

- New shared store `lesion_map` in `_STORE_GROUP` (genome group), typed `map[float]`, pre-seeded via `empty_lesion_map(n_bins)` (add to the map-init logic so it seeds all bins to 0.0, not `{}` — needed so additive deltas land).
- DNADamage/DNARepair `lesion_map` ports auto-wire by name to that store.
- Optionally emit `lesion_positions` — SKIP for phase 1 (keep the emit set stable); the positions helper is available for a viz/study later.

### 4. Tests (`tests/test_chromosome_state.py`)

- Unit: `empty_lesion_map` length; `add_lesions` returns `count` total additive weight over valid bins; `repair_sites` never repairs more than present and returns matching negative delta; `n_lesions` = Σ; round-trip: apply add then repair deltas to an empty map → non-negative counts, total = added − repaired.
- Process: DNADamage with `damaging_agent>0` populates `lesion_map` (Σ ≈ its `lesions` count); with `damaging_agent=0` emits an empty/zero delta (no behavior change).
- Process: DNARepair with a seeded `lesion_map` clears sites and its `lesions` delta equals −(sites cleared); never drives a bin negative.
- Integration: composite builds + runs 3 s; with default `damaging_agent=0` the `lesions` observable stays 0 and results are unchanged (non-regression); `lesion_map` store exists and is all-zero.

## Edge cases

- Repair when `lesion_map` empty ⇒ repairs 0, emits empty delta (no negative bins).
- `add_lesions(count=0)` ⇒ empty delta.
- Bins are addressed by `str(int_bin)` to match `map[float]` string keys; helper handles int/str consistently.

## Out of scope (phase 1)

- Migrating replication/supercoiling/condensation/segregation onto the shared
  structure (staged phases — figure-critical, each needs re-tuning).
- Folding `ChromosomeDynamics`' private occupancy into the shared structure.
- Per-region linking-number and bound-protein layers (documented stubs only).

## Fidelity impact

Introduces the **shared per-site chromosome structure** (bin-indexed) that Karr's
`CircularSparseMat` embodies, and puts the damage/repair subsystem on it: lesions
are now located sites that repair clears from the same structure, not a lumped
count. Zero baseline behavior change (no damage stimulus by default). The DNADamage/
DNARepair `description` fidelity lines are updated to say lesions are now tracked
per-site on the shared chromosome structure, with the remaining submodels staged.
