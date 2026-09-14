# Emergent mass + whole-cell composition for viva-Mgen (closing gap #2)

Date: 2026-09-14
Status: design (autonomous ruling — see "Key ruling"), pre-implementation

## Problem

`MassGrowthReproductionProcess` grows dry mass phenomenologically:
`Δmass = mass·(exp(μ·growth_fraction·Δt)−1)`, then splits it by the fixed fitted
dry-weight fractions (`DRY_MASS_FRACTIONS`). In the Karr 2012 model, total mass is
the **emergent** Σ(species count × molecular weight) over the actual molecular
inventory; the dry-weight fractions are an *output* of that inventory, not an
input. So today mass is not derived from what the cell actually contains, and the
composition is asserted rather than emergent.

## Key ruling (autonomous, documented)

The full form — make the primary `mass` the emergent Σ(species×MW) and fire
division at 2× birth of THAT — is **not safe for v1** and is deliberately scoped
out. Reason: the reduced, ParCa-fitted synthesis/decay balance drives the
proteome toward a **steady state**, not the sustained net doubling a real cell
cycle needs. Firing division off emergent total mass would therefore likely stall
the cell cycle and regress the already-shipped fig2/3/5 studies. Recalibrating the
whole synthesis network for net doubling is a large, figure-risking effort that
belongs to a later gap.

**v1 delivers exactly the FIDELITY_GAPS.md deliverable:** a mass process that
computes the emergent molecular mass Σ(species×MW) from the real species stores and
real MWs, exposes it and its emergent composition as observables, and **validates**
the emergent composition against the Karr fitted fractions. The fitted fractions
become a validation target, not the source of the reported composition. The
calibrated growth law continues to drive division (unchanged), so no figure
regresses. A follow-up (noted in FIDELITY_GAPS.md) can switch division to emergent
mass once synthesis is calibrated for doubling.

## What v1 builds

### 1. Molecular-weight helper (`viva_mgen/mass_composition.py`, new)

Pure functions computing per-species dry mass in grams from counts + real MWs:

- `rna_mass_g(rna_counts: dict) -> float` — Σ_g count·length_nt[g]·`_RNA_NT_MW`(340) / N_A. Uses `expression_defaults.gene_lengths()` (nt per gene) already keyed the same way as `rna_counts`.
- `protein_mass_g(protein_counts: dict) -> float` — Σ_g count·(length_nt[g]/3)·`_AA_MW`(≈110) / N_A.
- `dna_mass_g(chromosome_copy: float) -> float` — `chromosome_copy`·`GENOME_LENGTH_BP`·`_BP_MW`(≈660, both strands) / N_A.
- `metabolite_mass_g(pools: dict) -> float` — Σ pool_count·MW / N_A for the tracked metabolite pools (atp/gtp/ntp/amino_acid) with representative MWs (ATP 507, GTP 523, NTP≈500 avg, amino acid≈110).
- `emergent_composition(rna_counts, protein_counts, chromosome_copy, pools) -> dict` — returns `{RNA, protein, DNA, metabolite: grams}` and an inclusive `total` (RNA+protein+DNA+metabolite). Callers validating against Karr's macromolecule fractions sum RNA+protein+DNA themselves rather than using `total` — see ruling below.

Constants (`_AA_MW`, `_BP_MW`, `_N_A=6.022e23`) live in this module or `constants.py`; `_RNA_NT_MW` is reused from where it already is (parca) or duplicated as a named constant here — one source, cited.

### 2. Emergent outputs on the mass process

`MassGrowthReproductionProcess` gains:
- inputs (read-only sensors): `rna_counts` (map), `protein_counts` (map), `chromosome_copy` (float), and the metabolite pools `atp`/`gtp`/`ntp`/`amino_acid` (floats).
- outputs (overwrite observables): `emergent_mass` (float, fg; RNA+protein+DNA only) and `emergent_mass_fractions` (map[float]: RNA/protein/DNA fraction of `emergent_mass`, summing to 1). The free-metabolite pool is reported separately as `metabolite_mass` (float, fg) — a non-validated diagnostic, NOT included in `emergent_mass`/`emergent_mass_fractions`, because current pool sizes are uncalibrated allocation placeholders rather than a fitted target.
- The existing `mass`/`volume`/`division`/`mass_fractions` behavior is UNCHANGED (growth-law still drives division).
- `update()` computes the emergent composition from the sensor inputs and emits the three new observables each tick: `emergent_mass`/`emergent_mass_fractions` from RNA+protein+DNA, and `metabolite_mass` separately.

### 3. Composite wiring + emit

- Add `emergent_mass`, `emergent_mass_fractions`, and `metabolite_mass` to `_EMIT` and to the mass node's ports (they resolve to `physiology`/`metabolism` group stores; add `_STORE_GROUP` entries: `emergent_mass`→physiology, `emergent_mass_fractions`→physiology, `metabolite_mass`→physiology, as `_EMPTY_MAP_STORES` for the map).
- The mass node already exists; it just gains the sensor inputs + observable outputs. The sensor inputs wire to the existing `rna_counts`/`protein_counts`/`chromosome_copy`/pool stores (read-only; mass emits no delta to them).

### 4. Validation test

`tests/test_emergent_mass.py`:
- Unit: `rna_mass_g`/`protein_mass_g`/`dna_mass_g` give known values for hand-set counts (e.g. one 300-nt gene at count 1 → 300·340/N_A g).
- Composition: `emergent_composition` on a representative inventory yields fractions (over RNA/protein/DNA only, excluding metabolite) whose ORDERING and rough magnitudes match Karr (protein dominant, then RNA/DNA) — assert protein fraction is the largest.
- DNA sanity: `dna_mass_g(1.0)` ≈ 0.6–0.7 fg (matches the Karr DNA dry-mass fraction 0.1688 × 3.93 fg ≈ 0.66 fg), confirming `_BP_MW` is right.
- Integration: build the composite, run a few ticks, assert `emergent_mass` > 0 (RNA+protein+DNA only) and `emergent_mass_fractions` — keyed only `{RNA, protein, DNA}`, no `metabolite` — sums to ≈1.0 with protein the largest component; `metabolite_mass` is emitted separately and is not part of the validated total/fractions.
- Zero-total: `emergent_composition` with an all-empty/zero inventory returns all components 0.0 and `total` 0.0, documenting the input the process's divide-by-zero fraction guard handles.

## Data flow

```
rna_counts, protein_counts, chromosome_copy, atp/gtp/ntp/amino_acid
        │ (read-only sensors)
        ▼
MassGrowthReproductionProcess.update():
   emergent_mass = Σ(species × MW)/N_A  (fg)
   emergent_mass_fractions = per-component / total
   (mass/volume/division: growth-law, unchanged)
```

## Edge cases

- Empty inventory (tick 0): `emergent_mass` (RNA+protein+DNA) ≈ DNA-only (chromosome_copy=1); `emergent_mass_fractions` still sum to 1 (guard divide-by-zero → all-zero fractions if `emergent_mass` is 0).
- A gene in `rna_counts`/`protein_counts` with no length in `gene_lengths()`: fall back to an average length (1000 nt) rather than dropping it.

## Out of scope (v1)

- Driving division off emergent mass (needs synthesis calibrated for net doubling — figure-risking; deferred).
- Strict per-atom conservation of byproducts (water/Pi/PPi) — the metabolite pools already track the energy carriers; full atom-level balance is a larger effort tracked separately.

## Fidelity impact

Mass composition is now **computed from the real molecular inventory with real
MWs** and validated against Karr's fitted fractions, rather than asserted. The
validated `emergent_mass`/`emergent_mass_fractions` cover RNA+protein+DNA only —
the macromolecules Karr's fitted fractions actually target. The free-metabolite
pool is reported separately as `metabolite_mass`, a non-validated diagnostic
excluded from the total, since current pool sizes are uncalibrated allocation
placeholders rather than a fitted target. The `MassGrowth` `description`
fidelity line is updated to say the emergent macromolecule composition is
reported + validated, and that the growth-law total still drives division
pending synthesis calibration.
