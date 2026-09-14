# Dynamic metabolism ↔ proteome coupling for viva-Mgen (closing gap #6)

Date: 2026-09-14
Status: design (autonomous ruling — see "Key ruling"), pre-implementation

## Problem

`MetabolismFbaReproductionProcess` runs FBA over the real iPS189 network but with
**static** flux bounds. In the Karr 2012 model each reaction's flux capacity is
updated from the **current enzyme copy numbers** every tick (v_max = kcat·[enzyme]),
so proteome composition gates metabolism — a depleted or knocked-down enzyme
throttles its reaction. viva-Mgen's metabolism ignores the proteome (only
`nutrient_scale` + discrete `disrupted_genes` knockouts affect it).

## Key ruling (autonomous, documented)

Two constraints shape v1:

1. **Absolute kcats are unavailable.** The KB kcat tables are mostly empty (that's
   gap #5). So the coupling cannot be the literal `v_max = kcat·[enzyme]`. Instead
   use **relative enzyme gating**: scale each reaction's default bound by
   `min over its enzyme genes of clip(count_g / reference_g, floor, cap)`. At the
   reference (steady-state) proteome the factor ≈ 1 (default behavior); a
   knockdown drives it < 1 (throttle). This is faithful in mechanism (proteome
   gates flux) without fabricating kcats.

2. **The composite starts with a near-empty proteome** (protein_counts build from
   0). Gating metabolism ON by default would throttle the early cell cycle and
   risk regressing the already-shipped fig2/3/5 studies (whose behavior tests are
   tuned to the current build-from-zero dynamics). Making enzyme coupling the
   default therefore requires **birth-proteome seeding** (cells born with a
   proteome) — a separate faithfulness change with its own figure-retuning.

**v1 delivers the coupling MECHANISM, opt-in (`enzyme_coupling`, default False).**
Default-off means every existing figure/study is byte-for-byte unchanged (zero
regression). Turned on, metabolism's reaction bounds are gated by the live
proteome relative to a steady-state reference, so a metabolic-enzyme knockdown
graded-throttles its reaction and lowers growth. A dedicated test proves both
directions. `FIDELITY_GAPS.md` records that flipping it on by default is gated on
birth-proteome seeding + figure re-tuning.

## What v1 builds

### 1. Reference steady-state proteome (`viva_mgen/expression_defaults.py`)

Add `reference_protein_counts() -> dict` (gene-key → expected steady-state protein
count), from the existing panel:
`mRNA_ss = synthesis_rate / mrna_decay_rate`;
`protein_ss = translation_rate · mRNA_ss / protein_decay_rate`.
Keyed the same way as `protein_counts` (gene symbol/id). This is the model's own
self-consistent steady state, so scaling relative to it means "a gene above/below
its expected level speeds/slows its reactions."

### 2. Enzyme-gated bounds in metabolism (`viva_mgen/processes/metabolism.py`)

- New config: `enzyme_coupling` (bool, default **False**), `reference_protein_counts`
  (map[float], default `{}`), `enzyme_coupling_floor` (float, default 0.0),
  `enzyme_coupling_cap` (float, default 1.0).
- New input: `protein_counts` (map[float], read-only sensor).
- In `update()`, inside the existing reverting `with model:` context and BEFORE
  `optimize()`, when `enzyme_coupling` is on and a reference is present:
  - Build a normalized reference: `ref[norm_id] = count` via `normalize_gene_id`.
  - Build normalized live counts from `protein_counts` the same way.
  - For each reaction with `reaction.genes`, factor =
    `clip(min over g in reaction.genes of (live[norm(g.id)] / ref[norm(g.id)]),
    floor, cap)` — genes missing from the reference don't gate (contribute 1).
  - Scale the reaction's positive `upper_bound` and negative `lower_bound` by factor.
  - Reactions with no genes, or where no gene has a reference, are ungated (factor 1).
- Default off ⇒ this whole block is skipped ⇒ identical FBA to today. The existing
  `disrupted_genes` (hard knockout) and `reaction_bound_scale` paths are unchanged
  and compose with gating (knockout wins).

### 3. Composite wiring (`viva_mgen/composites/mgen.py`)

- Metabolism gains a read-only `protein_counts` input wired to the existing
  `proteome/protein_counts` store (metabolism emits NO delta to it).
- `build_mgen` gains an `enzyme_coupling: bool = False` parameter, passed to the
  metabolism node config, plus the default `reference_protein_counts` (from the
  helper). Default False ⇒ no behavior change; a study/test can build with True.
- The `mycoplasma_genitalium` generator exposes `enzyme_coupling` as a parameter
  (default False) so it's toggleable in the dashboard without code change.

### 4. Tests (`tests/test_dynamic_metabolism.py`)

- Unit: `reference_protein_counts()` returns positive counts for expressed genes,
  keyed compatibly with `protein_counts`.
- Coupling OFF (default): metabolism growth == the current wild-type growth
  (non-regression) regardless of `protein_counts`.
- Coupling ON, reference proteome: growth ≈ wild-type (factor ≈ 1).
- Coupling ON, a metabolic enzyme's protein driven to ~0: its reaction's flux
  drops and growth_fraction is strictly lower than the reference case — the
  graded generalization of a knockout.
- Composite builds + runs with `enzyme_coupling=True` (3 s), pools/growth sane.

## Edge cases

- Reference count 0 or missing for a gene ⇒ that gene doesn't gate (factor 1) — never divide by zero.
- `protein_counts` empty (cold start) with coupling ON ⇒ every gated reaction hits `floor`; with floor 0 that stalls growth — which is why coupling is default-OFF and callers enabling it should also seed a birth proteome (documented).
- A reaction gene not in the metabolic model's gene set is simply absent from `reaction.genes` — no effect.

## Out of scope (v1)

- Absolute `kcat·[enzyme]` bounds (needs kcats — gap #5).
- Making coupling the default (needs birth-proteome seeding + figure re-tuning).
- Isozyme (OR-rule) handling beyond the conservative `min` over all reaction genes;
  documented as a reduced approximation.

## Fidelity impact

Metabolism can now be **dynamically gated by the live proteome** (relative
enzyme availability), the graded generalization of the discrete knockout it
already supports — the core of Karr's per-tick enzyme→flux-bound coupling, minus
the absolute kcats the KB lacks. Delivered opt-in so no shipped figure regresses;
the metabolism `description` fidelity line is updated to state the coupling is
available (`enzyme_coupling`) and what enabling it by default awaits.
