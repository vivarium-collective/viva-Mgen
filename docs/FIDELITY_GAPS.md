# viva-Mgen fidelity gaps vs. Karr 2012 — completeness backlog

Tracked list of the remaining gaps between this reproduction and the original
Karr et al. 2012 *M. genitalium* whole-cell model. Submodel **coverage** is
complete (all 28 cell-process submodels + mass); the gaps are in the whole-cell
**integration algorithm**, a few **state representations**, and a handful of
**constants the knowledge base never stored**.

Each gap is worked as its own brainstorm → spec → plan → subagent-driven
implementation cycle (specs under `docs/superpowers/specs/`, plans under
`docs/superpowers/plans/`). Worked one at a time, in the priority order below.

## Status

- [~] **Gap 1 — Resource-allocation layer.** IN PROGRESS (branch
  `feat/resource-allocation`). Finite metabolite pools replenished by metabolism,
  partitioned among the submodels each tick by a central allocator
  (demand → allocate → run). Spec: `2026-09-14-resource-allocation-design.md`.

- [ ] **Gap 2 — Emergent mass + whole-cell conservation.**
  Today `MassGrowth` grows mass phenomenologically (exp growth × `growth_fraction`,
  split by fixed dry-weight fractions). Make total mass the **emergent**
  Σ(species count × molecular weight) over the actual molecular inventory
  (RNA, protein, metabolites, DNA), and enforce whole-cell mass/atom balance so
  byproducts (water, Pi, PPi, GDP, formate) are conserved rather than "delegated
  to the pools". Deliverable: a mass process that reads the real species stores
  and MWs; the fitted fractions become an emergent *output* to validate against,
  not an input. Depends on Gap 1 (finite pools) being in. Size: L.

- [ ] **Gap 6 — Dynamic metabolism ↔ proteome coupling.**
  FBA runs over the real iPS189 network but with static bounds; the original
  updates each reaction's flux bound from **current enzyme copy numbers** every
  tick. Wire metabolic-enzyme protein counts → per-reaction `v_max` bounds each
  step (kcat × [enzyme]), so proteome composition actually gates metabolism.
  Needs the gene→reaction→enzyme map (already partly in `genes.csv`
  associated_reactions). Independent of Gap 2. Size: M.

- [ ] **Gap 3 — Unified chromosome representation.**
  Replace the split (aggregate DNA submodels + a separate coordinate-resolved
  `ChromosomeDynamics`) with ONE shared per-site chromosome structure
  (Karr's `CircularSparseMat`: bound-protein footprints, damaged sites,
  per-region linking numbers) that replication, transcription, supercoiling,
  condensation, segregation, and repair all read/write. Largest structural
  change. Depends on nothing but touches many submodels. Size: XL.

- [ ] **Gap 4 — Within-submodel state de-reductions (cluster; splittable).**
  Where the mechanism is faithful but the state is lumped:
  - 4a. Macromolecular complexation: replace greedy limiting-subunit assembly
    with the network steady-state solve (real stoichiometry is already wired).
  - 4b. RNA processing/modification: carry the species→RNase/enzyme assignment
    so the real per-RNase kcats apply per substrate (not one pooled rate).
  - 4c. tRNA aminoacylation: model the 20 AA × 37 tRNA reactions with
    per-synthetase kcats instead of one lumped pool.
  - 4d. Protein processing II: carry the per-protein lipoprotein classification
    so the Lgt transferase step applies only to lipoproteins.
  Each sub-item is its own S/M task. Some (4b, 4c) overlap Gap 3's data.

- [ ] **Gap 5 — Constants the KB never stored (research task).**
  ProteinFolding / ProteinModification rate matrices, per-monomer protein
  half-lives, and the empty RNAModification / tRNAAminoacylation / Transcriptional-
  Regulation kinetic dicts are order-of-magnitude values today. Source real values
  from the Karr supplement / primary literature where they exist; where they
  genuinely don't, document each as irreducible with its citation. Deliverable:
  a sourced constants table + wiring, or a documented "irreducible" note per
  constant. Size: M (mostly research).

- [ ] **Gap 7 — Single-cell ensembles.**
  The mechanisms are stochastic but figures use a single representative run.
  Add an ensemble runner (N independent seeded cells) and aggregate to reproduce
  Karr's cell-to-cell variation (e.g. Fig 2 distributions). Runner/infra, not
  model change. Independent of the others. Size: M.

## Recommended order & rationale

1. **Gap 1** (in progress) — the integration algorithm; everything downstream
   reads truer pools once it lands.
2. **Gap 2** — emergent mass is the other half of "whole-cell"; highest
   scientific value after allocation, and benefits from finite pools.
3. **Gap 6** — closes the metabolism↔proteome loop; medium effort, high payoff.
4. **Gap 3** — the big structural unification; do it before Gap 4b/4c which
   reuse its per-site data.
5. **Gap 4** — de-reduce the remaining lumped submodels (sub-items in order 4a,
   4d, 4b, 4c).
6. **Gap 5** — sourced constants / irreducibility documentation.
7. **Gap 7** — ensembles for variation (can run in parallel anytime; ordered
   last as it's infra, not fidelity of a mechanism).

Each process's own `description` fidelity line remains the authoritative
statement of that submodel's current gap; this file is the program-level backlog.
