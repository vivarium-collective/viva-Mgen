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

- [x] **Gap 1 — Resource-allocation layer.** DONE (merged #13) (branch
  `feat/resource-allocation`). Finite metabolite pools replenished by metabolism,
  partitioned among the submodels each tick by a central allocator
  (demand → allocate → run). Spec: `2026-09-14-resource-allocation-design.md`.

- [x] **Gap 2 — Emergent mass + whole-cell conservation.** DONE.
  Macromolecule composition (merged #14): total mass is the **emergent**
  Σ(species count × molecular weight) over the actual molecular inventory
  (RNA, protein, metabolites, DNA); the fitted dry-weight fractions became an
  emergent *output* validated against Karr's values, not an input.
  **Synthesis/mass calibration (this round):** the just-added Karr-comparison
  report card exposed that the emergent inventory was grossly protein-deficient
  (emergent protein fraction ≈ 0.002 vs Karr's ~0.70). Root cause was two
  coupled defects, both now fixed:
  1. *Energy-carrier units.* Metabolism supplied ATP/GTP to the finite pools as
     raw FBA flux (~0.2–0.8 molecules/tick) while NTP/AA were molecule-scaled, so
     GTP drained in the first coarse tick and translation stalled after ~1600
     proteins. Metabolism now emits per-second precursor SUPPLY rates
     (`<pool>_supply` = base·growth_fraction), the allocator replenishes
     `rate·interval` (timestep-consistent at dt=1 and dt=300), and the raw FBA
     `atp/gtp_production` flux is retained as the Fig 5 signal.
  2. *Greedy per-gene GTP.* Translation consumed its GTP grant first-come in dict
     order, so one high-demand gene (MG471, anomalously stable mRNA × max rate)
     ate 96% of the budget and starved the rest (only 45 genes translated). It
     now shares the budget **proportionally** across genes (no gene monopolizes;
     ~273 genes translate).
  3. *Length-proportional GTP cost.* Translation charged a flat `gtp_per_protein`
     (600) per chain; it now charges `gtp_per_peptide_bond` (2) · (length_g/3),
     the faithful ~2 GTP/peptide-bond mechanism — long proteins draw more of the
     shared budget than short ones. fig2 fractions are unchanged (re-fitting
     `gtp_base_supply` to 3500 holds them at ~0.69:0.13:0.18) and fig5's energy
     budget becomes length-accurate (transcription share 0.001→0.006, translation
     0.997→0.985, both toward Karr).
  4. *Stable-RNA level.* The emergent RNA mass was over-represented (fraction
     ~0.18 vs Karr ~0.11), which dragged DNA low (~0.13 vs ~0.19). Splitting the
     inventory by RNA class showed the excess is entirely STABLE RNA — rRNA (~73%
     of RNA mass) and tRNA (~19%); real mRNA is only ~9% and at its proper level.
     Those stable species have no packaging/charging sink here and half-lives
     (rRNA ~20 h, tRNA ~45 min) exceeding the 9 h cycle, so they accumulate
     ~linearly instead of reaching steady state. `expression_defaults`
     .STABLE_RNA_SYNTHESIS_SCALE (0.37) rescales rRNA/tRNA/SRP-RNA synthesis
     (a proxy for the missing sink) — mRNA synthesis is untouched, so gene-
     expression timing (fig3 t50/t90) and the translation-feeding mRNA pool are
     unchanged. Stable-RNA genes are identified by `parca.rna_type`.
  5. *RNA-type classifier + protein-coding translation.* `parca.rna_type` had
     substring-matched "trna"/"rrna", so tRNA-synthetase / rRNA-methyltransferase
     PROTEINS were misclassified as non-coding RNA — which (a) gave those protein
     mRNAs stable-RNA half-lives in the fitted panel and (b) let the tRNA/rRNA
     GENES themselves be "translated" into spurious proteins: **92% of the
     emergent protein COUNT (~286k of ~311k copies) was non-coding genes**, led by
     76-nt tRNAs turned into 25-aa "proteins". `rna_type` now matches the RNA
     PRODUCTS (tRNA-<AA> / "ribosomal rRNA" / SRP 4.5S RNA), never the proteins
     that process them, and `translation_rates()` translates only mRNA genes. The
     GTP freed from the fake proteins flows to real ones, so the dry-mass fractions
     are unchanged (still ~0.70:0.18:0.11) while the protein population becomes
     realistic: **count ~311k → ~135k, mean length ~80 → ~198 aa**, zero non-coding
     "proteins".
  With gtp_base_supply (2500) and STABLE_RNA_SYNTHESIS_SCALE (0.37) fitted
  jointly, the emergent protein:DNA:RNA fractions land on ~0.703:0.185:0.113
  across seeds — essentially Karr's 0.703:0.192:0.105 (all three fig2 bands
  green; was 0.002:0.42:0.57). fig3 (t50/t90, exploration, collisions) and fig5
  (energy shares) stay green; the parca-fitting closure stays feasible.
  DEFERRED / investigated follow-ups: (a) protein length — NOT a defect on
  inspection: the by-gene mean of translated proteins is ~369 aa (matches real
  M. genitalium ~330–360; all 482 protein-coding genes translate); the ~198 aa
  copy-weighted mean is just abundant proteins skewing it, as expected. (b)
  whole-cell mass/atom balance (conserve water, Pi, PPi, GDP, formate) and
  division-driving from emergent mass — genuinely open, own task. (c) fig5's
  `transcription_energy_share` gate [0.04,0.12] (Karr's 7.1% of ATP+GTP): the
  transcription energy is now counted at the FULL stable-RNA synthesis rate (the
  pool-shaping STABLE_RNA_SYNTHESIS_SCALE is a degradation-sink proxy, not a real
  transcription-rate cut), raising it ~2%→2.5%. The residual gap is the reduced
  model translating each transcript ~3× more than Karr (higher translation-per-
  mRNA); the study itself documents it won't hit Karr's exact split, and closing
  it would perturb the fig2 protein calibration — a genuine reduced-model limit,
  not a defect.

- [x] **Gap 6 — Dynamic metabolism ↔ proteome coupling.** DONE (default-on).
  FBA runs over the real iPS189 network; each gene-associated reaction's flux
  bound is now scaled every step by the live enzyme abundance over a reference
  (`clip(live/reference, floor, cap)`), so proteome composition gates metabolism.
  The two blockers are resolved:
  - *Birth proteome.* The cell is seeded with `expression_defaults.birth_proteome()`
    (~half the emergent division proteome, emergent-distributed), so enzyme_coupling
    has live enzymes at t=0 (was: empty proteome → every reaction clipped to the
    floor → growth collapsed on tick 1) and the cell cycle is birth→double rather
    than accumulate-from-zero. `gtp_base_supply` re-fitted 2500→1900 for the halved
    per-cycle protein synthesis; fig2 fractions stay on Karr (~0.68:0.20:0.12).
  - *Consistent reference.* enzyme_coupling's reference is the emergent DIVISION
    proteome (`coupling_reference()` = 2×birth ≈ 135k), NOT the analytic
    `reference_protein_counts` (~3.0M, 22× too high, which made live/ref ~0.05 and
    starved every reaction). live/ref runs ~0.5 at birth → ~1 near division.
  - *Graceful floor.* `enzyme_coupling_floor` default 0.0→0.1 so a stochastic dip
    in one essential enzyme can't zero its reaction and collapse growth. Hard gene
    KNOCKOUTS still use the separate `disrupted_genes` path, so Fig 6 essentiality
    is unchanged (verified identical to pre-coupling).
  HONEST NOTE: at the healthy calibrated proteome the coupling is ACTIVE but
  NON-BINDING — iPS189 has flux slack (growth is flat down to ~4% of reference,
  per the Fig 7 kcat sweep), and the normal cycle only ranges live/ref 0.5→1.0, so
  growth is unaffected at health; the coupling bites only under severe enzyme
  depletion. So the mechanism is present and correct (the defining whole-cell
  feature) but latent in the healthy cell. The absolute kcat·[enzyme] form still
  awaits the KB kcats (gap #5). All 7 studies verified; 133 tests pass.

- [~] **Gap 3 — Unified chromosome representation.** IN PROGRESS.
  Replace the split (aggregate DNA submodels + a separate coordinate-resolved
  `ChromosomeDynamics`) with ONE shared per-site chromosome structure
  (Karr's `CircularSparseMat`: bound-protein footprints, damaged sites,
  per-region linking numbers) that replication, transcription, supercoiling,
  condensation, segregation, and repair all read/write. Largest structural
  change. Depends on nothing but touches many submodels. Size: XL.
  PHASE 1 landed: shared per-site chromosome structure (`chromosome_state.py`)
  with DNADamage/DNARepair migrated to per-site lesions on a shared
  `lesion_map`. STAGED for later phases: replication polymerized-regions,
  DNASupercoiling per-region linking number, condensation/segregation
  per-site, and folding ChromosomeDynamics' occupancy in.

- [ ] **Gap 4 — Within-submodel state de-reductions (cluster; splittable).**
  Where the mechanism is faithful but the state is lumped:
  - 4a. Macromolecular complexation: replace greedy limiting-subunit assembly
    with the network steady-state solve (real stoichiometry is already wired).
  - 4b. RNA processing/modification: carry the species→RNase/enzyme assignment
    so the real per-RNase kcats apply per substrate (not one pooled rate).
  - 4c. tRNA aminoacylation: model the 20 AA × 37 tRNA reactions with
    per-synthetase kcats instead of one lumped pool.
  - 4d. Protein processing II: carry the per-protein lipoprotein classification
    so the Lgt transferase step applies only to lipoproteins. **DONE (phase 1)** —
    per-protein maturation classification decoded from KB (14 lipoproteins, 20 secretory,
    35 N-terminal-Met-cleavage) now routes ProcessingII (Lgt only), ProcessingI
    (Met-cleavage subset), and Translocation (secretory+lipoprotein only).
  
  Each sub-item is its own S/M task. Some (4b, 4c) overlap Gap 3's data.

  **Known follow-up:** the read-only dashboard's `scripts/regen_composite_state.py`
  docstring→`_contract` splitter garbles `description`/`math` fields for processes
  with multi-line indented equation blocks (pre-existing; affects protein_folding,
  complexation, maturation). Cosmetic loom-rendering only; authoritative `description`
  attributes are correct. Track as regen/parser fix.

- [x] **Gap 5 — Constants the KB never stored (research task).** DONE
  (provenance audit + supplement sourcing; irreducible constants documented).
  Every numeric kinetic constant across all process `config_schema`s (83
  total) is now classified by source tier via
  `viva_mgen.provenance.audit_constants()`, dumped to
  `datasets/constant_provenance.json`, and written up in
  `docs/CONSTANT_PROVENANCE.md`. Final breakdown: **real_kb 11**, **real_supplement
  1**, **order_of_magnitude 61**, **irreducible 10**. The only Karr 2012 SI file in
  this workspace (16 pages) is the *Cell* main text + two SI figures, not the
  Data S1 parameter appendix, so it yielded exactly one explicit unambiguous
  match (`MassGrowthReproductionProcess.cell_cycle_length_s` = 32,400 s / 9.0 h,
  Figure 2A). The 10 `irreducible` constants (ProteinFolding
  spontaneous_rate/chaperone_rate, ProteinModification
  modification_specific_rate, RNAModification modification_rate/enzyme_kcat,
  tRNAAminoacylation synthetase_kcat, DNADamage base_rate/agent_rate,
  ProteinActivation default_k, TerminalOrganelleAssembly threshold) each have a
  KB dict confirmed empty and an SI-silent note naming the missing Data S1
  table that would be needed to source a real value. The remaining
  order-of-magnitude constants are stoichiometric facts, initial conditions,
  or engineering choices, not missing rate measurements.

- [x] **Gap 7 — Single-cell ensembles.** DONE (ensemble runner + aggregation + demo CLI landed).
  The mechanisms are stochastic but figures use a single representative run.
  Add an ensemble runner (N independent seeded cells) and aggregate to reproduce
  Karr's cell-to-cell variation (e.g. Fig 2 distributions). Runner/infra, not
  model change. Independent of the others. Size: M. Delivered: `run_ensemble`, `aggregate`,
  `final_values` (viva_mgen.ensemble); `scripts/run_ensemble.py` CLI with JSON aggregates
  (mean±std + per-cell + end-of-cycle distributions). A dashboard variation-figure study
  consuming it is staged.

## Comparison harness

Karr 2012 reference values (`datasets/karr_reference_values.json`) and
report-card gating metrics have landed: `viva_mgen.validation` scores emergent
macromolecular composition and single-cell mRNA variation against the paper's
reported values, and `fig2-growth` gates on them as four new secondary
`behavior_tests` (emergent protein/DNA/RNA fraction, single-cell mRNA CV) —
mirrored in `expected_behavior`. This pattern is extensible to other figures
as their studies land. The emergent-composition gates (Gap 2's Σ(species×MW)
output vs. Karr's fitted fractions) are the current fidelity target for the
mass/synthesis work.

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
