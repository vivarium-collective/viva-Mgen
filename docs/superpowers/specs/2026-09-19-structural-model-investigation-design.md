# Design: `structural-model` investigation — a 3D whole-cell structural model of *M. genitalium* with viva-parsimony

- **Date:** 2026-09-19
- **Workspace:** viva-Mgen (`viva_mgen`)
- **Branch:** `investigation/structural-model` (worktree `~/code/viva-mGen--structural-model`, off `origin/main` @ f2f9d09)
- **Status:** approved design → implementation plan next

## 1. Goal

Build **as faithful a reproduction as possible of the Maritan et al. 2022 3D
whole-cell structural model of *Mycoplasma genitalium*** — but packed with
**viva-parsimony** (the `parsimony` cellPACK-style engine, via `pbg_parsimony`)
in place of Maritan's CellPACK — and make the packed cell viewable from this
workspace so the investigation's results can be seen.

Reference: Maritan, Autin, Karr, Covert, Olson & Goodsell, *"Building Structural
Models of a Whole Mycoplasma Cell"*, J. Mol. Biol. 434 (2022) 167351
(`workspace/references/papers/`, BibTeX `maritan2022structural`; supplement
tables S1/S2/S3 in `papers/Maritan2022_supplementary/`).

This is the *M. genitalium* analog of the existing **`3d-ecoli`** repo
(`v2ecoli` → `pbg_parsimony`); here it is **`viva_mgen` → `pbg_parsimony`**,
following the same Maritan method that `3d-ecoli` itself cites. `3d-ecoli`'s
`ecoli_3d/build.py` is the reference bridge pattern.

## 2. Research question & hypothesis

**Question.** Can the Maritan 2022 3D structural model of *M. genitalium* be
reproduced with viva-parsimony instead of CellPACK, and does the `viva_mgen`
WCM reproduction's *simulated* cell state yield a comparable 3D cell?

**Hypothesis.** The Maritan ingredient roster (S1) + genome geometry (S2),
resolved to structures and packed at true abundance with the parsimony octree
engine into a single-membrane capsule cell, reproduces Maritan's headline
mesoscale properties (species placed, crowding / volume occupancy, membrane
localization, a single supercoiled nucleoid). Swapping the abundances for a
`viva_mgen` simulated run yields a structurally comparable cell, differing only
where the reproduction's predicted proteome diverges from the published counts.

## 3. Integration (decided: import + bridge)

- Add `pbg_parsimony` to `workspace.yaml` `imports:` (mode `reference`) and as a
  git dependency in `pyproject.toml`. No vendoring; tracks upstream — mirrors
  `3d-ecoli`'s `workspace.yaml`.
- The `parsimony` binary is already built at `~/code/parsimony/target/release/parsimony`;
  resolved via `PARSIMONY_HOME=~/code/parsimony` (or `PARSIMONY_BIN`). Document
  this as a build-time requirement (not a Python dep).
- **Results are viewed** through `pbg_parsimony`'s bundled viewer
  (`pbg_parsimony/viewer/`), surfaced in the viva-Mgen dashboard/report as a
  Visualization — the "so we can see the results" piece. Mirrors how `3d-ecoli`
  publishes its pack + viewer.

## 4. Bridge package: `viva_mgen/structural/`

The organism-specific half of the pipeline (the generic engine stays in
`pbg_parsimony`). Modules:

### `maritan_tables.py`
Parse the supplement once into committed data artifacts under
`viva_mgen/structural/data/` (so runtime builds don't depend on `.xlsx`):

- **S1** (`S1_MG_proteins_ingredients.xlsx`, ~996 rows) columns:
  `ProtID, Name, Function, Compartment (c|m|…), Type (monomer/complex),
  Molecular Weight, Sequence Length, Confidence (HHpred), Complex Biosynthesis,
  DNA footprint, DNA binding, Structural Model (PDB id), Chain Selection,
  Template, Data Source, Quality (ModFOLD/VoroMQA), Comments`.
  → per-species ingredient record: id=`ProtID`, display_name=`Name`,
  category=`Function`, compartment (`c`→cytoplasm/interior, `m`→membrane/surface),
  and a **structure ref**: PDB/mmCIF id from "Structural Model" when present,
  else AlphaFold via the S2 UniProt accession.
- **S2** (`S2_MG_genes.xlsx`) columns:
  `GeneID, Type, Coordinates, Length, Direction, Transcription Unit, Essential,
  Name, UniProt, Protein Monomer`.
  → gene records for (a) the UniProt→AlphaFold fallback and (b) genome geometry
  (coordinates + direction) for the `Chromosome`.
- **S3** (membrane-protein compartment assignments) — optional refinement of the
  `m` compartment / periplasmic-vs-membrane placement; not required for v1.

### `counts.py`
Two abundance providers returning `{ProtID: count}`:

- `maritan_counts()` — the **baseline**: WC-MG static copy numbers.
  **⚠ Open dependency (§7):** S1 carries *no* copy-number column, so this needs
  a WC-MG counts source. Resolve at implementation time.
- `mgen_sim_counts()` — the **variant**: per-protein counts read from a
  `viva_mgen` simulated run's protein-monomer store (the reproduction's genome-
  scale expression panel).

### `build.py`
`build_mgen_pack(counts_source, *, out_dir, top_n=None, state="birth")`:

1. Load parsed S1/S2 records; select species (all, or the `top_n` most abundant).
2. Map each to a `pbg_parsimony.Ingredient` (structure ref, count from
   `counts_source`, region from compartment, category/display_name from S1).
3. Cell geometry: a single **`Capsule`** sized to the *M. genitalium* cell
   volume (~0.3 µm scale). **No** gram-negative two-membrane `envelope` — MG has
   one membrane, no outer membrane / periplasm. Shape simplified (capsule), and
   the terminal organelle omitted — exactly the simplifications Maritan made.
4. A single supercoiled **`Chromosome`** from S2: 580,076 bp, one circular
   chromosome (birth state), `genome_csv` from S2 coordinates; optional
   RNAP/ribosome placement deferred to a later iteration.
5. Call `pbg_parsimony.build_pack(...)` → `pack.json` + sidecar + recipe.

### `pack_step.py`
A `mgen-structural` process-bigraph `Step` wrapping `build_mgen_pack`, so the
pack is a first-class workspace composite (mirrors `ecoli_3d/composite.py` /
`ecoli_3d/pack_step.py`).

## 5. Investigation & studies

Create via `/viva-investigation new structural-model` (draft PR, `investigation:`
title prefix per `AGENTS.md`). Two studies:

- **`s01-maritan-baseline`** — Maritan roster + structures + S2 genome +
  WC-MG static counts. The faithful reproduction. Acceptance bands: species
  placed (vs S1 roster size), cytoplasm crowding / volume occupancy in Maritan's
  reported range, membrane proteins seated in the bilayer, one supercoiled
  chromosome present.
- **`s02-reproduction-driven`** — same roster/structures, abundances from a
  `viva_mgen` simulated run. Contrast study: "published WC-MG" vs "what our
  reproduction predicts" (occupancy + per-category species counts).

## 6. Viewer / report wiring

Publish each pack to the bundled `pbg_parsimony` viewer and register a
workspace **Visualization** pointing at the pack (`/viva-viz` or a direct
`workspace.yaml` entry), so both cells are browsable from the dashboard already
GitHub-Pages-published. A reproducible build/publish script under
`viva_mgen/structural/publish/` (thin analog of `ecoli_3d/publish/`).

## 7. Open dependency — baseline copy numbers

S1 has **no abundance column**. Maritan drew counts from WholeCellKB. At
implementation time, in priority order:

1. Reuse WC-MG per-protein static counts already embedded in `viva_mgen`
   (`constants.py` / `datasets/`), if present — faithful to Maritan's method.
2. Otherwise a WholeCellKB proteins export committed under `structural/data/`.
3. Otherwise a single `viva_mgen` steady-state run, documented as
   "WC-MG-derived" (blurs baseline toward variant — least preferred).

Resolved when `counts.py` is implemented; does not block the plan.

## 8. Tests

- S1/S2 parse → expected record counts / a few known rows (DNA_GYRASE PDB 6RKW).
- ProtID→structure mapping: PDB when S1 provides one, else AlphaFold via UniProt.
- Both counts providers return `{ProtID: int}` over the roster.
- Small end-to-end pack: a handful of species + chromosome through the real
  `parsimony` binary → a valid `pack.json` with the expected ingredient count.

## 9. Fidelity caveats (honest, per workspace convention)

- **Packer:** parsimony (cellPACK-*style*), not Maritan's CellPACK.
- **Shape:** capsule-simplified; terminal organelle omitted (as Maritan did).
- **Structures:** PDB where S1 gives one, else AlphaFold by UniProt; homology
  quality per S1's ModFOLD/VoroMQA columns.
- **Copy numbers:** WC-MG static counts (baseline) — §7.

## 10. Scope / YAGNI

- Birth-state cell first; a dividing (two-nucleoid, septum) cell is a later
  stretch, not v1.
- No new VR work — reuse the viewer's existing WebXR support.
- `top_n` species-richness knob like `3d-ecoli`.

## 11. PR topology (finish-line note)

Per `AGENTS.md`: the **bridge is reusable infrastructure** and should ultimately
ship via a **feature PR against `main`** (`feat(structural): …`), while the
**investigation branch carries content** (study YAMLs, packed cells, the Maritan
references, rendered report). To keep momentum this is prototyped whole on the
investigation branch first (as `3d-ecoli` was prototyped inside `v2ecoli/structural/`),
then `viva_mgen/structural/` is carved into a feature PR once stable. The
investigation branch is **draft, not a merge target**.
