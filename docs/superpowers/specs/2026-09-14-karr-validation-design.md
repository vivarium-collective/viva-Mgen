# Karr-comparison gating metrics in the study report cards

Date: 2026-09-14
Status: design, pre-implementation (revised per user: metrics live in study tests/report cards, not a standalone CLI)

## Problem

We've closed fidelity gaps but have no *quantitative, gating* measure of whether the
reproduction matches the numbers Karr 2012 reports. The figure-studies already gate
on some Karr metrics via `behavior_tests` (e.g. fig3 expression-coverage times), and
report cards evaluate pass/fail. But two problems:

1. **The most meaningful new observables aren't gated.** fig2 gates `protein_fraction`/
   `rna_fraction` off the **phenomenological** `mass_fractions` store — which is set to
   the fitted dry-weight fractions by construction, so those tests pass *tautologically*.
   The genuinely-predictive quantity is the **emergent** composition (`emergent_mass_fractions`,
   Σ species×MW, gap #2) — nothing gates it. Likewise single-cell variation (gap #7
   ensemble) is reported by Karr (Fig 2) but ungated.
2. **No single Karr-reference source.** Targets/bands are re-typed per study with prose
   provenance; there's no versioned reference the studies cite.

We cannot run the original MATLAB model here (agreed), so ground truth is the
**quantitative results the paper reports**, curated once and gated in the studies.

## Ruling

Put the comparison **in the study report cards** as new `behavior_tests` (with
Karr-cited `pass_if` bands), backed by a shared reference dataset + extraction helpers.
The new gates measure the *emergent* quantities the recent gaps produced, so they are
honest predictors (not tautologies) and will move as the remaining gap work lands —
some will FAIL at baseline (e.g. emergent composition before synthesis calibration),
which is the intended improvement signal, not a defect. Report-card behavior_tests are
evaluated in the dashboard, not the pytest CI, so failing gates don't break CI — they
show as the fidelity target.

## What v1 builds

### 1. Karr reference dataset (`datasets/karr_reference_values.json`)

Single versioned source of the paper's quantitative targets:
`{id: {description, target, band: [lo,hi], unit, source, observable}}` — e.g.
`mass_fraction_protein` (0.62, band [0.52,0.72], Karr cell composition),
`mass_fraction_dna` (0.169), `mass_fraction_rna` (0.093),
`single_cell_mrna_cv` (target ~0.3, band (0,1.5], Fig 2 single-cell distributions),
`expression_t50_min` (18, Fig 3C), `expression_t90_min` (143, Fig 3C),
`doubling_time_h` (9.0, Fig 2A). Studies read their targets from here.

### 2. Shared metric helpers (`viva_mgen/validation.py`)

Pure helpers the study `run.py` files call so metric extraction is defined once:
- `reference(id) -> dict` — load a reference entry.
- `emergent_macro_fractions(row) -> {protein,DNA,RNA}` — renormalize a row's
  `emergent_mass_fractions` over the three macromolecules (matches the Karr
  renormalized targets).
- `mrna_cv(final_totals) -> float` — coefficient of variation across an ensemble's
  per-cell final total-mRNA values.
- `score(observed, ref) -> {observed, target, band, pass, closeness}` — pass = in band;
  closeness ∈ [0,1] graded distance (so a report can show "how close" even when failing).

### 3. New gating behavior_tests + derived scalars in fig2-growth

`workspace/studies/fig2-growth/`:
- `sims/run.py` — additionally compute, from the existing full-cycle run:
  - `emergent_protein_fraction`, `emergent_dna_fraction`, `emergent_rna_fraction`
    (from `emergent_mass_fractions` of the final row, via `validation.emergent_macro_fractions`).
  - `single_cell_mrna_cv` — run a small ensemble (`viva_mgen.ensemble.run_ensemble`,
    n≈6) and compute CV of per-cell final total mRNA (`validation.mrna_cv`).
  Add these to the emitted metrics dict the report card reads.
- `study.yaml` — add `behavior_tests` (and mirrored `expected_behavior`) gating each:
  - `emergent-mass-protein-fraction` → `emergent_protein_fraction` in band from
    `karr_reference_values.mass_fraction_protein` (renormalized target ≈0.705, band
    [0.55,0.80]), provenance citing Karr composition. `classification: secondary`
    (it may fail at baseline; secondary keeps the study's primary gates meaningful
    while still evaluating + reporting the emergent metric).
  - `emergent-mass-dna-fraction`, `emergent-mass-rna-fraction` — analogous.
  - `single-cell-mrna-variation` → `single_cell_mrna_cv` in (0, 1.5], provenance
    Fig 2 single-cell distributions. `classification: secondary`.
  Each `behavior_test` has `measure: {kind: derived_scalar, field: <name>}`,
  `pass_if: {op: range, low, high, provenance: {kind: experiment, note: "Karr 2012 <fig/source>"}}`,
  `requires_simulation: baseline`.

The existing tautological `protein_fraction`/`rna_fraction` (off phenomenological
`mass_fractions`) stay as-is (they verify the growth-law composition wiring) but the
new emergent gates are the real Karr comparison.

## Tests (`tests/test_validation.py`)

- `score`: in-band pass, out-of-band fail, closeness monotone (1 at target).
- `emergent_macro_fractions`: renormalizes {protein,DNA,RNA} to sum 1 from a sample row; ignores metabolite.
- `mrna_cv`: known values → known CV; all-equal → 0.
- reference dataset: every entry has description/target/band/unit/source/observable; band lo≤hi.
- fig2 run.py derived scalars: import its metric-computation and assert the new fields are produced on a tiny hand-built rows fixture + fake ensemble (don't run the full 9h sim in the unit test).

## Edge cases

- Emergent fractions all-zero at t=0 (empty proteome) — the metric uses the FINAL row after a full cycle; if protein still ~0 the gate fails (honest signal). Guard divide-by-zero → fractions 0, gate fails with a clear note.
- Ensemble CV with n<2 or all-identical → 0 (gate fails the >0 lower bound — flags no variation).

## Out of scope (v1)

- Rewiring other studies (fig4/fig5) — same pattern, added later; v1 establishes the reference + helpers + the fig2 emergent/variation gates (the highest-value new ones).
- A standalone scorecard CLI (superseded by report-card integration; the dashboard report card IS the scorecard surface).
- Making the emergent gates `primary` — they're `secondary` until the synthesis/mass work brings them into band.

## Fidelity impact

The Karr comparison now lives where evaluation happens — the study report cards —
as gating `behavior_tests` on the **emergent** quantities (composition from gap #2,
single-cell variation from gap #7), each citing a versioned Karr reference. This
replaces tautological composition checks with genuine predictors and gives a
report-card pass/fail that the remaining gap work (#3/#4/#6, synthesis calibration)
is measured against.
