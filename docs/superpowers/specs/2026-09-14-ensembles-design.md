# Single-cell ensembles for viva-Mgen (closing gap #7)

Date: 2026-09-14
Status: design, pre-implementation

## Problem

The submodels are stochastic (Poisson synthesis/decay, multinomial partitioning,
random DNA binding), but every figure is produced from a SINGLE representative run
with a fixed seed. Karr 2012 characterizes **cell-to-cell variation** by running
ensembles of independent cells (e.g. the Fig 2 single-cell distributions). viva-Mgen
has no ensemble runner: `build_mgen(seed=...)` already seeds the stochastic
processes, so different seeds give independent trajectories — but nothing runs and
aggregates a population.

## What v1 builds

A small, reusable ensemble layer (infrastructure, no model change).

### 1. `viva_mgen/ensemble.py` (new)

- `run_ensemble(n_cells, duration, *, build_kwargs=None, seeds=None, core=None) -> dict`
  Runs `n_cells` independent composites (via `build_mgen`), one per seed
  (`seeds` or `range(n_cells)`), each for `duration` seconds, and gathers each
  cell's emitter rows. Returns `{"seeds": [...], "cells": [rows_per_cell]}` where
  each `rows_per_cell` is the `gather_emitter_results(sim)[("emitter",)]` list.
  Builds one core (via `build_core`) if none given; passes `seed=<seed>` into
  `build_kwargs` per cell so trajectories are independent.
- `aggregate(cells, observable) -> dict`
  For a SCALAR observable present in the emitter rows, aligns cells to the common
  (minimum) row count and returns `{"t_index": [...], "mean": [...], "std": [...],
  "per_cell": [[...], ...], "n": n_cells}` — mean and population std across cells
  at each timepoint. Ignores cells lacking the observable; raises a clear error if
  none have it.
- `final_values(cells, observable) -> list[float]` — the last value of `observable`
  per cell (the distribution at cycle end), for histogram-style single-cell
  variation plots.

Pure/deterministic given seeds; no global state.

### 2. `scripts/run_ensemble.py` (new) — demonstration CLI

Runs an ensemble of the whole-cell composite and writes an aggregate:
`python scripts/run_ensemble.py [--n 8] [--hours 1.0] [--observable mass] [--out reports/ensemble_<obs>.json]`.
Builds `n` cells, runs each `hours`, aggregates the chosen observable, and writes
`{observable, n, seeds, t_index, mean, std, per_cell, final_values}` to JSON.
Prints a one-line summary (mean±std of the final value across cells) so the
cell-to-cell spread is visible. Offline, no network.

### 3. Tests (`tests/test_ensemble.py`)

- `run_ensemble(n_cells=3, duration=small)` returns 3 cells with non-empty rows and 3 distinct seeds.
- Cells are genuinely independent: for a stochastic observable (e.g. `rna_counts` total or `mass`), at least two cells' trajectories differ (variation present) — proves seeds diverge.
- `aggregate(cells, "mass")` returns mean/std/per_cell aligned to the common length; `std` is non-negative and >0 somewhere when cells differ.
- `final_values(cells, "mass")` returns one value per cell.
- Missing-observable: `aggregate(cells, "does_not_exist")` raises a clear error.

## Edge cases

- Cells with unequal row counts (a cell that divided earlier) → align to min length in `aggregate`; `final_values` uses each cell's own last row.
- A non-scalar observable (map/list) → `aggregate` documents it handles scalars only; a caller wanting a map aggregates a derived scalar (e.g. sum) first. `run_ensemble` returns raw rows so map observables remain accessible.
- `n_cells=1` → runs one cell, std all zeros (degenerate but valid).

## Out of scope (v1)

- Parallel/Ray execution (sequential is fine for modest n; the runner is structured so a parallel backend could drop in later).
- A dashboard ensemble study/viz (the CLI writes JSON; a figure study consuming it is a follow-up).
- Independent-founder lineages / multi-generation ensembles.

## Fidelity impact

Adds the **ensemble capability** Karr uses to characterize single-cell variation:
independent seeded cells run and aggregated to mean±spread + end-of-cycle
distributions. The mechanisms were already stochastic; this exercises them as a
population. No model/behavior change to the single-cell composite. `FIDELITY_GAPS.md`
gap #7 marked done (ensemble infra; a dashboard variation-figure study staged).
