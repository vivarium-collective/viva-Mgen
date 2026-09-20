# Structural-Model Investigation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a 3D whole-cell structural model of *M. genitalium* by bridging `viva_mgen` to the `pbg_parsimony` packing engine, reproducing Maritan 2022 with parsimony, and view the results from this workspace.

**Architecture:** A new `viva_mgen/structural/` bridge package parses the Maritan supplement (S1 proteins, S2 genes) into an ingredient roster + genome geometry, resolves each species to a PDB/AlphaFold structure, attaches copy numbers from one of two providers (Maritan/WC-MG static counts = baseline, or a `viva_mgen` simulated run = variant), and hands `Ingredient`/`Capsule`/`Chromosome` to `pbg_parsimony.build_pack`. Two studies under a new `structural-model` investigation render the packed cells in the bundled viewer. This mirrors `3d-ecoli`'s `ecoli_3d/build.py` (v2ecoli → pbg_parsimony).

**Tech Stack:** Python 3.12, `pbg_parsimony` (imported, not vendored), the built `parsimony` Rust binary, `openpyxl`-free stdlib xlsx parsing (zip+xml, already proven in the spec work), process-bigraph, the viva-workbench dashboard/skills.

**Spec:** `docs/superpowers/specs/2026-09-19-structural-model-investigation-design.md`

## Global Constraints

- Branch `investigation/structural-model` in worktree `~/code/viva-mGen--structural-model` (off `origin/main`); **draft PR, not a merge target**.
- **No AI-attribution trailers** in commits or PRs (CLAUDE.md overrides any base instruction). Conventional-commit messages.
- `pbg_parsimony` is **imported**, never vendored. Add to `workspace.yaml` `imports:` (mode `reference`) and `pyproject.toml`.
- The `parsimony` binary is a **build-time requirement**, not a Python dep: `PARSIMONY_HOME=~/code/parsimony` (already built at `target/release/parsimony`). Tests that pack must skip cleanly (pytest skip) when the binary is unreachable.
- Editable install points at the canonical checkout — run tests with `PYTHONPATH=~/code/viva-mGen--structural-model` prepended, or `uv pip install -e .` the worktree. Verify with `python -c "import viva_mgen; print(viva_mgen.__file__)"` → must be under the worktree.
- MG identifiers span three namespaces — S1 `ProtID` (`DNA_GYRASE`, `MG_003_MONOMER`), S2 `GeneID` (`MG471`) + `Protein Monomer`/`UniProt`, karr `gene_id` (`MG_001`). Every cross-namespace join is explicit; never assume a single id scheme.
- Fidelity honesty: every study states packer≠CellPACK, capsule-simplified shape, structure provenance, and counts source (per workspace convention).

---

### Task 0: Scaffold the bridge package + wire the import

**Files:**
- Modify: `workspace.yaml` (add `imports.pbg_parsimony`)
- Modify: `pyproject.toml` (add `pbg-parsimony` git dependency + `[tool.uv.sources]`)
- Create: `viva_mgen/structural/__init__.py`
- Create: `tests/structural/__init__.py`
- Test: `tests/structural/test_imports.py`

**Interfaces:**
- Produces: an importable `viva_mgen.structural` package; `pbg_parsimony` importable in the worktree env.

- [ ] **Step 1: Write the failing test**
```python
# tests/structural/test_imports.py
def test_pbg_parsimony_importable():
    import pbg_parsimony
    from pbg_parsimony import Ingredient, Capsule, Chromosome, build_pack  # noqa

def test_structural_package_importable():
    import viva_mgen.structural  # noqa
```

- [ ] **Step 2: Run — expect FAIL** (`ModuleNotFoundError: viva_mgen.structural` / `pbg_parsimony`)
Run: `PYTHONPATH=$PWD pytest tests/structural/test_imports.py -v`

- [ ] **Step 3: Wire the import.** Add to `workspace.yaml` under `imports:`:
```yaml
  pbg_parsimony:
    source: https://github.com/vivarium-collective/pbg-parsimony
    ref: main
    mode: reference
    description: |
      cellPACK-style 3D packing engine (Ingredient/Capsule/Chromosome/
      build_pack) used by viva_mgen.structural to pack the whole-cell
      molecular roster into a single-membrane M. genitalium capsule.
```
Add to `pyproject.toml` dependencies: `"pbg-parsimony"`, and under `[tool.uv.sources]`: `pbg-parsimony = { git = "https://github.com/vivarium-collective/pbg-parsimony", branch = "main" }`. Then install into the worktree env: `uv pip install -e ~/code/viva-parsimony --no-deps` (local checkout is fine; it's the same package). Create empty `viva_mgen/structural/__init__.py` with a one-line docstring and `tests/structural/__init__.py`.

- [ ] **Step 4: Run — expect PASS.** Also verify `python -c "import viva_mgen; print(viva_mgen.__file__)"` resolves under the worktree.

- [ ] **Step 5: Commit** — `git commit -m "feat(structural): scaffold viva_mgen.structural + import pbg_parsimony"`

---

### Task 1: Parse the Maritan supplement into committed data artifacts

**Files:**
- Create: `viva_mgen/structural/maritan_tables.py`
- Create: `viva_mgen/structural/data/` (output dir)
- Create: `viva_mgen/structural/parse_supplement.py` (one-shot generator)
- Test: `tests/structural/test_maritan_tables.py`

**Interfaces:**
- Produces:
  - `ProteinRow` dataclass: `prot_id:str, name:str, function:str, compartment:str, kind:str ("monomer"|"complex"), pdb_id:str|None, biosynthesis:str, seq_length:int|None`
  - `GeneRow` dataclass: `gene_id:str, gtype:str, coord:int, length:int, direction:str, tu:str, essential:bool, name:str, uniprot:str|None, protein_monomer:str|None`
  - `load_proteins() -> list[ProteinRow]` (reads `data/mgen_proteins.csv`)
  - `load_genes() -> list[GeneRow]` (reads `data/mgen_genes.csv`)

- [ ] **Step 1: Write the failing test**
```python
# tests/structural/test_maritan_tables.py
from viva_mgen.structural.maritan_tables import load_proteins, load_genes

def test_proteins_roster_size_and_known_row():
    rows = load_proteins()
    assert len(rows) > 400                      # ~996 species in S1
    gyr = next(r for r in rows if r.prot_id == "DNA_GYRASE")
    assert gyr.kind == "complex"
    assert gyr.pdb_id == "6RKW"                 # S1 "Structural Model"
    assert gyr.compartment == "c"

def test_genes_known_row_and_uniprot_join():
    genes = load_genes()
    trna = next(g for g in genes if g.gene_id == "MG471")
    assert trna.gtype.lower() == "trna"
    assert trna.essential is True
    # at least some protein-coding genes carry a UniProt accession
    assert any(g.uniprot for g in genes)
```

- [ ] **Step 2: Run — expect FAIL** (module/data missing).

- [ ] **Step 3: Implement.** Write `parse_supplement.py` that reads the three
`.xlsx` in `workspace/references/papers/Maritan2022_supplementary/` via the
stdlib zip+xml reader (see the approach already used during spec authoring:
`xl/sharedStrings.xml` + `xl/worksheets/sheetN.xml`, namespace
`spreadsheetml/2006/main`). Header rows: S1 data starts at row index 4 with the
`ProtID` header at row 3; S2 header at row 2. Emit `data/mgen_proteins.csv` and
`data/mgen_genes.csv` with the dataclass fields as columns. Normalize:
`compartment` lowercased; `kind` from the "Type (monomer/complex)" column
(`ProteinComplex`→`complex`, `ProteinMonomer`/monomer→`monomer`); `pdb_id` =
"Structural Model" cell when it looks like a PDB id (4-char alnum) else `None`;
`essential` from S2 "Essential" (`yes`→True). Then `maritan_tables.py` exposes
the dataclasses + `load_*` readers over the committed CSVs (csv stdlib). Run the
generator once and commit its CSV outputs (so runtime never touches `.xlsx`).

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit** — `git commit -m "feat(structural): parse Maritan S1/S2 supplement into committed roster + gene tables"`

---

### Task 2: Resolve each protein to a structure ref

**Files:**
- Create: `viva_mgen/structural/structures.py`
- Test: `tests/structural/test_structures.py`

**Interfaces:**
- Consumes: `ProteinRow`, `GeneRow` (Task 1); `pbg_parsimony.StructureRef`.
- Produces: `structure_ref_for(protein: ProteinRow, uniprot_by_prot: dict[str,str]) -> StructureRef | None`; `uniprot_index(genes, proteins) -> dict[str,str]` mapping `ProtID → UniProt`.

- [ ] **Step 1: Write the failing test**
```python
# tests/structural/test_structures.py
from pbg_parsimony import StructureRef
from viva_mgen.structural.maritan_tables import load_proteins, load_genes
from viva_mgen.structural.structures import structure_ref_for, uniprot_index

def test_pdb_takes_priority():
    proteins, genes = load_proteins(), load_genes()
    idx = uniprot_index(genes, proteins)
    gyr = next(r for r in proteins if r.prot_id == "DNA_GYRASE")
    ref = structure_ref_for(gyr, idx)
    assert ref.kind in ("pdb", "cif") and ref.ref == "6RKW"

def test_alphafold_fallback_when_no_pdb_but_uniprot():
    from viva_mgen.structural.maritan_tables import ProteinRow
    p = ProteinRow(prot_id="MG_XXX_MONOMER", name="x", function="f",
                   compartment="c", kind="monomer", pdb_id=None,
                   biosynthesis="", seq_length=100)
    ref = structure_ref_for(p, {"MG_XXX_MONOMER": "P47000"})
    assert ref.kind == "alphafold" and ref.ref == "P47000"

def test_none_when_no_structure_source():
    from viva_mgen.structural.maritan_tables import ProteinRow
    p = ProteinRow("NO_STRUCT", "x", "f", "c", "monomer", None, "", None)
    assert structure_ref_for(p, {}) is None
```

- [ ] **Step 2: Run — expect FAIL.**

- [ ] **Step 3: Implement.** `uniprot_index` joins S1 monomer `ProtID` → S2
`Protein Monomer`/`UniProt` (match on the MG gene id embedded in the monomer id
or S2's `Protein Monomer` column). `structure_ref_for`: if `protein.pdb_id`,
return `StructureRef(kind="pdb"|"cif", ref=pdb_id)` (cif if the id is a modern
long accession); elif a UniProt exists, `StructureRef(kind="alphafold", ref=uniprot)`; else `None`. Use `StructureRef`'s actual constructor — confirm its
field names in `pbg_parsimony/structures.py` before writing.

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit** — `git commit -m "feat(structural): resolve Maritan proteins to PDB/AlphaFold structure refs"`

---

### Task 3: Copy-number providers (baseline + variant)

**Files:**
- Create: `viva_mgen/structural/counts.py`
- Test: `tests/structural/test_counts.py`

**Interfaces:**
- Consumes: `ProteinRow` (Task 1); `datasets/karr_gene_expression.csv`, `datasets/karr_reference_values.json`; a `viva_mgen` run's `protein_counts` store.
- Produces: `maritan_counts(proteins, genes) -> dict[str,int]` (ProtID→count); `mgen_sim_counts(run_dir_or_state, proteins, genes) -> dict[str,int]`; helper `complex_count(biosynthesis, monomer_counts) -> int`.

- [ ] **Step 0 (discovery, spelled out):** Confirm the static-count source.
Run: `python -c "import json; d=json.load(open('datasets/karr_reference_values.json')); print(json.dumps(d['values'], indent=0)[:800])"`. If it contains per-gene protein copy numbers, use them directly. Otherwise derive steady-state counts from `karr_gene_expression.csv` (`synthesis_rate`, `half_life_min`) as `count ≈ synthesis_rate × half_life`, normalized to the model's total protein number — document whichever is chosen in a module docstring. This choice is baseline-only; it does not change the interface.

- [ ] **Step 1: Write the failing test**
```python
# tests/structural/test_counts.py
from viva_mgen.structural.maritan_tables import load_proteins, load_genes
from viva_mgen.structural.counts import maritan_counts, complex_count

def test_maritan_counts_cover_monomers_positive_ints():
    proteins, genes = load_proteins(), load_genes()
    counts = maritan_counts(proteins, genes)
    mono = [p for p in proteins if p.kind == "monomer"]
    covered = [p for p in mono if counts.get(p.prot_id, 0) > 0]
    assert len(covered) > 0.5 * len(mono)          # majority of monomers get a count
    assert all(isinstance(v, int) and v >= 0 for v in counts.values())

def test_complex_count_is_limiting_subunit():
    # "(2.0)MG_003_MONOMER'+(2.0)MG_004_MONOMER'" with subunit stoich 2 each
    mono = {"MG_003_MONOMER": 100, "MG_004_MONOMER": 50}
    bio = "(2.0)MG_003_MONOMER'+(2.0)MG_004_MONOMER'"
    assert complex_count(bio, mono) == 25           # min(100//2, 50//2)
```

- [ ] **Step 2: Run — expect FAIL.**

- [ ] **Step 3: Implement.** `complex_count` parses the S1 "Complex Biosynthesis"
string (`(stoich)MONOMER_ID'` terms joined by `+`) and returns
`min(monomer_counts[m] // stoich)` over subunits (0 if any subunit missing).
`maritan_counts`: assign monomer counts from the chosen static source (Step 0),
then complexes via `complex_count`. `mgen_sim_counts`: read the `protein_counts`
store (`map[float]`, keyed by MG gene id — emitted by
`viva_mgen/processes/translation.py`) from a completed run's emitter output
(xarray/parquet under the study's run dir), map gene id → S1 monomer `ProtID`
via the Task-1 join, round to int, and lift to complexes with `complex_count`.
Both return `{ProtID: int}`.

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit** — `git commit -m "feat(structural): baseline (WC-MG) + variant (sim) copy-number providers"`

---

### Task 4: `build_mgen_pack` — assemble ingredients + geometry, pack the cell

**Files:**
- Create: `viva_mgen/structural/build.py`
- Test: `tests/structural/test_build.py`

**Interfaces:**
- Consumes: Tasks 1–3; `pbg_parsimony.{Ingredient, Capsule, Chromosome, build_pack}`; `viva_mgen.constants.GENOME_LENGTH_BP` (580070).
- Produces: `build_mgen_pack(counts, *, out_dir, top_n=None, name="mgen") -> dict` (the `build_pack` return dict with `pack_path`, `n_placed`, …); `mgen_ingredients(counts, top_n=None) -> list[Ingredient]`; `mgen_capsule() -> Capsule`; `mgen_chromosome() -> Chromosome`.

- [ ] **Step 1: Write the failing test** (pure-assembly parts run without the binary; the pack step skips if unreachable)
```python
# tests/structural/test_build.py
import os, pytest
from viva_mgen.structural.build import mgen_ingredients, mgen_capsule, mgen_chromosome, build_mgen_pack
from viva_mgen.structural.counts import maritan_counts
from viva_mgen.structural.maritan_tables import load_proteins, load_genes

def _counts():
    return maritan_counts(load_proteins(), load_genes())

def test_ingredients_have_structure_or_sphere_and_region():
    ings = mgen_ingredients(_counts(), top_n=20)
    assert len(ings) == 20
    for ing in ings:
        assert ing.count > 0
        assert ing.region in ("interior", "surface", "fiber")
        assert (ing.structure is not None) or (ing.sphere_radius is not None)

def test_capsule_and_chromosome_scale():
    cap = mgen_capsule()
    assert cap.radius > 0 and cap.half_len >= cap.radius
    chrom = mgen_chromosome()
    assert chrom.n_chromosomes == 1 and chrom.beads > 1000

@pytest.mark.skipif(not (os.environ.get("PARSIMONY_HOME") or os.environ.get("PARSIMONY_BIN")),
                    reason="parsimony binary not configured")
def test_small_end_to_end_pack(tmp_path):
    res = build_mgen_pack(_counts(), out_dir=tmp_path, top_n=8, name="mgen_test")
    assert res["n_placed"] > 0
    assert (tmp_path).exists()
```

- [ ] **Step 2: Run — expect FAIL** (build module missing; end-to-end skipped if no binary).

- [ ] **Step 3: Implement.**
`mgen_capsule()` — a single `Capsule.from_volume_fl(MGEN_VOLUME_FL, radius_um=MGEN_RADIUS_UM)` with module constants `MGEN_VOLUME_FL ≈ 0.067`, `MGEN_RADIUS_UM ≈ 0.20` (M. genitalium ~0.3–0.4 µm cell; document as Maritan/literature-derived, refine against the paper's stated dimensions). No `envelope=` (single membrane).
`mgen_chromosome()` — `Chromosome(beads=…, n_chromosomes=1, genome_csv=<S2-derived genome csv path>)`; genome length from `GENOME_LENGTH_BP`; pick `beads` so contour ≈ genome bp (follow `ecoli_3d`'s `GENOME_BEADS` scaling ~ bp/135).
`mgen_ingredients(counts, top_n)` — for each `ProteinRow` with `counts[prot_id] > 0` (top-N by count if given): `Ingredient(id=prot_id, count=counts[prot_id], structure=structure_ref_for(...), region="surface" if compartment=="m" else "interior", compartment="cytoplasm", display_name=name, category=function, color=<category palette>)`. Proteins with neither PDB nor UniProt get a `sphere_radius` proxy from MW (skip if MW missing).
`build_mgen_pack(...)` — assemble the three, call `pbg_parsimony.build_pack(ingredients, capsule, chromosome, out_dir=out_dir, name=name)`, return its dict.

- [ ] **Step 4: Run — expect PASS** (with `PARSIMONY_HOME=~/code/parsimony` set, the end-to-end packs; otherwise it skips and the assembly tests pass).

- [ ] **Step 5: Commit** — `git commit -m "feat(structural): build_mgen_pack — roster+geometry -> parsimony pack"`

---

### Task 5: `mgen-structural` process-bigraph Step

**Files:**
- Create: `viva_mgen/structural/pack_step.py`
- Test: `tests/structural/test_pack_step.py`

**Interfaces:**
- Consumes: `build_mgen_pack` (Task 4); the process-bigraph `Step` base + registry used elsewhere in `viva_mgen` (follow an existing process, e.g. `viva_mgen/processes/mass.py`, for the exact base class + `register_*` pattern).
- Produces: `MgenStructuralStep` (a `Step` whose config selects the counts source `"maritan"|"sim"` and `top_n`, and whose output is the pack path/sidecar); `register_structural(core)`.

- [ ] **Step 1: Write the failing test**
```python
# tests/structural/test_pack_step.py
def test_step_registers_and_declares_io():
    from viva_mgen.structural.pack_step import MgenStructuralStep
    step = MgenStructuralStep({"counts_source": "maritan", "top_n": 8})
    outs = step.outputs()
    assert "pack_path" in outs
```

- [ ] **Step 2: Run — expect FAIL.**

- [ ] **Step 3: Implement** the Step wrapping `build_mgen_pack`, matching the base
class/signature conventions of an existing `viva_mgen` process. `counts_source
== "sim"` reads a run's `protein_counts` via `mgen_sim_counts`; `"maritan"` uses
`maritan_counts`. Provide `register_structural(core)`.

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit** — `git commit -m "feat(structural): mgen-structural process-bigraph Step + registration"`

---

### Task 6: Create the investigation + two studies, run, and wire the viewer

**Files:**
- Create (via skills): `workspace/investigations/structural-model/…`, `workspace/studies/s01-maritan-baseline/…`, `workspace/studies/s02-reproduction-driven/…`
- Create: `viva_mgen/structural/publish/` (thin build/publish script, analog of `ecoli_3d/publish/`)
- Modify: `workspace.yaml` (`visualizations:` entry pointing at the packs)

**Interfaces:**
- Consumes: everything above; the `/viva-investigation`, `/viva-study`, `/viva-viz`, `/viva-report` skills.

- [ ] **Step 1:** `/viva-investigation new structural-model` — title
`investigation: structural-model — a 3D structural model of M. genitalium with viva-parsimony`; set the question/hypothesis from the spec; **draft PR** per AGENTS.md.
- [ ] **Step 2:** `/viva-study` create `s01-maritan-baseline` — baseline composite drives `MgenStructuralStep(counts_source="maritan")`; `expected_behavior`/`behavior_tests`: species-placed count vs roster, cytoplasm volume-occupancy in Maritan's reported range, membrane proteins in the bilayer, one supercoiled chromosome. Add fidelity caveats.
- [ ] **Step 3:** Run it: `PARSIMONY_HOME=~/code/parsimony` then `/viva-study run` (or the canonical `sims/run.py`) → pack recorded; record the run.
- [ ] **Step 4:** `/viva-study` create `s02-reproduction-driven` — same roster, `counts_source="sim"` from a `viva_mgen` steady-state run; acceptance compares occupancy + per-category species counts to s01. Run it.
- [ ] **Step 5:** `/viva-viz` — register a Visualization per study pointing at the pack + the bundled `pbg_parsimony` viewer; add to `workspace.yaml` `visualizations:`. Add a `publish/` script that builds both packs reproducibly and copies them where the viewer/report reads them.
- [ ] **Step 6:** `/viva-report` (Pass A + B) to fold both studies into the investigation report; verify the packed cells are browsable from the dashboard.
- [ ] **Step 7: Commit** — study YAMLs, runs, viz wiring, publish script: `git commit -m "feat(structural): structural-model investigation — maritan-baseline + reproduction-driven studies, packed + viewable"`

---

## Self-Review

**Spec coverage:** §3 integration → Task 0; §4 bridge (maritan_tables/counts/build/pack_step) → Tasks 1–5; §5 investigation+studies → Task 6; §6 viewer/report → Task 6 Steps 5–6; §7 copy-number open dependency → Task 3 Step 0 (resolved concretely); §8 tests → each task's tests; §9 caveats → Task 6 study YAMLs; §11 PR topology → Global Constraints + Task 6 Step 1. No gaps.

**Placeholder scan:** the only deferred decision (static-count source) is a spelled-out discovery step (Task 3 Step 0) with a concrete fallback formula, not a placeholder. Two geometry constants (`MGEN_VOLUME_FL`, `MGEN_RADIUS_UM`) have concrete starting values to refine against the paper. StructureRef/Step base-class field names are flagged "confirm before writing" because they belong to imported code — deliberate, not vague.

**Type consistency:** `{ProtID: int}` counts flow from Task 3 → Task 4 `mgen_ingredients` → Task 5 Step. `ProteinRow`/`GeneRow` dataclass fields are used identically in Tasks 2–4. `structure_ref_for(protein, uniprot_by_prot)` signature matches its Task-2 definition and Task-4 call.
