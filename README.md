# viva-mGen

<!-- BEGIN:dashboard -->
<!-- END:dashboard -->

A **viva-native (process-bigraph) clean-room reproduction** of the Karr et al.
2012 *Cell* whole-cell model of *Mycoplasma genitalium*
([A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://www.cell.com/fulltext/S0092-8674(12)00776-3)).
Each cellular submodel is re-expressed as a native process-bigraph `Process`,
composed through shared cell-variable stores, and exercised by a showcase
**investigation with one study per figure** of the paper.

> **This is a reproduction, not the original.** Every process class is named
> `*ReproductionProcess` and is a from-scratch Python reimplementation of the
> corresponding Karr 2012 submodel — it is **not** the original MATLAB code and
> does **not** bridge to it. Fidelity is honest and per-study (see below):
> some figures are reproduced at full fidelity, others at reduced/representative
> fidelity, and each study says which.

## What it reproduces

The original whole-cell model integrates 28 submodels over 16 cell variables,
requires a MySQL knowledge base, and runs ~8–10 h per cell on a cluster. This
reproduction implements a **reduced-but-genuine core** of that model in
process-bigraph:

| Submodel (reproduction) | Reproduces | Fidelity |
|---|---|---|
| `MetabolismFbaReproductionProcess` | FBA metabolism over the published Suthers 2009 **iPS189** network | full (the real reconstruction the WCM used) |
| `MassGrowthReproductionProcess` | mass accumulation, dry-mass composition, ~9 h doubling | full (fitted constants) |
| `TranscriptionReproductionProcess` | stochastic mRNA synthesis | reduced (representative gene panel) |
| `TranslationReproductionProcess` | stochastic protein synthesis | reduced |
| `RnaDecayReproductionProcess` / `ProteinDecayReproductionProcess` | Poisson decay | reduced |
| `ReplicationReproductionProcess` | dNTP-buffered three-phase cell cycle | reduced mechanism (genuine dynamics) |

### Data provenance

* **Metabolic network** — Suthers et al. 2009 *M. genitalium* reconstruction
  **iPS189** (BioModels `MODEL1507180052`, `datasets/ips189.sbml.xml`; 351
  reactions, 126 genes), the model the WCM's metabolism submodel was built on.
* **Gene list + reference essentiality** — extracted from the Karr 2012 repo
  (`datasets/genes.csv`, 525 genes).
* **Fitted constants** — mass fractions, cell-cycle length (32400 s ≈ 9 h),
  maintenance and exchange bounds, taken verbatim from the repo's
  `parameters.json` (`viva_mgen/constants.py`).

The original knowledge base's per-species / kinetic tables are locked inside
undecodable MATLAB MCOS objects; the iPS189 network and the paper stand in for
those.

## The investigation: one study per figure

`investigations/mgen-wholecell-showcase/` binds seven studies, one per figure:

| Study | Figure | What it shows |
|---|---|---|
| `fig1-architecture` | Fig 1 | the integrated cell as composed pbg processes over shared stores |
| `fig2-growth-and-composition` | Fig 2 | ~9 h doubling, dry-mass composition, growth curve |
| `fig3-gene-expression` | Fig 3 | bursty low-copy mRNA, accumulating protein |
| `fig4-cell-cycle-regulation` | Fig 4 | emergent inverse initiation↔replication duration control |
| `fig5-energy-allocation` | Fig 5 | ATP/GTP production and expression energy budget |
| `fig6-gene-essentiality` | Fig 6 | in-silico single-gene KO vs reference essentiality |
| `fig7-kinetic-parameters` | Fig 7 | growth vs kinetic-parameter (kcat proxy) sweeps |

## Installation

    # development (editable)
    uv venv .venv && source .venv/bin/activate
    uv pip install -e ".[dev]"
    pytest

Once installed, processes register via `viva_mgen.core.build_core()`.

## Quick start

```python
from process_bigraph import Composite, gather_emitter_results
from viva_mgen.core import build_core
from viva_mgen.composites import fig2_growth

core = build_core()
sim = Composite({"state": fig2_growth(core)}, core=core)
sim.run(600.0)
print(gather_emitter_results(sim)[("emitter",)][-1])
```

Run a study locally:

    python studies/fig2-growth-and-composition/sims/run.py

## API reference

| Class | Ports (in → out) |
|---|---|
| `MetabolismFbaReproductionProcess` | `nutrient_scale` → `growth_rate`, `growth_fraction`, `atp_production`, `gtp_production`, `feasible` |
| `MassGrowthReproductionProcess` | `growth_fraction`, `mass` → `mass`, `mass_fractions`, `volume`, `division` |
| `TranscriptionReproductionProcess` | `ntp`, `rna_pol` → `rna_counts`, `ntp` |
| `TranslationReproductionProcess` | `rna_counts`, `gtp` → `protein_counts`, `gtp` |
| `RnaDecayReproductionProcess` | `rna_counts` → `rna_counts` |
| `ProteinDecayReproductionProcess` | `protein_counts` → `protein_counts` |
| `ReplicationReproductionProcess` | `dntp_synthesis_scale` → replication phase/fraction/dNTP/durations |

## Limitations

* Reduced core, not the full 28-submodel model: chromosome occupancy (Fig 3
  original), 128-cell populations (Fig 3/4/5 distributions), and the full
  3,011-strain deletion scan (Fig 6/7) are represented at reduced fidelity —
  each study states its scope.
* Expression uses a representative gene panel, not all ~480 genes.
* Essentiality (Fig 6) is computed for the ~125 metabolic genes iPS189 covers
  (the paper notes metabolic genes are the most debilitating).

## License

MIT. The original whole-cell model and iPS189 are likewise MIT / CC-licensed;
see their upstreams.
