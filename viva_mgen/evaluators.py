"""Workspace derived-scalar computers for viva-mGen study behavior tests.

Study run scripts (sims/run.py) compute figure observables in-process and
persist them to <workspace>/.pbg/runs.jsonl via append_run_event(..., observables=...).
These computers read those already-computed values back so a study's behavior
tests grade through the framework's derived-scalar seam — no science is
re-derived here.

Every declared derived-scalar field across the mGen studies is registered. The
field names are unique across studies, so resolving by field name alone (the
latest completed run whose observables carry it) targets the right run. A field
that a study's run has not (yet) persisted raises honestly and shows as
ungraded rather than fabricating a value; it grades automatically once the
study's run records it.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_RUN_LOG_RELPATH = ".pbg/runs.jsonl"

# fig1-architecture
FIG1_FIELDS = (
    "n_processes",
    "n_stores_coupled",
    "mass_fold_change",
    "replicated_fraction_final",
    "divides",
    "mrna_species_produced",
    "protein_species_produced",
)

# fig2-growth. The five emergent_*/single_cell_mrna_cv fields are declared but
# not currently persisted by the study's run — they grade once run.py records them.
FIG2_FIELDS = (
    "doubling_time_h",
    "final_mass_ratio",
    "protein_fraction",
    "rna_fraction",
    "dna_fold_change",
    "mean_mrna_per_gene",
    "mrna_protein_abs_corr",
    "emergent_protein_fraction",
    "emergent_dna_fraction",
    "emergent_rna_fraction",
    "single_cell_mrna_cv",
    "emergent_doubling_time_h",
)

# fig3-expression
FIG3_FIELDS = (
    "pct_explored_at_6min",
    "pct_explored_at_20min",
    "rnap_90pct_time_min",
    "rna_t50_min",
    "rna_t90_min",
    "n_collisions_per_cycle",
    "frac_collisions_by_rnap",
    "frac_collisions_displacing_smc",
    "collisions_density_pearson_r",
)

# fig4-cell-cycle
FIG4_FIELDS = (
    "cell_cycle_h",
    "init_dur_h",
    "repl_dur_h",
    "cyto_dur_h",
    "cv_repl_pct",
    "total_cv_below_phase_cvs",
    "r_init_repl",
    "r_dnaA_init",
    "r_dntp_repl",
)

# fig5-energy
FIG5_FIELDS = (
    "atp_to_gtp_ratio",
    "atp_gtp_to_redox_ratio",
    "translation_energy_share",
    "aminoacylation_energy_share",
    "transcription_energy_share",
    "unaccounted_share",
    "ntp_use_cv_across_cells",
)

# fig6-gene-essentiality
FIG6_FIELDS = (
    "essentiality_accuracy",
    "sensitivity",
    "specificity",
    "growth_distribution_bimodal",
    "pathology_metabolic_non_growing",
    "pathology_rna_ko_stops_rna",
    "pathology_protein_ko_stops_protein",
    "pathology_dna_ko_non_replicative",
    "pathology_cytokinesis_ko_non_fissive",
)

# fig7-kinetic-parameters
FIG7_FIELDS = (
    "growth_is_monotonic_in_kcat",
    "growth_saturates_at_wt",
    "growth_is_sigmoidal",
    "growth_dynamic_range",
    "growth_at_max_kcat",
    "growth_at_min_kcat",
    "kcat_half_max",
)

# parca-parameter-fitting. Declared but not yet gradeable — this study has no
# recorded run in runs.jsonl; it grades once the study is run.
PARCA_FIELDS = (
    "n_genes",
    "mrna_halflife_min_mean",
    "rrna_halflife_min",
    "nmp_au_fraction",
    "closed_loop_feasible",
)

ALL_FIELDS = (
    FIG1_FIELDS + FIG2_FIELDS + FIG3_FIELDS + FIG4_FIELDS
    + FIG5_FIELDS + FIG6_FIELDS + FIG7_FIELDS + PARCA_FIELDS
)


def _latest_observable(ws_root: Any, field: str) -> float:
    """Return `field` from the most-recently-completed run whose observables carry it.

    Reads <ws_root>/.pbg/runs.jsonl (append-only JSONL). Raises ValueError if no
    completed run recorded this observable — the evaluator surfaces that honestly
    rather than fabricating a value.
    """
    log = Path(ws_root) / _RUN_LOG_RELPATH
    if not log.is_file():
        raise ValueError(f"no run log at {log}")
    best_t = None
    best_v = None
    with log.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                ev = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not isinstance(ev, dict):
                continue
            obs = ev.get("observables")
            if not isinstance(obs, dict) or field not in obs:
                continue
            t = ev.get("completed_at") or 0.0
            if best_t is None or t >= best_t:
                best_t = t
                best_v = obs[field]
    if best_v is None:
        raise ValueError(f"observable {field!r} not found in any completed run")
    try:
        return float(best_v)
    except (TypeError, ValueError):
        raise ValueError(f"observable {field!r} is not numeric: {best_v!r}")


def register_derived_scalars(reg: dict) -> None:
    """Register a reader for every declared mGen derived-scalar field.

    field -> fn(reader, test, ws_root) -> float. Field names are unique across
    studies, so the by-field lookup in `_latest_observable` targets the right run.
    """
    for field in ALL_FIELDS:
        reg[field] = lambda reader, test, ws_root, _f=field: _latest_observable(ws_root, _f)
