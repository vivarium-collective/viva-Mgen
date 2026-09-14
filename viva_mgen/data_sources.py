"""Workspace data-source provider for the dashboard Resources tab.

Declared in ``workspace.yaml`` as
``dashboard.data_sources.provider: viva_mgen.data_sources:list_sources``. The
workbench calls ``list_sources()`` (no args) in the env worker and renders the
returned rows in the Resources tab, so every model parameter / knowledge-base
file — plus the source paper — is browsable there.

Each row supplies a ``url``: an external hyperlink that works BOTH in the live
dashboard and in the published static snapshot. The read-only bundle on GitHub
Pages does not stage input binaries, and a row without a ``url`` falls back to
the server-only ``/api/data-source-file`` endpoint, which 404s on Pages — so we
point every file at its committed copy in the GitHub source repo (the ``raw``
endpoint) and every paper at its canonical open-access page.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

# Fallback repo coordinates (used when the git remote can't be read in the env
# worker); the live remote is parsed first so a fork/rename keeps working.
_FALLBACK_OWNER_REPO = "vivarium-collective/viva-Mgen"
_BRANCH = "main"

# (key, path relative to repo root, category, optional note)
_FILES = [
    ("mgen-parameters", "datasets/parameters.csv", "parameters"),
    ("fitted-constants", "datasets/fitted_constants.json", "parameters"),
    ("karr-kinetic-parameters", "datasets/karr_parameters.json", "parameters",
     "Karr 2012 kinetic parameters (knowledge base constants) — the real "
     "M. genitalium whole-cell model per-process kinetic constants and initial "
     "states, transcribed verbatim from the reference repo's data/parameters.json "
     "and consumed via viva_mgen.kb.load_karr_parameters()."),
    ("metabolic-demand", "datasets/karr_metabolic_demand.json", "parameters",
     "Aggregate NMP + amino-acid demand the fitted transcriptome/proteome imply — "
     "the metabolic requirement the ParCa hands to metabolism (Karr FitConstants "
     "forward coupling). Computed from the decoded genome + per-gene synthesis rates; "
     "M. genitalium's AT-rich genome shows in the ~62%% A+U ribonucleotide demand."),
    ("complex-stoichiometry", "datasets/karr_complexes.json", "knowledge-base",
     "Real macromolecular-complex subunit stoichiometry for all 201 M. genitalium "
     "protein complexes, decoded natively from the Karr 2012 knowledge-base MCOS .mat "
     "(see viva_mgen.kb_decode.decode_protein_complexes + scripts/extract_kb_complexes.py): "
     "each complex's true integer subunit counts (e.g. DNA_GYRASE = 2 GyrB + 2 GyrA; the "
     "30S ribosomal subunit = 20 r-proteins + 16S rRNA; 50S = 32 r-proteins + 23S + 5S). "
     "Consumed by the complexation and ribosome-assembly submodels."),
    ("tf-regulation", "datasets/karr_tf_regulation.json", "knowledge-base",
     "Real transcription-factor regulatory network decoded from the KB TranscriptionUnit "
     "objects (viva_mgen.kb_decode.decode_transcription_regulation): 52 genes regulated by "
     "the 5 genuine M. genitalium TFs (MG_127/MG_236/MG_101 monomers, MG_205/MG_428 dimers) "
     "with their true per-edge activity fold-changes (>1 activates, <1 represses). Consumed "
     "by the TranscriptionalRegulation submodel."),
    ("protein-maturation", "datasets/karr_protein_maturation.json", "knowledge-base",
     "Per-protein maturation classification decoded from the Karr 2012 KB: "
     "signal-sequence type (14 lipoproteins, 20 secretory), N-terminal Met-cleavage "
     "flag (35 proteins), signal length. Routes the protein-maturation submodels to "
     "the correct protein subsets."),
    ("ips189-network", "datasets/ips189.sbml.xml", "knowledge-base"),
    ("gene-list", "datasets/genes.csv", "knowledge-base"),
    ("observed-gene-expression", "datasets/karr_gene_expression.csv", "knowledge-base",
     "Genuine per-gene knowledge-base parameters for all 525 M. genitalium genes, "
     "decoded natively from the Karr 2012 knowledge-base MCOS .mat (no MATLAB/Octave; "
     "see viva_mgen.kb_decode + scripts/extract_kb_genes.py): RNA type, genome "
     "coordinates and length, transcription direction, mRNA half-life, Karr's fitted "
     "synthesis rate, and the Weiner et al. 2003 observed expression profile "
     "(32/37/43 C). Row order matches genes.csv; consumed by the ParCa "
     "(viva_mgen.parca)."),
    ("karr-2012-supplement", "workspace/references/papers/Karr2012_supplementary_information.pdf", "reference"),
]

# External references (no repo file — linked to their canonical open-access page)
# (key, url, category, note)
_REFERENCES = [
    ("karr-2012-paper",
     "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413483/",
     "reference",
     "Karr JR, Sanghvi JC, Macklin DN, et al. A whole-cell computational model "
     "predicts phenotype from genotype. Cell 150(2):389-401 (2012). "
     "doi:10.1016/j.cell.2012.05.044 — the model this workspace reproduces."),
    ("weiner-2003-expression",
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC275481/",
     "reference",
     "Weiner J III, Zimmerman C-U, Göhlmann HWH, Herrmann R. Transcription "
     "profiles of the bacterium Mycoplasma pneumoniae grown at different "
     "temperatures. Nucleic Acids Res 31(21):6306-20 (2003). doi:10.1093/nar/gkg841 "
     "— the source microarray transcription profiles behind the observed gene "
     "expression (see the observed-gene-expression dataset)."),
]


def _raw_base(root: Path) -> str:
    """``owner/repo`` from the git remote, falling back to the known coordinates
    when git is unavailable (e.g. in the env worker)."""
    owner_repo = _FALLBACK_OWNER_REPO
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "config", "--get", "remote.origin.url"],
            capture_output=True, text=True, timeout=5,
        ).stdout.strip()
        if out:
            slug = out.rsplit("github.com", 1)[-1].lstrip(":/")
            slug = slug[:-4] if slug.endswith(".git") else slug
            if slug.count("/") == 1:
                owner_repo = slug
    except Exception:  # noqa: BLE001 — git absent in the env worker is fine
        pass
    return owner_repo


def list_sources() -> list:
    root = Path(__file__).resolve().parent.parent
    owner_repo = _raw_base(root)
    # "open ↗" points at GitHub's BLOB view (syntax-highlighted / CSV-as-table /
    # line numbers + a Download-raw button) rather than the raw endpoint, which
    # just dumps unusable text (e.g. a 471 KB SBML XML) into the browser.
    blob = f"https://github.com/{owner_repo}/blob/{_BRANCH}"
    rows = []
    for entry in _FILES:
        key, rel, category = entry[0], entry[1], entry[2]
        note = entry[3] if len(entry) > 3 else None
        p = root / rel
        row = {
            "key": key,
            "path": rel,
            "category": category,
            "kind": "file",
            "size_bytes": (p.stat().st_size if p.exists() else 0),
            "url": f"{blob}/{rel}",
        }
        if note:
            row["note"] = note
        rows.append(row)
    for key, url, category, note in _REFERENCES:
        rows.append({
            "key": key,
            "path": "",
            "category": category,
            "kind": "reference",
            "size_bytes": 0,
            "url": url,
            "note": note,
        })
    return rows
