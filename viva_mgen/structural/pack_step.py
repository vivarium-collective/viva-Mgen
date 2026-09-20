"""``MgenStructuralStep`` — pack the structural cell (Task 4's
``build_mgen_pack``) as a process-bigraph ``Step``.

A one-shot, stateless computation (assemble ingredients + geometry, invoke
``pbg_parsimony``'s Rust packer, write the pack to disk) rather than a
per-tick dynamical submodel, so this wraps ``Step`` (see
``ReferenceDataStep`` in viva-biomodels for the sibling-workspace precedent)
rather than the ``Process`` base the 28 Karr-2012 submodel reproductions in
``viva_mgen.processes.*`` use (e.g. ``viva_mgen/processes/mass.py``). It
follows the same conventions as those processes where they apply: a class
``config_schema``, an ``inputs()``/``outputs()`` port-schema pair, and
registration into a process-bigraph ``core`` via ``core.register_link``
(mirroring ``viva_mgen.core.register_processes``) under both the bare class
name and the dotted import path.
"""

from __future__ import annotations

from typing import Any, ClassVar, Dict

from process_bigraph import Step

from viva_mgen.structural.build import build_mgen_pack
from viva_mgen.structural.counts import maritan_counts, mgen_sim_counts
from viva_mgen.structural.maritan_tables import load_genes, load_proteins


class MgenStructuralStep(Step):
    """Pack the mgen structural cell and emit the pack path + sidecar.

    Config
    ------
    counts_source : ``"maritan"`` | ``"sim"``
        ``"maritan"`` (default) uses the baseline (WC-MG-flavoured)
        expression-derived copy numbers (``maritan_counts``, Task 3).
        ``"sim"`` uses this repo's own simulation output
        (``mgen_sim_counts``) — the ``protein_counts`` input state (an
        already-accumulated ``{gene_id: count}`` mapping, or a completed
        run directory path) supplies the counts.
    top_n : int | None
        Only the ``top_n`` highest-count ingredients are placed (``None``
        places every ingredient with a resolvable structure/MW basis).
    out_dir : str
        Directory the pack is written to.
    name : str
        Pack basename (``build_mgen_pack``'s ``name=``).

    Inputs
    ------
    protein_counts : map[float]
        Per-gene protein copy numbers, read only when
        ``counts_source == "sim"`` (ignored for ``"maritan"``, which derives
        its own counts from the committed Maritan tables).

    Outputs
    -------
    pack_path : overwrite[string]
        Path to the packed cell.
    sidecar_path : overwrite[string]
        Path to the pack's sidecar (per-ingredient display metadata).
    n_placed : overwrite[float]
        Number of ingredients actually placed in the pack.
    """

    description = (
        "Pack the mgen structural cell (ingredient roster + capsule + chromosome) via "
        "pbg_parsimony, wrapping Task 4's build_mgen_pack as a process-bigraph Step. "
        "counts_source selects the copy-number provider: 'maritan' (baseline, "
        "expression-derived WC-MG-flavoured counts) or 'sim' (this repo's own simulation "
        "protein_counts).\n"
        "Contract — in: protein_counts (map[float], read only for counts_source='sim'). "
        "out: pack_path, sidecar_path (overwrite[string]), n_placed (overwrite[float])."
    )

    config_schema: ClassVar[Dict[str, Any]] = {
        "counts_source": {"_type": "string", "_default": "maritan"},
        "top_n": {"_type": "maybe[integer]", "_default": None},
        "out_dir": {"_type": "string", "_default": "out/structural"},
        "name": {"_type": "string", "_default": "mgen"},
    }

    def inputs(self) -> Dict[str, str]:
        return {"protein_counts": "map[float]"}

    def outputs(self) -> Dict[str, str]:
        return {
            "pack_path": "overwrite[string]",
            "sidecar_path": "overwrite[string]",
            "n_placed": "overwrite[float]",
        }

    def _counts(self, state: Dict[str, Any]) -> Dict[str, int]:
        source = self.config.get("counts_source", "maritan")
        proteins = load_proteins()
        genes = load_genes()
        if source == "sim":
            protein_counts = (state or {}).get("protein_counts") or {}
            return mgen_sim_counts(protein_counts, proteins, genes)
        if source == "maritan":
            return maritan_counts(proteins, genes)
        raise ValueError(
            f"MgenStructuralStep: unknown counts_source {source!r}; expected "
            "'maritan' or 'sim'"
        )

    def update(self, state: Dict[str, Any]) -> Dict[str, Any]:
        counts = self._counts(state)
        result = build_mgen_pack(
            counts,
            out_dir=self.config["out_dir"],
            top_n=self.config.get("top_n"),
            name=self.config.get("name", "mgen"),
        )
        return {
            "pack_path": str(result["pack_path"]),
            "sidecar_path": str(result["sidecar_path"]),
            "n_placed": float(result["n_placed"]),
        }


def register_structural(core):
    """Register ``MgenStructuralStep`` into ``core`` under its bare name and
    dotted import path, mirroring ``viva_mgen.core.register_processes``."""
    core.register_link("MgenStructuralStep", MgenStructuralStep)
    core.register_link(
        f"{MgenStructuralStep.__module__}.MgenStructuralStep", MgenStructuralStep
    )
    return core
