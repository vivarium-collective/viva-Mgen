"""Structural-model composite: a one-shot 3D whole-cell pack of M. genitalium.

Unlike the time-series ``mycoplasma_genitalium`` composite (28 submodels
integrated over the cell cycle), this composite is a single
:class:`~viva_mgen.structural.pack_step.MgenStructuralStep` that turns the
Maritan-2022 ingredient roster + a copy-number source into a packed 3D cell
(a ``pack.json`` + sidecar) via the imported ``pbg_parsimony`` engine. It is
the composite the ``structural-model`` investigation's studies drive:

* ``counts_source="maritan"`` — the faithful Maritan/WC-MG baseline.
* ``counts_source="sim"`` — abundances from a viva_mgen simulated run (variant).

See ``docs/superpowers/specs/2026-09-19-structural-model-investigation-design.md``.
"""
from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

from ..structural.pack_step import MgenStructuralStep as _Step


def build_mgen_structural(core=None, *, counts_source: str = "maritan",
                          top_n=None, out_dir: str | None = None):
    """Return a composite doc wrapping a single ``MgenStructuralStep``.

    The step packs the cell on update and exposes ``pack_path`` / ``n_placed`` /
    ``sidecar_path``. A RAMEmitter captures the scalar ``n_placed`` so the run is
    recordable; the 3D pack itself is the ``pack.json`` the viewer reads.
    """
    if core is None:
        from ..core import build_core
        core = build_core()
    # Register the structural step type on this core.
    from ..structural.pack_step import register_structural
    register_structural(core)

    config = {"counts_source": counts_source}
    if top_n is not None:
        config["top_n"] = top_n
    if out_dir is not None:
        config["out_dir"] = out_dir

    ports = list(_Step(config, core=core).outputs())
    doc = {
        "structural": {
            "_type": "step",
            "address": f"local:{_Step.__module__}.{_Step.__name__}",
            "config": config,
            "inputs": {},
            "outputs": {p: ["structural_out", p] for p in ports},
        },
        "structural_out": {p: "" for p in ports},
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {"structural_out": {"n_placed": "integer"}}},
            "inputs": {"structural_out": ["structural_out"]},
        },
    }
    return doc


@composite_generator(
    name="mgen_structural",
    description="One-shot 3D whole-cell structural pack of M. genitalium — the Maritan-2022 ingredient roster (S1/S2) resolved to PDB/AlphaFold structures and packed at true abundance into a single-membrane capsule via the pbg_parsimony (parsimony) engine. counts_source selects the faithful Maritan/WC-MG baseline vs a viva_mgen simulated variant.",
    parameters={
        "counts_source": {"type": "string", "default": "maritan",
                          "description": "Copy-number source: 'maritan' (WC-MG baseline) | 'sim' (viva_mgen run)"},
        "top_n": {"type": "integer", "default": 0,
                  "description": "Place only the top-N species by count (0/unset = every ingredient)"},
    },
    emitters=[{"address": "local:RAMEmitter", "paths": ["structural_out/n_placed"]}],
)
def mgen_structural(core=None, *, counts_source: str = "maritan", top_n=0):
    return build_mgen_structural(core, counts_source=counts_source,
                                 top_n=(top_n or None))
