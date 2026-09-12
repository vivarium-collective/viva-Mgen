"""Workspace core builder for viva-Mgen.

The dashboard and every ``studies/*/sims/run.py`` call ``build_core()`` to get a
core with viva-Mgen's own Process classes registered. Editable installs are not
always auto-discovered by ``allocate_core()``'s distribution-keyed discovery, so
``local:MetabolismFbaReproductionProcess`` (etc.) are registered explicitly
here — see viva-superpowers docs/conventions/discovery.md and the viva-fenics
precedent.
"""

from __future__ import annotations

from process_bigraph import allocate_core

from .processes import (
    MetabolismFbaReproductionProcess,
    MassGrowthReproductionProcess,
    TranscriptionReproductionProcess,
    TranslationReproductionProcess,
    RnaDecayReproductionProcess,
    ProteinDecayReproductionProcess,
    ReplicationReproductionProcess,
)

_PROCESSES = (
    ("MetabolismFbaReproductionProcess", MetabolismFbaReproductionProcess),
    ("MassGrowthReproductionProcess", MassGrowthReproductionProcess),
    ("TranscriptionReproductionProcess", TranscriptionReproductionProcess),
    ("TranslationReproductionProcess", TranslationReproductionProcess),
    ("RnaDecayReproductionProcess", RnaDecayReproductionProcess),
    ("ProteinDecayReproductionProcess", ProteinDecayReproductionProcess),
    ("ReplicationReproductionProcess", ReplicationReproductionProcess),
)


def register_processes(core):
    """Register viva-Mgen's own Process classes into ``core``."""
    for name, cls in _PROCESSES:
        core.register_link(name, cls)
    return core


def build_core(core=None):
    """Return a core with viva-Mgen's processes registered."""
    if core is None:
        core = allocate_core()
    register_processes(core)
    return core
