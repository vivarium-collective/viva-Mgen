"""Workspace core builder for viva-Mgen.

The dashboard and every ``studies/*/sims/run.py`` call ``build_core()`` to get a
core with viva-Mgen's own Process classes registered. Editable installs are not
always auto-discovered by ``allocate_core()``'s distribution-keyed discovery, so
every ``*ReproductionProcess`` (the 28 Karr-2012 submodel reproductions across
``viva_mgen.processes.*``) is registered explicitly here — see
viva-superpowers docs/conventions/discovery.md and the viva-fenics precedent.
"""

from __future__ import annotations

from process_bigraph import allocate_core

from .processes import all_process_classes


def register_processes(core):
    """Register every viva-Mgen ``*ReproductionProcess`` class into ``core``."""
    for name, cls in all_process_classes().items():
        core.register_link(name, cls)
    return core


def build_core(core=None):
    """Return a core with all viva-Mgen submodel processes registered."""
    if core is None:
        core = allocate_core()
    register_processes(core)
    return core
