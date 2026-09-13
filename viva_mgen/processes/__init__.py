"""viva-Mgen process submodels — clean-room reproductions of the Karr 2012
Mycoplasma genitalium whole-cell model's 28 cellular-process submodels.

Every process class is a ``*ReproductionProcess`` (a from-scratch reduced
reimplementation, not the original MATLAB). ``all_process_classes()`` collects
them for registration in :mod:`viva_mgen.core`.
"""

from __future__ import annotations

import inspect

from process_bigraph import Process

# core (original) submodels
from .metabolism import MetabolismFbaReproductionProcess
from .mass import MassGrowthReproductionProcess
from .transcription import TranscriptionReproductionProcess
from .translation import TranslationReproductionProcess
from .decay import RnaDecayReproductionProcess, ProteinDecayReproductionProcess
from .replication import ReplicationReproductionProcess
# migrated submodels (one module per category)
from . import dna, rna, protein, cytokinesis, chromosome

_CORE_MODS = None


def all_process_classes() -> dict:
    """Return ``{class_name: class}`` for every ``*ReproductionProcess`` defined
    across the process submodules (core + migrated)."""
    from . import metabolism, mass, transcription, translation, decay, replication
    mods = [metabolism, mass, transcription, translation, decay, replication,
            dna, rna, protein, cytokinesis, chromosome]
    out = {}
    for mod in mods:
        for name, obj in inspect.getmembers(mod, inspect.isclass):
            if (name.endswith("ReproductionProcess") and issubclass(obj, Process)
                    and obj.__module__ == mod.__name__):
                out[name] = obj
    return out


__all__ = sorted(all_process_classes().keys())
