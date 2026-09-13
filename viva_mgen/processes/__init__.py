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
from .parameter_calculator import ParameterCalculatorReproductionProcess

_CORE_MODS = None


def all_process_classes() -> dict:
    """Return ``{class_name: class}`` for every ``*ReproductionProcess`` defined
    across the process submodules (core + migrated)."""
    from . import metabolism, mass, transcription, translation, decay, replication
    from . import parameter_calculator
    mods = [metabolism, mass, transcription, translation, decay, replication,
            dna, rna, protein, cytokinesis, chromosome, parameter_calculator]
    out = {}
    for mod in mods:
        for name, obj in inspect.getmembers(mod, inspect.isclass):
            if (name.endswith("ReproductionProcess") and issubclass(obj, Process)
                    and obj.__module__ == mod.__name__):
                out[name] = obj
    return out


def _condense_card_descriptions() -> None:
    """Keep loom node cards readable.

    The loom renders a process's ``description`` on its node card with
    ``white-space: nowrap`` and a fixed card width, so a multi-paragraph
    description spills out of the card and over its neighbours. Each process is
    authored with a rich, multi-line description (summary + mechanism + contract
    + fidelity); we preserve that verbatim as ``full_description`` (still the
    class/module docstrings for developers) and expose only its concise first
    line as the card ``description`` that ``Edge.describe()`` — and therefore the
    loom — surfaces. Idempotent: only condenses a class's own multi-line
    ``description``.
    """
    for cls in all_process_classes().values():
        desc = cls.__dict__.get("description")
        if not isinstance(desc, str) or "\n" not in desc:
            continue
        cls.full_description = desc
        cls.description = next((ln.strip() for ln in desc.splitlines() if ln.strip()), desc)


_condense_card_descriptions()

__all__ = sorted(all_process_classes().keys())
