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
from .allocation import AllocatorProcess  # noqa: F401

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
    from . import allocation
    out["AllocatorProcess"] = allocation.AllocatorProcess
    return out


def _alias_full_descriptions() -> None:
    """Expose each process's full multi-line ``description`` verbatim.

    The loom's contract renderer (bigraph_schema.contract.resolve_contract →
    ProcessContract.from_description) already splits a multi-line ``description``
    into a one-line SUMMARY (shown on the card), the governing EQUATIONS (rendered
    as KaTeX at the contract zoom tier), and the full prose (revealed at the full
    zoom tier / on click). So the whole authored description — mechanism, the
    ``Contract — in:/out:`` line, and the fidelity note — must reach the loom
    intact; we do NOT truncate it. ``full_description`` is kept as an alias for
    any developer/tooling that wants the verbatim text in one field.
    """
    for cls in all_process_classes().values():
        desc = cls.__dict__.get("description")
        if isinstance(desc, str) and desc.strip():
            cls.full_description = desc


_alias_full_descriptions()

__all__ = sorted(all_process_classes().keys())
