"""viva-mGen: a viva-native (process-bigraph) clean-room reproduction of the
Karr et al. 2012 *Mycoplasma genitalium* whole-cell model.

Every process class is a ``*ReproductionProcess`` — a from-scratch Python
reimplementation of the corresponding Karr 2012 submodel, not the original
MATLAB and not a bridge to it. The metabolic FBA network is the published
Suthers 2009 iPS189 reconstruction (see :mod:`viva_mgen.kb`).
"""

from .processes import (
    MetabolismFbaReproductionProcess,
    MassGrowthReproductionProcess,
    TranscriptionReproductionProcess,
    TranslationReproductionProcess,
    RnaDecayReproductionProcess,
    ProteinDecayReproductionProcess,
    ReplicationReproductionProcess,
)
from .core import build_core, register_processes
from . import composites  # noqa: F401  (fires @composite_generator registration)

__all__ = [
    "MetabolismFbaReproductionProcess",
    "MassGrowthReproductionProcess",
    "TranscriptionReproductionProcess",
    "TranslationReproductionProcess",
    "RnaDecayReproductionProcess",
    "ProteinDecayReproductionProcess",
    "ReplicationReproductionProcess",
    "build_core",
    "register_processes",
]
