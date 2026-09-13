"""viva-Mgen: a viva-native (process-bigraph) clean-room reproduction of the
Karr et al. 2012 *Mycoplasma genitalium* whole-cell model.

Every cellular-process submodel is a ``*ReproductionProcess`` — a from-scratch
reduced Python reimplementation of the corresponding Karr 2012 submodel, not the
original MATLAB and not a bridge to it. The metabolic FBA network is the
published Suthers 2009 iPS189 reconstruction (see :mod:`viva_mgen.kb`).
"""

from .processes import all_process_classes
from .core import build_core, register_processes
from . import composites  # noqa: F401  (fires @composite_generator registration)

# Re-export every submodel Process class at the package top level.
_PROCESS_CLASSES = all_process_classes()
globals().update(_PROCESS_CLASSES)

__all__ = sorted(_PROCESS_CLASSES) + ["build_core", "register_processes"]
