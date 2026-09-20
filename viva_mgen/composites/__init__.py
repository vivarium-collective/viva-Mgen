"""viva-Mgen composites: one reusable Mycoplasma genitalium whole-cell composite."""

from . import mgen  # noqa: F401  (fires @composite_generator registration)
from . import structural  # noqa: F401  (fires @composite_generator registration)
from .mgen import mycoplasma_genitalium, build_mgen
from .structural import mgen_structural, build_mgen_structural

__all__ = ["mycoplasma_genitalium", "build_mgen",
           "mgen_structural", "build_mgen_structural"]
