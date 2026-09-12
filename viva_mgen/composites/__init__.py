"""viva-Mgen composites: one reusable Mycoplasma genitalium whole-cell composite."""

from . import mgen  # noqa: F401  (fires @composite_generator registration)
from .mgen import mycoplasma_genitalium, build_mgen

__all__ = ["mycoplasma_genitalium", "build_mgen"]
