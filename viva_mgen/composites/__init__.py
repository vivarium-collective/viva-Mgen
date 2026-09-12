"""viva-mGen composite generators (one per paper figure)."""

from . import whole_cell  # noqa: F401
from . import cell_cycle  # noqa: F401
from . import genetics  # noqa: F401

from .whole_cell import fig1_architecture, fig2_growth, fig3_expression, fig5_energy
from .cell_cycle import fig4_cell_cycle
from .genetics import fig6_gene_essentiality, fig7_kinetic_parameters

__all__ = [
    "fig1_architecture", "fig2_growth", "fig3_expression", "fig5_energy",
    "fig4_cell_cycle", "fig6_gene_essentiality", "fig7_kinetic_parameters",
]
