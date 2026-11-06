"""Independent molecular descriptor implementations.

These modules evaluate specified mathematical models. They do not imply that
force-field parameters or polymer properties have been experimentally validated.
"""

from .torsion import FourierTerm, TorsionPotential

__all__ = ["FourierTerm", "TorsionPotential"]
