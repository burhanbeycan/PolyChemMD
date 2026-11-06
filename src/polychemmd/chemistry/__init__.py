"""Explicit repeat-unit chemistry contracts; RDKit is imported only on use."""

from .contracts import (
    AtomLineage,
    ChemistryDependencyError,
    RepeatUnitSpec,
    UnsupportedChemistryError,
)

__all__ = [
    "AtomLineage", "ChemistryDependencyError", "RepeatUnitSpec",
    "UnsupportedChemistryError",
]
