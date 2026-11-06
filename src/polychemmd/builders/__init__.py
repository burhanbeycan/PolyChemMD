"""Independent molecular graph builders; no MD/QM engine is imported."""

from .linear import LinearChainResult, build_linear

__all__ = ["LinearChainResult", "build_linear"]
