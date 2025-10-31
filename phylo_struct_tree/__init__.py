"""
Phylogenetic Structure Trees

A package for creating phylogenetic trees based on protein structure alignment.
"""

__version__ = "0.1.0"

from .structure_alignment import StructureAligner
from .distance_matrix import DistanceMatrixBuilder
from .tree_builder import TreeBuilder

__all__ = ["StructureAligner", "DistanceMatrixBuilder", "TreeBuilder"]
