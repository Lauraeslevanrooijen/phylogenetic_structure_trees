"""
Phylogenetic tree construction module.

Implements distance-based phylogenetic tree construction algorithms:
- UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- Neighbor-Joining (NJ)
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Phylo.BaseTree import Tree, Clade
import io


class TreeBuilder:
    """
    Builds phylogenetic trees from distance matrices using various algorithms.
    """

    def __init__(self):
        """Initialize the TreeBuilder."""
        self.tree = None
        self.method = None

    def build_tree(self, distance_matrix: np.ndarray, labels: List[str],
                   method: str = 'nj') -> Tree:
        """
        Build a phylogenetic tree from a distance matrix.

        Args:
            distance_matrix: Square distance matrix (numpy array)
            labels: List of labels for the taxa
            method: Tree building method ('upgma' or 'nj')

        Returns:
            Bio.Phylo Tree object

        Raises:
            ValueError: If method is not supported or matrix is invalid
        """
        if method not in ['upgma', 'nj']:
            raise ValueError(f"Unsupported method: {method}. Use 'upgma' or 'nj'")

        # Validate inputs
        if len(labels) != distance_matrix.shape[0]:
            raise ValueError("Number of labels must match matrix dimensions")

        if distance_matrix.shape[0] != distance_matrix.shape[1]:
            raise ValueError("Distance matrix must be square")

        # Convert numpy array to Biopython DistanceMatrix format
        bio_matrix = self._numpy_to_bio_matrix(distance_matrix, labels)

        # Build tree using appropriate method
        constructor = DistanceTreeConstructor()

        if method == 'upgma':
            self.tree = constructor.upgma(bio_matrix)
        else:  # nj
            self.tree = constructor.nj(bio_matrix)

        self.method = method

        return self.tree

    def _numpy_to_bio_matrix(self, matrix: np.ndarray, labels: List[str]) -> DistanceMatrix:
        """
        Convert a numpy distance matrix to Biopython DistanceMatrix format.

        Args:
            matrix: Numpy distance matrix
            labels: List of labels

        Returns:
            Bio.Phylo.TreeConstruction.DistanceMatrix
        """
        n = len(labels)

        # Biopython DistanceMatrix uses lower triangular format INCLUDING diagonal
        # Row i has i+1 elements (indices 0 to i inclusive)
        # For example, with 3 taxa:
        #   Row 0: [0.0]           (just the diagonal)
        #   Row 1: [d(1,0), 0.0]   (distance to 0, and diagonal)
        #   Row 2: [d(2,0), d(2,1), 0.0]  (distances to 0, 1, and diagonal)
        matrix_list = []

        for i in range(n):
            row = []
            for j in range(i + 1):  # Include diagonal (j <= i)
                row.append(float(matrix[i, j]))
            matrix_list.append(row)

        return DistanceMatrix(labels, matrix_list)

    def build_upgma_tree(self, distance_matrix: np.ndarray, labels: List[str]) -> Tree:
        """
        Build a UPGMA tree.

        UPGMA (Unweighted Pair Group Method with Arithmetic Mean) assumes
        a constant molecular clock.

        Args:
            distance_matrix: Square distance matrix
            labels: List of labels for the taxa

        Returns:
            Bio.Phylo Tree object
        """
        return self.build_tree(distance_matrix, labels, method='upgma')

    def build_nj_tree(self, distance_matrix: np.ndarray, labels: List[str]) -> Tree:
        """
        Build a Neighbor-Joining tree.

        NJ does not assume a constant molecular clock and is generally
        preferred for evolutionary studies.

        Args:
            distance_matrix: Square distance matrix
            labels: List of labels for the taxa

        Returns:
            Bio.Phylo Tree object
        """
        return self.build_tree(distance_matrix, labels, method='nj')

    def save_tree(self, filepath: str, format: str = 'newick'):
        """
        Save the tree to a file.

        Args:
            filepath: Output file path
            format: Tree format ('newick', 'nexus', 'phyloxml', 'nexml')

        Raises:
            ValueError: If no tree has been built yet
        """
        if self.tree is None:
            raise ValueError("No tree to save. Build a tree first.")

        Phylo.write(self.tree, filepath, format)

    def get_newick_string(self) -> str:
        """
        Get the tree in Newick format as a string.

        Returns:
            Newick format string

        Raises:
            ValueError: If no tree has been built yet
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        # Use StringIO to get Newick string
        handle = io.StringIO()
        Phylo.write(self.tree, handle, 'newick')
        return handle.getvalue()

    def print_tree(self):
        """
        Print a text representation of the tree.

        Raises:
            ValueError: If no tree has been built yet
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        Phylo.draw_ascii(self.tree)

    def get_tree_info(self) -> Dict:
        """
        Get information about the tree.

        Returns:
            Dictionary with tree statistics

        Raises:
            ValueError: If no tree has been built yet
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        terminals = self.tree.get_terminals()
        nonterminals = self.tree.get_nonterminals()

        # Calculate total branch length
        total_length = self.tree.total_branch_length()

        # Get tree depth (max distance from root to terminal)
        depths = self.tree.depths()
        max_depth = max(depths.values()) if depths else 0

        return {
            'method': self.method,
            'num_terminals': len(terminals),
            'num_internal_nodes': len(nonterminals),
            'total_branch_length': total_length,
            'max_depth': max_depth,
            'is_bifurcating': self.tree.is_bifurcating()
        }

    def get_terminal_names(self) -> List[str]:
        """
        Get the names of all terminal nodes (leaves).

        Returns:
            List of terminal node names

        Raises:
            ValueError: If no tree has been built yet
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        return [terminal.name for terminal in self.tree.get_terminals()]

    def root_tree(self, outgroup: Optional[str] = None) -> Tree:
        """
        Root the tree at a specified outgroup or midpoint.

        Args:
            outgroup: Name of the outgroup terminal node.
                     If None, roots at midpoint.

        Returns:
            Rooted tree

        Raises:
            ValueError: If no tree has been built yet or outgroup not found
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        if outgroup is None:
            # Root at midpoint
            self.tree.root_at_midpoint()
        else:
            # Root at specified outgroup
            try:
                self.tree.root_with_outgroup(outgroup)
            except ValueError as e:
                raise ValueError(f"Could not root tree with outgroup '{outgroup}': {e}")

        return self.tree

    def ladderize(self, reverse: bool = False):
        """
        Ladderize the tree (sort branches by number of terminals).

        Args:
            reverse: If True, sort in reverse order

        Raises:
            ValueError: If no tree has been built yet
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        self.tree.ladderize(reverse=reverse)

    def get_distance(self, terminal1: str, terminal2: str) -> float:
        """
        Get the patristic distance between two terminal nodes.

        Args:
            terminal1: Name of first terminal
            terminal2: Name of second terminal

        Returns:
            Patristic distance (sum of branch lengths along path)

        Raises:
            ValueError: If no tree has been built or terminals not found
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        return self.tree.distance(terminal1, terminal2)

    def collapse_short_branches(self, min_length: float = 0.0001):
        """
        Collapse branches shorter than a threshold.

        Args:
            min_length: Minimum branch length to keep

        Raises:
            ValueError: If no tree has been built yet
        """
        if self.tree is None:
            raise ValueError("No tree available. Build a tree first.")

        for clade in self.tree.find_clades():
            if clade.branch_length is not None and clade.branch_length < min_length:
                clade.branch_length = 0
