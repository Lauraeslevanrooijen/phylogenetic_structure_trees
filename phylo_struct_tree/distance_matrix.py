"""
Distance matrix module for converting RMSD values to distance matrices.
"""

import numpy as np
from typing import Dict, Tuple, List
import csv


class DistanceMatrixBuilder:
    """
    Builds distance matrices from pairwise RMSD values.

    The distance matrix can be used as input for phylogenetic tree
    construction algorithms like UPGMA and Neighbor-Joining.
    """

    def __init__(self):
        """Initialize the DistanceMatrixBuilder."""
        self.matrix = None
        self.labels = None

    def build_matrix(self, rmsd_dict: Dict[Tuple[str, str], float],
                     structure_ids: List[str]) -> np.ndarray:
        """
        Build a distance matrix from pairwise RMSD values.

        Args:
            rmsd_dict: Dictionary mapping (id1, id2) to RMSD value
            structure_ids: List of structure IDs (defines matrix order)

        Returns:
            Distance matrix as numpy array
        """
        n = len(structure_ids)
        matrix = np.zeros((n, n))

        for i, id1 in enumerate(structure_ids):
            for j, id2 in enumerate(structure_ids):
                if i == j:
                    matrix[i, j] = 0.0
                else:
                    key = (id1, id2)
                    if key in rmsd_dict:
                        matrix[i, j] = rmsd_dict[key]
                    else:
                        # Try reverse key
                        key = (id2, id1)
                        if key in rmsd_dict:
                            matrix[i, j] = rmsd_dict[key]
                        else:
                            raise ValueError(f"Missing RMSD value for pair ({id1}, {id2})")

        self.matrix = matrix
        self.labels = structure_ids

        return matrix

    def normalize_matrix(self, method: str = 'none') -> np.ndarray:
        """
        Normalize the distance matrix.

        Args:
            method: Normalization method:
                   'none' - no normalization
                   'minmax' - scale to [0, 1]
                   'zscore' - standardize to mean=0, std=1

        Returns:
            Normalized distance matrix
        """
        if self.matrix is None:
            raise ValueError("No matrix to normalize. Call build_matrix first.")

        if method == 'none':
            return self.matrix.copy()

        elif method == 'minmax':
            # Get non-diagonal values for scaling
            non_diag = self.matrix[~np.eye(self.matrix.shape[0], dtype=bool)]
            min_val = non_diag.min()
            max_val = non_diag.max()

            if max_val - min_val == 0:
                return np.zeros_like(self.matrix)

            normalized = (self.matrix - min_val) / (max_val - min_val)
            # Ensure diagonal is still 0
            np.fill_diagonal(normalized, 0)
            return normalized

        elif method == 'zscore':
            # Get non-diagonal values for scaling
            non_diag = self.matrix[~np.eye(self.matrix.shape[0], dtype=bool)]
            mean_val = non_diag.mean()
            std_val = non_diag.std()

            if std_val == 0:
                return np.zeros_like(self.matrix)

            normalized = (self.matrix - mean_val) / std_val
            # Ensure diagonal is still 0
            np.fill_diagonal(normalized, 0)
            return normalized

        else:
            raise ValueError(f"Unknown normalization method: {method}")

    def to_phylip(self, filepath: str):
        """
        Save distance matrix in PHYLIP format.

        Args:
            filepath: Output file path
        """
        if self.matrix is None or self.labels is None:
            raise ValueError("No matrix to save. Call build_matrix first.")

        with open(filepath, 'w') as f:
            # Write header (number of taxa)
            f.write(f"{len(self.labels)}\n")

            # Write matrix
            for i, label in enumerate(self.labels):
                # PHYLIP format: label (10 chars) followed by distances
                f.write(f"{label[:10]:<10}")
                for j in range(len(self.labels)):
                    f.write(f" {self.matrix[i, j]:.6f}")
                f.write("\n")

    def to_csv(self, filepath: str):
        """
        Save distance matrix in CSV format.

        Args:
            filepath: Output file path
        """
        if self.matrix is None or self.labels is None:
            raise ValueError("No matrix to save. Call build_matrix first.")

        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)

            # Write header
            writer.writerow([''] + self.labels)

            # Write matrix rows
            for i, label in enumerate(self.labels):
                row = [label] + list(self.matrix[i, :])
                writer.writerow(row)

    def get_matrix(self) -> np.ndarray:
        """
        Get the current distance matrix.

        Returns:
            Distance matrix as numpy array

        Raises:
            ValueError: If no matrix has been built yet
        """
        if self.matrix is None:
            raise ValueError("No matrix available. Call build_matrix first.")
        return self.matrix.copy()

    def get_labels(self) -> List[str]:
        """
        Get the labels (structure IDs) for the matrix.

        Returns:
            List of structure IDs

        Raises:
            ValueError: If no matrix has been built yet
        """
        if self.labels is None:
            raise ValueError("No labels available. Call build_matrix first.")
        return self.labels.copy()

    def is_valid_distance_matrix(self) -> Tuple[bool, str]:
        """
        Check if the current matrix is a valid distance matrix.

        A valid distance matrix should be:
        - Square
        - Symmetric
        - Have zeros on diagonal
        - Have non-negative values

        Returns:
            Tuple of (is_valid, error_message)
        """
        if self.matrix is None:
            return False, "No matrix has been built"

        # Check if square
        if self.matrix.shape[0] != self.matrix.shape[1]:
            return False, "Matrix is not square"

        # Check if symmetric
        if not np.allclose(self.matrix, self.matrix.T):
            return False, "Matrix is not symmetric"

        # Check diagonal is zero
        if not np.allclose(np.diag(self.matrix), 0):
            return False, "Diagonal is not zero"

        # Check non-negative
        if np.any(self.matrix < 0):
            return False, "Matrix contains negative values"

        # Check for infinite values
        if np.any(np.isinf(self.matrix)):
            return False, "Matrix contains infinite values"

        return True, "Matrix is valid"

    def get_summary_stats(self) -> Dict:
        """
        Get summary statistics for the distance matrix.

        Returns:
            Dictionary with statistics (mean, median, min, max, std)
        """
        if self.matrix is None:
            raise ValueError("No matrix available. Call build_matrix first.")

        # Get non-diagonal values
        non_diag = self.matrix[~np.eye(self.matrix.shape[0], dtype=bool)]

        return {
            'mean': float(np.mean(non_diag)),
            'median': float(np.median(non_diag)),
            'min': float(np.min(non_diag)),
            'max': float(np.max(non_diag)),
            'std': float(np.std(non_diag)),
            'size': len(self.labels)
        }
