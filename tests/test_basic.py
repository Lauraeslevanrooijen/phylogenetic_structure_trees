"""
Basic tests for phylo_struct_tree package.

These tests verify that the package imports correctly and basic functionality works.
"""

import unittest
import numpy as np
from pathlib import Path
import tempfile
import os


class TestImports(unittest.TestCase):
    """Test that all modules import correctly."""

    def test_import_structure_alignment(self):
        """Test importing StructureAligner."""
        from phylo_struct_tree import StructureAligner
        aligner = StructureAligner()
        self.assertIsNotNone(aligner)

    def test_import_distance_matrix(self):
        """Test importing DistanceMatrixBuilder."""
        from phylo_struct_tree import DistanceMatrixBuilder
        builder = DistanceMatrixBuilder()
        self.assertIsNotNone(builder)

    def test_import_tree_builder(self):
        """Test importing TreeBuilder."""
        from phylo_struct_tree import TreeBuilder
        builder = TreeBuilder()
        self.assertIsNotNone(builder)


class TestDistanceMatrix(unittest.TestCase):
    """Test distance matrix functionality."""

    def setUp(self):
        """Set up test data."""
        from phylo_struct_tree import DistanceMatrixBuilder
        self.builder = DistanceMatrixBuilder()

    def test_build_simple_matrix(self):
        """Test building a simple distance matrix."""
        # Create simple RMSD dictionary
        rmsd_dict = {
            ('A', 'A'): 0.0,
            ('A', 'B'): 2.0,
            ('A', 'C'): 4.0,
            ('B', 'A'): 2.0,
            ('B', 'B'): 0.0,
            ('B', 'C'): 3.0,
            ('C', 'A'): 4.0,
            ('C', 'B'): 3.0,
            ('C', 'C'): 0.0,
        }
        structure_ids = ['A', 'B', 'C']

        matrix = self.builder.build_matrix(rmsd_dict, structure_ids)

        # Check shape
        self.assertEqual(matrix.shape, (3, 3))

        # Check specific values
        self.assertEqual(matrix[0, 0], 0.0)
        self.assertEqual(matrix[0, 1], 2.0)
        self.assertEqual(matrix[1, 2], 3.0)

        # Check symmetry
        np.testing.assert_array_almost_equal(matrix, matrix.T)

    def test_matrix_validation(self):
        """Test matrix validation."""
        rmsd_dict = {
            ('A', 'A'): 0.0,
            ('A', 'B'): 1.0,
            ('B', 'A'): 1.0,
            ('B', 'B'): 0.0,
        }
        structure_ids = ['A', 'B']

        self.builder.build_matrix(rmsd_dict, structure_ids)
        is_valid, msg = self.builder.is_valid_distance_matrix()

        self.assertTrue(is_valid)

    def test_get_summary_stats(self):
        """Test getting summary statistics."""
        rmsd_dict = {
            ('A', 'A'): 0.0,
            ('A', 'B'): 2.0,
            ('A', 'C'): 4.0,
            ('B', 'A'): 2.0,
            ('B', 'B'): 0.0,
            ('B', 'C'): 3.0,
            ('C', 'A'): 4.0,
            ('C', 'B'): 3.0,
            ('C', 'C'): 0.0,
        }
        structure_ids = ['A', 'B', 'C']

        self.builder.build_matrix(rmsd_dict, structure_ids)
        stats = self.builder.get_summary_stats()

        self.assertIn('mean', stats)
        self.assertIn('min', stats)
        self.assertIn('max', stats)
        self.assertEqual(stats['min'], 2.0)
        self.assertEqual(stats['max'], 4.0)

    def test_normalize_minmax(self):
        """Test min-max normalization."""
        rmsd_dict = {
            ('A', 'A'): 0.0,
            ('A', 'B'): 2.0,
            ('B', 'A'): 2.0,
            ('B', 'B'): 0.0,
        }
        structure_ids = ['A', 'B']

        self.builder.build_matrix(rmsd_dict, structure_ids)
        normalized = self.builder.normalize_matrix(method='minmax')

        # Min-max should give values in [0, 1]
        self.assertTrue(np.all(normalized >= 0))
        self.assertTrue(np.all(normalized <= 1))

    def test_save_csv(self):
        """Test saving matrix to CSV."""
        rmsd_dict = {
            ('A', 'A'): 0.0,
            ('A', 'B'): 1.5,
            ('B', 'A'): 1.5,
            ('B', 'B'): 0.0,
        }
        structure_ids = ['A', 'B']

        self.builder.build_matrix(rmsd_dict, structure_ids)

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            temp_file = f.name

        try:
            self.builder.to_csv(temp_file)
            self.assertTrue(os.path.exists(temp_file))

            # Read and verify
            with open(temp_file, 'r') as f:
                lines = f.readlines()
                self.assertTrue(len(lines) > 0)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


class TestTreeBuilder(unittest.TestCase):
    """Test tree building functionality."""

    def setUp(self):
        """Set up test data."""
        from phylo_struct_tree import TreeBuilder
        self.builder = TreeBuilder()

        # Create a simple distance matrix
        self.matrix = np.array([
            [0.0, 2.0, 4.0],
            [2.0, 0.0, 3.0],
            [4.0, 3.0, 0.0]
        ])
        self.labels = ['A', 'B', 'C']

    def test_build_nj_tree(self):
        """Test building a Neighbor-Joining tree."""
        tree = self.builder.build_nj_tree(self.matrix, self.labels)

        self.assertIsNotNone(tree)
        terminal_names = self.builder.get_terminal_names()
        self.assertEqual(len(terminal_names), 3)
        self.assertIn('A', terminal_names)
        self.assertIn('B', terminal_names)
        self.assertIn('C', terminal_names)

    def test_build_upgma_tree(self):
        """Test building a UPGMA tree."""
        tree = self.builder.build_upgma_tree(self.matrix, self.labels)

        self.assertIsNotNone(tree)
        terminal_names = self.builder.get_terminal_names()
        self.assertEqual(len(terminal_names), 3)

    def test_get_tree_info(self):
        """Test getting tree information."""
        self.builder.build_nj_tree(self.matrix, self.labels)
        info = self.builder.get_tree_info()

        self.assertIn('method', info)
        self.assertIn('num_terminals', info)
        self.assertEqual(info['num_terminals'], 3)
        self.assertEqual(info['method'], 'nj')

    def test_save_tree(self):
        """Test saving tree to file."""
        self.builder.build_nj_tree(self.matrix, self.labels)

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.newick') as f:
            temp_file = f.name

        try:
            self.builder.save_tree(temp_file)
            self.assertTrue(os.path.exists(temp_file))

            # Verify file has content
            with open(temp_file, 'r') as f:
                content = f.read()
                self.assertTrue(len(content) > 0)
                self.assertIn('A', content)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_get_newick_string(self):
        """Test getting Newick string."""
        self.builder.build_nj_tree(self.matrix, self.labels)
        newick = self.builder.get_newick_string()

        self.assertTrue(isinstance(newick, str))
        self.assertTrue(len(newick) > 0)
        self.assertIn('A', newick)


if __name__ == '__main__':
    unittest.main()
