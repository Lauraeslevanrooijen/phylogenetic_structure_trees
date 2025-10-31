#!/usr/bin/env python3
"""
Command-line interface for phylogenetic structure tree construction.
"""

import sys
import os
from pathlib import Path
from typing import List

import click

from .structure_alignment import StructureAligner
from .distance_matrix import DistanceMatrixBuilder
from .tree_builder import TreeBuilder


@click.group()
@click.version_option(version='0.1.0')
def main():
    """
    Phylogenetic Structure Trees - Build phylogenetic trees from protein structures.

    This tool performs structural alignment of proteins and constructs
    phylogenetic trees based on structural similarity (RMSD).
    """
    pass


@main.command()
@click.option('--input', '-i', 'input_path', required=True,
              help='Input PDB/mmCIF file(s) or directory containing structures')
@click.option('--output', '-o', 'output_file', required=True,
              help='Output tree file (Newick format)')
@click.option('--method', '-m', type=click.Choice(['upgma', 'nj']), default='nj',
              help='Tree building method (default: nj)')
@click.option('--matrix-out', type=click.Path(),
              help='Optional: Save distance matrix to CSV file')
@click.option('--normalize', type=click.Choice(['none', 'minmax', 'zscore']), default='none',
              help='Normalize distance matrix (default: none)')
@click.option('--verbose', '-v', is_flag=True,
              help='Print verbose output')
@click.option('--print-tree', is_flag=True,
              help='Print ASCII tree to console')
def align(input_path, output_file, method, matrix_out, normalize, verbose, print_tree):
    """
    Align protein structures and build a phylogenetic tree.

    INPUT can be:
    - A directory containing PDB/mmCIF files
    - Multiple PDB/mmCIF files (space-separated)

    Examples:

    \b
    # Build tree from all structures in a directory
    phylo-struct-tree align -i structures/ -o tree.newick

    \b
    # Use UPGMA method and save distance matrix
    phylo-struct-tree align -i structures/ -o tree.newick -m upgma --matrix-out distances.csv

    \b
    # Print the tree to console
    phylo-struct-tree align -i structures/ -o tree.newick --print-tree
    """
    try:
        # Initialize components
        aligner = StructureAligner(quiet=not verbose)
        matrix_builder = DistanceMatrixBuilder()
        tree_builder = TreeBuilder()

        # Load structures
        if verbose:
            click.echo(f"Loading structures from: {input_path}")

        input_path_obj = Path(input_path)

        if input_path_obj.is_dir():
            structures = aligner.load_structures_from_directory(input_path)
            structure_ids = list(structures.keys())
        elif input_path_obj.is_file():
            structure = aligner.load_structure(str(input_path_obj))
            structure_ids = [input_path_obj.stem]
            click.echo("Error: Need at least 2 structures to build a tree", err=True)
            sys.exit(1)
        else:
            click.echo(f"Error: Input path not found: {input_path}", err=True)
            sys.exit(1)

        if len(structure_ids) < 2:
            click.echo("Error: Need at least 2 structures to build a tree", err=True)
            sys.exit(1)

        if verbose:
            click.echo(f"Loaded {len(structure_ids)} structures: {', '.join(structure_ids)}")

        # Calculate pairwise RMSD
        if verbose:
            click.echo("Calculating pairwise structural alignments...")

        rmsd_dict = aligner.calculate_pairwise_rmsd(structure_ids)

        # Build distance matrix
        if verbose:
            click.echo("Building distance matrix...")

        distance_matrix = matrix_builder.build_matrix(rmsd_dict, structure_ids)

        # Normalize if requested
        if normalize != 'none':
            if verbose:
                click.echo(f"Normalizing matrix using {normalize} method...")
            distance_matrix = matrix_builder.normalize_matrix(method=normalize)

        # Validate matrix
        is_valid, error_msg = matrix_builder.is_valid_distance_matrix()
        if not is_valid:
            click.echo(f"Warning: Distance matrix validation failed: {error_msg}", err=True)

        # Print matrix statistics
        if verbose:
            stats = matrix_builder.get_summary_stats()
            click.echo("\nDistance Matrix Statistics:")
            click.echo(f"  Size: {stats['size']} x {stats['size']}")
            click.echo(f"  Mean RMSD: {stats['mean']:.4f}")
            click.echo(f"  Median RMSD: {stats['median']:.4f}")
            click.echo(f"  Min RMSD: {stats['min']:.4f}")
            click.echo(f"  Max RMSD: {stats['max']:.4f}")
            click.echo(f"  Std Dev: {stats['std']:.4f}")

        # Save matrix if requested
        if matrix_out:
            matrix_builder.to_csv(matrix_out)
            if verbose:
                click.echo(f"Distance matrix saved to: {matrix_out}")

        # Build tree
        if verbose:
            click.echo(f"\nBuilding phylogenetic tree using {method.upper()}...")

        tree = tree_builder.build_tree(distance_matrix, structure_ids, method=method)

        # Get tree info
        if verbose:
            info = tree_builder.get_tree_info()
            click.echo("\nTree Statistics:")
            click.echo(f"  Method: {info['method'].upper()}")
            click.echo(f"  Terminal nodes: {info['num_terminals']}")
            click.echo(f"  Internal nodes: {info['num_internal_nodes']}")
            click.echo(f"  Total branch length: {info['total_branch_length']:.4f}")
            click.echo(f"  Max depth: {info['max_depth']:.4f}")
            click.echo(f"  Bifurcating: {info['is_bifurcating']}")

        # Print tree if requested
        if print_tree:
            click.echo("\nPhylogenetic Tree:")
            tree_builder.print_tree()

        # Save tree
        tree_builder.save_tree(output_file, format='newick')

        click.echo(f"\n✓ Tree saved to: {output_file}")

    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Unexpected error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.argument('structure_files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('--output', '-o', 'output_file',
              help='Output CSV file for distance matrix')
@click.option('--verbose', '-v', is_flag=True,
              help='Print verbose output')
def distance(structure_files, output_file, verbose):
    """
    Calculate pairwise RMSD distances between structures.

    STRUCTURE_FILES: Two or more PDB/mmCIF files

    Example:

    \b
    phylo-struct-tree distance struct1.pdb struct2.pdb struct3.pdb -o distances.csv
    """
    try:
        if len(structure_files) < 2:
            click.echo("Error: Need at least 2 structure files", err=True)
            sys.exit(1)

        aligner = StructureAligner(quiet=not verbose)
        matrix_builder = DistanceMatrixBuilder()

        # Load structures
        structure_ids = []
        for filepath in structure_files:
            path = Path(filepath)
            structure_id = path.stem
            aligner.load_structure(str(path), structure_id)
            structure_ids.append(structure_id)

        if verbose:
            click.echo(f"Loaded {len(structure_ids)} structures")

        # Calculate distances
        if verbose:
            click.echo("Calculating pairwise RMSD...")

        rmsd_dict = aligner.calculate_pairwise_rmsd(structure_ids)

        # Build matrix
        distance_matrix = matrix_builder.build_matrix(rmsd_dict, structure_ids)

        # Print matrix
        click.echo("\nDistance Matrix (RMSD):")
        click.echo("       " + "  ".join(f"{sid:>8}" for sid in structure_ids))

        for i, sid in enumerate(structure_ids):
            row_str = f"{sid:>6} "
            for j in range(len(structure_ids)):
                row_str += f"{distance_matrix[i, j]:>8.4f}  "
            click.echo(row_str)

        # Print statistics
        stats = matrix_builder.get_summary_stats()
        click.echo(f"\nStatistics:")
        click.echo(f"  Mean: {stats['mean']:.4f}")
        click.echo(f"  Min:  {stats['min']:.4f}")
        click.echo(f"  Max:  {stats['max']:.4f}")

        # Save if output file specified
        if output_file:
            matrix_builder.to_csv(output_file)
            click.echo(f"\n✓ Distance matrix saved to: {output_file}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.argument('tree_file', type=click.Path(exists=True))
def view(tree_file):
    """
    View a phylogenetic tree in ASCII format.

    TREE_FILE: Path to Newick format tree file

    Example:

    \b
    phylo-struct-tree view tree.newick
    """
    try:
        from Bio import Phylo

        tree = Phylo.read(tree_file, 'newick')

        click.echo(f"\nTree from: {tree_file}\n")
        Phylo.draw_ascii(tree)

        # Print info
        terminals = tree.get_terminals()
        click.echo(f"\nTerminal nodes: {len(terminals)}")
        click.echo(f"Total branch length: {tree.total_branch_length():.4f}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
