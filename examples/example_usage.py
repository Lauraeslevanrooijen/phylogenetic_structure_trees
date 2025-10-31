#!/usr/bin/env python3
"""
Example usage of the phylo_struct_tree package.

This script demonstrates how to use the package programmatically
to build phylogenetic trees from protein structures.
"""

from phylo_struct_tree import StructureAligner, DistanceMatrixBuilder, TreeBuilder


def example_basic_usage():
    """
    Basic example: Load structures, align them, and build a tree.
    """
    print("=" * 60)
    print("Example 1: Basic Usage")
    print("=" * 60)

    # Initialize components
    aligner = StructureAligner(quiet=False)
    matrix_builder = DistanceMatrixBuilder()
    tree_builder = TreeBuilder()

    # Load structures from a directory
    # Replace with your actual structure directory
    structure_dir = "structures/"

    print(f"\n1. Loading structures from {structure_dir}...")
    try:
        structures = aligner.load_structures_from_directory(structure_dir)
        structure_ids = list(structures.keys())
        print(f"   Loaded {len(structure_ids)} structures: {', '.join(structure_ids)}")
    except Exception as e:
        print(f"   Error: {e}")
        print("   Please provide a directory with PDB/mmCIF files")
        return

    # Calculate pairwise RMSD
    print("\n2. Calculating pairwise RMSD values...")
    rmsd_dict = aligner.calculate_pairwise_rmsd(structure_ids)
    print(f"   Calculated {len(rmsd_dict)} pairwise distances")

    # Build distance matrix
    print("\n3. Building distance matrix...")
    distance_matrix = matrix_builder.build_matrix(rmsd_dict, structure_ids)

    # Print matrix statistics
    stats = matrix_builder.get_summary_stats()
    print(f"   Matrix size: {stats['size']} x {stats['size']}")
    print(f"   Mean RMSD: {stats['mean']:.4f} Å")
    print(f"   Range: {stats['min']:.4f} - {stats['max']:.4f} Å")

    # Build tree using Neighbor-Joining
    print("\n4. Building phylogenetic tree (Neighbor-Joining)...")
    tree = tree_builder.build_nj_tree(distance_matrix, structure_ids)

    # Display tree information
    info = tree_builder.get_tree_info()
    print(f"   Terminal nodes: {info['num_terminals']}")
    print(f"   Internal nodes: {info['num_internal_nodes']}")
    print(f"   Total branch length: {info['total_branch_length']:.4f}")

    # Print ASCII tree
    print("\n5. Tree visualization:")
    tree_builder.print_tree()

    # Save tree
    output_file = "example_tree.newick"
    tree_builder.save_tree(output_file)
    print(f"\n6. Tree saved to: {output_file}")

    # Save distance matrix
    matrix_file = "example_distances.csv"
    matrix_builder.to_csv(matrix_file)
    print(f"   Distance matrix saved to: {matrix_file}")


def example_load_individual_files():
    """
    Example: Load individual structure files.
    """
    print("\n" + "=" * 60)
    print("Example 2: Loading Individual Files")
    print("=" * 60)

    aligner = StructureAligner()

    # Load individual PDB files
    structure_files = [
        "structure1.pdb",
        "structure2.pdb",
        "structure3.pdb"
    ]

    print("\nLoading individual structure files...")
    structure_ids = []

    for filepath in structure_files:
        try:
            structure = aligner.load_structure(filepath)
            structure_ids.append(structure.get_id())
            print(f"  ✓ Loaded: {filepath}")
        except FileNotFoundError:
            print(f"  ✗ File not found: {filepath}")

    if len(structure_ids) >= 2:
        print(f"\nSuccessfully loaded {len(structure_ids)} structures")
        # Continue with alignment and tree building...
    else:
        print("\nNeed at least 2 structures to build a tree")


def example_alignment_details():
    """
    Example: Get detailed alignment information.
    """
    print("\n" + "=" * 60)
    print("Example 3: Detailed Alignment Information")
    print("=" * 60)

    aligner = StructureAligner()

    # Load two structures
    try:
        struct1 = aligner.load_structure("structure1.pdb", "protein1")
        struct2 = aligner.load_structure("structure2.pdb", "protein2")

        print("\nGetting alignment details...")
        alignment_info = aligner.get_alignment_info(struct1, struct2)

        print(f"\nAlignment Results:")
        print(f"  RMSD: {alignment_info['rmsd']:.4f} Å")
        print(f"  Protein 1 residues: {alignment_info['num_residues_1']}")
        print(f"  Protein 2 residues: {alignment_info['num_residues_2']}")
        print(f"  Aligned residues: {alignment_info['aligned_residues']}")
        print(f"  Rotation matrix shape: {alignment_info['rotation_matrix'].shape}")

    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("Please provide valid PDB files")


def example_upgma_vs_nj():
    """
    Example: Compare UPGMA and Neighbor-Joining trees.
    """
    print("\n" + "=" * 60)
    print("Example 4: Comparing UPGMA and NJ Methods")
    print("=" * 60)

    aligner = StructureAligner()
    matrix_builder = DistanceMatrixBuilder()

    try:
        # Load structures
        structures = aligner.load_structures_from_directory("structures/")
        structure_ids = list(structures.keys())

        # Calculate distances
        rmsd_dict = aligner.calculate_pairwise_rmsd(structure_ids)
        distance_matrix = matrix_builder.build_matrix(rmsd_dict, structure_ids)

        # Build UPGMA tree
        print("\n1. Building UPGMA tree...")
        upgma_builder = TreeBuilder()
        upgma_tree = upgma_builder.build_upgma_tree(distance_matrix, structure_ids)
        upgma_info = upgma_builder.get_tree_info()
        print(f"   Total branch length: {upgma_info['total_branch_length']:.4f}")

        # Build NJ tree
        print("\n2. Building Neighbor-Joining tree...")
        nj_builder = TreeBuilder()
        nj_tree = nj_builder.build_nj_tree(distance_matrix, structure_ids)
        nj_info = nj_builder.get_tree_info()
        print(f"   Total branch length: {nj_info['total_branch_length']:.4f}")

        # Save both trees
        upgma_builder.save_tree("upgma_tree.newick")
        nj_builder.save_tree("nj_tree.newick")
        print("\n✓ Trees saved: upgma_tree.newick, nj_tree.newick")

    except Exception as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    print("\nPhylogenetic Structure Trees - Example Usage\n")

    # Run the basic example
    example_basic_usage()

    # Uncomment to run other examples:
    # example_load_individual_files()
    # example_alignment_details()
    # example_upgma_vs_nj()

    print("\n" + "=" * 60)
    print("Examples completed!")
    print("=" * 60)
