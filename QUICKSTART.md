# Quick Start Guide

Get up and running with phylogenetic structure trees in 5 minutes!

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd phylogenetic_structure_trees
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. (Optional) Install the package:
```bash
pip install -e .
```

## Your First Tree

### Step 1: Get some protein structures

Create a directory and download a few related structures:

```bash
mkdir structures
cd structures

# Download three ubiquitin structures (small, related proteins)
wget https://files.rcsb.org/download/1UBQ.pdb
wget https://files.rcsb.org/download/1UBI.pdb
wget https://files.rcsb.org/download/3ONS.pdb

cd ..
```

### Step 2: Build a tree

```bash
python -m phylo_struct_tree align -i structures/ -o my_first_tree.newick --verbose --print-tree
```

You should see output like:

```
Loading structures from: structures/
Loaded 3 structures: 1UBQ, 1UBI, 3ONS

Calculating pairwise structural alignments...

Building distance matrix...

Distance Matrix Statistics:
  Size: 3 x 3
  Mean RMSD: 2.1234
  ...

Building phylogenetic tree using NJ...

Phylogenetic Tree:
     _________ 1UBQ
____|
    |    _____ 1UBI
    |___|
        |_____ 3ONS

✓ Tree saved to: my_first_tree.newick
```

### Step 3: View your tree

```bash
python -m phylo_struct_tree view my_first_tree.newick
```

## Next Steps

### Save the distance matrix

```bash
python -m phylo_struct_tree align -i structures/ -o tree.newick --matrix-out distances.csv
```

Open `distances.csv` in Excel or any spreadsheet program to see the pairwise RMSD values.

### Try different tree methods

**UPGMA** (assumes constant evolution rate):
```bash
python -m phylo_struct_tree align -i structures/ -o upgma_tree.newick -m upgma
```

**Neighbor-Joining** (default, doesn't assume constant rate):
```bash
python -m phylo_struct_tree align -i structures/ -o nj_tree.newick -m nj
```

### Calculate distances only

```bash
python -m phylo_struct_tree distance structures/*.pdb -o distances.csv
```

### Use in Python code

```python
from phylo_struct_tree import StructureAligner, DistanceMatrixBuilder, TreeBuilder

# Load structures
aligner = StructureAligner()
structures = aligner.load_structures_from_directory("structures/")
structure_ids = list(structures.keys())

# Calculate RMSD
rmsd_dict = aligner.calculate_pairwise_rmsd(structure_ids)

# Build matrix
matrix_builder = DistanceMatrixBuilder()
distance_matrix = matrix_builder.build_matrix(rmsd_dict, structure_ids)

# Build tree
tree_builder = TreeBuilder()
tree = tree_builder.build_nj_tree(distance_matrix, structure_ids)

# Save
tree_builder.save_tree("output.newick")
tree_builder.print_tree()
```

## Common Use Cases

### 1. Analyze protein family evolution

Collect structures from a protein family and build a tree to see evolutionary relationships.

### 2. Compare protein conformations

Use structures of the same protein in different conformational states to see how different they are.

### 3. Assess structural similarity

Quickly calculate RMSD between multiple structures to identify the most similar pairs.

### 4. Quality control

Compare experimental structures to predicted models (e.g., from AlphaFold) to assess quality.

## Tips

- **CA atoms**: The tool uses alpha carbon atoms for alignment, which works for proteins
- **File formats**: Supports both PDB (.pdb) and mmCIF (.cif) formats
- **RMSD units**: RMSD values are in Ångströms (Å)
- **Minimum structures**: You need at least 2 structures to build a tree
- **Related structures**: Trees are most meaningful when structures are related (same family, homologs, etc.)

## Troubleshooting

**Issue**: Package not found when running `python -m phylo_struct_tree`

**Solution**: Make sure you're in the repository directory, or install with `pip install -e .`

---

**Issue**: No structures found in directory

**Solution**: Make sure your PDB/CIF files have the correct extensions (.pdb, .ent, .cif, .mmcif)

---

**Issue**: Very high RMSD values (>10 Å)

**Solution**: Your structures might be very different or unrelated. Make sure they're from similar proteins.

## Getting Help

- Check the full README.md for detailed documentation
- Look at examples in the `examples/` directory
- Run with `--help` flag: `python -m phylo_struct_tree align --help`

## What's Next?

Check out the `examples/` directory for more advanced usage patterns and the full README.md for complete documentation.

Happy tree building! 🌳
