# Phylogenetic Structure Trees

A command-line tool for creating phylogenetic trees based on protein structure alignment.

## Features

- **Structure Alignment**: Align protein structures using superposition and RMSD calculation
- **Distance Matrix**: Generate distance matrices from structural alignments
- **Phylogenetic Trees**: Build phylogenetic trees using distance-based methods (UPGMA, Neighbor-Joining)
- **CLI Interface**: Easy-to-use command-line interface

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic usage - align structures and build tree:

```bash
python -m phylo_struct_tree align --input structures/ --output tree.newick
```

### Align structures from PDB files:

```bash
python -m phylo_struct_tree align --input structure1.pdb structure2.pdb structure3.pdb --output tree.newick
```

### Specify tree building method:

```bash
python -m phylo_struct_tree align --input structures/ --method nj --output tree.newick
```

Available methods:
- `upgma`: UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- `nj`: Neighbor-Joining (default)

## Input Format

The tool accepts:
- PDB files (.pdb)
- mmCIF files (.cif)
- Directories containing multiple structure files

## Output Format

- Newick format tree files
- Optional: distance matrix (CSV)
- Optional: alignment visualization

## Requirements

- Python 3.8+
- Biopython
- NumPy
- SciPy

## How It Works

1. **Structure Loading**: Reads protein structures from PDB/mmCIF files
2. **Pairwise Alignment**: Performs structural superposition for all pairs using:
   - Alpha carbon (CA) alignment
   - SVD-based superposition
   - RMSD calculation
3. **Distance Matrix**: Converts RMSD values to evolutionary distances
4. **Tree Construction**: Builds phylogenetic tree using clustering algorithms

## License

MIT License
