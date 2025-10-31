# Examples

This directory contains example scripts and data for using the phylogenetic structure trees package.

## Example Scripts

### example_usage.py

Demonstrates programmatic usage of the package with several examples:

1. **Basic Usage**: Load structures, calculate RMSD, and build a tree
2. **Loading Individual Files**: Load specific PDB files instead of a directory
3. **Detailed Alignment**: Get detailed information about pairwise alignments
4. **Method Comparison**: Compare UPGMA and Neighbor-Joining tree methods

Run the examples:

```bash
python example_usage.py
```

## Getting Example Data

To test the package, you'll need some protein structure files. Here are some ways to get example data:

### Option 1: Download from PDB

Download structures directly from the [Protein Data Bank](https://www.rcsb.org/):

```bash
# Create a structures directory
mkdir structures

# Download some example structures (requires wget or curl)
wget https://files.rcsb.org/download/1UBQ.pdb -O structures/1UBQ.pdb
wget https://files.rcsb.org/download/1UBI.pdb -O structures/1UBI.pdb
wget https://files.rcsb.org/download/1AAR.pdb -O structures/1AAR.pdb
```

### Option 2: Use Biopython to download structures

```python
from Bio.PDB import PDBList

pdbl = PDBList()
pdbl.retrieve_pdb_file('1UBQ', pdir='structures/', file_format='pdb')
pdbl.retrieve_pdb_file('1UBI', pdir='structures/', file_format='pdb')
pdbl.retrieve_pdb_file('1AAR', pdir='structures/', file_format='pdb')
```

### Option 3: Use your own structures

If you have your own PDB or mmCIF files, simply place them in a directory:

```bash
mkdir structures
cp /path/to/your/structures/*.pdb structures/
```

## Example Workflows

### Workflow 1: Quick tree from directory

```bash
# Build a tree from all structures in a directory
python -m phylo_struct_tree align -i structures/ -o tree.newick

# View the tree
python -m phylo_struct_tree view tree.newick
```

### Workflow 2: Detailed analysis with matrix output

```bash
# Build tree and save distance matrix
python -m phylo_struct_tree align \
    -i structures/ \
    -o tree.newick \
    --matrix-out distances.csv \
    --verbose \
    --print-tree
```

### Workflow 3: Compare methods

```bash
# Build UPGMA tree
python -m phylo_struct_tree align -i structures/ -o upgma_tree.newick -m upgma

# Build NJ tree
python -m phylo_struct_tree align -i structures/ -o nj_tree.newick -m nj

# View both
python -m phylo_struct_tree view upgma_tree.newick
python -m phylo_struct_tree view nj_tree.newick
```

### Workflow 4: Calculate distances only

```bash
# Just calculate and display pairwise RMSD
python -m phylo_struct_tree distance \
    structures/1UBQ.pdb \
    structures/1UBI.pdb \
    structures/1AAR.pdb \
    -o distances.csv
```

## Expected Output

### Tree File (Newick format)

The output tree file will be in Newick format, which looks like:

```
((structure1:0.5,structure2:0.3):0.2,(structure3:0.4,structure4:0.6):0.1);
```

### Distance Matrix (CSV format)

The distance matrix will be saved as a CSV file:

```csv
,1UBQ,1UBI,1AAR
1UBQ,0.0,2.3,5.1
1UBI,2.3,0.0,4.8
1AAR,5.1,4.8,0.0
```

### ASCII Tree Visualization

When using `--print-tree` or the `view` command:

```
     ________________ 1UBQ
____|
    |         _______ 1UBI
    |________|
             |_______ 1AAR
```

## Tips

1. **Structure Selection**: Choose structures that are related (e.g., same protein family) for meaningful trees
2. **File Formats**: Both PDB (.pdb) and mmCIF (.cif) formats are supported
3. **Alignment**: The tool uses CA (alpha carbon) atoms for alignment
4. **Tree Methods**:
   - Use NJ (Neighbor-Joining) for most cases - doesn't assume constant evolution rate
   - Use UPGMA if you expect constant molecular clock

## Troubleshooting

**Problem**: "No CA atoms found"
- Solution: Make sure your structures contain protein (not just DNA/RNA)

**Problem**: "Need at least 2 structures"
- Solution: Provide a directory with at least 2 structure files

**Problem**: "Matrix is not symmetric"
- Solution: This is usually a bug - please report it

**Problem**: Very large RMSD values
- Solution: Structures might be very different or from different proteins
