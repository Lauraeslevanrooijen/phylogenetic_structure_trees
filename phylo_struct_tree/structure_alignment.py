"""
Structure alignment module for protein structure superposition and RMSD calculation.
"""

import os
import warnings
from typing import List, Tuple, Dict, Optional
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Superimposer, Selection
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom


class StructureAligner:
    """
    Handles loading and aligning protein structures.

    Uses Biopython's PDB module to:
    - Load structures from PDB/mmCIF files
    - Extract CA atoms for alignment
    - Perform structural superposition
    - Calculate RMSD between structures
    """

    def __init__(self, quiet: bool = True):
        """
        Initialize the StructureAligner.

        Args:
            quiet: If True, suppress PDB parser warnings
        """
        self.pdb_parser = PDBParser(QUIET=quiet)
        self.cif_parser = MMCIFParser(QUIET=quiet)
        self.superimposer = Superimposer()
        self.structures: Dict[str, Structure] = {}

        if quiet:
            warnings.filterwarnings('ignore')

    def load_structure(self, filepath: str, structure_id: Optional[str] = None) -> Structure:
        """
        Load a structure from a PDB or mmCIF file.

        Args:
            filepath: Path to the structure file
            structure_id: Optional ID for the structure. If None, uses filename.

        Returns:
            Loaded Bio.PDB.Structure object

        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file format is not supported
        """
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"Structure file not found: {filepath}")

        if structure_id is None:
            structure_id = filepath.stem

        # Determine file format and use appropriate parser
        if filepath.suffix.lower() in ['.pdb', '.ent']:
            structure = self.pdb_parser.get_structure(structure_id, str(filepath))
        elif filepath.suffix.lower() in ['.cif', '.mmcif']:
            structure = self.cif_parser.get_structure(structure_id, str(filepath))
        else:
            raise ValueError(f"Unsupported file format: {filepath.suffix}")

        self.structures[structure_id] = structure
        return structure

    def load_structures_from_directory(self, directory: str) -> Dict[str, Structure]:
        """
        Load all PDB/mmCIF structures from a directory.

        Args:
            directory: Path to directory containing structure files

        Returns:
            Dictionary mapping structure IDs to Structure objects
        """
        directory = Path(directory)

        if not directory.is_dir():
            raise NotADirectoryError(f"Not a directory: {directory}")

        # Find all structure files
        patterns = ['*.pdb', '*.ent', '*.cif', '*.mmcif']
        structure_files = []
        for pattern in patterns:
            structure_files.extend(directory.glob(pattern))

        if not structure_files:
            raise ValueError(f"No structure files found in {directory}")

        # Load all structures
        for filepath in structure_files:
            self.load_structure(str(filepath))

        return self.structures

    def get_ca_atoms(self, structure: Structure) -> List[Atom]:
        """
        Extract all CA (alpha carbon) atoms from a structure.

        Args:
            structure: Bio.PDB.Structure object

        Returns:
            List of CA atoms
        """
        ca_atoms = []

        # Get all CA atoms from the first model
        model = structure[0]
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])

        return ca_atoms

    def align_structures(self, structure1: Structure, structure2: Structure,
                        match_by_id: bool = False) -> Tuple[float, np.ndarray]:
        """
        Align two structures and calculate RMSD.

        Uses CA atoms for alignment. If structures have different numbers of residues,
        aligns based on the minimum common length.

        Args:
            structure1: First structure (will be moved to align with second)
            structure2: Second structure (reference)
            match_by_id: If True, match residues by their ID, otherwise by position

        Returns:
            Tuple of (RMSD value, rotation matrix)
        """
        ca_atoms1 = self.get_ca_atoms(structure1)
        ca_atoms2 = self.get_ca_atoms(structure2)

        if not ca_atoms1 or not ca_atoms2:
            raise ValueError("One or both structures have no CA atoms")

        # Match atoms
        if match_by_id:
            # Match by residue ID
            atoms1, atoms2 = self._match_atoms_by_id(ca_atoms1, ca_atoms2)
        else:
            # Match by position (align first N residues)
            min_length = min(len(ca_atoms1), len(ca_atoms2))
            atoms1 = ca_atoms1[:min_length]
            atoms2 = ca_atoms2[:min_length]

        if len(atoms1) < 3:
            raise ValueError("Need at least 3 matching atoms for alignment")

        # Perform superposition
        self.superimposer.set_atoms(atoms2, atoms1)

        return self.superimposer.rms, self.superimposer.rotran[0]

    def _match_atoms_by_id(self, atoms1: List[Atom], atoms2: List[Atom]) -> Tuple[List[Atom], List[Atom]]:
        """
        Match atoms from two structures by residue ID.

        Args:
            atoms1: CA atoms from first structure
            atoms2: CA atoms from second structure

        Returns:
            Tuple of (matched_atoms1, matched_atoms2)
        """
        # Create dictionaries mapping residue ID to atom
        atom_dict1 = {atom.get_parent().get_id(): atom for atom in atoms1}
        atom_dict2 = {atom.get_parent().get_id(): atom for atom in atoms2}

        # Find common residue IDs
        common_ids = set(atom_dict1.keys()) & set(atom_dict2.keys())
        common_ids = sorted(common_ids)

        matched_atoms1 = [atom_dict1[rid] for rid in common_ids]
        matched_atoms2 = [atom_dict2[rid] for rid in common_ids]

        return matched_atoms1, matched_atoms2

    def calculate_pairwise_rmsd(self, structure_ids: Optional[List[str]] = None) -> Dict[Tuple[str, str], float]:
        """
        Calculate pairwise RMSD for all loaded structures.

        Args:
            structure_ids: Optional list of structure IDs to compare.
                         If None, uses all loaded structures.

        Returns:
            Dictionary mapping (structure_id1, structure_id2) to RMSD value
        """
        if structure_ids is None:
            structure_ids = list(self.structures.keys())

        if len(structure_ids) < 2:
            raise ValueError("Need at least 2 structures for pairwise comparison")

        rmsd_dict = {}

        # Calculate all pairwise RMSDs
        for i, id1 in enumerate(structure_ids):
            for id2 in structure_ids[i+1:]:
                try:
                    rmsd, _ = self.align_structures(
                        self.structures[id1],
                        self.structures[id2]
                    )
                    rmsd_dict[(id1, id2)] = rmsd
                    rmsd_dict[(id2, id1)] = rmsd  # Symmetric
                except Exception as e:
                    print(f"Warning: Could not align {id1} and {id2}: {e}")
                    # Use a large distance for failed alignments
                    rmsd_dict[(id1, id2)] = float('inf')
                    rmsd_dict[(id2, id1)] = float('inf')

        # Add self-comparisons (RMSD = 0)
        for id in structure_ids:
            rmsd_dict[(id, id)] = 0.0

        return rmsd_dict

    def get_alignment_info(self, structure1: Structure, structure2: Structure) -> Dict:
        """
        Get detailed alignment information for two structures.

        Args:
            structure1: First structure
            structure2: Second structure

        Returns:
            Dictionary with alignment details including RMSD, num_residues, etc.
        """
        ca_atoms1 = self.get_ca_atoms(structure1)
        ca_atoms2 = self.get_ca_atoms(structure2)

        min_length = min(len(ca_atoms1), len(ca_atoms2))

        rmsd, rotation = self.align_structures(structure1, structure2)

        return {
            'rmsd': rmsd,
            'num_residues_1': len(ca_atoms1),
            'num_residues_2': len(ca_atoms2),
            'aligned_residues': min_length,
            'rotation_matrix': rotation
        }
