#!/usr/bin/env python
from pymol import cmd
import argparse
import os
import numpy as np
from scipy.spatial import cKDTree

# Amino acid mapping
aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'
}

def analyze_trp_interactions(input_file, output_file, distance_cutoff=5.0, k=100):
    """
    Analyze protein interactions with Trp ligand (UNK residue) using KDTree method
    
    Parameters:
    input_file: Input PDB file path
    output_file: Output file path
    distance_cutoff: Distance threshold in Angstroms
    k: Maximum number of neighbors to consider
    """
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    
    try:
        # Load PDB file
        cmd.load(input_file, base_name)
        
        # Remove solvent
        cmd.remove("solvent")
        
        # Get protein coordinates and residue information
        protein = cmd.get_model("polymer")
        protein_coords = np.array([atom.coord for atom in protein.atom])
        residue_info = [
            f"{atom.chain}:{aa_codes.get(atom.resn, 'X')}{atom.resi}" 
            for atom in protein.atom
        ]
        
        # Get ligand coordinates
        ligand = cmd.get_model("resn UNK")
        ligand_coords = np.array([atom.coord for atom in ligand.atom])
        
        # Check if ligand exists
        if len(ligand_coords) == 0:
            print(f"  Warning: No Trp ligand (UNK) found in {base_name}")
            with open(output_file, "a") as f:
                f.write(f"{base_name}\tNo Trp ligand found\n")
            return
        
        # Build KDTree for protein atoms
        tree = cKDTree(protein_coords)
        
        # Find protein atoms within distance cutoff of ligand atoms
        mapping = set()
        for lig_atom in ligand_coords:
            # Query the tree for nearby protein atoms
            distances, indices = tree.query(
                lig_atom, 
                k=k, 
                distance_upper_bound=distance_cutoff
            )
            
            # Filter out infinite distances
            valid_indices = indices[distances != np.inf]
            mapping.update(valid_indices)
        
        # Get unique residue information
        interacting_residues = set()
        for idx in mapping:
            residue_id = residue_info[idx]
            interacting_residues.add(residue_id)
        
        # Sort residues by chain and residue number
        sorted_residues = sorted(
            interacting_residues, 
            key=lambda x: (x.split(':')[0], int(''.join(filter(str.isdigit, x))))
        )
        
        # Format output
        formatted_residues = ", ".join([f"'{res}'" for res in sorted_residues])
        
        # Write results
        with open(output_file, "a") as f:
            f.write(f"{base_name}\t{formatted_residues}\n")
        
        # Print statistics
        print(f"Processed: {base_name}")
        print(f"  Trp ligand atoms: {len(ligand_coords)}")
        print(f"  Protein atoms: {len(protein_coords)}")
        print(f"  Interacting residues: {len(sorted_residues)}")
        
    except Exception as e:
        print(f"Error processing {base_name}: {str(e)}")
        with open(output_file, "a") as f:
            f.write(f"{base_name}\tError: {str(e)}\n")
    finally:
        # Clean up PyMOL objects
        cmd.delete(base_name)

def batch_process_directory(input_dir, output_file, distance_cutoff=5.0, k=100):
    """
    Batch process all PDB files in a directory
    
    Parameters:
    input_dir: Input directory path
    output_file: Output file path
    distance_cutoff: Distance threshold in Angstroms
    k: Maximum number of neighbors to consider
    """
    # Create output file with header
    with open(output_file, "w") as f:
        f.write("# Protein\tInteracting_Residues\n")
    
    # Collect all PDB files
    pdb_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".pdb"):
                pdb_files.append(os.path.join(root, file))
    
    total = len(pdb_files)
    if total == 0:
        print(f"No PDB files found in {input_dir}")
        return
    
    print(f"Processing {total} PDB files with Trp ligands...")
    print(f"Using distance cutoff: {distance_cutoff}?")
    print(f"Max neighbors per atom: {k}")
    
    # Process each file
    for i, pdb_file in enumerate(pdb_files):
        print(f"\nProcessing file [{i+1}/{total}]: {os.path.basename(pdb_file)}")
        analyze_trp_interactions(pdb_file, output_file, distance_cutoff, k)
    
    print(f"\nProcessing complete! Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze protein interactions with Trp ligands (UNK residues) using KDTree method")
    parser.add_argument("-i", "--input", required=True, 
                        help="Input PDB file or directory path")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file path")
    parser.add_argument("-d", "--distance", type=float, default=5.0,
                        help="Distance cutoff in Angstroms (default: 5.0)")
    parser.add_argument("-k", type=int, default=100,
                        help="Maximum number of neighbors to consider (default: 100)")
    
    args = parser.parse_args()
    
    # Check if input is file or directory
    if os.path.isdir(args.input):
        batch_process_directory(args.input, args.output, args.distance, args.k)
    elif os.path.isfile(args.input) and args.input.endswith(".pdb"):
        # Process single file
        with open(args.output, "w") as f:
            f.write("# Protein\tInteracting_Residues\n")
        analyze_trp_interactions(args.input, args.output, args.distance, args.k)
    else:
        print(f"Error: Input path {args.input} is not a valid PDB file or directory")
        return
    
    print(f"Processing complete! Results saved to {args.output}")

if __name__ == "__main__":
    cmd.feedback("disable", "all", "everything")
    main()