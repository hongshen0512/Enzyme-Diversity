#!/usr/bin/env python
from pymol import cmd
import argparse
import os

aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'
}

def analyze_interactions(input_file, output_file):
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    
    cmd.load(input_file, base_name)
    
    cmd.remove("solvent")
    
    cmd.select("ligand", "hetatm")
    
    cmd.select("nearby_residues", "byres ligand around 3.5")
    
    nearby_atoms = cmd.get_model("nearby_residues")
    
    interacting_residues = set()
    
    for atom in nearby_atoms.atom:
        if atom.hetatm:
            continue
        
        if atom.resn in aa_codes:
            res_id = f"{atom.chain}:{aa_codes[atom.resn]}{atom.resi}"
            interacting_residues.add(res_id)
    
    sorted_residues = sorted(
        interacting_residues, 
        key=lambda x: (x.split(':')[0], int(''.join(filter(str.isdigit, x))))
    )
    
    formatted_residues = ", ".join([f"'{res}'" for res in sorted_residues])
    
    with open(output_file, "a") as f:
        f.write(f"{base_name}\t{formatted_residues}\n")
    
    cmd.delete(base_name)
    cmd.delete("ligand")
    cmd.delete("nearby_residues")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze protein ligand interactions and extract residues within the 3.5 ? range")
    parser.add_argument("-i", "--input", required=True, 
                        help="input PDB file path")
    parser.add_argument("-o", "--output", required=True,
                        help="output path")
    
    args = parser.parse_args()
    
    with open(args.output, "w") as f:
        f.write("# Protein\tInteracting_Residues\n")
    
    analyze_interactions(args.input, args.output)
    
    print(f"finished saved to {args.output}")

if __name__ == "__main__":
    cmd.feedback("disable", "all", "everything")
    main()