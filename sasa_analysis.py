# -*- coding: latin-1 -*-
#!/usr/bin/env python
from pymol import cmd
import os
import argparse
import sys
from collections import defaultdict, Counter

# Amino acid mapping
AA_CODES = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'
}

def calculate_sasa_and_classify(pdb_file, threshold):
    base_name = os.path.basename(pdb_file).replace(".pdb", "")

    try:
        cmd.load(pdb_file, base_name)
        cmd.remove("solvent")
        cmd.remove("hetatm and not resn HOH")
        cmd.select("protein", "polymer")

        if cmd.count_atoms("protein") == 0:
            print(f"   Warning: No protein atoms selected")
            return None

        cmd.get_area("protein", load_b=1)

        core_aa = Counter()
        surface_aa = Counter()

        residues = cmd.get_model("protein").atom
        residue_dict = defaultdict(list)
        for atom in residues:
            key = (atom.chain, atom.resi, atom.resn)
            residue_dict[key].append(atom)

        for (chain, resi, resn), atoms in residue_dict.items():
            if resn not in AA_CODES:
                continue
            total_sasa = sum(atom.b for atom in atoms)
            aa_code = AA_CODES[resn]

            if total_sasa < threshold:
                core_aa[aa_code] += 1
            else:
                surface_aa[aa_code] += 1

        return {"core": dict(core_aa), "surface": dict(surface_aa)}
    except Exception as e:
        print(f"Error processing {base_name}: {e}")
        return None
    finally:
        cmd.delete(base_name)
        cmd.delete("protein")

def format_aa_counts(aa_dict):
    if not aa_dict:
        return "NoResidues"
    return " ".join([f"{aa}:{count}" for aa, count in sorted(aa_dict.items())])

def batch_process_pdb_files(input_dir, core_output, surface_output, threshold):
    for output_file in [core_output, surface_output]:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

    pdb_files = [os.path.join(root, f)
                 for root, _, files in os.walk(input_dir)
                 for f in files if f.endswith(".pdb")]

    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return

    print(f"Processing {len(pdb_files)} PDB files...")
    print(f"Threshold for core (<) vs surface (¡Ý): {threshold}")

    with open(core_output, "w") as core_f, open(surface_output, "w") as surf_f:
        core_f.write("# Protein Core Amino Acid Composition\n")
        core_f.write("# Format: Protein AA1:count AA2:count ...\n")

        surf_f.write("# Protein Surface Amino Acid Composition\n")
        surf_f.write("# Format: Protein AA1:count AA2:count ...\n")

        processed = failed = 0
        for i, pdb_file in enumerate(pdb_files, 1):
            base_name = os.path.basename(pdb_file).replace(".pdb", "")
            print(f"[{i}/{len(pdb_files)}] Processing: {base_name}")
            aa_counts = calculate_sasa_and_classify(pdb_file, threshold)
            if aa_counts is None:
                failed += 1
                core_f.write(f"{base_name} NoResidues\n")
                surf_f.write(f"{base_name} NoResidues\n")
                continue
            core_f.write(f"{base_name} {format_aa_counts(aa_counts['core'])}\n")
            surf_f.write(f"{base_name} {format_aa_counts(aa_counts['surface'])}\n")
            processed += 1

    print("\nProcessing complete!")
    print(f"Successfully processed: {processed} files")
    print(f"Failed: {failed} files")
    print(f"Core composition: {core_output}")
    print(f"Surface composition: {surface_output}")

def main():
    parser = argparse.ArgumentParser(
        description="Calculate SASA and classify residues into core/surface")
    parser.add_argument("-i", "--input", required=True,
                        help="Input directory with PDB files")
    parser.add_argument("-oc", "--core_output", required=True,
                        help="Core residue output file")
    parser.add_argument("-os", "--surface_output", required=True,
                        help="Surface residue output file")
    parser.add_argument("-t", "--threshold", type=float, default=25.0,
                        help="SASA threshold (?2) between core and surface (default: 25)")

    args = parser.parse_args()

    if sys.version_info < (3, 0):
        print("Error: Requires Python 3 or higher")
        sys.exit(1)

    cmd.feedback("disable", "all", "everything")
    batch_process_pdb_files(args.input, args.core_output,
                            args.surface_output, args.threshold)

if __name__ == "__main__":
    main()