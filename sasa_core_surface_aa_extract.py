#!/usr/bin/env python
from pymol import cmd
import os
import argparse
import sys
import re
from collections import defaultdict
import traceback

# Amino acid mapping
AA_CODES = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'
}

# Pre-calculated max SASA values from helical reference model
HELIX_MAX_SASA = {
    'A': 121.5, 'R': 258.2, 'N': 183.7, 'D': 181.4,
    'C': 157.1, 'Q': 211.8, 'E': 209.5, 'G': 97.8,
    'H': 210.8, 'I': 185.3, 'L': 189.0, 'K': 222.1,
    'M': 210.8, 'F': 225.9, 'P': 149.6, 'S': 145.8,
    'T': 161.8, 'W': 268.1, 'Y': 247.3, 'V': 163.7
}

def calculate_sasa_and_classify(pdb_file, sasa_threshold):
    base_name = os.path.basename(pdb_file).replace(".pdb", "")
    
    try:
        cmd.load(pdb_file, base_name)
        cmd.remove("solvent")
        cmd.remove("hetatm and not resn HOH")
        cmd.select("protein", "polymer")
        
        atom_count = cmd.count_atoms("protein")
        if atom_count == 0:
            print(f"   Warning: No protein atoms selected")
            return None, None
        
        cmd.get_area("protein", load_b=1)
        
        residues = cmd.get_model("protein").atom
        residue_dict = defaultdict(list)
        for atom in residues:
            # Use residue number as string to handle insertion codes
            resi_str = str(atom.resi)
            key = (atom.chain, resi_str, atom.resn)
            residue_dict[key].append(atom)
        
        residue_classification = {}
        residue_info = {}
        
        for (chain, resi, resn), atoms in residue_dict.items():
            if resn not in AA_CODES:
                continue
                
            total_sasa = sum(atom.b for atom in atoms)
            aa_code = AA_CODES[resn]
            max_sasa = HELIX_MAX_SASA.get(aa_code, 150.0)
            
            if max_sasa > 0:
                relative_sasa = (total_sasa / max_sasa) * 100
            else:
                relative_sasa = 0
            
            # Classify residue using threshold
            if relative_sasa > sasa_threshold:
                classification = "surface"
            else:
                classification = "core"
            
            # Use residue number as string
            residue_classification[(chain, resi)] = classification
            residue_info[(chain, resi)] = aa_code
        
        return residue_classification, residue_info
    except Exception as e:
        print(f"Error processing {base_name}: {str(e)}")
        traceback.print_exc()
        return None, None
    finally:
        try:
            cmd.delete(base_name)
        except:
            pass
        try:
            cmd.delete("protein")
        except:
            pass

def extract_chain_id(align_file):
    """Extract chain ID from alignment file"""
    try:
        with open(align_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(" "):
                    continue
                
                # Match target protein line
                parts = line.split()
                if len(parts) > 3 and '_' in parts[-1]:
                    target_id = parts[-1]
                    chain_id = target_id.split('_')[-1]
                    return chain_id
        
        # If not found, try to extract from file name
        base_name = os.path.basename(align_file)
        if '_' in base_name:
            parts = base_name.split('_')
            if len(parts) > 1:
                return parts[-1].replace(".txt", "")
        
        print("  Warning: Chain ID not found, using default 'A'")
        return 'A'
    except Exception as e:
        print(f"Error extracting chain ID: {str(e)}")
        return 'A'

def parse_alignment(align_file):
    """Parse alignment file and extract residue mapping"""
    ref_seq = []
    tar_seq = []
    ref_positions = []
    tar_positions = []
    
    try:
        with open(align_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header lines and summary lines
        # Look for the start of alignment blocks
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Check if this is an alignment block
            if line and line[0].isdigit():
                # This is a reference protein line
                ref_parts = re.split(r'\s+', line)
                ref_parts = [p for p in ref_parts if p]
                
                # Skip if not enough parts
                if len(ref_parts) < 4:
                    i += 1
                    continue
                
                # The next line is the alignment line
                if i + 1 >= len(lines):
                    i += 1
                    continue
                ali_line = lines[i+1].strip()
                
                # The line after that is the target protein line
                if i + 2 >= len(lines):
                    i += 3
                    continue
                tar_line = lines[i+2].strip()
                tar_parts = re.split(r'\s+', tar_line)
                tar_parts = [p for p in tar_parts if p]
                
                # Skip if not enough parts
                if len(tar_parts) < 4:
                    i += 3
                    continue
                
                # Extract sequence segments
                ref_seq_str = ref_parts[1]
                tar_seq_str = tar_parts[1]
                
                # Extract start positions (handle formats like "1-10" or "1")
                try:
                    ref_start_str = ref_parts[0]
                    if '-' in ref_start_str:
                        ref_start = int(ref_start_str.split('-')[0])
                    else:
                        ref_start = int(ref_start_str)
                except ValueError:
                    ref_start = 1
                
                try:
                    tar_start_str = tar_parts[0]
                    if '-' in tar_start_str:
                        tar_start = int(tar_start_str.split('-')[0])
                    else:
                        tar_start = int(tar_start_str)
                except ValueError:
                    tar_start = 1
                
                # Process this alignment block
                ref_pos = ref_start
                tar_pos = tar_start
                for j in range(len(ref_seq_str)):
                    ref_aa = ref_seq_str[j]
                    tar_aa = tar_seq_str[j]
                    
                    if ref_aa != '-':
                        ref_seq.append(ref_aa)
                        ref_positions.append(ref_pos)
                        ref_pos += 1
                    else:
                        ref_seq.append('-')
                        ref_positions.append(None)
                        
                    if tar_aa != '-':
                        tar_seq.append(tar_aa)
                        tar_positions.append(tar_pos)
                        tar_pos += 1
                    else:
                        tar_seq.append('-')
                        tar_positions.append(None)
                
                # Skip to next block
                i += 3
            else:
                i += 1
    except Exception as e:
        print(f"  Error parsing alignment file: {str(e)}")
        traceback.print_exc()
    
    # Print some debug info
    print(f"  Parsed alignment: {len(ref_seq)} positions")
    print(f"  Reference sequence: {''.join(ref_seq[:20])}...")
    print(f"  Target sequence: {''.join(tar_seq[:20])}...")
    
    return ref_seq, tar_seq, ref_positions, tar_positions

def main(pdb_dir, align_dir, core_output, surface_output, sasa_threshold):
    core_results = []
    surface_results = []
    
    print(f"Processing PDBs in: {pdb_dir}")
    print(f"Using alignments from: {align_dir}")
    print(f"Relative SASA threshold: {sasa_threshold}%")
    
    pdb_files = []
    for root, dirs, files in os.walk(pdb_dir):
        for file in files:
            if file.endswith(".pdb"):
                pdb_files.append(os.path.join(root, file))
    
    total_count = len(pdb_files)
    if total_count == 0:
        print(f"No PDB files found in {pdb_dir}")
        return
    
    print(f"Processing {total_count} PDB files...")
    
    processed_count = 0
    failed_count = 0
    
    for i, pdb_file in enumerate(pdb_files):
        pdb_id = os.path.basename(pdb_file).replace(".pdb", "")
        align_file = os.path.join(align_dir, f"2c44_{pdb_id}.txt")
        
        if not os.path.exists(align_file):
            print(f"Warning: Alignment file not found for {pdb_id}")
            failed_count += 1
            continue
            
        print(f"[{i+1}/{total_count}] Processing: {pdb_id}")
        
        try:
            # Get chain ID from alignment file
            chain_id = extract_chain_id(align_file)
            print(f"  Using chain ID: {chain_id}")
            
            # Calculate SASA and classify residues
            residue_class, residue_info = calculate_sasa_and_classify(pdb_file, sasa_threshold)
            if residue_class is None or residue_info is None:
                print(f"   Error: Failed to classify residues for {pdb_id}")
                failed_count += 1
                continue
            
            # Parse alignment
            ref_seq, tar_seq, ref_pos, tar_pos = parse_alignment(align_file)
            
            # Check if we got valid alignment data
            if len(ref_seq) == 0 or len(tar_seq) == 0:
                print(f"  Warning: No alignment data found for {pdb_id}")
                failed_count += 1
                continue
            
            # Collect residue pairs
            core_pairs = []
            surface_pairs = []
            
            for idx in range(len(ref_seq)):
                r_aa = ref_seq[idx]
                t_aa = tar_seq[idx]
                r_pos = ref_pos[idx]
                t_pos = tar_pos[idx]
                
                # Skip gaps or missing positions
                if r_aa == '-' or t_aa == '-' or r_pos is None or t_pos is None:
                    continue
                
                # Create key for residue lookup
                key = (chain_id, str(t_pos))
                
                # Find residue classification
                classification = None
                if key in residue_class:
                    classification = residue_class[key]
                else:
                    # Try to find residue by position only
                    for (ch, resi), cls in residue_class.items():
                        if resi == str(t_pos):
                            classification = cls
                            break
                
                if classification is None:
                    continue
                
                # Get target amino acid from structure
                t_aa_struct = residue_info.get(key, None)
                if t_aa_struct is None:
                    # Try to find residue by position only
                    for (ch, resi), aa in residue_info.items():
                        if resi == str(t_pos):
                            t_aa_struct = aa
                            break
                
                if t_aa_struct is None:
                    continue
                
                # Add to appropriate list
                if classification == "core":
                    core_pairs.append((r_aa, t_aa_struct))
                else:
                    surface_pairs.append((r_aa, t_aa_struct))
            
            # Calculate similarities
            core_sim = 0.0
            if core_pairs:
                matches = sum(1 for r, t in core_pairs if r == t)
                core_sim = (matches / len(core_pairs)) * 100
            
            surface_sim = 0.0
            if surface_pairs:
                matches = sum(1 for r, t in surface_pairs if r == t)
                surface_sim = (matches / len(surface_pairs)) * 100
            
            core_results.append((pdb_id, core_sim))
            surface_results.append((pdb_id, surface_sim))
            processed_count += 1
            
            print(f"  Core similarity: {core_sim:.1f}% ({len(core_pairs)} residues)")
            print(f"  Surface similarity: {surface_sim:.1f}% ({len(surface_pairs)} residues)")
            
        except Exception as e:
            print(f"Error processing {pdb_id}: {str(e)}")
            traceback.print_exc()
            failed_count += 1
    
    # Write results
    with open(core_output, "w") as f:
        f.write("id\tcor_similarity%\n")
        for pdb_id, sim in core_results:
            f.write(f"{pdb_id}\t{sim:.1f}%\n")
    
    with open(surface_output, "w") as f:
        f.write("id\tsurface_similarity%\n")
        for pdb_id, sim in surface_results:
            f.write(f"{pdb_id}\t{sim:.1f}%\n")
    
    print(f"\nProcessing complete!")
    print(f"Core similarity results: {core_output}")
    print(f"Surface similarity results: {surface_output}")
    print(f"Total processed: {processed_count} of {total_count} files")
    print(f"Failed: {failed_count} files")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate core and surface sequence similarity")
    parser.add_argument("-i", "--input", required=True, dest="ipdb",
                        help="Input directory with PDB files")
    parser.add_argument("-ialn", required=True, 
                        help="Directory with alignment files")
    parser.add_argument("-oc", required=True,
                        help="Output file for core similarity")
    parser.add_argument("-os", required=True,
                        help="Output file for surface similarity")
    parser.add_argument("-t", "--threshold", type=float, default=60.0,
                        help="Relative SASA threshold (default 60%)")
    
    args = parser.parse_args()
    
    if sys.version_info < (3, 0):
        print("Error: Requires Python 3 or higher")
        sys.exit(1)
    
    cmd.feedback("disable", "all", "everything")
    
    main(
        pdb_dir=args.ipdb,
        align_dir=args.ialn,
        core_output=args.oc,
        surface_output=args.os,
        sasa_threshold=args.threshold
    )