#!/usr/bin/env python

import os
import argparse
import time
from pymol import cmd, stored

def calculate_pocket_volume(pdb_path, distance=5.0):
    """Calculate pocket volume around ligand"""
    try:
        # Load PDB file
        cmd.reinitialize()
        cmd.load(pdb_path, "complex")
        
        # Identify ligand
        cmd.select("ligand", "hetatm and not (solvent or inorganic)")
        ligand_count = cmd.count_atoms("ligand")
        
        if ligand_count == 0:
            print(f"Warning: No ligand detected in {os.path.basename(pdb_path)}")
            return 0.0, 0.0
        
        # Define pocket around ligand
        pocket_sel = f"byres ligand around {distance} and polymer"
        cmd.select("pocket", pocket_sel)
        
        # Create molecular surface for the pocket
        cmd.show("surface", "pocket")
        
        # Calculate pocket volume using PyMOL's measurement
        stored.volume = 0.0
        cmd.iterate("pocket", "stored.volume += vdw")
        approx_volume = stored.volume
        
        # Calculate solvent accessible surface area
        cmd.set("dot_solvent", 1)
        cmd.set("dot_density", 3)
        sasa = cmd.get_area("pocket")
        
        return approx_volume, sasa
    
    except Exception as e:
        print(f"Error processing {os.path.basename(pdb_path)}: {str(e)}")
        return 0.0, 0.0
    finally:
        cmd.delete("all")

def process_pdb_folder(input_dir, output_file, distance=5.0):
    """Process all PDB files in directory"""
    start_time = time.time()
    results = []
    processed = 0
    no_ligand_count = 0
    
    print(f"Processing directory: {input_dir}")
    print(f"Using pocket distance: {distance} ?")
    
    for filename in sorted(os.listdir(input_dir)):
        if filename.lower().endswith(".pdb"):
            pdb_path = os.path.join(input_dir, filename)
            approx_vol, sasa = calculate_pocket_volume(pdb_path, distance)
            results.append((filename, approx_vol, sasa))
            processed += 1
            
            if approx_vol == 0.0:
                no_ligand_count += 1
                
            if processed % 10 == 0:
                print(f"Processed {processed} files...")
    
    # Write results
    with open(output_file, "w") as f:
        f.write("PDB_File\tApprox_Volume(A3)\tSASA(A2)\n")
        for filename, approx_vol, sasa in results:
            f.write(f"{filename}\t{approx_vol:.2f}\t{sasa:.2f}\n")
    
    # Performance stats
    elapsed = time.time() - start_time
    avg_time = elapsed / processed if processed > 0 else 0
    
    print("\n" + "="*50)
    print(f"Processing complete! Total files: {processed}")
    print(f"Files with ligands: {processed - no_ligand_count}")
    print(f"Files without ligands: {no_ligand_count}")
    print(f"Total time: {elapsed:.1f} sec | Avg: {avg_time:.1f} sec/file")
    print(f"Results saved to: {output_file}")
    print("="*50)

if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(
        description="Batch calculate pocket volumes for ligand-bound PDB structures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, 
                        help="Input directory containing PDB files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output results file")
    parser.add_argument("-d", "--distance", type=float, default=5.0,
                        help="Distance around ligand to define pocket (?)")
    args = parser.parse_args()
    
    # Initialize PyMOL in headless mode
    cmd.feedback("disable", "all", "everything")
    
    # Process files
    process_pdb_folder(args.input, args.output, args.distance)