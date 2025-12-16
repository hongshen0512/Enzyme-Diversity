#!/bin/bash
source /gss1/env_os7/autodock.env
source /gss1/env_os7/rdkit.env
source /gss1/env_os7/pymol.env

INPUT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/11_clades_ref_pdb"
REFERENCE_PROTEIN="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/2c44.pdb" 
OUTPUT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/11_ref_2c44_pdb"
LOG_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/logs"
PYMOL_SCRIPT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/pymol_scripts"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$PYMOL_SCRIPT_DIR"

if [ ! -f "$REFERENCE_PROTEIN" ]; then
    echo "error: Reference protein $REFERENCE_PROTEIN not found!"
    exit 1
fi

ref_name=$(basename "$REFERENCE_PROTEIN" .pdb)
total_files=$(ls "$INPUT_DIR"/*.pdb | wc -l)
count=0
success_count=0
fail_count=0

for input_file in "$INPUT_DIR"/*.pdb; do
    ((count++))
    base_name=$(basename "$input_file" .pdb)
    output_pdb="$OUTPUT_DIR/$base_name.pdb"
    log_file="$LOG_DIR/${base_name}_convert.log"
    pymol_script="$PYMOL_SCRIPT_DIR/${base_name}_align.py"
    
    echo "Processing [$count/$total_files]: $base_name.pdb"
    echo "Aligning to reference: $ref_name"
    
    cat > "$pymol_script" <<EOF
import os
import sys

cmd.load("$REFERENCE_PROTEIN", "ref_protein")
cmd.load("$input_file", "current_structure")

alignment = cmd.align(
    "current_structure and name CA", 
    "ref_protein and name CA",
    cycles=5,
    cutoff=2.0
)

print(f"Alignment Results:")
print(f"RMSD = {alignment[0]:.4f}")
print(f"Aligned Atoms = {alignment[1]}")
print(f"Rotation Matrix = {alignment[2]}")

cmd.save("$output_pdb", "current_structure", format="pdb")
cmd.quit()
EOF
    
    echo "Running PyMOL alignment..."
    pymol -c -u "$pymol_script" > "$log_file" 2>&1
    
    if [ -s "$output_pdb" ]; then
        rmsd=$(grep "RMSD = " "$log_file" | awk '{print $3}')
        atoms_aligned=$(grep "Aligned Atoms = " "$log_file" | awk '{print $4}')
        
        echo "Success: Saved aligned structure to $base_name.pdb"
        echo "RMSD = ${rmsd}, Aligned Atoms = $atoms_aligned"
        ((success_count++))
    else
        echo "ERROR: Alignment failed! Check $log_file"
        ((fail_count++))
    fi
done

summary_file="$LOG_DIR/alignment_summary.csv"
echo "protein,status,RMSD,aligned_atoms" > "$summary_file"
for input_file in "$INPUT_DIR"/*.pdb; do
    base_name=$(basename "$input_file" .pdb)
    log_file="$LOG_DIR/${base_name}_convert.log"
    output_pdb="$OUTPUT_DIR/$base_name.pdb"
    
    if [ -f "$output_pdb" ]; then
        rmsd=$(grep "RMSD = " "$log_file" | awk '{print $3}')
        atoms=$(grep "Aligned Atoms = " "$log_file" | awk '{print $4}')
        echo "$base_name,success,$rmsd,$atoms" >> "$summary_file"
    else
        echo "$base_name,failed,NA,NA" >> "$summary_file"
    fi
done

echo ""
echo "Processing completed!"
echo "Success: $success_count, Failed: $fail_count, Total: $total_files"
echo "Summary report saved to: $summary_file"