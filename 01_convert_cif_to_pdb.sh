#!/bin/bash
source /gss1/env_os7/autodock.env
source /gss1/env_os7/rdkit.env
source /gss1/env_os7/pymol.env

INPUT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/361_cif"
REFERENCE_PROTEIN="/gss1/home/shenhong/AF3_TNA_IAD/autodock/input/reference/9blv1.pdbqt" 
OUTPUT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/pdb_aligned"
LOG_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/logs"
PYMOL_SCRIPT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/pymol_scripts"

mkdir -p $OUTPUT_DIR
mkdir -p $LOG_DIR

if [ ! -f "$REFERENCE_PROTEIN" ]; then
    echo "errpr: no $REFERENCE_PROTEIN"
    exit 1
fi

ref_name=$(basename "$REFERENCE_PROTEIN" .pdbqt)

total_files=$(ls $INPUT_DIR/*.cif | wc -l)
count=0
success_count=0
fail_count=0

for cif_file in $INPUT_DIR/*.cif; do
    ((count++))
    base_name=$(basename "$cif_file" .cif)
    pdb_file="$OUTPUT_DIR/$base_name.pdb"
    log_file="$LOG_DIR/${base_name}_convert.log"
    pymol_script="$PYMOL_SCRIPT_DIR/${base_name}_convert.py"
    
    echo "deal [$count/$total_files]: $base_name.cif"
    echo " align to: $ref_name"
    
    cat > "$pymol_script" <<EOF
import os
import sys

cmd.load("$REFERENCE_PROTEIN", "ref_protein")

cmd.load("$cif_file", "current_structure")

ref_selection = "ref_protein and name CA"
current_selection = "current_structure and name CA"

alignment = cmd.align(current_selection, ref_selection, cycles=5, cutoff=2.0)

print("align:")
print(f"RMSD = {alignment[0]:.3f} ?")
print(f"align atom no = {alignment[1]}")
print(f"matrix = {alignment[2]}")

cmd.save("$pdb_file", "current_structure", format="pdb")


cmd.quit()
EOF
    
    echo "  PyMOL covert.."
    pymol -c -u "$pymol_script" > "$log_file" 2>&1
    
    if [ ! -s "$pdb_file" ]; then
        echo "  error: failed check $log_file"
        ((fail_count++))
    else
        rmsd=$(grep "RMSD = " "$log_file" | awk '{print $3}')
        atoms_aligned=$(grep "align atom no = " "$log_file" | awk '{print $4}')
        
        echo "  ok: saved to $base_name.pdb"
        echo "  RMSD = ${rmsd}?, align atom no = $atoms_aligned"
        
        rm -f ${pdb_file%.pdb}.pse
        ((success_count++))
    fi
done

summary_file="$LOG_DIR/alignment_summary.csv"
echo "protein,state,RMSD, align atom no" > "$summary_file"
for cif_file in $INPUT_DIR/*.cif; do
    base_name=$(basename "$cif_file" .cif)
    log_file="$LOG_DIR/${base_name}_convert.log"
    
    if [ -f "$OUTPUT_DIR/$base_name.pdb" ]; then
        rmsd=$(grep "RMSD = " "$log_file" | awk '{print $3}')
        atoms_aligned=$(grep "align atom no = " "$log_file" | awk '{print $4}')
        echo "$base_name,ok,$rmsd,$atoms_aligned" >> "$summary_file"
    else
        echo "$base_name,failed,NA,NA" >> "$summary_file"
    fi
done

echo ""
echo "ok!"
echo "success: $success_count, failed: $fail_count, total: $total_files"
echo "text to: $summary_file"