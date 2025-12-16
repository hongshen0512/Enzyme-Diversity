#!/bin/bash

cd /gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/pdb_aligned
source /gss1/env_os7/python_v3.11.0.env
source /gss1/env_os7/pymol.env

BASE_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock"
SCRIPTS_DIR="${BASE_DIR}/scripts"
RESULTS_DIR="${BASE_DIR}/processed/autodock_results"
MERGED_DIR="${BASE_DIR}/processed/merged_pdb"
KEYAA_DIR="${BASE_DIR}/processed/key_aa"

find . -maxdepth 1 -type f -name "*.pdb" | while read -r pdb_file; do
    base_name=$(basename "$pdb_file" .pdb)
    
    result_file="${RESULTS_DIR}/${base_name}_result.pdb"
    output_file="${MERGED_DIR}/${base_name}_trp.pdb"
    
    if [ ! -f "$result_file" ]; then
        echo "  no - $result_file"
        continue
    fi
    
    echo " merge: $pdb_file and ${base_name}_result.pdb"
    python "${SCRIPTS_DIR}/merge_pdbqt_pdb.py" \
        "$result_file" \
        "$pdb_file" \
        --out "$output_file"
    
    if [ ! -f "$output_file" ]; then
        echo " failed - $output_file"
        continue
    fi
    
    keyaa_file="${KEYAA_DIR}/${base_name}_trp_key_aa.txt"
    
    echo " analyze: ${base_name}_trp.pdb"
    python "${SCRIPTS_DIR}/pymol_interaction_4.5_sum.py" \
        -i "$output_file" \
        -o "$keyaa_file"
    
    if [ ! -f "$keyaa_file" ]; then
        echo " failed - $keyaa_file"
    else
        residue_count=$(awk -F'\t' 'NR>1 {print $2}' "$keyaa_file" | tr ',' '\n' | wc -l)
        echo " ok find $residue_count aa"
    fi
    
    echo "----------------------------------------"
done

echo "finish"