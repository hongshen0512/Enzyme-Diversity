#!/bin/bash

REF_CONFIG="/gss1/home/shenhong/AF3_TNA_IAD/autodock/input/reference/9blv1.txt"
CONFIG_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/configs"
PDBQT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/pdbqt"

mkdir -p $CONFIG_DIR

center_x=$(grep 'center_x' $REF_CONFIG | cut -d'=' -f2 | tr -d ' ')
center_y=$(grep 'center_y' $REF_CONFIG | cut -d'=' -f2 | tr -d ' ')
center_z=$(grep 'center_z' $REF_CONFIG | cut -d'=' -f2 | tr -d ' ')
size_x=$(grep 'size_x' $REF_CONFIG | cut -d'=' -f2 | tr -d ' ')
size_y=$(grep 'size_y' $REF_CONFIG | cut -d'=' -f2 | tr -d ' ')
size_z=$(grep 'size_z' $REF_CONFIG | cut -d'=' -f2 | tr -d ' ')

for pdbqt_file in $PDBQT_DIR/*.pdbqt; do
    base_name=$(basename "$pdbqt_file" .pdbqt)
    config_file="$CONFIG_DIR/${base_name}.txt"
    
    cat > "$config_file" <<EOF
receptor = $pdbqt_file
ligand = /gss1/home/shenhong/AF3_TNA_IAD/autodock/input/reference/trp.pdbqt

center_x = $center_x
center_y = $center_y
center_z = $center_z

size_x = $size_x
size_y = $size_y
size_z = $size_z

exhaustiveness = 8
num_modes = 10
energy_range = 4
EOF

    echo "ok: $config_file"
done

echo "finished"