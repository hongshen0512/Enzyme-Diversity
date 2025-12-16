#!/bin/bash

source /gss1/env_os7/autodock.env

CONFIG_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/configs"
OUTPUT_DIR="/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/autodock_results"

mkdir -p $OUTPUT_DIR

for config_file in $CONFIG_DIR/*.txt; do
    base_name=$(basename "$config_file" .txt)
    output_file="$OUTPUT_DIR/${base_name}_result.pdbqt"
    
    echo "vina: $base_name"
    vina --config "$config_file" --out "$output_file"
    
    if [ -s "$output_file" ]; then
        echo "sucess: $base_name"
    else
        echo "failed: $base_name - check $log_file"
    fi
done

echo "ok"

