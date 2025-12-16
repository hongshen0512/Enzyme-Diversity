#!/usr/bin/env python
import os
import subprocess
from glob import glob

MGLTOOLS_PATH = "/gss1/App_os7/miniconda3/envs/autodock" 
PDB_DIR = "/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/pdb_aligned"
PDBQT_DIR = "/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/pdbqt"
LOG_DIR = "/gss1/home/shenhong/AF3_TNA_IAD/autodock/processed/logs"

os.makedirs(PDBQT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

def prepare_receptor(pdb_file):
    base_name = os.path.basename(pdb_file).replace(".pdb", "")
    log_file = os.path.join(LOG_DIR, f"{base_name}_prepare.log")
    pdbqt_file = os.path.join(PDBQT_DIR, f"{base_name}.pdbqt")
    
    cmd = [
        os.path.join(MGLTOOLS_PATH, "bin/pythonsh"),
        os.path.join(MGLTOOLS_PATH, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"),
        "-r", pdb_file,
        "-o", pdbqt_file,
        "-A", "bonds_hydrogens",
        "-U", "nphs_lps_waters_nonstdres",
        "-v" 
    ]
    
    try:
        with open(log_file, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
        print(f"ok: {base_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"failed: {base_name} - error {e.returncode}")
        return False

def batch_preprocess():
    pdb_files = glob(os.path.join(PDB_DIR, "*.pdb"))
    total = len(pdb_files)
    success_count = 0
    
    print(f"start {total} pdb file.")
    
    for i, pdb_file in enumerate(pdb_files):
        print(f"\ndeal [{i+1}/{total}]: {os.path.basename(pdb_file)}")
        if prepare_receptor(pdb_file):
            success_count += 1
    
    print(f"\nsuccess: {success_count}/{total}")
    return success_count == total

if __name__ == "__main__":
    batch_preprocess()