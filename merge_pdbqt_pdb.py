import argparse
from pymol import cmd

def merge_pdb_files(file1, file2, output_file):
    cmd.delete("all")
    cmd.load(file1, "3m1v")
    cmd.load(file2, "COM")
    cmd.create("3m1v_COM", "3m1v COM")
    cmd.save(output_file, "3m1v_COM")

# Create the argument parser
parser = argparse.ArgumentParser(description='Merge PDB files')
parser.add_argument('file1', help='Path to the first PDB file')
parser.add_argument('file2', help='Path to the second PDB file')
parser.add_argument('--out', dest='output_file', help='Output file name')

# Parse the command-line arguments
args = parser.parse_args()

# Call the merge_pdb_files function with the provided arguments
merge_pdb_files(args.file1, args.file2, args.output_file)

