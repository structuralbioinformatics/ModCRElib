import sys
import os
import configparser
import argparse  # 1. Imported argparse


# Helper function to prepend a valid PDB header Needed cause modern dssp assumes pdb otherwise
def prepend_pdb_header(file_path):
    """
    Prepends a dummy HEADER line to the top of a PDB file.
    This triggers DSSP 4's standard PDB parsing fallback engine instead of crashing.
    """
    dummy_header = "HEADER    REMODELED PROTEIN SYSTEM\n"
    
    # Read the existing renumbered content
    with open(file_path, 'r') as f:
        content = f.read()
        
    # Overwrite the file starting with the header followed by the original content
    with open(file_path, 'w') as f:
        f.write(dummy_header + content)



# Read configuration file #
config = configparser.ConfigParser()
# Get scripts path (i.e. ".") #
exe_path = os.path.abspath(os.path.dirname(__file__))
if os.path.exists(os.path.join(exe_path,"..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..")
elif os.path.exists(os.path.join(exe_path,"..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..")
elif os.path.exists(os.path.join(exe_path,"..","..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..","..")
else:
   scripts_path = os.path.join(exe_path)


config_file  = os.path.join(scripts_path,"ModCRElib","configure/config.ini")
config.read(config_file)

# Get sbilib path #
sbilib_path = config.get("Paths", "sbilib_path")

sys.path.append(sbilib_path)
from SBILib.structure import PDB

# 2. Replaced sys.argv manual parsing with argparse
parser = argparse.ArgumentParser(
    description="Renumber PDB Models: This script takes all PDB files from a source directory, "
                "renumbers their residues and atoms starting from 1 for each chain, "
                "and saves the modified files into a destination directory."
)
parser.add_argument(
    "sourcedir", 
    help="Path to the directory containing the original .pdb files to be renumbered."
)
parser.add_argument(
    "destdir", 
    help="Path to the directory where the renumbered .pdb files will be saved."
)

# Parse the arguments
args = parser.parse_args()

# Map parsed arguments to variables
sourcedir = args.sourcedir
destdir = args.destdir

# 3. Print information about what the script is about to do
print(f"\n--- Running PDB Model Renumbering ---")
print(f"Source Directory:      {os.path.abspath(sourcedir)}")
print(f"Destination Directory: {os.path.abspath(destdir)}")
print(f"-------------------------------------\n")

# Verify source directory exists before running
if not os.path.exists(sourcedir):
    print(f"Error: Source directory '{sourcedir}' does not exist.")
    sys.exit(1)

# Ensure destination directory exists, create it if it doesn't
os.makedirs(destdir, exist_ok=True)

models = [f for f in os.listdir(sourcedir) if f.endswith(".pdb")]

if not models:
    print(f"No .pdb files found in {sourcedir}")
    sys.exit(0)

# Process the files
for x in models:
    fulldir = os.path.join(sourcedir, x)
    renamed = os.path.join(destdir, x)
    
    print(f"Processing: {x}...")
    protein = PDB(fulldir)
    for chainID in protein.chain_identifiers:
        protein.get_chain_by_id(chainID).renumber_residues(1)
        protein.get_chain_by_id(chainID).renumber_atoms(1)
    protein.write(renamed)
    prepend_pdb_header(renamed) # adding the pdb header to the new file

print("\n✓ Renumbering complete!")
