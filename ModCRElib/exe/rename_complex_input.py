import os
import sys
import configparser
import argparse

# Read configuration file #
config = configparser.ConfigParser()
exe_path = os.path.abspath(os.path.dirname(__file__))

if os.path.exists(os.path.join(exe_path, "..", "ModCRElib")):
    scripts_path = os.path.join(exe_path, "..")
elif os.path.exists(os.path.join(exe_path, "..", "..", "ModCRElib")):
    scripts_path = os.path.join(exe_path, "..", "..")
elif os.path.exists(os.path.join(exe_path, "..", "..", "..", "ModCRElib")):
    scripts_path = os.path.join(exe_path, "..", "..", "..")
else:
    scripts_path = os.path.join(exe_path)

config_file = os.path.join(scripts_path, "ModCRElib", "configure/config.ini")
config.read(config_file)

# Get sbilib path #
sbilib_path = config.get("Paths", "sbilib_path")
sys.path.append(sbilib_path)

from SBILib.structure import PDB

def main():
    # Set up argparse interface
    parser = argparse.ArgumentParser(
        description="Rename a PDB file to match ModCRElib complex-building format and delete the original."
    )
    parser.add_argument(
        "pdb_file", 
        type=str, 
        help="Path to the input PDB file that needs to be renamed."
    )
    
    args = parser.parse_args()
    original_file = args.pdb_file

    # Verify the file actually exists before processing
    if not os.path.exists(original_file):
        print(f"ERROR: File not found at {original_file}")
        sys.exit(1)

    # Print initial input information for the CLI wrapper
    print(f"--- Processing Input File ---")
    print(f"Input Path: {os.path.abspath(original_file)}")
    print(f"File Size: {os.path.getsize(original_file)} bytes")

    # Load the structure
    protein = PDB(original_file)

    # Extract chain and residue boundaries
    begin = None
    end = None
    for x in protein.proteins:
        begin = x.first_aminoacid.number
        end = x.last_aminoacid.number

    # Prevent crash if sequence boundaries could not be parsed
    if begin is None or end is None:
        print("ERROR: Could not determine residue sequence boundaries from the PDB file.")
        sys.exit(1)

    print(f"Detected Sequence Bounds: Residues {begin} to {end}")

    # Construct the new structured file name
    new_file = original_file.replace(".pdb", f":{begin}:{end}_TF.pdb")
    print(f"Target Output Path: {os.path.abspath(new_file)}")

    try:
        # 1. Write out the newly formatted PDB file
        protein.write(new_file)
        print("SUCCESS: New PDB file successfully written.")
        
        # 2. Automatically delete the old file after a successful write
        os.remove(original_file)
        print(f"CLEANUP: Original file '{os.path.basename(original_file)}' has been deleted.")
        print("-----------------------------")

    except Exception as e:
        print(f"CRITICAL ERROR during file operations: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

