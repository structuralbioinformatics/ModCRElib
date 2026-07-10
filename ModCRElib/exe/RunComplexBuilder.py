import os
import sys
import subprocess
import argparse
import configparser

def run_command(command, stdout_file=None, stderr_to_stdout=False):
    """Safely executes a shell command and monitors for errors."""
    try:
        out_target = subprocess.PIPE
        err_target = subprocess.PIPE
        
        if stdout_file:
            out_target = open(stdout_file, 'w')
            if stderr_to_stdout:
                err_target = out_target

        print(f"Running: {' '.join(command)}")
        
        result = subprocess.run(
            command, 
            stdout=out_target, 
            stderr=err_target, 
            text=True, 
            check=True
        )
        
        if stdout_file:
            out_target.close()
        return result

    except subprocess.CalledProcessError as e:
        print(f"\nCRITICAL ERROR: Command failed with exit code {e.returncode}")
        if not stdout_file and e.stderr:
            print(f"Error details:\n{e.stderr}")
        sys.exit(1)

def get_fiber_binary_path():
    """Resolves the fiber binary path dynamically from the relative config.ini."""
    exe_path = os.path.dirname(os.path.abspath(__file__))
    config_file = os.path.normpath(os.path.join(exe_path, "..", "configure", "config.ini"))
    
    if not os.path.exists(config_file):
        print(f"ERROR: Configuration file not found at {config_file}")
        sys.exit(1)
        
    config = configparser.ConfigParser()
    config.read(config_file)
    
    try:
        # Assuming the config has a section like [Paths] containing fiber_path
        # Adjust section name ("Paths") or option name ("fiber_path") if needed
        return config.get("Paths", "x3dna_path")
    except (configparser.NoSectionError, configparser.NoOptionError):
        print(f"ERROR: Could not find 'fiber_path' under [Paths] in {config_file}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="ModCRE Structural Complex Pipeline")
    
    # User-defined Directories
    parser.add_argument(
        "-d", "--binary-dir", 
        type=str, 
        default="BinaryInteractions",
        help="Directory to store binary interactions and DNA sequence structure."
    )
    parser.add_argument(
        "-o", "--complex-dir", 
        type=str, 
        default="Complex",
        help="Directory to save the final modeled complex and logs."
    )
    parser.add_argument(
        "--pdb", type=str, default="./pdb", help="Path for the downloaded PDB database."
    )
    parser.add_argument(
        "--pbm", type=str, default="./pbm", help="Path for structural protein-binding matrix directory."
    )
    parser.add_argument(
        "--dummy", type=str, default="./dummy", help="Path for temporary dummy/placeholder directory."
    )
    
    # Input sequence
    parser.add_argument(
        "--seq", 
        type=str, 
        default="aaaaaaaaaaaaaCACGTCACCTTGaaaaaaaaaaaaa",
        help="DNA sequence string for x3dna fiber construction."
    )
    
    args = parser.parse_args()

    # 1. Dynamically read fiber binary location from relative configuration
    fiber_bin = os.path.join(get_fiber_binary_path(), "fiber") 

    # 2. Ensure all user-specified directories exist
    target_dirs = [args.binary_dir, args.complex_dir, args.pdb, args.pbm, args.dummy]
    for directory in target_dirs:
        os.makedirs(directory, exist_ok=True)

    # 3. Step 1: Run X3DNA Fiber Command
    # Outputs 'dnasequence.pdb' cleanly into the user's specific binary interactions folder
    dna_output = os.path.join(args.binary_dir, "dnasequence.pdb")
    fiber_cmd = [
        fiber_bin, 
        f"-seq={args.seq}", 
        "-b", dna_output
    ]
    
    print("\n[Step 1/2] Generating DNA structure using x3dna fiber...")
    run_command(fiber_cmd)
    print(f"SUCCESS: DNA PDB generated at {dna_output}")
    # 4. Step 2: Run Complex Builder Command
    # Dynamically find the absolute path to complexbuilder.py based on this script's location
    exe_dir = os.path.dirname(os.path.abspath(__file__))
    complexbuilder_abs_path = os.path.join(exe_dir, "complexbuilder.py")

    log_path = os.path.join(args.complex_dir, "complex.log")
    builder_cmd = [
        sys.executable, 
        complexbuilder_abs_path,  # Use the absolute path instead of a relative one
        "-d", args.binary_dir, 
        "-o", args.complex_dir,
    ]
    
    print("\n[Step 2/2] Running Complex Builder modeling...")
    run_command(builder_cmd, stdout_file=log_path, stderr_to_stdout=True)
    print(f"SUCCESS: Pipeline complete. Logs stored at {log_path}")

    

if __name__ == "__main__":
    main()

