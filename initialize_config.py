import os
import shutil
import configparser
from pathlib import Path

def find_x3dna_bottom_up(start_path):
    """
    Looks for X3DNA by searching downward from start_path, then walking up 
    through parent directories, and finally checking standard system directories.
    """
    # 1. Search downward from the current package directory first
    try:
        for path in start_path.glob("**/x3dna*"):
            if path.is_dir() and ((path / "bin" / "find_pair").exists() or (path / "bin" / "fiber").exists()):
                return path
    except (PermissionError, FileNotFoundError):
        pass

    # 2. Walk upwards through parent directories and scan inside them
    current = start_path.parent
    # Stop when we hit the root directory
    while current != current.parent:
        try:
            # Look for an x3dna folder inside this parent level (non-recursive to avoid huge sweeps)
            for path in current.glob("x3dna*"):
                if path.is_dir() and ((path / "bin" / "find_pair").exists() or (path / "bin" / "fiber").exists()):
                    return path
            
            # Also do a quick 1-level deep check in case it's in a sibling folder (e.g., ../software/x3dna)
            for path in current.glob("*/x3dna*"):
                if path.is_dir() and ((path / "bin" / "find_pair").exists() or (path / "bin" / "fiber").exists()):
                    return path
        except (PermissionError, FileNotFoundError):
            pass
        current = current.parent

    # 3. Global fallbacks if it's completely outside the user's project tree
    fallback_roots = [Path.home(), Path("/opt"), Path("/usr/local")]
    for root in fallback_roots:
        if not root.exists() or root == start_path or start_path in root.parents:
            continue  # Skip if we already covered it
        try:
            for path in root.glob("**/x3dna*"):
                if path.is_dir() and ((path / "bin" / "find_pair").exists() or (path / "bin" / "fiber").exists()):
                    return path
        except (PermissionError, FileNotFoundError):
            continue

    return None

def generate_dynamic_config():
    # 1. Locate the absolute path of this package deployment
    pkg_root = Path(__file__).resolve().parent
    config_path = pkg_root / "ModCRElib" / "configure" / "config.ini"
    
    if not config_path.exists():
        print(f"Error: Base configuration file not found at {config_path}")
        return

    # 2. Extract active Conda environment values
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        print("Warning: No active Conda environment detected ($CONDA_PREFIX is empty).")
        print("Falling back to standard system-wide PATH evaluation.\n")
    
    # 3. Read the existing config template
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(config_path)
    
    print(f"Automating configuration updates for: {config_path}")
    print(f"Package Root detected at: {pkg_root}\n")

    # 4. Automate static file asset locations
    if "Paths" not in config:
        config["Paths"] = {}

    config["Paths"]["scripts_path"] = str(pkg_root)
    config["Paths"]["src_path"] = str(pkg_root / "ModCRElib")
    config["Paths"]["pdb_dir"] = str(pkg_root / "pdb")
    config["Paths"]["pbm_dir"] = str(pkg_root / "pbm")
    config["Paths"]["files_path"] = str(pkg_root / "files")
    config["Paths"]["modpy_path"] = str(pkg_root / "modpy")
    config["Paths"]["sbilib_path"] = str(pkg_root / "SBILib")
    
    # Sub-files mapping
    config["Paths"]["TF_GOMF"] = str(pkg_root / "files" / "TF_molecular_function_w.txt")
    config["Paths"]["TF_GOBP"] = str(pkg_root / "files" / "TF_biological_process_w.txt")
    config["Paths"]["nTF_GOMF"] = str(pkg_root / "files" / "NonTF_molecular_function.txt")
    config["Paths"]["nTF_GOBP"] = str(pkg_root / "files" / "NonTF_biological_process.txt")
    config["Paths"]["posKW"] = str(pkg_root / "files" / "Positive_keywords.txt")
    config["Paths"]["negKW"] = str(pkg_root / "files" / "Negative_keywords.txt")
    config["Paths"]["species"] = str(pkg_root / "files" / "speclist.txt")
    
    # Historical paths
    config["Paths"]["jaspar_dir"] = str(pkg_root / "files")
    config["Paths"]["cisbp_dir"] = str(pkg_root / "files")

    # 5. Automatically locate external executables using PATH/Conda environment
    binary_mappings = {
        "python_path": "python",
        "blast_path": "blastp",
        "clustalo_path": "clustalo",
        "dssp_path": "mkdssp",
        "emboss_path": "matcher",
        "hmmer_path": "hmmscan",
        "meme_path": "meme",
        "tmalign_path": "TMalign",
        "weblogo_path": "weblogo",
        "cd-hit": "cd-hit",
        "ghostscript_path": "gs",
        "mmseqs": "mmseqs",
    }

    for config_key, binary_name in binary_mappings.items():
        found_path = shutil.which(binary_name)
        if found_path:
            found_path = str(Path(shutil.which(binary_name)).parent)
        
        if not found_path and binary_name == "mkdssp":
            found_path = shutil.which("dssp")

        
        if found_path:
            config["Paths"][config_key] = found_path
            print(f"✓ Found {binary_name} -> {found_path}")
        else:
            config["Paths"][config_key] = ""
            print(f"✗ Could not automatically locate binary: {binary_name}")

    # =========================================================
    # 6. DYNAMIC AND BOTTOM-UP X3DNA SETUP
    # =========================================================
    x3dna_env = os.environ.get("X3DNA")
    x3dna_detected_path = ""

    if x3dna_env:
        x3dna_detected_path = str(Path(x3dna_env) / "bin")
        print(f"✓ Found X3DNA via environment path -> {x3dna_detected_path}")
    else:
        # Check standard system execution paths
        find_pair_path = shutil.which("find_pair") or shutil.which("fiber")
        if find_pair_path:
            x3dna_detected_path = str(Path(find_pair_path).parent)
            print(f"✓ Found X3DNA via active system PATH -> {x3dna_detected_path}")
        else:
            # Deep proximity file scan strategy
            print("⚠ X3DNA not found in active PATH. Initiating proximity search mapping...")
            deep_search_result = find_x3dna_bottom_up(pkg_root)
            if deep_search_result:
                x3dna_detected_path = str(deep_search_result / "bin")
                print(f"✓ Dynamically discovered local X3DNA setup -> {x3dna_detected_path}")

    # Commit evaluated location back to config object
    if x3dna_detected_path:
        config["Paths"]["x3dna_path"] = x3dna_detected_path
    else:
        config["Paths"]["x3dna_path"] = ""
        print("✗ X3DNA path could not be automatically determined anywhere.")





    # =========================================================
    # ROBUST DYNAMIC MODELLER DETECTOR
    # =========================================================
    modeller_bin_path = ""
    root_path = None

    # =========================================================
    # 1. FIXED CUSTOM ENVIRONMENT VARIABLES (Resilient Check)
    # =========================================================
    custom_env_root = os.environ.get("MODELLER_ROOT") or os.environ.get("MODELLER_HOME")
    if custom_env_root:
        # Expand tokens like '~' if used in shell exports
        potential_root = Path(os.path.expanduser(custom_env_root)).resolve()
        
        if potential_root.exists():
            # Check standard subdirectories AND the root directory itself case-insensitively
            possible_bins = [
                potential_root / "bin" / "modpy.sh",
                potential_root / "Bin" / "modpy.sh",
                potential_root / "modpy.sh"
            ]
            
            for candidate_bin in possible_bins:
                if candidate_bin.exists():
                    root_path = potential_root
                    modeller_bin_path = str(candidate_bin.parent)
                    print(f"✓ Found Modeller via custom environment variable -> {root_path}")
                    break

    # =========================================================
    # 2. CHECK CONDA ENVIRONMENT
    # =========================================================
    if not modeller_bin_path:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            candidate_bin = Path(conda_prefix) / "bin" / "modpy.sh"
            if candidate_bin.exists():
                root_path = Path(conda_prefix).resolve()
                modeller_bin_path = str(candidate_bin.parent)
                print(f"✓ Found Modeller via active Conda environment -> {root_path}")

    # Alternative Conda / Python path check
    if not modeller_bin_path:
        for path_str in sys.path:
            candidate_lib = Path(path_str) / "modlib"
            if candidate_lib.exists():
                for parent in candidate_lib.parents:
                    candidate_bin = parent / "bin" / "modpy.sh"
                    if candidate_bin.exists():
                        root_path = parent.resolve()
                        modeller_bin_path = str(candidate_bin.parent)
                        print(f"✓ Found Modeller via Python site-packages -> {root_path}")
                        break
            if modeller_bin_path:
                break

    # =========================================================
    # 3. CHECK ACTIVE SYSTEM PATH
    # =========================================================
    if not modeller_bin_path:
        modpy_executable = shutil.which("modpy.sh")
        if modpy_executable:
            resolved_bin = Path(modpy_executable).resolve()
            possible_root = resolved_bin.parent.parent
            if (possible_root / "modlib").exists():
                root_path = possible_root.resolve()
                modeller_bin_path = str(resolved_bin.parent)
                print(f"✓ Found Modeller via active system PATH -> {root_path}")

    # =========================================================
    # 4. SCAN STANDARD SYSTEM + USER INSTALLATION DIRECTORIES
    # =========================================================
    if not modeller_bin_path:
        current_os = platform.system()
        search_roots = []
        home = Path.home()

        if current_os == "Linux":
            # Added common unprivileged user binary folders (~/bin, ~/.local/lib)
            search_roots = [
                home / "bin", 
                home / ".local" / "lib", 
                Path("/usr/lib"), 
                Path("/usr/local/lib"), 
                Path("/opt")
            ]
        elif current_os == "Darwin":  # macOS
            search_roots = [home / "bin", Path("/opt/local/lib"), Path("/usr/local/lib"), Path("/usr/local")]
        elif current_os == "Windows":
            prog_files = os.environ.get("ProgramFiles", "C:\\Program Files")
            prog_files_x86 = os.environ.get("ProgramFiles(x86)", "C:\\Program Files (x86)")
            search_roots = [Path(prog_files), Path(prog_files_x86)]

        # Search for any folder starting with "Modeller" (case-insensitive)
        for root in search_roots:
            if not root.exists():
                continue
            try:
                # Direct check: Is the root directory itself a "modeller" directory?
                folders_to_check = [root] if "modeller" in root.name.lower() else list(root.iterdir())
                
                for item in folders_to_check:
                    if item.is_dir() and "modeller" in item.name.lower():
                        bin_file = "modpy.bat" if current_os == "Windows" else "modpy.sh"
                        
                        # Look inside standard structures inside this found directory
                        possible_bins = [item / "bin" / bin_file, item / bin_file]
                        for candidate_bin in possible_bins:
                            if candidate_bin.exists():
                                root_path = item.resolve()
                                modeller_bin_path = str(candidate_bin.parent)
                                print(f"✓ Found manual Modeller installation -> {root_path}")
                                break
                    if modeller_bin_path:
                        break
            except PermissionError:
                continue 
            
            if modeller_bin_path:
                break




    if not modeller_bin_path:
        print("✗ Modeller installation could not be detected automatically. Please run: export MODELLER_ROOT='path/to/install/modeller10.3'")
       



    # =========================================================
    # 7. INTERACTIVE CLUSTER CONFIGURATION POPULATION
    # =========================================================
    print("\n--- Cluster Configuration Setup ---")
    print("Please provide the following details. Press Enter to keep the default value [in brackets].\n")

    if "Cluster" not in config:
        config["Cluster"] = {}

    cluster_fields = [
        "cluster_name", "cluster_queue", "cluster_submit", "cluster_qstat",
        "max_jobs_in_queue", "min_jobs_in_queue", "command_queue"#,
        #"server_host", "server_user", "server_passwd", "server_directory", "server_python"
    ]

    for field in cluster_fields:
        current_val = config["Cluster"].get(field, "None")
        user_input = input(f"Enter {field} [{current_val}]: ").strip()
        if user_input:
            config["Cluster"][field] = user_input
        else:
            config["Cluster"][field] = current_val

    # 8. Write the populated data back to the file
    with open(config_path, "w") as configfile:
        config.write(configfile)
    print("\n[SUCCESS] config.ini successfully updated with local environment configuration.")

if __name__ == "__main__":
    generate_dynamic_config()
