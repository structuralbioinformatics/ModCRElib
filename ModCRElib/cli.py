import subprocess
import sys
import importlib
import shutil
from pathlib import Path
import os
import configparser

PACKAGE = Path(__file__).resolve().parent
EXE = PACKAGE / "exe"

COMMANDS = {
    "doctor": None,
    "setup": None,  
    "model": "model_multiple_proteins.py",
    "renumber": "renumberModels.py",
    "pwm": "pwm_pbm.py",
    "thread": "get_best_binding.py",
    "score": "scorer.py",
    "profile": "xprofiler.py",
    "get_json": "get_json.py",
    "aggregate": "aggregate_pwms.py",
    "prep_scan": "make_scan_ready.py",
    "scan": "scan.py",
    "prep_build": "rename_complex_input.py",
    "build_complex": "RunComplexBuilder.py",
}

REQUIRED_MODULES = [
    "numpy",
    "pandas",
    "scipy",
    "Bio",       # Biopython
    "SBILib",    # Make sure this is installed in site-packages or in your path
    "ModCRElib",
    "configparser",
]

REQUIRED_EXECUTABLES = [
    "fiber",     # x3dna provides 'fiber'
    "blastp",    # provided by blast
    "hmmscan",   # provided by hmmer
    "mkdssp",    # provided by dssp 
    "TMalign",   # provided by tmalign
    "matcher",   # provided by Modeller
    "mmseqs",    # rovided by mmseq2
]

REQUIRED_DIRECTORIES = [
    "pdb",
    "pbm",
    "PWMdatabase",
    "ExternalPWMs/CisBP_pwms",
    "ExternalPWMs/hocomoco_pwms",
    "ExternalPWMs/jaspar_pwms",
]

def run_script(script_name):
    """Run legacy executable scripts, forwarding arguments."""
    script = EXE / script_name
    if not script.exists():
        print(f"Error: cannot find {script}")
        sys.exit(1)

    cmd = [sys.executable, str(script), *sys.argv[2:]]
    subprocess.run(cmd, check=True)

def check_module(name):
    try:
        importlib.import_module(name)
        print(f"✓ {name}")
    except Exception as e:
        print(f"✗ {name}\n    {e}")

def check_program(name):
    path = shutil.which(name)
    if path:
        print(f"✓ {name} ({path})")
        return

    if name == "mkdssp" and shutil.which("dssp"):
        print(f"✓ mkdssp (found as dssp)")
        return

    if name == "fiber" and shutil.which("x3dna_fiber"):
        print(f"✓ fiber (found as x3dna_fiber)")
        return

    if name == "fiber":
        cli_dir = Path(__file__).resolve().parent
        config_path = cli_dir / "configure" / "config.ini"
        
        if config_path.exists():
            try:
                config = configparser.ConfigParser(allow_no_value=True)
                config.read(config_path)
                x3dna_path_raw = config.get("Paths", "x3dna_path", fallback="").strip()
                
                if x3dna_path_raw:
                    base_dir = Path(x3dna_path_raw)
                    bin_dir = base_dir if base_dir.name == "bin" else base_dir / "bin"
                    local_fiber = bin_dir / "fiber"
                    local_x3dna_fiber = bin_dir / "x3dna_fiber"
                    chosen_binary = local_fiber if local_fiber.exists() else local_x3dna_fiber

                    if chosen_binary.exists():
                        try:
                            result = subprocess.run(
                                [str(chosen_binary), "-h"], 
                                stdout=subprocess.DEVNULL, 
                                stderr=subprocess.DEVNULL,
                                timeout=2
                            )
                            print(f"✓ fiber (Validated execution via config.ini -> {chosen_binary})")
                            return
                        except (OSError, subprocess.TimeoutExpired):
                            pass 
            except Exception:
                pass

    print(f"✗ {name}")
    if name == "matcher":
        print("\n    [!] Modeller Installation Notice:")
        print("        Modeller needs to be installed with:")
        print("        wget https://salilab.org/modeller/10.3/modeller-10.3.tar.gz")
        print("        following the instructions.")
        print("        Additionally, export MODELLER_ROOT='/path/to/install' needs to be called")
        print("        before running 'modcrelib setup'.\n")
    elif name == "fiber":
        print("\n    [!] X3DNA Installation Notice:")
        print("        The 'fiber' utility belongs to the X3DNA package suite.")
        print("        Please download x3dna 2.4 from the official instructions page:")
        print("        -> http://forum.x3dna.org/site-announcements/download-instructions/")
        print("        following the instructions.\n")

def check_directory(base_package_path, relative_path_str):
    """Resolves the directory relative to the package root and checks its existence."""
    absolute_dir_path = (base_package_path / relative_path_str).resolve()
    if absolute_dir_path.exists() and absolute_dir_path.is_dir():
        print(f"✓ {relative_path_str} ({absolute_dir_path})")
    else:
        print(f"✗ {relative_path_str} (Expected at: {absolute_dir_path})")

def doctor():
    print("\nPython\n------")
    print(f"✓ Python {sys.version.split()[0]}")

    print("\nPython modules\n--------------")
    for module in REQUIRED_MODULES:
        check_module(module)

    print("\nExternal programs\n-----------------")
    for executable in REQUIRED_EXECUTABLES:
        check_program(executable)

    print("\nDirectories\n-----------")
    package_root = Path(__file__).resolve().parent.parent
    for directory in REQUIRED_DIRECTORIES:
        check_directory(package_root, directory)

def main():
    if len(sys.argv) < 2:
        print("""
ModCRElib

Usage:
    modcrelib <command> [options]

Commands:
    doctor            Diagnose missing dependencies
    setup             Configure paths, external files, and structural databases
    model             Build protein models
    pwm               Generate PWMs from a TF model
    thread            Generate thread files for a TF-DNA interaction
    score             Calculate the statistical potential score for TF-DNA interaction
    renumber          Renumber PDB files
    profile           Generate a scoring profile plot for a TF along a DNA sequence
    get_json          Will generate input (json files) for subsequent pwm aggregation steps
    aggregate         Aggregate PWM clusters for a TF
    prep_scan         Prepare pwm database file for scan
    scan              Scan dna sequence for TF binding sites
    prep_build        Prepare Binary interaction file for builder
    build_complex     Build a TF-DNA and TF-protein complex
""")
        sys.exit(0)

    command = sys.argv[1]

    if command == "doctor":
        doctor()
        return

    # =========================================================
    # CONSOLIDATED SETUP SEQUENCE PIPELINE
    # =========================================================
    if command == "setup":
        cli_dir = Path(__file__).resolve().parent
        parent_dir = cli_dir.parent
        
        print("==================================================")
        print(" Step 1: Initializing Environment & Configuration ")
        print("==================================================")
        init_script = parent_dir / "initialize_config.py"
        if init_script.exists():
            try:
                subprocess.run([sys.executable, str(init_script)], check=True)
            except subprocess.CalledProcessError as e:
                print(f"\n✗ Error executing initialization configuration setup: {e}")
                sys.exit(1)
        else:
            print(f"✗ Error: Cannot find 'initialize_config.py' at execution target: {init_script}")
            sys.exit(1)

        print("\n=============================================")
        print(" Step 2: Downloading External PWM Databases  ")
        print("=============================================")
        # Check if the external folders already exist and contain files
        ext_pwm_root = parent_dir / "ExternalPWMs"
        pwm_dirs = ["CisBP_pwms", "hocomoco_pwms", "jaspar_pwms"]
        
        # Check if any data exists in those folders
        pwms_exist = False
        if ext_pwm_root.exists():
            try:
                for sub in pwm_dirs:
                    sub_path = ext_pwm_root / sub
                    if sub_path.exists() and any(sub_path.iterdir()):
                        pwms_exist = True
                        break
            except Exception:
                pass

        if pwms_exist:
            print("→ External PWM databases already present inside ExternalPWMs/. Skipping download.")
        else:
            parent_dir_str = os.path.abspath(str(parent_dir))
            if parent_dir_str not in sys.path:
                sys.path.append(parent_dir_str)
                
            try:
                import download_pwms
                download_pwms.download_and_extract()
            except Exception as e:
                print(f"✗ Notice: Failure parsing standard external PWM downloads: {e}")

        print("\n=============================================")
        print(" Step 3: Fetching Structural Data Frameworks ")
        print("=============================================")
        # Verify if pdb and pbm already exist and are populated
        pdb_path = parent_dir / "pdb"
        pbm_path = parent_dir / "pbm"
        
        large_data_exists = False
        try:
            if pdb_path.exists() and any(pdb_path.iterdir()) and pbm_path.exists() and any(pbm_path.iterdir()):
                large_data_exists = True
        except Exception:
            pass

        if large_data_exists:
            print("→ Large structural directories ('pdb' and 'pbm') detected and populated. Skipping download.")
        else:
            cli_dir_str = os.path.abspath(str(cli_dir))
            if cli_dir_str not in sys.path:
                sys.path.insert(0, cli_dir_str)
                
            try:
                import download_large_data
                download_large_data.extract_and_move_assets(__file__)
            except ImportError as e:
                print(f"✗ Error: Could not locate 'download_large_data.py' adjacent to execution thread.")
                print(f"  Debug Details: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"✗ Structural frameworks deployment failed: {e}")
                sys.exit(1)

        print("\n[SUCCESS] ModCRElib environment workspace fully customized and online.")
        return

    script = COMMANDS.get(command)
    if script is None:
        print(f"Unknown command: {command}")
        sys.exit(1)

    run_script(script)

if __name__ == "__main__":
    main()
