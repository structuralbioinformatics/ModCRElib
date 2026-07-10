import os
import shutil
import sys
import tarfile
import urllib.request
import zipfile
from pathlib import Path

# Define paths relative to this script's position
ROOT_DIR = Path(__file__).resolve().parent
TARGET_BASE_DIR = ROOT_DIR / "ExternalPWMs"

# Updated multi-source configurations
DATABASES = {
    "jaspar_pwms": "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_meme.zip",
    "CisBP_pwms": "https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.27.tgz",
    "hocomoco_pwms": "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE/H14CORE_pwm.tar.gz"
}

def download_with_headers(url, output_path):
    """Downloads files safely using a browser user-agent string."""
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    req = urllib.request.Request(url, headers=headers)
    
    with urllib.request.urlopen(req) as response, open(output_path, 'wb') as out_file:
        total_size = int(response.info().get('Content-Length', 0))
        downloaded = 0
        block_size = 65536
        
        while True:
            buffer = response.read(block_size)
            if not buffer:
                break
            downloaded += len(buffer)
            out_file.write(buffer)
            
            if total_size > 0:
                percent = min(100, int(downloaded * 100 / total_size))
                sys.stdout.write(f"\r    Downloading... {percent}%")
                sys.stdout.flush()
        print("\n    Download complete.")

def flatten_directory(target_dir, sub_path_string):
    """Moves nested archive directories up to the root level of target_dir."""
    nested_dir = target_dir / sub_path_string
    if nested_dir.exists() and nested_dir.is_dir():
        print("    Flattening files from nested folder...")
        for item in nested_dir.iterdir():
            shutil.move(str(item), str(target_dir / item.name))
        
        top_nested_folder = target_dir / sub_path_string.split('/')[0]
        if top_nested_folder.exists():
            shutil.rmtree(top_nested_folder)

def split_multi_meme_file(target_dir):
    """Splits a multi-motif MEME file into individual motif .meme files."""
    source_file = target_dir / "motif_databases" / "CIS-BP_2.00" / "Homo_sapiens.meme"
    
    if not source_file.exists():
        print(f"    [Error] Could not locate multi-meme file at: {source_file}")
        return
        
    print("    Parsing and splitting multi-motif MEME file into individual entries...")
    
    with open(source_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    header_lines = []
    motifs_list = []
    current_motif = []
    
    # Separate global metadata headers from individual motif blocks
    for line in lines:
        if line.startswith("MOTIF"):
            if current_motif:
                motifs_list.append(current_motif)
                current_motif = []
            current_motif.append(line)
        elif current_motif:
            current_motif.append(line)
        else:
            header_lines.append(line)
            
    if current_motif:
        motifs_list.append(current_motif)

    # Write each motif out to its own file using its unique motif ID
    for motif_block in motifs_list:
        first_line_parts = motif_block[0].split()
        if len(first_line_parts) > 1:
            motif_id = first_line_parts[1] # e.g., M0123_2.00
            output_file = target_dir / f"{motif_id}.meme"
            
            with open(output_file, 'w', encoding='utf-8') as out:
                out.writelines(header_lines)
                if not header_lines[-1].endswith('\n'):
                    out.write('\n')
                out.writelines(motif_block)

    # Remove the massive parent directory uncompressed by MEME Suite to save space
    temp_unzipped_root = target_dir / "motif_databases"
    if temp_unzipped_root.exists():
        shutil.rmtree(temp_unzipped_root)
        
    print(f"    Successfully isolated {len(motifs_list)} individual motif profiles.")

def download_and_extract():
    print("Initializing ModCRElib Reference Database Sync Setup...")
    
    for folder_name, url in DATABASES.items():
        target_dir = TARGET_BASE_DIR / folder_name
        print(f"\nProcessing Database Layer [{folder_name}]...")
        
        if target_dir.exists() and any(target_dir.iterdir()):
            print(f"  -> Path '{target_dir.relative_to(ROOT_DIR)}' is already populated. Skipping.")
            continue
            
        target_dir.mkdir(parents=True, exist_ok=True)
        file_extension = ".tar.gz" if url.endswith(".tgz") or url.endswith(".tar.gz") else ".zip"
        temp_archive = TARGET_BASE_DIR / f"temp_{folder_name}{file_extension}"
        
        try:
            download_with_headers(url, temp_archive)
            
            with open(temp_archive, 'rb') as f:
                start_bytes = f.read(100)
                if b'<!DOCTYPE html>' in start_bytes or b'<html' in start_bytes:
                    raise ValueError("The host rejected your query or served an unexpected HTML page.")

            print(f"    Unpacking contents to {target_dir.relative_to(ROOT_DIR)}...")
            if file_extension == ".zip":
                with zipfile.ZipFile(temp_archive, 'r') as zip_ref:
                    zip_ref.extractall(target_dir)
            elif file_extension in [".tar.gz", ".tgz"]:
                with tarfile.open(temp_archive, 'r:gz') as tar_ref:
                    tar_ref.extractall(target_dir)
            print("    Archive expansion completed.")
            
            # Post-extraction normalization tasks
            if folder_name == "hocomoco_pwms":
                flatten_directory(target_dir, "H14CORE/pwm")
            elif folder_name == "CisBP_pwms":
                split_multi_meme_file(target_dir)
                
        except Exception as e:
            print(f"    [Processing Error]: {e}")
        finally:
            if temp_archive.exists():
                temp_archive.unlink()

    print("\nAll database layers have been normalized, isolated, and formatted successfully.")

if __name__ == "__main__":
    download_and_extract()

