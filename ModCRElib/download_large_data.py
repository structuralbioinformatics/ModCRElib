import os
import sys
import tarfile
import urllib.request
import ssl  # 1. Import the ssl module
import shutil
from pathlib import Path

def download_with_progress(url, output_path):
    """Downloads a file from a URL with an interactive terminal progress bar."""
    print(f"Downloading: {url}")
    try:
        req = urllib.request.Request(
            url, 
            headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
        )
        
        # 2. Create a bypass context for local issuer certification validation slips
        context = ssl._create_unverified_context()
        
        # 3. Pass the unverified context directly to urlopen
        with urllib.request.urlopen(req, context=context) as response:
            total_size = int(response.getheader('Content-Length', 0))
            block_size = 1024 * 1024  # 1 MB chunks
            downloaded = 0

            with open(output_path, 'wb') as f:
                while True:
                    buffer = response.read(block_size)
                    if not buffer:
                        break
                    f.write(buffer)
                    downloaded += len(buffer)
                    
                    if total_size > 0:
                        percent = int(downloaded * 100 / total_size)
                        sys.stdout.write(f"\r    [{'=' * (percent // 4)}{' ' * (25 - percent // 4)}] {percent}% ({downloaded / (1024*1024):.1f} MB / {total_size / (1024*1024):.1f} MB)")
                    else:
                        sys.stdout.write(f"\r    Downloaded: {downloaded / (1024*1024):.1f} MB")
                    sys.stdout.flush()
        print("\n    ✓ Download complete.")
    except Exception as e:
        print(f"\n    ✗ Failed to download {url}: {e}")
        raise

def extract_and_move_assets(cli_file_path):
    """Handles downloading, extracting, and rewriting the structural datasets."""
    cli_dir = Path(cli_file_path).resolve().parent
    target_destination = cli_dir.parent  # The folder ABOVE cli.py (../)

    urls = {
        "pdb": "https://sbi.upf.edu/modcre/views/images/pdb.tgz",
        "pbm": "https://sbi.upf.edu/modcre/views/images/pbm.tgz"
    }

    target_destination.mkdir(parents=True, exist_ok=True)

    for name, url in urls.items():
        archive_path = target_destination / f"{name}.tgz"
        
        try:
            download_with_progress(url, archive_path)
        except Exception:
            if archive_path.exists():
                archive_path.unlink()
            continue

        print(f"Extracting {name}.tgz archive package directly to {target_destination}...")
        try:
            with tarfile.open(archive_path, "r:gz") as tar:
                tar.extractall(path=target_destination)
            print(f"    ✓ Extracted {name} structure matrix successfully.")
        except Exception as e:
            print(f"    ✗ Decompression framework failure for {name}.tgz: {e}")
        finally:
            if archive_path.exists():
                archive_path.unlink()

    # =========================================================
    # 3. INTER-REPOSITORY MERGE & SWAP LOGIC
    # =========================================================
    source_pwms_dir = (cli_dir.parent / "PWMdatabase" / "pwms").resolve()
    target_pbm_dir = (target_destination / "pbm").resolve()
    destination_pwms_slot = target_pbm_dir / "pwms"

    print("\nProcessing PWM Database workspace consolidation...")
    
    if not source_pwms_dir.exists():
        print(f"⚠ Warning: Source cluster directory not found at: {source_pwms_dir}")
        print("  Skipping step: local folder swap sequence omitted.")
        return

    if not target_pbm_dir.exists():
        print(f"✗ Error: Target structure path missing at {target_pbm_dir}. Download may have failed.")
        return

    if destination_pwms_slot.exists():
        print(f"Removing pre-packaged reference tree framework at: {destination_pwms_slot}")
        if destination_pwms_slot.is_symlink() or destination_pwms_slot.is_file():
            destination_pwms_slot.unlink()
        else:
            shutil.rmtree(destination_pwms_slot)

    print(f"Injecting localized package data array to: {destination_pwms_slot}")
    shutil.copytree(source_pwms_dir, destination_pwms_slot)
    print("✓ Operational structure update step finalized successfully.\n")
