#!/bin/sh

echo "[INFO] Setting up environment..."

if command -v module >/dev/null 2>&1; then
    echo "[INFO] Detected module system → loading HPC modules"

    module purge

    module load Python/3.8.6-GCCcore-10.2.0
    module load ClustalW2/2.1-foss-2020b
    module load Clustal-Omega/1.2.4-foss-2020b
    module load EMBOSS/6.6.0-foss-2020b
    module load Modeller
    module load MEME/5.1.1-GCCcore-10.2.0-Python-3.8.6
    module load BLAST+
    module load HMMER/3.3.2-foss-2020b
    module load TMalign
    module load x3dna/2.3
    module load CD-HIT/4.8.1-GCC-10.2.0
    module load Ghostscript/9.53.3-GCCcore-10.2.0

else
    echo "[INFO] No module system detected → assuming conda/local environment"

    # Optional: activate conda if available
    if [ -n "$CONDA_PREFIX" ]; then
        echo "[INFO] Using active conda environment: $CONDA_PREFIX"
    else
        echo "[WARNING] No conda environment detected"
    fi
fi




