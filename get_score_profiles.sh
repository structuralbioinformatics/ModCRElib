#!/bin/sh
##$BATCH -l h_vmem=200G
#$BATCH --mem=80G
#purge previous loads
module purge
#start loading the path
#module load modulepath
#load the modules
module load Python/3.8.6-GCCcore-10.2.0
module load ClustalW2/2.1-foss-2020b
module load Clustal-Omega/1.2.4-foss-2020b
#module load dssp/1.0 # /usr/bin/dssp
module load EMBOSS/6.6.0-foss-2020b
module load Modeller
module load MEME/5.1.1-GCCcore-10.2.0-Python-3.8.6
module load BLAST+
module load HMMER/3.3.2-foss-2020b
module load TMalign
#module load hbplus/3.2
module load x3dna
module load CD-HIT/4.8.1-GCC-10.2.0
#module load Rosetta/2016.20
#module load imp/2.21.0-GCC-10.2.0
module load Ghostscript/9.53.3-GCCcore-10.2.0


pdb="/home/pgohl/ModCRE_Package/../TempFileStorage_ModCRE_Package/pdb"
pbm="/home/pgohl/ModCRE_Package/../TempFileStorage_ModCRE_Package/pbm"
dummy="./dummy"

python exe/xprofiler.py --dummy $dummy -i profiling/ASBprofilerinput.txt  -d profiling/profile.fa  --pbm $pbm --pdb $pdb -v --auto --known --refine 2 --plot --html > profiling/profiler.log 2>&1
#python exe/profiler.py --dummy $dummy -i profiling  -d profiling/profile.fa  --pbm $pbm --pdb $pdb -v  > profiling/profiler.log 2>&1

