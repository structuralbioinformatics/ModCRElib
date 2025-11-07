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
pdb="./pdb"


#python exe/model_multiple_proteins.py -v -d ./dummy -i NonPBMTFs.fa -o ./models --pdb=${pdb} --best --renumerate  --n-temp=1 --n-total=1 --n-model=2

python exe/model_multiple_proteins.py -v -d ./dummy -i Threads_list.txt -t -o ./scanned_models --pdb=${pdb} --best --renumerate  --n-temp=1 --n-total=1 --n-model=2
