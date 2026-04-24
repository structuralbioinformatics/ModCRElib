#!/bin/sh
source env_setup.sh

pdb="./pdb"

#python exe/model_multiple_proteins.py -v -d ./dummy -i NonPBMTFs.fa -o ./models --pdb=${pdb} --best --renumerate  --n-temp=1 --n-total=1 --n-model=2

python exe/model_multiple_proteins.py -v -d ./dummy -i Threads_list.txt -t -o ./scanned_models --pdb=${pdb} --best --renumerate  --n-temp=1 --n-total=1 --n-model=2
