#!/bin/sh
source env_setup.sh

pdb="./pdb"
pbm="./pbm"
dummy="./dummy"
delta=1


#python get_best_bindings_threads.py  dna_table.txt

#python exe/get_best_binding.py --standard --dummy $dummy -i remodels/MODEL_sp_P35869_AHR_HUMAN:34:272_5nj8_A_1.pdb  -o threads --pdb $pdb --pbm $pbm --seq Ref1.fa --dna agctggc --delta 1 --pwm pwm/MODEL_sp_P35869_AHR_HUMAN:34:272_5nj8_A_1.meme --verbose --auto --known 
python exe/get_best_binding.py --standard --dummy $dummy -i remodels/MODEL_sp_P35869_AHR_HUMAN:34:272_5nj8_A_2.pdb  -o threads --pdb $pdb --pbm $pbm --seq Ref1.fa --dna agctggc --delta 1 --pwm pwm/MODEL_sp_P35869_AHR_HUMAN:34:272_5nj8_A_2.meme --verbose --auto --known
