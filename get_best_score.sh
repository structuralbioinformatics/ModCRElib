#!/bin/sh
source env_setup.sh

pdb="./pdb"
pbm="./pbm"
dummy="./dummy"

#python get_best_score.py 

python exe/scorer.py --threading --dummy $dummy -i threads/MODEL_SP_P35869_AHR_HUMAN:34:272_5NJ8_A_2_A_CAGCTGGCTGTG.txt  --norm -o scores/MODEL_SP_P35869_AHR_HUMAN:34:272_5NJ8_A_2_A_CAGCTGGCTGTG --pdb $pdb --pbm $pbm  --verbose  --auto --known






