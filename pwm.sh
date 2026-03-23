#!/bin/sh
source env_setup.sh

pdb="./pdb"
pbm="./pbm"


python exe/pwm_pbm.py -i remodels -o pwm --pdb=${pdb} --pbm=${pbm}  --refine 2  --dummy=./dummy --auto --known


