#!/bin/sh
source env_setup.sh

pdb="./pdb"
pbm="./pbm"
dummy="./dummy"

#/soft/system/software/x3dna/2.3/bin/fiber  -seq=GAGCTGCCCACCCCGCTACAGCAGCTGGCTGTGCAAGGTGGCTGGACCACA -b BinaryInteractions/dnasequence.pdb
/soft/system/software/x3dna/2.3/bin/fiber -seq=aaaaaaaaaaaaatttaggaggcccaaaaaaaaaaaaa -b BinaryInteractions/dnasequence.pdb

python exe/complexbuilder.py -d BinaryInteractions -o Complex > Complex/complex.log 2>&1 

