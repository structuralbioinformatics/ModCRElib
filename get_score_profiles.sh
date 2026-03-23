#!/bin/sh
source env_setup.sh

pdb="./pdb"
pbm="./pbm"
dummy="./dummy"

python exe/xprofiler.py --dummy $dummy -i profiling/ASBprofilerinput.txt  -d profiling/profile.fa  --pbm $pbm --pdb $pdb -v --auto --known --refine 2 --plot --html > profiling/profiler.log 2>&1
#python exe/profiler.py --dummy $dummy -i profiling  -d profiling/profile.fa  --pbm $pbm --pdb $pdb -v  > profiling/profiler.log 2>&1

