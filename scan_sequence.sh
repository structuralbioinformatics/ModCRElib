#!/bin/sh
source env_setup.sh

pdb="./pdb"
pbm="./pbm"
dumm="./dummy"

python exe/scan.py --dummy=$dummy -i profiling/profile.fa -o scanning --pbm=$pbm --pdb=$pdb -s 9606 -v --db ./pwm/scanning_database.txt -c --ft 0.05 > scanning/Scan_output.log 2>&1
