"""
Convert one mmCIF structure file to PDB format.

Usage:
    python mmcif_to_pdb.py input.cif output.pdb
"""

import sys,os
from Bio.PDB import PDBIO, Superimposer, PDBParser,MMCIFParser,MMCIFIO

cifparser     = MMCIFParser()
x=cifparser.get_structure("structure",sys.argv[1])
io=PDBIO()
io.set_structure(x)
io.save(sys.argv[2])
print("Done")


