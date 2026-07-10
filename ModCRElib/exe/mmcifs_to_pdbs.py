"""
Batch-convert a folder of mmCIF files into PDB files with optional renaming.

Usage:
    python mmcifs_to_pdbs.py INPUT_DIR OUTPUT_DIR [LABEL] [SPECIE]
"""

import sys,os
import configparser

# Get scripts path (i.e. ".") #
exe_path = os.path.abspath(os.path.dirname(__file__))
if os.path.exists(os.path.join(exe_path,"..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..")
elif os.path.exists(os.path.join(exe_path,"..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..")
elif os.path.exists(os.path.join(exe_path,"..","..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..","..")
else:
   scripts_path = os.path.join(exe_path)

config_path  = os.path.join(scripts_path,"ModCRElib","configure")

# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = configparser.ConfigParser()
config_file = os.path.join(config_path, "config.ini")
config.read(config_file)

# Imports my functions #
from SBILib.structure import PDB

input_folder=os.path.abspath(sys.argv[1])
output_folder=os.path.abspath(sys.argv[2])
print("Output folder: %s"%output_folder)

if len(sys.argv)>3: 
   label=sys.argv[3]+"_"
   print("Label: %s"%label)
else:
   label=None
if len(sys.argv)>4: 
   specie=sys.argv[4]
   print("Specie: %s"%specie)
else:
   specie=""




cwd = os.getcwd()
os.chdir(input_folder)
for f in os.listdir(input_folder):
    if f.endswith("cif"):
       print("Read CIF %s"%f)
       pdb=PDB(cif_file=f)
       if label is not None:
        min=10000
        max=0
        chain=""
        for c in pdb.chains:
           if c.chaintype=="P":
              chain=chain+c.chain
              idx=c.protein_idx.split(";")
              if int(idx[0]) < min: min=int(idx[0])
              if int(idx[-1]) > max: max=int(idx[-1])
        frg=specie+":"+str(min)+":"+str(max)+"_alf3_"+chain+"_"
        tf_name=f.split("_")[1]
        tf_base=f.lstrip("fold_").rstrip(".cif").split("_model_")[0]  
        tf_sub="".join([x for x in tf_base.upper().split("_") if x!=tf_name.upper()])
        number=f.rstrip(".cif").split("_")[-1]   
        tf=label+tf_name.upper()+"_"+tf_name.upper()+"_"+tf_sub.upper()+frg+number+".pdb"
       else:
        tf = ".".join([x for x in f.split(".")[0:-1]])+".pdb"
       print("Write PDB %s"%tf)
       new_pdb=os.path.join(output_folder,tf)
       pdb.write(new_pdb,force=True)
os.chdir(cwd)
print("Done")

