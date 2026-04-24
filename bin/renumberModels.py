import sys,os
import configparser

# Read configuration file #
config = configparser.ConfigParser()
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


config_file  = os.path.join(scripts_path,"ModCRElib","configure/config.ini")
config.read(config_file)

# Get sbilib path #
sbilib_path = config.get("Paths", "sbilib_path")


sys.path.append(sbilib_path)
from SBILib.structure import PDB
sourcedir = sys.argv[1]
destdir = sys.argv[2]

models = [f for f in os.listdir(sourcedir) if f.endswith(".pdb")]



for x in models:
    fulldir = sourcedir +"/"+ x
    renamed = destdir +"/"+  x
    protein = PDB(fulldir)
    for chainID in protein.chain_identifiers:
        protein.get_chain_by_id(chainID).renumber_residues(1)
        protein.get_chain_by_id(chainID).renumber_atoms(1)
    protein.write(renamed)
    
    
