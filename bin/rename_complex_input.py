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

file=sys.argv[1]

protein = PDB(file)

for x in protein.proteins:
    begin = x.first_aminoacid.number
    end = x.last_aminoacid.number

file = file.replace(".pdb",f":{begin}:{end}_TF.pdb")
protein.write(file)

