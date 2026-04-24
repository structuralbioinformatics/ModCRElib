import os, sys, re
import configparser
import optparse
import shutil
import subprocess
import socket
import time

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
from ModCRElib.beans import functions

# Get python path #
python = os.path.join(config.get("Paths", "python_path"), "python")


# Imports jbonet's module #
from SBILib.structure import PDB
from SBILib.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBILib.structure.residue import ResidueOfNucleotide
from SBILib.structure.atom import AtomOfNucleotide


#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for renaming model files for complex building.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options controlling input/output folders
        and naming conventions.
    """

    parser = optparse.OptionParser("python models4complex.py -i INPUT_DIRECTORY -o OUTPUT_DIRECTORY [--code PDB_CODE --name NAME --start START --dummy DUMMY_DIR]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_folder", help="Input directory of models in PDB format (suffixed as pdb)", metavar="{directory}")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output directory of models in PDB format with the name modified to be used by complexbuilder", metavar="{directory}")
    parser.add_option("--code", action="store",default="1mdl", type="string", dest="code", help="PDB code of the closest template of the model or any other alphanumeric code of 4 digits (Chains/s are taken from the structure, default is 1mdl)", metavar="{string}")  
    parser.add_option("--name", action="store",default="TF", type="string", dest="name", help="Protein root-name of the proteins of the directory (a list is produced by numbering the structures of the folder, default is TF)", metavar="{string}")  
    parser.add_option("--start", default=1,action="store", type="int", dest="start", help="Start position in the DNA (last position is calculated with the DNA sequence of the model (default is 1)", metavar="{integer}")  
    parser.add_option("-j", "--join", default=False, action="store_true", dest="join", help="Joint all chains in each PDB to a single one (default = False)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
 


    (options, args) = parser.parse_args()

    if options.input_folder is None or options.output_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


def features(pdb_file):
    """
        Extract chain and boundary metadata from one PDB model.
    
        Args:
            pdb_file (str): Path to a PDB file.
    
        Returns:
            tuple: `(chain, first, last, dna_size)` where `chain` is the protein
            chain-id string, `first/last` are residue identifiers, and `dna_size`
            is the maximum nucleotide-chain length.
    
    Args:
        pdb_file (Any): Path to the input/output file."""
    pdb_obj=PDB(pdb_file)
    chain=""
    dna_size=0
    for chain_id in pdb_obj.chain_identifiers:
        chain_obj=pdb_obj.get_chain_by_id(chain_id)
        if chain_obj.chaintype == "N": 
           dna_size = max(dna_size,len(list(chain_obj.nucleotide_sequence())))
        if chain_obj.chaintype == "P":
           chain+=chain_id
           aa=chain_obj.first_aminoacid
           first=aa.identifier.replace(" ","")
           aa=chain_obj.last_aminoacid
           last=aa.identifier.replace(" ","")

    return chain,first,last,dna_size

def main():
    """
    Run the model-renaming workflow used by `complexbuilder.py`.

    Workflow:
        1. Parse runtime options and collect input PDB models.
        2. Extract chain/boundary metadata from each model.
        3. Rewrite model names and optional chain-joined variants.
        4. Save a table mapping original files to generated model names.

    Returns:
        None. Renamed model files are written to the output directory.
    """
    # Arguments & Options #
    options   = parse_options()
    folder    = options.input_folder
    outdir    = os.path.abspath(options.output_dir)
    if not os.path.exists(outdir): os.makedirs(outdir)
    name      = options.name
    code      = options.code
    start     = int(options.start)
    verbose   = options.verbose
    join      = options.join

    models    = []
    for pdb_file in os.listdir(folder):
        if not pdb_file.endswith("pdb"): continue
        models.append(os.path.join(folder,pdb_file))

    correspondence=dict()
    for i in range(len(models)):
        pdb_file     = models[i]
        protein_name = name+"_"+str(i+1)
        chain,first,last,dna_size = features(pdb_file)
        if len(chain)>1 and join:
           chain_base = chain[0]
           new_name = protein_name+"."+code+"_"+chain_base+"."+str(start)+"-"+str(start+dna_size-1)+":"+first+":"+last+"_"+chain_base+".pdb"
           new_file = os.path.join(outdir,new_name)
           pdb_obj=PDB(pdb_file)
           new_pdb=PDB()
           m=0
           n=0
           for chain_id in pdb_obj.chain_identifiers:
               chain_obj=pdb_obj.get_chain_by_id(chain_id)
               if chain_obj.chaintype == "N": 
                  new_chain=ChainOfNucleotide(new_name,chain_obj.chain)
                  for residue in chain_obj.nucleotides:
                      m = m + 1
                      residue.number =str(m)
                      residue.version=residue.version
                      new_chain.add_residue(residue)
                  new_pdb.add_chain(new_chain)
           new_chain=pdb_obj.get_chain_by_id(chain_base)
           for chain_id in list(chain):
               if chain_id==chain_base:continue
               chain_obj=pdb_obj.get_chain_by_id(chain_id)
               if not chain_obj.chaintype == "P":continue
               new_chain2=ChainOfProtein(new_name,chain_base)
               for residue in chain_obj.aminoacids:
                      n = n +1
                      residue.number =str(n)
                      residue.version=residue.version
                      new_chain2.add_residue(residue)
               new_chain.fuse(new_chain2)
           new_pdb.add_chain(new_chain)
           if verbose: print("Rewrite file %s to %s with single chain %s for proteins"%(os.path.basename(pdb_file),os.path.basename(new_file),chain_base))
           new_pdb.clean()
           new_pdb.write(new_file,force=True)    
        else:
           new_name = protein_name+"."+code+"_"+chain+"."+str(start)+"-"+str(start+dna_size-1)+":"+first+":"+last+"_["+chain+"].pdb"
           for chain_base in list(chain):
               new_name_chain = protein_name+"."+code+"_"+chain_base+"."+str(start)+"-"+str(start+dna_size-1)+":"+first+":"+last+"_"+chain_base+".pdb"
               new_file = os.path.join(outdir,new_name_chain)
               pdb_obj=PDB(pdb_file)
               new_pdb=PDB()
               m=0
               for chain_id in pdb_obj.chain_identifiers:
                   chain_obj=pdb_obj.get_chain_by_id(chain_id)
                   if chain_obj.chaintype == "N": 
                      new_chain=ChainOfNucleotide(new_name,chain_obj.chain)
                      for residue in chain_obj.nucleotides:
                          m = m + 1
                          residue.number =str(m)
                          residue.version=residue.version
                          new_chain.add_residue(residue)
                      new_pdb.add_chain(new_chain)
               for chain_id in pdb_obj.chain_identifiers:
                   if chain_obj.chaintype == "P" and chain_id==chain_base:
                      new_pdb.add_chain(chain_obj)
               if verbose: print("Write file %s for specific chain %s"%(os.path.basename(new_file),chain_base))
               new_pdb.clean()
               new_pdb.write(new_file,force=True)
        correspondence.setdefault(os.path.basename(pdb_file),new_name)

    fo=open(os.path.join(outdir,"table_of_models.txt"),"w")
    for i in range(len(models)):
        pdb_file = os.path.basename(models[i])
        fo.write("%d\t%s\t%s\n"%(i+1,pdb_file,correspondence[pdb_file]))
    fo.close()


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()


