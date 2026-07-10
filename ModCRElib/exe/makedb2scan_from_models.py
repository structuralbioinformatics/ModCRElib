"""
Build a PWM scanning database from structural models and PWM files.

The script extracts protein sequences from a directory of PDB models, prepares
the FASTA files expected by the ModCRE scanning pipeline, concatenates MEME PWM
files into a single database file, and then delegates the final database
assembly to ``ModCRElib.web.md2p.get_userdb``.
"""

import os, sys, re
import configparser
import json
import numpy as np
import optparse
import shutil
import subprocess
from time import time
import pickle
import pickle

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

# Imports jbonet's module #
from SBILib.structure import PDB
from SBILib.data      import aminoacids1to3

# Imports my functions #
from ModCRElib.profile import fimo, scorer, tomtom
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.structure.threading import threader
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.protein import dssp, model_protein
from ModCRElib.sequence import nr
from ModCRElib.msa import  pwm_pbm as PWM
from ModCRElib.sequence import homologs as HOMO
from ModCRElib.web import md2p as MD2P

# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Directories #
pbm_dir = config.get("Paths", "pbm_dir")
pdb_dir = config.get("Paths", "pdb_dir")


#---------------#
# Functions     #
#---------------#




def parse_options():
    """
    Parse command-line options for database generation from structural models.

    How to run:
        python makedb2scan_from_models.py -i PWMS_DIR -p MODELS_DIR
            [-o OUTPUT_DIR -v]

    Example:
        python makedb2scan_from_models.py -i pwm_models -p pdb_models
            -o scan_db -v

    The parser configures:
        - Input PWM directory and directory with structural models.
        - Optional output directory for generated scan database artifacts.
        - Verbose logging mode.

    Args:
        None.

    Returns:
        optparse.Values: Parsed options with the input PWM directory, the input
        models directory, the optional output directory, and the verbosity
        flag.

    Raises:
        SystemExit: Triggered by ``OptionParser`` when required arguments are
        missing or invalid.
    """

    parser = optparse.OptionParser("Usage: makedb2scan_from_models.py -i PWMS_DIR -p MODELS_DIR [-o OUTPUT_DIR -v]")

    parser.add_option("-i", "--pwm-dir", default=None, action="store", type="string", dest="pwms_dir", help="PWMs directory ", metavar="PWMS_DIR")
    parser.add_option("-p", "--pdb", default=None, action="store", type="string", dest="models_dir", help="Folder with protein 3D models ", metavar="MODELS_DIR")
    parser.add_option("-o", "--output-dir", default=None, action="store", type="string", dest="output_dir", help="Output directory (default = PWMS_DIR)", metavar="OUTPUT_DIR")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    
    (options, args) = parser.parse_args()

    if  options.models_dir is None or options.pwms_dir is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options


def main():
    """
    Generate a scan-ready PWM database from model structures and motif files.

    Workflow:
        1. Parse command-line arguments and resolve the working directories.
        2. Read each input PDB model and extract concatenated protein-chain
           sequences into a FASTA file when it does not already exist.
        3. Create per-protein FASTA files plus sequence lookup tables in the
           output ``sequences`` directory.
        4. Merge all MEME-formatted PWM files into ``pwm_database.txt``.
        5. Call ``MD2P.get_userdb`` to build the final PWM database,
           association file, and homologs directory used downstream.

    Returns:
        None. The function writes its results to the selected output directory
        and reports the generated paths to standard output.
    """

    # Arguments & Options #
    try:
      options = parse_options()
    except Exception as e:
      print(e)
      exit()

    uid_set      = set() #PWMs manually given, should be none
    pdb_dir      = os.path.abspath(options.models_dir)
    fasta_file   = os.path.join(os.path.dirname(pdb_dir),"sequences.fasta")
    database_dir = os.path.abspath(options.pwms_dir)

    pdb_set=set()
    for pdb_file in os.listdir(pdb_dir):
        if pdb_file.endswith("pdb"):
           pdb_set.add(os.path.join(pdb_dir,pdb_file))

    if not os.path.exists(fasta_file):
       sequence_dict={}
       for pdb_file in pdb_set:
        pdb_obj=PDB(pdb_file)
        name=pdb_obj.id
        sequence=""
        for chain_id in pdb_obj.chain_identifiers:
            chain=pdb_obj.get_chain_by_id(chain_id)
            if chain.chaintype == "P":
               sequence+=chain.protein_sequence
        sequence_dict.setdefault(name,sequence)
       fseq=open(fasta_file,"w")
       for name in sequence_dict:
        fseq.write(">%s\n%s\n"%(name,sequence_dict[name]))
       fseq.close()

    if options.output_dir is not None:
       output_dir= os.path.abspath(options.output_dir)
    else:
       output_dir= database_dir
    if output_dir != database_dir:
       if not os.path.exists(output_dir): 
          shutil.copytree(database_dir,output_dir)
       else:
          shutil.rmtree(output_dir)
          shutil.copytree(database_dir,output_dir)
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    folder       = os.path.abspath(os.path.join(output_dir))
    database_file= os.path.join(folder,"pwm_database.txt")
    if not os.path.exists(os.path.join(folder,"sequences")):os.makedirs(os.path.join(folder,"sequences"))
     
    if not os.path.exists(database_dir): 
              print("PWMS folder  not found")
              exit()
     

    sequence_db = open(os.path.join(folder, "sequences", "sequences.fasta"), "w")
    model_db    = open(os.path.join(folder, "sequences", "models.fasta"), "w")
    mdl2fasta={}
    for name,sequence in functions.parse_fasta_file(fasta_file):
        m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)",name)
        if m:
                        uid=m.group(2)
                        #uid_set.add(uid)
                        sequence_db.write(">%s\n%s\n"%(name,sequence))
                        model_db.write(">%s\n%s\n"%(m.group(2),sequence))
                        mdl2fasta.setdefault(m.group(2),name)
                        mdl2fasta.setdefault(name,m.group(2))
                        fasta = open(os.path.join(folder, "sequences", uid + ".fasta"), "w")
                        fasta.write(">%s\n%s\n"%(m.group(2),sequence))
                        fasta.close()
        else:
                        uid=name.split()[0].rstrip().replace("\r", "").replace("\n", "").replace(" ", "").replace("|","_")
                        #uid_set.add(uid)
                        sequence_db.write(">%s\n%s\n"%(name,sequence))
                        model_db.write(">%s\n%s\n"%(uid,sequence))
                        mdl2fasta.setdefault(uid,name)
                        mdl2fasta.setdefault(name,uid)
                        fasta = open(os.path.join(folder, "sequences", uid + ".fasta"), "w")
                        fasta.write(">%s\n%s\n"%(uid,sequence))
                        fasta.close()
                        print(("SKIP SEQUENCE %s :  Wrong format sp|ACCESSION|GENE_SPECIE"%name))
    sequence_db.close()
    model_db.close()
    if not  os.path.exists(database_file):
       for pwm_file in os.listdir(database_dir):
           if pwm_file.endswith("meme"):
              #print("Add ",pwm_file)
              os.system("cat %s >> %s"%(os.path.join(database_dir,pwm_file),database_file))
    #print("UID_SET",uid_set)
    #print("mdl2fasta",mdl2fasta)
    try:
              database_pwm, association_pwm, homologs_dir = MD2P.get_userdb(folder,fasta_file,mdl2fasta,uid_set,database_file) 
              print("DATABASE FILE:\t\t%s"%(database_pwm))
              print("DATABASE DIRECTORY:\t%s"%(os.path.dirname(database_pwm)))
              print("ASSOCIATION FILE:\t%s"%(association_pwm))
              print("HOMOLOGS DIRECTORY:\t%s"%(homologs_dir))
              print("\nDone")
    except Exception as e:
              print("Fail when constructing database of PWMs: check the format of the protein sequences")
              print(e)


#---------------#
# Main          #
#---------------#

if __name__ == "__main__":
    main()
