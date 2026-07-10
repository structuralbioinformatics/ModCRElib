import os, sys, re
import configparser
import optparse
import subprocess
import shutil

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
from SBILib.data import dna_complementary
from SBILib.structure import PDB
from SBILib.structure.chain import ChainOfNucleotide
from SBILib.structure.residue import ResidueOfNucleotide

# Import my modules #
from ModCRElib.structure.contacts import interface
from ModCRElib.structure.dna import x3dna
from ModCRElib.web import md2p as MD2P

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for standalone DNA structure generation.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for FASTA input and output settings.

    """

    parser = optparse.OptionParser("python build_dna.py -i FASTA_FILE [-o OUTPUT_FILE --conformation DNA_TYPE --dummy DUMMY_DIR]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", default=None, type="string", dest="fasta_file", help="DNA FastA file", metavar="{filename}")
    parser.add_option("-o", "--output", default=None, action="store", type="string", dest="output", help="Output file with the structure (default is the same as FastA name with pdb extension)", metavar="{filename}")
    parser.add_option("-c", "--conformation", default="B", action="store", type="string", dest="conformation", help="DNA conformation is either A, B, C, Z or B_bent (for nuclesome). Default is B.", metavar="{string}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
  

    (options, args) = parser.parse_args()

    if (options.fasta_file is None ):
        parser.error("missing arguments: type option \"-h\" for help")

    return options




def main():
    """
    Run the command-line DNA-building workflow from FASTA entries.

    Workflow:
        1. Parse command-line options.
        2. Iterate over FASTA sequences.
        3. Build one DNA structure per sequence with selected conformation.
        4. Write output PDB files.

    Args:
        None.

    Returns:
        None. Generated DNA PDB files are written to disk.
    """


    # Arguments & Options #
    options = parse_options()

    verbose    = options.verbose
    dummy_dir  = options.dummy_dir
    dna_str    = options.conformation
    fasta_file = options.fasta_file
    output     = options.output
    begining   = None
    ending     = None

    counter = 0
    for name,dna_seq in functions.parse_fasta_file(fasta_file):
        dna_file = MD2P.build_dna_structure(dna_seq, dna_str, dummy_dir,  begining, ending, verbose)
        if output is None:
           output_file = name.replace(" ","_").replace("|","_")+".pdb"
           if verbose: print("PDB Structure: \t %s"%( output_file ))
           shutil.copy(dna_file,output_file)
        else:
           counter += 1
           single_file = output.rstrip(".pdb") + ".pdb"
           output_file = output.rstrip(".pdb") + "_" +str(counter) + ".pdb"
           if verbose: print("PDB Structure: \t %s"%( output_file ))
           shutil.copy(dna_file,output_file)

    if counter==1: shutil.move(output_file,single_file)
    if verbose: print("Done")



#---------------#
# Main          #
#---------------#

if __name__ == "__main__":
    main()
