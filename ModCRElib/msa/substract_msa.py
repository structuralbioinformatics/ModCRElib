"""
Subtract one motif from another and write the resulting differential motif.

This script is a small command-line wrapper around :mod:`ModCRElib.msa.pwm_pbm`.
It loads two motifs in a supported format, computes their difference using the
``MSA.difference`` method, regenerates a synthetic sequence set, and writes the
result in MEME, PWM, and MSA formats together with a sequence logo.
"""

import os, sys, re
from collections import Counter
import configparser
import itertools
import numpy
import optparse
import subprocess
import time
import random

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
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases,dna_complementary
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.msa import pwm_pbm as PWM



# Defined maximum size of nMSA, minimum must be 1000
try:
  maxsize= int(config.get("Parameters", "max_sequences_in_MSA"))
except:
  maxsize=5000
if maxsize < 1000: maxsize=1000
  

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse the command line arguments for motif subtraction.

    Returns:
        optparse.Values: Parsed options describing the input motifs, output
        prefix, input format, overlap, mode (DNA/protein), and logo settings.
    """

    parser = optparse.OptionParser("python substract_msa.py -i main_pwm -j background_pwm  [-o output_file -m [MAXSIZE] --overlap [OVERLAP] --enhance --dummy dummy_dir --verbose --protein ]  ")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", "--main", action="store", default=None, type="string", dest="main_pwm", help="Main PWM ", metavar="{filename}")
    parser.add_option("-j", "--background", action="store", default=None, type="string", dest="background_pwm", help="Background PWM to substract in the Main PWM", metavar="{filename}")
    parser.add_option("-m", "--maxsize", action="store", default=None, type="string", dest="maxsize", help="Maximum number of sequences in a MSA file (default uses max_sequences_in_MSA of configuration) ", metavar="{int}")
    parser.add_option("--format", action="store", default="msa", type="string", dest="format", help="Format of the PWM can only be msa/meme/pwm/txt (default is msa)", metavar="{string}")
    parser.add_option("--overlap",action="store", default=None,type="string", dest="overlap", help="Overlap between the Main PWM and the background PWM (default is None and it overlaps the wole Background PWM on the Main PWM) ", metavar="{int}")
    parser.add_option("-o", action="store", default="output_msa", type="string", dest="output_file", help="Output MSA (default = 'output_msa')", metavar="{filename}")
    parser.add_option("-e", "--enhance", default=False, action="store_true", dest="enhance", help="Emhance mode to highlight the residues with maximum score in the PWM (default = False)", metavar="{boolean}")
    parser.add_option("-p", "--protein", default=False, action="store_true", dest="protein", help="DNA/protein mode. By default the PWMs are for DNA (default = False)", metavar="{boolean}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
  
 
    (options, args) = parser.parse_args()

    if options.main_pwm is None or options.background_pwm is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options


def main():
    """
    Run the command-line workflow for subtracting two motifs.

    Workflow:
        1. Read the main and background motifs in the selected format.
        2. Compute a differential motif with ``MSA.difference``.
        3. Rebuild a synthetic sequence set from the resulting PWM.
        4. Write the result as ``.meme``, ``.pwm``, ``.msa``, and a logo image.

    Returns:
        None. Files are written next to the requested output prefix.
    """
    # Arguments & Options #
    options = parse_options()
    main        =options.main_pwm
    background  =options.background_pwm
    output      =options.output_file
    verbose     =options.verbose
    pwm_format  =options.format
    protein     =options.protein
    enhance     =options.enhance
    dummy_dir   =options.dummy_dir

    if options.maxsize is not None:
       maxsize=int(options.maxsize)

    overlap=None
    if options.overlap is not None:
       overlap=int(options.overlap)

    if pwm_format != "meme" and pwm_format != "msa" and pwm_format != "pwm":
       pwm_format="msa"

    if protein:
       main_pwm        = PWM.pMSA(main,"main",option=pwm_format)
       background_pwm  = PWM.pMSA(background,"background",option=pwm_format)
    else:
       main_pwm        = PWM.nMSA(main,"main",option=pwm_format)
       background_pwm  = PWM.nMSA(background,"background",option=pwm_format)

    difference = main_pwm.difference(background_pwm,overlap=overlap,enhance=enhance)

    difference.set_motif(os.path.basename(output))
    difference.clean_sequences()
    difference.set_sequences(maxsi=maxsize)

    if verbose: sys.stdout.write("-- Write meme file %s...\n"%(output.rstrip(".msa")+".meme"))
    difference.write(output.rstrip(".msa")+".meme",option="meme")
    if verbose: sys.stdout.write("-- Write PWM file %s...\n"%(output.rstrip(".msa")+".pwm"))
    difference.write(output.rstrip(".msa")+".pwm",option="pwm")
    if verbose: sys.stdout.write("-- Write MSA file %s...\n"%(output.rstrip(".msa")+".msa"))
    difference.write(output.rstrip(".msa")+".msa",option="msa")
    if verbose: sys.stdout.write("-- Write LOGOS %s...\n"%(output.rstrip(".msa")+".logo"))
    PWM.write_protein_logo(difference,output.rstrip(".msa")+".logo",dummy_dir)
    if verbose: print("Done")    


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()

    


