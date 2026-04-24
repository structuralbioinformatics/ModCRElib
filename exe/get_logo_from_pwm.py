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
elif os.path.exists(os.path.join(exe_path,"..","..","..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..","..","..")
elif os.path.exists(os.path.join(exe_path,"..","..","..","..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..","..","..","..")
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


#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.
    """

    parser = optparse.OptionParser("python get_logo_from_pwm.py -i input_pwm_or_msa [-o output_directory --dummy dummy_dir]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", "--input", action="store", default=None, type="string", dest="input_pwm", help="Input motif file or directory. Accepted formats: MSA, PWM, MEME and TXT", metavar="{filename/directory}")
    parser.add_option("-o", "--outdir", action="store", default=None, type="string", dest="outdir", help="Directory for outputs (input rootname is used as output prefix)", metavar="{directory}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")

    (options, args) = parser.parse_args()

    if options.input_pwm is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


def get_motif_name(motif_file):
    """Get the motif name from the input file name."""
    name = os.path.basename(motif_file)
    if name.endswith(".meme.s"):
       return name[:-7]
    if name.endswith(".meme"):
       return name[:-5]
    if name.endswith(".msa"):
       return name[:-4]
    if name.endswith(".pwm"):
       return name[:-4]
    if name.endswith(".txt"):
       return name[:-4]
    return ".".join(name.split(".")[:-1]) if "." in name else name


def guess_input_option(motif_file):
    """Guess the motif format from the file extension."""
    name = os.path.basename(motif_file)
    if name.endswith(".meme") or name.endswith(".meme.s"):
       return "meme"
    if name.endswith(".pwm"):
       return "pwm"
    if name.endswith(".txt"):
       return "txt"
    return "msa"


def get_logo_from_pwm(motif_file,outdir,dummy_dir,verbose=False):
    """Read a motif file and create the corresponding logo files."""
    if not os.path.exists(motif_file):
       if verbose: sys.stdout.write("Input file %s does not exist\n"%(motif_file))
       return False

    motif_name = get_motif_name(motif_file)
    input_option = guess_input_option(motif_file)

    try:
       msa_obj = PWM.nMSA(motif_file,motif_name,input_option)
       msa_obj.set_sequences()
    except Exception as e:
       sys.stdout.write("Failed to read motif %s with format %s: %s\n"%(motif_file,input_option,e))
       return False

    if outdir is None:
       outdir = os.path.dirname(os.path.abspath(motif_file))
    if outdir == "":
       outdir = "."
    if not os.path.exists(outdir):
       os.makedirs(outdir)

    logo_file = os.path.join(outdir,motif_name+".logo")
    if verbose: sys.stdout.write("-- Write LOGOS %s...\n"%(logo_file))
    PWM.write_logo(msa_obj,logo_file,dummy_dir)
    return True


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

  # Arguments & Options #
  options = parse_options()
  input_file=options.input_pwm
  outdir    =options.outdir
  verbose   =options.verbose
  dummy_dir =options.dummy_dir

  if os.path.isdir(input_file):
    done=set()
    for pwm_file in sorted(os.listdir(input_file)):
        if not pwm_file.endswith(".msa") and not pwm_file.endswith(".pwm") and not pwm_file.endswith(".meme") and not pwm_file.endswith(".meme.s") and not pwm_file.endswith(".txt"):
           continue
        motif_name = get_motif_name(pwm_file)
        if motif_name in done:
           continue
        full_path = os.path.join(input_file,pwm_file)
        if verbose: sys.stdout.write("Get LOGO from: %s\n"%(full_path))
        if get_logo_from_pwm(full_path,outdir,dummy_dir,verbose):
           done.add(motif_name)
           sys.stdout.write("Done for %s\n"%(full_path))
        else:
           sys.stdout.write("Failed for %s\n"%(full_path))
  else:
    if get_logo_from_pwm(input_file,outdir,dummy_dir,verbose):
       sys.stdout.write("Done for %s\n"%(input_file))
    else:
       sys.stdout.write("Failed for %s\n"%(input_file))
