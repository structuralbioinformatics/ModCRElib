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
    This function parses the command line arguments and returns an optparse
    object.
    
    """

    parser = optparse.OptionParser("python get_logo_msa.py -i input_msa  [-o output_directory --dummy dummy_dir ]  ")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", "--msa", action="store", default=None, type="string", dest="input_msa", help="Alignment of modelled TF sequences associated with structural models. Format MSA or MEME or PWM. If the input is a folder it uses the MSA  or the MEME files ", metavar="{filename}")
    parser.add_option("-m", "--meme", default=False, action="store_true", dest="meme", help="Write/rewrite MEME file (default = False)", metavar="{boolean}")
    parser.add_option("-p", "--pwm", default=False, action="store_true", dest="pwm", help="Write/rewrite PWM (default = False)", metavar="{boolean}")
    parser.add_option("-o", "--outdir", action="store", default=None, type="string", dest="outdir", help="Directory for outputs (MSA file is used as a rootname) ", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
  
 
    (options, args) = parser.parse_args()

    if options.input_msa is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_logo_msa(msa_file,outdir,meme_do,pwm_do):
    name_motif=".".join(os.path.basename(msa_file).split(".")[:-1])
    try:
       msa_obj = PWM.nMSA(msa_file)
    except:
       try:
          msa_obj = PWM.nMSA(msa_file,name_motif,"meme")
       except:
          try:
            msa_obj = PWM.nMSA(msa_file,name_motif,"pwm")
          except:
            print("Failed to read input, check format is MSA, PWM or MEME")
            return False 
    if not os.path.exists(outdir): os.makedirs(outdir)
    meme_file = os.path.join(outdir,name_motif+".meme")
    pwm_file  = os.path.join(outdir,name_motif+".pwm")
    logo_file = os.path.join(outdir,name_motif+".logo")
    if not os.path.exists(meme_file) and meme_do:
       if verbose: sys.stdout.write("-- Write meme file %s...\n"%(meme_file))
       msa_obj.write(meme_file,option="meme")
    if not os.path.exists(pwm_file) and pwm_do:
       if verbose: sys.stdout.write("-- Write PWM file %s...\n"%(pwm_file))
       msa_obj.write(pwm_file,option="pwm")
    if verbose: sys.stdout.write("-- Write LOGOS %s...\n"%(logo_file))
    PWM.write_logo(msa_obj,logo_file,dummy_dir)
    return True

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

  # Arguments & Options #
  options = parse_options()
  input_file=options.input_msa
  outdir    =options.outdir
  verbose   =options.verbose
  pwm_do    =options.pwm
  meme_do   =options.meme
  dummy_dir =options.dummy_dir
  pwms=set()
  pwm_types=["msa","meme","pwm"]
  done=set()
  if os.path.isdir(input_file):
    for msa_file in os.listdir(input_file):
        if msa_file.endswith("msa"): pwms.add(".".join(os.path.basename(msa_file).split(".")[:-1]))
    for motif in pwms:
        if motif in done: continue
        for pwm_type in pwm_types:
            if motif in done: continue
            msa_file = motif+"."+pwm_type
            if verbose: print("Get LOGO from: ",msa_file)
            if get_logo_msa(os.path.join(input_file,msa_file),outdir,meme_do,pwm_do): 
               done.add(motif)
               print("Done for %s"%msa_file)
            else:
               print("Failed for %s"%msa_file)
  else:
    if get_logo_msa(input_file,outdir,meme_do,pwm_do):
       print("Done for %s"%input_file)
    else:
       print("Failed for %s"%input_file)



    
