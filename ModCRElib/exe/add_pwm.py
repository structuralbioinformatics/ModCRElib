"""
Concatenate two PWMs/MSAs into a single motif and write the result.

This script is a lightweight wrapper around :mod:`ModCRElib.msa.pwm_pbm`. It
loads two input motifs, concatenates them end-to-end with ``MSA.__add__``, and
writes the resulting motif together with an optional logo.
"""

import os, sys, re
import configparser
import optparse

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

# Import my modules #
from ModCRElib.msa import pwm_pbm as PWM

def parse_options():
    """
    Parse the command line options for concatenating two motifs.

    Returns:
        optparse.Values: Parsed CLI options describing the two input files,
        their formats, the output format/prefix, and whether protein mode is used.
    """

    parser = optparse.OptionParser("python add_pwm.py -a input_pwm_file -b  input_pwm_file  [-o output_pwm_file --fa format file A --fb format file B --fmt format output --name motif_name]")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-a", action="store", default=None, type="string", dest="input_pwm_A", help="PWM file (A) in any format 'pwm' 'msa' or 'meme')", metavar="{filename}")
    parser.add_option("-b", action="store", default=None, type="string", dest="input_pwm_B", help="PWM file (B) in any format 'pwm' 'msa' or 'meme')", metavar="{filename}")
    parser.add_option("-o", action="store", default="output_pwm", type="string", dest="output_file", help="Output file (default = 'output_pwm')", metavar="{filename}")
    parser.add_option("--name", action="store", default=None, type="string", dest="name_motif", help="Name of output motif (default is None and it merges the input)", metavar="{string}")
    parser.add_option("--fmt", action="store", default="meme", type="string", dest="format_pwm", help="Format of the output PWM (default is 'meme')", metavar="{string}")
    parser.add_option("--fa", action="store", default="meme", type="string", dest="format_A", help="Format of the input file (A) PWM (default is 'meme')", metavar="{string}")
    parser.add_option("--fb", action="store", default="meme", type="string", dest="format_B", help="Format of the input file (B) PWM (default is 'meme')", metavar="{string}")
    parser.add_option("-p", "--protein", default=False, action="store_true", dest="protein", help="PWMs are from protein sequences (default = False, therefore sequences are polymers of nucleotides)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
 
    (options, args) = parser.parse_args()

    if (options.input_pwm_A is None or options.input_pwm_A is None):
        parser.error("missing arguments: type option \"-h\" for help")
     
    return options

def main():
    """
    Run the command-line workflow for concatenating two motifs.

    Workflow:
        1. Read motif A and motif B in the requested formats.
        2. Concatenate them with ``msa_obj_a + msa_obj_b``.
        3. Write the merged motif and create a logo when possible.

    Returns:
        None. Output files are written next to the requested output prefix.
    """
    # Arguments & Options #
    options = parse_options()

    # Read input files
    if  options.verbose:sys.stdout.write("\t--Get PWM file A %s ...\n"% options.input_pwm_A)
    if  options.protein:
        msa_obj_a = PWM.pMSA(options.input_pwm_A,None,options.format_A)
    else:
        msa_obj_a = PWM.nMSA(options.input_pwm_A,None,options.format_A)
    if  options.verbose:sys.stdout.write("\t--Get PWM file B %s ...\n"% options.input_pwm_B)
    if  options.protein:
        msa_obj_b = PWM.pMSA(options.input_pwm_B,None,options.format_B)
    else:
        msa_obj_b = PWM.nMSA(options.input_pwm_B,None,options.format_B)

    # Add PWMs
    msa_obj = msa_obj_a + msa_obj_b
    if options.name_motif is not None:
       msa_obj.set_motif(options.name_motif)

    # Write PWM #
    if  options.verbose:sys.stdout.write("\t--Write PWM ...\n")
    format_file = options.format_pwm
    pwm_file    = options.output_file + "." + format_file
    logo_file   = options.output_file + ".logo"
    logo_gapped_file   = options.output_file + ".gapped.logo"
    if not os.path.exists(pwm_file):
       if  options.verbose:sys.stdout.write("\t\t--PWM in %s format...\n"%format_file)
       msa_obj.write(pwm_file, option=format_file)
    if not os.path.exists(logo_file+".fwd.png") or  not os.path.exists(logo_file+".rev.png"):
       if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
       if  options.verbose:sys.stdout.write("\t\t--Logos...\n")
       if  options.protein:
         PWM.write_protein_logo(msa_obj,logo_file, logo_gapped_file, options.dummy_dir)
       else:
         PWM.write_logo(msa_obj,logo_file, options.dummy_dir)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()


    
