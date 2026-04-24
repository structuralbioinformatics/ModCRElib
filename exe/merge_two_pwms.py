"""
Merge two DNA PWMs using TOMTOM-guided alignment and optional structural context.

This script takes a main PWM and a second PWM, aligns them with the helper
logic in :mod:`ModCRElib.msa.merge_pwms`, and writes one or more merged motifs
covering the full, main-only, or core-overlap regions.
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
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.protein import dssp,tmalign
from ModCRElib.structure.threading import threader, threading_to_triads
from ModCRElib.potential import spotentials
from ModCRElib.profile import tomtom
from ModCRElib.msa import merge_pwms
from ModCRElib.msa import pwm_pbm as PWM


# Defined maximum size of nMSA, minimum must be 1000
try:
  maxsize= int(config.get("Parameters", "max_sequences_in_MSA"))
except:
  maxsize=5000
if maxsize < 1000: maxsize=1000


#-------------#
# Options     #
#-------------#

def parse_options():
    """
    Parse the command line arguments for two-PWM merging.

    TOMTOM is used internally with the main PWM as the database/query anchor
    and the second PWM as the motif to align and merge.

    Returns:
        optparse.Values: Parsed options describing the input PWMs, structural
        model, output mode, and statistical-potential settings (required to execute funcion merge_two_pwm)
    """

    parser = optparse.OptionParser("python merge_pwm.py -i input_main_pwm -j input_add_pwm -w weight [-o output_file --dummy dummy_dir -w weight --pbm pbm_dir --pdb pdb_dir  --refine refine ]  ")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", default=None, type="string", dest="main_pwm", help="MAIN PWM associated with a model structure ", metavar="{filename}")
    parser.add_option("-j", action="store", default=None, type="string", dest="second_pwm", help="SECOND PWM  ", metavar="{filename}")
    parser.add_option("--model",action="store", default=None, type="string", dest="model", help="Model structure corresponding to the MAIN PWM (it may be in the format of thread)", metavar="{filename}")
    parser.add_option("-w",action="store", default=1,type="float",dest="weight", help="Weight of the second PWM, the weigh for the MAIN PWM is always 1 (defaults = 1)")
    parser.add_option("--sig",action="store", default=0.05,type="float",dest="significance",help="Limit of significance (p-value) required in TOMTOM to merge the PWMs (defaults is 0.05)")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", default=None, help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("--main",default=False, action="store_true", dest="main", help="Generate a PWM for the positions of the MAIN PWM associated with structure (default is False)", metavar="{boolean}")
    parser.add_option("--full",default=False, action="store_true", dest="full", help="Generate a PWM to cover all positions of the MAIN and SECOND PWM  (default is False)", metavar="{boolean}")
    parser.add_option("--core",default=False, action="store_true", dest="core", help="Generate a PWM for the overlapping positins of the MAIN and SECOND PWM  (default is False)", metavar="{boolean}")
    parser.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    parser.add_option("--meme", default=False, action="store_true", dest="meme", help="Use 'uniprobe2meme' to calculate the PWM matrix for 'FIMO' (default = False)")
    parser.add_option("--threading", default=False, action="store_true", dest="threading", help="Input file is a threading file of a PDB structure that (MUST!) exist in the PDB folder of ModCRE (default = False)")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment_restrict", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. Example: '45_A-48_A;50_A-51_A' (Default is None it applies to all amino-acids)")
    parser.add_option("--binding", default=None, action="store", type="string", dest="binding_restrict", help="Binding site of DNA to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d of the forward chain (first in PDB). (Default is None it applies to all nucleotides)")
    parser.add_option("-o", action="store", default="output_pwm", type="string", dest="output_file", help="Output pwm (default = 'output_pwm')", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")

    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("-t", action="store", default=None, type="float", dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default is taken from config file)", metavar="{float}")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="computation", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' or family-specific radius from configuration", metavar="{float}")
    group.add_option("--methylation",default=False, action="store_true", dest="methylation", help="Flag to use methylated cytosines with binding/non-binding specificity (default=False)")
    group.add_option("--refine",default=0, action="store", type="int", dest="refine", help="Level to refine the MSA and PWM scoring DNA binding sequences with full length (default=0, 1 refine and trim, 2 refine, rescale and cut-off)", metavar="{integer}")

    parser.add_option_group(group)
    (options, args) = parser.parse_args()

    if options.main_pwm is None or options.second_pwm is None or options.model is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


def main():
    """
    Run the command-line workflow for merging two PWMs.

    Workflow:
        1. Parse the input PWMs, structural model, and merge settings.
        2. Load optional family, fragment, and binding restrictions.
        3. Align and merge the motifs through ``merge_pwms.merge_two_pwm``.
        4. Write the requested output motifs and print their consensus sequences.

    Returns:
        None. Output files are written next to the requested output prefix.
    """
    # Arguments & Options #
    options = parse_options()

    # Initialize #

    # Dummy/temporary folder
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    dummy_dir=os.path.abspath(options.dummy_dir)

    # input structures and folders of ModCRE
    pbm_dir=None
    if options.pbm_dir is not None: pbm_dir = os.path.abspath(options.pbm_dir)
    pdb_dir=None
    threading = options.threading
    if options.threading:
       input_threading_file=options.model
       input_pdb_file=None
       if not input_threading_file.startswith("/"): input_threading_file=os.path.abspath(options.model)
    else:
       input_pdb_file=options.model
       input_threading_file=None
       if not input_pdb_file.startswith("/"): input_pdb_file=os.path.abspath(options.model)
    if options.pdb_dir is not None: pdb_dir = os.path.abspath(options.pdb_dir)

    #input PWMs
    main_pwm=os.path.abspath(options.main_pwm)
    second_pwm=os.path.abspath(options.second_pwm)
    weight=float(options.weight)
    meme = options.meme
    significance = float(options.significance)

    # Input potentials
    radius = options.radius
    potential_file = options.potential_file
    split_potential= options.split_potential
    auto_mode = options.auto_mode
    family_potentials = options.family_potentials
    pbm_potentials = options.pbm_potentials
    score_threshold = options.score_threshold
    taylor_approach = options.taylor_approach
    pmf = options.pmf
    bins_approach = options.computation
    known = options.known

    # Special input features
    methylation = options.methylation
    refine      = options.refine
    fragment_restrict=None
    if options.fragment_restrict is not None:
       fragment_restrict={}
       segment_fragment=options.fragment_restrict.split(";")
       for interval in segment_fragment:
           ac,bc = interval.split("-")
           chain1=chain2=" "
           if len(ac.split("_"))==1: a=ac
           if len(bc.split("_"))==1: b=bc
           if len(ac.split("_"))>1: a,chain1=ac.split("_")
           if len(bc.split("_"))>1: b,chain2=bc.split("_")
           if chain1==chain2:
              fragment_restrict.setdefault(chain1,[]).append((int(a),int(b)))

    binding_restrict=None
    if options.binding_restrict is not None:
       binding_restrict=[]
       segment_binding=options.binding_restrict.split(";")
       for interval in segment_binding:
           a,b = interval.split("-")
           aa = int(a)
           bb = int(b)
           binding_restrict.append((aa,bb))

    # Get families from input
    families = {}
    if options.verbose:sys.stdout.write("Check families...\n")
    if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
       sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
    else:
      for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family

    #  ouput names 
    output_pwm = os.path.basename(options.output_file)
    output_dir = os.path.dirname(options.output_file)
    output_file = options.output_file
    verbose = options.verbose
    main=options.main
    full=options.full
    core=options.core
    n=0
    if full:   n=n+1
    if core:   n=n+1
    if n>1:
         output_full_name   = output_pwm + ".full"
         output_core_name   = output_pwm + ".core"
         output_main_name   = output_pwm + ".main"
         output_name = None
    else:
         output_name = msa_obj.get_motif()
    

    result_full_pwm, msa_obj, over5, overlap = merge_pwms.merge_two_pwm(main_pwm,second_pwm, weight, refine, significance, input_threading_file, input_pdb_file, pdb_dir, pbm_dir, threading, fragment_restrict, binding_restrict, output_file, families, radius, potential_file, split_potential,  auto_mode,  family_potentials,  pbm_potentials,  score_threshold,  taylor_approach,  pmf, bins_approach, known,   methylation, dummy_dir, verbose)

    if result_full_pwm is None or msa_obj is None:
        print("The alignment with TOTOM was not significant")
        print("Done")
        exit()

    # Write Output files
    if output_name is not None:
        if full:
            result_pwm = result_full_pwm
            result_pwm.set_motif(output_name)
            merge_pwms.write_pwm(result_pwm,output_dir,output_name,meme,dummy_dir,verbose)
            print("Sequence FULL ",result_pwm.get_main_sequence())
        if main:
            result_pwm = msa_obj
            result_pwm.set_motif(output_name)
            print("Sequence MAIN ",result_pwm.get_main_sequence())
            merge_pwms.write_pwm(result_pwm,output_dir,output_name,meme,dummy_dir,verbose)
        if core:
            result_pwm = result_full_pwm.section(over5,(over5+overlap))
            result_pwm.set_motif(output_name)
            print("Sequence CORE ",result_pwm.get_main_sequence())
            merge_pwms.write_pwm(result_pwm,output_dir,output_name,meme,dummy_dir,verbose)
    else:
         result_main_pwm=msa_obj
         result_core_pwm=result_full_pwm.section(over5,(over5+overlap))
         result_full_pwm.set_motif(output_full_name)
         result_main_pwm.set_motif(output_main_name)
         result_core_pwm.set_motif(output_core_name)
         if full: merge_pwms.write_pwm(result_full_pwm,output_dir,output_full_name,meme,dummy_dir,verbose)
         if main: merge_pwms.write_pwm(result_main_pwm,output_dir,output_main_name,meme,dummy_dir,verbose)
         if core: merge_pwms.write_pwm(result_core_pwm,output_dir,output_core_name,meme,dummy_dir,verbose)
         if full: print("Sequence FULL ",result_full_pwm.get_main_sequence())
         if main: print("Sequence MAIN ",result_main_pwm.get_main_sequence())
         if core: print("Sequence CORE ",result_core_pwm.get_main_sequence())
    print("Done")


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
