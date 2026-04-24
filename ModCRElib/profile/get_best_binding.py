"""
Locate the best PWM-supported binding region in a DNA sequence for a TF-DNA structure.

This script scans an input DNA sequence with a PWM, selects the best supported
binding region, optionally rescoring nearby alternatives with structural
potentials, and finally generates threaded TF-DNA models for the selected site.
"""

import os, sys, re
import configparser
import optparse
import shutil
import subprocess
import gzip
import traceback
import configparser
import argparse
import shutil
import pandas as pd
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

# Get python path #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Import my functions #
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.structure.dna import x3dna, model_dna
from ModCRElib.structure.protein import dssp
from ModCRElib.structure.threading import pdb2thread
from ModCRElib.potential import spotentials
from ModCRElib.profile import fimo, scorer
from ModCRElib.msa import  pwm_pbm as PWM

# Imports jbonet's module #
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases,dna_complementary
from SBILib.structure import PDB

#-------------#
# Functions   #
#-------------#

def get_threads_with_best_binding(binding,delta,standard,pwm_file,dna_sequence,pdb_file,input_dir,pdb_dir,pbm_dir,radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins_approach, known, p_val_max, verbose=False,dummy_dir="/tmp",fragment_restrict=None, binding_restrict=None,clean=False,wobble=False,methylation=False):
    """
    Find the best-supported binding site and build threaded models for it.

    Args:
        binding (str | None): Optional binding sequence that must be matched.
        delta (int): Window extension around the initial FIMO match to rescore.
        standard (bool): If True, use only standard nucleotides while rescoring.
        pwm_file (str): PWM file in MEME format.
        dna_sequence (str): DNA sequence to scan.
        pdb_file (str): Input TF-DNA PDB structure.
        input_dir (str): Output/input directory for generated threading files.
        pdb_dir, pbm_dir (str): ModCRE data directories.
        radius, potential_file, split_potential, auto_mode, family_potentials,
        pbm_potentials, score_threshold, taylor_approach, pmf, bins_approach,
        known, p_val_max, verbose, dummy_dir, fragment_restrict,
        binding_restrict, clean, wobble, methylation: Workflow parameters
        forwarded to the structural scoring and threading helpers.

    Returns:
        list: Threading objects/files returned by ``pdb2thread.pdb2thread`` for
        the selected best binding region.
    """

    #Get msa_objs
    if verbose: sys.stdout.write("\t-- Assign MSA objs to %s \n"%pwm_file)
    msa_objs={}
    msa_obj=PWM.nMSA(pwm_file,option="meme")
    for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
        msa_objs.setdefault(sp,msa_obj)

    if binding is not None: 
       binding   = binding.upper()
       reverse   = functions.reverse_dna(binding)
    else:
       reverse   = None

    #Get families
    families = {}
    if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
         if verbose: sys.stdout.write("\t-- Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
    else:
         for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
           if line.startswith("#"): continue
           pdb_chain, family = line.split(";")
           families[pdb_chain] = family

    # Get FIMO and select binding fragment
    dummy_dna_file = os.path.join(dummy_dir,"dna_sequence.fa")
    #if methylation:
    #  shutil.copy(dna_file,dummy_dna_file)
    #else:
    f=open(dummy_dna_file,"w")
    f.write(">%s\n%s\n"%("DNA",dna_sequence.upper().replace("X","C").replace("O","C").replace("J","G").replace("Q","G")))
    f.close()
    if verbose:print("Run FIMO ")
    fimo_obj = fimo.get_fimo_obj(pwm_file, dummy_dna_file)
    found_sequence = None
    start=None
    for hit_obj in fimo_obj.get_hits(sort=True):
           p_val = hit_obj.get_p_value()
           if p_val>p_val_max: continue
           strand= hit_obj.get_strand()
           if strand == "+": check_seq   = hit_obj.get_sequence()
           if strand == "-": check_seq   = functions.reverse_dna(hit_obj.get_sequence())
           seq = hit_obj.get_sequence()
           start = int(hit_obj.get_start())
           test_sequence =dna_sequence[start-1:start+len(seq)-1]
           if binding is not None:
              if binding in check_seq or check_seq in binding:
                 if (has_methylation(dna_sequence) and has_methylation(test_sequence)) or (not has_methylation(dna_sequence) and not has_methylation(test_sequence)):
                    found_sequence=seq
                    bind_sequence = dna_sequence[start-1:start+len(seq)-1]
                    if verbose: print(("\t-- FIMO RESULT %s (%d) original seq %s strand %s (check %s ) P-value %s "%(hit_obj.get_sequence(),start,bind_sequence,strand,check_seq,str(hit_obj.get_p_value()))))
                    break
              elif reverse in check_seq or check_seq in reverse:
                 if (has_methylation(dna_sequence) and has_methylation(test_sequence)) or (not has_methylation(dna_sequence) and not has_methylation(test_sequence)):
                    found_sequence=seq
                    bind_sequence = dna_sequence[start-1:start+len(seq)-1]
                    if verbose: print(("\t-- FIMO RESULT %s (%d) original seq %s strand %s (check %s ) P-value %s "%(hit_obj.get_sequence(),start,bind_sequence,strand,check_seq,str(hit_obj.get_p_value()))))
                    break
              else:
                 continue
    if found_sequence is None:
       print("No matches found in the DNA sequence")
       raise ValueError("No matches found in the DNA sequence")

    # Get general data of structure
    if verbose: print("Get structure data")
    # Get PDB object
    if verbose: sys.stdout.write("\t-- Reading PDB file %s ...\n"%pdb_file)
    pdb_obj = PDB(pdb_file)
    if clean:
       pdb_obj.clean()
       dummy_pdb_file = os.path.join(dummy_dir,"pdb_file.pdb")
       pdb_obj.write(dummy_pdb_file)
    else:
       dummy_pdb_file = pdb_file
    pdb_name= pdb_obj.id
    # Get DSSP object #
    if verbose: sys.stdout.write("\t\t-- Get DSSP ...\n")
    dssp_obj = dssp.get_dssp_obj(dummy_pdb_file)
    # Get X3DNA object #
    if verbose: sys.stdout.write("\t\t-- Get DNA object ...\n")
    x3dna_obj = x3dna.get_x3dna_obj(dummy_pdb_file, dummy_dir)
    # Get contacts object #
    if verbose: sys.stdout.write("\t\t-- Get contacts ...\n")
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
    # Skip if no contacts #
    if len(contacts_obj._contacts) == 0: exit("No protein-DNA contacts found!")
    # Get interface object #
    if verbose: sys.stdout.write("\t\t-- Get interface ...\n")
    interface_obj = interface.get_interface_obj(pdb_obj, x3dna_obj, contacts_obj, os.path.abspath(dummy_dir))
    # Get DNA template
    dna_template_seq=[]
    for chain_id in pdb_obj.chain_identifiers:
        chain=pdb_obj.get_chain_by_id(chain_id)
        if chain.chaintype=="N": 
            dna_template_seq.append(chain.gapped_nucleotide_sequence())
    size = len(dna_template_seq[0])
    if verbose:sys.stdout.write("\t-- DNA template (%3d)  : %s\n"%(size,dna_template_seq[0]))
    if verbose:sys.stdout.write("\t-- Found sequence (%3d): %s\n"%(len(found_sequence),found_sequence))
    if verbose:sys.stdout.write("\t-- Bind sequence (%3d) : %s\n"%(len(bind_sequence),bind_sequence))
    #size = interface_obj.get_interface_length()
    # Check DNA sequence length in PDB
    # If the length of DNA in PDB matches the binding, return solution
    if len(found_sequence) == interface_obj.get_interface_length() and size==len(found_sequence) and not wobble:
       try:
          if verbose:sys.stdout.write("\t-- Model DNA ...\n")
          pdb_obj_new = model_dna.get_dna_model_pdb_obj(dummy_pdb_file,found_sequence, x3dna_obj, interface_obj, interface_range=None, dummy_dir=dummy_dir)
       except Exception as e:
          print(("Error: %s"%e))
          raise ValueError(e)
       try:
          if verbose:sys.stdout.write("\t-- Make threaded object ...\n")
          thread_objs = pdb2thread.pdb2thread(pdb_obj_new,pdb_name,x3dna_obj,input_dir,bind_sequence,verbose)
       except Exception as e:
          print(("Error: %s"%e))
          raise ValueError(e)
    # If the length of the binding differs from PDB
    else:
       # Load statistical potential
       if verbose: sys.stdout.write("\t-- Load potentials...\n")
       try:
         potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins_approach, known, None, dummy_dir,verbose)
       except Exception as e:
          raise Exception("Failed to get potentials")

       # Prepare the starting sequences
       max_score=0
       start_position = start
       if start + size -1 <= len(dna_sequence):
           bind_sequence  = dna_sequence[start-1:start+size-1]
       else:
           bind_sequence  = dna_sequence[len(dna_sequence)-size:len(dna_sequence)]
       found_sequence = bind_sequence.upper().replace("X","C").replace("O","C").replace("J","G").replace("Q","G")
       if verbose:sys.stdout.write("\t-- Found sequence (%3d - %3d): %s\n"%(start_position,start_position+len(found_sequence),found_sequence))
       if verbose:sys.stdout.write("\t-- Bind sequence (%3d - %3d) : %s\n"%(start_position,start_position+len(bind_sequence),bind_sequence))
       # Calculate energy scores #
       if verbose: sys.stdout.write("\t-- Calculate scores...\n")
       scr=scorer.score(potential=split_potential)
       # Loop to select DNA sequences and get best scoring sequence
       for i in range(max(0,start-delta),min((start+delta),(len(dna_sequence)-size-1))):
           test_sequence = dna_sequence[i:i+size]
           if has_methylation(dna_sequence) and not has_methylation(test_sequence): continue
           temp_sequence = test_sequence.upper().replace("X","C").replace("O","C").replace("J","G").replace("Q","G")
           pdb_test      = model_dna.get_dna_model_pdb_obj(dummy_pdb_file,temp_sequence, x3dna_obj, interface_obj, interface_range=None, dummy_dir=dummy_dir)
           if clean: pdb_test.clean()
           dummy_file    = os.path.join(dummy_dir,"test.pdb")
           pdb_test.write(dummy_file,force=True)
           # Get DSSP object #
           if verbose: sys.stdout.write("\t\t-- Get DSSP ...\n")
           dummy_dssp_obj   = dssp.get_dssp_obj(dummy_file)
           # Get X3DNA object #
           if verbose: sys.stdout.write("\t\t-- Get DNA object ...\n")
           dummy_x3dna_obj = x3dna.get_x3dna_obj(dummy_file, dummy_dir)
           # Get contacts object #
           if verbose: sys.stdout.write("\t\t-- Get contacts ...\n")
           dummy_contacts_obj = contacts.get_contacts_obj(pdb_test, dummy_x3dna_obj)
           # Get triads object #
           if verbose: sys.stdout.write("\t\t-- Get triads ...\n")
           if standard:
              # Get triads from standard pdb #
              if verbose: sys.stdout.write("\t\t\t-- Standard triads ...\n")
              triads_obj = triads.get_triads_obj(pdb_test,dummy_dssp_obj,dummy_x3dna_obj,dummy_contacts_obj)
           else:
              # Get triads, pdb etc from threading files #
              if verbose: sys.stdout.write("\t\t\t-- Non-Standard triads using threads...\n")
              thread_objs = pdb2thread.pdb2thread(pdb_test,pdb_name,dummy_x3dna_obj,dummy_dir,test_sequence,verbose)
              triads_obj = triads.Triads()
              for thread_obj in thread_objs:
                  if verbose: sys.stdout.write("\t\t\t-- Add triads from template %s chain %s ...\n"%(thread_obj._pdb_name.lower(),thread_obj._pdb_chain))
                  thread_dummy_dir = os.path.join(dummy_dir,thread_obj._pdb_name.lower()+"_"+thread_obj._pdb_chain)
                  if not os.path.exists(thread_dummy_dir): os.makedirs(thread_dummy_dir)
                  thread_triads_obj, thread_pdb_obj, thread_x3dna_obj = threading_to_triads.get_triads(thread_obj,dummy_dir,pdb_dir,verbose,thread_dummy_dir)
                  for triad_obj in thread_triads_obj.get_triads():
                      triads_obj.add_triad(triad_obj)
           if verbose: sys.stdout.write("\t\t-- Calculate scores ...\n")
           try:
               scr.calculate_energies_and_binding(triads_obj, dummy_x3dna_obj, msa_objs, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir,methylation)
           except Exception as e:
               print("Failed to get energy scores %s"%(e))
               raise Exception("Failed to get energy scores %s\n"%(e))
           test_score = float(scr.get_score(normal=True,potential=split_potential))
           if verbose: sys.stdout.write("\t\t-- Check score in [ %d - %d ] %f ...\n"%(i,i+size+1,test_score))
           if test_score > max_score:
              max_score      = test_score
              start_position = i
              bind_sequence  = test_sequence
              found_sequence = temp_sequence
       # Get the TEMPLATE and Threads
       if verbose: sys.stdout.write("\t-- Selected [ %d - %d ] %f ...\n"%(start_position,start_position+size+1,max_score))
       if verbose: sys.stdout.write("\t-- Selected sequence %s \n"%(bind_sequence))
       try:
          pdb_obj_new = model_dna.get_dna_model_pdb_obj(pdb_file,found_sequence, x3dna_obj, interface_obj, interface_range=None, dummy_dir=dummy_dir)
       except Exception as e:
          print(("Error: %s"%e))
          raise Exception("Failed to model DNA ",e)
       try:
          thread_objs = pdb2thread.pdb2thread(pdb_obj_new,pdb_name,x3dna_obj,input_dir,bind_sequence,verbose)
       except Exception as e:
          print(("Error: %s"%e))
          raise Exception("Failed to get threads ",e)

    return  thread_objs

     
def parse_options():
    """
    Parse the command line options for best-binding detection and threading.

    How to run:
        python get_best_binding.py -i PDB_FILE --pwm MEME_FILE --seq FASTA_FILE
            --pbm PBM_DIR --pdb PDB_DIR -o OUTPUT_DIR
            [--dna BINDING_SEQ --delta N --pval PVALUE --standard -v]
            [pwm options] [statistical potential options]

    Example:
        python get_best_binding.py -i complex.pdb --pwm motif.meme --seq dna.fa
            --pbm pbm_data --pdb pdb_data -o best_binding -v

    The parser configures:
        - Input structure, PWM motif, and DNA-sequence scan target.
        - PBM/PDB resources and output/logging controls.
        - Binding-window and wobble/clean behavior.
        - PWM refinement and statistical-potential configuration.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options describing the structure, PWM, DNA
        sequence, output paths, and scoring/refinement settings.
    """

    parser = optparse.OptionParser("get_best_binding -i PDB_FILE --pwm MEME_FILE -seq FASTA_FILE -o OUTPUT_NAME --pbm PBM_FOLDER --pdb PDB_FOLDER [--dna SPECIFIC_BINDING --dummy DUMMY_DIR --standard --delta DELTA --pval P-VALUE][PWM OPTIONS][STATISTICAL POTENTIAL OPTIONS]")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR") 
    parser.add_option('-i', '--input_file', dest = 'input_file', action = 'store', metavar="PDB_FILE",
                        help = 'PDB file of TF with DNA')
    parser.add_option('--seq', dest = 'sequence_file', action = 'store', default=None,metavar="FASTA_FILE",
                        help = 'Input file in FASTA format with DNA sequence to scan')
    parser.add_option("--pwm", action="store",  dest="pwm",  default=None,metavar="MEME_FILE",
                        help="PWM file in MEME format" )
    parser.add_option("--dna", action="store",  dest="binding_site",  default=None,metavar="{string}",
                        help="DNA binding site sequence that must be found by the PWM" )
    parser.add_option('--pval', dest = 'significance', action = 'store', default = 1.0,
                        help = 'Threshold of p-value significance (default is 1.0 to use all)')
    parser.add_option('--delta', dest = 'delta', action = 'store', default = 2,
                        help = 'Increase of the sequence interval around binding to find the complete size of the binding interface in DNA (default is +2 at both ends)')
    parser.add_option("--pdb", action="store",  dest="pdb_dir",  default=None, metavar="FOLDER",
                        help="PDB directory" )
    parser.add_option("--pbm", action="store",  dest="pbm_dir",  default=None, metavar="FOLDER",
                        help="PBM directory" )
    parser.add_option('--standard', default=False, dest = 'standard', action = 'store_true',metavar="{boolean}",
                        help = 'Flag to use standard unmethylated nucleotides for the best binding region. The output will have the right sequence (default is False)')
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of PDB files that have failed and have been completed")
    parser.add_option('-o', '--output_dir', dest = 'output_dir', action = 'store', default = None, metavar='{string}',
                        help = 'Output folder for PDB and thread files')
    parser.add_option('-v', '--verbose', dest = 'verbose', action = 'store_true',metavar="{boolean}",
                        help = 'Flag for verbose mode (default is False)')
    parser.add_option('-c', '--clean', dest = 'clean', action = 'store_true',metavar="{boolean}",
                        help = 'Flag to clean the PDB files and renumber staring at 1 (default is False)')
    parser.add_option('-w', '--wobble', dest = 'wobble', action = 'store_true',metavar="{boolean}",
                        help = 'Flag to wobble around binding and check variations of +/- delta around best binding (default is False, unless the best binding differs from the interface)')

    pwm_options = optparse.OptionGroup(parser,"Data to calculate redo the PWM")

    pwm_options.add_option("--meme", default=False, action="store_true", dest="meme", help="Use 'uniprobe2meme' to calculate the PWM matrix for 'FIMO' (default = False)")
    pwm_options.add_option("-r","--reset", default=False, action="store_true", dest="reset", help="Clean the sequences of the original MSA and reset them by a random selection in accordance with the PWM (default = False)")
    pwm_options.add_option("--refine",default=0, action="store", type="int", dest="refine", help="Level to refine the MSA and PWM scoring DNA binding sequences with full nucleotides length and without dummy N (default=0, 0 no-refinement keeps dummy N  nucleotides, 1 refine and trim, 2 refine with rescaling and cut-off)", metavar="{integer}")


    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=True, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = True)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-t", action="store", type="float", dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    group.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' from configuration and it uses default radius to load potentials automatically, otherwise it uses radius for both", metavar="{string}")

    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment_restrict", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. Example: '45_A-48_A;50_A-51_A' (Default is None it applies to all amino-acids)")
    parser.add_option("--binding", default=None, action="store", type="string", dest="binding_restrict", help="Binding site of DNA to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d of the forward chain (first in PDB). (Default is None it applies to all nucleotides)")

    group.add_option("--methylation", default=False, action="store_true", dest="methylation", help="Use methylated cytosine specificities of binding/non-binding to calculate PWM motifs and binding sites (default = False)")

    parser.add_option_group(pwm_options)
    parser.add_option_group(group)


    (options, args) = parser.parse_args()

    if options.input_file is None or options.pdb_dir is None  or options.pbm_dir is None or  options.sequence_file is None  :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def has_methylation(sequence):
    """
    Check whether a DNA sequence contains methylation-specific symbols.

    Args:
        sequence (str): DNA sequence to inspect.

    Returns:
        bool: True when the sequence contains any of ``X``, ``J``, ``O``, or ``Q``.
    """
    for x in list("XJOQ"):
        if x in sequence:
            return True
    return False

def main():
    """
    Run the command-line workflow for selecting the best PWM-supported binding site.

    Workflow:
        1. Read the input structure, PWM, and DNA sequence.
        2. Use FIMO and optional rescoring to choose the best binding region.
        3. Generate threaded TF-DNA models for the chosen region.

    Returns:
        None. Output models/files are written into the requested output directory.
    """
    # Arguments & Options #
    options = parse_options()
    if not os.path.exists(options.dummy_dir): 
        os.makedirs(options.dummy_dir)
    dummy_dir    = os.path.abspath(os.path.join(options.dummy_dir, str(os.getpid())))
    if not os.path.exists(dummy_dir): 
        os.makedirs(dummy_dir)
    verbose      = options.verbose
    pwm_file     = options.pwm
    binding      = options.binding_site
    if binding is not None: 
       binding   = binding.upper()
       reverse   = functions.reverse_dna(binding)
    else:
       reverse   = None
    dna_file     = options.sequence_file
    input_file   = os.path.abspath(options.input_file)
    output_dir   = os.path.abspath(options.output_dir)
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    pbm_dir      = os.path.abspath(options.pbm_dir)
    pdb_dir      = os.path.abspath(options.pdb_dir)
    p_val_max    = float(options.significance)
    delta        = int(options.delta)
    standard     = options.standard
    input_dir    = os.path.dirname(input_file)
    input_pdb_file       = input_file
    input_threading_file = None
    threading            = False
    potential_file       = options.potential_file
    radius               = options.radius
    fragment_restrict    = options.fragment_restrict
    binding_restrict     = options.binding_restrict
    split_potential      = options.split_potential
    taylor_approach      = options.taylor_approach
    score_threshold      = options.score_threshold
    auto_mode            = options.auto_mode
    bins_approach        = options.bins
    pmf                  = options.pmf
    known                = options.known
    pbm_potentials       = options.pbm_potentials
    family_potentials    = options.family_potentials
    methylation          = options.methylation
    meme                 = options.meme
    reset                = options.reset
    refine               = options.refine
    info                 = options.info
    clean                = options.clean
    wobble               = options.wobble


    if info is None:
        info = os.path.join(output_dir,"get_best_binding.execution.log")

    # Get full dna sequence
    if verbose: sys.stdout.write("Get DNA sequence %s ....\n"%dna_file)
    try:
      for header_title, sequence in functions.parse_fasta_file(dna_file):
        if verbose: sys.stdout.write("\t- %s \n"%header_title)
        dna_sequence = sequence
      if verbose: sys.stdout.write("\t- Sequence %s \n"%dna_sequence)
    except:
      log_file=open(info,"a")
      log_file.write("%s\tFAIL\n"%(os.path.basename(input_file)))
      log_file.close()
      exit()



    #Get data for scores, potentials and PWM
    try:
      families = {}
      if options.verbose:sys.stdout.write("Check families...\n")
      if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
          sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
      else:
         for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
           if line.startswith("#"): continue
           pdb_chain, family = line.split(";")
           families[pdb_chain] = family
    except:
      log_file=open(info,"a")
      log_file.write("%s\tFAIL\n"%(os.path.basename(input_file)))
      log_file.close()
      exit()


    # Get the PWM
    if pwm_file is None:
       pwm_name = input_file.rstrip(".pdb")
       pwm_file = input_file.rstrip(".pdb")+".meme"
    else:
       pwm_name = input_file.rstrip(".txt").rstrip("meme")
    if not os.path.exists(pwm_file):
       input_pdb_file = input_file
       #Select the format of PWM when methylation
       if meme: methylation_pwm = False
       else:    methylation_pwm = methylation

       try:
              if verbose: sys.stdout.write("\t-- Getting PWM %s ...\n"%(pwm_file))
              PWM.get_single_pwm(input_pdb_file,input_threading_file,threading,pwm_name,pbm_dir,pdb_dir,families,potential_file, radius, fragment_restrict, binding_restrict, "s3dc_dd",auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins_approach,known,meme,reset,refine,options.dummy_dir,verbose,methylation_pwm)
              if verbose: sys.stdout.write("\t-- Done PWM %s ...\n"%pwm_file)
       except:
              print(("%s\tFAIL\n"%(pwm_file)))
              log_file=open(info,"a")
              log_file.write("%s\tFAIL\n"%(os.path.basename(input_file)))
              log_file.close()
              exit()
    #Get Threads with best binding sequence
    try:
         if verbose: sys.stdout.write("\t-- Getting Threads ...\n")
         thread_objs= get_threads_with_best_binding(binding,delta,standard,pwm_file,dna_sequence,input_pdb_file,input_dir,pdb_dir,pbm_dir,radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins_approach, known,p_val_max,verbose,dummy_dir,fragment_restrict, binding_restrict,clean,wobble,methylation)
    except:
         log_file=open(info,"a")
         log_file.write("%s\tFAIL\n"%(os.path.basename(input_file)))
         log_file.close()
         exit()

    # Write PDB templates and threads with the best DNA sequence
    if verbose: sys.stdout.write("Write threading files\n")
    try:
      threading_pdb_name = thread_objs[0].get_pdb_name()
      output_complex = os.path.join(output_dir,threading_pdb_name+"_")
      thread_complex = []
      n=0
      for thr_obj in thread_objs:
        n=n+1
        threading_pdb_name  = thr_obj.get_pdb_name()
        threading_pdb_chain = thr_obj.get_pdb_chain()
        template_pdb        = os.path.join(input_dir,threading_pdb_name+"_"+threading_pdb_chain+".pdb")
        if verbose: sys.stdout.write("\t-- Write template file %s\n"%(os.path.basename(template_pdb)))
        pdb_obj = PDB(template_pdb)
        pdb_obj.clean()
        pdb_obj.write(os.path.join(output_dir,os.path.basename(template_pdb)))
        #shutil.move(template_pdb,os.path.join(output_dir,os.path.basename(template_pdb)))
        output_file = os.path.join(output_dir,threading_pdb_name+"_"+threading_pdb_chain)
        output_complex = output_complex + threading_pdb_chain
        kmers=thr_obj.get_kmers()
        for kmr in kmers.keys():
            dna_seq = kmr
        output_file = output_file+"_"+dna_seq+".txt"
        if os.path.exists(output_file):
            if verbose: print(("\t-- Reuse threading file %s\n"%(output_file)))
            thr_obj.write(output_file)
        else:
            if verbose: sys.stdout.write("\t-- Write threading file %s\n"%(output_file))
            thr_obj.write(output_file)
        thread_complex.append(output_file)
      log_file=open(info,"a")
      log_file.write("%s\tDONE\n"%(os.path.basename(input_file)))
      log_file.close()
      if n>1:
        output_complex = output_complex+"_"+dna_seq+".txt"
        if verbose: sys.stdout.write("\n\t-- Write thread complex  %s\n"%(output_complex))
        fo=open(output_complex,"w")
        for output_file in thread_complex:
          fo.write("%s\n"%output_file)
        fo.close()
    except:
      log_file=open(info,"a")
      log_file.write("%s\tFAIL\n"%(os.path.basename(input_file)))
      log_file.close()


    # Clean files #
    try:
       shutil.rmtree(dummy_dir)
    except:
       print("Remaining DUMMY files, not cleanded")

    if verbose: print("Done")


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
