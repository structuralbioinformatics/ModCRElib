"""
Merge a main DNA PWM with one or more additional PWMs using TOMTOM alignment.

This module contains the core implementation used by the PWM-merging command
line scripts. It can align two motifs, optionally refine the structure-anchored
region with ModCRE statistical potentials, reconstruct full merged motifs, and
write the resulting representations to disk.
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
    Parse the command line options for PWM-merging workflows.

    The main PWM acts as the TOMTOM anchor motif, while the PWM list supplies
    the motifs that will be aligned and merged against it.

    Returns:
        optparse.Values: Parsed CLI options describing input motifs, structural
        model, output mode, queueing options, and statistical-potential
        settings.
    """

    parser = optparse.OptionParser("python merge_pwm.py -i input_main_pwm -j input_add_pwm -w weight [-o output_file --dummy dummy_dir -w weight --pbm pbm_dir --pdb pdb_dir  --refine refine ]  ")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", default=None, type="string", dest="main_pwm", help="MAIN PWM associated with a model structure ", metavar="{filename}")
    parser.add_option("-j", action="store", default=None, type="string", dest="pwm_list", help="List of PWMs (format path/file  weigth) ", metavar="{filename}")
    parser.add_option("--factor",action="store", default=1.0, type="float", dest="factor", help="Factor that multiplies the weight defined in the List of PWMs (option j). Default is 1", metavar="{float}")
    parser.add_option("--model",action="store", default=None, type="string", dest="model", help="Model structure corresponding to the MAIN PWM (it may be in the format of thread)", metavar="{filename}")
    parser.add_option("--sig",action="store", default=0.05,type="float",dest="significance",help="Limit of significance (p-value) required in TOMTOM to merge the PWMs (defaults is 0.05)")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", default=None, help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("--main",default=False, action="store_true", dest="main", help="Generate a PWM for the positions of the MAIN PWM associated with structure (default is False)", metavar="{boolean}")
    parser.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    parser.add_option("--meme", default=False, action="store_true", dest="meme", help="Use 'uniprobe2meme' to calculate the PWM matrix for 'FIMO' (default = False)")
    parser.add_option("--threading", default=False, action="store_true", dest="threading", help="Input file is a threading file of a PDB structure that (MUST!) exist in the PDB folder of ModCRE (default = False)")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment_restrict", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. Example: '45_A-48_A;50_A-51_A' (Default is None it applies to all amino-acids)")
    parser.add_option("--binding", default=None, action="store", type="string", dest="binding_restrict", help="Binding site of DNA to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d of the forward chain (first in PDB). (Default is None it applies to all nucleotides)")
    parser.add_option("-o", action="store", default="output_pwm", type="string", dest="output_file", help="Output pwm (default = 'output_pwm')", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of PWMs that have failed and have been completed")
    parser.add_option("--reuse",default=False, action="store_true", dest="reuse", help="Reuse the information files. If the flag is used then profiles that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")
    parser.add_option("--parallel", default=False, action="store_true", dest="parallel", help="Run in parallel if the input is a directory (default = False)")
    parser.add_option("--complete",default=1.00, action="store", type="float", dest="complete", help="Ratio of completness over the total number of profiles top be done(default= 0.95). This is useful in the server to stop the profiler when the time of execution exceeds more than 48 hours ", metavar="RATIO")
  

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

    if options.main_pwm is None or options.pwm_list is None or options.model is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Functions   #
#-------------#

def merge_two_pwm(main_pwm, second_pwm, weight, refine, significance, input_threading_file, input_pdb_file, pdb_dir, pbm_dir, threading, fragment_restrict, binding_restrict, output_file, families, radius, potential_file, split_potential,  auto_mode,  family_potentials,  pbm_potentials,  score_threshold,  taylor_approach,  pmf, bins_approach, known,   methylation, dummy_dir, verbose):
    """
    Align and merge one secondary PWM into a main PWM.

    The two motifs are first compared with TOMTOM to obtain the best offset,
    strand, overlap, and significance. If the hit passes the p-value cutoff,
    both PWMs are merged into a full-length motif. Optionally, the region tied
    to the structural model is refined with the ModCRE scoring machinery.

    Args:
        main_pwm (str): Path to the main PWM in MEME format.
        second_pwm (str): Path to the second PWM in MEME format.
        weight (float): Relative weight assigned to the second PWM.
        refine (int): Refinement level for the structure-anchored region.
        significance (float): Maximum TOMTOM p-value allowed for merging.
        input_threading_file, input_pdb_file (str | None): Structural input
            describing the main PWM context.
        pdb_dir, pbm_dir (str | None): ModCRE data directories.
        threading (bool): Whether the structural input is a threading file.
        fragment_restrict, binding_restrict: Optional structural restrictions.
        output_file (str): Output prefix used for temporary TOMTOM output and
            derived motif names.
        families (dict): Mapping used for family-specific potentials.
        radius, potential_file, split_potential, auto_mode, family_potentials,
        pbm_potentials, score_threshold, taylor_approach, pmf, bins_approach,
        known, methylation, dummy_dir, verbose: Additional ModCRE workflow
        parameters forwarded to refinement helpers.

    Returns:
        tuple: ``(result_full_pwm, msa_obj, over5, overlap)`` where
        ``result_full_pwm`` is the full merged motif, ``msa_obj`` is the main
        structure-anchored section, ``over5`` is the 5' overlap start in the
        merged motif, and ``overlap`` is the TOMTOM overlap length.
        If TOMTOM significance is not sufficient, returns ``(None, None, offset, overlap)``.
    """

    # Compare the PWMs
    tomtom_results=tomtom.get_tomtom_obj(main_pwm,second_pwm,dummy_dir)
    tomtom_obj=tomtom_results.get_hits()[0]
    offset  = tomtom_obj.get_offset()
    strand  = tomtom_obj.get_strand()
    overlap = tomtom_obj.get_overlap()
    p_value = tomtom_obj.get_p_value()

    if verbose:sys.stdout.write("EXECUTE TOTOM %s\n"%(output_file+".tomtom"))
    tomtom_results.write(output_file+".tomtom")

    if p_value > significance: return None,None,offset,overlap

    if strand == "-":
       work_main_pwm = PWM.nMSA(main_pwm,"main","meme").get_complementary()
    else:
       work_main_pwm = PWM.nMSA(main_pwm,"main","meme")


    work_offset = abs(offset)
    if offset>0:
       reference_pwm = PWM.nMSA(second_pwm,"second","meme")
       target_pwm = work_main_pwm
       reference_weight = weight
       target_weight = 1.0
       prime5 = 0
       if work_offset + reference_pwm.get_binding_site_length() <= target_pwm.get_binding_site_length(): prime3 = 0
       if work_offset + reference_pwm.get_binding_site_length() > target_pwm.get_binding_site_length():  prime3 = work_offset + reference_pwm.get_binding_site_length() - target_pwm.get_binding_site_length() 
    else:
       reference_pwm = work_main_pwm
       target_pwm = PWM.nMSA(second_pwm,"second","meme")
       target_weight = weight
       reference_weight =  1.0
       prime5 = work_offset
       if work_offset + reference_pwm.get_binding_site_length() >= target_pwm.get_binding_site_length(): prime3 = 0
       if work_offset + reference_pwm.get_binding_site_length() < target_pwm.get_binding_site_length():  prime3 = target_pwm.get_binding_site_length() - work_offset - reference_pwm.get_binding_site_length() 
    lenght_main = work_main_pwm.get_binding_site_length()
    length_reference = reference_pwm.get_binding_site_length()
    length_target = target_pwm.get_binding_site_length()
    length_alignment = prime5 + lenght_main + prime3
    over5=work_offset
    over3=length_alignment-work_offset-overlap

    #Merge the PWMs
    merged_name =  "new_"+work_main_pwm.get_motif()
    if verbose: sys.stdout.write("Reference query sequence %s\n"%(reference_pwm.get_main_sequence()))
    if verbose: sys.stdout.write("Target fixed sequence    %s\n"%(target_pwm.get_main_sequence()))

    merged_pwm  = reference_pwm.merge(target_pwm, work_offset, reference_weight, target_weight, merged_name)
    if verbose: sys.stdout.write("Merged sequence          %s\n"%(merged_pwm.get_main_sequence()))

    full_length = merged_pwm.get_binding_site_length()
    if full_length != length_alignment:
        sys.stdout.write("ERROR\n")
        sys.stdout.write(("ref_length %d tgt_len %d offset %d overlap %d prime5 %d prime3 %d ali %d full %d\n"%(length_reference,length_target,offset,overlap,prime5,prime3, length_alignment,full_length)))
        exit()

    if strand == "-":
        full_pwm = merged_pwm.get_complementary()
        cp5 = prime3
        cp3 = prime5
        prime5 = cp3
        prime3 = cp5
        co5 = over3
        co3 = over5
        over5 = co3
        over3 = co5
    else:
        full_pwm = merged_pwm


    if prime5>0: msa_prime5 = full_pwm.section(0,prime5)
#    msa_obj = full_pwm.section(prime5,(full_length-prime3)+1)
#    if prime3>0: msa_prime3 = full_pwm.section((full_length-prime3)+1,full_length)
    msa_obj = full_pwm.section(prime5,(full_length-prime3))
    if prime3>0: msa_prime3 = full_pwm.section((full_length-prime3),full_length)

    if verbose:
       if prime5>0:sys.stdout.write("5p %s %d\n"%(msa_prime5.get_main_sequence(),prime5))
       if prime3>0:sys.stdout.write("3p %s %d\n"%(msa_prime3.get_main_sequence(),prime3))
       sys.stdout.write("MSA_OBJ %s %d %d\n"%(msa_obj.get_main_sequence(),len(msa_obj.get_main_sequence()),full_length-prime3-prime3))

    # Refine msa_obj if requested
    if refine>0:
      if threading:
         #Get triads, pdb etc from threading file #
         if  verbose:sys.stdout.write("\t--Get thread %s ...\n"% input_threading_file)
         triads_obj, pdb_obj, x3dna_obj = threading_to_triads.threading_triads(threading_file=input_threading_file, pdb_dir= pdb_dir)

      else:

        try:
         # Get PDB object #
         if  verbose:sys.stdout.write("\t--Get protein %s ...\n"% input_pdb_file)
         pdb_obj = PDB( input_pdb_file)

         # Get DSSP object #
         if  verbose:sys.stdout.write("\t\t-- calculate secondary structure ...\n")
         dssp_obj = dssp.get_dssp_obj(os.path.abspath( input_pdb_file))

         # Get X3DNA object #
         if  verbose:sys.stdout.write("\t--Get DNA %s ...\n"% input_pdb_file)
         x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath( input_pdb_file),dummy_dir=dummy_dir)
         if len(list(x3dna_obj.get_dinucleotides().keys()) )<1:
            sys.stdout.write("Missing DNA ...\n")
        
         # Get contacts object #
         if  verbose:sys.stdout.write("\t\t-- calculate contacts ...\n")
         contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
         if len(contacts_obj.get_contacts())<1:
            sys.stdout.write("Missing Protein-DNA contacts ...\n")

         # Get triads object #
         if  verbose:sys.stdout.write("\t\t-- calculate protein-dna pairs ...\n")
         triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)

        except Exception as e:
         raise Exception("Failed to get protein features")

      # Load statistical potential #
      if  verbose:sys.stdout.write("\t--Load potentials ...\n")
      try:
        potentials, thresholds , radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj,  pdb_dir, pbm_dir, families, radius, potential_file, split_potential,  auto_mode,  family_potentials,  pbm_potentials,  score_threshold,  taylor_approach,  pmf, bins_approach, known, None,   dummy_dir, verbose)
      except Exception as e:
         raise Exception("Failed to get potentials %s"%(e))

      #  Apply refinement
      if  verbose:sys.stdout.write("\t--Refine MSA object ...\n")
      try:
          msa_obj = PWM.refine_msa_obj(refine, msa_obj, triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, split_potential, thresholds, dummy_dir, verbose, methylation)
      except Exception as e:  
          raise Exception("Failed to get refined MSA %s"%(e))


    # construct the final full PWM
    if prime5>0: 
        result_full_pwm = msa_prime5
        result_full_pwm = result_full_pwm + msa_obj
    else:
        result_full_pwm = msa_obj
    if prime3>0:
        result_full_pwm = result_full_pwm + msa_prime3

    sys.stdout.flush()
    result_full_pwm.set_motif("full_"+os.path.basename(output_file))
    msa_obj.set_motif("main_"+os.path.basename(output_file))

    return result_full_pwm, msa_obj, over5, overlap




def write_pwm(msa_obj,folder,output_name,meme,dummy_dir,verbose):
    """
    Write a merged motif to the standard ModCRE output formats.

    Args:
        msa_obj (PWM.nMSA): Motif object to write.
        folder (str): Output directory.
        output_name (str): Output prefix without extension.
        meme (bool): If True, also export the MEME-suite ``.meme.s`` variant.
        dummy_dir (str): Temporary directory for helper files.
        verbose (bool): If True, print progress messages.

    Returns:
        None. Writes ``.pwm``, ``.meme``, ``.msa``, optional ``.meme.s``, and
        logo PNG files.
    """

    pwm_file = os.path.join(folder,output_name+".pwm")
    meme_file = os.path.join(folder,output_name+".meme.s")
    pwm_meme = os.path.join(folder,output_name+".meme")
    msa_file = os.path.join(folder,output_name+".msa")
    logo_file= os.path.join(folder,output_name+".logo")

    # Write PWM #
    if  verbose:sys.stdout.write("\t--Write PWM ...\n")

    if not os.path.exists(pwm_file):
     if  verbose:sys.stdout.write("\t\t--PWM...\n")
     try:
       msa_obj.write(pwm_file, option="pwm")
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)
    if not os.path.exists(pwm_meme):
     if  verbose:sys.stdout.write("\t\t--PWM in MEME format...\n")
     try:
       msa_obj.write(pwm_meme, option="meme")
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)
    if not os.path.exists(msa_file):
     if  verbose:sys.stdout.write("\t\t--MSA...\n")
     try:
       msa_obj.write(msa_file, option="msa")
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)

    if not os.path.exists(meme_file) and  meme:
     if  verbose:sys.stdout.write("\t\t--PWM by MEME...\n")
     try:
       PWM.write_pwm_by_meme(msa_obj,meme_file, dummy_dir)
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)

    if not os.path.exists(logo_file+".fwd.png") or  not os.path.exists(logo_file+".rev.png"):
     if  verbose:sys.stdout.write("\t\t--Logos...\n")
     try:
       PWM.write_logo(msa_obj,logo_file, dummy_dir)
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)

    sys.stdout.flush()

 

def get_single_pwm(options, input_file=None, output_file=None):
    """
    Merge one main PWM against all PWMs listed in the provided manifest.

    For each candidate PWM, the function attempts a TOMTOM-guided merge through
    :func:`merge_two_pwm`. The significant merged motifs are then combined into
    a consensus output PWM and written to disk.

    Args:
        options (optparse.Values): Parsed CLI options.
        input_file (str | None): Optional path overriding ``options.main_pwm``.
        output_file (str | None): Optional output prefix overriding
            ``options.output_file``.

    Returns:
        None. Writes the combined merged motif to disk when enough significant
        merges are found.
    """

    # Initialize #

    # Dummy/temporary folder
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    dummy_dir=os.path.abspath(options.dummy_dir)

    # input structures and folders of ModCRE
    pbm_dir=None
    if options.pbm_dir is not None: pbm_dir = os.path.abspath(options.pbm_dir)
    pdb_dir=None
    if options.pdb_dir is not None: pdb_dir = os.path.abspath(options.pdb_dir)
    threading = options.threading
    if input_file is None:
      output_file = options.output_file
      if options.threading:
        input_threading_file=options.model
        input_pdb_file=None
        if not input_threading_file.startswith("/"): input_threading_file=os.path.abspath(options.model)
      else:
        input_pdb_file=options.model
        input_threading_file=None
        if not input_pdb_file.startswith("/"): input_pdb_file=os.path.abspath(options.model)
      #input PWMs
      main_pwm=os.path.abspath(options.main_pwm)
    else:
      main_pwm=input_file
      if options.threading:
        input_threading_file=os.path.join(os.path.dirname(main_pwm),".".join([x for x in  os.path.basename(main_pwm).split(".")[:-1]])+".txt")
        input_pdb_file=None
      else:
        input_pdb_file=os.path.join(os.path.dirname(main_pwm),".".join([x for x in os.path.basename(main_pwm).split(".")[:-1]])+".pdb")
        input_threading_file=None


    #input PWMs
    pwm_list=os.path.abspath(options.pwm_list)
    meme = options.meme
    significance = float(options.significance)
    factor = float(options.factor)

    pwms=[]
    pwm_files = open(pwm_list,"r")
    for line in pwm_files:
        data = line.strip().split()
        if len(data)>1: pwms.append((data[0],factor*float(data[1])))
    pwm_files.close()


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
    output_pwm = os.path.basename(output_file)
    output_dir = os.path.dirname(output_file)
    output_name= output_pwm + ".combined"
    verbose = options.verbose
 
    # Do list of merged pairs
    merged_pwms=[]
    number=0
    for second_pwm,weight in pwms: 
        number=number+1
        #Compare a pwm with main pwm
        output_file = os.path.join(output_dir,output_pwm+"_"+os.path.basename(second_pwm))
        if verbose: sys.stdout.write("Compare %s vs %s \n"%(os.path.basename(main_pwm),os.path.basename(second_pwm)))
        try:
          result_full_pwm, msa_obj, over5, overlap = merge_two_pwm(main_pwm,second_pwm, weight, refine, significance, input_threading_file, input_pdb_file, pdb_dir, pbm_dir, threading, fragment_restrict, binding_restrict, output_file, families, radius, potential_file, split_potential,  auto_mode,  family_potentials,  pbm_potentials,  score_threshold,  taylor_approach,  pmf, bins_approach, known,   methylation, dummy_dir, verbose)
          if msa_obj is None: 
            if verbose: sys.stdout.write("Skip %s\n"%(second_pwm))
            continue
          msa_obj.set_motif(str(number))
          merged_pwms.append(msa_obj)
          print("Number of PWMS passing the threshold %d"%len(merged_pwms))
        except:
          print("Skip PWM %s"%second_pwm)
          continue

    # Combine merged pairs
    if len(merged_pwms) >1:
       main_pwm_obj = merged_pwms[0]
       rest_pwm_obj = merged_pwms[1:]
       if len(os.path.basename(pwm_list).split("."))>1:
          combined_name=".".join(os.path.basename(pwm_list).split(".")[:-1])+":"+"_".join([x.get_motif() for x in merged_pwms ])
       else:
          combined_name=os.path.basename(pwm_list).split(".")[0]+":"+"_".join([x.get_motif() for x in merged_pwms ])
       result_pwm   = main_pwm_obj.combine(rest_pwm_obj)
       result_pwm.set_motif(output_name+"_"+combined_name)
       if verbose: sys.stdout.write("Merged %s sequence:\t %s\n"%(main_pwm_obj.get_motif(),main_pwm_obj.get_main_sequence()))
       for compared_pwm in rest_pwm_obj:
           if verbose: sys.stdout.write("Merged %s sequence:\t %s\n"%(compared_pwm.get_motif(),compared_pwm.get_main_sequence()))
       if verbose: sys.stdout.write("Combined sequence:\t %s\n"%(result_pwm.get_main_sequence()))
       # Write outputs
       write_pwm(result_pwm,output_dir,output_name+"_"+combined_name,meme,dummy_dir,verbose)
    elif len(merged_pwms)==1:
        result_pwm   = merged_pwms[0]
        if len(os.path.basename(pwm_list).split("."))>1:
           combined_name=".".join(os.path.basename(pwm_list).split(".")[:-1])+":"+result_pwm.get_motif()
        else:
           combined_name=os.path.basename(pwm_list).split(".")[0]+":"+result_pwm.get_motif()
        result_pwm.set_motif(output_name+"_"+combined_name)
        if verbose: sys.stdout.write("Combined sequence:\t %s\n"%(result_pwm.get_main_sequence()))
        write_pwm(result_pwm,output_dir,output_name+"_"+combined_name,meme,dummy_dir,verbose)
    else:
        sys.stdout.write("Not enough files to be merged\n")

def main():
    """
    Run the command-line PWM-merging workflow.

    The script supports both single-run mode and batch/parallel processing over
    directories of models and PWMs. In both cases it delegates the actual
    motif-merging logic to :func:`get_single_pwm`.

    Returns:
        None. Output motifs and execution logs are written to disk.
    """
    # Arguments & Options #
    options = parse_options()
    if os.path.isdir(options.main_pwm) and os.path.isdir(options.model):
       
       if options.verbose: print("Running over all PWM files, PDB or THREAD must have the same name with suffix 'pdb or 'txt'")
       if not os.path.exists(options.output_file): os.makedirs(options.output_file)
       complete=float(options.complete)
       output_folder  = os.path.abspath(options.output_file)
       input_folder   = os.path.abspath(options.model)
       if options.threading:
         pdb_files = [x for x in os.listdir(input_folder) if x.endswith(".txt")]
       else:
         pdb_files = [x for x in os.listdir(input_folder) if x.endswith(".pdb")]
       name_of_pwms=set()
       for pdb_file in pdb_files:
         if options.threading: output_file = os.path.join( output_folder,pdb_file.rstrip(".txt") )
         else:                 output_file = os.path.join( output_folder,pdb_file.rstrip(".pdb") )
         name_of_pwms.add(os.path.basename(output_file)+".meme")
       # Iterate untill all protein profiles are done
       submitted=set()
       n_done=0
       if options.verbose: print("Start iteration to check and run profiles")
       info_file=options.info
       if info_file is None: info_file=os.path.join(output_folder,"pwm_list_execution.log")
       if options.reuse and functions.fileExist(info_file):
         if options.verbose: print(("Reuse previous information of runs from file %s"%info_file))
         for pdb_file in pdb_files:
           if options.threading:
              input_pdb_file = None
              input_threading_file = os.path.join(input_folder,pdb_file)
              output_file    = os.path.join( output_folder,pdb_file.rstrip(".txt") )
           else:
              input_pdb_file = os.path.join(input_folder,pdb_file)
              input_threading_file = None
              output_file    = os.path.join( output_folder,pdb_file.rstrip(".pdb") )
         if functions.fileExist(output_file+".meme"):
              if options.verbose:print(("\t-- Found PWM %s"%os.path.basename(output_file)))
              if options.verbose:print("\t\t-- rewrite the logos...")
              msa_obj=PWM.nMSA(output_file+".meme",output_file,"meme")
              write_logo(msa_obj,output_file+".logo")
       else:
         if options.verbose: print(("Open to write %s"%info_file))
         log_file = open(info_file,"w")
         log_file.write("#List of PWMs\n")
         log_file.close()
       done=functions.done_jobs(info_file)
       iterate = functions.check_done(done,name_of_pwms)
       maxtime = 3600 * 3
       start_time = d_time = 0
       while( iterate ):
         for pdb_file in pdb_files:
           if options.threading:
              input_pdb_file = None
              input_threading_file = os.path.join(input_folder,pdb_file)
              output_file    = os.path.join( output_folder,pdb_file.rstrip(".txt") )
              main_pwm       = os.path.join(input_folder,pdb_file.rstrip(".txt")+".meme")
           else:
              input_pdb_file = os.path.join(input_folder,pdb_file)
              input_threading_file = None
              output_file    = os.path.join( output_folder,pdb_file.rstrip(".pdb") )
              main_pwm       = os.path.join(input_folder,pdb_file.rstrip(".pdb")+".meme")
           if pdb_file in submitted: continue
           submitted.add(pdb_file)
           if options.verbose:sys.stdout.write("Generate PWM %s ...\n"%(os.path.basename(output_file)))
           if functions.fileExist(output_file+".meme"):
              if options.verbose:print(("\t-- Found PWM %s"%os.path.basename(output_file)))
              if options.reuse:
                 if options.verbose:print("\t\t-- rewrite the logos...")
                 msa_obj=PWM.nMSA(output_file+".meme",output_file,"meme")
                 write_logo(msa_obj,output_file+".logo")
              log_file=open(info_file,"r")
              skip_adding=False
              for line in log_file:
                  if os.path.basename(output_file)+".meme" in line.split(): skip_adding=True
              if skip_adding and options.verbose: print("\t\t-- Already in the information file as DONE")
              submitted.add(pdb_file)
              if skip_adding: continue
              if options.verbose: print("\t\t-- Add in the list of done")
              log_file.close()
              log_file=open(info_file,"a")
              log_file.write("%s\tDONE\n"%(os.path.basename(output_file)+".meme"))
              log_file.flush()
              log_file.close()
              continue
           if options.reuse:
              if options.verbose:print(("\t-- Not found protein profile %s"%os.path.basename(output_file)))
              log_file = open(info_file,"r")
              skip_adding=False
              for line in log_file:
                       if os.path.basename(output_file)+".meme" in line.split(): skip_adding=True
              if skip_adding and options.verbose: print("\t-- Already in the information file but FAILED")
              if skip_adding: submitted.add(pdb_file)
              if skip_adding: continue
              log_file.close()
           if options.parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              program=os.path.join(exe_path,"merge_pwms.py")
              python = os.path.join(config.get("Paths", "python_path"), "python")
              parameters = " -i %s "%main_pwm
              parameters = parameters + " -j %s "%options.pwm_list
              parameters = parameters + " --sig %f "%options.significance
              if options.threading:
                  parameters =   parameters + " --model %s "%input_threading_file
                  parameters =   parameters + " --threading "
              else:
                  parameters =   parameters + " --model %s "%input_pdb_file
              parameters =  parameters + " -o %s "%output_file
              parameters =  parameters + " --pdb=%s "%options.pdb_dir
              parameters =  parameters + " --pbm=%s "%options.pbm_dir
              parameters =  parameters + " --dummy=%s "%options.dummy_dir
              parameters =  parameters + " -s %s "%options.split_potential
              parameters =  parameters + "--info %s "%info_file
              parameters =  parameters + "--factor %f "%float(options.factor)
              if options.potential_file  is not None: parameters = parameters + " --file %s "%options.potential_file
              if options.score_threshold is not None: parameters =  parameters + " -t %s "%options.score_threshold
              if binding_restrict is not None  : parameters = parameters + " --binding %s "%options.binding_restrict
              if fragment_restrict is  not None: parameters = parameters + " --fragment %s "%options.fragment_restrict
              if options.radius > 0       :      parameters = parameters + " --radius %f "%options.radius
              if options.refine > 0       :      parameters = parameters + " --refine %d "%options.refine
              if options.verbose          :      parameters = parameters + " --verbose "
              if options.auto_mode        :      parameters = parameters + " --auto "
              if options.family_potentials:      parameters = parameters + " --family "
              if options.pmf              :      parameters = parameters + " --pmf "
              if options.pbm_potentials   :      parameters = parameters + " -p "
              if options.known            :      parameters = parameters + " --known "
              if options.taylor_approach  :      parameters = parameters + " --taylor "
              if options.computation      :      parameters = parameters + " --bins "
              if options.reset            :      parameters = parameters + " --reset "
              if options.methylation      :      parameters = parameters + " --methylation "
              if options.verbose:print(("\t-- Submmit %s %s"%(program,parameters)))
              functions.submit_command_to_queue("%s %s %s" % (python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
              submitted.add(pdb_file)
           else:
            try:
              if options.verbose: print("Run merging/refine for %s"%(output_file+".meme"))
              get_single_pwm(options,input_file=main_pwm,output_file=output_file+".meme")
              log_file=open(info_file,"a")
              log_file.write("%s\tDONE\n"%(os.path.basename(output_file)+".meme"))
              log_file.flush()
              log_file.close()
            except:
              log_file=open(info_file,"a")
              log_file.write("%s\tFAIL\n"%(os.path.basename(output_file)+".meme"))
              log_file.flush()
              log_file.close()
            submitted.add(pdb_file)

         #Check next iteration, profiles submitted and profiles done
         done=functions.done_jobs(info_file)
         iterate= functions.check_done(done,name_of_pwms) 
         if len(done) > n_done:
           n_done=len(done)
           if n_done>1 and start_time==0: start_time = float(time.time())
           if options.verbose: 
                sys.stdout.write("Number of profiles already done %d\n"%n_done)
                if d_time>0: sys.stdout.write("Time: %f\n"%d_time)
                sys.stdout.write("\t-- Check files done ...\n")
                log_file = open(info_file,"r")
                for line in log_file:
                   print(("\t\t-- %s"%line.strip()))
                log_file.flush()
                log_file.close()
                sys.stdout.write("\t-- Still running protein profiles  %s ...\n"%functions.check_done(done,name_of_pwms))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush() 
         #Check next iteration, if exceeding time and enough profiles stop iteration
         if start_time>0:
          current_time = float(time.time())
          d_time = current_time - start_time
          if float(  len(done) ) / len(name_of_pwms) > complete  and d_time > maxtime: iterate=False
          if d_time > maxtime and iterate and complete < 1.0: 
               complete = complete - 0.001
               maxtime  = maxtime  + 600
               if options.verbose: sys.stdout.write("Time: %f Done: %d (%d) Ratio to end: %f\n"%(d_time,len(done),len(name_of_pwms),complete))


    else:
       try:
          get_single_pwm(options, input_file=None,  output_file=None )
          if options.info is not None:
             log_file=open(options.info,"a")
             log_file.write("%s\tDONE\n"%(os.path.basename(options.output_file).rstrip(".meme")+".meme"))
             log_file.close()
       except Exception as e:
         print(("Failed PWM: %s\n"%e))
         if options.info is not None:
             log_file=open(options.info,"a")
             log_file.write("%s\tFAIL\n"%(os.path.basename(options.output_file).rstrip(".meme")+".meme"))
             log_file.close()

    sys.stdout.write("Done")


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()


