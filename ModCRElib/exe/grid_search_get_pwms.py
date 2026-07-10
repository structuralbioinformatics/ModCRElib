"""
Generate PWMs for grid-search combinations and evaluate TOMTOM consistency.
"""

import os, sys, re
from collections import Counter
import configparser
import itertools
import numpy
import optparse
import subprocess
import copy
import datetime
import math

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

config_path  = os.path.join(scripts_path,"configure")


# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = configparser.ConfigParser()
config_file = os.path.join(config_path, "config.ini")
config.read(config_file)

# Imports my functions #
from ModCRElib.beans import functions

# Imports jbonet's module #
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.protein import dssp,tmalign
from ModCRElib.potential import spotentials
from ModCRElib.builder import pbm
from ModCRElib.profile import tomtom
from ModCRElib.msa import pwm_pbm as PWM


def parse_options():
    """
    Parse CLI options for PWM generation/evaluation in grid-search mode.

    How to run:
        ``python grid_search_get_pwms.py --pbm pbm_dir --pdb pdb_dir -p model.pdb``

    Returns:
        optparse.Values: Namespace with PBM/PDB roots, optional precomputed
        files, threshold selection, and verbosity controls.
    """

    parser = optparse.OptionParser("python %prog --pbm=<pbm_dir> --pdb=<pdb_dir> --pwm=<pwm_dir> [--dummy=<dummy_dir> -o <output_dir> --start=<start_step> --stop=<stop_step> -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (output directory from pbm.py)", metavar="<pbm_dir>")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="<pdb_dir>")
    parser.add_option("-p", action="store", type="string", default=None, dest="pdb_file", help="PDB input file", metavar="<pdb_file>")
    parser.add_option("-d", action="store", type="string", default=None, dest="dssp_file", help="DSSP input file", metavar="<dssp_file>")
    parser.add_option("-x", action="store", type="string", default=None, dest="x3dna_file", help="X3DNA input file", metavar="<x3dna_file>")
    parser.add_option("-c", action="store", type="string", default=None, dest="contacts_file", help="contacts input file", metavar="<contacts_file>")
    parser.add_option("-t", action="store", type="string", default=None, dest="triads_file", help="triads input file", metavar="<triads_file>")
    parser.add_option("-r","--reduce_thresholds", default=False, action="store_true", dest="reduce_thresholds" , help="Reduce the number of thresholds to speed up the TOMTOM comparisons (default = False)") 
    parser.add_option("-s", default=None, action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    parser.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    


    (options, args) = parser.parse_args()

    if options.pbm_dir is None or options.pdb_dir is None or options.pdb_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


if __name__ == "__main__":

    # Step 1) Parse options and initialize output subdirectories.
    options = parse_options()
    if options.reduce_thresholds:
           selected_thresholds=[70,74,78,82,86,90,94,98]
    else:
           selected_thresholds=list(range(70, 101, 2))
    pdb_dir = os.path.join(options.pdb_dir)
    pbm_dir = os.path.join(options.pbm_dir)
    pdb_file= os.path.abspath(options.pdb_file)
    if options.verbose: sys.stdout.write("\t-- read PDB %s ...\n"%pdb_file)
    pdb_obj = PDB(pdb_file)
    tf_id   = os.path.basename(pdb_file).split("_mL_")[0]
    motif   = os.path.basename(pdb_file).split("_mL_")[1]
    family  = os.path.basename(pdb_file).split("_mL_")[2]
    template= os.path.basename(pdb_file).split("_")[-3] + "_" + os.path.basename(pdb_file).split("_")[-2]
    # Define the split potential that will be used, Es3dc_dd by default #
    pmf=options.pmf
    split_potential = options.split_potential
    if split_potential is None: 
       split_potential = config.get("Parameters", "split_potential")
       if split_potential is None: 
          split_potential = "s3dc_dd"
    if options.verbose: sys.stdout.write("\t-- potential %s ...\n"%split_potential)
    # Get dssp object #
    if options.dssp_file is None:
        if options.verbose: sys.stdout.write("\t-- caclulate dssp....\n")
        dssp_obj = dssp.get_dssp_obj(pdb_file)
    else:
        if options.verbose: sys.stdout.write("\t-- Use DSSP %s ...\n"%(os.path.basename(options.dssp_file)))
        dssp_obj = dssp.DSSP(os.path.abspath(options.dssp_file))
    # Get x3dna object #
    if options.x3dna_file is None:
        if options.verbose: sys.stdout.write("\t-- caclulate x3dna....\n")
        x3dna_obj = x3dna.get_x3dna_obj(pdb_file)
    else:
        if options.verbose: sys.stdout.write("\t-- Use X3DNA %s ...\n"%(os.path.basename(options.x3dna_file)))
        x3dna_obj = x3dna.X3DNA(os.path.abspath(options.x3dna_file))
    # Get contacts object #
    if options.contacts_file is None:
        if options.verbose: sys.stdout.write("\t-- calculate contacts....\n")
        contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
    else:
        if options.verbose: sys.stdout.write("\t-- Use contacts %s ...\n"%(os.path.basename(options.contacts_file)))
        contacts_obj = contacts.Contacts(os.path.abspath(options.contacts_file))
    # Get triads object #
    if options.triads_file is None:
        if options.verbose: sys.stdout.write("\t-- calculate triads....\n")
        triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)
    else:
        if options.verbose: sys.stdout.write("\t-- Use triads %s ...\n"%(os.path.basename(options.triads_file)))
        triads_obj = triads.Triads(os.path.abspath(options.triads_file))
    if options.verbose: sys.stdout.write("\t-- Get PWMs of %s %s %s %s ...\n"%(tf_id,motif,family,template))
    # Get structural homologs and find out which is the family of the input pdb #
    potential_pdb_chain = template
    # Iterate over data soruces #
    for data in ["pdb", "pbm"]:
        # Use both general or family potentials #
        #for fam in ["family"]:   
        for fam in ["general", "family"]:
            # Skip working with family if the input pdb doesn't has a clear pdb chain #
            if (fam == "family") and (family == "Unknown"):
                continue
            # Iterate over potentials #
            for taylor in [True, False]:
                for comp in ["acc", "bins"]:
                    # Select the potentials file #
                    if data == "pdb":
                        potentials_file = os.path.join(pdb_dir, "potentials", "")
                    elif data == "pbm":
                        potentials_file = os.path.join(pbm_dir, "potentials", "")
                    if fam == "general":
                        potentials_file += fam
                    elif (fam == "family") and (potential_pdb_chain != None):
                        potentials_file += potential_pdb_chain
                    if pmf == True:
                        potentials_file += ".pmf"
                    if taylor == True:
                        potentials_file += ".taylor"
                    if comp == "acc":
                        potentials_file += ".acc"
                    elif comp == "bins":
                        potentials_file += ".bins"
                    potentials_file += ".txt"
                    # Initialize potentials, binding sites and scores #
                    for prot in pdb_obj.proteins:
                      pdb_chain = prot.chain
                      if options.verbose: sys.stdout.write("\t\t-- chain %s...\n"%(pdb_chain))
                      for dist in [15.0, 22.0, 30.0]:
                        potentials={}
                        potentials.setdefault(pdb_chain,spotentials.Potentials(file_name=potentials_file, select_potential=split_potential))
                        # For each PDB chain... -> there is no need of splitting by pdb chain, we work for separated pdb chains #
                        # Get dinucleotide raw scores and binding site region #
                        radii={}
                        radii.setdefault(pdb_chain,dist)
                        for i in selected_thresholds:
                            thresholds={}
                            thresholds.setdefault(pdb_chain,float(i)/100.0)
                            # Define output files #
                            label = tf_id + "_mL_" + motif  + "_mL_" + family + "_mL_" + potential_pdb_chain + "_mL_" + os.path.basename(potentials_file).replace(".txt", "") + "_mL_" + data + "_mL_" + str(dist) + "_mL_" + str(i) + "_mL_" + pdb_chain
                            msa_file = os.path.join(os.path.abspath(options.output_dir), "msa", label + ".msa")
                            meme_file = os.path.join(os.path.abspath(options.output_dir), "pwms", label + ".meme.s")
                            #Check if they exist, otherwise make them
                            #if  os.path.exists(meme_file) and os.path.exists(msa_file): 
                            if  os.path.exists(meme_file) : 
                                if options.verbose: sys.stdout.write("\t\t-- Use PWM of %s ...\n"%meme_file)
                                continue
                            elif os.path.exists(msa_file):
                                if options.verbose: sys.stdout.write("\t\t-- Use PWM of %s ...\n"%msa_file)
                                msa_obj= PWM.nMSA(msa_file,label,"msa")
                            elif os.path.exists(msa_file+".Z") or os.path.exists(msa_file+".gz"):
                                if options.verbose: sys.stdout.write("\t\t-- Use PWM of compressed %s ...\n"%msa_file)
                                msa_obj= PWM.nMSA(msa_file,label,"msa",gz=True)
                            else:
                                if options.verbose: sys.stdout.write("\t\t-- Calculate PWM of %s ...\n"%msa_file)
                                msa_obj= PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, None, None, split_potential,thresholds)
                            # Write output files #
                            if options.verbose: sys.stdout.write("\t\t-- write PWM of %s ...\n"%label)
                            #msa_obj.write(file_name=msa_file, option="msa")
                            msa_obj.write(file_name=meme_file, option="meme")
                            #os.system("gzip %s"%(msa_file))
                            

                            

