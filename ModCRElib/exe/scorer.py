"""
Compute protein-DNA interaction scores from structures or threading models.

This module prepares the structural context required for scoring, evaluates
several statistical potentials over protein-DNA triads, and exposes a
container class for exporting, normalizing, serializing, and plotting the
resulting energy profiles.
"""

import os, sys, re
from collections import Counter
import configparser
import itertools
import numpy as np
import optparse
import shutil
import pickle
import random as RND

#import graphics
import matplotlib as mpl 
# agg backend is used to create plot as a .png file
mpl.use('Agg')
import matplotlib.pyplot as plt 

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
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.threading import pdb2thread, threading_to_triads, threader
from ModCRElib.potential import spotentials
from ModCRElib.structure.protein import dssp, tmalign
from ModCRElib.msa import  pwm_pbm as PWM

#-------------#
# Functions   #
#-------------#


def scale_energies(energies, split_potential, min_score, max_score):
    """
    Min-max scale one energy term according to its potential type.

    Args:
        energies (dict): Raw energies keyed by split-potential name.
        split_potential (str): Potential to scale.
        min_score (float): Lower bound used for scaling.
        max_score (float): Upper bound used for scaling.

    Returns:
        float: Scaled value for the selected potential.
    """

    norm_e = 0.0

    if split_potential == "3d":
        norm_e = PWM.scale(energies["3d"], max_score, min_score)
    if split_potential == "3dc":
        norm_e = PWM.scale(energies["3dc"], max_score, min_score)
    if split_potential == "local":
        norm_e = PWM.scale(-energies["local"], max_score, min_score)
    if split_potential == "pair":
        norm_e = PWM.scale(-energies["pair"], max_score, min_score)
    if split_potential == "s3dc":
        norm_e = PWM.scale(-energies["s3dc"], max_score, min_score)
    if split_potential == "s3dc_dd":
        norm_e = PWM.scale(-energies["s3dc_dd"], max_score, min_score)
    if split_potential == "s3dc_di":
        norm_e = PWM.scale(energies["s3dc_di"], max_score, min_score)


    return norm_e



def parse_data_scoring(pdb_dir,pbm_dir,dummy_dir,fragment_restrict,binding_restrict,input_pdb_file,input_threading_file,threading,template,random,dna_background,radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, norm,  verbose, methylation=False):
    """
    Prepare structural data and compute the main score plus random controls.

    Depending on the input mode, this function reads either a PDB structure or
    one or more threading files, derives the triads required for scoring,
    loads statistical potentials, calculates the observed score, and
    optionally builds random-score controls from background DNA sequences.

    Args:
        pdb_dir (str): Directory with structural data and optional family file.
        pbm_dir (str): Directory containing PBM-derived PWM data.
        dummy_dir (str): Working directory for temporary intermediate files.
        input_pdb_file (str): PDB file to score when ``threading`` is false.
        input_threading_file (str): Threading input file when ``threading`` is
            true.
        threading (bool): Whether the input is a threading description.
        random (int): Number of random control scores to generate.
        split_potential (str): Requested statistical potential or ``"all"``.
        verbose (bool): Whether to emit progress messages.
        methylation (bool): Whether methylated nucleotides are preserved.

    Returns:
        tuple: ``(scr, random_scores)`` where ``scr`` is the main ``score``
        object and ``random_scores`` is a list of control ``score`` objects.

    Raises:
        Exception: If structural parsing, potential loading, or scoring fails.
    """
    # Initialize #
    families = {}
    dummy_dir=dummy_dir+"_"+str(os.getpid())
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    # For family potentials... #
    if pdb_dir is not None:
     if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
       sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
     else:
      for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family
    if threading:
     try:
      try:
        pdb_name = None
        thr_obj = threader.Threaded(threading_file=input_threading_file)
        thr_obj._check_parsing()
        pdb_name =  thr_obj.get_pdb_name()
        if verbose: sys.stdout.write("\t-- Read single threading file on PDB %s ...\n"%pdb_name)
        if pdb_name is not None: 
           thread_complex=False
        else:
           fi=open(input_threading_file,"r")
           n=len([x for x in fi])
           fi.close()
           if n>50:
             if verbose: sys.stdout.write("\t-- Check your threading file: TOO LARGE COMPLEX...\n")
           if verbose: sys.stdout.write("\t-- Read complex threading file...\n")
           thread_complex=True
      except:
        fi=open(input_threading_file,"r")
        n=len([x for x in fi])
        fi.close()
        if n>50:
           if verbose: sys.stdout.write("\t-- Check your threading file: TOO LARGE COMPLEX...\n")
        if verbose: sys.stdout.write("\t-- Read complex threading file...\n")
        thread_complex=True
      if thread_complex:
       fi=open(input_threading_file,"r")
       if verbose: sys.stdout.write("\t-- Open %s..."%input_threading_file)
       threading_files = []
       for line in fi:
           threading_files.append(line.strip())
       fi.close()
       triads_obj=triads.Triads()
       pdb_obj=PDB()
       for threading_file in threading_files:
           # Get triads, pdb etc from threading file #
           if verbose: sys.stdout.write("\t-- Get triads from %s ...\n"%threading_file)
           try:
               triads_obj_file, pdb_obj_file, x3dna_obj = threading_to_triads.threading_triads(threading_file=threading_file, pdb_dir= pdb_dir, template=template, verbose=verbose, dummy_dir=dummy_dir)
           except Exception as e:
               print("Failed to read triads (%s)"%e)
           for chain_obj in pdb_obj_file.chains:
               if chain_obj.chaintype == "P": pdb_obj.add_chain(chain_obj)
           for triad_obj in triads_obj_file.get_triads():
               #print("Add triad ",triad_obj.return_as_string())
               triads_obj.add_triad(triad_obj)
       #Get x3dna obj with full DNA contacts (unnecessary, X3DNA works indenpendtly)
       #if template is not None:
       #   x3dna_obj = x3dna.get_x3dna_obj(template, dummy_dir)
       #else:
       #   threading_file = threading_files[0]
       #   thr_obj = threader.Threaded(threading_file=threading_files[0])
       #   thr_obj._check_parsing()
       #   pdb_name =  thr_obj.get_pdb_name()
       #   if verbose: sys.stdout.write("\t-- Get DNA binding ...\n")
       #   dummy_triads_obj_file, pdb_obj_file, dummy_x3dna_obj = threading_to_triads.threading_triads(threading_file=threading_file, pdb_dir= pdb_dir, template=template, verbose=verbose, dummy_dir=dummy_dir)
       #   for chain_obj in pdb_obj_file.chains:
       #        if chain_obj.chaintype == "N": pdb_obj.add_chain(chain_obj)
       #   chain_ids= "".join([c for c in pdb_obj.chain_identifiers])
       #   template_pdb = os.path.join(dummy_dir,pdb_name+"_"+chain_ids+".pdb")
       #   if not os.path.exists(template_pdb): pdb_obj.write(template_pdb)
       #   x3dna_obj = x3dna.get_x3dna_obj(template_pdb, dummy_dir)
       # Get data for random scores
       if random >0:
          if verbose: sys.stdout.write("\t-- Generate random DNA sequences ...\n")
          random_triads=[]
          background_sequence=""
          for line in functions.parse_file(dna_background):
              if line.startswith(">"):continue
              background_sequence += line.strip()          
          # Generate n random DNA sequences
          thr_obj = threader.Threaded(threading_file=threading_files[0])
          thr_obj._check_parsing()
          dna_random = []
          for nran in range(random):
              dna={}
              kmers=thr_obj.get_kmers()
              for dna_seq,binding in kmers.items():
                start=RND.randint(0,len(background_sequence)-len(dna_seq))
                dna.setdefault(background_sequence[start:start+len(dna_seq)],binding)
              dna_random.append(dna)
          # Generate random threading files
          for nran in range(random):
              dna=dna_random[nran]
              random_triads_obj=triads.Triads()
              for threading_file in threading_files:
                  thr_obj = threader.Threaded(threading_file=threading_file) 
                  thr_obj._check_parsing()
                  random_threading_file=os.path.abspath(os.path.join(os.path.dirname(threading_file),"random_threading_file_"+str(nran)+"_from_"+os.path.basename(threading_file)))
                  if not os.path.exists(random_threading_file):
                     thr_obj.set_dna(dna)
                     thr_obj.set_dna_fixed(dna)
                     thr_obj.write(random_threading_file)
                  if verbose: sys.stdout.write("\t-- %dth threading to get random triads ...\n"%(nran))
                  random_triads_obj_file,dummy_pdb_obj, dummy_x3dna_obj = threading_to_triads.threading_triads(threading_file=random_threading_file, pdb_dir= pdb_dir, template=template, verbose=verbose, dummy_dir=dummy_dir)
                  for triad_obj in random_triads_obj_file.get_triads():
                      random_triads_obj.add_triad(triad_obj)
                  # Clean dummy thread files
                  os.remove(random_threading_file)
              random_triads.append(random_triads_obj)
      else:
       # Get triads, pdb etc from threading file #
       triads_obj, pdb_obj, x3dna_obj = threading_to_triads.threading_triads(threading_file=input_threading_file, pdb_dir= pdb_dir, template=template, verbose=verbose, dummy_dir=dummy_dir)
       # Get data for random scores
       if random >0:
          random_triads=[]
          background_sequence=""
          for line in functions.parse_file(dna_background):
              if line.startswith(">"):continue
              background_sequence += line.strip()
          # Generate a random threading file
          thr_obj = threader.Threaded(threading_file=input_threading_file) 
          thr_obj._check_parsing()
          if verbose: sys.stdout.write("\t-- Generate random DNA sequences ...\n")
          for nran in range(random):
            random_threading_file=os.path.abspath(os.path.join(os.path.dirname(input_threading_file),"random_threading_file_"+str(nran)+"_from_"+os.path.basename(input_threading_file)))
            if not os.path.exists(random_threading_file):
              dna={}
              kmers=thr_obj.get_kmers()
              for dna_seq,binding in kmers.items():
                start=RND.randint(0,len(background_sequence)-len(dna_seq))
                dna.setdefault(background_sequence[start:start+len(dna_seq)],binding)
              thr_obj.set_dna(dna)
              thr_obj.set_dna_fixed(dna)
              random_threading_file=os.path.abspath(os.path.join(os.path.dirname(input_threading_file),"random_threading_file_"+str(nran)+"_from_"+os.path.basename(input_threading_file)))
              thr_obj.write(random_threading_file)
            if verbose: sys.stdout.write("\t-- %dth threading to get random triads ...\n"%(nran))
            random_triads_obj,dummy_pdb_obj, dummy_x3dna_obj = threading_to_triads.threading_triads(threading_file=random_threading_file, pdb_dir= pdb_dir, template=template, verbose=verbose, dummy_dir=dummy_dir)
            random_triads.append(random_triads_obj)
            # Clean dummy thread files
            os.remove(random_threading_file)
     except Exception as e:
       raise Exception("Failed to get structure data, error %s"%e)
              
    else:
     try:
       # Get PDB object #
       if verbose: sys.stdout.write("\t-- Reading PDB file %s ...\n"%(os.path.basename(input_pdb_file)))
       pdb_obj = PDB(input_pdb_file)
       # Get DSSP object #
       if verbose: sys.stdout.write("\t\t-- Get DSSP ...\n")
       dssp_obj = dssp.get_dssp_obj(input_pdb_file)
       # Get X3DNA object #
       if verbose: sys.stdout.write("\t\t-- Get DNA object ...\n")
       x3dna_obj = x3dna.get_x3dna_obj(input_pdb_file, dummy_dir)
       # Get contacts object #
       if verbose: sys.stdout.write("\t\t-- Get contacts ...\n")
       contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
       # Skip if no contacts #
       if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")
       # Get triads object #
       triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)
       #for triad_obj in triads_obj.get_triads():
       #        print("Add triad ",triad_obj.return_as_string())
     except Exception as e:
       raise Exception("Failed to get structure data")
     try:
       # Get data for random scores
       if random >0:
          random_triads=[]
          output_dir = os.path.abspath(os.path.join(dummy_dir,"threads"))
          if not os.path.exists(output_dir): os.makedirs(output_dir)
          name = pdb_obj.id
          pdb_file_in_thread = os.path.join(output_dir,os.path.basename(input_pdb_file))
          if not os.path.exists(pdb_file_in_thread): shutil.copy(input_pdb_file,pdb_file_in_thread)
          background_sequence=""
          for line in functions.parse_file(dna_background):
              if line.startswith(">"):continue
              background_sequence += line.strip()
          if verbose: sys.stdout.write("\t-- Generate threading files ...\n")
          try:
             thread_objs = pdb2thread.pdb2thread(pdb_obj,name,x3dna_obj,output_dir, None,verbose)
          except Exception as e:
             print(("Error: %s"%e))
             raise Exception("Failed to get threading files, error %s"%e)
          for nran in range(random):
              dna_size=0
              for thr_obj in thread_objs:
                  kmers=thr_obj.get_kmers()
                  for dna_seq in list(kmers.keys()):
                      if len(dna_seq)>dna_size: dna_size=len(dna_seq)
              start=RND.randint(0,len(background_sequence)-dna_size)
              random_triads_obj = triads.Triads()
              for thr_obj in thread_objs:
                thr_id = thr_obj.get_pdb_name()+"_"+thr_obj.get_pdb_chain()
                random_threading_file=os.path.abspath(os.path.join(output_dir,"random_threading_file_"+str(nran)+"_from_"+thr_id+".txt"))
                if not os.path.exists(random_threading_file):
                  dna={}
                  kmers=thr_obj.get_kmers()
                  for dna_seq,binding in kmers.items():
                    dna.setdefault(background_sequence[start:start+len(dna_seq)],binding)
                  thr_obj.set_dna(dna)
                  thr_obj.set_dna_fixed(dna)
                  thr_obj.write(random_threading_file)
                if verbose: sys.stdout.write("\t\t-- %dth threading to get random triads ...\n"%(nran))
                if verbose: sys.stdout.write("\t\t-- DNA sequence %s\n"%(background_sequence[start:start+len(dna_seq)]))
                thread_triads_obj,dummy_pdb_obj, dummy_x3dna_obj = threading_to_triads.threading_triads(threading_file=random_threading_file, pdb_dir= pdb_dir, template=None, verbose=verbose, dummy_dir=dummy_dir)
                for tri_obj in thread_triads_obj.get_triads():
                    random_triads_obj.add_triad(tri_obj)
              random_triads.append(random_triads_obj)
          # Clean dummy thread files
          shutil.rmtree(output_dir)
     except Exception as e:
       raise Exception("Failed to get random threads")


    # Load statistical potential
    if verbose: sys.stdout.write("\t-- Load potentials...\n")
    try:
      potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, None, dummy_dir,verbose)
    except Exception as e:
       raise Exception("Failed to get potentials")
    # Calculate all energy scores #
    if verbose: sys.stdout.write("\t-- Calculate scores...\n")
    scr=score()
    # Define all msa_objs with the same potential "s3dc_dd"
    msa_objs={}
    msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict,"s3dc_dd", thresholds, methylation)
    for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
        msa_objs.setdefault(sp,msa_obj)
    try:
      scr.calculate_energies_and_binding(triads_obj, x3dna_obj, msa_objs, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir,methylation)
      #bs=scr.get_binding_site()
      #print("NUCLEOTIDES",bs.keys())
      #for k,u in bs.iteritems():
      #  for tro in u:
      #    print("Nuc",k,"triad",tro.return_as_string())
    except Exception as e:
       raise Exception("Failed to get energy scores")
    # Generate random scores
    try:
      random_scores=[]
      if random>0:
       if verbose: sys.stdout.write("\t-- Calculate random energies...\n")
       for nran in range(random):
           random_triads_obj=random_triads[nran]
           if verbose: sys.stdout.write("\t\t-- %dth threaded triads...\n"%(nran))
           random_scr = score()
           random_scr.calculate_energies_and_binding(random_triads_obj, x3dna_obj, msa_objs, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir,methylation)
           random_scores.append(random_scr)
      if not verbose:
       shutil.rmtree(dummy_dir)

    except Exception as e:
       raise Exception("Failed to get random scores")

    return  scr,random_scores


#-------------#
# Classes     #
#-------------#

class score(object):
    """
    Store interaction energies and derived score tracks for one complex.

    The object keeps total energies, normalized energies, per-nucleotide and
    per-amino-acid contributions, plus the binding-site mapping used to
    summarize or visualize a scored protein-DNA complex.

    Object features:
        - Global raw energies per split-potential (`_energies`).
        - Global normalized energies per split-potential (`_norm_energies`).
        - Per-nucleotide energy contribution maps
          (`_energies_per_nucleotide`).
        - Per-amino-acid energy contribution maps
          (`_energies_per_aminoacid`).
        - Binding-site index to supporting triad assignments (`_binding_site`).
        - Default split-potential selector used by getters (`_potential`);
          this can be set to a specific potential, while `"all"` preserves all
          tracks and defaults score retrieval to `"s3dc_dd"`.
    """

    def __init__(self,binding_site=None, energies=None, energies_per_nucleotide=None, energies_per_aminoacid=None, norm_energies=None, potential=None):
        """
        Initialize a score container with optional precomputed tracks.

        Args:
            binding_site (dict, optional): Binding-site map keyed by nucleotide
                or dinucleotide positions, with supporting triad assignments.
            energies (dict, optional): Total raw energies per split-potential.
            energies_per_nucleotide (dict, optional): Per-position energy
                contributions per split-potential.
            energies_per_aminoacid (dict, optional): Per-residue energy
                contributions per split-potential.
            norm_energies (dict, optional): Normalized total energies per
                split-potential.
            potential (str, optional): Default split-potential selector used by
                score getters. If not provided, defaults to `"all"`.

        Returns:
            None.
        """

        self._energies={}
        self._energies_per_nucleotide={}
        self._energies_per_aminoacid={}
        self._norm_energies={}
        self._binding_site={}
        #Use all energies by default, but retrieves s3dc_dd with get_score by default
        self._potential="all"

        if energies is not None:
           self._energies=energies
        if energies_per_nucleotide is not None:
           self._energies_per_nucleotide=energies_per_nucleotide
        if energies_per_aminoacid is not None:
           self._energies_per_aminoacid=energies_per_aminoacid
        if norm_energies is not None:
           self._norm_energies=norm_energies
        if binding_site is not None:
           self._binding_site=binding_site
        if potential is not None:
           self._potential=potential

    def set_potential(self, potential=None):
        """Set the default split-potential returned by score accessors."""
        if potential is None: 
           self._potential="all"
        else:
           self._potential=potential

    def set_binding_site(self,binding_site):
        """Store the binding-site triads keyed by nucleotide or dinucleotide."""
        self._binding_site=binding_site

    def set_energies(self,energies):
        """Store raw total energies for each split-potential."""
        self._energies=energies

    def set_norm_energies(self,energies):
        """Store normalized total energies for each split-potential."""
        self._norm_energies=energies

    def set_energies_per_nucleotide(self,energies_per_nucleotide):
        """Store per-nucleotide energy contributions."""
        self._energies_per_nucleotide=energies_per_nucleotide

    def set_energies_per_aminoacid(self,energies_per_aminoacid):
        """Store per-amino-acid energy contributions."""
        self._energies_per_aminoacid=energies_per_aminoacid

    def get_energies(self):
        return self._energies

    def get_energies_per_nucleotide(self):
        return self._energies_per_nucleotide

    def get_energies_per_aminoacid(self):
        return self._energies_per_aminoacid

    def get_score_per_nucleotide(self, normal=False, potential=None):
        """Return the selected per-nucleotide score track."""
        #Select a score_per_nucleotide, default is defined by class or "s3dc_dd" if "all"
        if potential is None: split_potential = self._potential
        else:                 split_potential = potential
        if split_potential == "all": split_potential = "s3dc_dd"
        skip=False
        if potential not in self._energies: skip=True
        if normal and potential in self._energies:
            if abs(float(self._energies.get(potential))) < 1.0e-10: skip=True
        score_per_nucleotide={}
        for nuc in sorted(self._energies_per_nucleotide.keys()):
            if normal: 
                       if skip: score_per_nucleotide.setdefault(nuc,0.0)
                       else:    score_per_nucleotide.setdefault(nuc,float(100*self._energies_per_nucleotide[nuc][potential]/self._energies[potential]))
            else:      score_per_nucleotide.setdefault(nuc,float(self._energies_per_nucleotide[nuc][potential]))
        return score_per_nucleotide

    def get_score_per_aminoacid(self, normal=False, potential=None):
        """Return the selected per-amino-acid score track."""
        #Select a score_per_nucleotide, default is defined by class or "s3dc_dd" if "all"
        if potential is None: split_potential = self._potential
        else:                 split_potential = potential
        if split_potential == "all": split_potential = "s3dc_dd"
        skip=False
        if potential not in self._energies: skip=True
        if normal and potential in self._energies:
            if abs(float(self._energies.get(potential))) < 1.0e-10: skip=True
        score_per_aminoacid={}
        for aminoacid in sorted(self._energies_per_aminoacid.keys()):
            if normal: 
                       if skip: score_per_aminoacid.setdefault(aminoacid,0.0)
                       else:    score_per_aminoacid.setdefault(amimnoacid,float(100*self._energies_per_aminoacid[aminoacid][potential]/self._energies[potential]))
            else:      score_per_aminoacid.setdefault(aminoacid,float(self._energies_per_amimnoacid[aminoacid][potential]))
        return score_per_aminoacid

    def get_score(self, normal=False, potential=None):
        """Return one total score for the selected split-potential."""
        #Select a single score, default is defined by class or "s3dc_dd" if "all"
        if potential is None: split_potential = self._potential
        else:                 split_potential = potential
        if split_potential == "all": split_potential = "s3dc_dd"
        if normal:
            if split_potential in self._norm_energies: return  self._norm_energies[split_potential]
            else: return  0.0
        else:
            if split_potential in self._energies: return  self._energies[split_potential]
            else: return  0.0

    def get_norm_energies(self):
        return self._norm_energies

    def get_binding_site(self):
        return self._binding_site

    def get_potential(self):
        return self._potential

    def calculate_energies_and_binding(self, triads_obj, x3dna_obj, msa_objs, potentials, thresholds, radii, radius=0, fragment_restrict=None, binding_restrict=None, dummy_dir="/tmp",methylation=False):
        """
        Compute total and per-position energies from a set of triads.

        Args:
            triads_obj: Triads describing protein-DNA contacts.
            x3dna_obj: X3DNA object with base-pair indexing information.
            msa_objs (dict): MSA objects used for potential normalization.
            potentials (dict): Statistical potentials keyed by chain.
            thresholds: Thresholds associated with the loaded potentials.
            radii: Distance metadata returned with the potentials.
            radius (float): Maximum contact distance to include.
            fragment_restrict (dict, optional): Protein residue restrictions.
            binding_restrict (list, optional): DNA nucleotide restrictions.
            dummy_dir (str): Temporary directory for helper calculations.
            methylation (bool): Whether methylated bases are preserved.

        Returns:
            None. Results are stored on the current ``score`` instance.
        """
        # Initialize #
        split_potential = self._potential
        E_3d = 0
        E_3dc = 0
        E_local = 0
        E_pair = 0
        E_s3dc = 0
        E_s3dc_dd = 0
        E_s3dc_di = 0
        binding_site = {}
        dna_idx={}
        energies={}
        energies_per_nucleotide={}
        energies_per_aminoacid={}
        if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
        # Get original DNA sequence indexes
        basepairs=x3dna_obj.get_basepairs()
        for basepair in basepairs.keys():
              (fwd_pdb_chain, fwd_residue_num), (rev_pdb_chain, rev_residue_num) = basepairs.get(basepair)
              dna_idx.setdefault((fwd_pdb_chain,fwd_residue_num),int(basepair))
        # For each PDB chain... #
        twonucleotides={}
        for pdb_chain in sorted(potentials):
              # For each contact... #
              for triad_obj in triads_obj.get_triads():
                  # Initialize #
                  a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(";")
                  chain, residue_num = residue_A.split("-")
                  aminoacid = (chain,residue_num)
                  #check restriction in TF
                  if fragment_restrict is not None:
                   if pdb_chain in fragment_restrict:
                     belongs=False
                     for interval in fragment_restrict.get(pdb_chain):
                         if int(residue_num) >= int(interval[0]) and int(residue_num) <= int(interval[1]): belongs=True
                         if belongs: break
                     if not belongs: continue
                  dinucleotide = PWM.get_triad_dinucleotide(triad_obj, x3dna_obj)
                  bp_1F,bp_1R,bp_2F,bp_2R=residue_B.split(",")
                  # Skip if not the right chain... #
                  if chain != pdb_chain: continue
                  # Add dinucleotide to binding site and PDB numbering to list of twonucleotides#
                  chain_dna,base_number_1=bp_1F.split("-")
                  chain_dna,base_number_2=bp_2F.split("-")
                  try:
                     binucleotide = (dna_idx[(chain_dna,int(base_number_1))],dna_idx[(chain_dna,int(base_number_2))])
                  except:
                    sys.stdout.write("\t\t-Skip dinucleotide chain %s %s-%s unindexed\n"%(chain_dna,base_number_1,base_number_2))
                    continue
                  #Check restriction on DNA site
                  if binding_restrict is not None:
                     belongs=False
                     for interval in binding_restrict:
                         if int(binucleotide[0]) >= int(interval[0]) and int(binucleotide[1]) <= int(interval[1]): belongs=True
                         if belongs: break
                     if not belongs: continue
                  # confirm base_numbers from 3xdna object
                  confirm_bases=False
                  if x3dna_obj.has_dinucleotide(dinucleotide):
                     bp1,bp2=x3dna_obj.get_dinucleotide(dinucleotide)
                     if x3dna_obj.has_basepair(bp1) and x3dna_obj.has_basepair(bp2):
                        ((fwd_pdb_chain1, fwd_residue_num1), (rev_pdb_chain1, rev_residue_num1)) = x3dna_obj.get_basepair(bp1)
                        ((fwd_pdb_chain2, fwd_residue_num2), (rev_pdb_chain2, rev_residue_num2)) = x3dna_obj.get_basepair(bp2)
                        if fwd_pdb_chain1==chain_dna and fwd_pdb_chain2==chain_dna and int(fwd_residue_num1)==int(base_number_1) and int(fwd_residue_num2)==int(base_number_2):
                           confirm_bases=True
                  if not confirm_bases:
                     sys.stdout.write("Error: PDB residue numbers of nucleotides are different from 3XDNA definition of basepairs. Check your files\n")
                     raise ValueError("Error in SCORE: different basepairs")
                  #Check restriction on distance
                  dab = np.floor(float(distance))
                  if dab > radius: continue
                  #selected for binding site
                  try:
                   twonucleotides.setdefault(dinucleotide,(int(base_number_1),int(base_number_2)))
                   binding_site.setdefault(dinucleotide, []).append(triad_obj)
                   energies_per_nucleotide.setdefault(dinucleotide, {"3d":0.0, "3dc":0.0, "local":0.0, "pair":0.0, "s3dc":0.0, "s3dc_dd":0.0, "s3dc_di":0.0})
                   energies_per_aminoacid.setdefault(aminoacid, {"3d":0.0, "3dc":0.0, "local":0.0, "pair":0.0, "s3dc":0.0, "s3dc_dd":0.0, "s3dc_di":0.0})
                   # Skip if environment could not be obtained for amino acid #
                   if "None" in a_oa: continue
                   # Skip if environment could not be obtained for dinucleotide #
                   if "None" in b_ob: continue
                   # The following arguments arguments are defined as in the papers #
                   a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                   oa = "%s-%s-%s" % (hydrophobicity, degree_of_exposure, secondary_structure)
                   b, nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group = b_ob.split("-")
                   ob = "%s-%s-%s-%s" % (nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group)
                   a_b = "%s;%s" % (a, b)
                   a_b_oa_ob = "%s;%s" % (a_oa, b_ob)
                   oa_ob = "%s;%s" % (oa, ob)
                   ### Compute potentials for our input structure ### 
                   try:
                    if (split_potential == "3d") or (split_potential == "all"):
                      if potentials[chain].get_score("3d", distance=dab) is not None:
                          E_3d += potentials[chain].get_score("3d", distance=dab)
                          energies_per_nucleotide[dinucleotide]["3d"] += potentials[chain].get_score("3d", distance=dab)
                          energies_per_aminoacid[aminoacid]["3d"] += potentials[chain].get_score("3d", distance=dab)
                   except Exception as e:
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    print("3d",aminoacid,dinucleotide,a_b_oa_ob,a_b,dab,potentials[chain].get_score("3d",distance=dab))
                    continue  
                   try:
                    if (split_potential == "3dc") or (split_potential == "all"):
                      if potentials[chain].get_score("3dc", key=oa_ob, distance=dab) is not None:
                          E_3dc += potentials[chain].get_score("3dc", key=oa_ob, distance=dab)
                          energies_per_nucleotide[dinucleotide]["3dc"] += potentials[chain].get_score("3dc", key=oa_ob, distance=dab)
                          energies_per_aminoacid[aminoacid]["3dc"] += potentials[chain].get_score("3dc", key=oa_ob, distance=dab)
                   except Exception as e:
                    print("3dc",aminoacid,dinucleotide,a_b_oa_ob,a_b,dab,potentials[chain].get_score("3dc", key=oa_ob, distance=dab))
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    continue  
                   try:
                    if (split_potential == "local") or (split_potential == "all"):
                      if potentials[chain].get_score("local", key=a_oa) is not None:
                          E_local += potentials[chain].get_score("local", key=a_oa)
                          energies_per_nucleotide[dinucleotide]["local"] += potentials[chain].get_score("local", key=a_oa, distance=dab)
                          energies_per_aminoacid[aminoacid]["local"] += potentials[chain].get_score("local", key=a_oa, distance=dab)
                   except Exception as e:
                    print("local",aminoacid,dinucleotide,a_b_oa_ob,a_b,dab,potentials[chain].get_score("local", key=a_oa, distance=dab))
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    continue  
                   try:
                    if (split_potential == "pair") or (split_potential == "all"):
                      if potentials[chain].get_score("pair", key=a_b, distance=dab) is not None:
                          E_pair += potentials[chain].get_score("pair", key=a_b, distance=dab)
                          energies_per_nucleotide[dinucleotide]["pair"] += potentials[chain].get_score("pair", key=a_b, distance=dab)
                          energies_per_aminoacid[aminoacid]["pair"] += potentials[chain].get_score("pair", key=a_b, distance=dab)
                   except Exception as e:
                    print("pair",aminoacid,dinucleotide,a_b,dab,potentials[chain].get_score("pair", key=a_b, distance=dab))
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    continue  
                   try:
                    if (split_potential == "s3dc") or (split_potential == "all"):
                      if potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab) is not None:
                          E_s3dc += potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab)
                          energies_per_nucleotide[dinucleotide]["s3dc"] += potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab)
                          energies_per_aminoacid[aminoacid]["s3dc"] += potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab)
                   except Exception as e:
                    print("s3dc",aminoacid,dinucleotide,a_b_oa_ob,a_b,dab,potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab))
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    continue  
                   try:
                    if (split_potential == "s3dc_dd") or (split_potential == "all"):
                      if potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab) is not None:
                          E_s3dc_dd += potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab)
                          energies_per_nucleotide[dinucleotide]["s3dc_dd"] += potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab)
                          energies_per_aminoacid[aminoacid]["s3dc_dd"] += potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab)
                   except Exception as e:
                    print("s3dc_dd",aminoacid,dinucleotide,a_b_oa_ob,a_b,dab,potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab))
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    continue  
                   try:
                    if (split_potential == "s3dc_di") or (split_potential == "all"):
                      if potentials[chain].get_score("s3dc_di", key=a_b_oa_ob) is not None:
                          E_s3dc_di += potentials[chain].get_score("s3dc_di", key=a_b_oa_ob)
                          energies_per_nucleotide[dinucleotide]["s3dc_di"] += potentials[chain].get_score("s3dc_di", key=a_b_oa_ob)
                          energies_per_aminoacid[aminoacid]["s3dc_di"] += potentials[chain].get_score("s3dc_di", key=a_b_oa_ob)
                   except Exception as e:
                    print("s3dc_di",aminoacid,dinucleotide,a_b_oa_ob,a_b,dab,potentials[chain].get_score("s3dc_di", key=a_b_oa_ob, distance=dab))
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    continue  
                  except Exception as e:
                    sys.stdout.write("Skip triad due to a exception %s\n"%(e))
                    continue  
        # Get profile of energies per nucleotide (numbered as in PDB)#
        energies_per_nucleotide_processed = {}
        nucleotides=[]
        for dinucleotide,bps in twonucleotides.items():
              base1,base2=bps
              for i in (base1,base2):
                 try:
                    nucleotides.append(i)
                    if i not in energies_per_nucleotide_processed:   
                      energies_per_nucleotide_processed[i]={"3d":0.0,"3dc":0.0,"local":0.0,"pair":0.0,"s3dc":0.0,"s3dc_dd":0.0,"s3dc_di":0.0}
                    energies_per_nucleotide_processed[i]["3d"] += energies_per_nucleotide[dinucleotide]["3d"]
                    energies_per_nucleotide_processed[i]["3dc"] += energies_per_nucleotide[dinucleotide]["3dc"]
                    energies_per_nucleotide_processed[i]["local"] += energies_per_nucleotide[dinucleotide]["local"]
                    energies_per_nucleotide_processed[i]["pair"] += energies_per_nucleotide[dinucleotide]["pair"]
                    energies_per_nucleotide_processed[i]["s3dc"] += energies_per_nucleotide[dinucleotide]["s3dc"]
                    energies_per_nucleotide_processed[i]["s3dc_dd"] += energies_per_nucleotide[dinucleotide]["s3dc_dd"]
                    energies_per_nucleotide_processed[i]["s3dc_di"] += energies_per_nucleotide[dinucleotide]["s3dc_di"]
                 except Exception as e:
                    sys.stdout.write("Skip nucelotide %d due to a exception %s"%(i,e))
                    continue 
        # Energies per Nucleotide
        repeats=Counter(nucleotides)
        for i,s in repeats.items():
              for sp in energies_per_nucleotide_processed[i].keys():
                  energies_per_nucleotide_processed[i][sp]=energies_per_nucleotide_processed[i][sp]/s
        # Energies
        energies      = {"3d":E_3d, "3dc":E_3dc, "local":E_local, "pair":E_pair, "s3dc":E_s3dc, "s3dc_dd":E_s3dc_dd, "s3dc_di":E_s3dc_di}
        # Normalized energies
        norm_energies = {"3d":0.0, "3dc":0.0, "local":0.0, "pair":0.0, "s3dc":0.0, "s3dc_dd":0.0, "s3dc_di":0.0}
        # Get Msa_objs if they are not defined
        if split_potential == "all":
            for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
              try:
                if sp in msa_objs:
                   msa_obj = msa_objs.get(sp)
                else:
                   msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, sp, thresholds, methylation)
              except Exception as e:
                 sys.stdout.write("Error %s\n"%(e))
                 sys.stdout.write("\t-Skip MSA potential %s "%(sp))
                 continue
              try:
                min_score, max_score = PWM.get_min_max_scores(msa_obj, binding_site, x3dna_obj,  potentials,  sp, dummy_dir,methylation)
                norm_e = scale_energies(energies, sp, min_score, max_score)
                norm_energies[sp] = norm_e
              except Exception as e:
                sys.stdout.write("Error %s\n"%(e))
                sys.stdout.write("\t-Skip potential %s scaling [%f , %f] \n"%(sp,min_score, max_score))
                continue
        else:
          try:
            if split_potential in msa_objs:
               msa_obj = msa_objs.get(split_potential)
            else:
               msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, split_potential, thresholds,methylation)
          except Exception as e:
            sys.stdout.write("Error %s\n"%(e))
            sys.stdout.write("\t-Skip MSA potential %s \n"%(sp))
          try:
            min_score, max_score = PWM.get_min_max_scores(msa_obj, binding_site, x3dna_obj, potentials,  split_potential, dummy_dir,methylation)
            norm_e = scale_energies(energies, split_potential, min_score, max_score)
            norm_energies[split_potential] = norm_e
          except Exception as e:
            sys.stdout.write("Error %s\n"%(e))
            sys.stdout.write("\t-Skip potential %s scaling [%f , %f] \n"%(sp,min_score, max_score))
        
        # Constructor
        self._energies                = energies
        self._binding_site            = binding_site
        self._energies_per_nucleotide = energies_per_nucleotide_processed
        self._energies_per_aminoacid  = energies_per_aminoacid
        self._norm_energies           = norm_energies


    def to_pickle(self, output_file,  overwrite=False):
        """Serialize the score object to a pickle-friendly dictionary file."""
        # Create output file #
        if os.path.exists(output_file) and overwrite:
            os.remove(output_file)
        energies        = self.get_energies()
        norm_energies   = self.get_norm_energies()
        dummy={}
        for n in list(self.get_binding_site().keys()):
            dummy.setdefault(n,[])
        binding_site    = dummy
        potential       = self.get_potential()
        energies_per_nucleotide = self.get_energies_per_nucleotide()
        energies_per_aminoacid  = self.get_energies_per_aminoacid()
        x={}
        x.setdefault("energies",energies)
        x.setdefault("norm_energies",norm_energies)
        x.setdefault("binding_site",binding_site)
        x.setdefault("potential",potential)
        x.setdefault("energies_per_nucleotide",energies_per_nucleotide)
        x.setdefault("energies_per_aminoacid",energies_per_aminoacid)
        out = open(output_file,"wb")
        pickle.dump(x,out)
        out.close()

    def from_pickle(self, input_file):
        """Populate the score object from a pickle file created by ``to_pickle``."""
        # Read input file #
        inp = open(input_file,"rb")
        x=pickle.load(inp)
        inp.close()
        energies        = x["energies"]
        norm_energies   = x["norm_energies"]
        binding_site    = x["binding_site"]
        potential       = x["potential"]
        energies_per_nucleotide = x["energies_per_nucleotide"]
        energies_per_aminoacid  = x["energies_per_aminoacid"]
        self.set_energies(energies)
        self.set_norm_energies(norm_energies)
        self.set_binding_site(binding_site)
        self.set_potential(potential)
        self.set_energies_per_nucleotide(energies_per_nucleotide)
        self.set_energies_per_aminoacid(energies_per_aminoacid)

    def write(self, output_file, protein_name=None, normal=False, overwrite=False):
        """Write a human-readable text report with total and per-position scores."""

        # Create output file #
        if os.path.exists(output_file) and overwrite:
            os.remove(output_file)
        out = open(output_file, "a")

        # Get variables to write
        split_potential         = self._potential
        binding_site            = self._binding_site
        energies                = self._energies
        energies_per_nucleotide = self._energies_per_nucleotide
        energies_per_aminoacid  = self._energies_per_aminoacid
        norm_energies           = self._norm_energies
        if protein_name is None:  protein_name=".".join([x for x in os.path.basename(ouput_file).split(".")[0:-1]])

        # Write energies
        out.write("\n##############################################################################################################################################\n")
        out.write("##############################################################################################################################################\n")
        out.write("\n#%30s\t"%"Protein_scores")
        out.write("\t%20s\t%20s\t"%("Binding_site_length","Starting_site"))
        if split_potential == "3d" or split_potential == "all":       out.write("%15s\t"%"E_3d")
        if split_potential == "3dc" or split_potential == "all":      out.write("%15s\t"%"E_3dc")
        if split_potential == "local" or split_potential == "all":    out.write("%15s\t"%"E_local")
        if split_potential == "pair" or split_potential == "all":     out.write("%15s\t"%"E_pair")
        if split_potential == "s3dc" or split_potential == "all":     out.write("%15s\t"%"E_s3dc")
        if split_potential == "s3dc_dd" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_dd")
        if split_potential == "s3dc_di" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_di")
        out.write("\n")
        out.write("%30s\t"%protein_name)
        out.write("\t%20d\t%20d\t"%(len(binding_site)+1,min([int(x) for x in binding_site.keys()])))
        if split_potential == "3d" or split_potential == "all": out.write("%15.3f\t" % (float(energies["3d"])))
        if split_potential == "3dc" or split_potential == "all": out.write("%15.3f\t" % (float(energies["3dc"])))
        if split_potential == "local" or split_potential == "all": out.write("%15.3f\t" % (float(energies["local"])))
        if split_potential == "pair" or split_potential == "all": out.write("%15.3f\t" % (float(energies["pair"])))
        if split_potential == "s3dc" or split_potential == "all": out.write("%15.3f\t" % (float(energies["s3dc"])))
        if split_potential == "s3dc_dd" or split_potential == "all": out.write("%15.3f\t" % (float(energies["s3dc_dd"])))
        if split_potential == "s3dc_di" or split_potential == "all": out.write("%15.3f\t" % (float(energies["s3dc_di"])))
        out.write("\n")

        # Write Normalized energies
        if norm_energies is not None and len(list(norm_energies.keys()))>0 and normal:
          out.write("#%29s\t"%"Normalized (min-max scale)")
          if split_potential == "3d" or split_potential == "all":       out.write("%15s\t"%"E_3d")
          if split_potential == "3dc" or split_potential == "all":      out.write("%15s\t"%"E_3dc")
          if split_potential == "local" or split_potential == "all":    out.write("%15s\t"%"E_local")
          if split_potential == "pair" or split_potential == "all":     out.write("%15s\t"%"E_pair")
          if split_potential == "s3dc" or split_potential == "all":     out.write("%15s\t"%"E_s3dc")
          if split_potential == "s3dc_dd" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_dd")
          if split_potential == "s3dc_di" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_di")
          out.write("\n")
          out.write("%30s\t"%protein_name)
          if split_potential == "3d" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["3d"]))
          if split_potential == "3dc" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["3dc"]))
          if split_potential == "local" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["local"]))
          if split_potential == "pair" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["pair"]))
          if split_potential == "s3dc" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["s3dc"]))
          if split_potential == "s3dc_dd" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["s3dc_dd"]))
          if split_potential == "s3dc_di" or split_potential == "all": out.write("%15.3f\n" % float(norm_energies["s3dc_di"]))

        # Write Energies per nucleotide
        out.write("\n\n")
        out.write("#Nucleotide scores\n")
        out.write("%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n"%("#Nucleotide","E_3d","E_3dc","E_local","E_pair","E_s3dc","E_s3dc_dd","E_s3dc_di"))
        for nuc in sorted(energies_per_nucleotide.keys()):
           out.write("%15s\t" % str(nuc))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["3d"]) )
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["3dc"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["local"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["pair"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc_dd"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc_di"]))
           out.write("\n")

        # Write Energies per aminoacid
        residues={}
        for chain,residue in (list(energies_per_aminoacid.keys())):
            residues.setdefault(chain,set()).add(int(residue))
        skip={}
        for chain in sorted(residues.keys()):
            for residue in sorted([int(res_num) for res_num in residues[chain]]):
                 aminoacid=(str(chain),str(residue))
                 check=False
                 for  sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                     if  energies_per_aminoacid[aminoacid][sp] !=0: check=True
                 skip.setdefault(aminoacid,not check)
        out.write("\n")
        for chain in sorted(residues.keys()):
            out.write("#Amino-acid scores for chain %s\n"%str(chain))
            out.write("%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n"%("#Aminoacid","E_3d","E_3dc","E_local","E_pair","E_s3dc","E_s3dc_dd","E_s3dc_di"))
            for residue in sorted(residues[chain]):
                aminoacid=(str(chain),str(residue))
                if skip[aminoacid]:continue
                out.write("%15s\t" % (str(residue)+"-"+str(chain)))
                out.write("%15.3f\t" % float(energies_per_aminoacid[aminoacid]["3d"]) )
                out.write("%15.3f\t" % float(energies_per_aminoacid[aminoacid]["3dc"]))
                out.write("%15.3f\t" % float(energies_per_aminoacid[aminoacid]["local"]))
                out.write("%15.3f\t" % float(energies_per_aminoacid[aminoacid]["pair"]))
                out.write("%15.3f\t" % float(energies_per_aminoacid[aminoacid]["s3dc"]))
                out.write("%15.3f\t" % float(energies_per_aminoacid[aminoacid]["s3dc_dd"]))
                out.write("%15.3f\t" % float(energies_per_aminoacid[aminoacid]["s3dc_di"]))
                out.write("\n")

        # Write Normalized energies per nucleotide
        # Define the sign of energy for normalization
        sign={}
        for sp in ["3d", "3dc", "s3dc_di"]:
                sign[sp]=1.0
        for sp in [ "local", "pair", "s3dc", "s3dc_dd"]:
                sign[sp]=-1.0

        out.write("\n")
        if  normal:
            out.write("#Normalized scores per nucleotide \n")
            out.write("%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n"%("#Nucleotide","E_3d","E_3dc","E_local","E_pair","E_s3dc","E_s3dc_dd","E_s3dc_di"))
            maxx={}
            minx={}
            for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                x=[sign[sp]*float(energies_per_nucleotide[nuc][sp]) for nuc in sorted(energies_per_nucleotide.keys())]
                maxx.setdefault(sp,max(x))
                minx.setdefault(sp,min(x))
            for nuc in sorted(energies_per_nucleotide.keys()):
                out.write("%15s\t" % str(nuc))
                for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                  if split_potential == sp or  split_potential == "all":
                    if (maxx[sp] - minx[sp]) >0: out.write("%15.3f\t"%((sign[sp]*float(energies_per_nucleotide[nuc]["3d"])-minx[sp])/(maxx[sp] - minx[sp])))
                    else:out.write("%15.3f\t" %(0.0))
                  else:
                    out.write("%15.3f\t" %(0.0))
                out.write("\n")

        # Write Normalized energies per aminoacid
        out.write("\n")
        if normal:
           for chain in sorted(residues.keys()):
             out.write("#Normalized scores per aminoacid (ratio scaling) for chain %s\n"%str(chain))
             out.write("%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n"%("#Aminoacid","E_3d","E_3dc","E_local","E_pair","E_s3dc","E_s3dc_dd","E_s3dc_di"))
             maxx={}
             minx={}
             skip={}
             for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                 x=[sign[sp]*float(energies_per_aminoacid[(str(chain),str(residue))][sp]) for residue in sorted([int(res_num) for res_num in residues[chain]])]
                 maxx.setdefault(sp,max(x))
                 minx.setdefault(sp,min(x))
             for residue in sorted([int(res_num) for res_num in residues[chain]]):
                 aminoacid=(str(chain),str(residue))
                 check=False
                 for  sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                     if  energies_per_aminoacid[aminoacid][sp] !=0: check=True
                 skip.setdefault(aminoacid,not check)
             for residue in sorted([int(res_num) for res_num in residues[chain]]):
                aminoacid=(str(chain),str(residue))
                if skip[aminoacid]:continue
                out.write("%15s\t" % (str(residue)+"-"+str(chain)))
                for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                  if split_potential == sp or split_potential == "all":
                    if (maxx[sp] - minx[sp]) >0:out.write("%15.3f\t" % ((sign[sp]*float(energies_per_aminoacid[aminoacid][sp])-minx[sp])/(maxx[sp] - minx[sp])))
                    else:out.write("%15.3f\t" %(0.0))
                  else:
                    out.write("%15.3f\t" %(0.0))
                out.write("\n")


        out.close()



    def get_zscore(self,random):
        """
        Convert this score into Z-scores relative to random control scores.

        Args:
            random (list): List of ``score`` objects generated from controls.

        Returns:
            score: New score object containing Z-scored totals and profiles.
        """
    
        binding_site = self.get_binding_site()
        other = self.__class__(binding_site=binding_site)

        # list of split_potentials
        split_potentials  = ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]

        #Get energies
        energies={}
        norm_energies={}
        for normal in [True,False]:
            for sp in split_potentials:
                x=[r.get_score(normal,sp) for r in random]
                x.append(self.get_score(normal,sp))
                x=np.array(x) 
                z=(self.get_score(normal,sp)-x.mean())/x.std()
                if np.isnan(z): z=0.0
                if np.isinf(z): z=0.0
                if normal:
                    norm_energies.setdefault(sp,z)
                else:
                    energies.setdefault(sp,z)
        other.set_energies(energies)
        other.set_norm_energies(norm_energies)

        #Get zscores per nucleotide
        zscore_per_nucleotide={}
        energies_per_nucleotide = self.get_energies_per_nucleotide()
        for nuc in sorted(energies_per_nucleotide.keys()):
            zscore_per_nucleotide[nuc]={}
        for sp in split_potentials:
              ene = []
              for nuc in sorted(energies_per_nucleotide.keys()):
                    if sp in energies_per_nucleotide[nuc]:
                       ene.append(float(energies_per_nucleotide[nuc][sp]))
                    else:
                       ene.append(0.0)
              ene=np.array(ene)
              y=ene
              zscr = np.zeros(len(list(energies_per_nucleotide.keys())))
              for r in random:
                  random_per_nucleotide=r.get_energies_per_nucleotide()
                  x=[]
                  for nuc in sorted(energies_per_nucleotide.keys()):
                      if sp in random_per_nucleotide[nuc]:
                         x.append(float(random_per_nucleotide[nuc][sp]))
                      else:
                         x.append(0.0)
                  y = np.vstack((y,np.array(x)))
              zscr = (ene - y.mean(0))/y.std(0)
              zscr[np.isnan(zscr)]=0.
              zscr[np.isinf(zscr)]=0.
              n=0
              for nuc in sorted(energies_per_nucleotide.keys()):
                  zscore_per_nucleotide[nuc][sp]=zscr[n]
                  n=n+1
        other.set_energies_per_nucleotide(zscore_per_nucleotide)


        #Get zscores per amino-acid
        zscore_per_aminoacid={}
        energies_per_aminoacid = self.get_energies_per_aminoacid()
        residues={}
        for chain,residue in list(energies_per_aminoacid.keys()):
            residues.setdefault(chain,set()).add(int(residue))
        #Get length of amino-acids
        maxlen={}
        for chain in list(residues.keys()):
            maxlen.setdefault(chain,max(residues[chain]))
        for chain in list(maxlen.keys()):
            size   = maxlen[chain] + 1
            for residue in range(1,size):
                aminoacid=(str(chain),str(residue))
                zscore_per_aminoacid[aminoacid]={"3d":0.0, "3dc":0.0, "local":0.0, "pair":0.0, "s3dc":0.0, "s3dc_dd":0.0, "s3dc_di":0.0}
        for chain in list(maxlen.keys()):
            size   = maxlen[chain] + 1
            for sp in split_potentials:
              ene_aa = []
              for residue in range(1,size):
                aminoacid=(str(chain),str(residue))
                if aminoacid in energies_per_aminoacid:
                    if sp in energies_per_aminoacid[aminoacid]:
                       ene_aa.append(float(energies_per_aminoacid[aminoacid][sp]))
                    else:
                       ene_aa.append(0.0)
                else:
                   ene_aa.append(0.0)
              ene_aa=np.array(ene_aa)
              y=ene_aa
              zscr_aa = np.zeros(size)
              for r in random:
                  random_per_aminoacid=r.get_energies_per_aminoacid()
                  x=[]
                  for residue in range(1,size):
                    aminoacid=(str(chain),str(residue))
                    if aminoacid in random_per_aminoacid:
                      if sp in random_per_aminoacid[aminoacid]:
                         x.append(float(random_per_aminoacid[aminoacid][sp]))
                      else:
                         x.append(0.0)
                    else:
                       x.append(0.0)
                  y = np.vstack((y,np.array(x)))
              zscr_aa = (ene_aa - y.mean(0))/y.std(0)
              zscr_aa[np.isnan(zscr_aa)]=0.
              zscr_aa[np.isinf(zscr_aa)]=0.
              for residue in range(1,size):
                  aminoacid=(str(chain),str(residue))
                  zscore_per_aminoacid[aminoacid][sp]=zscr_aa[residue-1]
        other.set_energies_per_aminoacid(zscore_per_aminoacid)

        return other
            
    def get_protein_plot(self,output_name):
        """Save per-amino-acid score plots for each split-potential."""
        split_potentials  = ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]
        energies_per_aminoacid = self.get_energies_per_aminoacid()
        #print("PLOT ",energies_per_aminoacid)
        residues={}
        for chain,residue in list(energies_per_aminoacid.keys()):
            residues.setdefault(chain,set()).add(int(residue))
        #print("RESIDUES DICT",residues)
        #Define the sign of energy for normalization
        sign={}
        for sp in ["3d", "3dc", "s3dc_di"]:
                sign[sp]=1.0
        for sp in [ "local", "pair", "s3dc", "s3dc_dd"]:
                sign[sp]=-1.0

        #Get length of amino-acids
        maxlen={}
        for chain in list(residues.keys()):
            maxlen.setdefault(chain,max(residues[chain]))
        for norm in [True,False]:
         for chain in list(residues.keys()):
          #print("RESIDUES ",sorted(residues[chain]))
          if norm:
             maxx={}
             minx={}
             skip_zero={}
             for sp in split_potentials:
                 x=[sign[sp]*float(energies_per_aminoacid[(str(chain),str(residue))][sp]) for residue in sorted([int(res_num) for res_num in residues[chain]])]
                 maxx.setdefault(sp,max(x))
                 minx.setdefault(sp,min(x))
             for residue in sorted([int(res_num) for res_num in residues[chain]]):
                 aminoacid=(str(chain),str(residue))
                 check=False
                 for  sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                     if not check and energies_per_aminoacid[aminoacid][sp] !=0: check=True
                 skip_zero.setdefault(aminoacid,not check)
          for sp in split_potentials:
            label=chain+"_"+sp+"_x_aminoacid"
            if norm: label=label+"_normalized"
            out = output_name+"_"+label+".png"
            title = "Graph of energy %s"%label
            x=[]
            y=[]
            skip=False
            if norm:
               if abs(float(maxx[sp]-minx[sp])) < 1.0e-10: skip=True
            for residue in range(1,maxlen[chain]):
                aminoacid=(str(chain),str(residue))
                x.append(residue)
                if aminoacid in energies_per_aminoacid:
                   if sp in energies_per_aminoacid[aminoacid]:
                      if norm:
                         if skip:y.append(0.0)
                         elif skip_zero[aminoacid]:y.append(0.0)
                         else:   y.append((sign[sp]*float(energies_per_aminoacid[aminoacid][sp])-minx[sp])/(maxx[sp]-minx[sp]))
                      else:
                         y.append(energies_per_aminoacid[aminoacid][sp])
                   else:
                      y.append(0.0)
                else:
                   y.append(0.0)
            #print("PLOT ",sp,chain,norm,x,y)
            plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7, mfc='cyan')
            plt.xlabel('residue ')
            plt.ylabel('score ')
            plt.savefig(out,format="png")
            plt.close()

    def get_DNA_plot(self,output_name):
        """Save per-nucleotide score plots for each split-potential."""
        split_potentials  = ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]
        #Define the sign of energy for normalization
        sign={}
        for sp in ["3d", "3dc", "s3dc_di"]:
                sign[sp]=1.0
        for sp in [ "local", "pair", "s3dc", "s3dc_dd"]:
                sign[sp]=-1.0
        energies_per_nucleotide = self.get_energies_per_nucleotide()
        for norm in [True,False]:
          if norm:
             maxx={}
             minx={}
             for sp in split_potentials:
                x=[sign[sp]*float(energies_per_nucleotide[nuc][sp]) for nuc in sorted(energies_per_nucleotide.keys())]
                maxx.setdefault(sp,max(x))
                minx.setdefault(sp,min(x))
          for sp in split_potentials:
            label=sp+"_x_nucleotide"
            if norm: label=label+"_normalized"
            out = output_name+"_"+label+".png"
            title = "Graph of energy %s"%label
            y=[]
            x=[]
            skip=False
            if norm:
               if abs(float(maxx[sp]-minx[sp])) < 1.0e-10: skip=True
            for nucleotide in list(energies_per_nucleotide.keys()):
                x.append(nucleotide)
                if sp in energies_per_nucleotide[nucleotide]:
                   if norm:
                      if skip:y.append(0.0)
                      else:   y.append((sign[sp]*float(energies_per_nucleotide[nucleotide][sp])-minx[sp])/(maxx[sp]-minx[sp]))
                   else:
                      y.append(energies_per_nucleotide[nucleotide][sp])
            plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7, mfc='cyan')
            plt.xlabel('nucleotide ')
            plt.ylabel('score ')
            plt.savefig(out,format="png")
            plt.close()
 






#-------------#
# Options     #
#-------------#

def parse_options():
    """
    Parse command-line options for standalone score calculation.

    How to run:
        python scorer.py [--dummy DUMMY_DIR] -i INPUT_FILE
            [--pdb PDB_DIR --pbm PBM_DIR -o OUTPUT_ROOT]
            [--threading --template TEMPLATE --random N --plot -v]
            [statistical potential options]

    Example:
        python scorer.py -i input.pdb --pdb pdb_data --pbm pbm_data
            -o SCORE --plot -v

    The parser configures:
        - Input structure/threading source and optional template override.
        - PBM/PDB resources, output naming, random/background controls.
        - Fragment/binding restrictions and runtime logging options.
        - Statistical-potential controls used by energy and score computation.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options covering the input structure,
        scoring setup, random controls, and output preferences.

    Raises:
        SystemExit: Triggered by ``OptionParser`` when arguments are missing or
        invalid.
    """

    parser = optparse.OptionParser("Usage: scorer.py [--dummy=DUMMY_DIR] -i INPUT_FILE [-l LABEL -o OUTPUT_DIR --pbm=PBM_dir] --pdb=PDB_DIR [-m -v --threading] [-a -f -p -s SPLIT_POTENTIAL -t THRESHOLD -k -b --taylor --file POTENTIAL --radius RADIUS]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input PDB|THREADED file. Mandatory", metavar="INPUT_FILE")
    parser.add_option("-l", action="store", default=None, type="string", dest="label", help="Label to organize the output files as 'label.energies.txt' ", metavar="LABEL")
    parser.add_option("-o", "--output", default="SCORE", action="store", type="string", dest="output_name", help="Output ROOTNAME it uses the folder selected by the rootname if any (default = SCORE)", metavar="OUTPUT_ROOTNAME")
    parser.add_option("--pbm", action="store", type="string", default=None, dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py). This is ,mandatory unless using option --file", metavar="PBM_DIR")
    parser.add_option("--pdb", action="store", type="string", default=None, dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py). This is ,mandatory unless using option --template", metavar="PDB_DIR")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. (Default is None it applies to all amino-acids)")
    parser.add_option("-n", "--norm", default=False, action="store_true", dest="norm", help="Normalize scores using a min-max scaling")
    parser.add_option("--binding", default=None, action="store", type="string", dest="binding_restrict", help="Binding site of DNA to apply the potential (uses basepair order from 3XDNA). Format is 'a-b;c-d': two regions between nucleotides a-b and c-d of the forward chain (first in PDB). (Default is None it applies to all nucleotides)")
    parser.add_option("--threading", default=False, action="store_true", dest="threading", help="Input file is a threading file of a PDB structure that either exist in the PDB folder of ModCRE or it is given with option '--template' or it is in the same folder as the threading file(default = False)")
    parser.add_option("--random", default=0,type=int,action="store", dest="random",help="Number of random DNA sequences to generate a distribution of random scores ", metavar="INTEGER")
    parser.add_option("--background",default=None, action="store", type="string", dest="dna_background", help="File of a DNA sequence in FASTA format. This is used as background to generate a distribution of random scores", metavar="FASTA")
    parser.add_option("--template", action="store", default=None, type="string", dest="template", help="PDB structure of the template on threading. Option threading must be active", metavar="TEMPLATE")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("--plot", default=False, action="store_true", dest="plot", help="Plot the scores per nucleotide and amino-acid (default = False)")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of PDBThread files that have failed and have been completed")

    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="all", action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = all)", metavar="{string}")
    group.add_option("-t", action="store", type="float", dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    group.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' from configuration and it uses default radius to load potentials automatically, otherwise it uses radius for both", metavar="{string}")
    group.add_option("--methylation", default=False, action="store_true", dest="methylation", help="Use methylated cytosine specificities of binding/non-binding (default = False)")


    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()

    if options.input_file is None or (options.pdb_dir is None and options.template is None) or (options.pbm_dir is None and options.potential_file is None):
        parser.error("missing arguments: type option \"-h\" for help")
    if options.random > 0: 
       if options.dna_background is None:
        parser.error("missing DNA background file: type option \"-h\" for help")
       if not os.path.exists(options.dna_background) :
        parser.error("empty DNA background file : type option \"-h\" for help")
    if not re.search("^3d$|^3dc$|^local$|^pair$|^s3dc$|^s3dc_dd$|^s3dc_di$|^all$", options.split_potential):
        parser.error("incorrect value for -s argument: type option \"-h\" for help")
        raise ValueError("Error in SCORE: Wrong input")
    return options


#-------------#
# Main        #
#-------------#

def main():
    """
    Run the standalone scoring workflow for one structure or threading input.

    Workflow:
        1. Parse command-line options and resolve input/output/resource paths.
        2. Load or compute scores (and optional random/background distributions).
        3. Write score summaries, pickles, and optional plots to disk.
        4. Persist execution status in the selected log destination.

    Args:
        None.

    Returns:
        None. Score outputs are written to the selected output paths.
    """

    # Arguments & Options #
    options = parse_options()
    dummy_dir = options.dummy_dir
    verbose = options.verbose
    methylation = options.methylation
    info = options.info
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    pdb_dir = options.pdb_dir
    if pdb_dir is not None:
       if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)
    pbm_dir = options.pbm_dir
    if pbm_dir is not None:
       if not pbm_dir.startswith("/"): pbm_dir = os.path.abspath(options.pbm_dir)
    if options.threading:
       input_threading_file=options.input_file
       input_pdb_file=None
       if not input_threading_file.startswith("/"): input_threading_file=os.path.abspath(options.input_file)
    else:
       input_pdb_file=options.input_file
       input_threading_file=None
       if not input_pdb_file.startswith("/"): input_pdb_file=os.path.abspath(options.input_file)
    output_name = options.output_name
    if not output_name.startswith("/"): output_name = os.path.abspath(options.output_name)
    output_dir=os.path.dirname(output_name)
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    output_file = output_name + ".txt"
    if options.label is not None:
        output_file = output_name +"."+ options.label + ".txt"

    if info is None:
        info = os.path.join(output_dir,"scoring.log")




    fragment_restrict=None
    if options.fragment is not None:
       fragment_restrict={}
       segment_fragment=options.fragment.split(";")
       for interval in segment_fragment:
           ac,bc = interval.split("-")
           if len(ac.split("_"))>0: a,chain1=ac.split("_")
           if len(bc.split("_"))>0: b,chain2=bc.split("_")
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

    try:
      pickle_file=output_file.rstrip(".txt")+".pickle"
      if options.random>0: 
        zscore_pickle_file=output_file.rstrip(".txt")+".zscore.pickle"
        if os.path.exists(pickle_file) and os.path.exists(zscore_pickle_file):
           if verbose:sys.stdout.write("\t-- Open %s...\n"%pickle_file)
           scr=score()
           scr.from_pickle(pickle_file)
           if verbose:sys.stdout.write("\t-- Open %s...\n"%zscore_pickle_file)
           zscr=score()
           zscr.from_pickle(zscore_pickle_file)
        else:
           scr, random_scr = parse_data_scoring(pdb_dir,pbm_dir,dummy_dir,fragment_restrict,binding_restrict,input_pdb_file,input_threading_file,options.threading,options.template,options.random,options.dna_background,options.radius, options.potential_file, options.split_potential, options.auto_mode, options.family_potentials, options.pbm_potentials, options.score_threshold, options.taylor_approach, options.pmf , options.bins, options.known,options.norm,verbose,methylation)
           sys.stdout.flush()
           #Get Z-scores
           if verbose:sys.stdout.write("\t-- Calculate Z-scores... \n")
           zscr = scr.get_zscore(random_scr)
      else:
        if os.path.exists(pickle_file):
           if verbose:sys.stdout.write("\t-- Open %s...\n"%pickle_file)
           scr=score()
           scr.from_pickle(pickle_file)
        else:
           scr, random_scr = parse_data_scoring(pdb_dir,pbm_dir,dummy_dir,fragment_restrict,binding_restrict,input_pdb_file,input_threading_file,options.threading,options.template,options.random,options.dna_background,options.radius, options.potential_file, options.split_potential, options.auto_mode, options.family_potentials, options.pbm_potentials, options.score_threshold, options.taylor_approach, options.pmf , options.bins, options.known,options.norm,verbose,methylation)
           sys.stdout.flush()
    except Exception as e:
        print(("Failed scoring by %s"%e))
        log_file=open(info,"a")
        log_file.write("%s\tFAIL\n"%(os.path.basename(options.input_file)))
        log_file.close()
        exit()

    
    
    # Write the output #
    # Original
    try:
      protein_name=".".join([x for x in os.path.basename(options.input_file).split(".")[0:-1]])
      scr.write(output_file, protein_name=protein_name, normal=options.norm, overwrite=True)
      pickle_file=output_file.rstrip(".txt")+".pickle"
      if not os.path.exists(pickle_file): scr.to_pickle(pickle_file)

      if options.plot:
        scr.get_protein_plot(output_name)
        scr.get_DNA_plot(output_name)

    except Exception as e:
        print(("Failed scoring by %s"%e))
        log_file=open(info,"a")
        log_file.write("%s\tFAIL\n"%(os.path.basename(options.input_file)))
        log_file.close()
        exit()

    try:

      # Zscore/Random energies
      if options.random>0:
        zscore_output_name=output_file.rstrip(".txt")+".zscore"
        zscore_output_file=output_file.rstrip(".txt")+".zscore.txt"
        zscore_pickle_file=output_file.rstrip(".txt")+".zscore.pickle"
        zscr.write(zscore_output_file, protein_name=protein_name, normal=options.norm, overwrite=True)
        if not os.path.exists(zscore_pickle_file): zscr.to_pickle(zscore_pickle_file)
        if options.plot:
           zscr.get_protein_plot(zscore_output_name)
           zscr.get_DNA_plot(zscore_output_name)

    except Exception as e:
        print(("Failed scoring by %s"%e))
        log_file=open(info,"a")
        log_file.write("%s\tFAIL\n"%(os.path.basename(options.input_file)))
        log_file.close()
        exit()

       
     # random_output_file=output_file.strip("txt")+".random"
     # if os.path.exists(random_output_file): os.remove(random_output_file)
     # for random_score in random_scr:
     #     random_score.write(random_output_file, protein_name=protein_name, normal=options.norm, overwrite=True)
    log_file=open(info,"a")
    log_file.write("%s\tDONE\n"%(os.path.basename(options.input_file)))
    log_file.close()

    print("Done")


if __name__ == "__main__":
    main()

