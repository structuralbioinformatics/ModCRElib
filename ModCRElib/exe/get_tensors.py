"""
Generate training/testing tensor datasets from aligned PWMs.

Builds randomized PWM combinations grouped by similarity intervals and exports
pickled tensor sets for downstream ML workflows.
"""

import sys,os
import pandas as pd
import json
import random
import configparser
import pickle
from collections import defaultdict
import itertools
import numpy as np
import optparse
import time

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
from ModCRElib.msa import aggregate_pwms as APWM


#-------------#
# Options     #
#-------------#

def parse_options():
    """
    Parse CLI options for train/test PWM tensor generation.

    The command reads one JSON (or directory of JSON files) describing PWM/model
    similarities, then combines aligned PWMs from ``aggregate_pwms`` outputs into
    randomized training/testing tensor pickles.

    Returns:
        optparse.Values: Parsed CLI namespace.
    """
    parser = optparse.OptionParser("python get_tensors.py -i input_json -a aggregate_folder[ -o output_name   ]  ")
    parser.add_option("-i","--json", action="store", default=None, type="string", dest="input_json", help="Input JSON file or directory of JSON files with information of PWMs", metavar="{filename}")
    parser.add_option("--jaspar", action="store", default=None, type="string", dest="jaspar", help="Address of JASPAR PWMs", metavar="{filename}")
    parser.add_option("--cisbp", action="store", default=None, type="string", dest="cisbp", help="Address of CisBP PWMs", metavar="{filename}")
    parser.add_option("--hocomoco", action="store", default=None, type="string", dest="hocomoco", help="Address of HOCOMOCO PWMs", metavar="{filename}")
    parser.add_option("--a","--aggregate", action="store", default=None, type="string", dest="input_folder", help="Directory with aligned  PWMs as obtained by aggregate_pwms script", metavar="{filename}")
    parser.add_option("-n","--n_train",action="store", default=500, type="int",dest="num_train",help="Number of combinations of training profiles per target (default is 500)",metavar="{integer}")
    parser.add_option("-m","--n_test",action="store", default=100, type="int",dest="num_test",help="Number of combinations of testing profiles per target (default is 100)",metavar="{integer}")
    parser.add_option("-r", "--repetition", default=False, action="store_true", dest="repetition", help="Apply repeating combinations of testing profiles to mantain the balance of tests (default = False)", metavar="{boolean}")
    parser.add_option("-o", "--output", action="store", default="train_and_test", type="string", dest="outdir", help="Output directory with training and testing tensors (default is train_and_test)", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
 
    (options, args) = parser.parse_args()

    if options.input_json is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Functions   #
#-------------#

def bin_files_by_interval(json_file, restriction=None, start=30, end=100, step=5):
    """
    Bin PWM identifiers into similarity intervals.

    Args:
        json_file (str): JSON file with fields ``pwm``, ``database``, ``similarity``.
        restriction (set | list | None): Optional whitelist of PWM ids.
        start (int): Lower bound for interval generation.
        end (int): Upper bound for interval generation.
        step (int): Interval width.

    Returns:
        dict: ``(low, high) -> [pwm_ids]`` plus ``(0,0)`` bucket for ``modcre`` models.
    """
    intervals = {(i, i + step): [] for i in range(start, end, step)}
    intervals[(0,0)]=[]
    with open(json_file,"r") as json_fp:
         json_obj=json.load(json_fp)
    pwms=[]
    for model in json_obj:
        pwm_file = model["pwm"].rstrip(".meme")
        if restriction is not None:
           if pwm_file not in  restriction: continue
        if model["database"]== "modcre":
            pwms.append((0.0,pwm_file))
        else:
            pwms.append((100*float(model["similarity"]),pwm_file))
    for f in pwms:
        number = f[0]
        if number>0:
          for (low, high) in intervals:
            if low < number <= high:
                intervals[(low, high)].append(f[1])
                break
        else:
          intervals[(0,0)].append(f[1])
    return intervals

def random_configuration(intervals, allow_empty=True, empty_prob=0.5, modcre_pick_n=10):
    """
    Draw one random profile configuration across interval bins.

    Args:
        intervals (dict): Output from ``bin_files_by_interval``.
        allow_empty (bool): Whether non-modcre intervals may be left empty.
        empty_prob (float): Probability of leaving a non-modcre interval empty.
        modcre_pick_n (int): Number of entries to sample from the ``(0,0)`` bucket.

    Returns:
        dict: ``interval -> selected_pwm | None | [selected_modcre_pwms]``.
    """
    config = {}
    for interval, pwms in intervals.items():
        if not pwms:  # no files available, must be empty
            config[interval] = None
        elif interval==(0,0):
            if len(pwms)<modcre_pick_n: select=len(pwms)
            else: select=modcre_pick_n
            config[interval] = random.sample(pwms, select)
            if select<modcre_pick_n:
               config[interval].extend((modcre_pick_n-select)*[None])
            random.shuffle(config[interval])
        else:
            if allow_empty and random.random() < empty_prob:
                # leave interval empty with probability `empty_prob`
                config[interval] = None
            else:
                config[interval] = random.choice(pwms)
    return config

def generate_random_possibilities(json_file, restriction=None, n_possibilities=100, allow_empty=True, empty_prob=0.5, modcre_pick_n=10,repetition=True):
    """
    Generate many random interval configurations.

    Args:
        json_file (str): JSON file used to build interval bins.
        restriction (set | list | None): Optional whitelist of PWM ids.
        n_possibilities (int): Number of configurations to generate.
        allow_empty (bool): Allow empty intervals.
        empty_prob (float): Empty-interval probability when allowed.
        modcre_pick_n (int): Number sampled from ``(0,0)`` bucket.
        repetition (bool): If False, skip duplicate configurations.

    Returns:
        list[dict]: Random configurations as produced by ``random_configuration``.
    """
    intervals = bin_files_by_interval(json_file, restriction)
    possibilities = []
    for _ in range(n_possibilities):
        configuration = random_configuration(intervals, allow_empty, empty_prob,modcre_pick_n)
        if not repetition:
           if configuration in possibilities: continue
        possibilities.append(configuration)
    return possibilities

def generate_tensors(combinations,cluster_pwms,modcre_pick_n):
    """
    Convert random PWM configurations into numeric tensor combinations.

    Args:
        combinations (list[dict]): Output of ``generate_random_possibilities``.
        cluster_pwms (dict): ``motif_id -> aligned_meme_path``.
        modcre_pick_n (int): Expected leading slots for ``modcre`` members.

    Returns:
        list[tuple]: ``(max_interval, tensor_list)`` for non-empty combinations.
    """
    generated_list=[]
    for comb in combinations:
       cl_tensor=[]
       not_empty=False
       modcre_list = comb[(0,0)]
       if modcre_list is None:
          cl_tensor=modcre_pick_n*[None]
       else:
          for mdcr in modcre_list:
              if mdcr is None:
                 cl_tensor.append(None)
              else:
                 if mdcr not in cluster_pwms:
                    print("Skip ",mdcr)
                    cl_tensor.append(None)
                    continue
                 msa_obj = PWM.nMSA(cluster_pwms[mdcr],None,"meme")
                 input_tensor = APWM.pwm2tensor(msa_obj)
                 not_empty=True
                 cl_tensor.append(input_tensor)
       max_interval = (0,0)
       for interval in comb:
           if interval==(0,0): continue
           if comb[interval] is None or comb[interval] not in cluster_pwms: 
              cl_tensor.append(None)
           else:
              if interval[1] > max_interval[1]: max_interval=interval
              msa_obj = PWM.nMSA(cluster_pwms[comb[interval]],None,"meme")
              input_tensor = APWM.pwm2tensor(msa_obj)
              not_empty=True
              cl_tensor.append(input_tensor)
       if not_empty: generated_list.append( (max_interval,cl_tensor) )
    return generated_list

def get_modcre_pwms(json_file):
    """
    Extract PWM identifiers belonging to the ``modcre`` database from JSON.

    Args:
        json_file (str): Input model/PWM JSON.

    Returns:
        list[str]: PWM ids (without ``.meme`` suffix).
    """
    pwms=[]
    with open(json_file,"r") as json_fp:
         json_obj=json.load(json_fp)
    for model in json_obj:
        db=model["database"]
        if db=="modcre":
           pwms.append(model["pwm"].rstrip(".meme"))
    return pwms



# Program 
#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

  # Arguments & Options #
  options = parse_options()
  input_json   = options.input_json
  input_folder = options.input_folder
  outdir       = options.outdir
  jaspar       = options.jaspar
  cisbp        = options.cisbp
  hocomoco     = options.hocomoco
  num_train    = int(options.num_train)
  num_test     = int(options.num_test )
  repetition   = not options.repetition
  verbose      = options.verbose

  #Defintions
  ratio = 0.99 #ratio of similaritythat defines a PWM as reference
  modcre_pick_n = 6 #Number of ModCRE PWMs to include in the train/test tensor

  #Prepare the output
  if not os.path.exists(outdir): os.makedirs(outdir)

  #Reading files
  list_of_jsons=[]
  if os.path.isdir(input_json):
    for json_file in os.listdir(input_json):
        if json_file.endswith("json"): 
           list_of_jsons.append(os.path.join(input_json,json_file))
  else:
    list_of_jsons.append(input_json)


  for json_file in list_of_jsons:
    tf_name= os.path.basename(json_file).rstrip(".json")
    modcre_pwms = get_modcre_pwms(json_file)
    if verbose: print("- Making train/test of %s (%d models from ModCRE) ..."%(tf_name,len(modcre_pwms)))
    if os.path.exists(os.path.join(input_folder,tf_name)):
       output = os.path.join(outdir,tf_name)
       train  = os.path.join(outdir,tf_name,"train")
       test   = os.path.join(outdir,tf_name,"test")
       if not os.path.exists(output): os.makedirs(output)
       if not os.path.exists(train): os.makedirs(train)
       if not os.path.exists(test): os.makedirs(test)
       #Make the training sets
       reference_pwms = APWM.get_reference(json_file,jaspar,cisbp,hocomoco,ratio,True)
       if verbose: print("\t TRAINING ...")
       for pwm in reference_pwms:
           motif = os.path.basename(pwm).rstrip(".meme")
           if verbose: print("\t\t- Reference target: %s ..."%(motif)) 
           if not os.path.exists(os.path.join(input_folder,tf_name,"reference_"+motif)): continue
           msa_obj = PWM.nMSA(pwm,None,"meme")
           target_tensor = APWM.pwm2tensor(msa_obj)
           out=open(os.path.join(train,"target_"+motif+".pickle"),"wb")
           pickle.dump(target_tensor,out)
           out.close()
           clusters = os.listdir(os.path.join(input_folder,tf_name,"reference_"+motif))
           num_cl_train = int(float(num_train)/len(clusters))
           list_of_tensors=[]
           for cluster in clusters:
               cluster_pwms={}
               for cl_pwm in os.listdir(os.path.join(input_folder,tf_name,"reference_"+motif,cluster)):
                   if not cl_pwm.endswith("meme"): continue
                   if not "member_of_cluster" in cl_pwm: continue
                   cl_motif = cl_pwm.split("_member_of_cluster_")[0]
                   cluster_pwms.setdefault(cl_motif,os.path.join(input_folder,tf_name,"reference_"+motif,cluster,cl_pwm))
               restriction=list(cluster_pwms.keys())
               combinations = generate_random_possibilities(json_file, restriction, n_possibilities=num_cl_train, allow_empty=False, empty_prob=0.5, modcre_pick_n=modcre_pick_n)  
               combined_tensors = generate_tensors(combinations,cluster_pwms,modcre_pick_n)
               for interval, combi in combined_tensors:
                   list_of_tensors.append(combi)
           if verbose: print("\t\t- Saving training data %s "%("data_"+motif+".pickle"))
           out=open(os.path.join(train,"data_"+motif+".pickle"),"wb")
           pickle.dump(list_of_tensors,out)
           out.close()
       #Make the testing sets
       test_cluster=set()
       if len(modcre_pwms)<=0 : continue
       for cluster in os.listdir(os.path.join(input_folder,tf_name)):
           if not cluster.startswith("cluster"):continue
           if len([cl_pwm.split("_member_of_cluster_")[0] for cl_pwm in os.listdir(os.path.join(input_folder,tf_name,cluster)) if cl_pwm.split("_member_of_cluster_")[0] in modcre_pwms]):
              test_cluster.add(cluster)
       if verbose: print("\t TESTING ...")
       for cluster in  test_cluster:
           #if verbose: print("\t\t- Testing data for cluster %s..."%(cluster))
           cluster_pwms={}
           for cl_pwm in os.listdir(os.path.join(input_folder,tf_name,cluster)):
               if not cl_pwm.endswith("meme"): continue
               if not "member_of_cluster" in cl_pwm: continue
               cl_motif = cl_pwm.split("_member_of_cluster_")[0]
               cluster_pwms.setdefault(cl_motif,os.path.join(input_folder,tf_name,cluster,cl_pwm))
           restriction=list(cluster_pwms.keys())
           combinations = generate_random_possibilities(json_file, restriction, n_possibilities=num_test, allow_empty=False, empty_prob=0.5, modcre_pick_n=modcre_pick_n, repetition=repetition)  
           combined_tensors = generate_tensors(combinations,cluster_pwms,modcre_pick_n)
           test_tensors={}
           for interval,combi in combined_tensors:
               test_tensors.setdefault(interval,[]).append(combi)
           for interval in test_tensors:
              if verbose: print("\t\t- Saving testing data for cluster %s (maximum ID= %s)"%(cluster,str(interval[0])))
              out=open(os.path.join(test,"data_ID"+str(interval[0])+"_"+tf_name+"_"+cluster+".pickle"),"wb")
              pickle.dump(test_tensors[interval],out)
              out.close()

  #Exit    
  print("Done")




