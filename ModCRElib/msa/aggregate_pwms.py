import os, sys, re
from collections import Counter
import configparser
import itertools
import numpy as np
import optparse
import subprocess
import time
import random
import shutil
import hashlib
import json

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
    merge different types of PWMs, cluster and align them to produce a new set of combined PWMS
    the input are json files with information of the PWMs (i.e. name, features and location)
    the output are json files that combine the aligned PWMs, each json file correspond to a cluster
    of similar and aligned PWMs with the specific length selected and the location is an output folder
    with the modified PWMs
    """
    parser = optparse.OptionParser("python aggregate_pwms.py -i input_json [ -l length -o output_name --dummy dummy_dir  ]  ")
    parser.add_option("-i", action="store", default=None, type="string", dest="input_json", help="Input JSON file or directory of JSON files with information of PWMs", metavar="{filename}")
    parser.add_option("-l","--binding_site", action="store", default=50, type="int", dest="length", help="Length of the binding site of the output PWMs", metavar="{integer}")
    parser.add_option("--pvalue", action="store", default=0.005,type="float", dest="pvalue", help="P-value threshold of TOMTOM similarity between two PWMs (default 0.005)", metavar="{float}")
    #parser.add_option("--threshold", action="store", default=0.25,type="float", dest="threshold", help="Ratio of the maximum distance between dissimilar PWMs to use as distance threshold of the agglomerative clustering (default is 0.25)", metavar="{float}")
    parser.add_option("--threshold", action="store", default=0.01,type="float", dest="threshold", help="Distance threshold of the agglomerative clustering (default is 0.01)", metavar="{float}")
    parser.add_option("--gap_penalty", action="store", default=0.05,type="float", dest="gap_penalty", help="Gap penalty for the alignment of profiles (default is 0.05, but it is affected by the method)", metavar="{float}")
    parser.add_option("--reference", default=None, action="store", dest="reference", help="Threshold ratio of similarity to select PWMS of experimental databases as reference (default is None and it uses only centers of each cluster as reference)", metavar="{float}")
    parser.add_option("--jaspar", action="store", default=None, type="string", dest="jaspar", help="Address of JASPAR PWMs", metavar="{filename}")
    parser.add_option("--cisbp", action="store", default=None, type="string", dest="cisbp", help="Address of CisBP PWMs", metavar="{filename}")
    parser.add_option("--hocomoco", action="store", default=None, type="string", dest="hocomoco", help="Address of HOCOMOCO PWMs", metavar="{filename}")
    parser.add_option("--modcre", action="store", default=None, type="string", dest="modcre", help="Address of predicted PWMs with ModCRE", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("--kullback_leibler", default=False, action="store_true", dest="kullback_leibler", help="Flag to select  symmetric Kullback Leibler divergence (a.k.a. Jeffreys divergence) to compare profiles (default is scalar product)", metavar="{boolean}")
    parser.add_option("--jensen_shanon", default=False, action="store_true", dest="jensen_shanon", help="Flag to select Jensen Shanon divergence to compare profiles (default is scalar product)", metavar="{boolean}")
    parser.add_option("--pearson", default=False, action="store_true", dest="pearson", help="Flag to select Pearson correlation to compare profiles (default is scalar product)", metavar="{boolean}")
    parser.add_option("--trim", default=False, action="store_true", dest="trim", help="Flag to trim the alignment by the PWM of reference (default uses selected binding-site length)", metavar="{boolean}")
    parser.add_option("--keep_msa", default=False, action="store_true", dest="keep", help="Flag to keep the MSA files (default removes them after being used)", metavar="{boolean}")
    parser.add_option("-o", action="store", default="MAGG", type="string", dest="output", help="Output directory for the JSON files and the aligned new PWMs (default = 'MAGG')", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
    parser.add_option("--parallel", default=False, action="store_true", dest="parallel", help="Run in parallel if the input is a directory (default = False)")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of JSON files that have failed and have been completed")
    parser.add_option("--complete",default=1.00, action="store", type="float", dest="complete", help="Ratio of completness over the total number of profiles top be done(default= 0.95). This is useful in the server to stop the profiler when the time of execution exceeds more than 48 hours ", metavar="float")
 
    (options, args) = parser.parse_args()

    if options.input_json is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Functions   #
#-------------#

def is_padding(column, threshold=0.0001):
    """
    Determine if a column corresponds to padding (all zeros or near-zero).
    """
    return column.sum() < threshold


def get_trim(x, y):
    """
    Trims tensor x according to padding detected in y.
    x and y are 2D arrays from numpy
    Both x and y are expected to have shape (batch, length).
    Zeros on the left and right of y are treated as padding.
    """
    # Collapse across batch to find positions where at least one element is nonzero
    mask = (y != 0).any(axis=0)  # shape (length,)
    # Find left and right boundaries of nonzero region
    nonzero_indices = np.where(mask)[0]
    if nonzero_indices.size == 0:
        return x[:, :0]  # all zero -> return empty slice
    left, right = nonzero_indices[0], nonzero_indices[-1] + 1
    # Trim x
    return x[:, left:right]


def padding2D(arr,final_length=50,mode="right"):
    """
    Pads or cuts a 2D numpy array (B, L) along the last dimension.
    Args:
        arr (np.ndarray): Input array of shape (B, L).
        final_length (int): Desired length after padding/cutting.
        mode (str): Padding mode: 'left', 'right', or 'center'.
    Returns:
        np.ndarray: Array of shape (B, final_length).
    """
    if arr.ndim != 2:
        raise ValueError("Input must be a 2D array of shape (B, L)")
    batch, length = arr.shape
    # Case 1: Cut if too long
    if length > final_length:
        if mode == "right":  # keep left part
            return arr[:, :final_length]
        elif mode == "left":  # keep right part
            return arr[:, -final_length:]
        elif mode == "center":  # cut equally from both sides
            start = (length - final_length) // 2
            end = start + final_length
            return arr[:, start:end]
        else:
            raise ValueError("mode must be 'left', 'right', or 'center'")
    # Case 2: Pad if too short
    pad_size = final_length - length
    if mode == "right":
        pad_left, pad_right = 0, pad_size
    elif mode == "left":
        pad_left, pad_right = pad_size, 0
    elif mode == "center":
        pad_left = pad_size // 2
        pad_right = pad_size - pad_left
    else:
        raise ValueError("mode must be 'left', 'right', or 'center'")
    return np.pad(arr, ((0, 0), (pad_left, pad_right)), mode="constant", constant_values=0)


def pearson_correlation(a, b,gap_penalty,threshold=0.0001):
    """
    Calculate Pearson correlation between two 1D-tensors.
    """
    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       a_mean = a.mean()
       b_mean = b.mean()
       numerator =((a - a_mean) * (b - b_mean)).sum()
       denominator =np.sqrt(((a - a_mean) ** 2).sum()) * np.sqrt(((b - b_mean) ** 2).sum()) + 1e-8
       return numerator / denominator

def product(a,b,gap_penalty,threshold=0.0001):
    """
    Calculate dot product
    """

    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       return (a * b).sum()

def kullback_leibler(a,b,gap_penalty,threshold=0.0001):
    """
    Calculate Kullback Leibler divergence
    """
    epsilon = 1.0e-12
    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       p = (a + epsilon)/(a + epsilon).sum()
       q = (b + epsilon)/(b + epsilon).sum()
       return 0.5 * (p * (np.log(p) - np.log(q))).sum() + 0.5 * (q * (np.log(q) - np.log(p))).sum()

def jensen_shanon(a,b,gap_penalty,threshold=0.0001):
    """
    Calculate Kullback Leibler divergence
    """
    epsilon = 1.0e-12
    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       p = (a + epsilon)/(a + epsilon).sum()
       q = (b + epsilon)/(b + epsilon).sum()
       m = 0.5 * m
       return 0.5 * (p * (np.log(p) - np.log(m))).sum() + 0.5 * (q * (np.log(q) - np.log(m))).sum()
     
   
def get_metric(choice):
    if choice == "pearson":
        return pearson_correlation
    elif choice == "kullback_leibler":
        return kullback_leibler
    elif choice == "jensen_shanon":
        return jensen_shanon
    elif choice == "product":
        return product
    else:
        raise ValueError(f"Unknown metric: {choice}")


def align_without_gaps(profile1, profile2, metric, gap_penalty=0.5, threshold=0.0001):
    """
    Align two tensor profiles without gaps by attempting all possible offsets.

    Args:
        profile1 (Tensor): first profile of shape (4, L1).
        profile2 (Tensor): Second profile of shape (4, L2).
        gap_penalty (float): Gap penalty for padding.
        threshold (float): minimum sum to be considered non-padding.

    Returns:
        max_score (float): The maximum alignment score.
        max_offset (int): The offset yielding maximum score.
        alignment (list of tuples): List of pairs (i, j) with matching positions.
    """
    L1 = profile1.shape[1]
    L2 = profile2.shape[1]
    max_score = float('-inf')
    max_offset = 0
    max_alignment = []
    # Try all possible offsets
    for offset in range(-L2+1, L1):
        score = 0
        alignment = []
        for i in range(L1):
            j = i - offset
            if j < 0 or j >= L2:
                continue
            col1 = profile1[:, i]
            col2 = profile2[:, j]
            corr  = metric(col1, col2,gap_penalty,threshold)
            score += corr
            alignment.append((i, j))
        if score > max_score:
            max_score = score
            max_offset = offset
            max_alignment = alignment.copy()
    return max_score, max_offset, max_alignment


def build_aligned_profile(profile1, profile2, offset):
    """
    Shift and pad profile2 so it aligns with profile1 based on the offset.

    Args:
        profile1 (np.ndarray): shape (4, L1).
        profile2 (np.ndarray): shape (4, L2).
        offset (int): best offset from align_without_gaps.

    Returns:
        aligned_profile2 (np.ndarray): shape (4, L1),
            profile2 aligned and padded with zeros.
    """
    L1 = profile1.shape[1]
    L2 = profile2.shape[1]
    aligned = np.zeros((profile2.shape[0], L1), dtype=profile2.dtype)
    for i in range(L1):
        j = i - offset
        if 0 <= j < L2:
            aligned[:, i] = profile2[:, j]
    return aligned

def get_pwms(json_file,jaspar,cisbp,hocomoco,modcre,verbose=False,selected_path=None):
    pwms=[]
    models={}
    with open(json_file,"r") as json_fp:
         json_obj=json.load(json_fp)
    for model in json_obj:
        if selected_path is None:
           db=model["database"]
           if db=="jaspar":   path=jaspar
           if db=="cisbp" :   path=cisbp
           if db=="hocomoco": path=hocomoco
           if db=="modcre":   path=modcre
        else:
           path=selected_path
        if os.path.exists(os.path.join(path,model["pwm"])):
           try:
               msa=PWM.nMSA(os.path.join(path,model["pwm"]),None,"meme")
               binding_site_length = msa.get_binding_site_length()
           except Exception as e:
               if verbose: print("\t - Skip PWM %s (error %s)"%(model["pwm"],e))
               continue
           if binding_site_length >0:
               if verbose: print("\t - Add  PWM %s"%(model["pwm"]))
               pwms.append(os.path.join(path,model["pwm"]))
               models.setdefault(os.path.join(path,model["pwm"]),model)
           else:
               if verbose: print("\t - Skip PWM %s (error empty profile )"%(model["pwm"]))
               continue
        else:
           if verbose: print("\t - Skip PWM %s"%(model["pwm"]))
    if len(pwms)>300:
        reduced_pwms=[]
        db_pwms=[]
        for pwm in pwms:
            model=models[pwm]
            if model["database"]=="modcre": 
               reduced_pwms.append(pwm)
            elif float(model["similarity"])>0.50:
               reduced_pwms.append(pwm)
            else:     
               db_pwms.append(pwm)
        size=min(300,len(db_pwms),max(250,300-len(reduced_pwms)))
        reduced_pwms.extend(random.sample(db_pwms,size))
        pwms=reduced_pwms

    return pwms


def get_reference(json_file,jaspar,cisbp,hocomoco,modcre,ratio, verbose=False,selected_path=None):
    pwms=[]
    with open(json_file,"r") as json_fp:
         json_obj=json.load(json_fp)
    for model in json_obj:
        if selected_path is None:
           db=model["database"]
           if db=="jaspar":   path=jaspar
           if db=="cisbp" :   path=cisbp
           if db=="hocomoco": path=hocomoco
           if db=="modcre":   path=modcre
        else:
           path=selected_path
        similarity = float(model["similarity"])
        if  similarity >= ratio and db != "modcre" :
         if os.path.exists(os.path.join(path,model["pwm"])):
           try:
               msa=PWM.nMSA(os.path.join(path,model["pwm"]),None,"meme")
               binding_site_length = msa.get_binding_site_length()
           except Exception as e:
               if verbose: print("\t - Skip PWM %s (error %s)"%(model["pwm"],e))
               continue
           if binding_site_length >0:
               if verbose: print("\t - Add  PWM %s"%(model["pwm"]))
               pwms.append(os.path.join(path,model["pwm"]))
           else:
               if verbose: print("\t - Skip PWM %s (error empty profile )"%(model["pwm"]))
               continue
         else:
           if verbose: print("\t - Skip PWM %s"%(model["pwm"]))
    return pwms



def pwm2tensor(msa_obj):
    x=np.array(msa_obj.get_pwm())
    return x.astype(float).transpose()

def tensor2msa(tensor,name="MOTIF",sequences=False):
    #with torch
    #binding_site_length=tensor.size(-1)
    #pwm=np.array(tensor.transpose(-2,-1))
    #with numpy
    binding_site_length = tensor.shape[1]
    pwm=tensor.transpose()
    msa_obj=PWM.nMSA()
    msa_obj.set_motif(name)
    msa_obj.set_binding_site_length(binding_site_length)
    #print([["%.6f"%(freq) for freq in vector] for vector in pwm])
    msa_obj.set_pwm([["%.6f"%(freq) for freq in vector] for vector in pwm]) 
    if sequences:
       msa_obj.enhance
       msa_obj.set_sequences()
    return msa_obj

def makedb(pwms,folder):
    name=""
    for pwm in pwms:
        meme=os.path.basename(pwm).rstrip("meme").split(":")[0]
        name+=meme[0:len(meme)//2]
    db=os.path.join(folder,"database_"+hashlib.sha224(name.encode()).hexdigest()+".dat")
    if os.path.exists(db): os.remove(db)
    for pwm in pwms:
        os.system("cat %s >> %s"%(pwm,db))
    return db

def get_clusters(matrix,pvalue_threshold=0.005, distance_threshold=0.20):
    from sklearn.cluster import AgglomerativeClustering
    matrix[ np.abs(matrix) <= pvalue_threshold ] = 0
    #maximum=max(max([x for x in row]) for row in matrix)
    #minimum=min(min([x for x in row if x>0]) for row in D if len([x for x in row if x>0])>0)
    #distance_threshold = distance_threshold * max((maximum-minimum),0.1)
    cluster={}
    clustering = AgglomerativeClustering(
      n_clusters=None,  # no fixed number of clusters
      metric='precomputed',  # we provide distance matrix directly (instead of a similarity matrix)
      #affinity='precomputed', # Metric used to compute the linkage. Can be "euclidean", "l1", "l2", "manhattan", "cosine", or "precomputed (precomputed is a distance matrix instead of a similarity matrix)
      linkage='average',
      distance_threshold=distance_threshold  # tune this threshold to get clusters
    )
    labels = clustering.fit_predict(matrix)
    n_clusters = labels.max() + 1
    for cluster_id in range(n_clusters):
        members = np.where(labels == cluster_id)[0]
        min_sum_dist = float('inf')
        central_element = None
        for i in members:
            # Sum distance to other cluster members
            sum_dist = sum(matrix[i, j] for j in members if i != j)
            if sum_dist < min_sum_dist:
               min_sum_dist = sum_dist
               central_element = i
        cluster.setdefault(central_element,members)
    return cluster

def reverse_tensor(tensor,motif_name):
    msa_obj = tensor2msa(tensor,motif_name,False)
    msa_rev = msa_obj.get_complementary()
    reverse = pwm2tensor(msa_rev)
    return reverse

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    input_json   = options.input_json
    binding_site = int(options.length)
    jaspar       = options.jaspar
    cisbp        = options.cisbp
    hocomoco     = options.hocomoco
    modcre       = options.modcre
    dummy        = options.dummy_dir
    outdir       = options.output
    pv_threshold = float(options.pvalue)
    threshold    = float(options.threshold)
    gap_penalty  = float(options.gap_penalty)
    reference    = options.reference
    trim         = options.trim     
    keep         = options.keep
    verbose      = options.verbose
    parallel     = options.parallel
    info_file    = options.info
    complete     = float(options.complete)
   
    metric_choice = "product"
    if options.pearson: metric_choice="pearson"
    if options.kullback_leibler: metric_choice="kullback_leibler"
    if options.jensen_shanon: metric_choice="jensen_shanon"

    if not os.path.exists(dummy):  os.makedirs(dummy)
    if not os.path.exists(outdir): os.makedirs(outdir)
    list_of_jsons=[]
    name_of_tfs=[]
    if os.path.isdir(input_json):
        for json_file in os.listdir(input_json):
            if json_file.endswith("json"): 
               list_of_jsons.append(os.path.join(input_json,json_file))
               name_of_tfs.append(json_file.rstrip(".json"))
    else:
        list_of_jsons.append(input_json)
        name_of_tfs.append(input_json.rstrip(".json"))
    # Iterate untill all protein profiles are done
    submitted=set()
    n_done=0
    if verbose: print("Start iteration to check and run profiles")
    if info_file is None: info_file=os.path.join(outdir,"json_list_execution.log")
    if verbose: print(("Open to write Information file %s"%info_file))
    if not os.path.exists(info_file):
       log_file = open(info_file,"w")
       log_file.write("#List of JSONs\n")
       log_file.close()
    done=functions.done_jobs(info_file)
    iterate = functions.check_done(done,name_of_tfs)
    maxtime = 3600 * 3
    start_time = d_time = 0
    while( iterate ):
     #Analyze each JSON file
     for json_file in list_of_jsons:
       if json_file in submitted: continue
       if parallel:
            skip=False
            if os.path.exists(os.path.join(outdir,os.path.basename(json_file.rstrip(".json")))):
               if os.path.exists(os.path.join(outdir,os.path.basename(json_file.rstrip(".json")),"matrix_of_pvalues.dat")):
                  if len(os.listdir(os.path.join(outdir,os.path.basename(json_file.rstrip(".json"))))) > 1:
                     skip=True
            if skip: 
               if verbose: print("\t -JSON file %s is already done...."%(os.path.basename(json_file)))
               submitted.add(json_file)
               log_file=open(info_file,"a")
               log_file.write("%s\tDONE\n"%(json_file.rstrip(".json")))
               log_file.flush()
               log_file.close()
               continue
            #Make input
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            program    =  os.path.join(exe_path,"aggregate_pwms.py")
            python     =  os.path.join(config.get("Paths", "python_path"), "python")
            parameters =  " -i %s "%json_file 
            parameters =  parameters + " -o %s "%outdir
            parameters =  parameters + " --dummy=%s "%dummy
            parameters =  parameters + " --info=%s "%info_file
            parameters =  parameters + " --binding_site=%s "%str(binding_site)
            parameters =  parameters + " --pvalue=%s "%str(pv_threshold)
            parameters =  parameters + " --threshold=%s "%str(threshold)
            parameters =  parameters + " --gap_penalty=%s "%str(gap_penalty)
            if trim:      parameters = parameters +  " --trim "
            if keep:      parameters = parameters +  " --keep "
            if verbose:   parameters = parameters +  " --verbose "
            if options.kullback_leibler: parameters = parameters +  " --kullback_leibler "
            if options.pearson:          parameters = parameters +  " --pearson "
            if options.jensen_shanon:    parameters = parameters +  " --jensen_shanon "
            if jaspar is not None: parameters =  parameters + " --jaspar=%s "%jaspar
            if cisbp is not None: parameters =  parameters + " --cisbp=%s "%cisbp
            if hocomoco is not None: parameters =  parameters + " --hocomoco=%s "%hocomoco
            if modcre is not None: parameters =  parameters + " --modcre=%s "%modcre
            if reference is not None: parameters =  parameters + "--reference=%s "%str(reference)
            # Execute in cluster
            if verbose: print("\t %s %s %s\n" % (python,program,parameters))
            functions.submit_command_to_queue("%s %s %s" % (python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            submitted.add(json_file)

       else:
         try:
            output = os.path.join(outdir,os.path.basename(json_file).rstrip(".json"))
            if not os.path.exists(output): os.makedirs(output)

            #0) Get PWMs
            name = ".".join(os.path.basename(json_file).split(".")[:-1])
            if verbose: print("\t - Get PWMs of JSON %s"%(name))
            pwms = get_pwms(json_file,jaspar,cisbp,hocomoco,modcre,verbose)
        
            #1) Make DB
            if verbose:print("\t - Make DB of PWMs") 
            db=makedb(pwms,dummy)

            #2) Compare PWM vs DB
            if verbose:print("\t - Get pairing p-values") 
            adjacency_pvalues={}
            adjacency_strand={}
            pwm_names=[]
            pwm_names_dict={}
            for pwm in pwms:
              if verbose:print("\t\tGet tomtom of ",pwm)
              try:
                tomtom_obj = tomtom.get_tomtom_obj(db,pwm,dummy)
                #get p-values 
                query=tomtom_obj.get_query()
                if query is not None:
                   pwm_names.append(query.rstrip(".txt"))
                   pwm_names_dict.setdefault(pwm,len(pwm_names)-1)
                   for hit_obj in tomtom_obj.get_hits():
                      hit    = hit_obj.get_hit()
                      pvalue = hit_obj.get_p_value()
                      strand = hit_obj.get_strand()
                      if pvalue<pv_threshold and verbose: print("\t\t   Add  hit ", hit, " with p-value = ",pvalue)
                      adjacency_pvalues.setdefault((query.rstrip(".txt"),hit.rstrip(".txt")),pvalue)
                      adjacency_pvalues.setdefault((hit.rstrip(".txt"),query.rstrip(".txt")),pvalue)
                      adjacency_strand.setdefault((query.rstrip(".txt"),hit.rstrip(".txt")),strand)
                      adjacency_strand.setdefault((hit.rstrip(".txt"),query.rstrip(".txt")),strand)
                else:
                   print("\t\t\t- TOMTOM Failed query=None") 
                   pwms.remove(pwm)
              except Exception as e:
                print("\t\t\t- TOMTOM Failed ")
                pwms.remove(pwm)
            #3) Get comparison matrix
            if verbose:print("\t - Get Similarity Matrix") 
            comparison_matrix=[]
            for i in range(len(pwm_names)):
              vector=[]
              pwm_a = pwm_names[i]
              for j in range(len(pwm_names)):
                pwm_b = pwm_names[j]
                if (pwm_a,pwm_b) in adjacency_pvalues:
                   vector.append(adjacency_pvalues[(pwm_a,pwm_b)])
                else:
                   vector.append(1.0) 
              comparison_matrix.append(vector)
            comparison_matrix = np.array(comparison_matrix)

            #4) Get clusters
            if verbose:print("\t - Get Clusters of PWMs") 
            cluster = get_clusters(comparison_matrix,pv_threshold,threshold)
            membership={}
            cluster_number={}
            n=0
            for center in cluster:
              n=n+1
              if verbose:print("\t\t - Cluster %d ( PWM %d ) : %s"%(n, center,str(cluster[center])))
              for i in cluster[center]:
                membership.setdefault(i,center)
              cluster_number.setdefault(center,str(n))

            #5) Make output of TOMTOM p-values
            if verbose: print("\t - Generate output of TOMTOM matrix p-values")
            fo=open(os.path.join(output,"matrix_of_pvalues.dat"),"w")
            for i in range(len(pwm_names)):
               fo.write(" %d => %s\n"%(i,pwm_names[i]))
            fo.write("====  MATRIX  =====\n")
            for row in comparison_matrix:
               for pv in row:
                fo.write("\t%.6f "%(pv))
               fo.write("\n")
            fo.write("==== CLUSTERS =====\n")
            for center in cluster:
               fo.write(" Cluster %s  %s \n"%(cluster_number[center],str(cluster[center])))
            fo.close()

            #6) Get tensors of PWMs
            if verbose:print("\t - Get tensor-profiles for each PWM") 
            pwm_tensor=[]
            for pwm in pwms:
               msa_obj = PWM.nMSA(pwm,None,"meme")
               pwm_tensor.append(pwm2tensor(msa_obj))

            #7) Get aligned tensors per cluster
            aligned_pwm_tensor={}
            if verbose:print("\t - Align the tensors of each cluster") 
            for center in cluster:
              if verbose:print("\t\t - Alignment of cluster %s"%cluster_number[center]) 
              profile1= padding2D(pwm_tensor[center],final_length=binding_site,mode="center")
              aligned_pwm_tensor.setdefault(center,profile1)
              pwm_a = pwm_names[center]
              for i in cluster[center]:
                try:
                  pwm_b = pwm_names[i]
                  if (pwm_a,pwm_b) in adjacency_strand:
                    if adjacency_strand[(pwm_a,pwm_b)] == "-": profile2 = reverse_tensor(pwm_tensor[i],pwm_b)
                    if adjacency_strand[(pwm_a,pwm_b)] == "+": profile2 = pwm_tensor[i]
                  else:
                    if (pwm_b,pwm_a) in adjacency_strand:
                      if adjacency_strand[(pwm_b,pwm_a)] == "-": profile2 = reverse_tensor(pwm_tensor[i],pwm_b)
                      if adjacency_strand[(pwm_b,pwm_a)] == "+": profile2 = pwm_tensor[i]
                    else:
                      print("\t\t\t - Skip PROFILES %s %s without known orientation"%(pwm_a,pwm_b))
                      continue
                  metric=get_metric(metric_choice)
                  score, offset, alignment = align_without_gaps(profile1, profile2, metric, gap_penalty=gap_penalty)
                  if verbose: 
                      if (pwm_a,pwm_b) in adjacency_pvalues: print("\t\t\t - PROFILES %s - %s SCORE %.6f Offset %d Strand %s P-value %.6f (%.6f)"%(pwm_a,pwm_b,score,offset,adjacency_strand[(pwm_a,pwm_b)],adjacency_pvalues[(pwm_b,pwm_a)],comparison_matrix[center,i]))
                      elif (pwm_b,pwm_a) in adjacency_pvalues: print("\t\t\t - PROFILES %s - %s SCORE %.6f Offset %d Strand %s P-value %.6f (%.6f)"%(pwm_b,pwm_a,score,offset,adjacency_strand[(pwm_b,pwm_a)],adjacency_pvalues[(pwm_b,pwm_a)],comparison_matrix[center,i]))
                      else: print("\t\t\t - PROFILES %s - %s SCORE %.6f Offset %d P-value ?? (%.6f)"%(pwm_b,pwm_a,score,offset,comparison_matrix[center,i]))
                  aligned_profile2 = build_aligned_profile(profile1, profile2, offset)
                  aligned_pwm_tensor.setdefault(i,aligned_profile2)
                except Exception as e:
                  print("Failed to align %s vs %s (error %s)"%(pwm_a,pwm_b,e))
                 
            #8) Get PWMS from tensors
            msas={}
            if verbose:print("\t - Get PWMs of aligned tensors") 
            for i in aligned_pwm_tensor:
               motif_name = pwm_names[i]
               tensor     = aligned_pwm_tensor[i]
               msa_obj    = tensor2msa(tensor,motif_name,sequences=True)
               msas.setdefault(i,msa_obj)

            #9) Get the cluster-averaged PWM
            if verbose:print("\t - Get Averaged PWMs per cluster") 
            for center in cluster:
              if not os.path.exists(os.path.join(output,"cluster_"+cluster_number[center])):
                 os.makedirs(os.path.join(output,"cluster_"+cluster_number[center]))
              if verbose:print("\t\t - PWM average of cluster %s"%center) 
              pwm_list = []
              msa_obj  = msas[center]
              for i in cluster[center]:
                if i in msas: pwm_list.append(msas[i])
              msa_comb  = msa_obj.combine(pwm_list)
              msa_comb.set_motif(msa_obj.get_motif()+"_average_cluster_"+cluster_number[center])
              msa_comb.set_sequences()
              file_name = os.path.join(output,"cluster_"+cluster_number[center],msa_comb.get_motif())
              msa_comb.write(file_name+".meme","meme",True)
              if keep: msa_comb.write(file_name+".msa","msa",True)
              if keep: msa_comb.write(file_name+".pwm","pwm",True)
              PWM.write_logo(msa_comb,file_name+".logo",dummy)

            #10) Make output PWMs
            if verbose:print("\t - Generate PWMs and LOGOS of all aligned PWMs") 
            for i in msas:
              msa_obj = msas[i]
              if not os.path.exists(os.path.join(output,"cluster_"+cluster_number[membership[i]])):
                 os.makedirs(os.path.join(output,"cluster_"+cluster_number[membership[i]]))
              file_name = os.path.join(output,"cluster_"+cluster_number[membership[i]],msa_obj.get_motif())+"_member_of_cluster_"+cluster_number[membership[i]]
              msa_obj.write(file_name+".meme","meme",True)
              if keep: msa_obj.write(file_name+".msa","msa",True)
              if keep: msa_obj.write(file_name+".pwm","pwm",True)
              PWM.write_logo(msa_obj,file_name+".logo",dummy)

            #11) Select reference PWMs from JSON
            if reference is not None:
              if verbose: print("\t - Get REFERENCE PWMs of JSON %s"%(name))
              ratio = float(reference)
              reference_pwms = get_reference(json_file,jaspar,cisbp,hocomoco,modcre,ratio,verbose)
              reference_number=[]
              reference_names=[]
              for pwm in reference_pwms:
                reference_number.append(pwm_names_dict[pwm])
                reference_names.append(pwm_names[pwm_names_dict[pwm]])
              if verbose: print("\t - PWMS of reference are: %s"%(str(reference_names)))
               
            #12) Realign clusters with references 
            if reference is not None:
              realigned_pwm_tensor={}
              for pwm in reference_pwms:
                pwm_reference_number = pwm_names_dict[pwm]
                pwm_reference_name   = pwm_names[pwm_reference_number]
                if verbose:print("\t - Arrange alignments for reference %s "%(pwm_reference_name))
                profile1= padding2D(pwm_tensor[pwm_reference_number],final_length=binding_site,mode="center")
                if trim: realigned_pwm_tensor.setdefault((pwm_reference_number,pwm_reference_name),get_trim(profile1,profile1))
                else:    realigned_pwm_tensor.setdefault((pwm_reference_number,pwm_reference_name),profile1)
                for center in cluster:  
                   min_idx_cluster = int(np.argmin(np.array([ adjacency_pvalues[(pwm_reference_name,pwm_names[index])] for index in cluster[center]]) ))
                   min_idx         = cluster[center][min_idx_cluster]
                   pwm_min         = pwm_names[min_idx] 
                   if   adjacency_pvalues[(pwm_reference_name,pwm_min)] <  pv_threshold  :
                        if verbose:print("\t\t - Alignment of reference %s for cluster %s (best match %s)"%(pwm_reference_name,cluster_number[center],pwm_min)) 
                        try:
                          pwm_c = pwm_names[center]
                          pwm_a = pwm_reference_name
                          pwm_b = pwm_min
                          if (pwm_reference_name,pwm_min) in adjacency_strand and (pwm_min,pwm_c) in adjacency_strand:
                              pair_1 = (pwm_a,pwm_b)
                              pair_2 = (pwm_b,pwm_c)
                          else:
                              if verbose:print("\t\t - SKIP PROFILES %s %s %s without known orientation"%(pwm_a,pwm_b,pwm_c))
                              continue
                          reverse=False
                          if adjacency_strand[pair_1] == "-" and adjacency_strand[pair_2] == "-" : reverse=False
                          if adjacency_strand[pair_1] == "+" and adjacency_strand[pair_2] == "-" : reverse=True
                          if adjacency_strand[pair_1] == "-" and adjacency_strand[pair_2] == "+" : reverse=True
                          if adjacency_strand[pair_1] == "+" and adjacency_strand[pair_2] == "+" : reverse=False
                          if reverse: profile2 = reverse_tensor(aligned_pwm_tensor[min_idx],pwm_min)
                          else:       profile2 = aligned_pwm_tensor[min_idx]
                          metric=get_metric(metric_choice)
                          score, offset, alignment = align_without_gaps(profile1, profile2, metric, gap_penalty=gap_penalty)
                          if verbose : print("\t\t\t - PROFILES %s - %s SCORE %.6f Offset %d "%(pwm_reference_name,pwm_min,score,offset))
                          for i in cluster[center]:
                              if reverse: profile_i = reverse_tensor(aligned_pwm_tensor[i],pwm_names[i]) 
                              else:       profile_i = aligned_pwm_tensor[i]
                              realigned_profile_i = build_aligned_profile(profile1, profile_i, offset)
                              if trim:  realigned_pwm_tensor.setdefault((i,pwm_reference_name),get_trim(realigned_profile_i,profile1))
                              else:     realigned_pwm_tensor.setdefault((i,pwm_reference_name),realigned_profile_i)
                        except Exception as e:
                          print("Failed to align %s vs %s (error %s)"%(pwm_reference_name,pwm_min,e))


            #13) Get realigned PWMS from tensors
            if reference is not None:
              reference_msas={}
              if verbose:print("\t - Get PWMs of re-aligned tensors with respect to reference PWMs") 
              for i,ref in realigned_pwm_tensor:
                motif_name = pwm_names[i]
                tensor     = realigned_pwm_tensor[(i,ref)]
                msa_obj    = tensor2msa(tensor,motif_name,sequences=True)
                reference_msas.setdefault((i,ref),msa_obj)

            #14) Get the cluster-aggregate of realigned PWMs referred to a reference
            if reference is not None:
               if verbose:print("\t - Get Averaged realigned PWMs per reference and cluster") 
               pwm_list_by_reference={}
               for i,ref in reference_msas:
                   msa_obj = reference_msas[i,ref]
                   if pwm_names[i] != ref:
                      center=membership[i]
                      pwm_list_by_reference.setdefault((ref,center),[]).append(msa_obj)
               for ref,center in pwm_list_by_reference:
                   pwm_list =  pwm_list_by_reference[(ref,center)]
                   msa_obj  =  pwm_list[0]
                   msa_comb  = msa_obj.combine(pwm_list)
                   msa_comb.set_motif(pwm_names[center]+"_average_cluster_"+cluster_number[center]+"_referred_to_"+ref)
                   msa_comb.set_sequences()
                   if not os.path.exists(os.path.join(output,"reference_"+ref)):
                      os.makedirs(os.path.join(output,"reference_"+ref))
                   if not os.path.exists(os.path.join(output,"reference_"+ref,"cluster_"+cluster_number[center])):
                      os.makedirs(os.path.join(output,"reference_"+ref,"cluster_"+cluster_number[center]))
                   file_name = os.path.join(output,"reference_"+ref,"cluster_"+cluster_number[center],msa_comb.get_motif())
                   msa_comb.write(file_name+".meme","meme",True)
                   if keep: msa_comb.write(file_name+".msa","msa",True)
                   if keep: msa_comb.write(file_name+".pwm","pwm",True)
                   PWM.write_logo(msa_comb,file_name+".logo",dummy)

            #15) Make output of realigned PWMs
            if reference is not None:
              if verbose:print("\t - Generate PWMs and LOGOS of all re-aligned PWMs") 
              for i,ref in reference_msas:
                msa_obj = reference_msas[i,ref]
                if not os.path.exists(os.path.join(output,"reference_"+ref)):
                   os.makedirs(os.path.join(output,"reference_"+ref))
                if not os.path.exists(os.path.join(output,"reference_"+ref,"cluster_"+cluster_number[membership[i]])):
                   os.makedirs(os.path.join(output,"reference_"+ref,"cluster_"+cluster_number[membership[i]]))
                file_name = os.path.join(output,"reference_"+ref,"cluster_"+cluster_number[membership[i]],msa_obj.get_motif()+"_member_of_cluster_"+cluster_number[membership[i]]+"_referred_to_"+ref)
                msa_obj.write(file_name+".meme","meme",True)
                if keep: msa_obj.write(file_name+".msa","msa",True)
                if keep: msa_obj.write(file_name+".pwm","pwm",True)
                PWM.write_logo(msa_obj,file_name+".logo",dummy)

            #16) Done with a json file input  
            log_file=open(info_file,"a")
            log_file.write("%s\tDONE\n"%(json_file.rstrip(".json")))
            log_file.flush()
            log_file.close()
         except Exception as e:
            #16)Failed a json file input  
            if verbose:print("\t\t - Failed (Error %s)"%e) 
            log_file=open(info_file,"a")
            log_file.write("%s\tFAIL\n"%(json_file.rstrip(".json")))
            log_file.flush()
            log_file.close()

         submitted.add(json_file)
 

     #Check next iteration, profiles submitted and profiles done
     done=functions.done_jobs(info_file)
     iterate= functions.check_done(done,name_of_tfs) 
     if len(done) > n_done:
           n_done=len(done)
           #if n_done>1 and start_time==0: start_time = float(time.time())
           if len(submitted)==len(name_of_tfs) and start_time==0: start_time = float(time.time())
           if options.verbose: 
                sys.stdout.write("Number of files already done %d\n"%n_done)
                if d_time>0: sys.stdout.write("Time: %f\n"%d_time)
                sys.stdout.write("\t-- Check files done ...\n")
                log_file = open(info_file,"r")
                for line in log_file:
                   print(("\t\t-- %s"%line.strip()))
                log_file.flush()
                log_file.close()
                sys.stdout.write("\t-- Still running protein profiles  %s ...\n"%functions.check_done(done,name_of_tfs))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush() 
     #Check next iteration, if exceeding time and enough profiles stop iteration
     if start_time>0:
          current_time = float(time.time())
          d_time = current_time - start_time
          if float(  len(done) ) / len(name_of_tfs) > complete  and d_time > maxtime: iterate=False
          if d_time > maxtime and iterate and complete < 1.0: 
               complete = complete - 0.01
               maxtime  = maxtime  + 10
               sys.stdout.write("Time: %f Done: %d (%d) Ratio to end: %f\n"%(d_time,len(done),len(name_of_tfs),complete))

    
    #Exit    
    if not verbose:
        shutil.rmtree(dummy)
    print("Done")

