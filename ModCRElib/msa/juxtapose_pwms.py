"""
Cluster, align, and average related PWMs described in JSON manifests.

This script collects motif files from one or more PWM databases, compares them
with TOMTOM, clusters similar motifs, aligns the members of each cluster as
tensors, and writes both cluster-average motifs and aligned member motifs.
"""

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
    parser = optparse.OptionParser("python juxtapose_pwms.py -i input_json [ -l length -o output_name --dummy dummy_dir  ]  ")
    parser.add_option("-i", action="store", default=None, type="string", dest="input_json", help="Input JSON file or directory of JSON files with information of PWMs", metavar="{filename}")
    parser.add_option("-l","--binding_site", action="store", default=50, type="int", dest="length", help="Length of the binding site of the output PWMs", metavar="{integer}")
    parser.add_option("--pvalue", action="store", default=0.005,type="float", dest="pvalue", help="P-value threshold of TOMTOM similarity between two PWMs (default 0.005)", metavar="{float}")
    #parser.add_option("--threshold", action="store", default=0.25,type="float", dest="threshold", help="Ratio of the maximum distance between dissimilar PWMs to use as distance threshold of the agglomerative clustering (default is 0.25)", metavar="{float}")
    parser.add_option("--threshold", action="store", default=0.2,type="float", dest="threshold", help="Distance threshold of the agglomerative clustering (default is 0.2)", metavar="{float}")
    parser.add_option("--jaspar", action="store", default=None, type="string", dest="jaspar", help="Address of JASPAR PWMs", metavar="{filename}")
    parser.add_option("--cisbp", action="store", default=None, type="string", dest="cisbp", help="Address of CisBP PWMs", metavar="{filename}")
    parser.add_option("--hocomoco", action="store", default=None, type="string", dest="hocomoco", help="Address of HOCOMOCO PWMs", metavar="{filename}")
    parser.add_option("--modcre", action="store", default=None, type="string", dest="modcre", help="Address of predicted PWMs with ModCRE", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("--kullback_leibler", default=False, action="store_true", dest="kullback_leibler", help="Flag to select  Kullback Leibler divergence to compare profiles (default is scalar product)", metavar="{boolean}")
    parser.add_option("--pearson", default=False, action="store_true", dest="pearson", help="Flag to select Pearson correlation to compare profiles (default is scalar product)", metavar="{boolean}")
    parser.add_option("-o", action="store", default="MJUX", type="string", dest="output", help="Output directory for the JSON files and the aligned new PWMs (default = 'MJUX')", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
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

    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       return (a * (a.log() - b.log())).sum()
    
   
def get_metric(choice):
    """
    Select the column-comparison metric used during PWM alignment.

    Args:
        choice (str): Metric name. Supported values are ``pearson``,
            ``kullback_leibler``, and ``product``.

    Returns:
        callable: Function accepting two PWM columns and returning a similarity
        score.
    """
    if choice == "pearson":
        return pearson_correlation
    elif choice == "kullback_leibler":
        return kullback_leibler
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
    """
    Resolve motif files listed in a JSON manifest.

    Args:
        json_file (str): JSON file describing motifs and their source database.
        jaspar, cisbp, hocomoco, modcre (str | None): Base directories for each
            supported PWM collection.
        verbose (bool): If True, print skip/add messages.
        selected_path (str | None): Optional explicit directory overriding the
            database-specific path from the JSON metadata.

    Returns:
        list[str]: Existing PWM file paths that could be loaded successfully.
    """
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
        if os.path.exists(os.path.join(path,model["pwm"])):
           try:
               msa=PWM.nMSA(os.path.join(path,model["pwm"]),None,"meme")
               binding_site_length = msa.get_binding_site_length()
           except Exception as e:
               if verbose: print("\t - Skip PWM %s (error %s)"%(model["pwm"],e))
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
    """
    Convert an :class:`nMSA` PWM into a ``(alphabet, length)`` NumPy tensor.

    Args:
        msa_obj (PWM.nMSA): Motif object containing a PWM.

    Returns:
        np.ndarray: Float array with shape ``(4, L)`` for DNA motifs.
    """
    x=np.array(msa_obj.get_pwm())
    return x.astype(float).transpose()

def tensor2msa(tensor,name="MOTIF"):
    #with torch
    #binding_site_length=tensor.size(-1)
    #pwm=np.array(tensor.transpose(-2,-1))
    #with numpy
    """
    Convert a PWM tensor back into an :class:`nMSA` object.

    Args:
        tensor (np.ndarray): PWM-like array of shape ``(4, L)``.
        name (str): Motif name assigned to the returned object.

    Returns:
        PWM.nMSA: Motif object with the tensor stored as PWM and synthetic
        sequences generated from it.
    """
    binding_site_length = tensor.shape[1]
    pwm=tensor.transpose()
    msa_obj=PWM.nMSA()
    msa_obj.set_motif(name)
    msa_obj.set_binding_site_length(binding_site_length)
    msa_obj.set_pwm([["%.6f"%(freq) for freq in vector] for vector in pwm]) 
    msa_obj.enhance
    msa_obj.set_sequences()
    return msa_obj

def makedb(pwms,folder):
    """
    Concatenate multiple MEME files into a temporary TOMTOM database.

    Args:
        pwms (list[str]): Motif file paths to concatenate.
        folder (str): Directory where the temporary database file is created.

    Returns:
        str: Path to the generated TOMTOM database file.
    """
    name=""
    for pwm in pwms:
        meme=os.path.basename(pwm).rstrip("meme").split(":")[0]
        name+=meme[0:len(meme)//2]
    db=os.path.join(folder,"database_"+hashlib.sha224(name.encode()).hexdigest()+".dat")
    if os.path.exists(db): os.remove(db)
    for pwm in pwms:
        os.system("cat %s >> %s"%(pwm,db))
    return db

def get_clusters(matrix,pvalue_similarity=0.005, distance_threshold=0.20):
    """
    Cluster motifs from a TOMTOM-derived distance matrix.

    Args:
        matrix (np.ndarray): Pairwise distance/similarity matrix.
        pvalue_similarity (float): Values at or below this threshold are forced
            to zero before clustering.
        distance_threshold (float): Agglomerative clustering cut-off.

    Returns:
        dict: Mapping ``cluster_center_index -> member_indices``.
    """
    from sklearn.cluster import AgglomerativeClustering
    matrix[ np.abs(matrix) <= pvalue_similarity ] = 0
    #maximum=max(max([x for x in row]) for row in matrix)
    #minimum=min(min([x for x in row if x>0]) for row in D if len([x for x in row if x>0])>0)
    #distance_threshold = ratio_threshold * max((maximum-minimum),0.1)
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


def main():
    """
    Run the PWM juxtaposition workflow from the command line.

    Workflow:
        1. Load one or more JSON manifests describing PWM collections.
        2. Resolve motif files and build a temporary TOMTOM database.
        3. Compute pairwise motif similarities and cluster related PWMs.
        4. Align the PWMs in each cluster as tensors.
        5. Write aligned member motifs and cluster-average motifs to disk.

    Returns:
        None. Output files are written into the selected output directory.
    """
    # Arguments & Options #
    options = parse_options()
    input_json   = options.input_json
    binding_site = options.length
    jaspar       = options.jaspar
    cisbp        = options.cisbp
    hocomoco     = options.hocomoco
    modcre       = options.modcre
    dummy        = options.dummy_dir
    outdir       = options.output
    similarity   = options.pvalue
    threshold    = options.threshold
    verbose      = options.verbose
   
    metric_choice = "product"
    if options.pearson: metric_choice="pearson"
    if options.kullback_leibler: metric_choice="kullback_leibler"

    if not os.path.exists(dummy):  os.makedirs(dummy)
    if not os.path.exists(outdir): os.makedirs(outdir)
    list_of_jsons=[]
    if os.path.isdir(input_json):
        for json_file in os.listdir(input_json):
            if json_file.endswith("json"): list_of_jsons.append(os.path.join(input_json,json_file))
    else:
        list_of_jsons.append(input_json)

    for json_file in list_of_jsons:

        if len(list_of_jsons)>0:
           output = os.path.join(outdir,os.path.basename(json_file).rstrip(".json"))
           if not os.path.exists(output): os.makedirs(output)
        else:
           output=outdir

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
        for pwm in pwms:
            if verbose:print("\t\tGet tomtom of ",pwm)
            tomtom_obj = tomtom.get_tomtom_obj(db,pwm,dummy)
            #get p-values 
            query=tomtom_obj.get_query()
            pwm_names.append(query.rstrip(".txt"))
            for hit_obj in tomtom_obj.get_hits():
                hit    = hit_obj.get_hit()
                pvalue = hit_obj.get_p_value()
                strand = hit_obj.get_strand()
                if pvalue<similarity and verbose: print("\t\t   Add  hit ", hit, " with p-value = ",pvalue)
                adjacency_pvalues.setdefault((query.rstrip(".txt"),hit.rstrip(".txt")),pvalue)
                adjacency_pvalues.setdefault((hit.rstrip(".txt"),query.rstrip(".txt")),pvalue)
                adjacency_strand.setdefault((query.rstrip(".txt"),hit.rstrip(".txt")),strand)
                adjacency_strand.setdefault((hit.rstrip(".txt"),query.rstrip(".txt")),strand)
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
        cluster = get_clusters(comparison_matrix,similarity,threshold)
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
        pwm_tensor_reverse=[]
        for pwm in pwms:
            msa_obj = PWM.nMSA(pwm,None,"meme")
            pwm_tensor.append(pwm2tensor(msa_obj))
            msa_rev = msa_obj.get_complementary()
            pwm_tensor_reverse.append(pwm2tensor(msa_rev))

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
                    if adjacency_strand[(pwm_a,pwm_b)] == "-": profile2 = pwm_tensor_reverse[i]
                    if adjacency_strand[(pwm_a,pwm_b)] == "+": profile2 = pwm_tensor[i]
                  else:
                    if (pwm_b,pwm_a) in adjacency_strand:
                      if adjacency_strand[(pwm_b,pwm_a)] == "-": profile2 = pwm_tensor_reverse[i]
                      if adjacency_strand[(pwm_b,pwm_a)] == "+": profile2 = pwm_tensor[i]
                    else:
                      print("\t\t\t - Skip PROFILES %s %s without known orientation"%(pwm_a,pwm_b))
                      continue
                  metric=get_metric(metric_choice)
                  score, offset, alignment = align_without_gaps(profile1, profile2, metric, gap_penalty=0.05)
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
            msa_obj    = tensor2msa(tensor,motif_name)
            msas.setdefault(i,msa_obj)

        #9) Get the cluster-averaged PWM
        if verbose:print("\t - Get Averaged PWMs per cluster") 
        for center in cluster:
            if verbose:print("\t\t - PWM average of cluster %s"%center) 
            pwm_list = []
            msa_obj  = msas[center]
            for i in cluster[center]:
                pwm_list.append(msas[i])
            msa_comb  = msa_obj.combine(pwm_list)
            msa_comb.set_motif(msa_obj.get_motif()+"_average_cluster_"+cluster_number[center])
            msa_comb.set_sequences()
            file_name = os.path.join(output,msa_comb.get_motif())
            msa_comb.write(file_name+".meme","meme",True)
            msa_comb.write(file_name+".msa","msa",True)
            msa_comb.write(file_name+".pwm","pwm",True)
            PWM.write_logo(msa_comb,file_name+".logo",dummy)

        #10) Make output PWMs
        if verbose:print("\t - Generate PWMs and LOGOS of all aligned PWMs") 
        for i in msas:
            msa_obj = msas[i]
            file_name = os.path.join(output,msa_obj.get_motif())+"_member_of_cluster_"+cluster_number[membership[i]]
            msa_obj.write(file_name+".meme","meme",True)
            msa_obj.write(file_name+".msa","msa",True)
            msa_obj.write(file_name+".pwm","pwm",True)
            PWM.write_logo(msa_obj,file_name+".logo",dummy)
           
    
    #11) Exit    
    if not verbose:
        shutil.rmtree(dummy)
    print("Done")


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
