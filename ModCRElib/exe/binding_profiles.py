"""
Binding profiles for ModCRE: PWM + FIMO scans, statistical energies, and clustering.

Given a directory of TF–DNA models, the module builds position weight matrices with
``pwm_pbm``, runs FIMO on a target DNA FASTA, scores hits with statistical potentials,
clusters motifs, and exports per-model and per-cluster JSON (plus optional matplotlib
and Plotly HTML). Helper functions support merging tracks, windowing, and web-facing
plots.

How to run (CLI):
    ``python binding_profiles.py --dummy DUMMY_DIR -i MODEL_DIR -d DNA.fa \\
        --pdb PDB_DIR --pbm PBM_DIR [-o OUTPUT_DIR] [-l LABEL] [-v]``

    Entry point is ``main()``; the same helpers are importable for custom pipelines.

Outputs (under ``output_dir``):
    ``pwms/`` (MEME + MSA per motif and potential label), ``individual_profiles/``
    (JSON per model), ``cluster_profiles/`` (merged JSON per cluster), and
    ``clusters.txt`` (cluster membership). Optional Plotly/matplotlib helpers write
    beside or under ``plots/`` when those functions are used directly.

See also:
    ``get_potentials_and_thresholds()`` — table parser for optimal
    potential/threshold rows (for workflows that drive ``load_statistical_potentials``
    from a grid-search file rather than the fixed ``pair`` / ``s3dc_dd`` loop in
    ``main()``).
"""
import os, sys, re
import configparser
import optparse
import shutil
import subprocess
import numpy as np
import datetime
import pickle
#import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
import itertools
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes, text

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


# Imports jbonet's module #
from SBILib.structure import PDB
from SBILib.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBILib.structure.residue import ResidueOfNucleotide
from SBILib.structure.atom import AtomOfNucleotide

# Import my modules #
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.structure.protein import dssp
from ModCRElib.structure.dna import x3dna, model_dna
from ModCRElib.msa import  pwm_pbm
from ModCRElib.structure.threading import  threader
from ModCRElib.web import cluster_motifs
from ModCRElib.profile import fimo, scorer
from ModCRElib.sequence import blast



def get_potentials_and_thresholds(op_file):
    """
    Read split-potential names and score thresholds from a grid-search results file.

    Expects a header line plus data rows with whitespace-separated fields; field
    index 2 is the threshold and index 3 is the split-potential key (legacy layout).

    Args:
        op_file (str): Path to the optimal-potentials table.

    Returns:
        list[tuple]: ``(split_potential, threshold_token)`` for each data row after
        the header.

    Note:
        Rows with fewer than four whitespace-separated fields are skipped (blank or
        malformed lines).
    """

    potentials_and_thresholds = []
    optimal_file = open(op_file, "r").readlines()
    for line in optimal_file[1:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        potentials_and_thresholds.append((fields[3], fields[2]))

    return potentials_and_thresholds


def perform_match_analysis(msa_obj, fasta_file, motif_name, output_dir, dummy_dir, pval_cutoffs, pwm_meme=None, verbose=True):
    """
    Run FIMO with a motif PWM and count hits per nucleotide for several p-value views.

    Writes MEME/MSA under ``output_dir`` when ``pwm_meme`` is missing but ``msa_obj``
    is provided. Hit counts per position use strict ``p < cutoff`` (not ``<=``).

    Args:
        msa_obj: ``nMSA``-like object used to emit MEME if needed (or ``None`` if
            ``pwm_meme`` exists).
        fasta_file (str): DNA FASTA path.
        motif_name (str): Basename for generated ``.meme.s`` / ``.msa`` files.
        output_dir (str): Directory for PWM outputs.
        dummy_dir (str): FIMO scratch directory.
        pval_cutoffs (list): Keys (floats) for per-threshold count vectors.
        pwm_meme (str, optional): Existing MEME path; if absent, derived from ``msa_obj``.
        verbose (bool): Progress to stdout.

    Returns:
        tuple: ``(matches_per_nucleotide, fimo_obj)`` — dict keyed by cutoff to
        per-position integer counts, and the ``Fimo`` wrapper.
    """

    matches_per_nucleotide = {}
    if pwm_meme is None or not os.path.isfile(str(pwm_meme)):
        if msa_obj == None:
            print("If you dont provide a meme pwm you must provide a msa_obj")
            exit(0)
        msa_obj.set_motif(motif_name)
        # Write the PWM #
        if verbose:sys.stdout.write("\t\t--PWM in MEME format...\n")
        pwm_meme = os.path.join(output_dir, motif_name + ".meme.s")
        pwm_msa = os.path.join(output_dir, motif_name + ".msa")
        if not os.path.isfile(pwm_meme):
            msa_obj.write(pwm_meme, option="meme", overwrite=True)
            msa_obj.write(pwm_msa, option="msa", overwrite=True)

    # Count the length of the fasta sequence #
    seq_length = 0
    with open(fasta_file, 'r') as ff:
        for line in ff.readlines()[1:]:
            seq_length += len(str(line.rstrip()))
    ff.close()
    # Execute FIMO # 
    if verbose:sys.stdout.write("Starting FIMO computation ..." + str(datetime.datetime.now()) + "\n")
    stored_matches = int(seq_length*2)
    sys.stdout.write(pwm_meme + " " + fasta_file + " " + str(stored_matches) + " " + dummy_dir)
    fimo_obj = fimo.get_fimo_obj(pwm_meme, fasta_file, fimo_pvalue_threshold=1.0, max_stored_matches=stored_matches, dummy_dir=dummy_dir)
    # Process FIMO results #
    working_fimo_hits = fimo_obj.get_hits(sort=True)
    # Initialize the matches_per_nucleotide dictionary if needed #
    for i in pval_cutoffs:
        matches_per_nucleotide[i] = [0]*seq_length
    for i in pval_cutoffs:   
        # We iterate over the fimo hits #
        used_fimos = []
        for fh in working_fimo_hits:
            # The fimo hits are sorted. If we go below the pval cutoff we break the loop #
            if float(fh.get_p_value()) < float(i):
                start = fh.get_start()
                end = fh.get_end()
                sequence = fh.get_sequence()
                strand = fh.get_strand()
                hit = fh.get_hit()
                pval = fh.get_p_value()
                if not str(str(start) + "_" + str(end) + "_" + strand + "_" + sequence + "_" + hit + "_" + str(pval)) in used_fimos:
                    for k in range(fh.get_start()-1, fh.get_end()):
                        matches_per_nucleotide[i][k] += 1
                    used_fimos.append(str(str(start) + "_" + str(end) + "_" + strand + "_" + sequence + "_" + hit + "_" + str(pval)))
                else:
                    continue

    
    return matches_per_nucleotide, fimo_obj



def perform_energy_analysis(msa_obj, x3dna_obj, binding_site, fasta_file, triads_obj, potentials, dummy_dir, fimo_obj, pval_cutoffs, split_potential, verbose):
    """
    Score each FIMO hit (per p-value cutoff) with statistical potentials along the DNA.

    Uses ``pwm_pbm.get_min_max_scores`` and ``get_score_for_subseq``; scales raw scores
    with ``pwm_pbm.scale``. For ``split_potential`` in ``pair``, ``s3dc``, ``local``,
    scaling uses ``-subseq_score``; otherwise the raw sign is kept.

    Args:
        msa_obj: PWM/MSA used for score calibration.
        x3dna_obj: Structure context for dinucleotide mapping.
        binding_site (dict): Dinucleotide keys mapping to triad lists for subsequence
            scoring.
        fasta_file (str): Same DNA FASTA as FIMO.
        triads_obj: Triads container (reserved for callers; unused in the body).
        potentials: Loaded statistical potentials per chain.
        dummy_dir (str): Scratch path for scoring helpers.
        fimo_obj: ``Fimo`` instance or ``None`` (empty result dicts).
        pval_cutoffs (list): Hits with ``p <= cutoff`` are scored.
        split_potential (str): Potential family key.
        verbose (bool): Timestamped progress lines.

    Returns:
        tuple: ``(energies_per_nucleotide, energies_per_nucleotide_raw)`` — per-cutoff
        lists of per-position lists of scaled and raw scores.
    """

    if verbose:sys.stdout.write("Starting energy computation ..." + str(datetime.datetime.now()) + "\n")
    energies_per_nucleotide = {}
    energies_per_nucleotide_raw = {}
    # Count the length of the fasta sequence #
    dna_seq = ""
    with open(fasta_file, 'r') as ff:
        for line in ff.readlines()[1:]:
            dna_seq += str(line.rstrip())
    ff.close()
    if fimo_obj != None:
        # Score the sequences according to different fimo p-values #
        # Initialize the dictionary #
        for pc in pval_cutoffs:
            energies_per_nucleotide[pc] = [[] for _ in range(len(dna_seq))]
            energies_per_nucleotide_raw[pc] = [[] for _ in range(len(dna_seq))]
        # Include scores in the dictionary #
        working_fimo_hits = list(set(fimo_obj.get_hits(sort=True)))
        min_score, max_score = pwm_pbm.get_min_max_scores(msa_obj, binding_site, x3dna_obj, potentials, split_potential, dummy_dir)
        for pc in pval_cutoffs:
            used_fimos = []
            for fh in working_fimo_hits:
                if float(fh.get_p_value()) <= pc:
                    start = fh.get_start()
                    end = fh.get_end()
                    sequence = fh.get_sequence()
                    strand = fh.get_strand()
                    hit = fh.get_hit()
                    pval = fh.get_p_value()
                    if not str(str(start) + "_" + str(end) + "_" + strand + "_" + sequence + "_" + hit + "_" + str(pval)) in used_fimos:
                        subseq_score = pwm_pbm.get_score_for_subseq(binding_site, x3dna_obj, sequence, potentials, split_potential)
                        # Scale the score #
                        if float(max_score - min_score) != 0.0:
                            if split_potential in ["pair", "s3dc", "local"]:
                            #if split_potential in ["s3dc", "local"]:
                                scaled_score = pwm_pbm.scale(-subseq_score, max_score, min_score)
                            else:
                                scaled_score = pwm_pbm.scale(subseq_score, max_score, min_score)
                        else:
                            scaled_score = 0.0
                        # Include this score also in the lists that correspond to higher pvalue thresholds #
                        for k in range(start, end):
                            energies_per_nucleotide[pc][k].append(scaled_score)
                            energies_per_nucleotide_raw[pc][k].append(subseq_score)

                        used_fimos.append(str(str(start) + "_" + str(end) + "_" + strand + "_" + sequence + "_" + hit + "_" + str(pval)))
                    else:
                        continue

    if verbose:sys.stdout.write("Energy computation finished ..." + str(datetime.datetime.now()) + "\n")


    return energies_per_nucleotide, energies_per_nucleotide_raw


def perform_fimo_analysis(fimo_obj, fasta_file, pval_cutoffs, verbose):
    """
    Collect per-nucleotide FIMO scores and ``-log10(p)`` for hits under each cutoff.

    For each cutoff, every position inside a hit with ``p <= cutoff`` appends that
    hit's score and log p-value (deduplicated by hit identity string).

    Args:
        fimo_obj: ``Fimo`` wrapper with sorted hits.
        fasta_file (str): DNA FASTA path (sequence after header).
        pval_cutoffs (list): Threshold keys (typically floats).
        verbose (bool): Unused; kept for API compatibility.

    Returns:
        tuple: ``(scores_per_nucleotide, logpval_per_nucleotide)`` dicts keyed by
        cutoff, values are lists-of-lists per DNA position.
    """

    scores_per_nucleotide = {}
    logpval_per_nucleotide = {} 
    # Count the length of the fasta sequence #
    dna_seq = ""
    with open(fasta_file, 'r') as ff:
        for line in ff.readlines()[1:]:
            dna_seq += str(line.rstrip())
    ff.close()
    # Score the sequences according to different fimo p-values #
    working_fimo_hits = fimo_obj.get_hits(sort=True)
    # Initialize dictionaries #
    for pc in pval_cutoffs:
        scores_per_nucleotide[pc] = [[] for _ in range(len(dna_seq))]
        logpval_per_nucleotide[pc] = [[] for _ in range(len(dna_seq))]
        used_fimos = []
        for fh in working_fimo_hits:
            # The fimo hits are sorted. If we go below the pval cutoff we break the loop #
            if float(fh.get_p_value()) <= pc:
                start = fh.get_start()
                end = fh.get_end()
                sequence = fh.get_sequence()
                strand = fh.get_strand()
                hit = fh.get_hit()
                pval = fh.get_p_value()
                score = fh.get_score()
                logpval = np.log10(fh.get_p_value())*-1
                if not str(str(start) + "_" + str(end) + "_" + strand + "_" + sequence + "_" + hit + "_" + str(pval)) in used_fimos:
                    for k in range(start, end):
                        scores_per_nucleotide[pc][k].append(score)
                        logpval_per_nucleotide[pc][k].append(logpval)

                    used_fimos.append(str(str(start) + "_" + str(end) + "_" + strand + "_" + sequence + "_" + hit + "_" + str(pval)))
                else:
                    continue


    return scores_per_nucleotide, logpval_per_nucleotide


def dummy_profile(fasta_file, pval_cutoffs):
    """
    Build a neutral per-position profile (single zero score per nucleotide per cutoff).

    Args:
        fasta_file (str): DNA FASTA path.
        pval_cutoffs (list): Keys for the output dict.

    Returns:
        dict: ``cutoff -> [[0.0], ...]`` with length matching DNA sequence.
    """

    profile = {}
    # Get DNA length #
    dna_seq = ""
    with open(fasta_file, 'r') as ff:
        for line in ff.readlines()[1:]:
            dna_seq += str(line.rstrip())
    ff.close()
    dna_len = len(dna_seq)
    # Iterate pval_cutoffs #
    for pc in pval_cutoffs:
        profile[pc] = [[0.0] for _ in range(dna_len)]

    return profile


def merge_data_per_nucleotide(match_list, mode, pval_cutoffs):
    """
    Combine per-motif per-nucleotide vectors into one structure for a cluster.

    ``add`` sums integer match counts across motifs. ``merge`` unions list-valued
    tracks so stricter cutoffs inherit the broadest hit lists from cutoff ``1.0``.

    Args:
        match_list (list): One dict per motif, each keyed by ``pval_cutoffs``.
        mode (str): ``\"add\"`` or ``\"merge\"``.
        pval_cutoffs (list): Cutoff keys shared across dicts.

    Returns:
        dict: Merged structure keyed by cutoff; shape matches input mode.
    """

    final_dict = {}
    seq_length = len(match_list[0][pval_cutoffs[0]])
    if mode == "add":
        # Merge dictionaries by adding the values in the different lists #
        for pc in pval_cutoffs:
            final_dict[pc] = [0]*seq_length
            for md in match_list:
                for i in range(0, seq_length):
                    final_dict[pc][i] += md[pc][i]
        

    elif mode == "merge":
        
        final_dict[1.0] = [[] for _ in range(seq_length)]
        for md in match_list:
            for i in range(0, seq_length):
                final_dict[1.0][i] += md[1.0][i]
        for pc in pval_cutoffs:
            if pc == 1.0:
                continue
            final_dict[pc] = [[] for _ in range(seq_length)]
            for md in match_list:
                for i in range(0, seq_length):
                    # If any motif matches a position with a certain threshold then it's the same than the one with theshold 1.0 #
                    if len(md[pc][i]) > 0:
                        final_dict[pc][i] = final_dict[1.0][i]


    return final_dict


def scale_to_percentage(matches_per_nucleotide):
    """
    Normalize each cutoff's match-count vector to percent of the maximum count.

    Args:
        matches_per_nucleotide (dict): Cutoff keys to lists of numeric counts.

    Returns:
        dict: Same keys, values scaled in-place and returned.
    """

    for pc in sorted(matches_per_nucleotide.keys()):
        # Get the highest score #
        max_score = max(matches_per_nucleotide[pc])
        for k in range(0, len(matches_per_nucleotide[pc])):
            if max_score != 0.0:
                matches_per_nucleotide[pc][k] = (float(matches_per_nucleotide[pc][k])/float(max_score))*100.0
            else:
                matches_per_nucleotide[pc][k] = 0.0


    return matches_per_nucleotide


def average_scores(energies_per_nucleotide, pval_cutoffs):
    """
    Replace each position's list of scores with ``[mean, std]`` (or ``[0,0]`` if empty).

    Args:
        energies_per_nucleotide (dict): Cutoff -> list of per-position score lists.
        pval_cutoffs (list): Cutoffs to process.

    Returns:
        dict: Updated ``energies_per_nucleotide``.
    """

    # Iterate over the different pvalue cutoffs #
    for pc in pval_cutoffs:
        for k in range(0, len(energies_per_nucleotide[pc])):
            if len(energies_per_nucleotide[pc][k]) > 0:
                mean_std = []
                mean_std.append(np.mean(energies_per_nucleotide[pc][k]))
                mean_std.append(np.std(energies_per_nucleotide[pc][k]))
                energies_per_nucleotide[pc][k] = mean_std
            else:
                energies_per_nucleotide[pc][k] = [0.0, 0.0]
    

    return energies_per_nucleotide


def set_window_size(dict_per_nucleotide, winsize, mode, pval_cutoffs, window_spacing):
    """
    Downsample per-nucleotide series with sliding windows for plotting.

    ``simple`` averages scalar values per window. ``average`` averages mean and
    std-dev tracks from ``average_scores`` output.

    Args:
        dict_per_nucleotide (dict): Cutoff -> per-position data.
        winsize (int): Window width in nucleotides (1 = no smoothing).
        mode (str): ``\"simple\"`` or ``\"average\"``.
        pval_cutoffs (list): Cutoffs present in ``dict_per_nucleotide``.
        window_spacing (int): Step between window centers.

    Returns:
        dict: For each cutoff, ``value`` / optional ``std_dev`` and
        ``nucleotide_position`` lists (1-based positions in output).
    """

    nuc_len = len(dict_per_nucleotide[pval_cutoffs[0]])
    if mode == "simple":
        # The dictionary we are working with only contains an array of values #
        output_dict = {}
        for pc in pval_cutoffs:
            output_dict[pc] = {"value": [], "nucleotide_position": []}
            if winsize > 1:
                for k in range(int(winsize/2), nuc_len-int(winsize/2), window_spacing):
                    # Iterate all positions in the window and get their values #
                    w_val = []
                    for i in range(int(k)-int(winsize/2), int(k)+int(winsize/2)):
                        w_val.append(dict_per_nucleotide[pc][i])
                    output_dict[pc]["value"].append(np.mean(w_val))
                    output_dict[pc]["nucleotide_position"].append(k+1)
            else:
                for k in range(0, nuc_len, window_spacing):
                    output_dict[pc]["value"].append(dict_per_nucleotide[pc][k])
                    output_dict[pc]["nucleotide_position"].append(k+1)


    if mode == "average":
        # The dictionary we are working with is an array of lists containing the average and the standard deviation #
        output_dict = {}
        for pc in pval_cutoffs:
            output_dict[pc] = {"value": [], "std_dev": [], "nucleotide_position": []}
            if winsize > 1:
                for k in range(int(winsize/2), nuc_len-int(winsize/2), window_spacing):
                    # Iterate all positions in the window and get their values #
                    w_av = []
                    w_std = []
                    for i in range(int(k)-int(winsize/2), int(k)+int(winsize/2), window_spacing):
                        w_av.append(dict_per_nucleotide[pc][i][0])
                        w_std.append(dict_per_nucleotide[pc][i][1])
                    output_dict[pc]["value"].append(np.mean(w_av))
                    output_dict[pc]["std_dev"].append(np.mean(w_std))
                    output_dict[pc]["nucleotide_position"].append(k+1)
            else:
                for k in range(0, nuc_len, window_spacing):
                    output_dict[pc]["value"].append(dict_per_nucleotide[pc][k][0])
                    output_dict[pc]["std_dev"].append(dict_per_nucleotide[pc][k][1])
                    output_dict[pc]["nucleotide_position"].append(k+1)

    return output_dict


def write_fimo_output(matches_per_nucleotide, form, input_label, output_dir, verbose, plot=False):
    """
    Export FIMO match counts as BED (``form == \"bed\"``) with scores on a 0–1000 scale.

    Args:
        matches_per_nucleotide (dict): Cutoff -> per-position counts.
        form (str): Output format; only ``\"bed\"`` is implemented here.
        input_label (str): Track/chromosome label prefix.
        output_dir (str): Destination folder.
        verbose (bool): Log written path.
        plot (bool): Unused.

    Returns:
        None.
    """

    if form == "bed":
        # Write a bed file #
        feature_list = []
        highest_value = 0
        for pv in list(matches_per_nucleotide.keys()):
            start = None
            value = None
            previous_value = None
            for k in range(0, len(matches_per_nucleotide[pv])):
                value = matches_per_nucleotide[pv][k]
                if start == None:
                    start = k
                    previous_value = value
                if value != previous_value:
                    feature_list.append([start, k, pv, value])
                    start = k + 1
                    previous_value = value
                    if value > highest_value:
                        highest_value = value
        # Set the value of one PWM match in a 1 to 1000 scale (bed files) #
        match_value = int(round(float(1000.0/float(highest_value)), 0))
        # Write the file #
        bed_file = os.path.join(output_dir, input_label + "_fimo.bed")
        bf = open(bed_file, "w")
        bf.write("track name=TF_binding description='Statistical potentials based PWMs matches' useScore=1\n")
        for feat in sorted(feature_list, key=lambda x: int(x[0])):
            bf.write(input_label + "\t" + str(feat[0]) + "\t" + str(feat[1]) + "\t" + str(feat[2]) + "\t" + str(feat[3]*match_value) + "\n")
        bf.close()
        if verbose:sys.stdout.write("\t\t--Output BED file written in " + bed_file + "\n")
        

def get_positive_scores(scores):
    """
    Flip the sign of every value in a scores dict (in-place).

    Args:
        scores (dict): Numeric values by key.

    Returns:
        dict: Same object after mutation.
    """

    for key in list(scores.keys()):
        scores[key] = scores[key]*-1

    return scores


def cut_dna_sequence(dna_seq, dummy_dir, start, stop):
    """
    Write a FASTA slice of ``dna_seq`` into ``dummy_dir`` for downstream tools.

    Args:
        dna_seq (str): Path to a single-sequence FASTA.
        dummy_dir (str): Output directory.
        start (int, optional): 1-based inclusive start (``None`` = from beginning).
        stop (int, optional): 1-based inclusive end (``None`` = to end).

    Returns:
        str: Path to the written dummy FASTA fragment.
    """

    # Get the new sequence #
    header, sequence = functions.parse_fasta_file_single_sequence(dna_seq)
    header = os.path.basename(dna_seq)
    new_sequence = ""
    if (start != None) and (stop != None):
        new_sequence = sequence[int(start-1):int(stop)]
    elif (start != None) and (stop == None):
        new_sequence = sequence[int(start-1):]
    elif (start == None) and (stop != None):
        new_sequence = sequence[:int(stop)]
    # Save the sequence in a file in the dummy dir #
    dummy_file = os.path.join(dummy_dir, header + str(start) + "_" + str(stop) + "_dummy.fa")
    df = open(dummy_file, "w")
    df.write(">" + header.rstrip() + "\n" + new_sequence.rstrip())
    df.close()

    return dummy_file


def include_pwm_into_json(json_dict, msa_obj, pwm_meme, root_path, dummy_dir):
    """
    Materialize MEME/PWM/logo files and attach PWM metadata to the last motif in ``json_dict``.

    Expects ``json_dict[\"motifs\"][-1]`` to exist; fills ``pwm`` with paths,
    frequency matrix, and information content.

    Args:
        json_dict (dict): Experiment/motif bundle (mutated).
        msa_obj: MSA object with ``write`` for MEME/PWM.
        pwm_meme (str): Target MEME path (``.meme.s``).
        root_path (str): Stripped from paths for web-relative URIs.
        dummy_dir (str): Passed to ``pwm_pbm.write_logo``.

    Returns:
        dict: Updated ``json_dict``.
    """

    # Set the names of the files #
    meme_file = pwm_meme
    pwm_file = pwm_meme.replace(".meme.s", ".pwm")
    logo_file = pwm_meme.replace(".meme.s", ".logo")
    # Create the meme, raw_pwm and logo files #
    if not os.path.isfile(meme_file):
        msa_obj.write(meme_file, option="meme")
    msa_obj.write(pwm_file, option="pwm")
    pwm_pbm.write_logo(msa_obj, logo_file, dummy_dir)
    # Parse the output files and include them into the json dict #
    if not "pwm" in json_dict["motifs"][-1]:
        json_dict["motifs"][-1]["pwm"] = {}
    json_dict["motifs"][-1]["pwm"]["meme_file"] = meme_file.replace(root_path, "")
    json_dict["motifs"][-1]["pwm"]["pwm_file"] = pwm_file.replace(root_path, "")
    json_dict["motifs"][-1]["pwm"]["logos"] = [logo_file.replace(root_path, "") + ".fwd.png", logo_file.replace(root_path, "") + ".rev.png"] 
    json_dict["motifs"][-1]["pwm"]["frequencies"] = []
    raw_pwm = []
    for line in functions.parse_list_file(meme_file):
        m = re.search("(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)", line)
        if m:
            raw_pwm.append([float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))])
    for i in range(0, 4):
        row = []
        for j in range(0, len(raw_pwm)):
            row.append(raw_pwm[j][i])
        json_dict["motifs"][-1]["pwm"]["frequencies"].append(row)
    # Compute the information content of the pwm #
    pwm_informativity = pwm_pbm.get_pwm_informativity(meme_file)
    json_dict["motifs"][-1]["pwm"]["information_content"] = round(pwm_informativity, 3)

    return json_dict


def write_energies_output(energies_per_nucleotide, form, input_label, output_dir, verbose, plot=False):
    """
    Write per-position mean/std energy summaries to CSV or BED and pickle the series.

    Args:
        energies_per_nucleotide (list): Per-position lists of raw scores (iterable
            of positions; each position aggregated to mean/std in the loop).
        form (str): ``\"csv\"`` or ``\"bed\"``.
        input_label (str): Output basename prefix.
        output_dir (str): Destination directory.
        verbose (bool): Unused.
        plot (bool): Unused.

    Returns:
        None. Writes ``*_energies.{form}`` and pickles ``[[mean, std], ...]`` per
        position to ``energies_per_nucleotide.p`` (does not mutate the input list).
    """

    data = energies_per_nucleotide
    output_file = os.path.join(output_dir, input_label + "_energies." + form)
    row_summaries = []
    with open(output_file, "w") as e_file:
        if form == "csv":
            e_file.write("Nucleotide\tEnergy_average\tEnergy_std_dev\n")
        if form == "bed":
            e_file.write("track name=TF_binding description='Statistical potentials based scores' useScore=1\n")
        for i in range(0, len(data)):
            av = np.mean(data[i])
            std_dev = np.std(data[i])
            row_summaries.append([av, std_dev])
            if form == "csv":
                e_file.write(str(i) + "\t" + str(av) + "\t" + str(std_dev) + "\n")
            if form == "bed":
                e_file.write(input_label + "\t" + str(i) + "\t" + str(i + 1) + "\t" + "Es3dc_dd" + "\t" + str(av * 1000) + "\n")
    pickle_path = os.path.join(output_dir, "energies_per_nucleotide.p")
    with open(pickle_path, "wb") as pf:
        pickle.dump(row_summaries, pf)

#-----------#
#   Plots   #
#-----------#

def make_plot_per_nucleotide(matches_per_nucleotide, input_label, output_dir, pval_cutoffs, mode):
    """
    Save matplotlib PNGs for several window sizes/spacings over FIMO match profiles.

    Args:
        matches_per_nucleotide (dict): Cutoff -> per-position values (``simple`` mode)
            or ``average_scores`` output (``average`` mode).
        input_label (str): Basename prefix for PNG files.
        output_dir (str): Parent of ``plots/`` (``plots`` is sibling of this dirname).
        pval_cutoffs (list): Thresholds to draw.
        mode (str): ``\"simple\"`` or ``\"average\"`` (passed to ``set_window_size``).

    Returns:
        None.
    """

    plot_dir = os.path.join(os.path.dirname(output_dir), "plots")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    # Include the results for p-value = 1.0 #
    winsizes = [1, 10, 50, 100, 500, 1000]
    win_spacings = [1, 5, 25, 50, 250, 500]
    # Iterate over different window sizes and window spacings # 
    for ws in winsizes:
        for wsp in win_spacings:
            fig = figure(figsize=(22, 8))
            # Get the dictionary with the appropiate window size and spacing #
            dict_per_nucleotide = set_window_size(matches_per_nucleotide, ws, mode, pval_cutoffs, wsp)
            # Now iterate over the different p-value cutoffs #
            handle_list = []
            for pc in sorted(pval_cutoffs, reverse=True):
                # Se colors for the different thresholds #
                if pc == 1.0:
                    color = "black"
                if pc == 0.5:
                    color = "blue"
                if pc == 0.05:
                    color = "green"
                if pc == 0.005:
                    color = "orange"
                if pc == 0.0005:
                    color = "red"
                # Plot the values #
                l, = plt.plot(dict_per_nucleotide[pc]["nucleotide_position"], dict_per_nucleotide[pc]["value"], color=color, linewidth=2.0, linestyle="-")
                handle_list.append(l)
                if mode == "average":
                    # Include standard deviation #
                    top_std = []
                    bottom_std = []
                    for i in range(0, len(dict_per_nucleotide[pc]["value"])):
                        top_std.append(dict_per_nucleotide[pc]["value"][i] + dict_per_nucleotide[pc]["std_dev"][i])
                        bottom_std.append(dict_per_nucleotide[pc]["value"][i] - dict_per_nucleotide[pc]["std_dev"][i])
                    plt.fill_between(dict_per_nucleotide[pc]["nucleotide_position"], top_std, bottom_std, facecolor=color, interpolate=True, alpha=0.2)
            # Set the legend and save the plot #
            legend(handles=(handle_list), labels=("Thr=1.0", "Thr=0.5", "Thr=0.05", "Thr=0.005", "Thr=0.0005"), bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., prop={'size':14})
            plt.subplots_adjust(left=0.05, right=0.9, top=0.9, bottom=0.1)
            fig.savefig(os.path.join(plot_dir, input_label + "_winsize" + str(ws) + "_winspace" + str(wsp) + "_fimo_matches.png"))
            plt.close(fig)


# Include a variable wwith DNA sequence #
def plotly_plot(html_file, binding_profile_dict_list, binding_profile_variables, dna_seq, dna, motif, fullscreen=False):
    """
    Build an interactive Plotly HTML for one DNA: matches, log p-values, and energies.

    ``binding_profile_dict_list`` and ``binding_profile_variables`` are parallel:
    each dict is keyed by p-value cutoff; variables include ``\"matches\"``,
    ``\"logpval\"``, ``\"energies_s3dc_dd\"``.

    Args:
        html_file (str): Output HTML path (non-fullscreen).
        binding_profile_dict_list (list): Profile dicts to layer.
        binding_profile_variables (list): Variable name per dict (parallel).
        dna_seq (str): DNA sequence string (one letter per position).
        dna (str): Label for the x-axis title.
        motif (str): Plot title fragment.
        fullscreen (bool): If True, also write ``*_fullscreen.html``.

    Returns:
        None.
    """

    import plotly
    import plotly.graph_objs as go
    from plotly.offline import plot

    final_trace_list = []
    trace_list = []

    print(("html_file: " + str(html_file) + "    " + str(type(html_file))))
    print(("binding_profile_variables: " + str(binding_profile_variables) + "    " + str(type(binding_profile_variables))))
    print(("dna_seq: " + str(dna_seq) + "    " + str(type(dna_seq))))
    print(("dna: " + str(dna) + "    " + str(type(dna_seq))))
    print(("motif: " + str(motif) + "    " + str(type(dna_seq))))

    for k in range(0, len(binding_profile_dict_list)):
        # Define the binding profile dictionary and the variable to work with #
        binding_profile_dict = binding_profile_dict_list[k]
        var = binding_profile_variables[k]
        # Generate traces and append them to a trace list #
        pvalue_cutoffs = sorted(list(binding_profile_dict.keys()), key=float, reverse=True)
        # Iterate over the different p-value cutoffs #
        if var == "matches":
            x_nucleotides = []
            x_values = list(range(1, len(binding_profile_dict[1.0]) + 1))
            for x in x_values:
                #x_nucleotides.append(dna_seq[x-1] + "\n" + str(x))
                x_nucleotides.append(dna_seq[x-1] + " " + str(x))
            color_list = ["firebrick", "orangered", "orange", "gold"]
            for i in range(0, len(pvalue_cutoffs)):
                y_vals = []
                # Remove excesive decimals #
                for val in binding_profile_dict[pvalue_cutoffs[i]]:
                    y_vals.append(round(float(val), 2))
                if float(pvalue_cutoffs[i]) != float(10**-1):
                    trace = go.Scatter(
                        x = x_nucleotides,
                        y = y_vals,
                        mode = 'lines',
                        name = str(pvalue_cutoffs[i]),
                        line = dict(
                            color = (color_list[i]),
                            width = 3,),
                        visible = "legendonly"
                    )
                    trace_list.append(trace)
                elif float(pvalue_cutoffs[i]) == float(10**-1):
                    trace = go.Scatter(
                        x = x_nucleotides,
                        y = y_vals,
                        mode = 'lines',
                        name = str(pvalue_cutoffs[i]),
                        line = dict(
                            color = (color_list[i]),
                            width = 3),
                        visible = True
                    )
                    trace_list.append(trace)
        elif var in ["logpval", "energies_s3dc_dd"]:
            pvalue_cutoffs = sorted(list(binding_profile_dict.keys()), key=float, reverse=True)
            color_list = ["firebrick", "orangered", "orange", "gold"]
            #for i in range(0, len(pvalue_cutoffs)):
            av_values = []
            std_values = []
            for tup in binding_profile_dict[pvalue_cutoffs[0]]:
                av_values.append(round(float(tup[0]), 2))
                std_values.append(round(float(tup[1]), 2))
            upper_values = []
            lower_values = []
            filtered_av_values = []
            x_values = []
            
            for k in range(0, len(av_values)):
                x_values.append(k+1)
                upper_values.append((av_values[k]+std_values[k]))
                lower_values.append((av_values[k]-std_values[k]))
            x_nucleotides = []
            for x in x_values:
                x_nucleotides.append(dna_seq[x-1] + "\n" + str(x))
            trace = go.Scatter(
                x = x_nucleotides,
                y = av_values,
                mode = 'lines',
                name = str(pvalue_cutoffs[0]),
                legendgroup = str(pvalue_cutoffs[0]),
                line = dict(
                    color = (color_list[0]),
                    width = 3,),
                fill='tonexty',
                visible = False
            )
            upper_bound = go.Scatter(
                x = x_nucleotides,
                y = upper_values,
                mode = 'lines',
                name = str(pvalue_cutoffs[0]),
                legendgroup = str(pvalue_cutoffs[0]),
                line = dict(
                    color=(color_list[0]),
                    width = 0.8),
                fill='tonexty',
                showlegend = False,
                hoverinfo='none',
                visible = False
            )
            lower_bound = go.Scatter(
                x = x_nucleotides,
                y = lower_values,
                mode = 'lines',
                name = str(pvalue_cutoffs[0]),
                legendgroup = str(pvalue_cutoffs[0]),
                line = dict(
                    color=(color_list[0]),
                    width = 0.8),
                showlegend = False,
                hoverinfo='none',
                visible = False
            )
            trace_list.append(lower_bound)
            trace_list.append(trace)
            trace_list.append(upper_bound)
    # Define the layout and axis names #
    config = {'linkText': "", 'scrollZoom': True}
    # Determine a dropdown with the variables that will be shown in this plot #
    
    updatemenus = [
        dict(
            buttons = [
                dict(
                    label='Matches',
                    method='update',
                    args=[{'visible': ["legendonly", True, "legendonly", "legendonly", False, False, False, False, False, False]},
                        {'yaxis': {'title': 'Matches', 'titlefont': dict(size=16), 'range': [-10.0, 110.0]}}]
                    ),
                dict(
                    label='FIMO: Log10(P-value)',
                    method='update',
                    args=[{'visible': [False, False, False, False, True, True, True, False, False, False]},
                        {'yaxis': {'title': 'FIMO: Log10(P-value)', 'titlefont': dict(size=16)}}]
                    ),
                dict(
                    label='Statistical energy: s3dc_dd',
                    method='update',
                    args=[{'visible': [False, False, False, False, False, False, False, True, True, True]},
                        {'yaxis': {'title': 'Statistical energy: s3dc_dd', 'titlefont': dict(size=16), 'range': [0.0, 1.0]}}]
                    ),
                ],
            direction = 'down',
            pad = {'r': 10, 't': 40},
            showactive = True,
            x = 0.19,
            y = 1.40,
            bgcolor = "white"
        )
    ]
    
    annotations = list([
        dict(text='P-value\nthresholds:', x=1.085, y=1.20, xref='paper', yref='paper', showarrow=False, font=dict(size=14), align="center")
    ])
    layout = dict(
    
        title = 'Binding profile for motif: ' + str(motif),
        xaxis = dict(title = 'Nucleotide positions for DNA: ' + str(dna), rangeslider=dict(visible=True, bordercolor="grey", borderwidth=4), range=[1, len(binding_profile_dict_list[0][1.0])], titlefont=dict(size=16), tickformat=',d'),
        yaxis = dict(title = "Matches", titlefont=dict(size=16)),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        updatemenus = updatemenus,
        annotations = annotations
        )

    # Create the plot #
    fig = dict(data=trace_list, layout=layout)#, config=config)
    #plot(fig, filename=html_file, show_link=False, auto_open=False)
    plot(fig, filename=html_file, show_link=False)

    
    if fullscreen == True:

        updatemenus = [
            dict(
                buttons = [
                    dict(
                        label='Matches',
                        method='update',
                        args=[{'visible': ["legendonly", True, "legendonly", "legendonly", False, False, False, False, False, False]},
                            {'yaxis': {'title': 'Matches', 'titlefont': dict(size=20), 'range': [-10.0, 110.0]}}]
                        ),
                    dict(
                        label='FIMO: Log10(P-value)',
                        method='update',
                        args=[{'visible': [False, False, False, False, True, True, True, False, False, False]},
                            {'yaxis': {'title': 'FIMO: Log10(P-value)', 'titlefont': dict(size=20)}}]
                        ),
                    dict(
                        label='Statistical energy: s3dc_dd',
                        method='update',
                        args=[{'visible': [False, False, False, False, False, False, False, True, True, True]},
                            {'yaxis': {'title': 'Statistical energy: s3dc_dd', 'titlefont': dict(size=20), 'range': [0.0, 1.0]}}]
                        ),
                    ],
                direction = 'down',
                pad = {'r': 20, 't': 45},
                font = {'size': 16},
                showactive = True,
                x = 0.2,
                y = 1.29,
                bgcolor = "white"
            )
        ]

        annotations = list([
            dict(text='P-value\nthresholds:', x=1.085, y=1.05, xref='paper', yref='paper', showarrow=False, font=dict(size=18), align="center")
        ])
        layout = dict(
            title = 'Binding profile for motif: ' + str(motif),
            titlefont=dict(size=24),
            xaxis = dict(title = 'Nucleotide positions for DNA: ' + str(dna), rangeslider=dict(visible=True, bordercolor="grey", borderwidth=4), range=[1, len(binding_profile_dict_list[0][1.0])], titlefont=dict(size=20), tickformat=',d'),
            yaxis = dict(title = "Matches", titlefont=dict(size=20)),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            updatemenus = updatemenus,
            annotations = annotations
            )

        # Create the plot #
        fig = dict(data=trace_list, layout=layout)#, config=config)
        #plot(fig, filename=html_file, show_link=False, auto_open=False)
        html_file_fullscreen = html_file.replace(".html", "") + "_fullscreen.html"
        plot(fig, filename=html_file_fullscreen, show_link=False)
    

def plotly_plot_comparison(html_file, binding_profile_dict_list, binding_profile_variables, dna, motif, fullscreen=False):
    """
    Plotly overlay for two DNA sequences (allele or condition comparison).

    ``binding_profile_dict_list`` has shape ``[profile_bundle_dna0, profile_bundle_dna1]``
    where each bundle is a list of dicts paired with ``binding_profile_variables``.

    Args:
        html_file (str): Output HTML path.
        binding_profile_dict_list (list): Length-2 list of profile bundles.
        binding_profile_variables (list): Shared variable names (``matches``, ...).
        dna (list): Two DNA id strings for legend and axis titles.
        motif (str): Title fragment.
        fullscreen (bool): Emit fullscreen HTML variant when True.

    Returns:
        None.
    """

    import plotly
    import plotly.graph_objs as go
    from plotly.offline import plot

    color_list_complete = [
        # Red scale #
        ["firebrick", "orangered", "orange", "gold"],
        # Blue scale #
        ["navy", "blue", "dodgerblue", "cyan"]
    ]

    trace_list = []

    for j in range(0, len(dna)):
        dna_id = dna[j]
        binding_profile_data = binding_profile_dict_list[j]
        color_list = color_list_complete[j]

        for k in range(0, len(binding_profile_data)):
            # Define the binding profile dictionary and the variable to work with #
            binding_profile_dict = binding_profile_data[k]
            var = binding_profile_variables[k]
            # Generate traces and append them to a trace list #
            pvalue_cutoffs = sorted(list(binding_profile_dict.keys()), key=float, reverse=True)
            # Iterate over the different p-value cutoffs #
            if var == "matches":
                x_values = list(range(1, len(binding_profile_dict[1.0]) + 1))
                for i in range(0, len(pvalue_cutoffs)):
                    y_vals = []
                    # Remove excesive decimals #
                    for val in binding_profile_dict[pvalue_cutoffs[i]]:
                        y_vals.append(round(float(val), 2))
                    if float(pvalue_cutoffs[i]) != 0.1:
                        trace = go.Scatter(
                            x = x_values,
                            y = y_vals,
                            mode = 'lines',
                            name = dna_id + " P < " + str(pvalue_cutoffs[i]),
                            legendgroup = " P < " + str(pvalue_cutoffs[i]),
                            line = dict(
                                color = (color_list[i]),
                                width = 3,),
                            visible = "legendonly"
                        )
                        trace_list.append(trace)
                    elif float(pvalue_cutoffs[i]) == 0.1:
                        trace = go.Scatter(
                            x = x_values,
                            y = y_vals,
                            mode = 'lines',
                            name = dna_id + " P < " + str(pvalue_cutoffs[i]),
                            legendgroup = " P < " + str(pvalue_cutoffs[i]),
                            line = dict(
                                color = (color_list[i]),
                                width = 3),
                            visible = True
                        )
                        trace_list.append(trace)
            elif var in ['logpval', 'energies_s3dc_dd']:
                pvalue_cutoffs = sorted(list(binding_profile_dict.keys()), key=float, reverse=True)
                #for i in range(0, len(pvalue_cutoffs)):
                av_values = []
                std_values = []
                for tup in binding_profile_dict[pvalue_cutoffs[0]]:
                    av_values.append(round(float(tup[0]), 2))
                    std_values.append(round(float(tup[1]), 2))
                upper_values = []
                lower_values = []
                filtered_av_values = []
                x_values = []
                for k in range(0, len(av_values)):
                    x_values.append(k+1)
                    upper_values.append((av_values[k]+std_values[k]))
                    lower_values.append((av_values[k]-std_values[k]))
                trace = go.Scatter(
                    x = x_values,
                    y = av_values,
                    mode = 'lines',
                    name = dna_id + " P < " + str(pvalue_cutoffs[0]),
                    legendgroup = " P < " + str(pvalue_cutoffs[0]),
                    line = dict(
                        color = (color_list[0]),
                        width = 3,),
                    fill='tonexty',
                    visible = False
                )
                upper_bound = go.Scatter(
                    x = x_values,
                    y = upper_values,
                    mode = 'lines',
                    name = dna_id + " P < " + str(pvalue_cutoffs[0]),
                    legendgroup = " P < " + str(pvalue_cutoffs[0]),
                    line = dict(
                        color=(color_list[0]),
                        width = 0.8),
                    fill='tonexty',
                    showlegend = False,
                    hoverinfo='none',
                    visible = False
                )
                lower_bound = go.Scatter(
                    x = x_values,
                    y = lower_values,
                    mode = 'lines',
                    name = dna_id + " P < " + str(pvalue_cutoffs[0]),
                    legendgroup = " P < " + str(pvalue_cutoffs[0]),
                    line = dict(
                        color=(color_list[0]),
                        width = 0.8),
                    showlegend = False,
                    hoverinfo='none',
                    visible = False
                )
                trace_list.append(lower_bound)
                trace_list.append(trace)
                trace_list.append(upper_bound)
    # Define the layout and axis names #
    config = {'linkText': "", 'scrollZoom': True}
    # Determine a dropdown with the variables that will be shown in this plot #
    
    updatemenus = [
        dict(
            buttons = [
                dict(
                    label='Matches',
                    method='update',
                    args=[{'visible': ["legendonly", True, "legendonly", "legendonly", False, False, False, False, False, False, 
                    "legendonly", True, "legendonly", "legendonly", False, False, False, False, False, False]},
                        {'yaxis': {'title': 'Matches', 'titlefont': dict(size=16)}}]
                    ),
                dict(
                    label='FIMO: Log10(P-value)',
                    method='update',
                    args=[{'visible': [False, False, False, False, True, True, True, False, False, False, False, False, False, 
                    False, True, True, True, False, False, False]},
                        {'yaxis': {'title': 'FIMO: Log10(P-value)', 'titlefont': dict(size=16)}}]
                    ),
                dict(
                    label='Statistical energy: s3dc_dd',
                    method='update',
                    args=[{'visible': [False, False, False, False, False, False, False, True, True, True, False, False, False, 
                    False, False, False, False, True, True, True]},
                        {'yaxis': {'title': 'Statistical energy: s3dc_dd', 'titlefont': dict(size=16), 'range': [0.0, 1.0]}}]
                    ),
                ],
            direction = 'down',
            pad = {'r': 10, 't': 24},
            showactive = True,
            x = 0.20,
            y = 1.33,
            bgcolor = "white"
        )
    ]

    annotations = list([
        #dict(text='Variable:', x=-0.03, y=1.18, xref='paper', yref='paper', showarrow=False, font=dict(size=14)),
        dict(text='P-value\nthresholds:', x=1.15, y=1.20, xref='paper', yref='paper', showarrow=False, font=dict(size=14), align="center")
    ])
    layout = dict(
        title = 'Binding profile for motif: ' + str(motif),
        xaxis = dict(title = 'Nucleotide positions for DNA: ' + str(dna[0]) + " and " + str(dna[1]), rangeslider=dict(visible=True, bordercolor="grey", borderwidth=4), range=[1, len(binding_profile_dict_list[0][0][1.0])], titlefont=dict(size=16), tickformat=',d'),
        yaxis = dict(title = "Matches", titlefont=dict(size=16)),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        updatemenus = updatemenus,
        annotations = annotations
        )

    # Create the plot #
    fig = dict(data=trace_list, layout=layout)#, config=config)
    #plot(fig, filename=html_file, show_link=False, auto_open=False)
    plot(fig, filename=html_file, show_link=False)

    if fullscreen == True:

        updatemenus = [
            dict(
                buttons = [
                dict(
                    label='Matches',
                    method='update',
                    args=[{'visible': ["legendonly", True, "legendonly", "legendonly", False, False, False, False, False, False, 
                    "legendonly", True, "legendonly", "legendonly", False, False, False, False, False, False]},
                        {'yaxis': {'title': 'Matches', 'titlefont': dict(size=16)}}]
                    ),
                dict(
                    label='FIMO: Log10(P-value)',
                    method='update',
                    args=[{'visible': [False, False, False, False, True, True, True, False, False, False, False, False, False, 
                    False, True, True, True, False, False, False]},
                        {'yaxis': {'title': 'FIMO: Log10(P-value)', 'titlefont': dict(size=16)}}]
                    ),
                dict(
                    label='Statistical energy: s3dc_dd',
                    method='update',
                    args=[{'visible': [False, False, False, False, False, False, False, True, True, True, False, False, False, 
                    False, False, False, False, True, True, True]},
                        {'yaxis': {'title': 'Statistical energy: s3dc_dd', 'titlefont': dict(size=16), 'range': [0.0, 1.0]}}]
                    ),
                ],
                direction = 'down',
                pad = {'r': 20, 't': 20},
                font = {'size': 16},
                showactive = True,
                x = 0.18,
                y = 1.31,
                bgcolor = "white"
            )
        ]

        annotations = list([
            #dict(text='Variable:', x=-0.04, y=1.24, xref='paper', yref='paper', showarrow=False, font=dict(size=18)),
            dict(text='P-value\nthresholds:', x=1.12, y=1.08, xref='paper', yref='paper', showarrow=False, font=dict(size=18), align="center")
        ])
        layout = dict(
            title = 'Binding profile for motif: ' + str(motif),
            titlefont=dict(size=24),
            xaxis = dict(title = 'Nucleotide positions for DNA: ' + str(dna[0]) + " and " + str(dna[1]), rangeslider=dict(visible=True, bordercolor="grey", borderwidth=4), range=[1, len(binding_profile_dict_list[0][0][1.0])], titlefont=dict(size=20), tickformat=',d'),
            yaxis = dict(title = "Matches", titlefont=dict(size=20)),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            updatemenus = updatemenus,
            annotations = annotations
            )

        # Create the plot #
        fig = dict(data=trace_list, layout=layout)#, config=config)
        #plot(fig, filename=html_file, show_link=False, auto_open=False)
        html_file_fullscreen = html_file.replace(".html", "") + "_fullscreen.html"
        plot(fig, filename=html_file_fullscreen, show_link=False)



def get_energies_for_individual_nucleotides(triads_obj, fimo_obj, fasta_file, potentials, split_potential, pval_cutoffs, dummy_dir, verbose, binding_site, x3dna_obj):
    """
    Per-nucleotide energy lists inside FIMO hits (partial implementation).

    Initializes per-position lists for each cutoff. The detailed per-triad
    decomposition block is commented; the active path scores whole-hit sequences
    via ``pwm_pbm.get_score_for_subseq`` when ``pval < pv`` (strict).

    Args:
        triads_obj: Triads object (reserved for future per-residue breakdown).
        fimo_obj: FIMO hits provider.
        fasta_file (str): DNA FASTA path.
        potentials: Statistical potentials by chain.
        split_potential (str): Potential family key.
        pval_cutoffs (list): Cutoff keys.
        dummy_dir (str): Unused; API compatibility.
        verbose (bool): Unused.
        binding_site (dict): Dinucleotide -> triad lists for ``get_score_for_subseq``.
        x3dna_obj: X3DNA helper for dinucleotide mapping.

    Returns:
        tuple: ``(energies_per_nucleotide, contribution_per_nucleotide)`` keyed by
        cutoff (contribution remains empty unless the commented block is enabled).
    """

    # Initialize #
    energies_per_nucleotide = {}
    contribution_per_nucleotide = {}
    # Count the length of the fasta sequence #
    dna_seq = ""
    with open(fasta_file, 'r') as ff:
        for line in ff.readlines()[1:]:
            dna_seq += str(line.rstrip())
    ff.close()
    for pv in pval_cutoffs:
        energies_per_nucleotide[pv] = [[] for _ in range(len(dna_seq))]
        contribution_per_nucleotide[pv] = [[] for _ in range(len(dna_seq))]
    # Work using the binding sites identified by fimo #
    if fimo_obj != None:
        for pv in pval_cutoffs:
            working_fimo_hits = list(set(fimo_obj.get_hits(sort=True)))
            for fh in working_fimo_hits:
                start = fh.get_start()
                end = fh.get_end()
                sequence = fh.get_sequence()
                strand = fh.get_strand()
                hit = fh.get_hit()
                pval = fh.get_p_value()
                if pval < pv:
                    # Get the score for the entire binding site and for the individual nucleotides #
                    subseq_score = pwm_pbm.get_score_for_subseq(binding_site, x3dna_obj, sequence, potentials, split_potential)
                    """
                    for i in sorted(individual_scores.keys()):
                        if subseq_score != 0:
                            ind_score = individual_scores[i]
                            contribution = (ind_score/subseq_score)*100.0
                        else:
                            ind_score = 0.0
                            contribution = 0.0
                        nuc = i+start-1
                        energies_per_nucleotide[pv][nuc].append(ind_score)
                        contribution_per_nucleotide[pv][nuc].append(contribution)
                    """

    return energies_per_nucleotide, contribution_per_nucleotide


def parse_options():
    """
    Parse CLI options for ``binding_profiles.py``.

    How to run:
        python binding_profiles.py --dummy DUMMY_DIR -i MODEL_DIR -d DNA.fa \\
            --pdb PDB_DATA --pbm PBM_DATA [-o OUT] [-l LABEL] [-v]

    Configures:
        ``dummy_dir``, ``input_dir`` (``-i``), ``dna_seq`` (``-d`` FASTA),
        ``output_dir`` (``-o``), ``pdb_dir`` (``--pdb``), ``pbm_dir`` (``--pbm``),
        ``label`` (``-l``, suffix for output names), ``verbose`` (``-v``).

    Args:
        None.

    Returns:
        optparse.Values: Namespace with the attributes above. Callers must supply
        ``input_dir``, ``dna_seq``, ``pdb_dir``, and ``pbm_dir`` or
        ``parser.error`` is invoked.

    Raises:
        SystemExit: On missing required options (via ``optparse``).
    """

    parser = optparse.OptionParser(
        "Usage: binding_profiles.py [--dummy=DUMMY_DIR] -i MODEL_DIR -d DNA.fa "
        "--pdb PDB_DIR --pbm PBM_DIR [-o OUTPUT_DIR] [-l LABEL] [-v]"
    )

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Input directory containing protein models", metavar="{filename}")
    parser.add_option("-d", action="store", type="string", dest="dna_seq", help="DNA sequence in a fasta file.", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("-l", action="store", type="string", default="", dest="label", help="Label to include in the output models name", metavar="{str}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")

    (options, args) = parser.parse_args()

    if options.input_dir is None or options.pdb_dir is None or options.pbm_dir is None or options.dna_seq is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


def main():
    """
    CLI entry: score each model against DNA, cluster motifs, and write JSON profiles.

    Workflow:
        1. Parse options and create output folders (``pwms/``, ``individual_profiles/``,
           ``cluster_profiles/``).
        2. For each split-potential in ``pair`` and ``s3dc_dd`` and for each PMF
            flag (raw vs z-scored branch from ``pwm_pbm``), load statistical
            potentials per template, merge chain binding sites for scoring, build
            the MSA, run FIMO, statistical energy analysis, and per-hit energy lists.
        3. Cluster motifs with ``cluster_motifs.get_motif_clusters``, merge or
           average tracks across cluster members, write cluster JSON files, and
           append one line per cluster to ``clusters.txt``.

    Note:
        Structural homolog reuse (``structural_homologs_by_chain``) resets when the
        template id changes between successive models in the input directory.

    Returns:
        None.
    """

    # Arguments & Options #
    options = parse_options()
    # Get the uniprot IDs for the proteins under study #
    input_dir = os.path.abspath(options.input_dir)
    dna_seq = os.path.abspath(options.dna_seq)
    output_dir = os.path.abspath(options.output_dir)
    dummy_dir = os.path.abspath(options.dummy_dir)
    pdb_dir = os.path.abspath(options.pdb_dir)
    pbm_dir = os.path.abspath(options.pbm_dir)
    verbose = options.verbose
    if options.label == None:
        input_label = os.path.basename(dna_seq).replace(".fa", "")
    else:
        input_label = os.path.basename(dna_seq).replace(".fa", "") + "_" + options.label

    _label_for_files = (options.label or "").strip()

    # Create output directories #
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(os.path.join(output_dir, "pwms")):
        os.mkdir(os.path.join(output_dir, "pwms"))
    if not os.path.exists(os.path.join(output_dir, "individual_profiles")):
        os.mkdir(os.path.join(output_dir, "individual_profiles"))
    if not os.path.exists(os.path.join(output_dir, "cluster_profiles")):
        os.mkdir(os.path.join(output_dir, "cluster_profiles"))


    # We define the pvalue thresholds that we will use to define the differents levels of affinity for the TF-DNA interaction #
    pval_cutoffs = sorted([1.0, float(10**-1), float(10**-2), float(10**-3)], key=float, reverse=False)
    # Split-potential modes to evaluate (do not reuse name ``potentials``; loader returns a dict).
    split_potential_names = ["pair", "s3dc_dd"]
    pmfs = [True, False]
    # Load the families dictionary #
    families = {}
    for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family

    # Split-potential name ``pot`` selects the statistical potential and MSA mode for FIMO/PWM naming #
    structural_homologs_by_chain = None
    for pot in split_potential_names:
        for pmf in pmfs:
            # Handle the labeling of the output files #
            if pmf == True:
                pot_label = pot + "_pmf"
            else:
                pot_label = pot
            pot_tag = (pot_label + "_" + _label_for_files) if _label_for_files else pot_label

            # Create a PWM for each model #
            motif_list = []
            for model in sorted(os.listdir(input_dir)):
                if not model.endswith(".pdb"):
                    continue

                motif = {}
                motif["pdb_file"] = os.path.join(input_dir, model)
                motif["pdb_name"] = model.split(".pdb")[0]
                # Handle the domains #
                if motif["pdb_name"].count(":") == 2:
                    motif["domain"] = motif["pdb_name"].split(":")[1] + ":" + motif["pdb_name"].split(":")[2].split("_")[0]
                elif motif["pdb_name"].count(":") == 4:
                    motif["domain"] = motif["pdb_name"].split(":")[1] + ":" + motif["pdb_name"].split(":")[2] + "," + motif["pdb_name"].split(":")[3] + ":" + motif["pdb_name"].split(":")[4].split("_")[0]
                motif["template"] = motif["pdb_name"].split("_")[1] + "_" + motif["pdb_name"].split("_")[2]
                # If the model has the same template than the last motif we keep the structural homologs by chain #
                if len(motif_list) > 0:
                    if motif_list[-1]["template"] != motif["template"]:
                        structural_homologs_by_chain = None
                        
                # Get PDB object #
                pdb_path = os.path.join(input_dir, model)
                pdb_obj = PDB(pdb_path)
                motif_name = model.replace(".pdb", "")
            
                # Get DSSP object #
                if verbose:sys.stdout.write("\t\t-- calculate secondary structure ...\n")
                dssp_obj = dssp.get_dssp_obj(os.path.join(input_dir, model), dummy_dir)

                # Get X3DNA object #
                if verbose:sys.stdout.write("Get DNA %s ...\n"%os.path.join(input_dir, model))
                x3dna_obj = x3dna.get_x3dna_obj(os.path.join(input_dir, model))
                if len(list(x3dna_obj.get_dinucleotides().keys())) < 1:
                    sys.stdout.write("Missing DNA ...\n")
                    continue
                
                # Get contacts object #
                if verbose:sys.stdout.write("\t\t-- calculate contacts ...\n")
                contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj, contacts_type="pdi", distance_type="dinucleotides", dummy_dir=dummy_dir)
                if len(contacts_obj.get_contacts())<1:
                    sys.stdout.write("Missing Protein-DNA contacts ...\n")
                    continue

                # Get triads object #
                if verbose:sys.stdout.write("\t\t-- calculate protein-dna pairs ...\n")
                triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)

                # This means that these two models don't have the same template, new potentials must be computed #
                if verbose:sys.stdout.write("Load potentials ..." + str(datetime.datetime.now()) + "\n")
                auto_mode_eff = True if ((pot == "s3dc_dd") and (not pmf)) else False
                potentials, thresholds, radii, structural_homologs_by_chain = pwm_pbm.load_statistical_potentials(
                    pdb_obj,
                    os.path.abspath(options.pdb_dir),
                    pbm_dir,
                    families,
                    0.0,
                    potential_file_entry=None,
                    split_potential=pot,
                    auto_mode=auto_mode_eff,
                    family_potentials=False,
                    pbm_potentials=False,
                    score_threshold=False,
                    taylor_approach=False,
                    pmf_approach=pmf,
                    known_pdb=False,
                    structural_homologs_by_chain=structural_homologs_by_chain,
                    dummy_dir=dummy_dir,
                    verbose=verbose,
                )

                merged_binding_site = {}
                for pdb_chain in sorted(potentials.keys()):
                    _, bs = pwm_pbm.get_scores_and_binding_site(
                        triads_obj, x3dna_obj, potentials, radii, None, None, pot, pdb_chain
                    )
                    for dinuc, triad_list in bs.items():
                        merged_binding_site.setdefault(dinuc, []).extend(triad_list)

                if verbose:sys.stdout.write("Build MSA ...\n")
                msa_obj = pwm_pbm.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, None, None, pot, thresholds)

                # Get the data per nucleotide #
                motif_label = motif_name + "_" + pot_label
                matches_per_nucleotide, fimo_obj = perform_match_analysis(msa_obj, dna_seq, motif_label, output_dir, dummy_dir, pval_cutoffs, verbose=verbose)
                scores_per_nucleotide, logpval_per_nucleotide = perform_fimo_analysis(fimo_obj, dna_seq, pval_cutoffs, verbose)
                energies_per_nucleotide, energies_per_nucleotide_raw = perform_energy_analysis(msa_obj, x3dna_obj, merged_binding_site, dna_seq, triads_obj, potentials, dummy_dir, fimo_obj, pval_cutoffs, pot, verbose)
                energies_per_individual_nucleotide, contribution_per_individual_nucleotide = get_energies_for_individual_nucleotides(triads_obj, fimo_obj, dna_seq, potentials, pot, pval_cutoffs, dummy_dir, verbose, merged_binding_site, x3dna_obj)
                
                # Append individual dictionaries of nucleotide data to a list, they will be merged afterwards #
                motif["complete_meme_pwm"] = os.path.join(output_dir, "pwms", motif_label + ".meme.s")
                motif["msa_file"] = os.path.join(output_dir, "pwms", motif_label + ".msa")
                motif["matches_per_nucleotide"] = matches_per_nucleotide
                motif["scores_per_nucleotide"] = scores_per_nucleotide
                motif["logpval_per_nucleotide"] = logpval_per_nucleotide
                motif["energies_per_nucleotide"] = energies_per_nucleotide
                motif["energies_per_nucleotide_raw"] = energies_per_nucleotide_raw
                motif["energies_per_individual_nucleotide"] = energies_per_individual_nucleotide
                motif["contribution_per_individual_nucleotide"] = contribution_per_individual_nucleotide
                motif["IC"] = pwm_pbm.get_pwm_informativity(motif["complete_meme_pwm"])
                
                # Save the individual profiles as json files #
                functions.dumpJSON(os.path.join(output_dir, "individual_profiles", "FIMO_matches_" + motif["pdb_name"] + "_" + pot_tag + ".json"), matches_per_nucleotide)
                functions.dumpJSON(os.path.join(output_dir, "individual_profiles", "FIMO_scores_" + motif["pdb_name"] + "_" + pot_tag + ".json"), scores_per_nucleotide)
                functions.dumpJSON(os.path.join(output_dir, "individual_profiles", "FIMO_logpval_" + motif["pdb_name"] + "_" + pot_tag + ".json"), logpval_per_nucleotide)
                functions.dumpJSON(os.path.join(output_dir, "individual_profiles", "energies_" + motif["pdb_name"] + "_" + pot_tag + ".json"), energies_per_nucleotide)
                functions.dumpJSON(os.path.join(output_dir, "individual_profiles", "energies_raw_" + motif["pdb_name"] + "_" + pot_tag + ".json"), energies_per_nucleotide_raw)
                functions.dumpJSON(os.path.join(output_dir, "individual_profiles", "individual_energies_" + motif["pdb_name"] + "_" + pot_tag + ".json"), energies_per_individual_nucleotide)
                functions.dumpJSON(os.path.join(output_dir, "individual_profiles", "individual_energies_contribution_" + motif["pdb_name"] + "_" + pot_tag + ".json"), contribution_per_individual_nucleotide)
                
                # Append the motif to the motif_list #
                motif_list.append(motif)
            
            # Get the clusters, it is a list of lists where each cluster is a list of motif dictionaries #
            motifs_clustered = cluster_motifs.get_motif_clusters(motif_list, pval_cutoffs, equal_bs=False)
            for i in range(0, len(motifs_clustered)):
                cluster = motifs_clustered[i]
                
                # Initializing data lists #
                match_list = []
                score_list = []
                logpval_list = []
                energy_list = []
                energy_raw_list = []
                i_energy_list = []
                perc_energy_list = []
                
                # Iterate the motifs in the cluster #
                for motif in cluster:
                    match_list.append(motif["matches_per_nucleotide"])
                    score_list.append(motif["scores_per_nucleotide"])
                    logpval_list.append(motif["logpval_per_nucleotide"])
                    energy_list.append(motif["energies_per_nucleotide"])
                    energy_raw_list.append(motif["energies_per_nucleotide_raw"])
                    i_energy_list.append(motif["energies_per_individual_nucleotide"])
                    perc_energy_list.append(motif["contribution_per_individual_nucleotide"])
                
                # Merge data from different models #
                matches_per_nucleotide_global = merge_data_per_nucleotide(match_list, "add", pval_cutoffs)
                scores_per_nucleotide_global = merge_data_per_nucleotide(score_list, "merge", pval_cutoffs)
                logpval_per_nucleotide_global = merge_data_per_nucleotide(logpval_list, "merge", pval_cutoffs)
                energies_per_nucleotide_global = merge_data_per_nucleotide(energy_list, "merge", pval_cutoffs)
                energies_per_nucleotide_raw_global = merge_data_per_nucleotide(energy_raw_list, "merge", pval_cutoffs)
                energies_per_individual_nucleotide_global = merge_data_per_nucleotide(i_energy_list, "merge", pval_cutoffs)
                contribution_per_individual_nucleotide_global = merge_data_per_nucleotide(perc_energy_list, "merge", pval_cutoffs)
                
                # Average and scale data #
                matches_per_nucleotide_global = scale_to_percentage(matches_per_nucleotide_global)
                scores_per_nucleotide_global = average_scores(scores_per_nucleotide_global, pval_cutoffs)
                logpval_per_nucleotide_global = average_scores(logpval_per_nucleotide_global, pval_cutoffs)
                energies_per_nucleotide_global = average_scores(energies_per_nucleotide_global, pval_cutoffs)
                energies_per_nucleotide_raw_global = average_scores(energies_per_nucleotide_raw_global, pval_cutoffs)
                energies_per_individual_nucleotide_global = average_scores(energies_per_individual_nucleotide_global, pval_cutoffs)
                contribution_per_individual_nucleotide_global = average_scores(contribution_per_individual_nucleotide_global, pval_cutoffs)
                
                # Store the data in jsons #
                functions.dumpJSON(os.path.join(output_dir, "cluster_profiles", "FIMO_matches_cluster" + str(i) + "_" + pot_tag + ".json"), matches_per_nucleotide_global)
                functions.dumpJSON(os.path.join(output_dir, "cluster_profiles", "FIMO_scores_cluster" + str(i) + "_" + pot_tag + ".json"), scores_per_nucleotide_global)
                functions.dumpJSON(os.path.join(output_dir, "cluster_profiles", "FIMO_logpval_cluster" + str(i) + "_" + pot_tag + ".json"), logpval_per_nucleotide_global)
                functions.dumpJSON(os.path.join(output_dir, "cluster_profiles", "energies_cluster" + str(i) + "_" + pot_tag + ".json"), energies_per_nucleotide_global)
                functions.dumpJSON(os.path.join(output_dir, "cluster_profiles", "energies_raw_cluster" + str(i) + "_" + pot_tag + ".json"), energies_per_nucleotide_raw_global)
                functions.dumpJSON(os.path.join(output_dir, "cluster_profiles", "individual_energies_cluster" + str(i) + "_" + pot_tag + ".json"), energies_per_individual_nucleotide_global)
                functions.dumpJSON(os.path.join(output_dir, "cluster_profiles", "individual_energies_contribution_cluster" + str(i) + "_" + pot_tag + ".json"), contribution_per_individual_nucleotide_global)
                
                # Write the cluster's file #
                clusters_file = os.path.join(output_dir, "clusters.txt")
                cf = open(clusters_file, "ab")
                mot_in_cluster = []
                for motif in cluster:
                    mot_in_cluster.append(motif["pdb_name"])
                cf.write("Cluster " + str(i) + ": " + ";".join(mot_in_cluster) + "\n")
                cf.close()


if __name__ == "__main__":
    main()

