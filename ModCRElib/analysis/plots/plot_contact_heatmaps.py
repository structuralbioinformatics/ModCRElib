import os, sys, re
import configparser
import optparse
import shutil
import subprocess
import difflib
import collections
import math
import hashlib
import pickle
import scipy
from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, text
import pandas as pd
import seaborn as sns
import time
import matplotlib.ticker as tkr

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
from SBILib.structure         import PDB
from SBILib.structure.chain   import ChainOfProtein, ChainOfNucleotide
from SBILib.structure.residue import ResidueOfNucleotide
from SBILib.structure.atom    import AtomOfNucleotide
from SBILib.data              import nucleic1to3, dna_complementary
from SBILib.external.blast    import blast_parser
from SBILib.data import aminoacids_polarity_boolean, aminoacids3to1, nitrogenous_bases

# Imports my  modules #
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.potential import spotentials

# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")


def plot_contact_heatmap(frequencies, distances, label, output_dir, tope, detail):
    """
    Render contact-frequency heatmaps for amino-acid/dinucleotide environments.

    Args:
        frequencies (dict): Contact counts keyed by triad feature labels.
        distances (list[float]): Distance bins used to build contacts.
        label (str): Label used in output filenames.
        output_dir (str): Destination folder for generated plots.
        tope (int | None): Optional upper cap for color scaling.
        detail (bool): Whether to generate per-AA/per-dinucleotide detail plots.

    Returns:
        None.
    """

    # Get all possible ticks for amino acids and for dinucleotides #
    aa_ticks = []
    # Amino acids are sorted according to the hydrophobicity scale of Kyte and Doolittle. A simple method for displaying the hydropathic character of a protein. Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32. #
    aa = ["ILE", "VAL", "LEU", "PHE", "CYS", "MET", "ALA", "GLY", "THR", "SER", "TRP", "TYR", "PRO", "HIS", "GLU", "GLN", "ASP", "ASN", "LYS", "ARG"]
    exposure = ["E", "B"]
    sstructure = ["E", "H", "C"]

    dn_ticks = []
    nn = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    stra = ["F", "R"]
    grov = ["A", "I"]
    chem = ["B", "N"]

    for a in sorted(aa):
        num = 0
        for key in list(frequencies.keys()):
            if a in key:
                num += frequencies[key]
    
    for a in sorted(aa):
        if aminoacids_polarity_boolean[aminoacids3to1[a]] == True:
            hidro = "P"
        else:
            hidro = "N"
        for e in exposure:
            for s in sstructure:
                aa_ticks.append(a + "-" + hidro + "-" + e + "-" + s)

    for n in sorted(nn):
        nb = "".join(nitrogenous_bases[nucleotide] for nucleotide in n)
        for s in stra:
            for g in grov:
                for c in chem:
                    dn_ticks.append(n + "-" + nb + "-" + s + "-" + g + "-" + c)
    
    matrix = []
    for nnt in dn_ticks:
        row = []
        for aat in aa_ticks:
            if aat + ";" + nnt in list(frequencies.keys()):
                #
                if tope == None:
                    row.append(frequencies[aat + ";" + nnt])
                else:
                    if frequencies[aat + ";" + nnt] > tope:
                        row.append(np.log(tope))
                    else:
                        row.append(np.log(frequencies[aat + ";" + nnt]))
            else:
                row.append(np.nan)
        matrix.append(row)

    y_lines = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]
    x_lines = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120]

    reduced_aa_ticks = []
    reduced_dn_ticks = []
    for a in sorted(aa):
        reduced_aa_ticks += ["", "", " " + aminoacids3to1[a], "", "", ""]

    for n in sorted(nn):
        reduced_dn_ticks += ["", "  ", "   ", "  " + n, "", "", "", ""]

    df = pd.DataFrame(data=np.array(matrix), index=reduced_dn_ticks, columns=reduced_aa_ticks)
    #print(df.to_string())
    fig = figure(figsize=(50, 45))
    sns.set(font_scale=6)
    mask = df.isnull()

    formatter = tkr.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    if tope != None:
        ax = sns.heatmap(df, annot=False, cmap="Reds", cbar=True, annot_kws={"size": 18}, cbar_kws={"format": formatter}, mask=mask, vmin=0.0, vmax=np.log(tope))
    else:
        ax = sns.heatmap(df, annot=False, cmap="Reds", cbar=True, annot_kws={"size": 18}, cbar_kws={"format": formatter}, mask=mask, vmin=0.0)

    for item in ax.get_xticklabels():
        item.set_fontsize(80)
        item.set_fontweight('bold')
    for item in ax.get_yticklabels():
        item.set_fontsize(80)
        item.set_fontweight('bold')

    ax.hlines(y_lines, ax.get_xlim()[0], ax.get_xlim()[1])
    ax.vlines(x_lines, ax.get_ylim()[0], ax.get_ylim()[1])

    plt.subplots_adjust(left=0.15)
    fig.savefig(os.path.join(output_dir, "contacts_" + label.replace("/", "-") + "_" + str(tope) + ".png"))
    print(("plot created at: " + str(os.path.join(output_dir, "contacts_" + label.replace("/", "-") + "_" + str(tope) + ".png"))))
    plt.close(fig)

    if detail == True:

        for a in sorted(aa):
            for n in sorted(nn):
                matrix = []
                subplot_dt = []
                for nnt in dn_ticks:
                    if not nnt.startswith(n):
                        continue
                    row = []
                    subplot_dt.append(nnt)
                    subplot_at = []
                    for aat in aa_ticks:
                        if not aat.startswith(a):
                            continue
                        subplot_at.append(aat)
                        if aat + ";" + nnt in list(frequencies.keys()):
                            if tope == None:
                                row.append(frequencies[aat + ";" + nnt])
                            else:
                                if frequencies[aat + ";" + nnt] > tope:
                                    row.append(np.log(tope))
                                else:
                                    row.append(np.log(frequencies[aat + ";" + nnt]))
                        else:
                            row.append(np.nan)
                    matrix.append(row)

                df = pd.DataFrame(data=np.array(matrix), index=subplot_dt, columns=subplot_at)

                fig = figure(figsize=(20, 20))
                sns.set(font_scale=6)
                mask = df.isnull()

                formatter = tkr.ScalarFormatter(useMathText=True)
                formatter.set_scientific(True)
                formatter.set_powerlimits((-2, 2))
                if tope != None:
                    ax = sns.heatmap(df, annot=False, cmap="Reds", cbar=True, annot_kws={"size": 18}, cbar_kws={"format": formatter}, mask=mask, vmin=0.0, vmax=np.log(tope))
                else:
                    ax = sns.heatmap(df, annot=False, cmap="Reds", cbar=True, annot_kws={"size": 18}, cbar_kws={"format": formatter}, mask=mask, vmin=0.0)

                for item in ax.get_xticklabels():
                    item.set_rotation(90)
                    item.set_fontsize(70)
                    item.set_fontweight('bold')
                for item in ax.get_yticklabels():
                    item.set_rotation(0)
                    item.set_fontsize(70)
                    item.set_fontweight('bold')

                plt.subplots_adjust(left=0.35, bottom=0.35)
                fig.savefig(os.path.join(output_dir, a + "_" + n + "_" + label.replace("/", "-") + ".png"))
                print(("plot created at: " + str(os.path.join(output_dir, a + "_" + n + "_" + label.replace("/", "-") + ".png"))))
                plt.close(fig)
    




def parse_triads_files(input_file, only_pdb=False):
    """
    Parse triad potentials and collapse frequencies into scalar totals.

    Args:
        input_file (str): Triads input file.
        only_pdb (bool): Reserved compatibility flag (currently unused).

    Returns:
        dict: Feature-to-frequency mapping.
    """

    
    f_dab, f_a_dab, f_a_dab_oa, f_a_b_dab, f_dab_oa_ob, f_a_b_dab_oa_ob = spotentials.get_frequencies(input_file, computation=True,approach=False)
    new_frequencies = {}
    for key in f_a_b_dab_oa_ob:
        new_frequencies[key] = sum(f_a_b_dab_oa_ob[key])

    return new_frequencies

    


def parse_options():
    """
    Parse CLI options for contact heatmap generation.

    Returns:
        optparse.Values: Namespace with input/output controls and plot settings.
    """

    parser = optparse.OptionParser("tf_mutations_analyzer.py --input_tf= --mutation= --dummy= --pbm= --pdb=")

    # Directory arguments  
    parser.add_option("-o", "--output_dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="OUTPUT_DIR")
    parser.add_option("-i", "--input_dir", default="./", action="store", type="string", dest="input_dir", help="input_dir", metavar="INPUT_DIR")
    parser.add_option("--dummy", default="dummy/", action="store", type="string", dest="dummy_dir", help="Dummy_dir", metavar="DUMMY_DIR")    
    parser.add_option("-l", "--label", default="", action="store", type="string", dest="label", help="label", metavar="LABEL") 
    parser.add_option("-m", "--max", default=10000, action="store", type="int", dest="max", help="maximum number of contacts counted by cuadrant", metavar="MAX") 
    parser.add_option("-d", "--detail", default=False, action="store_true", dest="detail", help="make detailed plots for each amoni acid - dinucleotide combination", metavar="DETAIL") 
    parser.add_option("-f", "--families", default=False, action="store", type="string", dest="family_file", help="family_file", metavar="FAMILY_FILE") 
    parser.add_option("--frequencies", default=False, action="store_true", dest="frequencies", help="use precomputed frequencies pickles to make the plots", metavar="FREQUENCIES_PICKLE")
    
    (options, args) = parser.parse_args()

    return options


#-------------#
# Main        #
#-------------#
if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    output_dir = os.path.abspath(options.output_dir)
    input_dir = os.path.abspath(options.input_dir)
    dummy_dir = os.path.abspath(options.dummy_dir)
    bin_distance = float(config.get("Parameters", "bin_distance_bins"))
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    distances = list(np.arange(0, max_contact_distance + bin_distance, bin_distance))
    if not os.path.exists(output_dir):os.makedirs(output_dir)
    if not os.path.exists(dummy_dir):os.makedirs(dummy_dir)

    # Get contact frequencies #
    if options.frequencies == False:
        # Parse the families file #
        families = {}
        for line in functions.parse_file(os.path.abspath(options.family_file)):
            pdb_chain, family = line.split(";")
            families[pdb_chain] = family.replace("/","-")
        # Iterate nr files #
        for nr in os.listdir(input_dir):
            if nr.endswith(".txt"):
                input_file = os.path.join(input_dir, nr)
                frequencies = parse_triads_files(input_file)
                fold = nr.split(".")[0]
                try:
                    family = families[fold]
                except:
                    family = "Unknown"
                if options.label == None:
                    label = fold + "_" + family
                else:
                    label = fold + "_" + family + "_" + options.label

                pickle.dump(frequencies, open(os.path.join(output_dir, "frequencies_" + label + ".p"), "wb"))
                # Make plots #
                plot_contact_heatmap(frequencies=frequencies, distances=distances, label=label, output_dir=output_dir, tope=options.max, detail=options.detail)

    # Parse contact frequencies #
    if options.frequencies == True:
        for fr in os.listdir(input_dir):
            if fr.endswith(".p"):
                frequencies_pickle = os.path.join(input_dir, fr)
                frequencies = pickle.load(open(frequencies_pickle, "rb"))
                if options.label == "":
                    label = "_".join(fr.split(".")[0].split("_")[1:4])
                else:
                    label = "_".join(fr.split(".")[0].split("_")[1:4]) + "_" + label
                # Make plots #
                plot_contact_heatmap(frequencies=frequencies, distances=distances, label=label, output_dir=output_dir, tope=options.max, detail=options.detail)

