import os, sys, re
import configparser
import copy
import json
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import spotentials as potentials
import optparse
import pandas as pd
import datetime
import seaborn as sns
from scipy.stats import gaussian_kde
import itertools

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

# Import my modules #
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.potential import spotentials

# Import jbonet's module #
from SBILib.data import aminoacids3to1, aminoacids_polarity_boolean


def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python plot_potentials.py -i input_file [-a --dummy=dummy_dir -o output_file ]")
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Input: folder of potentials file", metavar="{filename}")
    parser.add_option("-f", action="store", type="string", dest="fold", help="Input: PDB chain to get the corresponding potentials files", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output file (default = stdout)", metavar="{filename}")
    parser.add_option("--methylation",default=False, action="store_true", dest="methylation", help="Flag to use methylated cytosines with binding/non-binding specificity (default=False)")
    parser.add_option("--pmf",default=False, action="store_true", dest="pmf", help="Flag to use Potentials of Mean Force without Zscore normalization (default=False)")



    (options, args) = parser.parse_args()
    return options



if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    import pandas as pd
    import numpy as np
    import datetime
    import pickle
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    import pylab as pl
    from pylab import plot, show, savefig, xlim, figure,  ylim, legend, boxplot, setp, axes, text
    from matplotlib import lines
    from scipy import stats
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    nucleotides = list("ACGT")
    if options.methylation:
       nucleotides.extend(list("XJOQ"))
    distances     = [i for i in range(50)]
    dinucleotides = ["".join(x) for x in itertools.product(nucleotides, repeat=2)]
    #dinucleotides = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    # Since there are more families (and general as well) that work better with accumulative potentials, we will make this plots with accumulative potentials #    
   
    if not os.path.exists(options.output_dir):
       os.makedirs(options.output_dir)
     
    if options.pmf:
      acc = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".pmf.acc.txt"))
      bins = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".pmf.bins.txt"))
      acc_taylor = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".pmf.taylor.acc.txt"))
      bins_taylor = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".pmf.taylor.bins.txt"))
    else:
      acc = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".acc.txt"))
      bins = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".bins.txt"))
      acc_taylor = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".taylor.acc.txt"))
      bins_taylor = spotentials.Potentials(file_name=os.path.join(os.path.abspath(options.input_dir), options.fold + ".taylor.bins.txt"))
    
    # By default we work with pair potentials, but this can be changed by selecting the type of pmf at this point of the code #
    acc_pair = acc._pmf_pair
    bins_pair = bins._pmf_pair
    acc_taylor_pair = acc_taylor._pmf_pair
    bins_taylor_pair = bins_taylor._pmf_pair
    
    acc_distances = sorted(acc._distances.keys())
    bins_distances = sorted(bins._distances.keys())
    
    for aa in aminoacids:
        
        for i in range(0, len(dinucleotides)):
            dn = dinucleotides[i]
            fig, s_plots = plt.subplots(nrows=1, ncols=1, figsize=(16, 10))
            s_plots.patch.set_facecolor("white")
            s_plots.grid(color="lightgray")
            # Get the data #
            
            if not str(aa+";"+dn) in list(acc_pair.keys()):
                acc_pair.setdefault(str(aa+";"+dn), [None]*len(distances))
                bins_pair.setdefault(str(aa+";"+dn), [None]*len(distances))
                acc_taylor_pair.setdefault(str(aa+";"+dn), [None]*len(distances))
                bins_taylor_pair.setdefault(str(aa+";"+dn), [None]*len(distances))

            pair_acc = acc_pair[aa+";"+dn]
            pair_bins = bins_pair[aa+";"+dn]
            pair_taylor_acc = acc_taylor_pair[aa+";"+dn]
            pair_taylor_bins = bins_taylor_pair[aa+";"+dn]
            # Plot the data #
            try:
              s_plots.plot(acc_distances, pair_acc, color="firebrick", linewidth="5")
              s_plots.plot(bins_distances, pair_bins, color="orange", linewidth="5")
              s_plots.plot(acc_distances, pair_taylor_acc, color="blueviolet", linewidth="5")
              s_plots.plot(bins_distances, pair_taylor_bins, color="dodgerblue", linewidth="5")
              # Scale axis #
              s_plots.set_xlim(0.0, 30.0)
              s_plots.set_ylim(-5.0, 5.0)
              plt.xticks(fontsize=40)
              plt.yticks(fontsize=40)
              plt.ylabel("Score", fontsize=40, fontweight="bold")
              plt.xlabel("Distance", fontsize=40, fontweight="bold")
              #plt.legend(handles=(handle_list), bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., prop={'size':25})
              plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)
              fig.suptitle(options.fold + "  " + aa + "  " + dn, fontsize=25)
              if options.pmf:
                fig.savefig(os.path.join(os.path.abspath(options.output_dir), options.fold + "_" + aa + "_" + dn + ".pmf.png")) 
                print(("plot created in: " + os.path.join(os.path.abspath(options.output_dir), options.fold + "_" + aa + "_" + dn + ".pmf.png")))
              else:
                fig.savefig(os.path.join(os.path.abspath(options.output_dir), options.fold + "_" + aa + "_" + dn + ".png")) 
                print(("plot created in: " + os.path.join(os.path.abspath(options.output_dir), options.fold + "_" + aa + "_" + dn + ".png")))
            except Exception as e:
                print(("Error: %s %s %s"%(aa,dn,e)))
                continue
