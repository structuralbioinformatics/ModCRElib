"""
Query TFinDit and infer TF/TF-homolog status per PDB entry.

Downloads (and caches) per-PDB TFinDit pages, then recursively follows homolog
links to classify entries as direct TF, TF-homolog, or unclassified.
"""

import os, sys, re
import configparser
import optparse
from time import sleep
import urllib.request, urllib.error, urllib.parse

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

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    Parse CLI options for TFinDit classification.

    How to run:
        ``python tfindit.py --pdb pdb_dir -o output_dir [-v]``

    Returns:
        optparse.Values: Namespace with PDB folder, output folder, and verbosity.
    """

    parser = optparse.OptionParser("python tfindit.py -p pdb_dir [--dummy dummy_dir -o output_dir -t tfindit_file -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (where all PDBs are placed; e.g. \"/db/rcsb/pdb-remediated.931/data/structures/all/pdb/\")", metavar="{directory}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")

    (options, args) = parser.parse_args()

    if options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_tfindit_file(pdb_name, filename):
    """
    Download one TFinDit page and persist it locally.

    Args:
        pdb_name (str): Four-letter PDB id.
        filename (str): Destination text file.
    """
    # Initialize #
    content = get_tfindit_content(pdb_name)

    if content is not None:
        functions.write(filename, content.read())

def get_tfindit_content(pdb_name):
    """
    Retrieve raw TFinDit HTTP content with retry logic.

    Args:
        pdb_name (str): Four-letter PDB id.

    Returns:
        HTTPResponse | None: Open URL response or ``None`` after retries.
    """
    # Initialize #
    content = None
    trials = 0

    while content is None:
        if trials == 5:
            break
        try:
            content = urllib.request.urlopen("http://bioinfozen.uncc.edu/tfindit/index.php?pdbid=%s&tab=pdb&seqid=95&evalue=1.0e-08&hitcov=0&querycov=0&seqid2=95&evalue2=1.0e-08&hitcov2=0&querycov2=0" % pdb_name, timeout=1)
        except:
            # Sleep and try again... #
            trials += 1
            sleep(5)

    return content

def recursive_search_tf_homolog_pdb_chains_in_tfindit(pdb_name, file_name, level, fathers=set()):
    """
    Recursively inspect TFinDit homolog links to assign TF status.

    Args:
        pdb_name (str): Current PDB id.
        file_name (str): Cached TFinDit file for ``pdb_name``.
        level (int): Recursion depth (0 = direct TF check).
        fathers (set): Already-visited PDB ids.

    Returns:
        str | None: ``"TF"``, ``"TF-homolog"``, or ``None``.
    """
    # Initialize #
    search = set()

    if level <= 1:
        # For each line... #
        for line in functions.parse_file(file_name):
            # This is a TF... #
            m = re.search("<td>chains (\S+): .+</td>", line)
            if m:
                if level == 0: return "TF"
                else: return "TF-homolog"
            if level == 0:
                for m in re.finditer("<td><font size='2'>(\S{4})\S</font></td><td><font size='2'>(\S{4})\S</font></td>", line):
                    search.add(m.group(1))
                    search.add(m.group(2))
        for search_pdb_name in search:
            if search_pdb_name == pdb_name: continue
            if search_pdb_name in fathers: continue
            search_file_name = file_name.replace(pdb_name, search_pdb_name)
            if not os.path.exists(search_file_name):
                get_tfindit_file(search_pdb_name, search_file_name)
            if recursive_search_tf_homolog_pdb_chains_in_tfindit(search_pdb_name, search_file_name, level + 1, fathers.union(search)) is not None:
                return "TF-homolog"

    return None

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Step 1) Parse options and ensure output folder exists.
    options = parse_options()

    # Create output directory #
    if options.output_dir is not None:
        if not os.path.exists(options.output_dir):
            os.makedirs(options.output_dir)
 
    # Step 2) Download/cache TFinDit pages for all PDB entries.
    for pdb_file in os.listdir(options.pdb_dir):
        # Initialize #
        pdb_name = pdb_file[:4].lower()
        # Skip if TFinDit file already exists #
        tfindit_file = os.path.abspath(os.path.join(options.output_dir, pdb_name + ".txt"))
        if os.path.exists(tfindit_file): continue
        # Verbose mode #
        if options.verbose: sys.stdout.write("%s...\n" % pdb_name)
        # Get TFinDit file #
        get_tfindit_file(pdb_name, tfindit_file)

    # Step 3) Classify each PDB recursively as TF / TF-homolog / none.
    tfindit = {}
    # For each PDB entry... #
    for pdb_file in os.listdir(options.pdb_dir):
        # Initialize #
        pdb_name = pdb_file[:4].lower()
        # Verbose mode #
        if options.verbose: sys.stdout.write("%s...\n" % pdb_name)
        # Get TFinDit file #
        tfindit_file = os.path.abspath(os.path.join(options.output_dir, pdb_name + ".txt"))
        # Recursively find if is TF #
        tfindit.setdefault(pdb_name, recursive_search_tf_homolog_pdb_chains_in_tfindit(pdb_name, tfindit_file, 0))

    # Step 4) Export non-empty classifications.
    output_file = os.path.join(os.path.abspath(options.output_dir), "tfindit.txt")
    functions.write(output_file, "#PDB;status")
    # For each PDB entry... #
    for pdb_name in sorted(tfindit):
        if tfindit[pdb_name] is not None:
            functions.write(output_file, "%s;%s" % (pdb_name, tfindit[pdb_name]))

