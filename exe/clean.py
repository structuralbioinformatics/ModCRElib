"""
PDB sanitization helpers for ModCRE preprocessing.

The script removes problematic/discontinuous DNA segments and overlapping
nucleotide chains while preserving protein chains, then writes a cleaned PDB.
"""

import os, sys, re
import configparser
import optparse

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

# Imports jbonet's module #
from SBILib.structure import PDB
from SBILib.structure.chain import ChainOfNucleotide

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse CLI arguments for cleaning a PDB file.

    How to run:
        ``python clean.py -i input.pdb [--dummy /tmp -o clean.pdb]``

    Returns:
        optparse.Values: Namespace with ``input_file``, ``dummy_dir``, ``output_file``.
    """

    parser = optparse.OptionParser("python clean.py -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", default="clean.pdb", action="store", type="string", dest="output_file", help="Output file (default = clean.pdb)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_clean_pdb_obj(pdb_file, dummy_dir="/tmp"):
    """
    Load and sanitize a PDB structure for downstream DNA/protein workflows.

    Args:
        pdb_file (str): Input structure file.
        dummy_dir (str): Reserved for compatibility (unused in current body).

    Returns:
        PDB: Cleaned structure object.

    Raises:
        ValueError: If the input cannot be parsed or cleaned.
    """

    try:
        # Get PDB object #
        pdb_obj = PDB(pdb_file)
        # Remove discontinued nucleotides and overlapped DNA chains; it interferes with 3DNA #
        pdb_obj = remove_discontinued_nucleotides_and_overlapped_chains(pdb_obj)
    except:
        raise ValueError("Could not clean PDB %s" % pdb_file)

    return pdb_obj

def remove_discontinued_nucleotides_and_overlapped_chains(pdb_obj):
    """
    Keep protein chains and one non-overlapping continuous block per DNA chain.

    DNA chains are split into continuous nucleotide blocks, ranked by block length,
    and accepted if they do not geometrically overlap previously accepted residues.

    Args:
        pdb_obj (PDB): Original structure.

    Returns:
        PDB: Sanitized structure preserving proteins and filtered DNA chains.
    """

    # Initialize #
    blocks = {}
    chains = {}
    nucleotides = []
    max_overlap_distance = float(config.get("Parameters", "max_overlap_distance"))
    clean_pdb_obj = PDB()

    # Add protein chains #
    for protein_chain in pdb_obj.proteins:
        clean_pdb_obj.add_chain(protein_chain)
    
    # For each DNA chain... #
    for dna_chain in pdb_obj.nucleotides:
        # Initialize #
        blocks.setdefault((dna_chain.pdb, dna_chain.chain), 0)
        chains.setdefault((dna_chain.pdb, dna_chain.chain), [])
        chains[(dna_chain.pdb, dna_chain.chain)].append([])
        chains[(dna_chain.pdb, dna_chain.chain)][-1].append(dna_chain.nucleotides[0])
        # For each residue... #
        for i in range(1, len(dna_chain.nucleotides)):
            if not chains[(dna_chain.pdb, dna_chain.chain)][-1][-1].is_followed(dna_chain.nucleotides[i]):
                # Chain is discontinued! #
                chains[(dna_chain.pdb, dna_chain.chain)].append([])
            # Add nucleotide to chain #
            chains[(dna_chain.pdb, dna_chain.chain)][-1].append(dna_chain.nucleotides[i])
            # Chain largest block length #
            if len(chains[(dna_chain.pdb, dna_chain.chain)][-1]) > blocks[(dna_chain.pdb, dna_chain.chain)]:
                blocks[(dna_chain.pdb, dna_chain.chain)] = len(chains[(dna_chain.pdb, dna_chain.chain)][-1])
    # For each DNA chain... #
    for pdb, pdb_chain in sorted(blocks, key=lambda x: (blocks[x] * -1, x)):
        # For continuous set of nucleotides... #
        for continuous_nucleotides in sorted(chains[(pdb, pdb_chain)], key=lambda x: len(x), reverse=True):
            # Initialize #
            is_overlapped = False
            # For continuous nucleotide... #
            for continuous_nucleotide in continuous_nucleotides:
                # Skip if is overlapped #
                if is_overlapped: continue
                # For nucleotide in DNA chain... #
                for nucleotide in nucleotides:
                    # If residue overlaps... #
                    a, b, distance = continuous_nucleotide.distance(nucleotide, "geometric")
                    if distance < max_overlap_distance:
                        is_overlapped = True
                        break
            # If non-overlapped continuous set of nucleotides... #
            if not is_overlapped:
                chain_obj = ChainOfNucleotide(pdb, pdb_chain)
                for continuous_nucleotide in continuous_nucleotides:
                    chain_obj.add_residue(continuous_nucleotide)
                    nucleotides.append(continuous_nucleotide)
                clean_pdb_obj.add_chain(chain_obj)
                break

    return clean_pdb_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get clean PDB object #
    pdb_obj = get_clean_pdb_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Create output #
    pdb_obj.write(os.path.abspath(options.output_file), force=True, clean=True)
