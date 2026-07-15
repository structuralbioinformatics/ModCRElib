import os, sys, re
import configparser
import optparse
import subprocess
import shutil
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Residue import Residue
from Bio.PDB import Superimposer
from copy import deepcopy

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
from SBILib.data import dna_complementary
from SBILib.structure import PDB
from SBILib.structure.chain import ChainOfNucleotide
from SBILib.structure.residue import ResidueOfNucleotide

# Import my modules #
from ModCRElib.structure.contacts import interface
from ModCRElib.structure.dna import x3dna



template_dir = os.path.join(
    scripts_path,
    "ModCRElib",
    "data",
    "dna_templates"
)



TEMPLATES = {}

parser = PDBParser(QUIET=True)

for name in ["DA","DC","DG","DT"]:

    s = parser.get_structure(
        name,
        os.path.join(template_dir, name + ".pdb")
    )

    TEMPLATES[name] = next(s.get_residues())

BACKBONE = {
    "P","OP1","OP2",
    "O5'","C5'",
    "C4'","O4'",
    "C3'","O3'",
    "C2'","C1'"
}



#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for protein-DNA complex DNA remodeling.
    """
    parser = optparse.OptionParser("python model_dna.py  -p pdb_file ( -i interface_file -p pdb_file and (-s dna_sequence or (-t threading_file --pdb pdb_dir ) ) ][--dummy=dummy_dir -o output_dir]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", default=None, type="string", dest="interface_file", help="Interface file (from interface.py)", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-p", action="store", type="string", dest="pdb_file", help="PDB file (e.g. from model_protein.py)", metavar="{filename}")
    parser.add_option("-s", action="store", type="string", dest="dna_sequence", help="DNA sequence ", metavar="{string}")
    parser.add_option("-t", action="store", default=None, type="string", dest="threading_file", help="Threading file (e.g. from threader.py)", metavar="{string}")
    parser.add_option("--pdb", action="store", default=None, type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="PDB_DIR")

    (options, args) = parser.parse_args()

    if (options.interface_file is None or options.pdb_file is None or (options.dna_sequence is None and (options.threading_file is None or options.pdb_dir is None)) ):
        parser.error("missing arguments: type option \"-h\" for help")

    return options


def mutate_residue(residue, new_base):
    """
    Replace the nucleobase of a DNA residue while keeping the original
    sugar-phosphate backbone coordinates.

    Parameters
    ----------
    residue : Bio.PDB.Residue
        Residue in the original structure to mutate.

    new_base : str
        One of "A", "C", "G", "T".
    """

    template_name = {
        "A": "DA",
        "C": "DC",
        "G": "DG",
        "T": "DT"
    }[new_base]

    # Make a fresh copy so templates never get modified
    template = TEMPLATES[template_name].copy()

    # Backbone atoms used for alignment
    ALIGN = [
        "C1'",
        "C2'",
        "C3'",
        "C4'",
        "O4'"
    ]

    fixed = [residue[a] for a in ALIGN]
    moving = [template[a] for a in ALIGN]

    sup = Superimposer()
    sup.set_atoms(fixed, moving)

    # Rotate/translate the template into the residue frame
    sup.apply(list(template.get_atoms()))

    # Remove every atom that is NOT backbone
    for atom in list(residue):
        if atom.get_name() not in BACKBONE:
            residue.detach_child(atom.id)

    # Copy in the new base atoms
    for atom in template:
        if atom.get_name() in BACKBONE:
            continue
        residue.add(atom.copy())

    residue.resname = template.resname

def get_dna_model_pdb_obj(pdb_file, dna_sequence, x3dna_obj, interface_obj, interface_range=None, dummy_dir="/tmp"):
    """
    Mutate interface DNA bases using a native implementation of the 
    find_pair -> analyze -> parameter modification -> rebuild pipeline.
    """
    basepairs = set()
    nucleobases = set()
    model_pdb_file = os.path.join(dummy_dir, str(os.getpid()) + "mutated.pdb")
    
    # Text mutation tracking dictionary: { basepair_index: "A-T" }
    target_mutations = {}

    # Correct DNA sequence #
    dna_sequence = re.sub("U", "T", dna_sequence)
    dna_sequence = re.sub("Q", "T", dna_sequence)
    dna_sequence = re.sub("J", "T", dna_sequence)
    dna_sequence = re.sub("X", "C", dna_sequence)
    dna_sequence = re.sub("O", "C", dna_sequence)
    dna_sequence = re.sub("[^acgtACGT]", "N", dna_sequence)
    
    # Get nucleotides list #
    nucleotides = list(dna_sequence)
    
    if interface_obj.get_interface_basepairs() is not None:
        interface_range = interface_obj.get_interface_basepairs()
    elif interface_range is None:
        sys.stdout.write("[WARNING] No DNA interface to model. Original DNA preserved.\n")
        return PDB(pdb_file)


    parser = PDBParser(QUIET=True)
    
    structure = parser.get_structure(
        "original",
        pdb_file
    )
    residue_lookup = {}
    
    for chain in structure[0]:
    
        for residue in chain:
    
            residue_lookup[(chain.id, residue.id[1])] = residue


    for basepair in interface_range:
    
        if not nucleotides:
            break
    
        nucleotide = nucleotides.pop(0).upper()
    
        ((chain1, res1), (chain2, res2)) = x3dna_obj.get_basepair(basepair)
    
        if nucleotide == "N":
            continue
    
        mutate_residue(
            residue_lookup[(chain1, res1)],
            nucleotide
        )
    
        mutate_residue(
            residue_lookup[(chain2, res2)],
            dna_complementary[nucleotide]
        )



    clean_pdb = os.path.join(
        dummy_dir,
        f"{os.getpid()}_mutated.pdb"
    )

    io = PDBIO()            
    io.set_structure(structure)
    io.save(clean_pdb)
    model_pdb_obj = PDB(clean_pdb)









    ###############################################################################
    # Cleanup
    ###############################################################################
    
    for f in (model_pdb_file, clean_pdb):
        if os.path.exists(f):
            os.remove(f)
    
    return model_pdb_obj


def main():
    options = parse_options()

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    if not os.path.exists(options.dummy_dir):
        os.makedirs(options.dummy_dir)
        
    pdb_file = options.pdb_file
    if not pdb_file.startswith("/"): pdb_file = os.path.abspath(pdb_file)
    
    pdb_dir = options.pdb_dir
    if pdb_dir is not None and not pdb_dir.startswith("/"): 
        pdb_dir = os.path.abspath(options.pdb_dir)
        
    if options.threading_file is not None:
       threading_file = options.threading_file
       if not threading_file.startswith("/"): threading_file = os.path.abspath(threading_file)
       
    if options.interface_file is not None:
       interface_file = options.interface_file
       if not interface_file.startswith("/"): interface_file = os.path.abspath(interface_file) 

    dna = set()
    kmers = {}
    x3dna_obj = x3dna.get_x3dna_obj(pdb_file)

    if options.threading_file is not None and options.pdb_dir is not None:
        if os.path.exists(threading_file) and os.path.exists(pdb_dir):
            import threader 
            thread_obj = threader.Threaded(threading_file)
            pdb_name = thread_obj.get_pdb_name()
            pdb_chain = thread_obj.get_pdb_chain()
            interface_file = os.path.join(pdb_dir, "interfaces", pdb_name + "_" + pdb_chain + ".txt")
            kmers = thread_obj.get_kmers()
            kmers_fixed = thread_obj.get_kmers_fixed()
            for dna_sequence in kmers.keys():
                if "N" not in dna_sequence: dna.add(dna_sequence)
            for dna_sequence in kmers_fixed.keys():
                dna.add(dna_sequence)
        else:
            sys.stderr.write("Threading file %s is not found\n" % (threading_file))
            exit(0)
    else:
        kmers.setdefault(options.dna_sequence, 0)
        dna.add(options.dna_sequence)
        
    if os.path.exists(interface_file):
        interface_obj = interface.Interface(interface_file)
    else:
        sys.stderr.write("Interface file %s is not found\n" % (interface_file))
        exit(0)

    for dna_sequence in dna:
        if len(dna_sequence) != interface_obj.get_interface_length():
            raise ValueError("DNA sequence does not match the interface length.")

        try:
            pdb_obj = get_dna_model_pdb_obj(pdb_file, dna_sequence, x3dna_obj, interface_obj, interface_range=None, dummy_dir=options.dummy_dir)
            pdb_obj.write(os.path.join(options.output_dir, "model." + dna_sequence + ".pdb"))
        except Exception as run_err:
            sys.stderr.write(f"\n[FATAL TERMINATION] Sequence {dna_sequence} failed processing: {str(run_err)}\n")
            raise run_err

if __name__ == "__main__":
    main()
