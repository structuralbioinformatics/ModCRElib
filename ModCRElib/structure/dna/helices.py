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

# Imports my functions #
from ModCRElib.beans import functions

# Import jbonet's module #
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.structure.contacts import contacts
from ModCRElib.structure.dna import x3dna

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for protein-chain to DNA-helix assignment.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for input/output paths.

    """

    parser = optparse.OptionParser("python helices.py -i INPUT_FILE [--dummy DUMMY_DIR -o OUTPUT_DIR]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_protein_chains_dna_helices(pdb_obj, x3dna_obj, contacts_obj, output_dir="./", dummy_dir="/tmp"):
    """
        Assign the most-contacted DNA helix to each protein chain.
    
        In this script, a DNA helix is defined as a contiguous nucleotide polymer
        segment together with its complementary strand, forming a double-stranded
        DNA helix.
    
        For every protein chain in `pdb_obj`, the function counts contacts against
        helix-associated dinucleotides and writes the top supported helix to
        `<pdbid>_<chain>.txt` in `output_dir`.
    
        Args:
            pdb_obj (PDB): Input PDB structure.
            x3dna_obj (X3DNA): Parsed DNA basepair/helix annotations.
            contacts_obj (Contacts): Precomputed protein-DNA contacts.
            output_dir (str, optional): Destination directory for helix files.
            dummy_dir (str, optional): Reserved temporary directory argument.
    
        Returns:
            None.
    
    Args:
        pdb_obj (Any): PDB object or structure file input.
        x3dna_obj (Any): DNA identifier, sequence, or DNA-related data.
        contacts_obj (Any): Contact object/data used by this routine.
        output_dir (Any): Directory path used by this operation.
        dummy_dir (Any): Directory path used by this operation."""

    # Initialize #
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    min_contacts = int(config.get("Parameters", "min_contacts"))

    # For protein chain... #
    for protein_chain_obj in pdb_obj.proteins:
        # Initialize #
        contacts = {}
        helix_file = os.path.join(output_dir, pdb_obj.id.lower() + "_" + protein_chain_obj.chain + ".txt")
        # For each helix... #
        for helix in x3dna_obj.get_dna_helices():
            # Initialize #
            contacts.setdefault(helix, set())
            # Get helix dinucleotides #
            dinucleotides = x3dna_obj.get_helix_dinucleotides(helix)
            # For each contact... #
            for contact_obj in contacts_obj.get_contacts():
                # If protein chain... #
                if contact_obj._A_chain == protein_chain_obj.chain:
                    # For each dinucleotide... #
                    for dinucleotide in x3dna_obj.get_basepair_dinucleotides(x3dna_obj.get_residue_basepair(contact_obj._B_chain[0], contact_obj._B_residue_obj[0].number)):
                        # If dinucleotide in helix dinucleotides... #
                        if dinucleotide in dinucleotides:
                            contacts[helix].add(contact_obj)
                            break
        # For each helix... #
        for helix in sorted(contacts, key=lambda x: contacts[x], reverse=True):
            # Skip if not enough contacts #
            if len(contacts[helix]) < min_contacts: continue
            # Erase helix file if already exists #
            if os.path.exists(helix_file): os.remove(helix_file)
            # Create output #
            functions.write(helix_file, str(helix))
            break

def main():
    """
    Run the command-line DNA-helix assignment workflow.

    Workflow:
        1. Parse command-line options.
        2. Build PDB, X3DNA, and contact objects.
        3. Write one best helix assignment per protein chain.

    Args:
        None.

    Returns:
        None. Helix assignment files are written to the output directory.
    """

    # Arguments & Options #
    options = parse_options()

    # Get PDB object #
    pdb_obj = PDB(os.path.abspath(options.input_file))

    # Get X3DNA object #
    x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Get contacts object #
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj, "pdi", "dinucleotides", os.path.abspath(options.dummy_dir))

    # Get protein chains DNA helices #
    get_protein_chains_dna_helices(pdb_obj, x3dna_obj, contacts_obj, os.path.abspath(options.output_dir), os.path.abspath(options.dummy_dir))


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
