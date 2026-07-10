import os, sys, re
import Bio.PDB
import configparser
import optparse
import subprocess

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
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.structure.contacts import contacts,interface
from ModCRElib.structure.dna import x3dna

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    Parse command-line options for DNA-based structural superimposition.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for static/mobile inputs and output.

    """

    parser = optparse.OptionParser("Usage: dna_superimposer.py -a pdb_static -b pdb_mobile [--dummy=dummy_dir -i interface_alignment -o output_dir -v]")

    parser.add_option("-a", action="store", type="string", dest="pdb_static", help="PDB file (static; e.g. from model_protein.py)", metavar="{finelame}")
    parser.add_option("-b", action="store", type="string", dest="pdb_mobile", help="PDB file (mobile; e.g. from model_protein.py)", metavar="{finelame}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="interface_alignment", help="Interface alignment file (from interfaces_aligner.py)", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    
    (options, args) = parser.parse_args()

    if options.pdb_file_a is None or options.pdb_file_b is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def dna_superimposition_2(pdb_file_A, pdb_file_B, extension, dummy_dir="/tmp"):
    """
        Superimpose a mobile DNA model onto a static DNA model by sequence offset.
    
        Args:
            pdb_file_A (str): Static/reference PDB file.
            pdb_file_B (str): Mobile PDB file to transform.
            extension (int): Basepair offset applied to static alignment positions.
            dummy_dir (str, optional): Temporary directory root.
    
        Returns:
            None.
    
    Args:
        pdb_file_A (Any): Value used by this routine.
        pdb_file_B (Any): Value used by this routine.
        extension (Any): Value used by this routine.
        dummy_dir (Any): Directory path used by this operation."""

    # Initialize #
    dummy_file = os.path.join(dummy_dir, "%s.pdb" % os.getpid())
    alignment_A = set()
    alignment_B = set()
    backbone_atoms = set(["P", "O1P", "O2P", "O3P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 
                               "OP1", "OP2", "OP3", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*"])
    pdb_obj_A = PDB(pdb_file_A)
    pdb_obj_B = PDB(pdb_file_B)
    a_chains=[]
    for c in pdb_obj_A.chains:
     identifier = c.chain
     if identifier not in a_chains and c.chaintype=="N": 
        a_chains.append(identifier)
    b_chains=[]
    for c in pdb_obj_B.chains:
     identifier = c.chain
     if identifier not in a_chains and c.chaintype=="N": 
        b_chains.append(identifier)
    sequence_A  = pdb_obj_A.get_chain_by_id(a_chains[0]).gapped_nucleotide_sequence()
    sequence_B  = pdb_obj_B.get_chain_by_id(b_chains[0]).gapped_nucleotide_sequence()
    for i in range(len(sequence_B)):
        alignment_B.add((b_chains[0],i+1))
        alignment_A.add((a_chains[0],i+extension))
    # Start the parser #
    pdb_parser = Bio.PDB.PDBParser(QUIET=True) 
    # Load the PDB files #
    static_pdb = pdb_parser.get_structure("static", pdb_file_A)
    mobile_pdb = pdb_parser.get_structure("mobile", pdb_file_B)
    # Use the first model in each PDB #
    static_model = static_pdb[0]
    mobile_model = mobile_pdb[0]
    # Get the residues to be aligned #
    static_residues = []
    mobile_residues = []
    # For each chain... #
    for chain in static_model:
        # For each residue... #
        for residue in chain:
            # Initialize #
            residue_num = residue.get_id()
            # If residue and chain in alignment... #
            if (chain.get_id(), residue_num[1]) in alignment_A:
                static_residues.append(residue)
    # For each chain... #
    for chain in mobile_model:
        # For each residue... #
        for residue in chain:
            # Initialize #
            residue_num = residue.get_id()
            # If residue and chain in alignment... #
            if (chain.get_id(), residue_num[1]) in alignment_B:
                mobile_residues.append(residue)
    # Get the atoms to be aligned #
    static_atoms = []
    mobile_atoms = []
    # For each residue... #
    for i in range(len(static_residues)):
        static_residue = static_residues[i]
        mobile_residue = mobile_residues[i]
        # For each backbone atom... #
        for backbone_atom in backbone_atoms:
            # Initialize #
            atom_in_static_residue = False
            atom_in_mobile_residue = False
            # For each atom... #
            for atom in static_residue:
                if atom.name == backbone_atom:
                    atom_in_static_residue = True
                    break
            # For each atom... #
            for atom in mobile_residue:
                if atom.name == backbone_atom:
                    atom_in_mobile_residue = True
                    break
            # If both residues have backbone atom... #
            if atom_in_static_residue and atom_in_mobile_residue:
                static_atoms.append(static_residue[backbone_atom])
                mobile_atoms.append(mobile_residue[backbone_atom])
    # Initiate the superimposer #
    superimposer = Bio.PDB.Superimposer() 
    superimposer.set_atoms(static_atoms, mobile_atoms)
    superimposer.apply(mobile_model.get_atoms()) 
    # Save the superposed structure #
    io = Bio.PDB.PDBIO()
    io.set_structure(mobile_pdb)
    io.save(dummy_file)

    # Get PDB object #
    pdb_obj = PDB(dummy_file)

    


def dna_superimposition(pdb_file_A, pdb_file_B, alignment, dummy_dir="/tmp"):
    """
        Superimpose two DNA-containing PDBs using explicit residue-pair alignment.
    
        Args:
            pdb_file_A (str): Static/reference PDB file.
            pdb_file_B (str): Mobile PDB file to transform.
            alignment (list): Residue alignment pairs as parsed chain/position tuples.
            dummy_dir (str, optional): Temporary directory root.
    
        Returns:
            PDB: Superimposed mobile structure as PDB object.
    
    Args:
        pdb_file_A (Any): Value used by this routine.
        pdb_file_B (Any): Value used by this routine.
        alignment (Any): Sequence/alignment content used in this step.
        dummy_dir (Any): Directory path used by this operation."""

    # Initialize #
    dummy_file = os.path.join(dummy_dir, "%s.pdb" % os.getpid())
    alignment_A = set()
    alignment_B = set()
    backbone_atoms = set(["P", "O1P", "O2P", "O3P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 
                               "OP1", "OP2", "OP3", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*"])

    # For each alignment position... #
    for i in alignment:
        alignment_A.add((i[0][0], int(i[0][1])))
        alignment_B.add((i[1][0], int(i[1][1])))
    # Start the parser #
    pdb_parser = Bio.PDB.PDBParser(QUIET=True) 
    # Load the PDB files #
    static_pdb = pdb_parser.get_structure("static", pdb_file_A)
    mobile_pdb = pdb_parser.get_structure("mobile", pdb_file_B)
    # Use the first model in each PDB #
    static_model = static_pdb[0]
    mobile_model = mobile_pdb[0]
    # Get the residues to be aligned #
    static_residues = []
    mobile_residues = []
    # For each chain... #
    for chain in static_model:
        # For each residue... #
        for residue in chain:
            # Initialize #
            residue_num = residue.get_id()
            # If residue and chain in alignment... #
            if (chain.get_id(), residue_num[1]) in alignment_A:
                static_residues.append(residue)
    # For each chain... #
    for chain in mobile_model:
        # For each residue... #
        for residue in chain:
            # Initialize #
            residue_num = residue.get_id()
            # If residue and chain in alignment... #
            if (chain.get_id(), residue_num[1]) in alignment_B:
                mobile_residues.append(residue)
    # Get the atoms to be aligned #
    static_atoms = []
    mobile_atoms = []
    # For each residue... #
    for i in range(len(static_residues)):
        static_residue = static_residues[i]
        mobile_residue = mobile_residues[i]
        # For each backbone atom... #
        for backbone_atom in backbone_atoms:
            # Initialize #
            atom_in_static_residue = False
            atom_in_mobile_residue = False
            # For each atom... #
            for atom in static_residue:
                if atom.name == backbone_atom:
                    atom_in_static_residue = True
                    break
            # For each atom... #
            for atom in mobile_residue:
                if atom.name == backbone_atom:
                    atom_in_mobile_residue = True
                    break
            # If both residues have backbone atom... #
            if atom_in_static_residue and atom_in_mobile_residue:
                static_atoms.append(static_residue[backbone_atom])
                mobile_atoms.append(mobile_residue[backbone_atom])
    # Initiate the superimposer #
    superimposer = Bio.PDB.Superimposer() 
    superimposer.set_atoms(static_atoms, mobile_atoms)
    superimposer.apply(mobile_model.get_atoms()) 
    # Save the superposed structure #
    io = Bio.PDB.PDBIO()
    io.set_structure(mobile_pdb)
    io.save(dummy_file)

    # Get PDB object #
    pdb_obj = PDB(dummy_file)

    # Remove files #
    os.remove(dummy_file)

    return pdb_obj

def main():
    """
    Run the command-line DNA superimposition workflow.

    Args:
        None.

    Returns:
        None. Superimposed structure is written to the output directory.
    """

    # Arguments & Options #
    options = parse_options()

    # Create output directory #
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # If interface alignment... #
    if os.path.exists(os.path.abspath(options.interface_alignment)):
        # Get PDB objects #
        pdb_obj_a = PDB(os.path.abspath(options.pdb_file_a))
        pdb_obj_b = PDB(os.path.abspath(options.pdb_file_b))

        # Get X3DNA objects #
        x3dna_obj_a = x3dna.get_x3dna_obj(os.path.abspath(options.pdb_file_a), os.path.abspath(options.dummy_dir))
        x3dna_obj_b = x3dna.get_x3dna_obj(os.path.abspath(options.pdb_file_b), os.path.abspath(options.dummy_dir))

        # Get contacts objects #
        contacts_obj_a = contacts.get_contacts_obj(pdb_obj_a, x3dna_obj_a, "pdi", options.distance_type, os.path.abspath(options.dummy_dir))
        contacts_obj_b = contacts.get_contacts_obj(pdb_obj_b, x3dna_obj_b, "pdi", options.distance_type, os.path.abspath(options.dummy_dir))

        # Get interface objects #
        interface_obj_a = interface.get_interface_obj(pdb_obj_a, x3dna_obj_a, contacts_obj_a, options.max_distance, os.path.abspath(options.dummy_dir))
        interface_obj_b = interface.get_interface_obj(pdb_obj_b, x3dna_obj_b, contacts_obj_b, options.max_distance, os.path.abspath(options.dummy_dir))

        # Get best interface alignment #
        alignment = align_interfaces(interface_obj_a, interface_obj_b, x3dna_obj_a, x3dna_obj_b, options.reverse)
    # Else... #
    else:
        # Initialize #
        alignment = []
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.interface_alignment)):
            if line.startswith("#"): continue
            line = line.split(";")
            alignment.append([(line[0], line[1]), (line[2], line[3])])

    # Get superimposed DNA object #
    pdb_obj = dna_superimposition(os.path.abspath(options.pdb_file_a), os.path.abspath(options.pdb_file_b), alignment, os.path.abspath(options.dummy_dir))

    # Initialize #
    superimposed_file = os.path.abspath(os.path.join(options.output_dir, "superimposed.pdb"))
    # Write PDB file #
    pdb_obj.write(output_file, force=True)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
