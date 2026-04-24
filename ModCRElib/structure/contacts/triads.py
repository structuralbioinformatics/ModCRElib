import os, sys, re
import configparser
import copy
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
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, aminoacids_surface, nitrogenous_bases, dna_complementary
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.contacts import contacts
from ModCRElib.structure.protein import dssp

#-------------#
# Classes     #
#-------------#

class Triads(object):
    """
    Store and manipulate protein-DNA triad descriptors.

    Triads are generated from SBILib-backed contacts (via `contacts.py`) plus
    DSSP/X3DNA annotations, then serialized as semicolon-delimited rows.

    Object features:
        - Stores triad records in `_triads` as `Triad` objects.
        - Supports loading/saving canonical triad text format.
        - Provides interface-distance filtering (`get_interface_triads`).
        - Provides overlap metrics across triad sets
          (`get_common_triads`, `get_percentage_common_triads`).
    """

    def __init__(self, file_name=None):
        """
                Initialize a triads container.
        
                Args:
                    file_name (str, optional): Path to a serialized triads file.
        
                Returns:
                    None.
        
        Args:
            file_name (Any): Value used by this routine."""
        self._file = file_name
        self._triads = []
        # Initialize #
        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        """
        Parse triads from a serialized text file.

        Expected format:
            - Optional header/comment lines starting with `#`.
            - One triad per row with five `;`-separated fields:
              `A_environment;B_environment;distance;residue_A;residue_B`.
            - Rows correspond to outputs produced by `Triad.return_as_string()`
              after building triads from SBILib-derived contacts.
        """
        for line in functions.parse_file(self._file):
            if line.startswith("#"): continue
            A_environment, B_environment, distance, residue_A, residue_B = line.split(";")
            self.add_triad(Triad(A_environment, B_environment, float(distance), residue_A, residue_B))

    def add_triad(self, triad_obj):
        """Add one triad object to the container.
        
        Args:
            triad_obj (Any): Triad object/data used for interaction modeling."""
        self._triads.append(triad_obj)

    def get_triads(self):
        """Return a copy of stored triads."""
        return copy.copy(self._triads)

    def get_interface_triads(self,d=None):
        """
                Return triads within the interface-distance threshold.
        
                Args:
                    d (float, optional): Distance threshold; inferred if None.
        
                Returns:
                    tuple: `(d, interface_triads)` with selected triads.
        
        Args:
            d (Any): Value used by this routine."""
        interface = []
        if d is None:
          d = float(config.get("Parameters", "interface_distance"))
          c = float(config.get("Parameters", "interface_minimum"))
          found = False
          while (d<=float(config.get("Parameters", "max_contact_distance")) and not found ):
             for triad in self.get_triads():
               if triad.get_contact_distance() <d:
                 interface.append(triad)
             if len(interface)>=c:
              found = True
             else:
              d += 5
        else:
          for triad in self.get_triads():
             if triad.get_contact_distance() <d:
                 interface.append(triad)
        return (d,interface)

    def get_common_triads(self, triads_obj):
        """
                Count common triads between two triad sets.
        
                Args:
                    triads_obj (Triads): Triads object to compare with.
        
                Returns:
                    int: Number of common triad environment keys.
        
        Args:
            triads_obj (Any): Triad object/data used for interaction modeling."""
        # INitialize #
        d = float(config.get("Parameters", "interface_distance"))
        c = float(config.get("Parameters", "interface_minimum"))
        found = False
        while (d<=float(config.get("Parameters", "max_contact_distance")) and not found ):
          a = {}
          b = {}
          for triad in self.get_triads():
            a.setdefault((triad._A_environment, triad._B_environment), 0)
            if triad.get_contact_distance() <d:
               a[(triad._A_environment, triad._B_environment)] += 1
          for triad in triads_obj.get_triads():
            b.setdefault((triad._A_environment, triad._B_environment), 0)
            if triad.get_contact_distance() <d:
               b[(triad._A_environment, triad._B_environment)] += 1
          if sum(a[i] for i in a.keys()) >=c and  sum(b[i] for i in b.keys()) >=c :
             found = True
          else:
             d += 5
        return sum(min(a[i], b[i]) for i in a if i in b)


    def get_percentage_common_triads(self, triads_obj):
        """
                Compute overlap ratio of common triads between two sets.
        
                Args:
                    triads_obj (Triads): Triads object to compare with.
        
                Returns:
                    float: Fraction of common triads over smaller interface set.
        
        Args:
            triads_obj (Any): Triad object/data used for interaction modeling."""
        # INitialize #
        d = float(config.get("Parameters", "interface_distance"))
        c = float(config.get("Parameters", "interface_minimum"))
        found = False
        while (d<=float(config.get("Parameters", "max_contact_distance")) and not found ):
          a = {}
          b = {}
          for triad in self.get_triads():
            a.setdefault((triad._A_environment, triad._B_environment), 0)
            if triad.get_contact_distance() <d:
               a[(triad._A_environment, triad._B_environment)] += 1
          for triad in triads_obj.get_triads():
            b.setdefault((triad._A_environment, triad._B_environment), 0)
            if triad.get_contact_distance() <d:
               b[(triad._A_environment, triad._B_environment)] += 1
          if sum(a[i] for i in a.keys()) >=c and  sum(b[i] for i in b.keys()) >=c :
             found = True
          else:
             d += 5

        return float(sum(min(a[i], b[i]) for i in a if i in b))/float(min(len(self.get_interface_triads(d)[1]),len(triads_obj.get_interface_triads(d)[1])))



    def write(self, file_name):
        """
                Write triads to the canonical serialized text format.
        
                Output format:
                    - Header:
                      `#aminoacid-hydrophobicity-degree_of_exposure-secondary_structure;`
                      `dinucleotide-nucleotide_types-dna_strand-dna_groove-dna_chemical_group;`
                      `distance;residue_A;residue_B`
                    - One semicolon-delimited triad row per line.
                    - Rows preserve descriptors derived from SBILib contact geometry.
        
                Args:
                    file_name (str): Output file path.
        
                Returns:
                    None.
        
        Args:
            file_name (Any): Value used by this routine."""
        # Initialize #
        if os.path.exists(file_name): os.remove(file_name)

        functions.write(file_name, "#aminoacid-hydrophobicity-degree_of_exposure-secondary_structure;dinucleotide-nucleotide_types-dna_strand-dna_groove-dna_chemical_group;distance;residue_A;residue_B")
        for triad_obj in self.get_triads():
            functions.write(file_name, triad_obj.return_as_string())

class Triad(object):
    """
    Represent one protein-DNA interaction triad descriptor.

    A triad stores:
        - `A_environment`: amino-acid environment descriptor
          (`AminoAcidEnvironment` or equivalent serialized text),
        - `B_environment`: dinucleotide environment descriptor
          (`DinucleotideEnvironment` or equivalent serialized text),
        - `distance`: contact distance between the protein residue and DNA
          context used to build the triad,
        - `residue_A`: protein residue identifier string,
        - `residue_B`: DNA residue identifier string(s), typically a
          dinucleotide/basepair context.

    Here, a dinucleotide means two contiguous nucleotide residues in
    sequence (adjacent positions on one DNA strand).

    Canonical serialized row format:
        `{A_environment};{B_environment};{distance};{residue_A};{residue_B}`
    """

    def __init__(self, A_environment, B_environment, distance, residue_A, residue_B):
        """Initialize one triad descriptor.
        
        Args:
            A_environment (Any): Value used by this routine.
            B_environment (Any): Value used by this routine.
            distance (Any): Value used by this routine.
            residue_A (Any): Value used by this routine.
            residue_B (Any): Value used by this routine."""
        self._A_environment = A_environment
        self._B_environment = B_environment
        self._distance = distance
        self._residue_A = residue_A
        self._residue_B = residue_B

    def get_triad(self):
        """Return triad environments paired with residue identifiers."""
        return copy.copy(((self._residue_A, self._A_environment), (self._residue_B, self._B_environment)))

    def get_contact_distance(self):
        """Return triad contact distance."""
        return copy.copy(self._distance)

    def get_residues(self):
        """Return residue identifiers involved in the triad."""
        return copy.copy(self._residue_A, self._residue_B)

    def return_as_string(self):
        """
        Serialize triad as a semicolon-delimited text row.

        Format:
            `{A_environment};{B_environment};{distance};{residue_A};{residue_B}`

        Where:
            - `A_environment` is serialized by `AminoAcidEnvironment.return_as_string()`.
            - `B_environment` is serialized by `DinucleotideEnvironment.return_as_string()`.
            - `distance` is the triad contact distance.
            - `residue_A` and `residue_B` are residue identifier strings.

        Returns:
            str: Serialized triad row.
        """
        try:
            return "%s;%s;%s;%s;%s" % (self._A_environment.return_as_string(), self._B_environment.return_as_string(), self._distance, self._residue_A, self._residue_B)
        except:
            return "%s;%s;%s;%s;%s" % (self._A_environment, self._B_environment, self._distance, self._residue_A, self._residue_B)


class AminoAcidEnvironment(object):
    """
    Describe amino-acid environment features used in triads.
    """
    
    def __init__(self, aminoacid_single_letter, accessible_surface_area, secondary_structure):
        """Initialize amino-acid environment descriptor.
        
        Args:
            aminoacid_single_letter (Any): Value used by this routine.
            accessible_surface_area (Any): Value used by this routine.
            secondary_structure (Any): Value used by this routine."""
        self._aminoacid_single_letter = aminoacid_single_letter
        self._accessible_surface_area = accessible_surface_area
        self._secondary_structure = secondary_structure

    def get_hydrophobicity(self):
        """Return amino-acid polarity class (`P`/`N`)."""
        if self._aminoacid_single_letter in aminoacids_polarity_boolean:
            if aminoacids_polarity_boolean[self._aminoacid_single_letter]:
                return "P"
            return "N"

        return None

    def get_degree_of_exposure(self):
        """Return exposure class (`B`/`E`) from relative accessibility."""
        min_exposure = float(config.get("Parameters", "min_exposure"))

        if self._aminoacid_single_letter in aminoacids_surface and self._accessible_surface_area != None:
            if self._accessible_surface_area / aminoacids_surface[self._aminoacid_single_letter] <= min_exposure:
                return "B"
            return "E"

        return None

    def get_secondary_structure(self):
        """Return normalized secondary-structure code."""
        allowed_secondary_structures = config.get("Parameters", "allowed_secondary_structures")

        if self._secondary_structure != None:
            if self._secondary_structure not in allowed_secondary_structures:
                return "C" # residue is in a "coil" region
            return self._secondary_structure

        return None

    def return_as_string(self):
        """Serialize amino-acid environment as text."""
        return "%s-%s-%s-%s" % (aminoacids1to3[self._aminoacid_single_letter], self.get_hydrophobicity(), self.get_degree_of_exposure(), self.get_secondary_structure())

class DinucleotideEnvironment(object):
    """
    Describe dinucleotide environment features used in triads.

    A dinucleotide is defined here as two contiguous nucleotide residues in
    sequence (adjacent positions on the same strand).
    """

    def __init__(self, dinucleotide_sequence, dna_strand, dna_groove, dna_chemical_group):
        """Initialize dinucleotide environment descriptor.
        
        Args:
            dinucleotide_sequence (Any): Sequence/alignment content used in this step.
            dna_strand (Any): DNA identifier, sequence, or DNA-related data.
            dna_groove (Any): DNA identifier, sequence, or DNA-related data.
            dna_chemical_group (Any): DNA identifier, sequence, or DNA-related data."""
        self._dinucleotide = dinucleotide_sequence
        self._dna_strand = dna_strand
        self._dna_groove = dna_groove
        self._dna_chemical_groove = dna_chemical_group

    def get_nitrogenous_bases(self):
        """
        Return purine/pyrimidine pattern for the dinucleotide.

        Args:
            None.

        Returns:
            str: Nitrogenous-base classes (`U`/`Y`) per nucleotide.
        """
        return "".join(nitrogenous_bases[nucleotide] for nucleotide in self._dinucleotide)

    def return_as_string(self):
        """Serialize dinucleotide environment as text."""
        return "%s-%s-%s-%s-%s" % (self._dinucleotide, self.get_nitrogenous_bases(), self._dna_strand, self._dna_groove, self._dna_chemical_groove)

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for triad extraction.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for triad generation.
    """

    parser = optparse.OptionParser("python triads.py -i INPUT_FILE [-c --dummy DUMMY_DIR -o OUTPUT_FILE -d DISTANCE_TYPE]")

    parser.add_option("-c", default=False, action="store_true", dest="complementary", help="Complementary triad (default = False)", metavar="{boolean}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")
    parser.add_option("-d", default="dinucleotides", action="store", type="string", dest="distance_type", help="Distance type (i.e. \"basepairs\", \"dinucleotides\" or \"mindist\"; default = dinucleotides)", metavar="{string}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
 

    (options, args) = parser.parse_args()
    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if options.distance_type != "basepairs" and options.distance_type != "dinucleotides" and options.distance_type != "mindist":
        parser.error("incorrect distance type: accepted values are \"basepairs\", \"dinucleotides\" or \"mindist\"")
    return options

def get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj, complementary=False, dummy_dir="/tmp"):
    """
        Translate protein-DNA contacts into triad descriptors.
    
        Args:
            pdb_obj (PDB): Input PDB structure.
            dssp_obj (DSSP): DSSP annotations.
            x3dna_obj (X3DNA): DNA annotation object.
            contacts_obj (Contacts): Protein-DNA contacts object, typically created
                by `contacts.get_contacts_obj()` from SBILib structure atoms.
            complementary (bool, optional): Add complementary triads.
            dummy_dir (str, optional): Reserved temporary directory argument.
    
        Returns:
            Triads: Computed triads object.
    
    Args:
        pdb_obj (Any): PDB object or structure file input.
        dssp_obj (Any): Value used by this routine.
        x3dna_obj (Any): DNA identifier, sequence, or DNA-related data.
        contacts_obj (Any): Contact object/data used by this routine.
        complementary (Any): Value used by this routine.
        dummy_dir (Any): Directory path used by this operation."""

    # Initialize #
    done = set() # Removes redundant triads
    aminoacid_environments = {}
    aminoacid_dna_grooves = {}
    triads_obj = Triads()

    # For each contact #
    for contact_obj in contacts_obj.get_contacts():
        if (contact_obj._A, tuple(contact_obj._B)) in done: continue
        if not pdb_obj.chain_exists(contact_obj._A_chain): continue
        if not pdb_obj.get_chain_by_id(contact_obj._A_chain).residue_exists(str(contact_obj._A_residue_obj.number)): continue
        # Initialize #
        proceed = False
        # For each dinucleotide... #
        for i in range(len(contact_obj._B)):
            proceed = True
            if not pdb_obj.chain_exists(contact_obj._B_chain[i]):
                proceed = False
                break
            if not pdb_obj.get_chain_by_id(contact_obj._B_chain[i]).residue_exists(str(contact_obj._B_residue_obj[i].number)):
                proceed = False
                break
        # This is a valid contact: proceed... #
        if proceed:
            # Get amino acid object #
            aminoacid_obj = pdb_obj.get_chain_by_id(contact_obj._A_chain).get_residue_by_identifier(str(contact_obj._A_residue_obj.number))
            # Get amino acid environment #
            if (contact_obj._A_chain, contact_obj._A_residue_obj.number) not in aminoacid_environments: aminoacid_environments[(contact_obj._A_chain, contact_obj._A_residue_obj.number)] = get_aminoacid_environment_obj(aminoacid_obj, contact_obj._A_chain, dssp_obj)
            aminoacid_environment = aminoacid_environments[(contact_obj._A_chain, contact_obj._A_residue_obj.number)]
            # Get amino acid DNA groove #
            if (contact_obj._A_chain, contact_obj._A_residue_obj.number) not in aminoacid_dna_grooves: aminoacid_dna_grooves[(contact_obj._A_chain, contact_obj._A_residue_obj.number)] = get_aminoacid_dna_groove(contact_obj._A, pdb_obj, x3dna_obj)
            dna_groove = aminoacid_dna_grooves[(contact_obj._A_chain, contact_obj._A_residue_obj.number)]
            # Get dinucleotide sequence #
            dinucleotide_sequence = contact_obj._B_residue_obj[0].single_letter + contact_obj._B_residue_obj[2].single_letter
            # Get DNA strand #
            dna_strand = get_closest_strand_to_aminoacid_cb(contact_obj._A, contact_obj._B)
            # Get DNA chemical group #
            dinucleotide_residues = []
            for i in range(len(contact_obj._B)):
                dinucleotide_residues.append(pdb_obj.get_chain_by_id(contact_obj._B_chain[i]).get_residue_by_identifier(str(contact_obj._B_residue_obj[i].number)))
            dna_chemical_group = get_closest_chemical_group(aminoacid_obj, dinucleotide_residues)
            # Get dinucleotide environment #
            dinucleotide_environment = DinucleotideEnvironment(dinucleotide_sequence, dna_strand, dna_groove, dna_chemical_group)
            # Add triad #
            triads_obj.add_triad(Triad(aminoacid_environment, dinucleotide_environment, contact_obj._distance, "%s-%s" % (contact_obj._A_chain, contact_obj._A_residue_obj.number), "%s-%s,%s-%s,%s-%s,%s-%s" % (contact_obj._B_chain[0], contact_obj._B_residue_obj[0].number, contact_obj._B_chain[1], contact_obj._B_residue_obj[1].number, contact_obj._B_chain[2], contact_obj._B_residue_obj[2].number, contact_obj._B_chain[3], contact_obj._B_residue_obj[3].number)))
            # Get complementary dinucleotide environment #
            if complementary:
                dinucleotide_environment = DinucleotideEnvironment(get_complementary_dna_sequence(dinucleotide_sequence), get_complementary_dna_strand(dna_strand), dna_groove, dna_chemical_group)
                # Add triad #
                triads_obj.add_triad(Triad(aminoacid_environment, dinucleotide_environment, contact_obj._distance, "%s-%s" % (contact_obj._A_chain, contact_obj._A_residue_obj.number), "%s-%s,%s-%s,%s-%s,%s-%s" % (contact_obj._B_chain[0], contact_obj._B_residue_obj[0].number, contact_obj._B_chain[1], contact_obj._B_residue_obj[1].number, contact_obj._B_chain[2], contact_obj._B_residue_obj[2].number, contact_obj._B_chain[3], contact_obj._B_residue_obj[3].number)))

    return triads_obj

def get_aminoacid_environment_obj(aminoacid_obj, pdb_chain, dssp_obj):
    """
        Build amino-acid environment descriptor from structural annotations.
    
        Args:
            aminoacid_obj (ResidueOfAminoAcid): Amino-acid residue.
            pdb_chain (str): PDB chain identifier.
            dssp_obj (DSSP): DSSP annotation object.
    
        Returns:
            AminoAcidEnvironment: Computed amino-acid environment object.
    
    Args:
        aminoacid_obj (Any): Value used by this routine.
        pdb_chain (Any): PDB structure identifier and chain selector.
        dssp_obj (Any): Value used by this routine."""

    accessible_surface_area = dssp_obj.get_accessible_surface_area(pdb_chain, aminoacid_obj.number)
    secondary_structure = dssp_obj.get_secondary_structure(pdb_chain, aminoacid_obj.number)

    return AminoAcidEnvironment(aminoacid_obj.single_letter, accessible_surface_area, secondary_structure)

def get_aminoacid_dna_groove(aminoacid_atom_obj, pdb_obj, x3dna_obj):
    """
        This function returns whether an {AminoAcid} is located in the major
        or the minor groove according to the definition of grooves from
        El Hassan & Calladine, 1998:
        
        Select the closest phosphate from each strand to the cb {Atom} of the
        {AminoAcid}; let this be at position "p" for strand "s" and at position
        "P" for strand "S".
                                  cb.
                                 /   \ 
           strand s . ,-"-.   ,-p-. ,-P-.   ,-"-. ,-"-.   ,-"-. ,-"-.   ,
                     X | | \ / | | X | | \ / | | X | | \ / | | X | | \ /
                    / \| | |X| | |/ \| | |X| | |/ \| | |X| | |/ \| | |X|
           strand S    `-!-' `-!-'   `-!-' `-!-'   `-!-' `-!-'   `-.-' `-
    
        If p < P, calculate the distances between the phosphate {Atom} at "p" in
        "s" and the phosphate {Atom}s from "p+2" to "p+5" in "S". The minimum
        defines the groove width where the amino acid is located (direct groove).
        Also, calculate the distances between the phosphate atom at "p" in "s"
        and the phosphate {Atom}s from "p-2" to "p-5" in "S". The minimum defines
        the groove width where the amino acid is not located (indirect groove).
    
        If direct groove width > indirect groove width, the amino acid is located
        in the major groove ("A"); otherwise, it is located in the minor groove
        ("I").
    
        * If in p > P, instead of calculating the distances from "p+2" to "p+5"
        in "S", use the distances from "p-2" to "p-5" in "S" for the direct groove
        width and viceversa for the indirect groove width, and apply the same
        criterion to select the DNA groove where the amino acid is located.
    
        ** If p = P, the amino acid location is dubious.
            Determine whether an amino-acid atom lies in major or minor DNA groove.
    
        The groove assignment follows the geometric criterion from
        El Hassan and Calladine (1998), using phosphates around the two closest
        basepair positions (one per strand):
            1. Find the closest phosphate index in each strand (`p` and `P`).
            2. For each strand, measure phosphate-phosphate distances to the
               opposite strand in the windows `+1..+5` and `-1..-5`.
            3. Depending on whether `p < P` or `p > P`, one direction is treated
               as the direct groove and the other as the indirect groove.
            4. If `min(direct) > min(indirect)`, classify as major groove (`"A"`),
               otherwise classify as minor groove (`"I"`).
    
        The function returns `None` when the geometry is ambiguous (for example,
        missing phosphates, unresolved basepairs, or insufficient distances).
    
        Args:
            aminoacid_atom_obj (AtomOfAminoAcid): Query amino-acid atom.
            pdb_obj (PDB): Input PDB structure.
            x3dna_obj (X3DNA): DNA annotation object.
    
        Returns:
            str or None: `"A"` for major groove, `"I"` for minor groove, or None.
    
    Args:
        aminoacid_atom_obj (Any): Value used by this routine.
        pdb_obj (Any): PDB object or structure file input.
        x3dna_obj (Any): DNA identifier, sequence, or DNA-related data."""

    # Initialize #
    closest_nucleotides = get_closest_phosphates_to_aminoacid_cb(aminoacid_atom_obj, pdb_obj, x3dna_obj)

    if closest_nucleotides[0] != closest_nucleotides[1]:
        # Initialize #
        direct = set()
        indirect = set()
        # For each strand... #
        for strand in [0, 1]:
            # Get complementary strand #
            if strand == 0:
                complementary_strand = 1
            else:
                complementary_strand = 0
            # Get phosphate #
            p = closest_nucleotides[strand]
            # Get basepair #
            basepair_p = x3dna_obj.get_basepair(p)
            # Get nucleotide p object #
            if not pdb_obj.chain_exists(basepair_p[strand][0]):
                continue
            if not pdb_obj.get_chain_by_id(basepair_p[strand][0]).residue_exists(str(basepair_p[strand][1])):
                continue
            nucleotide_p_obj = pdb_obj.get_chain_by_id(basepair_p[strand][0]).get_residue_by_identifier(str(basepair_p[strand][1]))
            # Get phosphate p atom object #
            dna_atom_p_obj = contacts.get_nucleotide_p_or_bb(nucleotide_p_obj)
            # For upstream/downstream i, j #
            for multiplier in [1, -1]:
                # The minimum is generally between +/-1 and +/-5
                for Pn in range(p + (1 * multiplier), p + ((5 + 1) * multiplier), multiplier):
                    if x3dna_obj.has_basepair(Pn):
                        # Get basepair #
                        basepair_Pn = x3dna_obj.get_basepair(Pn)
                        # Get nucleotide Pn object #
                        if not pdb_obj.chain_exists(basepair_Pn[complementary_strand][0]):
                            continue
                        if not pdb_obj.get_chain_by_id(basepair_Pn[complementary_strand][0]).residue_exists(str(basepair_Pn[complementary_strand][1])):
                            continue
                        nucleotide_Pn_obj = pdb_obj.get_chain_by_id(basepair_Pn[complementary_strand][0]).get_residue_by_identifier(str(basepair_Pn[complementary_strand][1]))
                        # Get phosphate Pn atom object #
                        dna_atom_Pn_obj = contacts.get_nucleotide_p_or_bb(nucleotide_Pn_obj)
                        # Get distance #
                        distance = dna_atom_p_obj.distance(dna_atom_Pn_obj)
                        # Skip if distance could not be calculated #
                        if distance == -1: continue
                        # Add distance #
                        if p < closest_nucleotides[complementary_strand]:
                            if multiplier == 1:
                                direct.add(distance)
                            if multiplier == -1:
                                indirect.add(distance)
                        if p > closest_nucleotides[complementary_strand]:
                            if multiplier == 1:
                                indirect.add(distance)
                            if multiplier == -1:
                                direct.add(distance)
        # Skip if direct/indirect grooves could not be calculated #
        if len(direct) == 0 or len(indirect) == 0: return None
        # Return major groove #
        if min(direct) > min(indirect):
            return "A"
        # Return minor groove #
        else:
            return "I"

    return None

def get_closest_phosphates_to_aminoacid_cb(aminoacid_atom_obj, pdb_obj, x3dna_obj):
    """
        Find nearest phosphate positions on forward and reverse strands.
    
        Args:
            aminoacid_atom_obj (AtomOfAminoAcid): Query amino-acid atom.
            pdb_obj (PDB): Input PDB structure.
            x3dna_obj (X3DNA): DNA annotation object.
    
        Returns:
            tuple: `(i, j)` closest basepair indices for each strand.
    
    Args:
        aminoacid_atom_obj (Any): Value used by this routine.
        pdb_obj (Any): PDB object or structure file input.
        x3dna_obj (Any): DNA identifier, sequence, or DNA-related data."""

    # Initialize #
    i = []
    j = []

    for basepair in x3dna_obj.get_basepairs():
        (pdb_chain_i, residue_num_i), (pdb_chain_j, residue_num_j) = x3dna_obj.get_basepair(basepair)
        if pdb_obj.chain_exists(pdb_chain_i) and pdb_obj.chain_exists(pdb_chain_j):
            if pdb_obj.get_chain_by_id(pdb_chain_i).residue_exists(str(residue_num_i)) and pdb_obj.get_chain_by_id(pdb_chain_j).residue_exists(str(residue_num_j)):
                dna_atom_obj = contacts.get_nucleotide_p_or_bb(pdb_obj.get_chain_by_id(pdb_chain_i).get_residue_by_identifier(str(residue_num_i)))
                distance = aminoacid_atom_obj.distance(dna_atom_obj)
                if distance != -1:
                    i.append((basepair, distance))
                dna_atom_obj = contacts.get_nucleotide_p_or_bb(pdb_obj.get_chain_by_id(pdb_chain_j).get_residue_by_identifier(str(residue_num_j)))
                distance = aminoacid_atom_obj.distance(dna_atom_obj)
                if distance != -1:
                    j.append((basepair, distance))

    i.sort(key=lambda x: x[1])
    j.sort(key=lambda x: x[1])

    return i[0][0], j[0][0]

def get_closest_strand_to_aminoacid_cb(aminoacid_atom_obj, dinucleotide_atoms):
    """
        Return closest DNA strand (`F`/`R`) to an amino-acid atom.
    
        Args:
            aminoacid_atom_obj (AtomOfAminoAcid): Query amino-acid atom.
            dinucleotide_atoms (list): Dinucleotide atom list.
    
        Returns:
            str or None: `F`, `R`, or None when unresolved.
    
    Args:
        aminoacid_atom_obj (Any): Value used by this routine.
        dinucleotide_atoms (Any): Value used by this routine."""

    # Initialize #
    distances = []

    for dna_atom_obj in dinucleotide_atoms:
        # Get distance #
        distance = aminoacid_atom_obj.distance(dna_atom_obj)
        if distance != -1:
            distances.append(distance)

    # If not enough distances... #
    if len(distances) != 4: return None

    # If min. distance index is even; it is a forward nucleotide #
    if distances.index(min(distances)) % 2 == 0:
            return "F"

    return "R"

def get_closest_chemical_group(aminoacid_obj, dinucleotide_residues):
    """
        Return closest DNA chemical group (`B`/`N`) to an amino acid.
    
        Args:
            aminoacid_obj (ResidueOfAminoAcid): Amino-acid residue.
            dinucleotide_residues (list): Dinucleotide residue list.
    
        Returns:
            str or None: `B`, `N`, or None when unresolved.
    
    Args:
        aminoacid_obj (Any): Value used by this routine.
        dinucleotide_residues (Any): Value used by this routine."""

    # Initialize #
    contacts = []

    for nucleotide_obj in dinucleotide_residues:
        contacts.append(aminoacid_obj.distance(nucleotide_obj, "min"))

    # For each contact... #
    for contact in sorted(contacts, key=lambda x: x[-1]):
        # Skip if not contact #
        if contact[1] == None: continue
        # If contact is through backbone... #
        if contact[1].is_backbone:
            return "B"
        return "N"

    return None

def get_complementary_dna_sequence(dna_sequence):
    """
        Return reverse-complement DNA sequence.
    
        Args:
            dna_sequence (str): Input DNA sequence.
    
        Returns:
            str: Reverse-complement sequence.
    
    Args:
        dna_sequence (Any): DNA identifier, sequence, or DNA-related data."""

    return "".join(dna_complementary[nucleotide] for nucleotide in dna_sequence[::-1])

def get_complementary_dna_strand(dna_strand):
    """
        Return complementary strand label.
    
        Args:
            dna_strand (str): Input strand label (`F` or `R`).
    
        Returns:
            str: Complementary strand label.
    
    Args:
        dna_strand (Any): DNA identifier, sequence, or DNA-related data."""

    if dna_strand == "F":
        return "R"

    return "F"

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    """
    Run the command-line triads extraction workflow.

    Args:
        None.

    Returns:
        None. Triads are written to output file or stdout.
    """

    # Arguments & Options #
    options = parse_options()

    # Get PDB object #
    pdb_obj = PDB(os.path.abspath(options.input_file))

    # Get X3DNA object #
    dssp_obj = dssp.get_dssp_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Get X3DNA object #
    x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Get contacts object #
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj, "pdi", options.distance_type, os.path.abspath(options.dummy_dir))

    # Get triads object #
    triads_obj = get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj, options.complementary, os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        triads_obj.write(os.path.abspath(options.output_file))
    else:
        sys.stdout.write("#aminoacid-hydrophobicity-degree_of_exposure-secondary_structure;dinucleotide-nucleotide_types-dna_strand-dna_groove-dna_chemical_group;distance\n")
        for triad_obj in triads_obj.get_triads():
            sys.stdout.write(triad_obj.return_as_string())
