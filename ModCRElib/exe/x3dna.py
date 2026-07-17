import os, sys, re
import configparser
import copy
import optparse
import shutil
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

#-------------#
# Classes     #
#-------------#

class X3DNA(object):
    """
    Parse and store basepair, helix, and dinucleotide annotations from X3DNA.

    Object features:
        - Preserves raw X3DNA output lines in `_file_content`.
        - Stores forward-strand sequence in `_sequence`.
        - Tracks residue-to-basepair assignments in `_residues`.
        - Stores basepair definitions in `_basepairs` as paired strand residues.
        - Stores contiguous dinucleotide windows in `_dinucleotides`, plus
          reverse lookup in `_inverse_dinucleotides`.
        - Stores helix-to-basepair assignments in `_helix` and derived helix
          dinucleotide indices in `_helix_dinucleotides`.
        - Exposes query helpers for residue/basepair/dinucleotide/helix
          membership and retrieval, including sequence slicing by basepair range.
    """

    def __init__(self, file_name):
        """
                Initialize an X3DNA container from a `find_pair` output file.
        
                Args:
                    file_name (str): Path to the X3DNA output file.
        
                Returns:
                    None.
        
        Args:
            file_name (Any): Value used by this routine."""
        self._file = file_name
        self._file_content = []
        self._sequence = ""
        self._residues = {}
        self._basepairs = {}
        self._dinucleotides = {}
        self._inverse_dinucleotides = {}
        self._helix = {}
        self._helix_dinucleotides = {}
        # Initialize #
        self._parse_file()
        self._initialize_dinucleotides()
        self._initialize_helix_dinucleotides()

    def _parse_file(self):
        """Parse basepairs/helices and cache raw output lines."""
        for line in functions.parse_file(self._file):
            # Get base-pair residues #
            m = re.search("(\d+)\s+\S\s\S{4}>(\S):\D*(\d+)_:\[\S{3}\](\w)\S{5}\w\[\S{3}\]:\D*(\d+)_:(\S)<\S{4}", line)
            if m:
                basepair = int(m.group(1))
                fwd_pdb_chain = m.group(2)
                fwd_residue_num = int(m.group(3))
                self._sequence += m.group(4)
                rev_pdb_chain = m.group(6)
                rev_residue_num = int(m.group(5))
                self._residues[(fwd_pdb_chain, fwd_residue_num)] = basepair
                self._residues[(rev_pdb_chain, rev_residue_num)] = basepair
                self._basepairs[basepair] = ((fwd_pdb_chain, fwd_residue_num), (rev_pdb_chain, rev_residue_num))
            m = re.search("^#####\s+Helix\s+#(\d+)\s+\(\d+\):\s+(\d+)\s+\-\s+(\d+)", line)
            if m:
                helix = m.group(1)
                first = int(m.group(2))
                last = int(m.group(3))
                self._helix.setdefault(helix, list(range(first, last + 1)))
            # Add line to file content #
            if line != "":
               self._file_content.append(line)
            

    def _initialize_dinucleotides(self):
        """Initialize consecutive basepair windows as dinucleotides."""
        basepairs = [key for key in sorted(self._basepairs.keys())]
        while len(basepairs) > 1:
            i = basepairs.pop(0)
            if self.has_dinucleotide(i): continue
            self._dinucleotides.setdefault(i, (i, basepairs[0]))
            self._inverse_dinucleotides.setdefault((i, basepairs[0]),i)

    def _initialize_helix_dinucleotides(self):
        """Map helix identifiers to their constituent dinucleotide starts."""
        for helix in self._helix:
            self._helix_dinucleotides.setdefault(helix, [])
            basepairs = copy.copy(self.get_helix_basepairs(helix))

            while len(basepairs) > 1:
                i = basepairs.pop(0)
                self._helix_dinucleotides[helix].append(i)

    def has_residue(self, pdb_chain, residue_num):
        """Return whether a residue is mapped to a basepair.
        
        Args:
            pdb_chain (Any): PDB structure identifier and chain selector.
            residue_num (Any): Numeric identifier/index used by this routine."""
        return (pdb_chain, residue_num) in self._residues

    def has_basepair(self, basepair):
        """Return whether a basepair index exists.
        
        Args:
            basepair (Any): Value used by this routine."""
        return basepair in self._basepairs

    def has_dinucleotide(self, dinucleotide):
        """Return whether a dinucleotide index exists.
        
        Args:
            dinucleotide (Any): Value used by this routine."""
        return dinucleotide in self._dinucleotides

    def has_helix(self, helix):
        """Return whether a DNA helix identifier exists.
        
        Args:
            helix (Any): Value used by this routine."""
        return helix in self._helix

    def helix_has_residue(self, helix, pdb_chain, residue_num):
        """Return whether a residue belongs to a specific helix.
        
        Args:
            helix (Any): Value used by this routine.
            pdb_chain (Any): PDB structure identifier and chain selector.
            residue_num (Any): Numeric identifier/index used by this routine."""
        if self.has_residue(pdb_chain, residue_num) and self.has_helix(helix):
            return self._residues[(pdb_chain, residue_num)] in self.get_helix_basepairs(helix)

        return False

    def get_residue_basepair(self, pdb_chain, residue_num):
        """Return basepair index assigned to a residue.
        
        Args:
            pdb_chain (Any): PDB structure identifier and chain selector.
            residue_num (Any): Numeric identifier/index used by this routine."""
        if self.has_residue(pdb_chain, residue_num):
            return copy.copy(self._residues[(pdb_chain, residue_num)])

        return None

    def get_basepair(self, basepair):
        """Return residue tuples that define one basepair.
        
        Args:
            basepair (Any): Value used by this routine."""
        if self.has_basepair(basepair):
            return copy.copy(self._basepairs[basepair])

        return None

    def get_basepairs(self):
        """Return all basepair mappings."""
        return copy.copy(self._basepairs)

    def get_helix_basepairs(self, helix):
        """Return ordered basepairs assigned to one helix.
        
        Args:
            helix (Any): Value used by this routine."""
        if self.has_helix(helix):
            return copy.copy(self._helix[helix])

        return None

    def get_dinucleotide(self, dinucleotide):
        """Return the basepair pair represented by one dinucleotide index.
        
        Args:
            dinucleotide (Any): Value used by this routine."""
        if self.has_dinucleotide(dinucleotide):
            return copy.copy(self._dinucleotides[dinucleotide])

        return None
  
    def get_inverse_dinucleotide(self, basepair_1, basepair_2):
        """Return dinucleotide index for a basepair pair.
        
        Args:
            basepair_1 (Any): Value used by this routine.
            basepair_2 (Any): Value used by this routine."""
        if (basepair_1,basepair_2) in self._inverse_dinucleotides:
           return copy.copy(self._inverse_dinucleotides[(basepair_1,basepair_2)])
        return None
 
    def get_sequence(self):
        """Return forward-strand sequence parsed from X3DNA output."""
        return copy.copy(self._sequence)

    def get_dinucleotides(self):
        """Return all dinucleotide mappings."""
        return copy.copy(self._dinucleotides)

    def get_inverse_dinucleotides(self):
        """Return reverse mapping from basepair pairs to dinucleotide index."""
        return copy.copy(self._inverse_dinucleotides)


    def get_helix_dinucleotides(self, helix):
        """Return dinucleotide indices that belong to a helix.
        
        Args:
            helix (Any): Value used by this routine."""
        if self.has_helix(helix):
            return copy.copy(self._helix_dinucleotides[helix])

        return None

    def get_basepair_dinucleotides(self, basepair):
        """Return dinucleotide indices that include a given basepair.
        
        Args:
            basepair (Any): Value used by this routine."""
        dinucleotides = []

        for dinucleotide in self.get_dinucleotides():
            if basepair in self.get_dinucleotide(dinucleotide):
                dinucleotides.append(dinucleotide)
        
        return copy.copy(dinucleotides)

    def get_basepair_helix(self, basepair):
        """Return helix identifier that contains a basepair.
        
        Args:
            basepair (Any): Value used by this routine."""
        if self.has_basepair(basepair):
            for helix in self._helix:
                if basepair in self.get_helix_basepairs(helix):
                    return helix

        return None

    def get_dna_helices(self):
        """Return all parsed DNA helices."""
        return copy.copy(self._helix)

    def get_nucleotide_sequence(self, A, B):
        """Return sequence segment between basepair indices `A` and `B`.
        
        Args:
            A (Any): Value used by this routine.
            B (Any): Value used by this routine."""
        if self.has_basepair(A) and self.has_basepair(B):
            return self._sequence[A - 1:B]

        return None

    def write(self, file_name):
        """Write cached raw X3DNA output lines to file.
        
        Args:
            file_name (Any): Value used by this routine."""
        for line in self._file_content:
            functions.write(file_name, line)

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for X3DNA annotation extraction.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for input, output, and dummy paths.

    """

    parser = optparse.OptionParser("python x3dna.py -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_x3dna_obj(pdb_file, dummy_dir="/tmp"):
    """
        Run X3DNA `find_pair` and parse its output into an `X3DNA` object.
    
        Args:
            pdb_file (str): Input PDB file path.
            dummy_dir (str, optional): Temporary working directory root.
    
        Returns:
            X3DNA: Parsed X3DNA annotation object.
    
    Args:
        pdb_file (Any): Path to the input/output file.
        dummy_dir (Any): Directory path used by this operation."""

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        #x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        x3dna_path = config.get("Paths", "x3dna_path")
        os.environ['X3DNA'] = x3dna_path[:-4]
        # Get current working directory #
        cwd = os.getcwd()
        # Create tmp directory #
        tmp = os.path.join(dummy_dir, str(os.getpid()))
        if not os.path.exists(tmp): os.makedirs(tmp)
        # Change directory #
        os.chdir(tmp)
        # Exec process #
        process = subprocess.check_output([os.path.join(x3dna_path, "find_pair"), pdb_file, "3dna.out"], stderr=subprocess.STDOUT, env=os.environ)
        # Get X3DNA object #
        x3dna_obj = X3DNA("3dna.out")
        # Return to original directory #
        os.chdir(cwd)
        # Erase tmp directory #
        shutil.rmtree(tmp)
    except Exception as e:
        raise ValueError("Could not exec X3DNA for %s with error %s" % (pdb_file,e))

    return x3dna_obj

def main():
    """
    Run the command-line X3DNA extraction workflow.

    Workflow:
        1. Parse command-line options.
        2. Run X3DNA `find_pair` on the input PDB.
        3. Write raw X3DNA output to file or stdout.

    Args:
        None.

    Returns:
        None. X3DNA lines are written to file or stdout.
    """

    # Arguments & Options #
    options = parse_options()

    # Get X3DNA object #
    x3dna_obj = get_x3dna_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        x3dna_obj.write(os.path.abspath(options.output_file))
    else:
        for line in x3dna_obj._file_content:
            sys.stdout.write("%s\n" % line)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
