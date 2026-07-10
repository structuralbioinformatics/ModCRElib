import os, sys, re
import configparser
import copy
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

#-------------#
# Classes     #
#-------------#

class DSSP(object):
    """
    Parse and store per-residue DSSP accessibility and secondary structure.

    Object features:
        - Caches raw DSSP file lines in `_file_content` for round-trip writing.
        - Stores residue annotations in `_residues` keyed by
          `(pdb_chain, residue_number)`.
        - Exposes residue-level getters for:
            - accessible surface area (`get_accessible_surface_area`)
            - secondary-structure code (`get_secondary_structure`)
        - Supports residue presence checks (`has_residue`) and full dictionary
          retrieval (`get_dictionary_of_residues`).
    """

    def __init__(self, file_name):
        """
                Initialize a DSSP container from a DSSP output file.
        
                Args:
                    file_name (str): Path to a DSSP output file.
        
                Returns:
                    None.
        
        Args:
            file_name (Any): Value used by this routine."""
        self._file = file_name
        self._file_content = []
        self._residues = {}
        # Initialize #
        self._parse_file()

    def _parse_file(self):
        """Parse DSSP rows and cache relevant residue-level annotations."""
        for line in functions.parse_file(self._file):
            if not line.endswith(".") and "#" not in line and line != "":
                m = re.search("(\d+)\s\S\s[ACDEFGHIKLMNPQRSTVWY]", line)
                # Get DSSP info #
                if m:
                    pdb_chain = line[11:12]
                    residue_num = int(line[5:10])
                    accessible_surface_area = float(line[35:38])
                    secondary_structure = line[16:17]
                    self._residues[(pdb_chain, residue_num)] = (accessible_surface_area, secondary_structure)
            # Add line to file content #
            if line != "":
               self._file_content.append(line)

    def has_residue(self, pdb_chain, residue_num):
        """Return whether DSSP data exists for a residue.
        
        Args:
            pdb_chain (Any): PDB structure identifier and chain selector.
            residue_num (Any): Numeric identifier/index used by this routine."""
        if (pdb_chain, residue_num) in self._residues:
            return True

        return False
            

    def get_dictionary_of_residues(self):
        """Return all parsed residue annotations."""
        return copy.copy(self._residues)


    def get_accessible_surface_area(self, pdb_chain, residue_num):
        """Return DSSP accessible-surface area for one residue.
        
        Args:
            pdb_chain (Any): PDB structure identifier and chain selector.
            residue_num (Any): Numeric identifier/index used by this routine."""
        if self.has_residue(pdb_chain, residue_num):
            return copy.copy(self._residues[(pdb_chain, residue_num)][0])

        return None

    def get_secondary_structure(self, pdb_chain, residue_num):
        """Return DSSP secondary-structure code for one residue.
        
        Args:
            pdb_chain (Any): PDB structure identifier and chain selector.
            residue_num (Any): Numeric identifier/index used by this routine."""
        if self.has_residue(pdb_chain, residue_num):
            return copy.copy(self._residues[(pdb_chain, residue_num)][1])

        return None

    def write(self, file_name):
        """Write cached DSSP output lines to file.
        
        Args:
            file_name (Any): Value used by this routine."""
        for line in self._file_content:
            functions.write(file_name, line)

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for DSSP annotation extraction.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for input, output, and dummy paths.

    """

    parser = optparse.OptionParser("python dssp.py -i INPUT_FILE [--dummy DUMMY_DIR -o OUTPUT_FILE]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_dssp_obj(pdb_file, dummy_dir="/tmp"):
    """
        Run DSSP on a PDB file and parse results into a `DSSP` object.
    
        Args:
            pdb_file (str): Input PDB file path.
            dummy_dir (str, optional): Temporary working directory root.
    
        Returns:
            DSSP: Parsed DSSP object.
    
    Args:
        pdb_file (Any): Path to the input/output file.
        dummy_dir (Any): Directory path used by this operation."""

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        dssp_path = os.path.join(src_path, config.get("Paths", "dssp_path"))
        # Exec process #
        dssp_out="dssp_"+str(os.getpid())+".out"
        #print("Execute ",os.path.join(dssp_path, "dssp"), pdb_file, os.path.join(dummy_dir,dssp_out))
        process = subprocess.check_output([os.path.join(dssp_path, "mkdssp"), pdb_file, os.path.join(dummy_dir,dssp_out)], stderr=subprocess.STDOUT)
        # Get DSSP object #
        dssp_obj = DSSP(os.path.join(dummy_dir,dssp_out))
        # Remove DSSP file #
        os.remove(os.path.join(dummy_dir,dssp_out))
    except:
        print("Failed ",os.path.join(dssp_path, "mkdssp"),pdb_file,os.path.join(dummy_dir,dssp_out))
        raise ValueError("Could not exec DSSP for %s" % pdb_file)

    return dssp_obj

def main():
    """
    Run the command-line DSSP workflow.

    Workflow:
        1. Parse command-line options.
        2. Execute DSSP for the input PDB.
        3. Write DSSP output to file or stdout.

    Args:
        None.

    Returns:
        None. DSSP output lines are written to file or stdout.
    """

    # Arguments & Options #
    options = parse_options()

    # Get DSSP object #
    dssp_obj = get_dssp_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        dssp_obj.write(os.path.abspath(options.output_file))
    else:
        for line in dssp_obj._file_content:
            sys.stdout.write("%s\n" % line)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
