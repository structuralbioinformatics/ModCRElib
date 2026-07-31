import os, sys, re
import configparser
import copy
import numpy
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

class TMalign(object):
    """
    Parse TM-align output, including scores, transform, and alignments.

    Object features:
        - Stores raw TM-align stdout lines in `_file_content`.
        - Stores one or more TM-scores in `_tm_scores`.
        - Stores rigid-body superposition parameters as:
            - rotation matrix (`_matrix`)
            - translation vector (`_vector`)
        - Stores parsed alignment strings in `_query_alignment` and
          `_hit_alignment`.
    """

    def __init__(self, file_content):
        """
                Initialize a TM-align result container from command output lines.
        
                Args:
                    file_content (list): Raw TM-align stdout split by lines.
        
                Returns:
                    None.
        
        Args:
            file_content (Any): Value used by this routine."""
        self._file_content = file_content
        self._tm_scores = []
        self._matrix = numpy.identity(3, float)
        self._vector = numpy.zeros(3, float)
        self._query_alignment = None
        self._hit_alignment = None
        # Initialize #
        self._parse_file()

    def _parse_file(self):
        """Parse TM-scores, rigid transform, and alignment strings."""
        for line in self._file_content:
            # Capture TM-scores #
            m = re.search("^TM-score=\s*(\S+)", line)
            if m:
                self._tm_scores.append(float(m.group(1)))
            # Capture vector and matrix #
            m = re.search("^(0|1|2)\t\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$", line)
            if m:
                self._matrix[int(m.group(1))][0] = float(m.group(3))
                self._matrix[int(m.group(1))][1] = float(m.group(4))
                self._matrix[int(m.group(1))][2] = float(m.group(5))
                self._vector[int(m.group(1))] = float(m.group(2))
            # Capture alignment #
            m = re.search("^([\w-]+)$", line)
            if m:
                if self._query_alignment == None:
                    self._query_alignment = m.group(1)
                else:
                    self._hit_alignment = m.group(1)

    def get_tm_scores(self):
        """Return parsed TM-scores."""
        return self._tm_scores

    def get_matrix(self):
        """Return rotation matrix from TM-align superposition."""
        return self._matrix

    def get_vector(self):
        """Return translation vector from TM-align superposition."""
        return self._vector

    def get_query_alignment(self):
        """Return query alignment string."""
        return self._query_alignment

    def get_hit_alignment(self):
        """Return hit alignment string."""
        return self._hit_alignment

    def write(self, file_name):
        """Write raw TM-align output lines to file.
        
        Args:
            file_name (Any): Value used by this routine."""
        for line in self._file_content:
            functions.write(file_name, line)

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for TM-align structure comparison.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for input structures and outputs.

    """

    parser = optparse.OptionParser("python tmalign.py -a PDB_FILE_A -b PDB_FILE_B [--dummy DUMMY_DIR -o OUTPUT_FILE]")

    parser.add_option("-a", action="store", type="string", dest="pdb_file_a", help="PDB file", metavar="{filename}")
    parser.add_option("-b", action="store", type="string", dest="pdb_file_b", help="PDB file", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.pdb_file_a is None or options.pdb_file_b is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_tmalign_obj(pdb_file_a, pdb_file_b, dummy_dir="/tmp"):
    """
        Run TM-align and parse its output.
    
        Note:
            `pdb_file_b` is superimposed onto `pdb_file_a`.
    
        Args:
            pdb_file_a (str): Reference PDB file.
            pdb_file_b (str): Mobile PDB file.
            dummy_dir (str, optional): Reserved temporary directory argument.
    
        Returns:
            TMalign: Parsed TM-align result object.
    
    Args:
        pdb_file_a (Any): Value used by this routine.
        pdb_file_b (Any): Value used by this routine.
        dummy_dir (Any): Directory path used by this operation."""

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        tmalign_path = os.path.join(src_path, config.get("Paths", "tmalign_path"))
        # If Mac OS X... #
        try: 
            sys.platform == "darwin"
            process = subprocess.check_output([os.path.join(tmalign_path, "TMalign"), pdb_file_a, pdb_file_b], stderr=subprocess.STDOUT)
        # Else... #
        except:
            process = subprocess.check_output([os.path.join(tmalign_path, "TMalign"), pdb_file_a, pdb_file_b], stderr=subprocess.STDOUT)
        # Get TMalign object #
        tmalign_obj = TMalign(process.decode().split("\n"))
    except:
        raise ValueError("Could not exec TMalign for %s %s" % (pdb_file_a, pdb_file_b))
    return tmalign_obj

def main():
    """
    Run the command-line TM-align workflow.

    Workflow:
        1. Parse command-line options.
        2. Run TM-align on input structures.
        3. Write TM-align output to file or stdout.

    Args:
        None.

    Returns:
        None. TM-align output is written to file or stdout.
    """

    # Arguments & Options #
    options = parse_options()

    # Get TMalign object #
    tmalign_obj = get_tmalign_obj(os.path.abspath(options.pdb_file_a), os.path.abspath(options.pdb_file_b), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        tmalign_obj.write(os.path.abspath(options.output_file))
    else:
        for line in tmalign_obj._file_content:
            sys.stdout.write("%s\n" % line)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
