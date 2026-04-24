import os, sys, re
import configparser
import copy
import json
import numpy
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
from SBILib.data import aminoacids3to1, aminoacids_polarity_boolean

#ratios
# 102939     methylated binding files
# 233423     methylated non-binding files
# 2211745    unspecific bindings
# Total  2548107
# 980211 PBM triads standard modcre
# 9315   PDB triads standard modcre
# Total 989526
# ratio 0.3883377

#ratio =  0.3883377
ratio =  1.0

#-------------#
# Class       #
#-------------#

class Potentials(object):
    """
    Hold and query statistical potential tables loaded from disk.

    The object stores multiple PMF variants and the corresponding distance bins
    used to index per-distance scores.

    Object features:
        - Serialized source path and optional selective-loading mode
          (`_file`, `_potential`).
        - Distance-dependent PMF tables for global and context-aware models
          (`_pmf_3d`, `_pmf_3dc`, `_pmf_s3dc`, `_pmf_s3dc_dd`, `_pmf_pair`).
        - Single-value/context PMF tables for local and distance-independent
          variants (`_pmf_local`, `_pmf_s3dc_di`).
        - Distance-bin index map used to convert raw distances into PMF-array
          positions (`_distances`).
        - Accessors that resolve distance bins and return potential values for
          requested keys/potential families.
        - Serialization helper to write the complete PMF state back to disk.
    """

    def __init__(self, file_name=None, select_potential=None):
        """
        Build a potentials container and optionally load it from a file.

        Args:
            file_name (str, optional): Path to a serialized potentials file.
            select_potential (str, optional): Name of a specific potential to
            load (for example ``3d`` or ``s3dc``), or ``all``.

        Returns:
            None.
        """
        self._file = file_name
        self._pmf_3d = None
        self._pmf_3dc = None
        self._pmf_s3dc = None
        self._pmf_local = None
        self._pmf_pair = None
        self._pmf_s3dc_dd = None
        self._pmf_s3dc_di = None
        self._distances = {}
        self._potential = select_potential

        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        """
        Parse a serialized potentials file into in-memory PMF structures.

        Returns:
            None.
        """
        for line in functions.parse_file(self._file):
            potential, json_obj = line.strip("\n").split("\t")
            if self._potential is None or self._potential=="all" or self._potential == potential[4:]:
                if potential == "pmf_3d" :
                    self._pmf_3d = json.loads(json_obj)
                if potential == "pmf_3dc"  :
                    self._pmf_3dc = json.loads(json_obj)
                if potential == "pmf_s3dc"  :
                    self._pmf_s3dc = json.loads(json_obj)
                if potential == "pmf_s3dc_dd"  :
                    self._pmf_s3dc_dd = json.loads(json_obj)
                if potential == "pmf_s3dc_di"  :
                    self._pmf_s3dc_di = json.loads(json_obj)
                if potential == "pmf_local"  :
                    self._pmf_local = json.loads(json_obj)
                if potential == "pmf_pair"  :
                    self._pmf_pair = json.loads(json_obj)
            if potential == "distances":
                distances = json.loads(json_obj)
                for position in range(len(distances)):
                    self._distances[distances[position]] = position

    def _get_distance(self, distance):
        """
        Return the first stored distance bin upper bound above a value.

        Args:
            distance (float): Query distance in angstroms.

        Returns:
            float or None: Matching bin edge, or ``None`` if out of range.
        """
        for d in self._distances:
            if distance < d:
                return d

        return None

    def _get_distance_position(self, distance):
        """
        Get the index position associated with a distance bin edge.

        Args:
            distance (float): Distance bin edge value.

        Returns:
            int or None: Position in PMF vectors, or ``None`` if not present.
        """
        if distance in self._distances:
            return self._distances[distance]

        return None

    def get_score(self, potential, key=None, distance=None):
        """
        Retrieve a potential score for a key and distance.

        Args:
            potential (str): Potential type (``3d``, ``3dc``, ``s3dc``,
            ``s3dc_dd``, ``s3dc_di``, ``local``, or ``pair``).
            key (str, optional): Residue/environment key used by keyed PMFs.
            distance (float, optional): Query distance used to select a bin.

        Returns:
            float or None: Score value if available; otherwise ``None``.
        """
        distance = numpy.floor(self._get_distance(distance))
        position = self._get_distance_position(distance)
        try:
          if position is not None:
            if potential == "3d":
                return self._pmf_3d[position]
            if potential == "3dc":
                if key in self._pmf_3dc:
                    return self._pmf_3dc[key][position]
            if potential == "s3dc":
                if key in self._pmf_s3dc:
                    return self._pmf_s3dc[key][position]
            if potential == "s3dc_dd":
                if key in self._pmf_s3dc_dd:
                    return self._pmf_s3dc_dd[key][position]
            if potential == "pair":
                if key in self._pmf_pair:
                    return self._pmf_pair[key][position]
          if potential == "local":
            if key in self._pmf_local:
                return self._pmf_local[key][0]
          if potential == "s3dc_di":
            if key in self._pmf_s3dc_di:
                return self._pmf_s3dc_di[key][0]
        except:
          return None

    def write(self, file_name):
        """
        Write all stored PMF tables and distance bins to disk.

        Args:
            file_name (str): Output file path.

        Returns:
            None.
        """
        functions.write(file_name, "pmf_3d\t%s" % json.dumps(self._pmf_3d, separators=(",", ":")))
        functions.write(file_name, "pmf_3dc\t%s" % json.dumps(self._pmf_3dc, separators=(",", ":")))
        functions.write(file_name, "pmf_s3dc\t%s" % json.dumps(self._pmf_s3dc, separators=(",", ":")))
        functions.write(file_name, "pmf_s3dc_dd\t%s" % json.dumps(self._pmf_s3dc_dd, separators=(",", ":")))
        functions.write(file_name, "pmf_s3dc_di\t%s" % json.dumps(self._pmf_s3dc_di, separators=(",", ":")))
        functions.write(file_name, "pmf_local\t%s" % json.dumps(self._pmf_local, separators=(",", ":")))
        functions.write(file_name, "pmf_pair\t%s" % json.dumps(self._pmf_pair, separators=(",", ":")))
        functions.write(file_name, "distances\t%s" % json.dumps(self._distances, separators=(",", ":")))

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    Parse command-line options for methylation-aware potential generation.

    The parser configures the triad manifest to process, output destination,
    and optional post-processing steps (Taylor approximation, smoothing,
    Z-scoring, and binning mode).

    Returns:
        optparse.Values: Parsed CLI options describing input/output paths and
        potential-computation behavior.
    """

    parser = optparse.OptionParser("python spotentials.py -i input_file [-a --dummy=dummy_dir -o output_file -s -v -z -b]")

    parser.add_option("-a", default=False, action="store_true", dest="approach", help="Approach missing contacts with evolutionary transitions from BLOSUM62 (i.e. Taylor's approach; default = False)", metavar="{boolean}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (full path to triads files to be used to derive the potentials)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")
    parser.add_option("-s", default=False, action="store_true", dest="smooth", help="Smooth potentials (default = False)", metavar="{boolean}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
    parser.add_option("-z", default=False, action="store_true", dest="zscores", help="Calculate scanning-Zscores (default = False)", metavar="{boolean}")
    parser.add_option("-b", "--bins", default=False, action="store_true",  dest="computation", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error('missing arguments: type option \'-h\' for help')

    return options

def get_statistical_potentials(file_name, approach=False, smooth=False, zscores=False, computation=False, dummy_dir="/tmp"):
    """
    Derive statistical PMFs from triad-contact input files.

    Args:
        file_name (str): File listing triad files to aggregate.
        approach (bool, optional): Apply Taylor-like amino-acid transition
        approximation to fill sparse entries.
        smooth (bool, optional): Smooth PMFs by local bin averaging.
        zscores (bool, optional): Convert PMFs to per-context Z-scores.
        computation (bool, optional): Use binned mode instead of cumulative
        mode for frequency counting.
        dummy_dir (str, optional): Unused compatibility argument.

    Returns:
        tuple: ``(pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di,
        pmf_local, pmf_pair, distances)``.
    """

    # Initialize #
    if computation:
       bin_distance = float(config.get("Parameters", "bin_distance_bins"))
       smooth_sample= int(config.get("Parameters", "smooth_potential_bins"))
    else:
       bin_distance = float(config.get("Parameters", "bin_distance_accu"))
       smooth_sample= int(config.get("Parameters", "smooth_potential_accu"))
       
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    distances = list(numpy.arange(0, max_contact_distance + bin_distance, bin_distance)) # +bin_distance because of how range works
    kB = 8.31441e-3 # Boltzmann constant
    T = 298.0       # standard temperature

    # Get frequencies #
    
    f_dab, f_a_dab, f_a_dab_oa, f_a_b_dab, f_dab_oa_ob, f_a_b_dab_oa_ob = get_frequencies(file_name, computation=computation, approach=approach)
    # Calculate potentials #
    pmf_3d = [None] * len(distances)
    pmf_3dc = {}
    pmf_s3dc = {}
    pmf_s3dc_dd = {}
    pmf_s3dc_di = {}
    pmf_local = {}
    pmf_pair = {}
    # PMF 3D #
    for dr in range(len(distances)):
        if f_dab[dr] > 0:
            pmf_3d[dr] = kB * T * numpy.log(f_dab[dr] / float(sum(f_dab)))
    # PMF 3DC #
    for oa_ob in f_dab_oa_ob:
        pmf_3dc[oa_ob] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_dab_oa_ob[oa_ob][dr] > 0:
                pmf_3dc[oa_ob][dr] = (kB * T * numpy.log(f_dab_oa_ob[oa_ob][dr] / float(sum(f_dab_oa_ob[oa_ob])))) - pmf_3d[dr]
    # PMF S3DC #
    for a_b_oa_ob in f_a_b_dab_oa_ob:
        a_oa, b_ob = a_b_oa_ob.split(";")
        a_list = a_oa.split("-")
        b_list = b_ob.split("-")
        a = a_list.pop(0)
        b = b_list.pop(0)
        oa = "-".join(a_list)
        ob = "-".join(b_list)
        a_b = a + ";" + b
        oa_ob = oa + ";" + ob
        pmf_s3dc[a_b_oa_ob] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_a_b_dab_oa_ob[a_b_oa_ob][dr] > 0:
                pmf_s3dc[a_b_oa_ob][dr] = pmf_3d[dr] + pmf_3dc[oa_ob][dr] - (kB * T * numpy.log(f_a_b_dab_oa_ob[a_b_oa_ob][dr] / float(sum(f_a_b_dab_oa_ob[a_b_oa_ob]))))
    # PMF S3DC dd #
    for a_b_oa_ob in f_a_b_dab_oa_ob:
        pmf_s3dc_dd[a_b_oa_ob] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_a_b_dab_oa_ob[a_b_oa_ob][dr] > 0:
                pmf_s3dc_dd[a_b_oa_ob][dr] = - (kB * T * numpy.log(f_a_b_dab_oa_ob[a_b_oa_ob][dr] / float(f_dab[dr])))
    # PMF S3DC di #
    for a_b_oa_ob in f_a_b_dab_oa_ob:
        pmf_s3dc_di[a_b_oa_ob] = [None]
        if sum(f_a_b_dab_oa_ob[a_b_oa_ob]) > 0:
            pmf_s3dc_di[a_b_oa_ob] = [kB * T * numpy.log(sum(f_a_b_dab_oa_ob[a_b_oa_ob]) / float(sum(f_dab)))]
    # PMF local #
    for a_oa in f_a_dab_oa:
        pmf_local[a_oa] = [None]
        if sum(f_a_dab_oa[a_oa]) > 0:
            a_list = a_oa.split("-")
            a = a_list.pop(0)
            pmf_local[a_oa] = [(kB * T * numpy.log(sum(f_a_dab_oa[a_oa]) / float(sum(f_dab)))) - (kB * T * numpy.log(sum(f_a_dab[a]) / float(sum(f_dab))))]
    # PMF pair #
    for a_b in f_a_b_dab:
        pmf_pair[a_b] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_a_b_dab[a_b][dr] > 0:
                pmf_pair[a_b][dr] = pmf_3d[dr] - (kB * T * numpy.log(f_a_b_dab[a_b][dr] / float(sum(f_a_b_dab[a_b]))))

    # Approach potentials #
    if approach:
        # PMF S3DC #
        approached = []
        for a_b_oa_ob in pmf_s3dc:
            approached.append([a_b_oa_ob, approach_pmf(a_b_oa_ob, pmf_s3dc, "-")])
        for a_b_oa_ob, approached_pmf in approached:
            pmf_s3dc[a_b_oa_ob] = approached_pmf
        # PMF S3DC dd #
        approached = []
        for a_b_oa_ob in pmf_s3dc_dd:
            approached.append([a_b_oa_ob, approach_pmf(a_b_oa_ob, pmf_s3dc_dd, "-")])
        for a_b_oa_ob, approached_pmf in approached:
            pmf_s3dc_dd[a_b_oa_ob] = approached_pmf
        # PMF S3DC di #
        approached = []
        for a_b_oa_ob in pmf_s3dc_di:
            approached.append([a_b_oa_ob, approach_pmf(a_b_oa_ob, pmf_s3dc_di, "+")])
        for a_b_oa_ob, approached_pmf in approached:
            pmf_s3dc_di[a_b_oa_ob] = approached_pmf
        # PMF pair #
        approached = []
        for a_b in pmf_pair:
            approached.append([a_b, approach_pmf(a_b, pmf_pair, "-")])
        for a_b, approached_pmf in approached:
            pmf_pair[a_b] = approached_pmf

    # Calculate Z-scores #
    if zscores:
        pmf_3dc = calculate_zscores(pmf_3dc, "pmf_3dc", bin_distance)
        pmf_s3dc = calculate_zscores(pmf_s3dc, "pmf_s3dc", bin_distance)
        pmf_s3dc_dd = calculate_zscores(pmf_s3dc_dd, "pmf_s3dc_dd", bin_distance)
        pmf_s3dc_di = calculate_zscores(pmf_s3dc_di, "pmf_s3dc_di", bin_distance)
        pmf_local = calculate_zscores(pmf_local, "pmf_local", bin_distance)
        pmf_pair = calculate_zscores(pmf_pair, "pmf_pair", bin_distance)

    # Smooth potentials #
    if smooth:
        # PMF 3D #
        for dr in range(len(distances)):
            pmf_3d[dr] = smooth_bin(pmf_3d, dr,  smooth_sample)
        # PMF 3DC #
        for envpair in pmf_3dc:
            for dr in range(len(distances)):
                pmf_3dc[envpair][dr] = smooth_bin(pmf_3dc[envpair], dr, smooth_sample )
        # PMF S3DC #
        for pair in pmf_s3dc:
            for dr in range(len(distances)):
                pmf_s3dc[pair][dr] = smooth_bin(pmf_s3dc[pair], dr, smooth_sample )
        # PMF S3DC dd #
        for pair in pmf_s3dc_dd:
            for dr in range(len(distances)):
                pmf_s3dc_dd[pair][dr] = smooth_bin(pmf_s3dc_dd[pair], dr, smooth_sample ) 
        # PMF pair #
        for pair in pmf_pair:
            for dr in range(len(distances)):
                pmf_pair[pair][dr] = smooth_bin(pmf_pair[pair], dr, smooth_sample )


    return pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances

def get_frequencies(file_name, computation=False, approach=False):
    """
    Extract frequency tables required to compute all statistical PMFs.

    Args:
        file_name (str): File listing triad files to parse.
        computation (bool, optional): Use discrete bins (True) or cumulative
        counting (False).
        approach (bool, optional): Initialize extra amino-acid substitutions
        for Taylor-like approximation.

    Returns:
        tuple: ``(f_dab, f_a_dab, f_a_dab_oa, f_a_b_dab, f_dab_oa_ob,
        f_a_b_dab_oa_ob)`` frequency structures.
    """

    # Initialize #
    if computation:
       bin_distance = float(config.get("Parameters", "bin_distance_bins"))
    else:
       bin_distance = float(config.get("Parameters", "bin_distance_accu"))
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    distances = list(numpy.arange(0, max_contact_distance + bin_distance, bin_distance)) # +bin_distance because of how range works
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    path = file_name.split("/")
    pathbase = os.path.dirname(scripts_path)

    # The following arguments arguments are defined as in the papers:
    # "a" and "b" are the residues "a" and "b"
    # "oa" and "ob" are the environments of residues "a" and "b"
    # "a_oa" and "b_ob" are the residues "a" and "b" in their respective environments "oa" and "ob"
    # "dab" is the distance between residues "a" and "b"
    # "a_b_oa_ob" is a contact between residue "a_oa" and residue "b_ob"
    f_dab = [0] * len(distances) # contacts at a max. distance "dab"
    f_a_dab = {}                 # contacts at a max. distance "dab" involving residue "a"
    f_a_dab_oa = {}              # contacts at a max. distance "dab" involving residue "a" in environment "oa"
    f_a_b_dab = {}               # contacts at a max. distance "dab" involving residues "a" and "b"
    f_dab_oa_ob = {}             # contacts at a max. distance "dab" involving environments "oa" and "ob"
    f_a_b_dab_oa_ob = {}         # contacts at a max. distance "dab" involving residues "a" and "b" in their respective environments "oa" and "ob"

    # For each line... #
    sys.stderr.write("\t\t--Frequencies: Read File %s\n"%(file_name))

    for line in functions.parse_file(file_name):
        # Skip line if not a file name #
        if line.startswith("#"): continue
        # Try... #
        if "methylation" in line:
           factor = ratio
        else:
           factor = 1.0
        try:
            # For each line... #
            if pathbase in line: 
                triad_file=line
            else:
                triad_file = os.path.join(pathbase, line)
            sys.stderr.write("\t\t--Frequencies: add triad %s\n"%triad_file)
            for line in functions.parse_file(triad_file):
                # Skip line if not a triad #
                if line.startswith("#"): continue
                # Get triad #
                line = line.split(";")
                dab = line[2]
                # Adjust distance #
                dab = adjust_distance_to_bin(float(dab),bin_distance)
                # Skip if distance is too large #
                if dab > max_contact_distance: continue
                a_list = line[0].split("-")
                b_list = line[1].split("-")
                a = a_list.pop(0)
                b = b_list.pop(0)
                oa = "-".join(a_list)
                ob = "-".join(b_list)
                a_oa = a + "-" + oa
                b_ob = b + "-" + ob
                a_b_oa_ob = a_oa + ";" + b_ob
                a_b = a + ";" + b
                oa_ob = oa + ";" + ob
                if approach:
                    # Initialize #
                    for aa in aminoacids:
                        hydrophobicity = "N"
                        if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                            hydrophobicity = "P"
                        oaa = hydrophobicity + "-" + a_list[1] + "-" + a_list[2]
                        aa_oaa = aa + "-" + oaa
                        aa_b = aa + ";" + b
                        oaa_ob = oaa + ";" + ob
                        aa_b_oaa_ob = aa_oaa + ";" + b_ob
                        # Initialize frequencies #
                        if aa not in f_a_dab:
                            f_a_dab[aa] = [0] * len(distances)
                        # For local only do amino acids #
                        if aa_oaa not in f_a_dab_oa:
                            f_a_dab_oa[aa_oaa] = [0] * len(distances)
                        if aa_b not in f_a_b_dab:
                            f_a_b_dab[aa_b] = [0] * len(distances)
                        if oaa_ob not in f_dab_oa_ob:
                            f_dab_oa_ob[oaa_ob] = [0] * len(distances)
                        if aa_b_oaa_ob not in f_a_b_dab_oa_ob:
                            f_a_b_dab_oa_ob[aa_b_oaa_ob] = [0] * len(distances)
                # Initialize frequencies #
                if a not in f_a_dab:
                    f_a_dab[a] = [0] * len(distances)
                if b not in f_a_dab:
                    f_a_dab[b] = [0] * len(distances)
                # For local only do amino acids #
                if a_oa not in f_a_dab_oa:
                    f_a_dab_oa[a_oa] = [0] * len(distances)
                if a_b not in f_a_b_dab:
                    f_a_b_dab[a_b] = [0] * len(distances)
                if oa_ob not in f_dab_oa_ob:
                    f_dab_oa_ob[oa_ob] = [0] * len(distances)
                if a_b_oa_ob not in f_a_b_dab_oa_ob:
                    f_a_b_dab_oa_ob[a_b_oa_ob] = [0] * len(distances)
                # Update #
                for dr in range(1,len(distances)):
                    if computation :
                        if distances[dr-1] < dab <= distances[dr]:
                            f_dab[dr] += 1*factor
                            f_a_dab[a][dr] += 1*factor
                            f_a_dab[b][dr] += 1*factor
                            f_a_b_dab[a_b][dr] += 1*factor
                            # For local only do amino acids #
                            f_a_dab_oa[a_oa][dr] += 1*factor
                            f_a_b_dab_oa_ob[a_b_oa_ob][dr] += 1*factor
                            f_dab_oa_ob[oa_ob][dr] += 1*factor
                    else:
                        if dab <= distances[dr]:
                            f_dab[dr] += 1*factor
                            f_a_dab[a][dr] += 1*factor
                            f_a_dab[b][dr] += 1*factor
                            f_a_b_dab[a_b][dr] += 1*factor
                            # For local only do amino acids #
                            f_a_dab_oa[a_oa][dr] += 1*factor
                            f_a_b_dab_oa_ob[a_b_oa_ob][dr] += 1*factor
                            f_dab_oa_ob[oa_ob][dr] += 1*factor
        # Except... #
        except: pass

    return f_dab, f_a_dab, f_a_dab_oa, f_a_b_dab, f_dab_oa_ob, f_a_b_dab_oa_ob

def adjust_distance_to_bin(distance,bin_distance):
    """
    Snap a raw distance to the lower edge of its bin.

    Args:
        distance (float): Input distance.
        bin_distance (float): Width of each distance bin.

    Returns:
        float: Bin-aligned distance value.
    """

    # Initialize #
    distances = []

    for i in numpy.arange(numpy.floor(distance), numpy.floor(distance) + 1, bin_distance):
        if i <= distance:
            distances.append(i)

    return distances[-1]


def approach_pmf(a_b_oa_ob, pmf, symbol):
    """
    Estimate missing PMF values with amino-acid transition probabilities.

    Args:
        a_b_oa_ob (str): Residue/environment key to approach.
        pmf (dict): PMF table keyed by residue/environment pair.
        symbol (str): ``"+"`` for attractive terms, ``"-"`` for repulsive
        terms in the approximation formula.

    Returns:
        list: Approached PMF profile for the requested key.
    """

    # Initialize #
    kB = 8.31441e-3 # Boltzmann constant
    T = 298.0       # standard temperature
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    approached_pmf = copy.copy(pmf[a_b_oa_ob])
    a_oa, b_ob = a_b_oa_ob.split(";")
    if len(a_oa.split("-"))>1:
       a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
       env=1
    else:
       a= a_oa
       b= b_ob
       env=0

    # For each distance bin... #
    for dr in range(len(approached_pmf)):
            # Approach PMF for amino acid at distance bin #
            # Function to determine best transition amino acid #
            f_max = []
            for aa in aminoacids:
                # Skip if amino acid is a... #
                if aa == a: continue
                if (env):
                  hydrophobicity = "N"
                  if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                    hydrophobicity = "P"
                  oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                  aa_oaa = aa + "-" + oaa
                  aa_b_oaa_ob = aa_oaa + ";" + b_ob
                else:
                  aa_b_oaa_ob = aa + ";" + b
                if aa_b_oaa_ob in pmf:
                    if pmf[aa_b_oaa_ob][dr] is not None:
                        if symbol == "+":
                            f_max.append((aa, get_transition_probability(A=a, B=aa), pmf[aa_b_oaa_ob][dr], get_transition_probability(A=a, B=aa) * numpy.exp(pmf[aa_b_oaa_ob][dr] / (kB * T))))
                        else:
                            f_max.append((aa, get_transition_probability(A=a, B=aa), pmf[aa_b_oaa_ob][dr], get_transition_probability(A=a, B=aa) * numpy.exp(-pmf[aa_b_oaa_ob][dr] / (kB * T))))
            # Skip if no contacts were found #
            if len(f_max) > 0:
                # Get transition probability and energy of amino acid maximizing previous function #
                f_max.sort(key=lambda x: x[-1], reverse=True)
                max_aa = f_max[0][0]
                max_probability = f_max[0][1]
                max_energy = f_max[0][2]
                # Fill the summation using Taylor's approach #
                summation = 0.0
                for aa in aminoacids:
           	         # Skip if amino acid is a... #
                    if aa == a: continue
                    # Skip if amino acid is b max... #
                    if aa == max_aa: continue
                    if (env):
                      hydrophobicity = "N"
                      if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                        hydrophobicity = "P"
                      oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                      aa_oaa = aa + "-" + oaa
                      aa_b_oaa_ob = aa_oaa + ";" + b_ob
                    else:
                      aa_b_oaa_ob = aa + ";" + b
                    if aa_b_oaa_ob in pmf:
                        if pmf[aa_b_oaa_ob][dr] is not None:
                            probability = get_transition_probability(A=a, B=aa)
                            energy = pmf[aa_b_oaa_ob][dr]
                            if symbol == "+":
                                summation += (probability / max_probability) * numpy.exp((energy - max_energy) / (kB * T))
                            else:
                                summation += (probability / max_probability) * numpy.exp(-((energy - max_energy) / (kB * T)))
                if approached_pmf[dr] is None:
                   if symbol == "+":
                    approached_pmf[dr] = kB * T * ((numpy.log(max_probability) + (max_energy / (kB * T)) + summation))
                   else:
                    approached_pmf[dr] = -kB * T * ((numpy.log(max_probability) - (max_energy / (kB * T)) + summation))
                else:
                   if symbol == "+":
                    #if (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp( (pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) < 1 :
                       #approached_pmf[dr] +=  kB * T * ( (numpy.log(get_transition_probability(A=a,B=a)) + (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp( (pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )
                    if (max_probability * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) )) < 1 :
                       approached_pmf[dr] +=  kB * T * (  max_probability * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )
                   else:
                    #if (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) < 1 :
                       #approached_pmf[dr] += -kB * T * ( (numpy.log(get_transition_probability(A=a,B=a)) + (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )
                    if (max_probability * numpy.exp(+(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) )) < 1 :
                       approached_pmf[dr] += -kB * T * (  max_probability * numpy.exp(+(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )


    return approached_pmf

def get_transition_probability(A, B):
    """
    Compute transition weight between two amino acids using BLOSUM62.

    Args:
        A (str): Source amino acid in 3-letter code.
        B (str): Target amino acid in 3-letter code.

    Returns:
        float: Frequency-weighted transition probability-like score.
    """
    
    l = 0.347
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    frequences = [0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.034, 0.073]
    blosum62   = [['4', '-1', '-2', '-2', '0', '-1', '-1', '0', '-2', '-1', '-1', '-1', '-1', '-2', '-1', '1', '0', '-3', '-2', '0'],
                 ['-1', '5', '0', '-2', '-3', '1', '0', '-2', '0', '-3', '-2', '2', '-1', '-3', '-2', '-1', '-1', '-3', '-2', '-3'],
                 ['-2', '0', '6', '1', '-3', '0', '0', '0', '1', '-3', '-3', '0', '-2', '-3', '-2', '1', '0', '-4', '-2', '-3'],
                 ['-2', '-2', '1', '6', '-3', '0', '2', '-1', '-1', '-3', '-4', '-1', '-3', '-3', '-1', '0', '-1', '-4', '-3', '-3'],
                 ['0', '-3', '-3', '-3', '9', '-3', '-4', '-3', '-3', '-1', '-1', '-3', '-1', '-2', '-3', '-1', '-1', '-2', '-2', '-1'],
                 ['-1', '1', '0', '0', '-3', '5', '2', '-2', '0', '-3', '-2', '1', '0', '-3', '-1', '0', '-1', '-2', '-1', '-2'],
                 ['-1', '0', '0', '2', '-4', '2', '5', '-2', '0', '-3', '-3', '1', '-2', '-3', '-1', '0', '-1', '-3', '-2', '-2'],
                 ['0', '-2', '0', '-1', '-3', '-2', '-2', '6', '-2', '-4', '-4', '-2', '-3', '-3', '-2', '0', '-2', '-2', '-3', '-3'],
                 ['-2', '0', '1', '-1', '-3', '0', '0', '-2', '8', '-3', '-3', '-1', '-2', '-1', '-2', '-1', '-2', '-2', '2', '-3'],
                 ['-1', '-3', '-3', '-3', '-1', '-3', '-3', '-4', '-3', '4', '2', '-3', '1', '0', '-3', '-2', '-1', '-3', '-1', '3'],
                 ['-1', '-2', '-3', '-4', '-1', '-2', '-3', '-4', '-3', '2', '4', '-2', '2', '0', '-3', '-2', '-1', '-2', '-1', '1'],
                 ['-1', '2', '0', '-1', '-3', '1', '1', '-2', '-1', '-3', '-2', '5', '-1', '-3', '-1', '0', '-1', '-3', '-2', '-2'],
                 ['-1', '-1', '-2', '-3', '-1', '0', '-2', '-3', '-2', '1', '2', '-1', '5', '0', '-2', '-1', '-1', '-1', '-1', '1'],
                 ['-2', '-3', '-3', '-3', '-2', '-3', '-3', '-3', '-1', '0', '0', '-3', '0', '6', '-4', '-2', '-2', '1', '3', '-1'],
                 ['-1', '-2', '-2', '-1', '-3', '-1', '-1', '-2', '-2', '-3', '-3', '-1', '-2', '-4', '7', '-1', '-1', '-4', '-3', '-2'],
                 ['1', '-1', '1', '0', '-1', '0', '0', '0', '-1', '-2', '-2', '0', '-1', '-2', '-1', '4', '1', '-3', '-2', '-2'],
                 ['0', '-1', '0', '-1', '-1', '-1', '-1', '-2', '-2', '-1', '-1', '-1', '-1', '-2', '-1', '1', '5', '-2', '-2', '0'],
                 ['-3', '-3', '-4', '-4', '-2', '-2', '-3', '-2', '-2', '-3', '-2', '-3', '-1', '1', '-4', '-3', '-2', '11', '2', '-3'],
                 ['-2', '-2', '-2', '-3', '-2', '-1', '-2', '-3', '2', '-1', '-1', '-2', '-1', '3', '-3', '-2', '-2', '2', '7', '-1'],
                 ['0', '-3', '-3', '-3', '-1', '-2', '-2', '-3', '-3', '3', '1', '-2', '1', '-1', '-2', '-2', '0', '-3', '-1', '4']]

    return frequences[aminoacids.index(A)] * frequences[aminoacids.index(B)] * numpy.exp(l * int(blosum62[aminoacids.index(A)][aminoacids.index(B)]))


def smooth_bin(array, position, n):
    """
    Smooth one PMF bin by averaging values in a local window.

    Args:
        array (list): PMF values along distance bins.
        position (int): Index of the bin to smooth.
        n (int): Number of bins to include on each side.

    Returns:
        float or None: Local mean ignoring ``None`` entries.
    """

    values = []
    first = max(0,position-n)
    last  = min(position+n+1,len(array))
    for i in range(first,last):
        if array[i] != None: values.append(array[i])
    if len(values)>0: 
       return numpy.mean(values)
    else:
       return None

def calculate_zscores(pmf, pmf_name, bin_distance):
    """
    Convert PMF values into Z-scores within comparable amino-acid contexts.

    Args:
        pmf (dict): Potential table to standardize.
        pmf_name (str): Potential family name (for example ``pmf_3dc`` or
        ``pmf_pair``) controlling grouping logic.
        bin_distance (float): Distance-bin width used to rebuild bin indices.

    Returns:
        dict or None: Z-scored PMF table, or ``None`` for unknown PMF names.
    """

    # Initialize #
    zpmf = {}
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    distances = list(numpy.arange(0, max_contact_distance + bin_distance, bin_distance)) # +bin_distance because of how range works


    # Z-score PMF 3DC #
    if pmf_name == "pmf_3dc":
        # For each environment pair... #
        for oa_ob in pmf:
            oa, ob = oa_ob.split(";")
            hydrophobicity, exposure, secondary_structure = oa.split("-")
            zpmf[oa_ob] = [None] * len(distances)
            for dr in range(len(distances)):
                if pmf[oa_ob][dr] is not None:
                    energy = []
                    for aa in aminoacids:
                        hydrophobicity = "N"
                        if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                            hydrophobicity = "P"
                        oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                        oaa_ob = oaa + ";" + ob
                        if oaa_ob in pmf:
                            if pmf[oaa_ob][dr] is not None:
                                energy.append(pmf[oaa_ob][dr])
                    mean = numpy.mean(energy)
                    std = numpy.std(energy)
                    if std != 0.0:
                        zpmf[oa_ob][dr] = (pmf[oa_ob][dr] - mean) / std
        return zpmf
    # Z-score PMF S3DC or S3DC dd #
    if pmf_name == "pmf_s3dc" or pmf_name == "pmf_s3dc_dd":
        # For each residue-environment pair... #
        for a_b_oa_ob in pmf:
            a_oa, b_ob = a_b_oa_ob.split(";")
            a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
            zpmf[a_b_oa_ob] = [None] * len(distances)
            for dr in range(len(distances)):
                if pmf[a_b_oa_ob][dr] is not None:
                    energy = []
                    for aa in aminoacids:
                        hydrophobicity = "N"
                        if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                            hydrophobicity = "P"
                        oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                        aa_oaa = aa + "-" + oaa
                        aa_b_oaa_ob = aa_oaa + ";" + b_ob
                        if aa_b_oaa_ob in pmf:
                            if pmf[aa_b_oaa_ob][dr] is not None:
                                energy.append(pmf[aa_b_oaa_ob][dr])
                    mean = numpy.mean(energy)
                    std = numpy.std(energy)
                    if std != 0.0:
                        zpmf[a_b_oa_ob][dr] = (pmf[a_b_oa_ob][dr] - mean) / std
        return zpmf
    # Z-score PMF S3DC di #
    if pmf_name == "pmf_s3dc_di":
        # For each residue-environment pair... #
        for a_b_oa_ob in pmf:
            a_oa, b_ob = a_b_oa_ob.split(";")
            a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
            zpmf[a_b_oa_ob] = [None]
            if pmf[a_b_oa_ob][0] is not None:
                energy = []
                for aa in aminoacids:
                    hydrophobicity = "N"
                    if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                        hydrophobicity = "P"
                    oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                    aa_oaa = aa + "-" + oaa
                    aa_b_oaa_ob = aa_oaa + ";" + b_ob
                    if aa_b_oaa_ob in pmf:
                        if pmf[aa_b_oaa_ob][0] is not None:
                            energy.append(pmf[aa_b_oaa_ob])
                mean = numpy.mean(energy)
                std = numpy.std(energy)
                if std != 0.0:
                    zpmf[a_b_oa_ob][0] = (pmf[a_b_oa_ob][0] - mean) / std
        return zpmf
    # Z-score PMF local #
    if pmf_name == "pmf_local":
        # For each residue-environment... #
        for a_oa in pmf:
            a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
            zpmf[a_oa] = [None]
            if pmf[a_oa][0] is not None:
                energy = []
                for aa in aminoacids:
                    hydrophobicity = "N"
                    if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                        hydrophobicity = "P"
                    oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                    aa_oaa = aa + "-" + oaa
                    if aa_oaa in pmf:
                        if pmf[aa_oaa][0] is not None:
                            energy.append(pmf[aa_oaa])
                mean = numpy.mean(energy)
                std = numpy.std(energy)
                if std != 0.0:
                    zpmf[a_oa][0] = (pmf[a_oa][0] - mean) / std
        return zpmf
    # Z-score PMF pair #
    if pmf_name == "pmf_pair":
        # For each residue pair... #
        for a_b in pmf:
            a, b = a_b.split(";")
            zpmf[a_b] = [None] * len(distances)
            for dr in range(len(distances)):
                if pmf[a_b][dr] is not None:
                    energy = []
                    for aa in aminoacids:
                        aa_b = aa + ";" + b
                        if aa_b in pmf:
                            if pmf[aa_b][dr] is not None:
                                energy.append(pmf[aa_b][dr])
                    mean = numpy.mean(energy)
                    std = numpy.std(energy)
                    if std != 0.0:
                        zpmf[a_b][dr] = (pmf[a_b][dr] - mean) / std
        return zpmf

    return None

def main():
    """
    Run the methylation-aware statistical-potential generation workflow.

    Workflow:
        1. Parse runtime options and resolve absolute input/output paths.
        2. Read the triad manifest and compute all PMF potential families.
        3. Populate a :class:`Potentials` container with the computed tables.
        4. Serialize the resulting potentials to file or print them to stdout.

    Returns:
        None. Potentials are written to disk or emitted to standard output.
    """
    # Arguments & Options #
    options = parse_options()

    # Get statistical potentials #
    pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = get_statistical_potentials(os.path.abspath(options.input_file), options.approach, options.smooth, options.zscores, options.computation, os.path.abspath(options.dummy_dir))

    # Initialize #
    potentials_obj = Potentials()
    # Assign the potentials #
    potentials_obj._pmf_3d = pmf_3d
    potentials_obj._pmf_3dc = pmf_3dc
    potentials_obj._pmf_s3dc = pmf_s3dc
    potentials_obj._pmf_s3dc_dd = pmf_s3dc_dd
    potentials_obj._pmf_s3dc_di = pmf_s3dc_di
    potentials_obj._pmf_local = pmf_local
    potentials_obj._pmf_pair = pmf_pair
    potentials_obj._distances = distances

    # Write output #
    if options.output_file is not None:
        potentials_obj.write(os.path.abspath(options.output_file))
    else:
        sys.stdout.write("pmf_3d\t%s\n" % json.dumps(potentials_obj._pmf_3d, separators=(",", ":")))
        sys.stdout.write("pmf_3dc\t%s\n" % json.dumps(potentials_obj._pmf_3dc, separators=(",", ":")))
        sys.stdout.write("pmf_s3dc\t%s\n" % json.dumps(potentials_obj._pmf_s3dc, separators=(",", ":")))
        sys.stdout.write("pmf_s3dc_dd\t%s\n" % json.dumps(potentials_obj._pmf_s3dc_dd, separators=(",", ":")))
        sys.stdout.write("pmf_s3dc_di\t%s\n" % json.dumps(potentials_obj._pmf_s3dc_di, separators=(",", ":")))
        sys.stdout.write("pmf_local\t%s\n" % json.dumps(potentials_obj._pmf_local, separators=(",", ":")))
        sys.stdout.write("pmf_pair\t%s\n" % json.dumps(potentials_obj._pmf_pair, separators=(",", ":")))
        sys.stdout.write("distances\t%s\n" % json.dumps(potentials_obj._distances, separators=(",", ":")))


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()

