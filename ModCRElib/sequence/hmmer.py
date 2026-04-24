import os, sys, re
import configparser
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


# Imports jbonet's module #
from SBILib.data import aminoacids1to3

#-------------#
# Classes     #
#-------------#

class Hmmer(object):
    """
    Parse and store HMMER text output for one query.

    Object features:
        - Raw HMMER output lines preserved from one `hmmsearch` execution
          (`_file_content`).
        - Query identifier and query-model length parsed from header lines
          (`_query`, `_query_length`).
        - Ordered hit collection (`_hits`) where each hit is a `HmmerHit`
          instance containing one or more domain segments (`HmmerHSP`).
        - Parsing logic for inclusion markers, domain boundaries, aligned
          query/hit sequences, and similarity counts.
        - Export helper that writes compact homolog rows compatible with the
          downstream homolog-processing pipeline.
    """

    def __init__(self, file_content):
        """
        Initialize and parse HMMER plain-text output lines.

        Args:
            file_content (list[str]): Line-wise `hmmsearch` output text.

        Returns:
            None.
        """
        self._file_content = file_content
        self._query = None
        self._query_length = None
        self._hits = []
        # Initialize #
        self._parse_file()

    def _get_file(self):
        """Return the source file path when available."""
        return self._file

    def _parse_file(self):
        """Parse HMMER output lines and populate hit/HSP containers."""
        for line in self._file_content:
            if line.startswith("#"): continue
            if "[No individual domains that satisfy reporting thresholds (although complete target did)]" in line:
                self._hits.pop(-1)
                continue
            m = re.search("Query:\s+(\S+)\s+\[M=(\d+)\]", line)
            if m:
                self._query = m.group(1)
                self._query_length = int(m.group(2))
            m = re.search(">>\s+(.+)", line)
            if m:
                self._hits.append(HmmerHit(m.group(1)))
            m = re.search("(\d+)\s+(\?|\!)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+.*\s+(\d+)\s+(\d+)\s+.*\s+\d+\s+\d+\S+", line)
            if m:
                domain = int(m.group(1))
                inclusion = m.group(2)
                score = float(m.group(3))
                bias = float(m.group(4))
                e_value = float(m.group(5))
                query_from = int(m.group(6))
                query_to = int(m.group(7))
                hit_from = int(m.group(8))
                hit_to = int(m.group(9))
                self._hits[-1].add_hsp(HmmerHSP(domain, inclusion, score, bias, e_value, query_from, query_to, hit_from, hit_to))
            m = re.search("== domain (\d+)", line)
            if m:
                domain = int(m.group(1))
            m = re.search("\S+\s+(\d+)\s([\w\-\.]+)\s(\d+)", line)
            if m:
                if len(re.findall("\w", m.group(2))) == int(m.group(3)) - int(m.group(1)) + 1:
                    if len(self._hits) > 0:
                        if len(self._hits[-1]._hsps) > 0:
                            self._hits[-1]._hsps[domain - 1].add_sequence(m.group(2).upper().replace(".", "-"))
            if "+" in line:
                if len(self._hits) > 0:
                    if len(self._hits[-1]._hsps) > 0:
                        if len(self._hits[-1]._hsps[domain - 1]._query_sequence) > len(self._hits[-1]._hsps[domain - 1]._hit_sequence):
                            self._hits[-1]._hsps[domain - 1].add_similarities(len(re.findall("\+", line)))

    def get_hits(self):
        """Return hits sorted by the best-domain E-value."""
        return sorted(self._hits, key=lambda x: x.get_best_hsp().get_e_value())

    def write(self, file_name, filter_hits=False):
        """Write compacted hit records to file, optionally filtering hits."""
        for hmmer_hit_obj in self.get_hits():
            for hmmer_hsp_obj in hmmer_hit_obj.get_hsps():
                try:
                    string = "%s\t%s\t%s\t-1\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self._query, self._query_length, hmmer_hit_obj._hit, hmmer_hsp_obj.get_identities(), hmmer_hsp_obj.get_similarities(), hmmer_hsp_obj.get_gaps(), hmmer_hsp_obj.get_e_value(), hmmer_hsp_obj._query_sequence, hmmer_hsp_obj._hit_sequence, hmmer_hsp_obj.get_alignment_blocks())
                    if filter_hits:
                        if hmmer_hsp_obj._inclusion == "?": continue
                    functions.write(file_name, string)
                except:
                    pass

class HmmerHit(object):
    """
    Store one HMMER hit and all associated HSP/domain segments.

    Object features:
        - Hit identifier/name as reported by HMMER (`_hit`).
        - List of parsed HSP/domain alignment segments (`_hsps`) represented
          by `HmmerHSP` objects.
        - Access helpers for sorted HSP retrieval and best-segment selection.
    """

    def __init__(self, hit):
        """
        Initialize one hit with its identifier string.

        Args:
            hit (str): HMMER hit identifier/description line content.

        Returns:
            None.
        """
        self._hit = hit
        self._hsps = []
        
    def get_hit(self):
        """Return the hit identifier."""
        return self._hit

    def add_hsp(self, hmmer_hsp_obj):
        """Append one parsed HSP/domain object to the hit."""
        self._hsps.append(hmmer_hsp_obj)

    def get_hsps(self):
        """Return HSPs sorted by E-value."""
        return sorted(self._hsps, key=lambda x: x.get_e_value())

    def get_best_hsp(self):
        """Return the lowest-E-value HSP for this hit."""
        for hmmer_hsp_obj in self.get_hsps():
            return hmmer_hsp_obj

class HmmerHSP(object):
    """
    Store one aligned HMMER domain segment (HSP-like record).

    Object features:
        - Domain-level statistics and coordinates (`_domain`, `_inclusion`,
          `_score`, `_bias`, `_e_value`, `_query_from`, `_query_to`,
          `_hit_from`, `_hit_to`).
        - Reconstructed aligned query/hit sequences (`_query_sequence`,
          `_hit_sequence`) assembled from HMMER text blocks.
        - Similarity counter for non-identity conservative matches
          (`_similarities`).
        - Helper methods that derive identities, similarities, gaps, and
          contiguous alignment-block coordinate mappings.
    """

    def __init__(self, domain, inclusion, score, bias, e_value, query_from, query_to, hit_from, hit_to):
        """
        Initialize one HMMER domain segment.

        Args:
            domain (int): Domain index within the hit.
            inclusion (str): Inclusion marker (`!` included, `?` marginal).
            score (float): Domain score reported by HMMER.
            bias (float): Domain bias term from HMMER output.
            e_value (float): Domain E-value.
            query_from (int): Query alignment start position.
            query_to (int): Query alignment end position.
            hit_from (int): Hit alignment start position.
            hit_to (int): Hit alignment end position.

        Returns:
            None.
        """
        self._domain = domain
        self._inclusion = inclusion
        self._score = score
        self._bias = bias
        self._e_value = e_value
        self._query_from = query_from
        self._query_to = query_to
        self._hit_from = hit_from
        self._hit_to = hit_to
        self._query_sequence = ""
        self._hit_sequence = ""
        self._similarities = 0

    def add_sequence(self, sequence):
        """Append sequence text to query or hit side depending on state."""
        if len(self._query_sequence) == len(self._hit_sequence):
            self.add_query_sequence(sequence)
        else:
            self.add_hit_sequence(sequence)
        
    def add_query_sequence(self, sequence):
        """Append residues to the aligned query sequence."""
        self._query_sequence += sequence

    def add_hit_sequence(self, sequence):
        """Append residues to the aligned hit sequence."""
        self._hit_sequence += sequence

    def add_similarities(self, similarities):
        """Accumulate HMMER similarity markers for this segment."""
        self._similarities += similarities

    def get_e_value(self):
        """Return the domain E-value."""
        return self._e_value

    def get_query_sequence(self):
        """Return the aligned query sequence string."""
        return self._query_sequence

    def get_hit_sequence(self):
        """Return the aligned hit sequence string."""
        return self._hit_sequence

    def get_identities(self):
        """Count exact residue identities across aligned positions."""
        identities = 0

        for i in range(len(self._query_sequence)):
            if self._query_sequence[i] == self._hit_sequence[i]:
                identities += 1

        return identities

    def get_similarities(self):
        """Return identities plus HMMER similarity-only matches."""
        return self.get_identities() + self._similarities

    def get_gaps(self):
        """Return the total number of gap characters in both alignments."""
        return len(re.findall("\-", self._query_sequence)) + len(re.findall("\-", self._hit_sequence))

    def get_alignment_blocks(self):
        """Build compact contiguous query:hit alignment block coordinates."""
        query_blocks = [[]]
        hit_blocks = [[]]
        query_positions = list(range(self._query_from, self._query_to + 1))
        hit_positions = list(range(self._hit_from, self._hit_to + 1))
        for i in range(len(self._query_sequence)):
            if self._query_sequence[i] in aminoacids1to3 and self._hit_sequence[i] in aminoacids1to3:
                query_blocks[-1].append(query_positions[0])
                hit_blocks[-1].append(hit_positions[0])
            else:
                if len(query_blocks[-1]) > 0:
                    query_blocks.append([])
                if len(hit_blocks[-1]) > 0:
                    hit_blocks.append([])
            if self._query_sequence[i] in aminoacids1to3:
                query_positions.pop(0)
            if self._hit_sequence[i] in aminoacids1to3:
                hit_positions.pop(0)

        return ";".join(["%s:%s,%s:%s" % (query_blocks[i][0], hit_blocks[i][0], query_blocks[i][-1], hit_blocks[i][-1]) for i in range(len(query_blocks))])

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for HMMER searches.

    How to run:
        python hmmer.py -d DATABASE_FASTA -i INPUT_HMM
            [--dummy DUMMY_DIR --filter -o OUTPUT_FILE]

    Example:
        python hmmer.py -d proteins.fa -i query.hmm --filter -o hits.tsv

    The parser configures:
        - Input HMM profile and target sequence database.
        - Optional filtering by inclusion threshold markers.
        - Output destination and temporary working directory.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options describing HMMER execution.
    """

    parser = optparse.OptionParser("python hmmer.py -d DATABASE_FILE -i INPUT_FILE [--dummy DUMMY_DIR -f -o OUTPUT_FILE]")

    parser.add_option("-d", action="store", type="string", dest="database_file", help="Database file (in FASTA format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-f", "--filter", default=False, action="store_true", dest="filter_hits", help="Filter hits under inclusion threshold (default = False)", metavar="{boolean}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in HMM format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.database_file == None or options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_pfam(database_file, input_file, dummy_dir="/tmp"):
    """
    Run `hmmscan` and collect matching PFAM family identifiers.

    Args:
        database_file (str): HMM database file.
        input_file (str): Query FASTA file.
        dummy_dir (str, optional): Reserved temporary directory argument.

    Returns:
        set: Family identifiers reported by `hmmscan`.
    """
    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        hmmer_path = os.path.join(src_path, config.get("Paths", "hmmer_path"))
        # Exec process #
        process = subprocess.check_output([os.path.join(hmmer_path, "hmmscan"), database_file, input_file], stderr=subprocess.STDOUT)
        # Get BLAST object #
        families=set()
        for line in process.decode().split("\n"):
            if line.startswith(">>"):
                families.add(line.split()[1])
    except:
        os.system("%s %s %s"%(os.path.join(hmmer_path, "hmmscan"), database_file, input_file))
        raise ValueError("Could not exec hmmscan for %s" % input_file)

    return families


def get_hmmer_obj(database_file, input_file, dummy_dir="/tmp"):
    """
    Run `hmmsearch` and parse its output into a `Hmmer` object.

    Args:
        database_file (str): Target sequence database in FASTA format.
        input_file (str): Query profile/input accepted by `hmmsearch`.
        dummy_dir (str, optional): Reserved temporary directory argument.

    Returns:
        Hmmer: Parsed HMMER result object.
    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        hmmer_path = os.path.join(src_path, config.get("Paths", "hmmer_path"))
        # Exec process #
        process = subprocess.check_output([os.path.join(hmmer_path, "hmmsearch"), input_file, database_file], stderr=subprocess.STDOUT)
        # Get BLAST object #
        hmmer_obj = Hmmer(process.decode().split("\n"))
    except:
        raise ValueError("Could not exec hmmsearch for %s" % input_file)

    return hmmer_obj

#-------------#
# Main        #
#-------------#

def main():
    """
    Run the command-line HMMER workflow for profile-based homolog search.

    Workflow:
        1. Parse runtime options and create the dummy directory if needed.
        2. Execute `hmmsearch` with the requested query/database pair.
        3. Optionally filter out hits below HMMER inclusion thresholds.
        4. Write compacted results to file or stdout.

    Returns:
        None. Parsed HMMER hits are serialized to disk or printed.
    """
    # Arguments & Options #
    options = parse_options()
    if not os.path.exists(os.path.abspath(options.dummy_dir)):
       os.makedirs(os.path.abspath(options.dummy_dir))

    # Get HMMER object #
    hmmer_obj = get_hmmer_obj(os.path.abspath(options.database_file), os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        # Write output #
        dummy_file = os.path.abspath(os.path.join(options.dummy_dir, str(os.getpid()) + ".txt"))
        hmmer_obj.write(dummy_file, options.filter_hits)
        shutil.copy(dummy_file, options.output_file)
        os.remove(dummy_file)
    else:
        hmmer_obj.write(options.output_file, options.filter_hits)


if __name__ == "__main__":
    main()
