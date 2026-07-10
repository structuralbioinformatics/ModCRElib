"""
Run MEME/FIMO motif scanning and expose the results through small helper classes.

This module provides a lightweight wrapper around the FIMO command-line tool,
plus in-memory containers for parsed FIMO hits that are reused by several
scanning workflows in ModCRE.
"""

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


# Import my modules #
from ModCRElib.structure.contacts import triads


#-------------#
# Classes     #
#-------------#

class Fimo(object):
    """
    Container for the parsed output of one FIMO execution.

    The object stores the original text lines and exposes the detected motif
    hits as :class:`FimoHit` instances.

    Object features:
        - Raw FIMO text output preserved line-by-line (`_file_content`).
        - Parsed hit collection as `FimoHit` objects (`_hits`).
        - Query/sequence identifier extracted from FIMO lines (`_query`).
        - Parsing logic that supports alternative FIMO output layouts and
          normalizes reverse-strand hits to complementary DNA sequence context.
        - Accessors and filters for query metadata, hit presence checks, and
          iteration over stored hits.
        - Serialization helper to write the original FIMO output back to disk.
    """

    def __init__(self, file_content):
        """
        Create a parsed FIMO-result object.

        Args:
            file_content (list[str]): Lines of FIMO output.
        """
        self._file_content = file_content
        self._hits = []
        # Initialize #
        self._parse_file()

    def _parse_file(self):
        """
        Parse the raw FIMO output lines into :class:`FimoHit` objects.

        Returns:
            None. Populates ``self._hits`` and ``self._query`` in place.
        """

        #if config.get("Parameters", "fimo_pvalue_threshold") is not None:
        #   fimo_pvalue_threshold=float(config.get("Parameters", "fimo_pvalue_threshold"))
        #else:
        #   fimo_pvalue_threshold=1.0
        # We collect all p-values selected by fimo without exception
        fimo_pvalue_threshold=1.0
        for linefile in self._file_content:
          if linefile.startswith("#"): continue
          try:
             line = linefile.split("\t")
             #print("LENGTH LINE "+str(len(line))+ " line -3 "+line[-3] + " Threshold "+ str(fimo_pvalue_threshold)+" LINE " + str(line))
             if len(line) == 9:
                 read_condition=False
                 if line[-3] != "":
                    if float(line[-3]) < fimo_pvalue_threshold: read_condition=True
                 if read_condition:
                      strand=line[4]
                      if strand == "+": sequence=line[-1]
                      else: sequence = triads.get_complementary_dna_sequence(line[-1])
                      self._query = line[1]
                      self._hits.append(FimoHit(line[0], int(line[2]), int(line[3]), strand, float(line[5]), float(line[6]), sequence))
             if len(line) == 10:
                 read_condition=False
                 if line[-2] != "":
                    if float(line[-2]) < fimo_pvalue_threshold: read_condition=True
                 if line[-3] != "":
                    if float(line[-3]) < fimo_pvalue_threshold: read_condition=True
                 if read_condition:
                      strand=line[5]
                      if strand == "+": sequence=line[-1]
                      else: sequence = triads.get_complementary_dna_sequence(line[-1])
                      self._query = line[1]
                      self._hits.append(FimoHit(line[0], int(line[3]), int(line[4]), strand, float(line[-4]), float(line[-3]), sequence))
          except:
             sys.stderr.write("WARNING FIMO: skip %s\n"%(linefile))

    def get_query(self):
        """Return the query/sequence identifier associated with the FIMO run."""
        return self._query

    def has_hit(self, hit_name):
        """
        Check whether a motif hit with the given name is present.

        Args:
            hit_name (str): Motif identifier.

        Returns:
            bool: True if at least one hit matches ``hit_name``.
        """
        for hit_obj in self._hits:
            if hit_name == hit_obj.get_hit():
                return True

        return False

    def get_hits(self, sort=False):
        """
        Return the list of FIMO hits.

        Args:
            sort (bool): If True, sort hits by ascending p-value.

        Returns:
            list[FimoHit]: Parsed motif hits.
        """
        if sort:
            return sorted(self._hits, key=lambda x: x.get_p_value())

        return self._hits

    def get_hit(self, hit_name):
        """
        Return the first hit with the given motif name.

        Args:
            hit_name (str): Motif identifier.

        Returns:
            FimoHit | None: Matching hit or ``None`` if absent.
        """
        if self.has_hit(hit_name):
            for hit_obj in self._hits:
                if hit_name == hit_obj.get_hit():
                    return hit_obj

        return None

    def write(self, file_name):
        """
        Write the raw FIMO output back to disk.

        Args:
            file_name (str): Output path.
        """
        for line in self._file_content:
            functions.write(file_name, line)

class FimoHit(object):
    """
    Representation of one FIMO motif match on a DNA sequence.

    Object features:
        - Motif identifier reported by FIMO (`_hit`).
        - Match interval on the scanned sequence (`_start`, `_end`), using
          FIMO's coordinate convention from the parsed output.
        - Match orientation (`_strand`), typically `+` or `-`.
        - Match significance/intensity metrics (`_score`, `_p_value`).
        - Matched sequence string (`_sequence`); for reverse-strand hits this
          is already normalized by the parser to the complementary sequence
          context used by this codebase.
        - Lightweight accessor methods used by profile/scoring workflows to
          build per-position binding and score tracks.
    """

    def __init__(self, hit, start, end, strand, score, p_value, sequence):
        """
        Create one FIMO hit record.

        Args:
            hit (str): Motif identifier.
            start (int): Match start position.
            end (int): Match end position.
            strand (str): Match strand, usually ``+`` or ``-``.
            score (float): FIMO score.
            p_value (float): FIMO p-value.
            sequence (str): Matched DNA sequence.
        """
        self._hit = hit
        self._start = start
        self._end = end
        self._strand = strand
        self._score = score
        self._p_value = p_value
        self._sequence = sequence
        
    def get_hit(self):
            """Return the motif identifier for this hit."""
            return self._hit

    def get_start(self):
            """Return the start coordinate of this hit."""
            return self._start

    def get_strand(self):
            """Return the strand of this hit."""
            return self._strand

    def get_end(self):
            """Return the end coordinate of this hit."""
            return self._end

    def get_score(self):
            """Return the FIMO score for this hit."""
            return self._score

    def get_p_value(self):
            """Return the FIMO p-value for this hit."""
            return self._p_value

    def get_sequence(self):
            """Return the matched DNA sequence for this hit."""
            return self._sequence

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse the command line options for running FIMO.

    How to run:
        python fimo.py -d DATABASE_MEME -i INPUT_FASTA
            [--dummy DUMMY_DIR --ft P_VALUE --max MAX_MATCHES -o OUTPUT_FILE]

    Example:
        python fimo.py -d motifs.meme -i dna.fa --ft 1e-4 -o fimo_hits.txt

    The parser configures:
        - Motif database and FASTA input required by FIMO.
        - Optional p-value and match-cap controls.
        - Temporary/output path configuration.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options describing the motif database,
        input FASTA file, thresholds, and output path.
    """

    parser = optparse.OptionParser("python fimo.py -d database_file -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("-d", action="store", type="string", dest="database_file",  default=None, help="Database file (in MEME format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file",  default=None, help="Input file (in FASTA format)", metavar="{filename}")
    parser.add_option("--ft", action="store", type="float", dest="fimo_pvalue_threshold",  default=None, help="P-value threhold for fimo matches", metavar="{float}")
    parser.add_option("--max", action="store", type="int",dest="max_stored_matches",  default=None, help="Maximum number of matches stored", metavar="{integer}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None or options.database_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_fimo_obj(database_file, fasta_file,fimo_pvalue_threshold=None, max_stored_matches=None, dummy_dir="/tmp"):
    """
    Execute FIMO and return the parsed results.

    Args:
        database_file (str): MEME-format motif database.
        fasta_file (str): FASTA file to scan.
        fimo_pvalue_threshold (float | None): Optional p-value threshold.
        max_stored_matches (int | None): Optional FIMO storage cap.
        dummy_dir (str): Temporary directory used for FIMO output.

    Returns:
        Fimo: Parsed FIMO result object.

    Raises:
        ValueError: If FIMO execution fails.
    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        meme_path =  config.get("Paths", "meme_path")
        if not os.path.exists(config.get("Paths", "meme_path")):
           meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
        if fimo_pvalue_threshold is None:
           if config.get("Parameters", "fimo_pvalue_threshold") is not None:
              fimo_pvalue_threshold=float(config.get("Parameters", "fimo_pvalue_threshold"))
           else:
              fimo_pvalue_threshold=1
        if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
        # Exec process #
        log_file=os.path.join(dummy_dir,"fimo"+str(os.getpid())+".log")
        output_fimo=os.path.join(dummy_dir,"fimo"+str(os.getpid()))
        n=1
        while os.path.exists(output_fimo):
              output_fimo=output_fimo+str(n)
              n = n + 1

        try:
           fimo_file=os.path.join(output_fimo,"fimo.txt")
           if max_stored_matches is None:
              process = subprocess.check_output([os.path.join(meme_path, "fimo"), "-o",output_fimo, "--text --thresh", str(fimo_pvalue_threshold) , database_file, fasta_file], stderr=subprocess.STDOUT)
              #print("%s -o %s --text --thresh %s %s %s > %s\n"%(os.path.join(meme_path, "fimo"), output_fimo, str(fimo_pvalue_threshold),database_file, fasta_file,fimo_file))
           else:
              process = subprocess.check_output([os.path.join(meme_path, "fimo"), "-o",output_fimo, "--thresh", str(fimo_pvalue_threshold) ,"--max-stored-scores", str(max_stored_matches), database_file, fasta_file], stderr=subprocess.STDOUT)
           fimo_output=open(fimo_file,"r")
           # Get Fimo object #
           fimo_obj = Fimo(fimo_output.read().split("\n"))
           shutil.rmtree(output_fimo)
           return fimo_obj
        except:
         try:
           if not os.path.exists(output_fimo): os.makedirs(output_fimo)
           fimo_file=os.path.join(output_fimo,"fimo.txt")
           sys.stdout.write("\t-- execute system '--text' fimo option\n")
           if max_stored_matches is None:
              #print("%s -o %s --text --thresh %s %s %s > %s\n"%(os.path.join(meme_path, "fimo"), output_fimo, str(fimo_pvalue_threshold),database_file, fasta_file,fimo_file))
              os.system("%s -o %s --text --thresh %s %s %s > %s\n"%(os.path.join(meme_path, "fimo"), output_fimo, str(fimo_pvalue_threshold),database_file, fasta_file,fimo_file))
           else:
              #print("%s -o %s --text --thresh %s --max-stored-scores %s %s %s > %s\n"%(os.path.join(meme_path, "fimo"), output_fimo, str(fimo_pvalue_threshold),str(max_stored_matches),database_file, fasta_file,fimo_file))
              os.system("%s -o %s --text --thresh %s --max-stored-scores %s %s %s > %s\n"%(os.path.join(meme_path, "fimo"), output_fimo, str(fimo_pvalue_threshold),str(max_stored_matches),database_file, fasta_file,fimo_file))
           fimo_output=open(fimo_file,"r")
           # Get Fimo object #
           fimo_obj = Fimo(fimo_output.read().split("\n"))
           shutil.rmtree(output_fimo)
           return fimo_obj
         except:
           raise ValueError("Could not exec FIMO for %s" % fasta_file)
    except:
        raise ValueError("Could not exec FIMO for %s" % fasta_file)

def main():
    """
    Run the command-line FIMO wrapper.

    Workflow:
        1. Parse command-line options and resolve absolute paths.
        2. Execute FIMO and parse the resulting motif hits.
        3. Write parsed output to file or standard output.

    Args:
        None.

    Returns:
        None. Writes raw FIMO output to a file or stdout.
    """
    # Arguments & Options #
    options = parse_options()

    # Get FIMO object #
    database_file          =os.path.abspath(options.database_file)
    input_file             =os.path.abspath(options.input_file)
    if options.fimo_pvalue_threshold is not None: 
       fimo_pvalue_threshold  =float(options.fimo_pvalue_threshold)
    else:
       fimo_pvalue_threshold  = None
    if options.max_stored_matches is not None:
       max_stored_matches     =int(options.max_stored_matches)
    else:
       max_stored_matches     =None
    dummy_dir              =os.path.abspath(options.dummy_dir)
    if options.output_file is not None:
       output_file            =os.path.abspath(options.output_file)
    fimo_obj = get_fimo_obj(database_file, input_file, fimo_pvalue_threshold,max_stored_matches, dummy_dir)

    # Output #
    if options.output_file is not None:
        fimo_obj.write(output_file)
    else:
        for line in fimo_obj._file_content:
            sys.stdout.write("%s\n" % line)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
