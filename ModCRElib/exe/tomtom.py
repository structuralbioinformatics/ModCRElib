"""
Run TOMTOM motif comparisons and parse the resulting hit table.

This module provides lightweight wrappers around the MEME-suite ``tomtom``
tool, together with container classes that expose the parsed query and hit
information for downstream profile workflows.
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

class Tomtom(object):
    """
    Represent the parsed TOMTOM output for one query motif.

    The object stores the raw TOMTOM output lines, the inferred query motif
    identifier, and the ordered list of parsed ``TomtomHit`` matches.

    Object features:
        - Raw TOMTOM text output preserved line-by-line (`_file_content`).
        - Parsed query motif identifier (`_query`).
        - Ordered list of parsed hit objects (`_hits`), each represented as a
          `TomtomHit` instance with rank/statistical fields.
        - Cached hit count (`_size`) for fast reporting.
        - Helper methods to query existence, retrieve one hit by name, iterate
          all hits, and optionally sort hits by significance.
        - Serialization helper to write the original TOMTOM output back to
          disk for reproducible downstream processing.
    """

    def __init__(self, file_content):
        """
        Initialize a TOMTOM result object from raw output lines.

        Args:
            file_content (list[str]): TOMTOM text output split into lines.
                These lines are stored verbatim and parsed into `TomtomHit`
                objects by `_parse_file()`.

        Returns:
            None.
        """
        self._file_content = file_content
        self._hits = []
        self._query = None
        self._size  = 0
        # Initialize #
        self._parse_file()

    def _parse_file(self):
        """Parse the raw TOMTOM text output into ``TomtomHit`` objects."""
        rank=0
        for line in self._file_content:
            if line.startswith("Warning:"): continue
            if line.startswith("Query_ID"): continue
            if line.startswith("Processing"): continue
            if line.startswith("Estimat"): continue
            if line.startswith("Warning"): continue
            if line.startswith("#"): continue
            line = line.split("\t")
            if len(line) == 10:
                rank = rank + 1
                self._query = line.pop(0)
                hit, offset, p_value, e_value, q_value, overlap, query_sequence, hit_sequence, strand = line
                self._hits.append(TomtomHit(hit, offset, p_value, e_value, q_value, overlap, query_sequence, hit_sequence, strand, rank))
            else:
              if rank==0:
                sys.stderr.write("Could not execute TOMTOM first line %s\n"%(str(line)))
        self._size = rank

    def get_query(self):
        """Return the query motif identifier parsed from the TOMTOM output."""
        return self._query

    def get_size(self):
        """Return the number of parsed TOMTOM hits."""
        return self._size

    def has_hit(self, hit_name):
        """Check whether the result contains a hit with the given name."""
        for hit_obj in self._hits:
            if hit_name == hit_obj.get_hit():
                return True

        return False

    def get_hits(self, sort=False):
        """Return all parsed hits, optionally sorted by p-value."""
        if sort:
            return sorted(self._hits, key=lambda x: x.get_p_value())

        return self._hits

    def get_hit(self, hit_name):
        """Return one parsed hit by motif name, or ``None`` if absent."""
        if self.has_hit(hit_name):
            for hit_obj in self._hits:
                if hit_name == hit_obj.get_hit():
                    return hit_obj

        return None

    def write(self, file_name):
        """Write the original TOMTOM output lines to a file."""
        for line in self._file_content:
            functions.write(file_name, line)

class TomtomHit(object):
    """
    Store one TOMTOM motif-match result.

    Each hit contains the matched motif identifier plus the TOMTOM-reported
    alignment, significance, overlap, strand, and rank metadata.

    Object features:
        - Target motif identifier matched against the query (`_hit`).
        - Alignment offset between query and target motifs (`_offset`).
        - TOMTOM significance metrics (`_p_value`, `_e_value`, `_q_value`).
        - Number of aligned overlapping positions (`_overlap`).
        - Aligned query and target consensus sequence fragments
          (`_query_sequence`, `_hit_sequence`).
        - Strand/orientation annotation and rank order in the parsed result
          (`_strand`, `_rank`).
        - Accessor helpers used by profile-comparison and filtering workflows.
    """

    def __init__(self, hit, offset, p_value, e_value, q_value, overlap, query_sequence, hit_sequence, strand, rank):
        """
        Initialize one TOMTOM hit from parsed tabular fields.

        Args:
            hit (str): Identifier of the matched motif in the TOMTOM database.
            offset (int | str): Alignment offset reported by TOMTOM.
            p_value (float | str): TOMTOM p-value for this match.
            e_value (float | str): TOMTOM E-value for this match.
            q_value (float | str): TOMTOM q-value for this match.
            overlap (int | str): Number of overlapping aligned positions.
            query_sequence (str): Aligned query consensus sequence fragment.
            hit_sequence (str): Aligned target consensus sequence fragment.
            strand (str): Strand/orientation annotation for the match.
            rank (int | str): Rank order assigned while parsing output.

        Returns:
            None.
        """
        self._hit = hit
        self._offset = int(offset)
        try:
          self._p_value = float(p_value)
        except:
          self._p_value = p_value
        try:
          self._e_value = float(e_value)
        except:
          self._e_value = e_value
        try:
          self._q_value = float(q_value)
        except:
          self._q_value = q_value
        self._overlap = int(overlap)
        self._query_sequence = query_sequence
        self._hit_sequence = hit_sequence
        self._strand = strand
        self._rank   = int(rank)
        
    def get_hit(self):
            """Return the matched motif identifier."""
            return self._hit

    def get_offset(self):
            """Return the alignment offset reported by TOMTOM."""
            return self._offset

    def get_p_value(self):
           """Return the TOMTOM p-value for this hit."""
           return self._p_value

    def get_e_value(self):
            """Return the TOMTOM E-value for this hit."""
            return self._e_value

    def get_q_value(self):
           """Return the TOMTOM q-value for this hit."""
           return self._q_value

    def get_overlap(self):
            """Return the motif overlap length for this hit."""
            return self._overlap

    def get_query_sequence(self):
            """Return the aligned query consensus sequence."""
            return self._query_sequence

    def get_hit_sequence(self):
            """Return the aligned target consensus sequence."""
            return self._hit_sequence

    def get_strand(self):
            """Return the strand annotation for the motif match."""
            return self._strand

    def get_rank(self):
            """Return the rank assigned during TOMTOM parsing."""
            return self._rank

    def write(self,file_name):
            """Write a one-line summary of the hit to a file."""
            functions.write(file_name,"Hit %s Offset %d pValue %e eValue %e qValue %e Overlap %d Query sequence %s Hit sequence %s Strand %s Rank %d\n"%(self._hit,self._offset,self._p_value,self._e_value,self._q_value,self._overlap,self._query_sequence,self._hit_sequence,self._strand, self._rank))

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for TOMTOM execution.

    How to run:
        python tomtom.py -d DATABASE_MEME -i QUERY_MEME
            [--dummy DUMMY_DIR -o OUTPUT_FILE]

    Example:
        python tomtom.py -d motifs.meme -i query_pwm.meme -o tomtom_hits.txt

    The parser configures:
        - Query motif and target motif database.
        - Temporary directory used by TOMTOM fallbacks.
        - Optional output-file destination.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options with the query motif, motif
        database, dummy directory, and optional output file.

    Raises:
        SystemExit: Triggered by ``OptionParser`` when required arguments are
        missing or invalid.
    """

    parser = optparse.OptionParser("python tomtom.py -d database_file -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("-d", action="store", type="string", dest="database_file", help="Database file (in MEME format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (pwm in MEME format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_tomtom_obj(database_file, pwm_file, dummy_dir="/tmp"):
    """
    Execute TOMTOM for one PWM against a motif database and parse the result.

    Args:
        database_file (str): MEME-format motif database used as target.
        pwm_file (str): Query motif file in MEME format.
        dummy_dir (str): Temporary directory for fallback file-based execution.

    Returns:
        Tomtom: Parsed TOMTOM result object.

    Raises:
        ValueError: If TOMTOM execution fails or its output cannot be parsed.
    """
    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        meme_path =  config.get("Paths", "meme_path")
        if not os.path.exists(config.get("Paths", "meme_path")):
           meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
        # Always define output_tomtom upfront
        output_tomtom = os.path.join(dummy_dir, "tomtom" + str(os.getpid()))
        if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
        if not os.path.exists(output_tomtom): os.makedirs(output_tomtom)
        tomtom_file = os.path.join(output_tomtom, "tomtom.txt")
        try:
         # Attempt 1: subprocess directly (no shell redirection needed)
         tomtom_file=os.path.join(output_tomtom,"tomtom.txt")
         sys.stdout.write("\t\t\t-- %s -text -thresh 1.0 %s %s >& %s \n"%(os.path.join(meme_path, "tomtom"),pwm_file,database_file,tomtom_file))
         process = subprocess.check_output([os.path.join(meme_path, "tomtom"), "-text", "-thresh", "1.0", pwm_file, database_file], stderr=subprocess.STDOUT)
         # Get Fimo object #
         tomtom_obj = Tomtom(process.decode().split("\n"))
         return tomtom_obj
        except Exception as e:
         sys.stdout.write("\t\t\t-- subprocess failed (%s), trying os.system\n" % str(e))
         try:
           # Attempt 2: os.system with POSIX-compatible redirection (> file 2>&1)
           cmd = "%s -text -thresh 1.0 %s %s > %s 2>&1" % (
           os.path.join(meme_path, "tomtom"), pwm_file, database_file, tomtom_file)
           sys.stdout.write("\t\t\t-- %s\n" % cmd)
           os.system(cmd)
           #sys.stdout.write("\t\t\t-- %s -text -thresh 1.0 %s %s >& %s \n"%(os.path.join(meme_path, "tomtom"),pwm_file,database_file,tomtom_file))
           #os.system("%s -text -thresh 1.0 %s %s >& %s"%(os.path.join(meme_path, "tomtom"),pwm_file,database_file,tomtom_file))
           #tomtom_output=open(tomtom_file,"r")
           if not os.path.exists(tomtom_file):
                    raise FileNotFoundError("tomtom produced no output file: %s" % tomtom_file)
           with open(tomtom_file, "r") as tomtom_output:
                    tomtom_obj = Tomtom(tomtom_output.read().split("\n"))
           if 'output_tomtom' in locals() and os.path.exists(output_tomtom):
               shutil.rmtree(output_tomtom)

           return tomtom_obj
         except Exception as e:
              raise ValueError("Could not exec Tomtom for  %s, error %s" % (pwm_file,e))
    except Exception as e:
        raise ValueError("Could not exec Tomtom for %s, error %s" % (pwm_file,e))





#-------------#
# Main        #
#-------------#

def main():
    """
    Run the TOMTOM command-line workflow for one query motif.

    Workflow:
        1. Parse command-line arguments and resolve the selected paths.
        2. Execute TOMTOM against the requested motif database.
        3. Write the parsed TOMTOM output either to a file or to standard
           output.

    Args:
        None.

    Returns:
        None. Output is written to the chosen destination.
    """

    # Arguments & Options #
    options = parse_options()

    # Get Tomtom object #
    tomtom_obj = get_tomtom_obj(os.path.abspath(options.database_file), os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        tomtom_obj.write(os.path.abspath(options.output_file))
    else:
        for line in tomtom_obj._file_content:
            sys.stdout.write("%s\n" % line)


if __name__ == "__main__":
    main()
