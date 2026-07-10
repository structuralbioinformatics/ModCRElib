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
from SBILib.external.blast import blast_parser

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for running a BLASTP search.

    How to run:
        python blast.py -d DATABASE_FASTA -i QUERY_FASTA
            [--dummy DUMMY_DIR --filter -o OUTPUT_FILE]

    Example:
        python blast.py -d proteins.fa -i query.fa --filter -o hits.tsv

    The parser configures:
        - Query sequence and searchable protein database.
        - Optional twilight-zone filtering of BLAST hits.
        - Output destination and temporary working directory.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options describing input/output paths and
        optional filtering behavior.
    """

    parser = optparse.OptionParser("python blast.py -d DATABASE_FILE -i INPUT_FILE [--dummy DUMMY_DIR -f -o OUTPUT_FILE]")

    parser.add_option("-d", action="store", type="string", dest="database_file", help="Database file (in FASTA format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-f", "--filter", default=False, action="store_true", dest="filter_hits", help="Filter twilight zone hits (default = False)", metavar="{boolean}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in FASTA format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.database_file is None or options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_blast_obj(database_file, input_file, dummy_dir="/tmp"):
    """
    Run `blastp` and parse the XML output into a BLAST result object.

    Args:
        database_file (str): Target protein database in FASTA format.
        input_file (str): Query FASTA file.
        dummy_dir (str, optional): Directory used for temporary BLAST output.

    Returns:
        BlastOutput: Parsed BLAST result object from `blast_parser`.
    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        blast_path = os.path.join(src_path, config.get("Paths", "blast_path"))
        blast_out = os.path.join(dummy_dir, os.path.basename(input_file).strip()+str(os.getpid()) + ".txt")
        # Exec process #
        #print("%s -query %s -db %s -out %s -outfmt 5 > %s"%(os.path.join(blast_path, "blastp"),input_file,database_file,blast_out,blast_out+".log"))
        process = subprocess.check_output([os.path.join(blast_path, "blastp"), "-query", input_file, "-db", database_file, "-out", blast_out, "-outfmt", "5"], stderr=subprocess.STDOUT)
        # Get BLAST object #
        blast_obj = blast_parser.parse_blast(query_sequence=None, blast_output_file=blast_out, selfHit=True, hitIDformat="all")
        # Remove BLAST file #
        #os.remove(blast_out)
    except:
        raise ValueError("Could not exec blastp for %s" % input_file)

    return blast_obj

#-------------#
# Main        #
#-------------#

def main():
    """
    Run the command-line BLAST workflow for homolog discovery.

    Workflow:
        1. Parse runtime options and resolve absolute paths.
        2. Execute BLASTP against the selected sequence database.
        3. Optionally apply twilight-zone filtering thresholds.
        4. Write compacted hits to file or print them to stdout.

    Returns:
        None. Compacted BLAST output is written to disk or printed.
    """
    # Arguments & Options #
    options = parse_options()

    # Get BLAST object #
    blast_obj = get_blast_obj(os.path.abspath(options.database_file), os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Twilight zone #
    tz_parameter = 0
    tz_type = None
    if options.filter_hits:
        tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
        tz_type = config.get("Parameters", "twilight_zone_type")

    # Output #
    if options.output_file is not None:
        # Write output #
        dummy_file = os.path.abspath(os.path.join(options.dummy_dir, str(os.getpid()) + ".txt"))
        functions.write(dummy_file, blast_obj.str_compacted_blast(tz_parameter=tz_parameter, tz_type=tz_type))
        shutil.copy(dummy_file, options.output_file)
        os.remove(dummy_file)
    else:
        print((blast_obj.str_compacted_blast(tz_parameter=tz_parameter, tz_type=tz_type)))


if __name__ == "__main__":
    main()
