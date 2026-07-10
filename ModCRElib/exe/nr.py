import os, sys, re
import configparser
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


# Import my modules #
from ModCRElib.structure.contacts import triads
from ModCRElib.structure.threading import threader

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    Parse command-line options for non-redundant triad selection.

    How to run:
        python nr.py --pdb PDB_DIR -r ROOT_DIR
            [--pbm PBM_DIR -i PDB_CHAIN -t THRESHOLD -o OUTPUT_FILE -v]

    Example:
        python nr.py --pdb pdb_data -r /path/to/project --pbm pbm_data
            -i general -t 0.35 -o nr_triads.list -v

    The parser configures:
        - PDB/PBM resource paths and root workspace directory.
        - Chain scope and optional TF inclusion/exclusion filters.
        - Redundancy mode/threshold and output destination.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options describing the NR workflow.
    """

    parser = optparse.OptionParser("python nr.py --pdb PDB_DIR -r ROOT_DIR [-f FILTER_TF_ID -i PDB_CHAIN -l LIST_FILE -n -o OUTPUT_FILE --pbm PBM_DIR -t THRESHOLD -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-f", action="store", default=None, type="string", dest="filter_tf_id", help="Filter TF id (this TF will not be considered). More than one TF or PDB_chain code can be included in a list separated by ',' (i.e. 2a07_G,3qrf_F,4wk8_G). Default is None and it uses the list of the configuration file dfined as 'tf_rejected' ", metavar="{string}")
    parser.add_option("-i", default="general", action="store", type="string", dest="pdb_chain", help="Input PDB chain (either \"general\" or a PDB chain as in \"1ahd_P\"; default=\"general\")", metavar="{string}")
    parser.add_option("-l", action="store", type="string", dest="list_file", help="List file of allowed TF ids (only these TF will be considered)", metavar="{string}")
    parser.add_option("-n", "--nr", default=False, action="store_true", dest="nr", help="High Non-redundant mode (default = False). If applied: all vs all TFs are compared; otherwise, it uses non-redudant PDB files and filters out similar KMERS without testing")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (output directory from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="{directory}")
    parser.add_option("-r", "--root", action="store", type="string", dest="root_dir", help="Root directory (i.e. protein-DNA directory)", metavar="{directory}")
    parser.add_option("-t", default=0.35, action="store", type="float", dest="threshold", help="Redundancy threshold (default = 0.35)", metavar="{float}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.pdb_dir is None or options.root_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def equivalent_kmers(kmer_a,kmer_b,min_redundant):
    """
    Compare two k-mers and test whether they are equivalent by overlap.

    The comparison allows shifts left/right and considers two k-mers
    equivalent when at least `min_redundant` nucleotide identities can be
    aligned.

    Args:
        kmer_a (str): First k-mer.
        kmer_b (str): Second k-mer.
        min_redundant (int): Minimum identity count required.

    Returns:
        bool: True when the k-mers are considered equivalent.
    """
    shift = max(0, (min(len(kmer_a), len(kmer_b)) - min_redundant))
    if len(kmer_a) >=  len(kmer_b):
       large = list(kmer_a)
       short = list(kmer_b)
    else:
       large = list(kmer_b)
       short = list(kmer_a)
    difference=len(large)-len(short)
    maxim_shift= difference+shift
    is_equivalent = False
    for shi in range(-maxim_shift,maxim_shift+1,1):
        equivalent=0
        for i in range(len(large)):
          pos = i + shi
          if pos in range(len(short)):
             if large[i] == short[pos]:
                equivalent += 1
        if equivalent >= min_redundant:
           is_equivalent = True
           break
    return  is_equivalent      

         
    

def get_nr_triads(pdb_dir, pdb_chain="general", threshold=0.50, pbm_dir=None, filter_tf_id=None, list_file=None, non_redundant=False):
    """
    Build a non-redundant list of triad files for potential derivation.

    Args:
        pdb_dir (str): PDB workspace directory generated by `pdb.py`.
        pdb_chain (str, optional): Chain scope (`general` or `PDB_CHAIN`).
        threshold (float, optional): Redundancy threshold for triad overlap.
        pbm_dir (str, optional): PBM workspace directory from `pbm.py`.
        filter_tf_id (str, optional): Comma-separated TF identifiers to skip.
        list_file (str, optional): File listing allowed TF identifiers.
        non_redundant (bool, optional): Use precomputed NR complexes to speed
            up processing.

    Returns:
        list: Paths to triad files selected as non-redundant.
    """

    # Initialize #
    done = set()
    allowed_tfs = set()

    # If list of allowed TFs... #
    if list_file is not None:
        # For each line... #
        for line in functions.parse_file(list_file):
            allowed_tfs.add(line)
    
    if pdb_chain=="general":
      folds_file = os.path.join(pdb_dir, "folds", pdb_chain + ".txt")
      if not os.path.exists(folds_file):
       fold_general={}
       for folds in os.listdir(os.path.join(pdb_dir, "folds")):
        for line in functions.parse_file(os.path.join(pdb_dir,"folds", folds.strip())): 
         if line.startswith("#"): continue
         pdb_tm_chain, tm_score = line.split(";")
         fold_general.setdefault(pdb_tm_chain,tm_score)
       fd=open(folds_file,"w")
       fd.write("#pdb_chain;tm-score\n")
       for pdb_tm_chain,tm_score in fold_general.items():
        fd.write("%s;%s\n" % (pdb_tm_chain, tm_score))
       fd.close()
    
    if pbm_dir is None:
        # Initialize #
        fold = []
        nr_pdb = []
        triads_objects = []
        folds_file = os.path.join(pdb_dir, "folds", pdb_chain + ".txt")
        # If folds file exists... #
        if os.path.exists(folds_file):
            # For each line... #
            for line in functions.parse_file(folds_file):
                if line.startswith("#"): continue
                pdb_chain, tm_score = line.strip().split(";")
                # Add PDB chain to fold
                fold.append(pdb_chain)
        # For each triads file... #
        for triads_file in sorted(os.listdir(os.path.join(pdb_dir, "triads"))):
            # Skip if PDB chain not in fold #
            if "index.txt" in triads_file:continue
            if len(fold) > 0 :
                m = re.search("(\S{4}\_\S).txt", triads_file)
                if not m: 
                  sys.stderr.write("\t\tSkip triad %s\n"%triads_file)
                  continue
                if m.group(1) not in fold: 
                  sys.stderr.write("\t\tSkip %s not found in pdb/fold of %s \n"% (m.group(1),pdb_chain))
                  continue
            sys.stdout.write("\t\t-- cluster similar folds of %s \n"%m.group(1))
            triads_objects.append(triads.Triads(os.path.join(pdb_dir, "triads", triads_file)))
        # For each triads object... #
        for triads_obj in sorted(triads_objects, key=lambda x: len(x.get_triads()), reverse=True):
            sys.stdout.write("\t\t\t-- check redundancy of %s\n"%(triads_obj._file))
            # Store and skip if threshold >= 1.0
            if threshold >= 1.0:
               nr_pdb.append(triads_obj._file)
               continue
            # Initialize #
            is_nr = True
            # For each nr triads object... #
            for nr_triads_file in nr_pdb:
                nr_triads_obj=triads.Triads(nr_triads_file)
                if triads_obj.get_percentage_common_triads(nr_triads_obj) > threshold:
                    is_nr = False
                    break
            if is_nr:
                nr_pdb.append(triads_obj._file)
        return nr_pdb
    else:
        # Initialize #
        fold = []
        nr_pbm = []
        nr_pdb = []
        threading_files = []
        triads_files = {}
        triads_dummy = {}
        folds_file = os.path.join(pdb_dir, "folds", pdb_chain + ".txt")
        # If folds file exists... #
        sys.stdout.write("\t\t-- cluster similar folds of %s \n"%pdb_chain)
        if os.path.exists(folds_file):
            # For each line... #
            for line in functions.parse_file(folds_file):
                if line.startswith("#"): continue
                pdb_chain, tm_score = line.split(";")
                # Add PDB chain to fold
                fold.append(pdb_chain)
        sys.stdout.write("\t\t\t-- Folds to collect : %s\n"%(",".join(fold)) )
        if not non_redundant:
            # For each line... #
            for line in functions.parse_file(os.path.join(pdb_dir, "nr", pdb_chain + ".txt")):
                if line.startswith("#"): continue
                m = re.search("\/(\S{4}\_\S).txt", line)
                if not m: 
                  sys.stderr.write("\t\tSkip file %s\n"%line)
                  continue
                nr_pdb.append(os.path.join(pdb_dir, "triads", m.group(1) + ".txt"))
        else:
            # For each line... #
            for line in functions.parse_file(os.path.join(pdb_dir, "nr", pdb_chain + ".txt")):
                if line.startswith("#"): continue
                m = re.search("\/(\S{4}\_\S).txt", line)
                if not m: 
                  sys.stderr.write("\t\tSkip file %s\n"%line)
                  continue
                nr_pbm.append(os.path.join(pdb_dir, "triads", m.group(1) + ".txt"))

        # For each threading file... #
        for threading_file in sorted(os.listdir(os.path.join(pbm_dir, "threading"))):
            if "index.txt" in threading_file: continue
            m = re.search("(\S+).\d+.(\S{4}\_\S).txt", threading_file)
            # Skip if PDB chain not in fold #
            if not m:
              sys.stderr.write("\t\tWrong threading file %s\n"%threading_file)
              continue
            if len(fold) > 0:
                if m.group(2) not in fold: continue
            # Skip if threading file belongs to filtered TF #
            if filter_tf_id is not None:
               list_tf_rejected = filter_tf_id.split(",")
               if len(list_tf_rejected) <= 0: list_tf_rejected=None
            else:
               list_tf_rejected = config.get("Parameters", "tf_rejected")
            if list_tf_rejected is not None:
               skip = False
               for tf_out in list_tf_rejected:
                   if m.group(1) == tf_out: skip = True
                   if m.group(2) == tf_out: skip = True
               if skip: continue
            # Skip if threading file belongs to non-allowed TF #
            if list_file is not None:
                if m.group(1) not in allowed_tfs: continue
            # For each line... #
            for line in functions.parse_file(os.path.join(pbm_dir, "threading", threading_file)):
                m = re.search("identity  : (.+)", line)
                if m:
                    sys.stdout.write("\t\t-- collect thread %s\n"%threading_file)
                    sys.stdout.flush()
                    threading_files.append((threading_file, float(m.group(1))))
                    break
        # For each triads file... #
        number_of_triads = len(os.listdir(os.path.join(pbm_dir, "triads")))
        number_check = 0
        for triads_file in sorted(os.listdir(os.path.join(pbm_dir, "triads"))):
            if "index.txt" in triads_file:continue
            number_check = number_check + 1
            sys.stdout.write("\t\t-- collect triad %s (%d out of %d)\n"%(triads_file,number_check,number_of_triads))
            sys.stdout.flush()
            m = re.search("(\S+).(\S{4}\_\S).([ACGT]{1,8}).\d+\-\d+.txt", triads_file)
            if not m: 
              sys.stderr.write("\t\tWrong triad file %s\n"%triads_file)
              continue
            # Skip if PDB chain not in fold #
            if len(fold) > 0:
                if m.group(2) not in fold: continue
            triads_dummy.setdefault(m.group(1), {})
            triads_dummy[m.group(1)].setdefault(m.group(2), {})
            triads_dummy[m.group(1)][m.group(2)].setdefault(m.group(3), triads_file)
        number_of_triads = len(triads_dummy.keys())
        number_check = 0
        for tf in triads_dummy.keys():
          triads_files.setdefault(tf,{})
          number_check = number_check + 1
          sys.stdout.write("\t\t-- check similar KMERS of %s (%d out of %d)\n"%(tf,number_check,number_of_triads))
          sys.stdout.flush()
          for pdb in triads_dummy[tf].keys():
            triads_files[tf].setdefault(pdb,{})
            kmer_set  = set([kmer for kmer in triads_dummy[tf][pdb].keys()])
            kmer_list = []
            for kmer in kmer_set:
                add = True
                for k in range(len(kmer_list)):
                    kmer_selected=kmer_list[k]
                    if kmer_selected in kmer:
                       kmer_list.pop(k)
                       kmer_list.append(kmer)
                       add=False
                    if kmer in kmer_selected:
                       add=False
                if add: kmer_list.append(kmer)
            kmer_selected_set=set(kmer_list)
            for kmer in kmer_selected_set:
                triads_files[tf][pdb].setdefault(kmer, triads_dummy[tf][pdb][kmer])          
        
            
        # Define the set of done... #
        done_tf_pdb={}
        # For each threading file... #
        number_of_threads = len(threading_files)
        number_check = 0
        for threading_file, percentage_identity in sorted(threading_files, key=lambda x: x[-1], reverse=True):
            number_check = number_check + 1
            sys.stdout.write("\t\t-- checking threading file %s (%d out of %d)\n"%(threading_file,number_check,number_of_threads))
            sys.stdout.flush()
            if non_redundant:
                m = re.search("T\d+_\d+.\d+.\d+.(\S{6}).txt", threading_file)
                if m:
                    pdb_chain = m.group(1)
                    nr_file = os.path.join(pbm_dir, "nr", pdb_chain + ".txt")
                    if pdb_chain in done: continue
                    # For each line... #
                    for line in functions.parse_file(os.path.join(nr_file)):
                        if os.path.basename(pdb_dir)+"/triads/" in line: continue
                        line = line.split("/")
                        if line[-1] in done: continue
                        # Define triads_obj
                        triads_obj = triads.Triads(os.path.join(pbm_dir, "triads", line[-1]))
                        # Store and skip if threshold >= 1.0
                        if threshold >= 1.0:
                           nr_pbm.append(triads_obj._file)
                           done.add(line[-1])
                           continue
                        # Initialize #
                        is_nr = True
                        # For each nr triads object... #
                        for nr_triads_file in nr_pbm:
                            nr_triads_obj=triads.Triads(nr_triads_file)
                            if triads_obj.get_percentage_common_triads(nr_triads_obj) > threshold:
                                is_nr = False
                                break
                        if is_nr:
                            nr_pbm.append(triads_obj._file)
                        done.add(line[-1])
                    done.add(pdb_chain)
                    sys.stdout.write("\t\t\t-- done with %s (all done %d non-redundant %d)\n"%(pdb_chain,len(done), len(nr_pbm))) 
            else:
                # Initialize #
                kmers = []
                m = re.search("(\S+).(\S{4}\_\S).txt", threading_file)
                min_redundant_kmer = int(config.get("Parameters","min_redundant_kmer"))
                # For each line... #
                if not m:
                    sys.stderr.write("\t\tWrong threading file %s\n"%threading_file)
                    continue
                thread_obj=threader.Threaded(os.path.join(pbm_dir, "threading", threading_file))
                kmers=[kmer for kmer in thread_obj.get_kmers().keys()]
                tf_pdb=(m.group(1),m.group(2))
                # For each k-mer... #
                for kmer in kmers:
                 for kmer_interface in triads_files[m.group(1)][m.group(2)].keys():
                    # Initialize #
                    is_nr = True
                    if kmer_interface not in kmer: continue
                    sys.stdout.write("\t\t\t-- analyze kmer %s of triad %s in %s from %s\n"%(kmer_interface,triads_files[m.group(1)][m.group(2)][kmer_interface],kmer,threading_file))
                    try:
                      triads_obj = triads.Triads(os.path.join(pbm_dir, "triads", triads_files[m.group(1)][m.group(2)][kmer_interface]))
                    except:
                      sys.stdout.write("Failed to open Triads %s %s %s  for Threading %s\n"%(m.group(1),m.group(2),kmer_interface,threading_file))
                      exit(0)
                    # Skip if no triads #
                    if len(triads_obj.get_triads()) == 0: continue
                    # Store and skip if threshold >= 1.0
                    if threshold >= 1.0:
                        nr_pbm.append(triads_obj._file)
                        continue
                    else:
                      # For each nr triads object... #
                      for nr_triads_file in nr_pdb:
                        nr_triads_obj= triads.Triads(nr_triads_file)
                        if triads_obj.get_percentage_common_triads(nr_triads_obj) > threshold:
                            is_nr = False
                            sys.stdout.write("\t\t\t\t-- skip: found pdb triad %s too similar with triad %s \n"%(nr_triads_file,triads_files[m.group(1)][m.group(2)][kmer_interface]))
                            sys.stdout.flush()
                            break
                      if is_nr:
                         if tf_pdb not in done_tf_pdb:
                            done_tf_pdb.setdefault(tf_pdb,set()).add(kmer_interface)
                         else:
                            for done_kmer in done_tf_pdb.get(tf_pdb):
                                is_nr = not equivalent_kmers(done_kmer,kmer_interface,min_redundant_kmer)
                                if not is_nr:
                                   sys.stdout.write("\t\t\t\t-- skip (redundant with kmer %s): kmer %s of triad %s in %s from %s\n"%(done_kmer,kmer_interface,triads_files[m.group(1)][m.group(2)][kmer_interface],kmer,threading_file))
                                   sys.stdout.flush()
                                   break
                            if is_nr:
                                done_tf_pdb.setdefault(tf_pdb,set()).add(kmer_interface)
                      if is_nr:
                        # For each nr triads object... #
                        if len(nr_pbm)>0:
                          for nr_triads_file in nr_pbm:
                            #sys.stdout.write("\t\t\t\t-- test %s common with %s \n"%(nr_triads_file,triads_obj._file))
                            nr_triads_obj= triads.Triads(nr_triads_file)
                            #sys.stdout.write("\t\t\t\t-- test %s common is %f > %f ?\n"%(nr_triads_obj._file,triads_obj.get_percentage_common_triads(nr_triads_obj),threshold))
                            if triads_obj.get_percentage_common_triads(nr_triads_obj) > threshold:
                                is_nr = False
                                break
                      if is_nr:
                        nr_pbm.append(triads_obj._file)

        return nr_pdb + nr_pbm


#-------------#
# Main        #
#-------------#

def main():
    """
    Run the command-line workflow that builds non-redundant triad sets.

    Workflow:
        1. Parse runtime options and move into the requested root directory.
        2. Optionally skip execution when the destination file already exists.
        3. Collect non-redundant triad files from PDB or PDB+PBM sources.
        4. Write the selected triad-file list to the requested output target.

    Returns:
        None. Selected triad paths are written to the output stream/file.
    """
    # Arguments & Options #
    options = parse_options()

    # Change directory #
    os.chdir(options.root_dir)

    # Exit if output file already exists #
    if options.output_file is not None:
        if os.path.exists(options.output_file): exit(0)

    # Get non-redundant triads #
    if options.verbose:
       print("PDB_CHAIN     ",options.pdb_chain)
       print("PDB FOLDER    ",options.pdb_dir) 
       print("PBM FOLDER    ",options.pbm_dir) 
       print("LIST TFS      ",options.list_file)
       print("NON-REDUNDANT ",options.nr)

    nr_triads_files = get_nr_triads(options.pdb_dir, options.pdb_chain, float(options.threshold), options.pbm_dir, options.filter_tf_id, options.list_file, options.nr)

    # For each nr triads object... #
    for nr_triads_file in nr_triads_files:
        functions.write(options.output_file, nr_triads_file)


if __name__ == "__main__":
    main()


