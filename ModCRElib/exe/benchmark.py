import os, sys, re
import configparser
import optparse
import matplotlib.patches as mpatches
import matplotlib.pyplot as plot
import numpy
import pandas
import seaborn
from scipy.stats import wilcoxon, rankdata
import socket
import subprocess
from sklearn import metrics


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
from ModCRElib.sequence import nr
from ModCRElib.potentials import spotentials
from ModCRElib.structure.contacts import triads
from ModCRElib.structure.dna import x3dna


#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse CLI options for benchmark data generation and evaluation.

    How to run:
        ``python benchmark.py --pbm pbm_dir --pdb pdb_dir [options]``

    Returns:
        optparse.Values: Namespace with benchmark inputs, fold controls, and
        execution flags.
    """

    parser = optparse.OptionParser("python benchmark.py --pbm=pbm_dir --pdb=pdb_dir [--dummy=dummy_dir -f folds -o output_dir -r randoms --start=start_step --stop=stop_step -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-f", default=5, action="store", type="int", dest="folds", help="Number of folds to split the data (default = 5)", metavar="{int}")
    parser.add_option("-n", default=10, action="store", type="int", dest="negatives", help="Number of times to randomly select negatives (for statistical purposes; default = 10)", metavar="{int}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (output directory from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="{directory}")
    parser.add_option("-r", default=100, action="store", type="int", dest="randoms", help="Number of randomly selected negative k-mers (default = 100)", metavar="{int}")    
    parser.add_option("--start", default=1, action="store", type="int", dest="start_step", help="Start at a given step (default = 1; first)", metavar="{int}")
    parser.add_option("--stop", default=7, action="store", type="int", dest="stop_step", help="Stop at a given step (default = 7; last)", metavar="{int}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("-j","--parallel", default=False, action="store_true", dest="parallel", help="Submit JOBS to Queues in parallel. The program stops, you need to re-start with the requested value (default = False)", metavar="{boolean}")
    parser.add_option("-g","--GC-content" , default=False, action="store_true", dest="gc_content", help="Negatives should have similar GC content than positives (Default is not )", metavar="{boolean}")
    parser.add_option("--nr", default=False,action="store_true", dest="non_redundant", help="High Non-redundant mode (default = False). If applied: all vs all TFs are compared; otherwise, it uses non-redudant PDB files and filters out similar KMERS without testing", metavar="{boolean}")
    parser.add_option("--tfs", action="store", type="string", dest="tfs_file", help="TFs file (from CIS-BP; i.e. cisbp_1.02.tfs.sql)", metavar="{filename}")
    parser.add_option("--family", action="store", type="string", dest="families_file", help="Families file (from CIS-BP; i.e. cisbp_1.02.tf_families.sql)", metavar="{filename}")
    parser.add_option("--motifs", action="store", type="string", dest="motifs_file", help="Motifs file (from CIS-BP; i.e. cisbp_1.02.motifs.sql)", metavar="{filename}")
   
 

    (options, args) = parser.parse_args()

    if options.pbm_dir is None or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def split_list_into_n_sublists_old(list_of_items, n=5):
    """
    Split a list of items into ``n`` approximately balanced sublists.

    Args:
        list_of_items (list): Input items.
        n (int): Number of output sublists.

    Returns:
        list[list]: Partitioned sublists.
    """

    # Initialize #
    list_of_lists = []
    # For each sublist... #
    for i in range(n):
        list_of_lists.append([])
    # For each item... #
    for item in list_of_items:
        # For each sublist... #
        for i in range(n):
            # Initialize #
            add = True
            # For each remaining sublist... #
            for j in range(i, n):
                if len(list_of_lists[i]) != len(list_of_lists[j]):
                    add = False
                    break
            # If add... #
            if add:
                list_of_lists[i].append(item)
                break

    return list_of_lists


def split_list_into_n_sublists(list_of_items, n=5):
    """
    Shuffle items and distribute them across ``n`` sublists.

    Args:
        list_of_items (list): Input items.
        n (int): Number of output sublists.

    Returns:
        list[list]: Shuffled and balanced sublists.
    """
    from random import shuffle
    # Initialize #
    list_of_lists = []
    # For each sublist... #
    for i in range(n):
        list_of_lists.append([])
    if len(list_of_items) <=0: return list_of_lists
    shuffled = list_of_items
    shuffle(shuffled)
    i=0
    for j in range(len(shuffled)):
        list_of_lists[i].append(shuffled[j])
        i = i+1
        if i == n: i=0
    return list_of_lists

def get_data_frame(file_name):
    """
    Load parsed benchmark rows into a typed pandas DataFrame.

    Args:
        file_name (str): Input benchmark table path.

    Returns:
        pandas.DataFrame: DataFrame with columns ``fold``, ``family``, ``score``,
        and ``class``.
    """

    # Initialize #
    data = []
    families_filter = set(config.get("Parameters", "families_filter").split(","))
    
    # For each line... #
    for line in functions.parse_file(file_name):
        if line.startswith("#"): continue
        line = line.split(";")
        if line[7] in families_filter: continue
        data.append(line[6:])
        data[-1][0] = int(data[-1][0])
        data[-1][2] = float(data[-1][2])
        data[-1][3] = int(data[-1][3])

    # Convert list of lists to 2D dataframe #
    data_frame = pandas.DataFrame(data)
    data_frame.columns = ["fold", "family", "score", "class"]

    return data_frame

def get_aucpr_data(data_frame, families, folds, file_name):
    """
    Compute and write AUC-PR values for each family/fold.

    Args:
        data_frame (pandas.DataFrame): Benchmark observations.
        families (list[str]): Families to evaluate.
        folds (int): Number of cross-validation folds.
        file_name (str): Output file path.

    Returns:
        None.
    """
    
    # Write output #
    functions.write(file_name, "#family;%s" % ";".join(map(str, [i + 1 for i in range(folds)])))
    # Initialize #
    aucpr = []
    # For each family... #
    for family in families:
        # Initialize #
        aucpr = []
        df_family = data_frame[(data_frame["family"] == family)]
        # For each fold... #
        for i in range(folds):
            # Initialize #
            df_fold = df_family[(df_family["fold"] == i + 1)]
            if len(df_fold) > 0:
                aucpr.append(get_aucpr(df_fold["score"].tolist(), df_fold["class"].tolist()))
            else:
                aucpr.append(None)
        functions.write(file_name, "%s;%s" % (family, ";".join(map(str, aucpr))))

def get_auroc_data(data_frame, families, folds, file_name):
    """
    Compute and write AUROC values for each family/fold.

    Args:
        data_frame (pandas.DataFrame): Benchmark observations.
        families (list[str]): Families to evaluate.
        folds (int): Number of cross-validation folds.
        file_name (str): Output file path.

    Returns:
        None.
    """
    
    # Write output #
    functions.write(file_name, "#family;%s" % ";".join(map(str, [i + 1 for i in range(folds)])))
    # Initialize #
    auroc = []
    # For each family... #
    for family in families:
        # Initialize #
        auroc = []
        df_family = data_frame[(data_frame["family"] == family)]
        # For each fold... #
        for i in range(folds):
            # Initialize #
            df_fold = df_family[(df_family["fold"] == i + 1)]
            if len(df_fold) > 0:
                auroc.append(get_auroc(df_fold["score"].tolist(), df_fold["class"].tolist()))
            else:
                auroc.append(None)
        functions.write(file_name, "%s;%s" % (family, ";".join(map(str, auroc))))


def get_auroc(scores, labels):
   """
   Compute AUROC from prediction scores and binary labels.

   Args:
       scores (list[float]): Prediction scores.
       labels (list[int]): Binary ground-truth labels.

   Returns:
       float | None: AUROC value, or ``None`` on computation failure.
   """
   pred = numpy.array(scores)
   y    = numpy.array(labels)
   auroc = None
   try:
     fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=1)
     auroc = metrics.auc(fpr, tpr)
   except Exception as e:
     print(("Failed %s"%e))
 
   return auroc

def get_aucpr(scores, labels):
    """
    Compute AUC-PR via trapezoidal integration on interpolated points.

    Args:
        scores (list[float]): Prediction scores.
        labels (list[int]): Binary ground-truth labels.

    Returns:
        float: Area under precision-recall curve.
    """

    # Initialize #
    TPA = 0
    TPB = 0
    FPA = 0
    FPB = 0
    points = []
    TP_dict = {}    
    paired_list = [x for x in zip(scores, labels)]
    paired_list.sort(key=lambda x: x[0], reverse=True)
    total_positives = sum(labels)

    for cutoff, label in paired_list:
        TP_dict.setdefault(cutoff, [0,0])
        if label:
            TP_dict[cutoff][0] += 1
        else:
            TP_dict[cutoff][1] += 1

    sorted_cutoffs = sorted(TP_dict.keys(), reverse=True)

    TPB = TP_dict[sorted_cutoffs[0]][0]
    FPB = TP_dict[sorted_cutoffs[0]][1]

    # Initialize #
    points.extend(interpolate(0, TPB, 0, FPB, total_positives))

    for cutoff in range(1, len(sorted_cutoffs)):
        TPA += TP_dict[sorted_cutoffs[cutoff - 1]][0]
        TPB = TPA + TP_dict[sorted_cutoffs[cutoff]][0]
        FPA += TP_dict[sorted_cutoffs[cutoff - 1]][1]
        FPB = FPA + TP_dict[sorted_cutoffs[cutoff]][1]
        p = interpolate(TPA, TPB, FPA, FPB, total_positives)
        points.extend(p)

    x, y = [x for x in zip(*points)]

    return numpy.trapz(x=x, y=y)

def interpolate(TPA, TPB, FPA, FPB, total_positives):
    """
    Interpolate precision-recall points between consecutive confusion states.

    Args:
        TPA (int | float): True positives at point A.
        TPB (int | float): True positives at point B.
        FPA (int | float): False positives at point A.
        FPB (int | float): False positives at point B.
        total_positives (int | float): Total positives in the dataset.

    Returns:
        list[tuple[float, float]]: Interpolated ``(recall, precision)`` points.
    """

    # Initialize #
    points = []
    TPA = float(TPA)
    TPB = float(TPB)
    FPA = float(FPA)
    FPB = float(FPB)

    if (TPA - TPB) != 0:
        skew = (FPB-FPA)/(TPB-TPA)
        for x in range(int(TPB) - int(TPA) + 1):
            if (TPA + x + FPA + skew * x) > 0:
                points.append(((TPA + x) / total_positives, (TPA + x) / (TPA + x + FPA + skew * x)))

    return points

def get_precision_recall_fpr_data(data_frame, families, scores, folds, ratio, file_name_a, file_name_b, file_name_c):
    """
    Write precision, recall, and false-positive-rate tables by threshold.

    Args:
        data_frame (pandas.DataFrame): Benchmark observations.
        families (list[str]): Families to evaluate.
        scores (list[float]): Score thresholds.
        folds (int): Number of cross-validation folds.
        ratio (float): Reserved compatibility parameter.
        file_name_a (str): Precision output table.
        file_name_b (str): Recall output table.
        file_name_c (str): FPR output table.

    Returns:
        None.
    """

    # Initialize #
    total_positives = {}
    total_negatives = {}
    # For each family... #
    for family in families:
        # Initialize #
        df_family = data_frame[(data_frame["family"] == family)]
        # For each fold... #
        for i in range(folds):
            total_positives.setdefault((family, i + 1), len(df_family[(df_family["fold"] == i + 1) & (df_family["class"] == 1)]))
            total_negatives.setdefault((family, i + 1), len(df_family[(df_family["fold"] == i + 1) & (df_family["class"] == 0)]))

    # Write output #
    functions.write(file_name_a, "#family;score;%s" % ";".join(map(str, [i + 1 for i in range(folds)])))
    functions.write(file_name_b, "#family;score;%s" % ";".join(map(str, [i + 1 for i in range(folds)])))
    functions.write(file_name_c, "#family;score;%s" % ";".join(map(str, [i + 1 for i in range(folds)])))

    # For each score... #
    for score in scores:
        # Initialize #
        df_scores = data_frame[(data_frame["score"] >= score)]
        # For each family... #
        for family in families:
            # Initialize #
            precision = []
            recall = []
            fpr = []
            df_family = df_scores[(df_scores["family"] == family)]
            # For each fold... #
            for i in range(folds):
                # Initialize #
                df_fold = df_family[(df_family["fold"] == i + 1)]
                fold = len(df_fold)
                precision.append(None)
                recall.append(None)
                fpr.append(None)
                if fold > 0:
                    df_positives = df_fold[(df_fold["class"] == 1)]
                    positives = float(len(df_positives))
                    negatives = fold - positives
                    precision[-1] = (positives / (positives + negatives)) * 100
                    if total_positives[(family, i + 1)]>0: recall[-1] = (positives / total_positives[(family, i + 1)]) * 100
                    if total_negatives[(family, i + 1)]>0: fpr[-1] = (negatives / total_negatives[(family, i + 1)]) * 100
            functions.write(file_name_a, "%s;%s;%s" % (family, score, ";".join(map(str, precision))))
            functions.write(file_name_b, "%s;%s;%s" % (family, score, ";".join(map(str, recall))))
            functions.write(file_name_c, "%s;%s;%s" % (family, score, ";".join(map(str, fpr))))



def get_precision_recall_data(data_frame, families, scores, folds, ratio, file_name_a, file_name_b):
    """
    Write precision and recall tables by score threshold.

    Args:
        data_frame (pandas.DataFrame): Benchmark observations.
        families (list[str]): Families to evaluate.
        scores (list[float]): Score thresholds.
        folds (int): Number of cross-validation folds.
        ratio (float): Reserved compatibility parameter.
        file_name_a (str): Precision output table.
        file_name_b (str): Recall output table.

    Returns:
        None.
    """

    # Initialize #
    total_positives = {}
    # For each family... #
    for family in families:
        # Initialize #
        df_family = data_frame[(data_frame["family"] == family)]
        # For each fold... #
        for i in range(folds):
            total_positives.setdefault((family, i + 1), len(df_family[(df_family["fold"] == i + 1) & (df_family["class"] == 1)]))

    # Write output #
    functions.write(file_name_a, "#family;score;%s" % ";".join(map(str, [i + 1 for i in range(folds)])))
    functions.write(file_name_b, "#family;score;%s" % ";".join(map(str, [i + 1 for i in range(folds)])))

    # For each score... #
    for score in scores:
        # Initialize #
        df_scores = data_frame[(data_frame["score"] >= score)]
        # For each family... #
        for family in families:
            # Initialize #
            precision = []
            recall = []
            df_family = df_scores[(df_scores["family"] == family)]
            # For each fold... #
            for i in range(folds):
                # Initialize #
                df_fold = df_family[(df_family["fold"] == i + 1)]
                fold = len(df_fold)
                precision.append(None)
                recall.append(None)
                if fold > 0:
                    df_positives = df_fold[(df_fold["class"] == 1)]
                    positives = float(len(df_positives))
                    negatives = fold - positives
                    precision[-1] = (positives / (positives + negatives)) * 100
                    recall[-1] = (positives / total_positives[(family, i + 1)]) * 100
            functions.write(file_name_a, "%s;%s;%s" % (family, score, ";".join(map(str, precision))))
            functions.write(file_name_b, "%s;%s;%s" % (family, score, ";".join(map(str, recall))))

def get_auprc_values(file_name, family=None):
    """
    Yield numeric AUPRC values from a benchmark result file.

    Args:
        file_name (str): Input table path.
        family (str | None): Optional family filter.

    Yields:
        float: Parsed AUPRC values.
    """

    # For each line... #
    for line in functions.parse_file(file_name):
        if line.startswith("#"): continue
        line = line.split(";")
        if family is not None:
            if family != line.pop(0): continue
        for value in line:
            try: yield float(value)
            except: pass

def one_sided_wilcoxon_test(a, b, bonferroni_correction=1):
    """
    Run one-sided Wilcoxon signed-rank test and return corrected significance.

    Args:
        a (array-like): First paired sample.
        b (array-like): Second paired sample.
        bonferroni_correction (int | float): Multiplicative correction factor.

    Returns:
        tuple[float, float, bool, float]: Median(a), median(b), ``b>a`` flag, and
        corrected one-sided p-value.
    """

    try:
      statistic, p_value = wilcoxon(a, b)
      a = numpy.median(a)
      b = numpy.median(b)
    except:
      print("Failed one sided wilcoxon test")
      p_value=1.0
      a = 0.0
      b = 0.0

    return a, b, b > a, (p_value / 2) * bonferroni_correction

def get_precision_recall_score(data_frame, thd_precision=None, thd_recall=None, thd_score=None):
    """
    Select best precision/recall operating point under optional constraints.

    Args:
        data_frame (pandas.DataFrame): Summary statistics by score.
        thd_precision (float | None): Optional minimum precision constraint.
        thd_recall (float | None): Reserved parameter.
        thd_score (float | None): Optional maximum score to inspect.

    Returns:
        tuple[float, float | None, float | None]: Selected precision, recall, score.
    """

    # Initialize #
    precision = 0
    recall = None
    score = None

    for i in reversed(numpy.arange(0, 1 + 0.01, 0.01)):
        if thd_score is None:
            df = data_frame[(data_frame["score"] == i) & (data_frame["statistic"] == "Precision")]
            mean = numpy.mean(df["value"])
            if mean > precision:
                precision = mean
                df = data_frame[(data_frame["score"] == i) & (data_frame["statistic"] == "Recall")]
                recall = numpy.mean(df["value"])
                score = i
        else:
            if thd_precision is None: 
                df = data_frame[(data_frame["score"] == i) & (data_frame["statistic"] == "Precision")]
                precision = numpy.mean(df["value"])
                df = data_frame[(data_frame["score"] == i) & (data_frame["statistic"] == "Recall")]
                recall = numpy.mean(df["value"])
                score = i
            else:
                df = data_frame[(data_frame["score"] == i) & (data_frame["statistic"] == "Precision")]
                mean = numpy.mean(df["value"])
                if mean > thd_precision:
                    precision = mean
                    df = data_frame[(data_frame["score"] == i) & (data_frame["statistic"] == "Recall")]
                    recall = numpy.mean(df["value"])
                    score = i
                elif mean > precision:
                    precision = mean
                    df = data_frame[(data_frame["score"] == i) & (data_frame["statistic"] == "Recall")]
                    recall = numpy.mean(df["value"])
                    score = i
            if i == thd_score:
                break

    return precision, recall, score

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()


    # Create output "main" subdirs #
    for subdir in ["nr", "folds", "potentials", "results", "parsed", "figures"]:
        if not os.path.exists(os.path.join(options.output_dir, subdir)):
            os.makedirs(os.path.join(options.output_dir, subdir))
    if not os.path.exists(os.path.join(options.dummy_dir)): os.makedirs(os.path.join(options.dummy_dir))

    # Verbose mode #
    if options.verbose: sys.stdout.write("\n")
    non_redundant = options.non_redundant

    ################################
    # 1. Create non-redundant sets #
    ################################
    if options.verbose: sys.stdout.write("Create non-redundant sets...\n\n")
    # Skip if starts later #
    if options.start_step <= 2:
      nr_file = os.path.join(options.output_dir, "nr", "family.txt")
      # Skip if nr file already exists or some redundancies are accpeted#
      if not os.path.exists(nr_file) :
        tfs_in_families = {}
        ##############################
        # 1.1 Parse motifs info      #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing motif information...\n")
        sys.stdout.flush()
        # Initialize #
        motifs = {}
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.motifs_file)):
            m = re.search("\('(.+)', '(.+)', '(.+)', '.+', 'PBM', '(.+)', '.+', '.+'\),*", line)
            if m:
                # Skip if motif file does not exist #
                motif_name=m.group(1)
                if not os.path.exists(os.path.join(options.pbm_dir, "motifs", motif_name + ".txt")): continue
                # Skip if not characterized protein sequence #
                if m.group(4) == "NULL": continue
                motifs.setdefault(m.group(1), [m.group(2), m.group(3), re.sub("[^A-Z]", "X", m.group(4).upper())])
        ##############################
        # 1.2 Parse TFs families     #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing TF families...\n")
        sys.stdout.flush()
        # Initialize #
        families = {}
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.families_file)):
            #m = re.search("\('(.+)', '(.+)', '.+', \d+, .+\),*", line)
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+'\),*", line)
            if m:
                families.setdefault(m.group(1), set(m.group(2).split(",")))
        ##############################
        # 1.3 Parse TFs info         #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing TF information...\n")
        sys.stdout.flush()
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '(.+)', '(.+)', '[DIN]'\),*", line)
            if m:
                # Initialize #
                tf_motifs = []
                # For each motif... #
                for motif in motifs:
                    if motifs[motif][0] == m.group(1):
                        tf_motifs.append(motif)
                if len(tf_motifs) > 0:
                    tf_id     =  m.group(1)
                    tf_family =  m.group(2)
                    #for i in range(len(tf_motifs)):
                    #    for fam in families[tf_family]:
                    #         tfs_in_families.setdefault(fam.replace("/",":").replace(" ","_"),set()).add(tf_id+"."+str(i))
                    for fam in families[tf_family]:
                        tfs_in_families.setdefault(fam.replace("/",":").replace(" ","_"),set()).add(tf_id)
        ##############################
        # 1.4 TF lists by family     #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parse non-redundant triads by families...\n")
        sys.stdout.flush()
        # For each family....
        nr_triads=set()
        if config.get("Parameters", "tf_rejected") is not None:
            list_tf_rejected = config.get("Parameters", "tf_rejected").split(",")
        else:
            list_tf_rejected = None
        for tf_family in tfs_in_families:
            list_of_tfs    = [x for x in tfs_in_families[tf_family]]
            #print("LIST",list_of_tfs)
            list_file      = os.path.join(options.output_dir, "nr", "TFs_family_"+tf_family+".txt")
            folds_file     = os.path.join(options.pdb_dir, "folds", "family_"+tf_family+".txt")
            nr_pdb_file    = os.path.join(options.pdb_dir, "nr", "family_"+tf_family+".txt")
            nr_pbm_file    = os.path.join(options.pbm_dir, "nr", "family_"+tf_family+".txt")
            folds_name     = "family_"+tf_family
            nr_triads_file = os.path.join(options.output_dir, "nr", "NR_family_"+tf_family+".txt")
            if not os.path.exists(list_file):
               if options.verbose: sys.stdout.write("\t\t\t--write %s...\n"%(list_file))
               sys.stdout.flush()
               for tf in list_of_tfs:
                   functions.write(list_file,tf)
            else:
               if options.verbose: sys.stdout.write("\t\t\t--use %s...\n"%(list_file))
               sys.stdout.flush()
            if not os.path.exists(folds_file):
               folds = set()
               if options.verbose: sys.stdout.write("\t\t\t--make list of folds %s\n"%(folds_name))
               for threading_file in sorted(os.listdir(os.path.join(options.pbm_dir, "threading"))):
                   if "index.txt" in threading_file: continue
                   m = re.search("(\S+).\d+.(\S{4}\_\S).txt", threading_file)
                   # Skip if PDB chain not in fold #
                   if not m:
                      sys.stderr.write("\t\tWrong threading file %s\n"%threading_file)
                      continue
                   #print("CHECK",m.group(1))
                   if m.group(1) in list_of_tfs:
                      folds.add(m.group(2))
               add_folds=set()
               for pdb_chain in folds:
                   pdb_nr_file = os.path.join(options.pdb_dir,"nr",pdb_chain+".txt")
                   for line in functions.parse_file(pdb_nr_file):
                       if line.startswith("#"): continue
                       add_folds.add(os.path.basename(line.rstrip(".txt")))
               folds.update(add_folds)
               #print("FOLDS",folds)
               functions.write(folds_file,"#pdb_chain;tm-score")
               for pdb_chain in folds:
                   functions.write(folds_file,pdb_chain+";1.0")
               for pdb_chain in folds:
                   functions.write(nr_pdb_file,os.path.join(options.pdb_dir,"triads",pdb_chain+".txt"))
            if os.path.exists(nr_triads_file):
               if options.verbose: sys.stdout.write("\t\t\t-- parse non-redundant triads of %s\n"%(nr_triads_file))
               for line in functions.parse_file(nr_triads_file):
                   if line.startswith("#"):continue
                   tf_triad = os.path.basename(line).rstrip(".txt")
                   if list_tf_rejected is not None:
                      skip = False
                      for tf_out in list_tf_rejected:
                          if tf_triad == tf_out: skip = True
                      if skip: continue
                   nr_triads.add(line)
            else:
               if options.parallel:
                   # Submit to queue #
                   if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                   else: cluster_queue=config.get("Cluster", "cluster_queue")
                   label=" "
                   if non_redundant: label=" --nr "
                   dummy_folder = os.path.join(options.dummy_dir,"dummy_nr")
                   if not os.path.exists(dummy_folder): os.makedirs(dummy_folder)
                   if options.verbose: sys.stdout.write("Submit %s %s %s --pdb=%s -r %s -i %s -o %s --pbm=%s -t %s -l %s --dummy %s -v\n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "nr.py"), label,options.pdb_dir, os.path.abspath(os.getcwd()), folds_name ,nr_triads_file , options.pbm_dir, config.get("Parameters", "max_redundancy_family"),list_file,dummy_folder))
                   functions.submit_command_to_queue("%s %s %s --pdb=%s -r %s -i %s -o %s --pbm=%s -t %s -l %s --dummy %s -v" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "nr.py"), label,options.pdb_dir, os.path.abspath(os.getcwd()), folds_name ,nr_triads_file , options.pbm_dir, config.get("Parameters", "max_redundancy_family"),list_file,dummy_folder), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                   sys.stdout.flush()
               else:
                   nr_triads_family = nr.get_nr_triads(options.pdb_dir, threshold=float(config.get("Parameters", "max_redundancy_family")), pbm_dir=options.pbm_dir, list_file=list_file, non_redundant=non_redundant)
                   nr_triads.update(nr_triads_family)
                   for triad in nr_triads_family:
                       functions.write(nr_triads_file,triad)
        ##############################
        for nr_triad_file in nr_triads:
            functions.write(nr_file, nr_triad_file)
      # Skip if nr file already exists #
      #elif not os.path.exists(nr_file):
        # Get redundant triads #
      #  nr_triads = nr.get_nr_triads(options.pdb_dir, threshold=float(config.get("Parameters", "max_redundancy_family")), pbm_dir=options.pbm_dir, non_redundant=non_redundant)
        # For each nr triads object... #
      #  for nr_triads_obj in nr_triads:
      #          functions.write(nr_file, nr_triads_obj._file)
    # Exit if stops here #
    if options.stop_step == 1:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 2. Split nr set into n folds #
    ################################
    if options.verbose: sys.stdout.write("Split non-redundant sets into folds...\n\n")
    # Skip if starts later #
    if options.start_step <= 2:

        ##############################
        # 2.1 Split general nr file  #
        ##############################
        if options.verbose: sys.stdout.write("\t-- split nr PBM general set into folds...\n")
        # Dont use split on PDB, it can not be compatred with PDB+PBM neither with e-scores
        folds_file = os.path.join(options.output_dir, "folds", "general.txt")
        # Skip if folds file already exists #
        if not os.path.exists(folds_file) :
            # Initialize #
            pbm = []
            pdb = []
            # For each line... #
            for line in functions.parse_file(os.path.join(options.pdb_dir,"nr","general.txt")):
                pdb.append(line)
            for line in functions.parse_file(os.path.join(options.pbm_dir, "nr", "general.txt")):
                pbm.append(line)
        if not os.path.exists(folds_file):
            # Get folds #
            folds_pbm = split_list_into_n_sublists(pbm, options.folds)
            # Write output #
            functions.write(folds_file, "#filename;fold")
            # For each fold... #
            for i in range(len(folds_pbm)):
                # Skip if nr file already exists #
                nr_file = os.path.join(options.output_dir, "nr", "general.%s.txt" % (i + 1))
                if not os.path.exists(nr_file):
                    for j in range(len(folds_pbm)):
                        if j != i:
                            for file_name in folds_pbm[j]:
                                functions.write(nr_file, file_name)
                    # For each file... #
                    for file_name in folds_pbm[i]:
                        functions.write(folds_file, "%s;%s" % (file_name, i + 1))
                else:
                    for file_name in folds_pbm[i]:
                        functions.write(folds_file, "%s;%s" % (file_name, i + 1))


        ##############################
        # 2.2 Split family nr file   #
        ##############################
        if options.verbose: sys.stdout.write("\t-- split nr PBM family set into folds...\n")
        # Skip if folds file already exists #
        tf_families=set()
        for textfile in os.listdir(os.path.join(options.output_dir,"nr")):
            if not textfile.startswith("NR_family"): continue
            tf_families.add(textfile.rstrip("txt").rstrip("\."))
        # get families of structures from PBM directory
        pdb_chains={}
        for line in functions.parse_file(os.path.join(options.pbm_dir, "families.txt")):
                if line.startswith("#"): continue
                line = line.split(";")
                pdb_chains.setdefault(line[0], set())
                pdb_chains[line[0]].add(line[-1].replace("/",":").replace(" ","_"))
        # Get NR folds  family-specific
        family_file = os.path.join(options.output_dir, "folds", "family.txt")
        functions.write(family_file, "#filename;fold")
        for tf_fam in tf_families:
          check_fam = "_".join(tf_fam.split("_")[2:])
          if options.verbose:sys.stdout.write("\t\t-- Family %s\n"%(tf_fam))
          sys.stdout.flush()
          folds_file = os.path.join(options.output_dir, "folds", tf_fam+".txt")
          if not os.path.exists(folds_file) :
            # Initialize #
            pbm = []
            pdb = []
            files = {}
            # For each line... #
            for line in functions.parse_file(os.path.join(options.output_dir, "nr",tf_fam+".txt")):
                # If PDB line... #
                if "/pdb/triads/" in line: pdb.append(line)
                pbm.append(line)
            # Get folds #
            folds_pbm = split_list_into_n_sublists(pbm, options.folds)
            # Write output #
            functions.write(folds_file, "#filename;fold")
            # For each fold... #
            for i in range(options.folds):
                # For each file... #
                for file_name in folds_pbm[i]:
                    functions.write(folds_file , "%s;%s" % (file_name, i + 1))
                    functions.write(family_file, "%s;%s" % (file_name, i + 1))
                    files.setdefault(os.path.basename(file_name), i)
          else:
            if options.verbose:sys.stdout.write("\t\tUse %s\n"%(folds_file))
            folds_f = open(folds_file,"r")
            # Initialize #
            pbm = []
            pdb = []
            files = {}
            folds_pbm = []
            for line in folds_f:
                if line.startswith("#"):continue
                file_name,ifold = line.split(";")
                files.setdefault(os.path.basename(file_name), int(ifold)-1)
            folds_f.close()
          # For each nr file... #
          for file_name in os.listdir(os.path.join(options.pbm_dir, "nr")):
                if options.verbose: sys.stdout.write("\t\t-- Check file %s\n"%(os.path.join(options.pbm_dir, "nr", file_name)))
                if file_name == "general.txt": continue
                if file_name.startswith("family"):continue
                try:
                 pdb_chain = file_name[:6]
                 # get family of pdb_chain
                 if pdb_chain not in pdb_chains:continue
                 family = pdb_chains[pdb_chain]
                 skip = True
                 for test_fam in check_fam.split(":"):
                   if skip:
                      if test_fam in pdb_chains[pdb_chain] or check_fam in pdb_chains[pdb_chain]: skip= False
                 if skip:
                    if options.verbose: sys.stdout.write("\t\t\t-- Skip, wrong family %s != %s \n"%(str(family),check_fam))
                    sys.stdout.flush()
                    continue
                 # For each fold... #
                 for i in range(options.folds):
                    # Skip if nr file already exists #
                    nr_file = os.path.join(options.output_dir, "nr", "%s.%s.txt" % (pdb_chain, i + 1))
                    if options.verbose: sys.stdout.write("\t\t\t-- prepare NR data %s.[%s].%s.txt for potentials\n"%(pdb_chain,tf_fam, i + 1))
                    sys.stdout.flush()
                    n=0
                    if not os.path.exists(nr_file):
                        # For each line... #
                        for line in functions.parse_file(os.path.join(options.pbm_dir, "nr", file_name)):
                            # If PDB line... #
                            tf_triad = os.path.basename(line)
                            if tf_triad in files:
                                if files[tf_triad] == i: continue
                                functions.write(nr_file, line)
                                n = n +1
                        if options.verbose: sys.stdout.write("\t\t\t-- Selected %d NR data %s.[%s].%s.txt for potentials using %s\n"%(n,pdb_chain,tf_fam, i + 1,file_name))
                    else:
                        if options.verbose: sys.stdout.write("\t\t\t-- Use existing %d NR data %s.[%s].%s.txt for potentials\n"%(n,pdb_chain,tf_fam, i + 1))
                except Exception as e:
                    if options.verbose: sys.stdout.write("\t\t\t-- Error in family %s for %s skip due to exception %s\n"%(tf_fam,pdb_chain,e))
          else:
             for line in functions.parse_file(folds_file):
                 if line.startswith("#"): continue
                 functions.write(family_file, line)
                 
             

    # Exit if stops here #
    if options.stop_step == 2:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 3. Derive stat. potentials   #
    ################################
    if options.verbose: sys.stdout.write("Derive statistical potentials...\n\n")
    # Skip if starts later #
    if options.start_step <= 3:
        # For each nr file... #
        for nr_file in os.listdir(os.path.join(options.output_dir, "nr")):
            # Skip if not a fold file... #
            if nr_file == "family.txt": continue
            if nr_file.startswith("TFs_family"): continue
            if nr_file.startswith("NR_family"): continue
            # Initialize #
            #m = re.search("(\S+).(\S+).(\d+).txt$", nr_file)
            m = re.search("(\S+).(\d+).txt$", nr_file)
            if not m: continue
            #if options.verbose: print("\tCheck potential ",m.group(1),  m.group(2),m.group(3))
            if options.verbose: print("\tCheck potential ",m.group(1),  m.group(2))
            # Skip if potentials file already exists #
            #potentials_file = os.path.join(options.output_dir, "potentials", m.group(1) + "." + m.group(2) + "." + m.group(3) + ".txt")
            potentials_file = os.path.join(options.output_dir, "potentials", m.group(1) + "." + m.group(2)  + ".txt")
            if not os.path.exists(potentials_file):
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    functions.submit_command_to_queue("%s %s -i %s -o %s -s -b -z" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                # Else... #
                else:
                    # Derive statistical potentials #
                    pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = spotentials.get_statistical_potentials(os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), approach=False, smooth=True, zscores=True,computation=True, dummy_dir=options.dummy_dir)
                    # Initialize #
                    potentials_obj = spotentials.Potentials()
                    # Add statistical potentials to object #
                    potentials_obj._pmf_3d = pmf_3d
                    potentials_obj._pmf_3dc = pmf_3dc
                    potentials_obj._pmf_s3dc = pmf_s3dc
                    potentials_obj._pmf_s3dc_dd = pmf_s3dc_dd
                    potentials_obj._pmf_s3dc_di = pmf_s3dc_di
                    potentials_obj._pmf_local = pmf_local
                    potentials_obj._pmf_pair = pmf_pair
                    potentials_obj._distances = distances
                    # Write statistical potentials #
                    potentials_obj.write(potentials_file)
            # Skip if potentials file already exists #
            #potentials_file = os.path.join(options.output_dir, "potentials", m.group(1) + "." + m.group(2) + "." + m.group(3) + ".taylor.txt")
            potentials_file = os.path.join(options.output_dir, "potentials", m.group(1) + "." + m.group(2)  + ".taylor.txt")
            if not os.path.exists(potentials_file):
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    functions.submit_command_to_queue("%s %s -i %s -o %s -a -s -b -z" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                # Else... #
                else:
                    # Derive statistical potentials #
                    pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = spotentials.get_statistical_potentials(os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), approach=True, smooth=True, zscores=True,computation=True, dummy_dir=options.dummy_dir)
                    # Initialize #
                    potentials_obj = spotentials.Potentials()
                    # Add statistical potentials to object #
                    potentials_obj._pmf_3d = pmf_3d
                    potentials_obj._pmf_3dc = pmf_3dc
                    potentials_obj._pmf_s3dc = pmf_s3dc
                    potentials_obj._pmf_s3dc_dd = pmf_s3dc_dd
                    potentials_obj._pmf_s3dc_di = pmf_s3dc_di
                    potentials_obj._pmf_local = pmf_local
                    potentials_obj._pmf_pair = pmf_pair
                    potentials_obj._distances = distances
                    # Write statistical potentials #
                    potentials_obj.write(potentials_file)
    # Exit if stops here #
    if options.stop_step == 3:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 4. Benchmark potentials      #
    ################################
    if options.verbose: sys.stdout.write("Benchmark potentials on non-redundant sets...\n\n")
    # Skip if starts later #
    if options.start_step <= 4:
        # For each subdir... #
        for subdir in ["pdb.general", "pdb.general.taylor", "pbm.general", "pbm.general.taylor", "pdb.family", "pdb.family.taylor", "pbm.family", "pbm.family.taylor"]:
            if not os.path.exists(os.path.join(options.output_dir, "results", subdir)):
                os.makedirs(os.path.join(options.output_dir, "results", subdir))
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "folds", "general.txt")):
            if line.startswith("#"): continue
            # Initialize #
            potentials = {}
            file_name, fold = line.split(";")
            m = re.search("/(T\S+.\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
            if not m: continue
            tf_id = m.group(1)
            pdb_chain = m.group(2)
            kmer = m.group(3)
            start = m.group(4)
            end = m.group(5)
            # Get triads file #
            triads_file = os.path.abspath(os.path.join(options.pbm_dir, "triads", "%s.%s.%s.%s-%s.txt" % (tf_id, pdb_chain, kmer, start, end)))
            # Get X3DNA file #
            x3dna_file = os.path.abspath(os.path.join(options.pdb_dir, "x3dna", pdb_chain[:4] + ".txt"))

            ##############################
            # 4.1 PDB general potentials #
            ##############################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pdb.general", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.pdb_dir, "potentials", "general.bins.txt"))
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s\n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))

            ###################################################
            # 4.2 PDB general potentials approached by Taylor #
            ###################################################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pdb.general.taylor", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.pdb_dir, "potentials", "general.taylor.bins.txt"))
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))

            ####################################
            # 4.3 PDB + PBM general potentials #
            ####################################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pbm.general", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.output_dir, "potentials", "general." + fold + ".txt"))
                #get the fold of the PDB file of the template
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))

            #########################################################
            # 4.4 PDB + PBM general potentials approached by Taylor #
            #########################################################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pbm.general.taylor", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.output_dir, "potentials", "general." + fold + ".taylor.txt"))
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))

        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "folds", "family.txt")):
            if line.startswith("#"): continue
            # Initialize #
            potentials = {}
            file_name, fold = line.split(";")
            m = re.search("/(T\S+.\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
            if not m: continue
            tf_id = m.group(1)
            pdb_chain = m.group(2)
            kmer = m.group(3)
            start = m.group(4)
            end = m.group(5)
            # Get triads file #
            triads_file = os.path.abspath(file_name)
            # Get X3DNA file #
            x3dna_file = os.path.abspath(os.path.join(options.pdb_dir, "x3dna", pdb_chain[:4] + ".txt"))

            #############################
            # 4.5 PDB family potentials #
            #############################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pdb.family", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.pdb_dir, "potentials", pdb_chain + ".bins.txt"))
                if not os.path.exists(potentials_file): 
                   sys.stdout.write("\t\t-- Not found potential %s\n"%(potentials_file))
                   continue
                sys.stdout.flush()
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))

            ##################################################
            # 4.6 PDB family potentials approached by Taylor #
            ##################################################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pdb.family.taylor", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.pdb_dir, "potentials", pdb_chain + ".taylor.bins.txt"))
                if not os.path.exists(potentials_file): 
                   sys.stdout.write("\t\t-- Not found potential %s\n"%(potentials_file))
                   continue
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))

            ###################################
            # 4.7 PDB + PBM family potentials #
            ###################################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pbm.family", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.output_dir, "potentials", pdb_chain  +"."+ fold + ".txt"))
                if not os.path.exists(potentials_file): 
                   sys.stdout.write("\t\t-- Not found potential %s\n"%(potentials_file))
                   continue
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))

            #########################################################
            # 4.8 PDB + PBM general potentials approached by Taylor #
            #########################################################
            results_file = os.path.abspath(os.path.join(options.output_dir, "results", "pbm.family.taylor", tf_id + "." + pdb_chain + "." + start + "-" + end + "." + fold + ".txt"))
            if options.verbose: sys.stdout.write("\t-- Check %s\n"%(results_file))
            sys.stdout.flush()
            # Skip if results file already exists #
            if not os.path.exists(results_file):
                # Statistical potentials file #
                potentials_file = os.path.abspath(os.path.join(options.output_dir, "potentials", pdb_chain + "."+ fold + ".taylor.txt"))
                if not os.path.exists(potentials_file): 
                   sys.stdout.write("\t\t-- Not found potential %s\n"%(potentials_file))
                   continue
                # If parallel... #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if options.verbose: sys.stdout.write("\t-- submit %s %s -c %s -o %s -p %s -t %s -x %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
                    sys.stdout.flush()
                    functions.submit_command_to_queue("%s %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("python %s -c %s -o %s -p %s -t %s -x %s" % (os.path.join(scripts_path, "benchmark.results.py"), pdb_chain[-1], results_file, potentials_file, triads_file, x3dna_file))
    # Exit if stops here #
    if options.stop_step == 4:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 5. Parse benchmarked files   #
    ################################
    if options.verbose: sys.stdout.write("Parse benchmark files...\n\n")
    label = " "
    if options.gc_content: label += " -g "
    label += " -r %d "%(int(options.randoms))
    # Skip if starts later #
    if options.start_step <= 5:
        # For each time... #
        for n in range(options.negatives):
            # If parallel... #
            if options.parallel:
                if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                else: cluster_queue=config.get("Cluster", "cluster_queue")
                # Submit to queue #
                if options.verbose: sys.stdout.write("\t-- submit %s %s -n %s -o %s --pbm=%s  %s \n" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.parsed.py"), str(n + 1), os.path.abspath(options.output_dir), os.path.abspath(options.pbm_dir),label))
                functions.submit_command_to_queue("%s %s -n %s -o %s --pbm=%s  %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "benchmark.parsed.py"), str(n + 1), os.path.abspath(options.output_dir), os.path.abspath(options.pbm_dir),label), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            else:
                os.system("python %s -n %s -o %s --pbm=%s %s " % (os.path.join(scripts_path, "benchmark.parsed.py"), str(n + 1), options.output_dir, options.pbm_dir,label))
    # Exit if stops here #
    if options.stop_step == 5:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 6. Create input files 4 figs #
    ################################
    if options.verbose: sys.stdout.write("Create input files for figures...\n\n")
    # Skip if starts later #
    if options.start_step <= 6:
        # Initialize #
        families = set()
        families_filter = set(config.get("Parameters", "families_filter").split(","))
        scores = numpy.arange(0, 1.01, step=0.01)
        # Read PDB family file #
        for line in functions.parse_file(os.path.join(options.output_dir, "parsed", "pdb.family.1.txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[7] not in families_filter: families.add(line[7])
        families = sorted(families, key=lambda x: x.upper())
        # For each time... #
        for n in range(options.negatives):
            # Skip if AUCPR/Precision/Coverage for PDB general file already exists #
            aucpr_pdb_general_file = os.path.join(options.output_dir, "parsed", "pdb.general.aucpr.%s.txt" % str(n + 1))
            auroc_pdb_general_file = os.path.join(options.output_dir, "parsed", "pdb.general.auroc.%s.txt" % str(n + 1))
            precision_pdb_general_file = os.path.join(options.output_dir, "parsed", "pdb.general.precision.%s.txt" % str(n + 1))
            recall_pdb_general_file = os.path.join(options.output_dir, "parsed", "pdb.general.recall.%s.txt" % str(n + 1))
            fpr_pdb_general_file = os.path.join(options.output_dir, "parsed", "pdb.general.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pdb_general_file) or not os.path.exists(auroc_pdb_general_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pdb.general.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pdb_general_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pdb_general_file, recall_pdb_general_file,fpr_pdb_general_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pdb_general_file)
            # Skip if AUCPR/Precision/Coverage for PDB general Taylor's file already exists #
            aucpr_pdb_general_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.general.taylor.aucpr.%s.txt" % str(n + 1))
            auroc_pdb_general_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.general.taylor.auroc.%s.txt" % str(n + 1))
            precision_pdb_general_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.general.taylor.precision.%s.txt" % str(n + 1))
            recall_pdb_general_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.general.taylor.recall.%s.txt" % str(n + 1))
            fpr_pdb_general_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.general.taylor.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pdb_general_taylor_file) or not os.path.exists(auroc_pdb_general_taylor_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pdb.general.taylor.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pdb_general_taylor_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pdb_general_taylor_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pdb_general_taylor_file, recall_pdb_general_taylor_file,fpr_pdb_general_taylor_file)
            # Skip if AUCPR/Precision/Coverage for PDB + PBM general file already exists #
            aucpr_pbm_general_file = os.path.join(options.output_dir, "parsed", "pbm.general.aucpr.%s.txt" % str(n + 1))
            auroc_pbm_general_file = os.path.join(options.output_dir, "parsed", "pbm.general.auroc.%s.txt" % str(n + 1))
            precision_pbm_general_file = os.path.join(options.output_dir, "parsed", "pbm.general.precision.%s.txt" % str(n + 1))
            recall_pbm_general_file = os.path.join(options.output_dir, "parsed", "pbm.general.recall.%s.txt" % str(n + 1))
            fpr_pbm_general_file = os.path.join(options.output_dir, "parsed", "pbm.general.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pbm_general_file) or not os.path.exists(auroc_pbm_general_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pbm.general.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pbm_general_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pbm_general_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pbm_general_file, recall_pbm_general_file,fpr_pbm_general_file)
            # Skip if AUCPR/Precision/Coverage for PDB + PBM general Taylor's file already exists #
            aucpr_pbm_general_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.general.taylor.aucpr.%s.txt" % str(n + 1))
            auroc_pbm_general_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.general.taylor.auroc.%s.txt" % str(n + 1))
            precision_pbm_general_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.general.taylor.precision.%s.txt" % str(n + 1))
            recall_pbm_general_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.general.taylor.recall.%s.txt" % str(n + 1))
            fpr_pbm_general_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.general.taylor.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pbm_general_taylor_file) or not os.path.exists(auroc_pbm_general_taylor_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pbm.general.taylor.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pbm_general_taylor_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pbm_general_taylor_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pbm_general_taylor_file, recall_pbm_general_taylor_file, fpr_pbm_general_taylor_file)
            # Skip if AUCPR/Precision/Coverage for PDB family file already exists #
            aucpr_pdb_family_file = os.path.join(options.output_dir, "parsed", "pdb.family.aucpr.%s.txt" % str(n + 1))
            auroc_pdb_family_file = os.path.join(options.output_dir, "parsed", "pdb.family.auroc.%s.txt" % str(n + 1))
            precision_pdb_family_file = os.path.join(options.output_dir, "parsed", "pdb.family.precision.%s.txt" % str(n + 1))
            recall_pdb_family_file = os.path.join(options.output_dir, "parsed", "pdb.family.recall.%s.txt" % str(n + 1))
            fpr_pdb_family_file = os.path.join(options.output_dir, "parsed", "pdb.family.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pdb_family_file) or not os.path.exists(auroc_pdb_family_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pdb.family.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pdb_family_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pdb_family_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pdb_family_file, recall_pdb_family_file,  fpr_pdb_family_file)
            # Skip if AUCPR/Precision/Coverage for PDB family Taylor's file already exists #
            aucpr_pdb_family_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.family.taylor.aucpr.%s.txt" % str(n + 1))
            auroc_pdb_family_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.family.taylor.auroc.%s.txt" % str(n + 1))
            precision_pdb_family_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.family.taylor.precision.%s.txt" % str(n + 1))
            recall_pdb_family_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.family.taylor.recall.%s.txt" % str(n + 1))
            fpr_pdb_family_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.family.taylor.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pdb_family_taylor_file) or not os.path.exists(auroc_pdb_family_taylor_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pdb.family.taylor.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pdb_family_taylor_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pdb_family_taylor_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pdb_family_taylor_file, recall_pdb_family_taylor_file, fpr_pdb_family_taylor_file)
            # Skip if AUCPR/Precision/Coverage for PDB + PBM family file already exists #
            aucpr_pbm_family_file = os.path.join(options.output_dir, "parsed", "pbm.family.aucpr.%s.txt" % str(n + 1))
            auroc_pbm_family_file = os.path.join(options.output_dir, "parsed", "pbm.family.auroc.%s.txt" % str(n + 1))
            precision_pbm_family_file = os.path.join(options.output_dir, "parsed", "pbm.family.precision.%s.txt" % str(n + 1))
            recall_pbm_family_file = os.path.join(options.output_dir, "parsed", "pbm.family.recall.%s.txt" % str(n + 1))
            fpr_pbm_family_file = os.path.join(options.output_dir, "parsed", "pbm.family.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pbm_family_file) or not os.path.exists(auroc_pbm_family_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pbm.family.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pbm_family_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pbm_family_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pbm_family_file, recall_pbm_family_file, fpr_pbm_family_file)
            # Skip if AUCPR/Precision/Coverage for PDB + PBM family Taylor's file already exists #
            aucpr_pbm_family_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.family.taylor.aucpr.%s.txt" % str(n + 1))
            auroc_pbm_family_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.family.taylor.auroc.%s.txt" % str(n + 1))
            precision_pbm_family_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.family.taylor.precision.%s.txt" % str(n + 1))
            recall_pbm_family_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.family.taylor.recall.%s.txt" % str(n + 1))
            fpr_pbm_family_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.family.taylor.fpr.%s.txt" % str(n + 1))
            if not os.path.exists(aucpr_pbm_family_taylor_file) or not os.path.exists(auroc_pbm_family_taylor_file):
                data_frame = get_data_frame(os.path.join(options.output_dir, "parsed", "pbm.family.taylor.%s.txt" % str(n + 1)))
                get_aucpr_data(data_frame, families, options.folds, aucpr_pbm_family_taylor_file)
                get_auroc_data(data_frame, families, options.folds, auroc_pbm_family_taylor_file)
                get_precision_recall_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pbm_family_taylor_file, recall_pbm_family_taylor_file)
                get_precision_recall_fpr_data(data_frame, families, scores, options.folds, float(1) / options.randoms, precision_pbm_family_taylor_file, recall_pbm_family_taylor_file, fpr_pbm_family_taylor_file)
        # Skip if Wilcoxon file already exists #
        wilcoxon_file = os.path.join(options.output_dir, "parsed", "wilcoxon.txt")
        if not os.path.exists(wilcoxon_file):
            # Initialize #
            wilcoxon_tests = []
            # PDB general potentials #
            auprc_values = []
            auprc_taylor_values = []
            # For each time... #
            for n in range(options.negatives):
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pdb.general.aucpr.%s.txt" % str(n + 1))):
                    auprc_values.append(value)
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pdb.general.taylor.aucpr.%s.txt" % str(n + 1))):
                    auprc_taylor_values.append(value)
            mean, mean_taylor, taylor, p_value = one_sided_wilcoxon_test(auprc_values, auprc_taylor_values)
            wilcoxon_tests.append(["pdb.general", mean, mean_taylor, taylor, p_value])
            # PDB + PBM general potentials #
            auprc_values = []
            auprc_taylor_values = []
            # For each time... #
            for n in range(options.negatives):
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pbm.general.aucpr.%s.txt" % str(n + 1))):
                    auprc_values.append(value)
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pbm.general.taylor.aucpr.%s.txt" % str(n + 1))):
                    auprc_taylor_values.append(value)
            mean, mean_taylor, taylor, p_value = one_sided_wilcoxon_test(auprc_values, auprc_taylor_values)
            wilcoxon_tests.append(["pbm.general", mean, mean_taylor, taylor, p_value])
            # PDB family potentials #
            auprc_values = []
            auprc_taylor_values = []
            # For each time... #           
            for n in range(options.negatives):
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pdb.family.aucpr.%s.txt" % str(n + 1))):
                    auprc_values.append(value)
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pdb.family.taylor.aucpr.%s.txt" % str(n + 1))):
                    auprc_taylor_values.append(value)
            mean, mean_taylor, taylor, p_value = one_sided_wilcoxon_test(auprc_values, auprc_taylor_values)
            wilcoxon_tests.append(["pdb.family", mean, mean_taylor, taylor, p_value])
            # PDB + PBM family potentials #
            auprc_values = []
            auprc_taylor_values = []
            # For each time... #
            for n in range(options.negatives):
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pbm.family.aucpr.%s.txt" % str(n + 1))):
                    auprc_values.append(value)
                for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pbm.family.taylor.aucpr.%s.txt" % str(n + 1))):
                    auprc_taylor_values.append(value)
            mean, mean_taylor, taylor, p_value = one_sided_wilcoxon_test(auprc_values, auprc_taylor_values)
            wilcoxon_tests.append(["pbm.family", mean, mean_taylor, taylor, p_value])
            # For each family... #
            for family in families:
                family
                # Initialize #
                auprc_values = []
                auprc_taylor_values = []
                for n in range(options.negatives):
                    for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pbm.family.aucpr.%s.txt" % str(n + 1)), family):
                        auprc_values.append(value)
                    for value in get_auprc_values(os.path.join(options.output_dir, "parsed", "pbm.family.taylor.aucpr.%s.txt" % str(n + 1)), family):
                        auprc_taylor_values.append(value)
                mean, mean_taylor, taylor, p_value = one_sided_wilcoxon_test(auprc_values, auprc_taylor_values)
                wilcoxon_tests.append([family, mean, mean_taylor, taylor, p_value])
            # Write output #
            functions.write(wilcoxon_file, "#family;mean;mean_taylor;taylor;p-value")
            # For each wilcoxon test... #
            for family, mean, mean_taylor, taylor, p_value in wilcoxon_tests:
                functions.write(wilcoxon_file, "%s;%s;%s;%s;%s" % (family, mean, mean_taylor, str(taylor), str(p_value)))
    # Exit if stops here #
    if options.stop_step == 6:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 7. Create input files 4 figs #
    ################################
    if options.verbose: sys.stdout.write("Create figures...\n\n")
    # Skip if starts later #
    if options.start_step <= 7:

        ####################
        # 7.0. Tables      #
        ####################
        file_table_aucpr =  os.path.join(options.output_dir, "figures", "aucpr.csv")
        file_table_auroc =  os.path.join(options.output_dir, "figures", "auroc.csv")
        if not os.path.exists(file_table_aucpr) or not os.path.exists(file_table_auroc):
            data_prc =[] 
            data_roc =[] 
            aucpr={}
            auroc={}
            families = set()
            # For each time... #
            for n in range(options.negatives):
                # Files #
                files = [os.path.join(options.output_dir, "parsed", "pdb.general.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pdb.general.taylor.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.general.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.general.taylor.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pdb.family.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pdb.family.taylor.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.family.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.family.taylor.aucpr.%s.txt" % str(n + 1))]
                # For each file... #
                for file_name in files:
                    word = os.path.basename(file_name).split(".")
                    feature = ".".join(word[:word.index("aucpr")])
                    if os.path.exists(file_name):
                     faucpr =  open(file_name,"r")
                     for line in faucpr:
                        if line.startswith("#"):continue
                        scores = line.split(";")
                        family = scores[0]
                        families.add(family)
                        for x in scores[1:]:
                            try:
                               if x is not None: aucpr.setdefault((family,feature),[]).append(float(x))
                            except:
                               continue
                     faucpr.close()
                    file_auroc = os.path.join(options.output_dir, "parsed",feature+".auroc.%s.txt" % str(n + 1))
                    if os.path.exists(file_auroc):
                     fauroc =  open(file_auroc,"r")
                     for line in fauroc:
                        if line.startswith("#"):continue
                        scores = line.split(";")
                        family = scores[0]
                        families.add(family)
                        for x in scores[1:]:
                            try:
                               if x is not None: auroc.setdefault((family,feature),[]).append(float(x))
                            except:
                               continue
                     fauroc.close()
            features = ["pdb.general","pdb.general.taylor","pbm.general","pbm.general.taylor","pdb.family","pdb.family.taylor","pbm.family","pbm.family.taylor"]
            
            for fam in sorted(families):
                data=[]
                data.append(fam)
                for feature in features:
                   try:
                    if aucpr[(fam,feature)]:
                       x    = numpy.array(aucpr[(fam,feature)])
                       AUC  = x.mean()
                       RMSD = x.std()
                       data.append(AUC)
                       data.append(RMSD)
                    else:
                       data.append(None)
                       data.append(None)
                   except:
                       data.append(None)
                       data.append(None)
                data_prc.append(data)
                data=[]
                data.append(fam)
                for feature in features:
                   try:
                    if auroc[(fam,feature)]:
                       x    = numpy.array(auroc[(fam,feature)])
                       AUC  = x.mean()
                       RMSD = x.std()
                       data.append(AUC)
                       data.append(RMSD)
                    else:
                       data.append(None)
                       data.append(None)
                   except:
                       data.append(None)
                       data.append(None)
                data_roc.append(data)
            columns = ["Family","pdb.general","STD","pdb.general.taylor","STD","pbm.general","STD","pbm.general.taylor","STD","pdb.family","STD","pdb.family.taylor","STD","pbm.family","STD","pbm.family.taylor","STD"]
            table_aucpr = pandas.DataFrame(data_prc)
            table_aucpr.columns=columns
            table_auroc = pandas.DataFrame(data_roc)
            table_auroc.columns=columns
            table_aucpr.to_csv(file_table_aucpr)
            table_auroc.to_csv(file_table_auroc)


        ####################
        # 7.1. Boxplots #
        ####################
        boxplot_file = os.path.join(options.output_dir, "figures", "aucpr.boxplot.svg")
        # Skip if tsplot file already exists #
        if not os.path.exists(boxplot_file):
            data = []
            # For each time... #
            for n in range(options.negatives):
                # Files #
                files = [os.path.join(options.output_dir, "parsed", "pdb.general.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pdb.general.taylor.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.general.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.general.taylor.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pdb.family.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pdb.family.taylor.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.family.aucpr.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.family.taylor.aucpr.%s.txt" % str(n + 1))]
                # For each file... #
                for file_name in files:
                    # Initialize #
                    file_name = file_name
                    data_type = "pdb"
                    if "pbm" in file_name: data_type = "pdb+pbm"
                    potential_type = "general"
                    if "family" in file_name: potential_type = "family"
                    taylors_approach = "No"
                    if "taylor" in file_name: taylors_approach = "Yes"
                    # For each line... #
                    for line in functions.parse_file(file_name):
                        if line.startswith("#"): continue
                        if line.startswith("All"): continue
                        line = line.split(";")
                        family = line.pop(0)
                        for i in range(len(line)):
                            data.append(numpy.array([file_name, data_type, potential_type, taylors_approach, family, (i + 1), None]))
                            if line[i] != "None": data[-1][-1] = float(line[i])
            # Get data frame #
            data_frame = pandas.DataFrame(data)
            data_frame.columns = ["file", "data", "potential", "taylor", "family", "fold", "auprc"]
            # Seaborn's context #
            seaborn.set(context="paper", style="ticks")
            # Initialize plot #
            fig, axes = plot.subplots(nrows=1, ncols=2, sharex=False, sharey=False, figsize=(2.4803133333333336, 2.4803133333333336))
            # For each ax... #
            for i, ax in enumerate(axes):
                # Initialize #
                if i == 0:
                    df = data_frame[(data_frame["potential"] == "general")]
                else:
                    df = data_frame[(data_frame["potential"] == "family")]
                # Make boxplots #
                seaborn.boxplot(x="data", y="auprc", data=df, hue="taylor", palette={"No": "white", "Yes": "#4477AA"}, saturation=1, ax=ax, flierprops=dict(marker="o", markersize=1.5))
                # Tweek figure #
                ax.set(ylim=(-0.05, 0.95))
                # Set axes #
                if i == 0:
                    ax.set_xlabel("general")
                    ax.set_ylabel("AUCPR")
                    ax.set_yticks(numpy.arange(0, 0.9 + 0.1, 0.1))
                    ax.legend(handles=[mpatches.Patch(facecolor="#4477AA", edgecolor="black", linewidth=0.5, label="Taylor's approach")], loc="upper left")
                else:
                    ax.set_xlabel("family")
                    ax.set_ylabel("")
                    ax.set_yticks([])
                    ax.spines["left"].set_visible(False)
                    ax.tick_params(left="off")
                    ax.legend().set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.tick_params(right="off")
                ax.spines["top"].set_visible(False)
                ax.tick_params(top="off")
                ax.set_xticklabels(["PDB", "PDB+PBM"])
                # Fixing titles #
                ax.set_title("")
                # # Make square plot #
                # ax.set(aspect=(2. / ax.get_data_ratio()))
            # Save figure #
            fig.savefig(boxplot_file, bbox_inches="tight", pad_inches=0.1, transparent=True)

      
        #########################
        # 7.2. Time-series plot #
        #########################
        data_frame_file = os.path.join(options.output_dir, "figures", "precision_recall.tsplot.csv")
        family_data_frame_file = os.path.join(options.output_dir, "figures", "family.precision_recall.tsplot.csv")
        max_precision_file = os.path.join(options.output_dir, "figures", "max_precision.txt")
        # Skip if data frame file already exists #
        if not os.path.exists(data_frame_file):
            # Initialize #
            data = []
            ignore = set()
            ir_pbm_family_file = os.path.join(options.output_dir, "figures", "pbm.family.ir.txt")
            cov_pbm_family_file = os.path.join(options.output_dir, "figures", "pbm.family.cov.txt")
            ir_pbm_family_taylor_file = os.path.join(options.output_dir, "figures", "pbm.family.taylor.ir.txt")
            cov_pbm_family_taylor_file = os.path.join(options.output_dir, "figures", "pbm.family.taylor.cov.txt")
            files = [ir_pbm_family_file, cov_pbm_family_file, ir_pbm_family_taylor_file, cov_pbm_family_taylor_file]        
            # For each time... #
            for n in range(options.negatives):
                # Files #
                files = [os.path.join(options.output_dir, "parsed", "pbm.family.precision.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.family.taylor.precision.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.family.recall.%s.txt" % str(n + 1)),
                         os.path.join(options.output_dir, "parsed", "pbm.family.taylor.recall.%s.txt" % str(n + 1))]
                # For each file... #
                for file_name in files:
                    statistic = "Precision"
                    if "recall" in file_name: statistic = "Recall"
                    taylor = False
                    if "taylor" in file_name: taylor = True
                    # For each line... #
                    for line in functions.parse_file(file_name):
                        if line.startswith("#"): continue
                        if line.startswith("All"): continue
                        line = line.split(";")
                        family = line.pop(0)
                        score = line.pop(0)
                        for i in range(len(line)):
                            if float(score) == float(0) and line[i] == "None": ignore.add((family, str(i + 1)))
                            if (family, str(i + 1)) in ignore: continue
                            try:
                                data.append([float(score), statistic, "%s;%s;%s" % (family, str(i + 1), n), float(line[i])])
                            except:
                                data.append([float(score), statistic, "%s;%s;%s" % (family, str(i + 1), n), float(0)])
            # Sort data #
            data.sort(key=lambda x: x[2])
            data.sort(key=lambda x: x[0])
            # Get data frame #
            data_frame = pandas.DataFrame(data)
            data_frame.columns = ["score", "statistic", "family+fold+n", "value"]
            data_frame.to_csv(data_frame_file)
        # Skip if family data frame file already exists #
        if not os.path.exists(family_data_frame_file):
            # Initialize #
            data = []
            # For each line... #
            for line in functions.parse_file(data_frame_file):
                line = line.split(",")
                if line[1] == "score": continue
                family, fold, n = line[3].split(";")
                data.append([float(line[1]), line[2], family, "%s;%s" % (fold, n), float(line[4])])
            # Get data frame #
            data_frame = pandas.DataFrame(data)
            data_frame.columns = ["score", "statistic", "family", "fold+n", "value"]
            data_frame.to_csv(family_data_frame_file)
        # Skip deprecated TSplots
        print("Done")
        exit()

"""
        # Skip if tsplot file already exists #
        tsplot_file = os.path.join(options.output_dir, "figures", "precision_recall.tsplot.svg")
        if not os.path.exists(tsplot_file):
            # Initialize #
            data_frame = pandas.read_csv(data_frame_file, index_col=0, parse_dates=True)
            # Seaborn's context #
            seaborn.set(context="paper", style="ticks")
            # Initialize plot #
            fig, ax = plot.subplots(figsize=(2.4803133333333336, 2.4803133333333336))
            # Twin the x-axis to make 2 independent y-axes #
            axes = [ax, ax.twinx()]
            # Make plots #
            seaborn.lineplot(data_frame[(data_frame["statistic"] == "Precision")].reset_index(), time="score", unit="family+fold+n", value="value", ci=100, color="#4477AA", legend=False, ax=axes[0])
            seaborn.lineplot(data_frame[(data_frame["statistic"] == "Recall")].reset_index(), time="score", unit="family+fold+n", value="value", ci=100, color="#CC6677", legend=False, ax=axes[1])
            # Make labels/spines/etc. #
            axes[0].set_xlabel("Score")
            axes[0].set_ylabel("Precision")
            axes[0].set_ylim([0, 75])
            axes[0].spines["top"].set_visible(False)
            axes[0].tick_params(top="off")
            axes[1].set_ylabel("Recall")
            axes[1].set_ylim([0, 106.6666666666666666667])
            axes[1].spines["top"].set_visible(False)
            axes[1].tick_params(top="off")
            axes[1].spines["left"].set_visible(False)
            axes[1].tick_params(left="off")
            axes[1].spines["bottom"].set_visible(False)
            axes[1].tick_params(bottom="off")
            # # Get precision, recall and score #
            # precision, recall, score = get_precision_recall_score(data_frame)
            # functions.write(max_precision_file, "%s;%s;%s" % (precision, recall, score))
            # Set score line #
            ax.axvline(0.95, color="black", linestyle=":", linewidth=2)
            # Set legend #
            axes[1].legend([axes[0].get_lines()[0], axes[1].get_lines()[0], axes[0].get_lines()[-1]], ["Precision", "Recall", "Score = 0.95"], loc="center left")
            # # Make square plot #
            # axes[0].set(adjustable="box-forced", aspect=(1. / axes[0].get_data_ratio()))
            # axes[1].set(adjustable="box-forced", aspect=(1. / axes[1].get_data_ratio()))
            # Save figure #
            fig.savefig(tsplot_file, bbox_inches="tight", pad_inches=0.1, transparent=True)
        # Initialize #
        families = [["AP2", "APSES", "ARID/BRIGHT", "bHLH", "bZIP", "C2H2 ZF", "CxxC", "DM", "E2F", "Ets", "Forkhead", "GATA"],
                    ["GCM", "Homeodomain", "IRF", "Myb/SANT", "NAC/NAM", "Nuclear receptor", "POU", "Paired box", "Rel", "Sox", "T-box", "Zinc cluster"]]
        complexes = [[75, 16, 24, 906, 350, 1010, 75, 19, 45, 1324, 594, 225],
                     [15, 7254, 283, 400, 36, 2263, 779, 59, 21, 159, 194, 40]]
        # # For each line... #
        # for line in functions.parse_file(max_precision_file):
        #      # Get threshold precision, recall and score #
        #     thd_precision, thd_recall, thd_score = line.split(";")
        # For each family... #
        for n in range(len(families)):
            # Skip if tsplot file already exists #
            tsplot_file = os.path.join(options.output_dir, "figures", "family.precision_recall.tsplot.%s.svg" % (n + 1))
            if not os.path.exists(tsplot_file):
                # Initialize #
                data_frame = pandas.read_csv(family_data_frame_file, index_col=0, parse_dates=True)
                # Seaborn's context #
                seaborn.set(context="paper", palette="colorblind", style="ticks")
                # Initialize plot #
                fig, axes = plot.subplots(nrows=4, ncols=3, sharex=False, sharey=False, figsize=(7.44094, 9.921253333333334))
                # For each row... #
                for i, row in enumerate(axes):
                    # For each ax... #
                    for j, ax in enumerate(row):
                        # Initialize #
                        family = i * 3 + j
                        print(families[n][family])
                        # If family plot... #
                        if family < 12:
                            df = data_frame[(data_frame["family"] == families[n][family])]
                            # Make plots #
                            seaborn.lineplot(df[(df["statistic"] == "Precision")].reset_index(), time="score", unit="fold+n", value="value", ci=95, color="#4477AA", legend=False, ax=ax)
                            seaborn.lineplot(df[(df["statistic"] == "Recall")].reset_index(), time="score", unit="fold+n", value="value", ci=95, color="#CC6677", legend=False, ax=ax)
                            # Tweek figure #
                            ax.set(ylim=(0, 105))
                            # Don't show anything... #
                            ax.spines["top"].set_visible(False)
                            ax.tick_params(top="off")
                            ax.tick_params(bottom="off")
                            ax.tick_params(left="off")
                            ax.tick_params(right="off")
                            ax.set_xlabel("")
                            ax.set_xticks([])
                            ax.set_ylabel("")
                            ax.set_yticks([])
                            # Set axes #
                            if j < 1:
                                ax.set_ylabel("Precision")
                                ax.set_yticks([0, 25, 50, 75, 100])
                                ax.tick_params(left="on")
                            if i == 3:
                                ax.set_xlabel("Score")
                                ax.set_xticks(numpy.arange(0, 1.5, 0.5))
                                ax.tick_params(bottom="on")
                            if j == 2:
                                ax.yaxis.tick_right()
                                ax.yaxis.set_label_position("right")
                                ax.set_ylabel("Recall")
                                ax.set_yticks([0, 25, 50, 75, 100])
                                ax.tick_params(right="on")
                            # # Get precision, recall and score at threshold #
                            # precision, recall, score = get_precision_recall_score(df, thd_score=float(thd_score))
                            # # If precision < 75%, get an appropriate threshold #
                            # if precision < 75:
                            #     precision, recall, score = get_precision_recall_score(df, thd_score=float(thd_score), thd_precision=75)
                            # Set score line #
                            ax.axvline(0.95, color="black", linestyle=":", linewidth=2)
                            # ax.plot([score, score], [0, max_precision], label="Score", color="gray", linestyle=":", linewidth=1)
                            # ax.plot([score], [precision], "o", markersize=5, label="Max. precision", color="black")
                            # Set title #
                            ax.set_title(families[n][family])
                            # Set # of complexes #
                            if families[n][family] == "Ndt80/PhoG":
                                ax.text(0.075, 3, complexes[n][family], fontsize=8)
                            else:
                                ax.text(0.015, 3, complexes[n][family], fontsize=8)
                            # Set legend #
                            if family == 14:
                                # legend = [ax.get_lines()[0], ax.get_lines()[1],  ax.get_lines()[-2], ax.get_lines()[-1]]
                                legend = [ax.get_lines()[0], ax.get_lines()[1],  ax.get_lines()[-1]]
                        # Else... #
                        else:
                            ax.axis("off")
                            # ax.legend(legend, ["Precision", "Recall", "75% Precision", "Precision @ 0.92"], bbox_to_anchor=(1.3, 0.75))
                            ax.legend(legend, ["Precision", "Recall", "Score = 0.95"], bbox_to_anchor=(1.3, 0.75))
                        # # Make square plot #
                        # ax.set(adjustable="box-forced", aspect=(1./ax.get_data_ratio()))
                        
                # Save figure #
                fig.savefig(tsplot_file, bbox_inches="tight", pad_inches=0.1, transparent=True)
"""
