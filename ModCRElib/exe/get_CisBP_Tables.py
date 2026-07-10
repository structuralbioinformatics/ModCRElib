"""
Generate normalized CIS-BP SQL tables from mixed SQL/TSV sources.
"""

import os, sys, re
import optparse
from Bio import motifs as mm
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import configparser
import pandas as pd
import numpy as np
from functools import reduce

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
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBILib.structure import PDB


def parse_options():
    """
    Parse CLI options for CIS-BP table reconstruction.

    How to run:
        ``python get_CisBP_Tables.py --sql_motif motifs.sql --tf tf_info.tsv --ps proteins.tsv``

    Returns:
        optparse.Values: Namespace with SQL/TSV input files and output root.
    """

    parser = optparse.OptionParser("python --sql_motif sql_motif_file  --tf tf_info_file --ps proteins_file [ -o rootname -v --sql_tf sql_tf_file ]")

    parser.add_option("--sql_motif", action="store", default=None, type="string", dest="sql_motif_file", help="SQL motifs file (from CIS-BP)", metavar="{file}")
    parser.add_option("--sql_tf", action="store", default=None, type="string", dest="sql_tf_file", help="SQL TF file (from CIS-BP)", metavar="{file}")
    parser.add_option("--tf", action="store", default=None, type="string", dest="tfs_file", help="TF info of all (from CIS-BP)", metavar="{file}")
    parser.add_option("--ps", action="store", default=None, type="string", dest="prot_file", help="Protein sequences (from CIS-BP)", metavar="{file}")
    parser.add_option("-o", action="store", type="string",default="CisBP",dest="root", help="CisBP file in SQL format with table of motifs (default CisBP)", metavar="{rootname}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.prot_file is None or options.tfs_file is None  or options.sql_motif_file is None:
         parser.error("missing arguments: type option \"-h\" for help")

    return options

def parse_sql_inserts(filename, table_name, columns):
    """
    Parse SQL INSERT statements for one table into a DataFrame.

    Args:
        filename (str): SQL file path.
        table_name (str): SQL table name to extract.
        columns (list[str]): Output DataFrame columns.

    Returns:
        pandas.DataFrame: Parsed table rows.
    """
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    # Find all INSERT statements for this table
    pattern = rf"INSERT INTO\s+`{table_name}`.*?VALUES\s*(.*?);"
    matches = re.findall(pattern, content, flags=re.S)
    all_rows = []
    for match in matches:
        # Split rows by '),(' and clean parentheses
        rows = re.findall(r"\((.*?)\)", match, flags=re.S)
        for row in rows:
            # Split values by commas not inside quotes
            parts = re.findall(r"(?:'[^']*'|[^,]+)", row)
            cleaned = [x.strip().strip("'") if x.strip() != "NULL" else None for x in parts]
            all_rows.append(cleaned)
    # Create DataFrame
    return pd.DataFrame(all_rows, columns=columns)


def pandas_to_sql(df, table_name, output_file, varchar_len=5000):
    """
    Write a SQL file to create and populate a table from a DataFrame.

    Args:
        df (pandas.DataFrame): Data to export.
        table_name (str): SQL table name to create.
        output_file (str): Output SQL file path.
        varchar_len (int): Default varchar length for text-like columns.

    Returns:
        None.
    """
    with open(output_file, "w", encoding="utf-8") as f:
        # --- CREATE TABLE statement
        f.write(f"CREATE TABLE `{table_name}` (\n")
        for col, dtype in zip(df.columns, df.dtypes):
            # Map pandas dtypes to SQL types
            if pd.api.types.is_integer_dtype(dtype):
                sql_type = "INT"
            elif pd.api.types.is_float_dtype(dtype):
                sql_type = "FLOAT"
            elif pd.api.types.is_bool_dtype(dtype):
                sql_type = "BOOLEAN"
            else:
                sql_type = f"VARCHAR({varchar_len})"
            f.write(f"  `{col}` {sql_type},\n")
        f.seek(f.tell() - 2)  # remove last comma
        f.write("\n);\n\n")
        # --- INSERT statements
        for _, row in df.iterrows():
            values = []
            for val in row:
                if pd.isna(val):
                    values.append("NULL")
                else:
                    # escape single quotes
                    sval = str(val).replace("'", "''")
                    values.append(f"'{sval}'")
            values_str = ", ".join(values)
            f.write(f"INSERT INTO `{table_name}` VALUES ({values_str});\n")


def merge_keep_same_join_diff(dfs, on, sep=", "):
    """
    Merge multiple DataFrames while preserving agreements and joining conflicts.

    Args:
        dfs (list[pandas.DataFrame]): Input DataFrames sharing the merge keys.
        on (list[str] | str): Merge key column(s).
        sep (str): Separator used when joining distinct values.

    Returns:
        pandas.DataFrame: Merged table with harmonized non-key columns.
    """
    if isinstance(on, str):
        on = [on]
    # 1) Rename non-key columns in each df to make them unique per source:
    renamed_dfs = []
    all_base_cols = set()
    for i, df in enumerate(dfs, start=1):
        df = df.copy()
        src_tag = f"__src{i}"
        rename_map = {}
        for col in df.columns:
            if col not in on:
                rename_map[col] = f"{col}{src_tag}"
                all_base_cols.add(col)
        df = df.rename(columns=rename_map)
        renamed_dfs.append(df)
    # 2) Outer-merge all renamed dfs on the key columns
    merged = reduce(lambda a, b: pd.merge(a, b, on=on, how="outer"), renamed_dfs)
    # 3) For each base column, find all source-specific columns and collapse them
    result = merged[on].copy()  # start with keys
    for base_col in sorted(all_base_cols):
        # collect the source-specific column names that correspond to this base_col
        src_cols = [c for c in merged.columns if c.startswith(base_col + "__src")]
        if not src_cols:
            # no occurrences of this column in any source -> create as NaN
            result[base_col] = np.nan
            continue
        def collapse_row(values):
            # values is a list-like of values from the source-specific columns for a single row
            cleaned = []
            for v in values:
                if pd.isna(v):
                    continue
                # normalize to str and strip whitespace
                s = str(v).strip()
                if s == "":
                    continue
                cleaned.append(s)
            # preserve insertion order but remove exact duplicates
            uniq = []
            for x in cleaned:
                if x not in uniq:
                    uniq.append(x)
            if not uniq:
                return np.nan
            if len(uniq) == 1:
                return uniq[0]
            return sep.join(uniq)
        # apply per-row collapse across these src_cols
        result[base_col] = merged[src_cols].apply(lambda row: collapse_row(row.tolist()), axis=1)
    # Reorder columns: keys first (in same order), then base_cols
    final_cols = list(on) + [c for c in result.columns if c not in on]
    return result[final_cols]

#-------------#
# Main        #
#-------------#

def main():
    """
    Build normalized CIS-BP SQL tables from motif/TF/protein sources.

    Workflow:
        1. Define required ModCRE schemas.
        2. Parse CLI options and load source tables.
        3. Merge/aggregate records by TF identifier.
        4. Export one SQL file per target table.

    Returns:
        None.
    """
    # Step 1) Define target table schemas required by ModCRE.
    cisbp_tables={}
    cisbp_3_tables={}
    cisbp_tables.setdefault("tf_families",["Family_ID", "Family_Name", "DBDs", "DBD_Count", "Cutoff"])
    cisbp_tables.setdefault("motifs",["Motif_ID", "TF_ID", "MSource_ID", "DBID", "Motif_Type", "Motif_Sequence", "IUPAC", "IUPAC_REV"])
    cisbp_tables.setdefault("motif_sources",["MSource_ID", "MSource_Identifier", "MSource_Type", "MSource_Author", "MSource_Year", "PMID", "MSource_Version"])
    cisbp_tables.setdefault("tfs",["TF_ID", "Family_ID", "TSource_ID", "DBID", "TF_Name", "TF_Species", "TF_Status"])
    cisbp_tables.setdefault("proteins",["Protein_ID", "TF_ID", "DBID", "TF_Species", "Protein_Sequence"])
    cisbp_3_tables.setdefault("proteins",["Protein_ID", "TF_ID", "DBID", "Protein_Type", "Protein_Sequence"])

    # Step 2) Parse options and read input SQL/TSV sources.
    options = parse_options()
    #read tables
    motif_table= parse_sql_inserts(os.path.abspath(options.sql_motif_file),"motifs",cisbp_tables["motifs"])
    tf_info    = pd.read_csv(options.tfs_file,sep="\t",dtype="str")
    prot_info  = pd.read_csv(options.prot_file,sep="\t",dtype="str")
    if options.sql_tf_file is not None: 
       tf_table   = parse_sql_inserts(os.path.abspath(options.sql_tf_file),"tfs",cisbp_tables["tfs"])
    # Step 3) Merge all sources on TF_ID and restrict to selected TFs.
    if options.sql_tf_file is not None: 
       merged_all=merge_keep_same_join_diff([prot_info,tf_info,motif_table,tf_table],on="TF_ID")
    else:
       merged_all=merge_keep_same_join_diff([prot_info,tf_info,motif_table],on="TF_ID")
    #select the TF_IDs
    selected_ids = prot_info["TF_ID"].values
    filtered = merged_all[merged_all["TF_ID"].isin(selected_ids)]
    # Step 4) Export one SQL file per normalized target table.
    for table_name,table_keys in cisbp_tables.items():
        result   = (filtered.groupby(table_keys[0])[table_keys[1:]].agg(lambda x: ", ".join(sorted(set(x.dropna().astype(str))))).reset_index())
        output_file=options.root+"."+table_name+".sql"
        pandas_to_sql(result, table_name, output_file)


if __name__ == "__main__":
   main()


        
              
