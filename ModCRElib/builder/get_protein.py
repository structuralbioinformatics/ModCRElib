"""
Merge CIS-BP SQL tables into a compact TF/protein feature table.
"""

import re,sys,os
import pandas as pd

def parse_sql_inserts(filename, table_name, columns):
    """
    Parse SQL INSERT statements for one table into a DataFrame.

    Args:
        filename (str): SQL file path.
        table_name (str): SQL table name to extract.
        columns (list[str]): Column labels to assign.

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

# --- Define column names for each table (as in your SQL schemas)
domains_cols = ["Domain_ID", "Domain_Name", "Pfam_Name", "Pfam_DBID", "Inter_DBID", "Domain_Type"]
prot_features_cols = ["PFeature_ID", "Protein_ID", "Domain_ID", "ProtFeature_FromPos", "ProtFeature_ToPos", "ProtFeature_Sequence"]
proteins_cols = ["Protein_ID", "TF_ID", "DBID", "Protein_Type", "Protein_Sequence"]
tf_families_cols = ["Family_ID", "Family_Name", "DBDs", "DBD_Count", "Cutoff", "SR_Model", "SR_NoThreshold"]
tfs_cols = ["TF_ID", "Family_ID", "TSource_ID", "DBID", "TF_Name", "TF_Species", "TF_Status"]
motifs_cols = ["Motif_ID","TF_ID","MSource_ID","DBID","Motif_Type","Motif_Sequence","IUPAC","IUPAC_REV"]

if __name__ == "__main__":

   # Step 1) Read SQL file paths and output TSV path.
   SQLdomains         = sys.argv[1]
   SQLprot_features   = sys.argv[2]
   SQLproteins        = sys.argv[3]
   SQLtf_families     = sys.argv[4]
   SQLtfs             = sys.argv[5]
   SQLmotifs          = sys.argv[6]
   output             = sys.argv[7]
   if len(sys.argv)>7:
      experiment      = sys.argv[8]
   else:
      experiment=None
   
   # Step 2) Parse all required SQL tables.
   domains = parse_sql_inserts(SQLdomains, "domains", domains_cols)
   prot_features = parse_sql_inserts(SQLprot_features, "prot_features", prot_features_cols)
   proteins = parse_sql_inserts(SQLproteins, "proteins", proteins_cols)
   tf_families = parse_sql_inserts(SQLtf_families, "tf_families", tf_families_cols)
   tfs = parse_sql_inserts(SQLtfs, "tfs", tfs_cols)
   motifs = parse_sql_inserts(SQLmotifs,"motifs",motifs_cols)
   # Step 3) Merge relational tables on their key fields.
   merged = (
    tfs
    .merge(tf_families, on="Family_ID", how="left")
    .merge(proteins, on="TF_ID", how="left")
    .merge(motifs, on="TF_ID", how="left")
    .merge(prot_features, on="Protein_ID", how="left")
    .merge(domains, on="Domain_ID", how="left")
   )

   # Step 4) Optionally filter by motif experiment label.
   if experiment is not None:
      merged = merged[merged["Motif_Type"]==experiment]
   # Step 5) Keep/rename the core feature columns.
   df = merged[
    ["TF_ID", "TF_Species", "TF_Name", "Protein_ID", "Pfam_DBID","Family_ID", "Family_Name", "DBID_x", "DBDs", "Motif_Type", "Protein_Sequence"]
    ].rename(columns={"DBID_x": "DBID"})

   # Step 6) Aggregate repeated entries by TF_ID and export TSV.
   result = (
    df.groupby("TF_ID")[["TF_Species", "TF_Name", "Protein_ID", "Pfam_DBID","Family_ID", "Family_Name", "DBID", "DBDs", "Motif_Type", "Protein_Sequence"]]
      .agg(lambda x: ", ".join(sorted(set(x.dropna().astype(str)))))
      .reset_index()
   )
   # --- Export to TSV
   result.to_csv(output, sep="\t", index=False)

   print("✅ Saved table as %s"%output)
   print(result.head())


