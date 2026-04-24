"""
Export feature-like FASTA entries from SQL INSERT tuples.
"""

import re
import sys,os

def sql_to_fasta(sql_file, fasta_out):
    """
    Convert SQL tuple entries into FASTA records.

    Args:
        sql_file (str): SQL INSERT dump file.
        fasta_out (str): Output FASTA path.

    Returns:
        None.
    """

    # Regex to capture tuples inside INSERT INTO ... VALUES (...),(...)
    tuple_pattern = re.compile(r"\((.*?)\)", re.DOTALL)

    with open(sql_file, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    # Remove the beginning "INSERT INTO ..." and trailing semicolon
    content = re.sub(r"^INSERT INTO.*?VALUES", "", content, flags=re.IGNORECASE | re.DOTALL)
    content = content.strip().rstrip(";")

    # Find all tuples of values
    entries = tuple_pattern.findall(content)

    with open(fasta_out, "w", encoding="utf-8") as fasta:
        for entry in entries:
            # Split by commas not inside quotes
            fields = re.findall(r"'(.*?)'", entry, re.DOTALL)
            if len(fields) < 4:
                continue  # skip malformed rows

            header = fields[0]+"/"+fields[1]+"/"+fields[2]  # second element
            seq = fields[3].replace("\n", "").replace("\r", "").replace("*","").strip()

            fasta.write(f">{header}\n")

            # wrap lines to 60 characters
            for i in range(0, len(seq), 60):
                fasta.write(seq[i:i+60] + "\n")

    print(f"✅ Extracted {len(entries)} sequences into {fasta_out}")


if __name__ == "__main__":
    # Step 1) Read positional input/output arguments.
    infile =sys.argv[1]
    outfile=sys.argv[2]
    # Step 2) Convert SQL dump records into FASTA.
    sql_to_fasta(infile, outfile)
