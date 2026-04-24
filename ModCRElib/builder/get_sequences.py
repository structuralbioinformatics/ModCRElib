"""
Extract protein sequences from SQL INSERT dumps and write FASTA output.
"""

import re
import sys,os

def sql_to_fasta(sql_file, fasta_out,choice=1):
    """
    Parse SQL value tuples and export selected identifier column as FASTA headers.

    Args:
        sql_file (str): SQL INSERT dump file.
        fasta_out (str): Output FASTA path.
        choice (int): Field index used as FASTA header.

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

            if len(fields) < 5:
                continue  # skip malformed rows

            header = fields[choice]  # second element
            seq = fields[4].replace("\n", "").replace("\r", "").replace("*","").strip()

            fasta.write(f">{header}\n")

            # wrap lines to 60 characters
            for i in range(0, len(seq), 60):
                fasta.write(seq[i:i+60] + "\n")

    print(f"✅ Extracted {len(entries)} sequences into {fasta_out}")


if __name__ == "__main__":
    # Step 1) Read CLI positional arguments.
    infile =sys.argv[1]
    outfile=sys.argv[2]
    choice=int(sys.argv[3])
    # Step 2) Convert SQL tuples into FASTA.
    sql_to_fasta(infile, outfile,choice)
