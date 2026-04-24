#!/usr/bin/env python3
"""
Split a SQL dump into one file per table, including DDL and INSERTs.

Writes files to output_dir/<table_name>.sql
"""

import sys,os
import re
from collections import defaultdict

# Regex to capture table names in CREATE and INSERT statements (handles backticks and simple schema.table)
RE_CREATE = re.compile(r'^\s*CREATE\s+TABLE\s+(?:IF\s+NOT\s+EXISTS\s+)?(?:`(?P<bt>[^`]+)`|(?P<plain>[A-Za-z0-9_\.]+))', re.IGNORECASE)
RE_INSERT = re.compile(r'^\s*INSERT\s+INTO\s+(?:`(?P<ibt>[^`]+)`|(?P<iplain>[A-Za-z0-9_\.]+))', re.IGNORECASE)
RE_DROP = re.compile(r'^\s*DROP\s+TABLE\s+(?:IF\s+EXISTS\s+)?(?:`(?P<dbt>[^`]+)`|(?P<dplain>[A-Za-z0-9_\.]+))', re.IGNORECASE)
RE_ALTER = re.compile(r'^\s*ALTER\s+TABLE\s+(?:`(?P<abt>[^`]+)`|(?P<aplain>[A-Za-z0-9_\.]+))', re.IGNORECASE)

def normalize_table_name(raw):
    """
    Normalize SQL table identifiers for filesystem-friendly naming.

    Args:
        raw (str | None): Raw table identifier, possibly with schema prefix.

    Returns:
        str | None: Table name without schema/backticks, or ``None`` for null input.
    """
    if raw is None:
        return None
    # if schema.table format, keep only the last part
    if '.' in raw:
        raw = raw.split('.')[-1]
    return raw.replace('`', '').strip()

def open_for_table(output_dir, table,rootname):
    """
    Open (append mode) the per-table SQL output file.

    Args:
        output_dir (str): Destination directory.
        table (str): Normalized table name.
        rootname (str): Prefix included in output filename.

    Returns:
        tuple: ``(file_handle, path)`` for the table SQL file.
    """
    safe = re.sub(r'[^\w\-\.]', '_', table)
    path = os.path.join(output_dir, f"{rootname}.{safe}.sql")
    # open in append mode so INSERTs can be added as we go
    return open(path, "a", encoding="utf-8"), path

def split_dump(input_file, output_dir="tables_sql",rootname="CisBP"):
    """
    Split a monolithic SQL dump into per-table SQL files.

    Captures ``CREATE TABLE``, ``INSERT INTO`` (including multiline inserts),
    and optional ``DROP/ALTER TABLE`` statements; ignores global directives.

    Args:
        input_file (str): Source SQL dump.
        output_dir (str): Destination folder for split table files.
        rootname (str): Prefix used in generated filenames.

    Returns:
        None.
    """
    os.makedirs(output_dir, exist_ok=True)

    current_table = None          # name while we're inside a CREATE TABLE statement
    capture_create = False        # True when we're collecting a CREATE TABLE statement
    create_buffer = []            # lines for DDL until its terminating semicolon
    open_files = {}               # map table -> open file handle
    file_paths = {}               # map table -> path
    counts = defaultdict(int)     # counts for DDL/INSERT lines written

    # helper to ensure file handle for a table exists (open in append)
    def ensure_file(table):
        if table not in open_files:
            fh, p = open_for_table(output_dir, table,rootname)
            open_files[table] = fh
            file_paths[table] = p
        return open_files[table]

    # stream the input
    with open(input_file, "r", encoding="utf-8", errors="ignore") as infile:
        for raw_line in infile:
            line = raw_line.rstrip("\n")

            # If currently capturing CREATE TABLE block:
            if capture_create:
                create_buffer.append(line + "\n")
                # Detect end of CREATE (a semicolon that terminates the statement).
                # We assume a semicolon at the end of a line ends the statement.
                if line.strip().endswith(";"):
                    # flush create_buffer to file for current_table
                    if current_table:
                        fh = ensure_file(current_table)
                        fh.writelines(create_buffer)
                        counts[(current_table, "ddl")] += len(create_buffer)
                        print(f"Saved DDL for table {current_table} -> {file_paths[current_table]}")
                    # reset
                    create_buffer = []
                    capture_create = False
                    current_table = None
                continue  # go to next line

            # Not currently inside a CREATE: check for CREATE start
            m_create = RE_CREATE.match(line)
            if m_create:
                raw_table = m_create.group("bt") or m_create.group("plain")
                tbl = normalize_table_name(raw_table)
                current_table = tbl
                capture_create = True
                create_buffer = [line + "\n"]
                # If the CREATE DDL is one-line and ends with ;, handle immediately
                if line.strip().endswith(";"):
                    fh = ensure_file(current_table)
                    fh.writelines(create_buffer)
                    counts[(current_table, "ddl")] += len(create_buffer)
                    print(f"Saved DDL for table {current_table} -> {file_paths[current_table]}")
                    create_buffer = []
                    capture_create = False
                    current_table = None
                continue

            # Check for INSERT INTO statements (single- or multi-line)
            m_ins = RE_INSERT.match(line)
            if m_ins:
                raw_table = m_ins.group("ibt") or m_ins.group("iplain")
                tbl = normalize_table_name(raw_table)
                fh = ensure_file(tbl)
                fh.write(line + "\n")
                counts[(tbl, "insert")] += 1
                # If the INSERT is multi-line, keep appending following lines until semicolon
                if not line.strip().endswith(";"):
                    # read and append subsequent lines until semicolon is reached
                    for follow in infile:
                        fh.write(follow)
                        counts[(tbl, "insert")] += 1
                        if follow.rstrip("\n").strip().endswith(";"):
                            break
                continue

            # Optionally capture DROP TABLE / ALTER TABLE lines into corresponding table files
            m_drop = RE_DROP.match(line)
            if m_drop:
                raw_table = m_drop.group("dbt") or m_drop.group("dplain")
                tbl = normalize_table_name(raw_table)
                fh = ensure_file(tbl)
                fh.write(line + "\n")
                counts[(tbl, "drop")] += 1
                continue

            m_alt = RE_ALTER.match(line)
            if m_alt:
                raw_table = m_alt.group("abt") or m_alt.group("aplain")
                tbl = normalize_table_name(raw_table)
                fh = ensure_file(tbl)
                fh.write(line + "\n")
                counts[(tbl, "alter")] += 1
                continue

            # Otherwise ignore (global statements like SET, USE, LOCK TABLES, etc.)
            # If you'd like to prepend a header with these global statements to every file,
            # you could capture them here and write to each open file.

    # Close files
    for fh in open_files.values():
        fh.close()

    # Print summary
    tables = sorted(file_paths.keys())
    print("\nDone. Summary:")
    print(f"Total tables written: {len(tables)}")
    for t in tables:
        ddl = counts.get((t, "ddl"), 0)
        ins = counts.get((t, "insert"), 0)
        dr = counts.get((t, "drop"), 0)
        alt = counts.get((t, "alter"), 0)
        print(f" - {t}: DDL-lines={ddl}, INSERT-lines={ins}, DROP={dr}, ALTER={alt}, file={file_paths[t]}")

if __name__ == "__main__":
    infile=sys.argv[1]
    outdir=sys.argv[2]
    rootname=sys.argv[3]
    split_dump(infile, outdir,rootname)

