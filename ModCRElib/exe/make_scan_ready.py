import re
import sys

sourceDatabaseTXT = sys.argv[1]
scanDatabaseTXT = sys.argv[2]

with open(sourceDatabaseTXT) as f_in, open(scanDatabaseTXT, "w") as f_out:
    for line in f_in:
        if line.startswith("MOTIF"):
            # Find both model strings
            parts = re.findall(r'MODEL_[^ ]+', line)
            new_parts = []
            for p in parts:
                # Example: MODEL_sp_P35869_AHR_HUMAN:34:272_5nj8_A_1
                m = re.search(r'_(\w{4})_([A-Z])_(\d+)', p)
                if m:
                    pdb, chain, num = m.groups()
                    new_parts.append(f"{pdb}_{num}{chain}")
            if len(new_parts) == 2:
                line = f"MOTIF {new_parts[0]} {new_parts[1]}\n"
        f_out.write(line)
