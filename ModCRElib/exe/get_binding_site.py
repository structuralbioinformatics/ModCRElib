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

# Import jbonet's module #
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.contacts import contacts
from ModCRElib.structure.contacts import interface

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    Parse command-line options for binding-site sequence extraction.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for PDB input, gap size, and output.
    """

    parser = optparse.OptionParser("python get_binding_site.py -i INPUT_FILE [-d DISTANCE_TYPE --dummy DUMMY_DIR -g GAP -o OUTPUT_FILE]")

    parser.add_option("-d", default="basepairs", action="store", type="string", dest="distance_type", help="Distance type (i.e. \"basepairs\", \"dinucleotides\" or \"mindist\"; default = dinucleotides)", metavar="{string}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-g", default=1,action="store", type="int", dest="gap", help="Gap of nucleotides between binding-sites (default=1)", metavar="{integer}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if options.distance_type != "basepairs" and options.distance_type != "dinucleotides" and options.distance_type != "mindist":
        parser.error("incorrect distance type: accepted values are \"basepairs\" (for contact distance between amino acid cb and geometric center of a basepair), \"dinucleotides\" (for contact distance between amino acid cb and geometric center of a dinucleotide) and \"mindist\" (for minimum distance)")

    return options

def main():
    """
    Run the binding-site sequence extraction workflow.

    Workflow:
        1. Parse command-line arguments.
        2. Load structure, DNA annotations, contacts, and interface tables.
        3. Merge nearby interface fragments and extract DNA sequences.
        4. Write auxiliary outputs and extracted sequences.

    Args:
        None.

    Returns:
        None. Binding-site sequences and optional helper outputs are written.
    """

    # Arguments & Options #
    options = parse_options()

    # Get PDB object #
    pdb_obj = PDB(os.path.abspath(options.input_file))

    # Get X3DNA object #
    x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Get contacts object #
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj, "pdi", options.distance_type, os.path.abspath(options.dummy_dir))

    # Get interface object #
    interface_obj = interface.get_interface_obj(pdb_obj, x3dna_obj, contacts_obj, os.path.abspath(options.dummy_dir))

    # Get interface sequence
    gap=options.gap
    min_motif_size=int(config.get("Parameters", "min_motif_size"))
    start=[]
    end=[]
    sequences=[]
    interface=interface_obj.get_interface()
    check=False
    bp_list = [basepair for basepair in sorted(interface)]
    for i in range(len(bp_list)):
        basepair = bp_list[i]
        if interface[basepair]['contacts'] > 0:
            if not check:
                check=True
                start.append(i)
        elif check:
            check=False
            for j in range(i,i+gap):
                if check: continue
                gap_basepair = bp_list[j]
                if interface[gap_basepair]['contacts'] > 0:
                    check=True
            if not check:
                end.append(i)
    if check: end.append(i)
    if options.output_file is not None: print("Starts ",start,"\nEnds ",end)
    fragments=zip(start,end)
    for start,end in fragments:
      if (end-start) > int(config.get("Parameters", "min_motif_size")):
        interface_start = bp_list[start]
        interface_end   = bp_list[end]
        sequence        = x3dna_obj.get_nucleotide_sequence(interface_start,interface_end)
        if options.output_file is not None:
           print("\tFragment sequence: ",sequence," located in ", interface_start,interface_end)
        sequences.append(sequence)

    # Output #
    if options.output_file is not None:
        outinterface =options.output_file+".interface.out"
        outx3dna     =options.output_file+".x3dna.out"
        outcontact   =options.output_file+".contact.out"
        interface_obj.write(os.path.abspath(outinterface))
        x3dna_obj.write(os.path.abspath(outx3dna))
        contacts_obj.write(os.path.abspath(outcontact))
        for sequence in sequences:
            functions.write(options.output_file,sequence)
    else:
        for sequence in sequences:
            sys.stdout.write("%s\n"%sequence)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()

