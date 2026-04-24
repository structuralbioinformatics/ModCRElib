import os, sys, re
import optparse
import configparser
import shutil
import itertools as it

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
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases, dna_complementary
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.structure.contacts import contacts,interface,triads
from ModCRElib.structure.dna import x3dna, model_dna
from ModCRElib.structure.protein import dssp, model_protein
from ModCRElib.structure.threading import threader
from ModCRElib.sequence import blast, hmmer, nr, homologs



dna_complementary["X"]="J"
dna_complementary["J"]="X"
dna_complementary["O"]="Q"
dna_complementary["Q"]="O"
dna_complementary["N"]="N"
nitrogenous_bases["X"]=nitrogenous_bases["C"]
nitrogenous_bases["O"]=nitrogenous_bases["C"]
nitrogenous_bases["J"]=nitrogenous_bases["G"]
nitrogenous_bases["Q"]=nitrogenous_bases["G"]


def  get_triads(thr_obj,thr_dir,pdb_dir,verbose=False,dummy_dir="/tmp",clean=False):
    """
        Resolve and load triads/x3dna/pdb objects for a threading target.
    
        Args:
            thr_obj (Threaded): Parsed threading object.
            thr_dir (str): Directory containing threading/template files.
            pdb_dir (str): PDB resources directory.
            verbose (bool, optional): Enable verbose logs.
            dummy_dir (str, optional): Temporary directory root.
            clean (bool, optional): Clean PDB before deriving triads.
    
        Returns:
            tuple: `(triads_obj, pdb_obj, x3dna_obj)`.
    
    Args:
        thr_obj (Any): Value used by this routine.
        thr_dir (Any): Directory path used by this operation.
        pdb_dir (Any): Directory path used by this operation.
        verbose (Any): Boolean flag controlling routine behavior.
        dummy_dir (Any): Directory path used by this operation.
        clean (Any): Boolean flag controlling routine behavior."""
    # get pdb_name and chain of the template#
    dna_sequence={}
    pdb_name, pdb_chain = thr_obj._pdb_name.lower(), thr_obj._pdb_chain
    pdb_structure = os.path.join(pdb_name+"_"+pdb_chain)
    if not os.path.exists(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")):
       template_chain = os.path.abspath(os.path.join(thr_dir,thr_obj._pdb_name+"_"+pdb_chain+".pdb"))
       if not os.path.exists(template_chain):
          template = os.path.abspath(os.path.join(thr_dir,thr_obj._pdb_name+".pdb"))
       else:
          template = template_chain
       if not os.path.exists(template):
          sys.stdout.write("ERROR: file %s is not found. The PDB folder MUST contain the PDB of the thread\n"%(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")))
          sys.stdout.write("ERROR: file %s is not found. The template is not in the folder of the thread\n"%(template_chain))
          sys.stdout.write("ERROR: file %s is not found. The template is not in the folder of the thread\n"%(template))
          raise ValueError("No structure reference\n")
       else:
          try:
             if verbose: print("\t-- Get TRIADS of template",template)
             triads_obj, pdb_obj, x3dna_obj = get_triads_from_template(template,verbose,dummy_dir,clean)
          except Exception as e:
             raise ValueError(e)
    else:
       try:
             if verbose: print("\t-- Get TRIADS from folder",pdb_dir)
             triads_obj, pdb_obj, x3dna_obj =  get_triads_from_folder(thr_obj,pdb_dir)
       except Exception as e:
             raise ValueError(e)
    return triads_obj, pdb_obj, x3dna_obj

def  get_triads_from_folder(thr_obj,pdb_dir):
    """
        Load template triads and structural context from prepared PDB folders.
    
        Args:
            thr_obj (Threaded): Parsed threading object.
            pdb_dir (str): PDB resources directory.
    
        Returns:
            tuple: `(triads_obj, pdb_obj, x3dna_obj)`.
    
    Args:
        thr_obj (Any): Value used by this routine.
        pdb_dir (Any): Directory path used by this operation."""
    # get pdb_name and chain of the template#
    dna_sequence={}
    pdb_name, pdb_chain = thr_obj._pdb_name.lower(), thr_obj._pdb_chain
    pdb_structure = os.path.join(pdb_name+"_"+pdb_chain)
    if not os.path.exists(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")):
       sys.stdout.write("ERROR: file %s is not found. The PDB folder MUST contain the PDB of the thread"%(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")))
    # get helix structure #
    for helix in functions.parse_file(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")):
        # Get PDB object #
        pdb_obj = PDB(os.path.join(pdb_dir, "split", pdb_name + ".dna." + helix + ".pdb"))
        break
    pdb_obj_prot = PDB(os.path.join(pdb_dir, "split", pdb_structure + ".pdb"))
    pdb_obj.add_chains(pdb_obj_prot.chains)
    for chain_id in pdb_obj.chain_identifiers:
        chain=pdb_obj.get_chain_by_id(chain_id)
        if chain.chaintype=="N":
            dna_sequence.setdefault(chain_id,chain.nucleotide_sequence())
    # Get X3DNA object #
    x3dna_obj = x3dna.X3DNA(os.path.join(pdb_dir, "x3dna", pdb_name + ".txt"))

    # Get contacts object #
    contacts_obj = contacts.Contacts(os.path.join(pdb_dir, "contacts",pdb_name + ".txt"))

    # Skip if no contacts #
    if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")


    # Get helix #
    dna_helices={}
    for helix in x3dna_obj.get_dna_helices():
        # Get helix dinucleotides #
        dinucleotides = x3dna_obj.get_helix_dinucleotides(helix)
        # For each contact... #
        for contact_obj in contacts_obj.get_contacts():
            if contact_obj._A_chain == pdb_chain:
                dna_helices.setdefault(helix,set()).add(tuple(contact_obj._B_chain))

    # Get interface object #
    #interface_obj = interface.Interface(os.path.join(pdb_dir, "interfaces", pdb_structure + ".txt"))
    #interface_basepairs = interface_obj.get_interface_basepairs()
    #for dna_helix in dna_helices.iterkeys():
    #    dna_sequence_interface = dna_sequence[list(dna_helices[dna_helix])[0][0]][interface_obj.get_start()-1:interface_obj.get_interface_length()+interface_obj.get_start()+1]
    # Get nucleotide positions #

    # Get triads object #
    triads_obj = triads.Triads(os.path.join(pdb_dir, "triads", pdb_name + "_" + pdb_chain + ".txt"))

    return triads_obj, pdb_obj, x3dna_obj

def get_triads_from_template(template,verbose,dummy_dir,clean=False):
    """
        Compute triads directly from a template PDB file.
    
        Args:
            template (str): Template PDB path.
            verbose (bool): Enable verbose logs.
            dummy_dir (str): Temporary directory root.
            clean (bool, optional): Clean PDB before processing.
    
        Returns:
            tuple: `(triads_obj, pdb_obj, x3dna_obj)`.
    
    Args:
        template (Any): Value used by this routine.
        verbose (Any): Boolean flag controlling routine behavior.
        dummy_dir (Any): Directory path used by this operation.
        clean (Any): Boolean flag controlling routine behavior."""

    # Get PDB object #
    if verbose: sys.stdout.write("\t-- Reading PDB file %s ...\n"%(os.path.basename(template)))
    pdb_obj = PDB(template)
    if clean: pdb_obj.clean()
    dummy_pdb=os.path.join(dummy_dir,os.path.basename(template))
    if not os.path.exists(dummy_pdb):
       pdb_obj.write(dummy_pdb)

    # Get DSSP object #
    if verbose: sys.stdout.write("\t\t-- Get DSSP ...\n")
    dssp_obj = dssp.get_dssp_obj(dummy_pdb)

    # Get X3DNA object #
    if verbose: sys.stdout.write("\t\t-- Get DNA object ...\n")
    x3dna_obj = x3dna.get_x3dna_obj(dummy_pdb, dummy_dir)

    # Get contacts object #
    if verbose: sys.stdout.write("\t\t-- Get contacts ...\n")
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)

    # Skip if no contacts #
    if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")

    # Get triads object #
    triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)
    

    return triads_obj, pdb_obj, x3dna_obj


def thread_kmer_triads(kmer,triads_obj):
    """
        Thread one k-mer sequence onto an input triads object.
    
        Args:
            kmer (str): DNA k-mer sequence.
            triads_obj (Triads): Input triads object.
    
        Returns:
            Triads: Triads object updated for the threaded k-mer.
    
    Args:
        kmer (Any): Value used by this routine.
        triads_obj (Any): Triad object/data used for interaction modeling."""
    # Initialize #
    full_dna_seq = list(it.islice(kmer,0,None))
    full_dna_pos = list(range(1,len(full_dna_seq)+1))
    basepairs    = dict(list(zip(full_dna_pos,full_dna_seq)))
    thread_triads_obj = triads.Triads()
    kmer = list(kmer)
    # For each triad object... #
    for triad_obj in triads_obj.get_triads():
   	    # Initialize #
        aminoacid_environment = []
        dinucleotide_environment = []
        a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(';')
        # Get chain, number #
        chain, number = residue_A.split("-")
        nucleotide_environment = residue_B.split(",")
        dinucleotide_environment = b_ob.split("-")
        # Information stored in b_ob  
        # b, twobases, dna_strand, dna_groove, dna_chemical_group = b_ob.split('-')
        # Get aa_environment and change its residues according to sequence on strand in F2_aa #
        # Forward strand #
        if dinucleotide_environment[2] == 'F':
           if int(nucleotide_environment[0].split("-")[1]) in list(basepairs.keys()) and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()):
              aminoacid_environment = a_oa.split("-")
              # Information stored in a_oa  
              # a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
              # is preserved
              dinucleotide_environment[0] = basepairs[int(nucleotide_environment[0].split("-")[1])] + basepairs[int(nucleotide_environment[2].split("-")[1])]
              dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
              triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
              thread_triads_obj.add_triad(triad_obj)
        # Reverse strand #
        if dinucleotide_environment[2] == 'R':
           if int(nucleotide_environment[0].split("-")[1]) in list(basepairs.keys()) and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()):
              aminoacid_environment = a_oa.split("-")
              # Information stored in a_oa  
              # a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
              # is preserved
              dinucleotide_environment[0] = basepairs[int(nucleotide_environment[0].split("-")[1])] + basepairs[int(nucleotide_environment[2].split("-")[1])]
              dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
              triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
              thread_triads_obj.add_triad(triad_obj)
    # Return  new threaded triads
    return thread_triads_obj


def threading_triads(threading_file, pdb_dir=None,  template=None, verbose=True, dummy_dir="/tmp",clean=False):
    """
        Build threaded triads from one threading file and template context.
    
        Args:
            threading_file (str): Threading file path.
            pdb_dir (str, optional): PDB resources directory.
            template (str, optional): Explicit template PDB path override.
            verbose (bool, optional): Enable verbose logs.
            dummy_dir (str, optional): Temporary directory root.
            clean (bool, optional): Clean PDB before triad derivation.
    
        Returns:
            tuple: `(thread_triads_obj, pdb_obj, x3dna_obj)`.
    
    Args:
        threading_file (Any): Path to the input/output file.
        pdb_dir (Any): Directory path used by this operation.
        template (Any): Value used by this routine.
        verbose (Any): Boolean flag controlling routine behavior.
        dummy_dir (Any): Directory path used by this operation.
        clean (Any): Boolean flag controlling routine behavior."""
    thr_dir = os.path.dirname(threading_file)
    thr_obj = threader.Threaded(threading_file=threading_file)
    # add threading data to object #
    thr_obj._check_parsing()
    pdb_name, pdb_chain = thr_obj._pdb_name.lower(), thr_obj._pdb_chain
    pdb_structure = os.path.join(pdb_name+"_"+pdb_chain)
    # Get nucleotide positions #
    full_dna_seq = list(it.islice(list(thr_obj._dna.keys())[0],0,None))
    position=1 + int(thr_obj.get_kmer(list(thr_obj._dna.keys())[0]))
    full_dna_pos = list(range(position,len(full_dna_seq)+position))
    dna_int_obj = dict(list(zip(full_dna_pos,full_dna_seq)))
    #if template is None, use pdb_dir
    if template is None:
       try:
         triads_obj, pdb_obj, x3dna_obj  = get_triads(thr_obj,thr_dir,pdb_dir,verbose,dummy_dir,clean)
       except ValueError as e:
         if verbose: sys.stderr.write("Failed get triads from folder %s\n"%e)
         exit(0)
    else:
       try:
         triads_obj, pdb_obj, x3dna_obj  = get_triads_from_template(os.path.abspath(template),verbose,dummy_dir,clean)
       except ValueError as e:
         if verbose: sys.stderr.write("Failed get triads from template %s\n"%e)
         exit(0)


    # Initialize #
    protein = {}
    kmers = []
    for line in functions.parse_file(threading_file):
        if line == "": continue
        if line.startswith("#") or line.startswith("//"): continue
        elif line.startswith(">"):
            threading = line[1:]
        elif threading == "protein":
            line = line.split(";")
            protein.setdefault((line[0], int(line[1])), line[2])
        elif threading == "dna":
            kmers.append(line.split(";"))
    # For each k-mer... #
    for kmer, position in kmers:
        # Initialize #
        thread_triads_obj = triads.Triads()
        basepairs = {}
        kmer = list(kmer)
        # For each contact... #
        for pos in list(dna_int_obj.keys()):
            basepairs.setdefault(pos, dna_int_obj[pos])
        # Skip if triads file already exists #
        # Initialize #
        done = set()
        # For each triad object... #
        for triad_obj in triads_obj.get_triads():
        	# Initialize #
            aminoacid_environment = []
            dinucleotide_environment = []
            a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(';')
            if residue_A.split("-")[0] == pdb_chain:
                # Get chain, number #
                chain, number = residue_A.split("-")
                nucleotide_environment = residue_B.split(",")
                dinucleotide_environment = b_ob.split("-")
                # Information stored in b_ob  
                # b, twobases, dna_strand, dna_groove, dna_chemical_group = b_ob.split('-')
                # Get aa_environment and change its residues according to sequence on strand in F2_aa #
                # Forward strand #
                if (pdb_chain,int(number)) not in list(protein.keys()): continue
                if protein[(pdb_chain,int(number))] not in list(aminoacids1to3.keys()): continue
                if int(nucleotide_environment[0].split("-")[1]) not in basepairs: continue
                if int(nucleotide_environment[2].split("-")[1]) not in basepairs: continue
                if dinucleotide_environment[2] == 'F':
                    if int(nucleotide_environment[0].split("-")[1]) in list(basepairs.keys()) and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()) and int(nucleotide_environment[2].split("-")[1]) in list(basepairs.keys()):
                            aminoacid_environment = a_oa.split("-")
                            # Information stored in a_oa  
                            #a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                            try:
                              aminoacid_environment[0] = aminoacids1to3[protein[(pdb_chain,int(number))]]
                            except:
                              print(("PDB_CHAIN %s N %d PROTEIN %s"%(pdb_chain,int(number),str(protein))))
                              exit()
                            if aminoacids_polarity_boolean[protein[(pdb_chain,int(number))]]:
                                    hydrophobicity=aminoacid_environment[1] = "P"
                            else:
                                    hydrophobicity=aminoacid_environment[1] = "N"
                            if basepairs[int(nucleotide_environment[0].split("-")[1])] not in nitrogenous_bases: continue
                            if basepairs[int(nucleotide_environment[2].split("-")[1])] not in nitrogenous_bases: continue
                            try:
                              dinucleotide_environment[0] = basepairs[int(nucleotide_environment[0].split("-")[1])] + basepairs[int(nucleotide_environment[2].split("-")[1])]
                            except:
                              print(("TRIAD %s"%triad_obj.return_as_string()))
                              print(("Basepairs %s"%basepairs))
                              exit()

                            dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
                            triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
                            thread_triads_obj.add_triad(triad_obj)
                            #print("Add triad",triad_obj.return_as_string())

                # Reverse strand #
                if dinucleotide_environment[2] == 'R':
                    #Conversion to complementary is unnecesary, R or F is only the bound strand, but sequence is always forward
                    if int(nucleotide_environment[0].split("-")[1]) in list(basepairs.keys()) and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()) and int(nucleotide_environment[2].split("-")[1]) in list(basepairs.keys()):
                            aminoacid_environment = a_oa.split("-")
                            # Information stored in a_oa  
                            #a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                            try:
                              aminoacid_environment[0] = aminoacids1to3[protein[(pdb_chain,int(number))]]
                            except:
                              print(("PDB_CHAIN %s N %d PROTEIN %s"%(pdb_chain,int(number),str(protein))))
                              exit()
                            if aminoacids_polarity_boolean[protein[(pdb_chain,int(number))]]:
                                    hydrophobicity=aminoacid_environment[1] = "P"
                            else:
                                    hydrophobicity=aminoacid_environment[1] = "N"
                            try:
                              #dinucleotide_environment[0] = dna_complementary[basepairs[int(nucleotide_environment[0].split("-")[1])]] + dna_complementary[basepairs[int(nucleotide_environment[2].split("-")[1])]]
                              dinucleotide_environment[0] = basepairs[int(nucleotide_environment[0].split("-")[1])] + basepairs[int(nucleotide_environment[2].split("-")[1])]
                            except:
                              print(("TRIAD %s"%triad_obj.return_as_string()))
                              print(("Basepairs %s"%basepairs))
                              exit()
                            #if not nitrogenous_bases.has_key(dna_complementary[basepairs[int(nucleotide_environment[0].split("-")[1])]]): continue
                            #if not nitrogenous_bases.has_key(dna_complementary[basepairs[int(nucleotide_environment[2].split("-")[1])]]): continue
                            if basepairs[int(nucleotide_environment[0].split("-")[1])] not in nitrogenous_bases: continue
                            if basepairs[int(nucleotide_environment[2].split("-")[1])] not in nitrogenous_bases: continue
                            dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
                            triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
                            thread_triads_obj.add_triad(triad_obj)
                            #print("Add triad",triad_obj.return_as_string())
    return thread_triads_obj, pdb_obj, x3dna_obj


def parse_options():
    """
    Parse command-line options for threading-to-triads conversion.

    How to run:
        python threading_to_triads.py --threading THREAD_FILE
            [-o OUTPUT_TRIADS --pdb PDB_DIR | --template TEMPLATE_PDB]
            [--dummy DUMMY_DIR --clean -v]

    Example:
        python threading_to_triads.py --threading sample_thread.txt --pdb pdb_data -o threaded_triads.dat -v

    The parser configures:
        - Threading input file (`--threading`) and output triads file (`-o`).
        - Structural context source (`--pdb` folder or `--template` file).
        - Runtime controls (`--dummy`, `--clean`, `--verbose`).

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for threading/paths/output.
    """
    parser = optparse.OptionParser(
        "python threading_to_triads.py --threading THREAD_FILE "
        "[-o OUTPUT_TRIADS --pdb PDB_DIR | --template TEMPLATE_PDB] "
        "[--dummy DUMMY_DIR --clean -v]"
    )
    parser.add_option("--threading", action="store", type="string", dest="threading_file",  default=None, help="path to threading file", metavar="THREADING_FILE")
    parser.add_option("-o", default="threaded_triads.dat", action="store", type="string", dest="output", help="Output file (default = threaded_triad.dat)", metavar="{file}")
    parser.add_option("--pdb", action="store", type="string", default="pdb", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("--template", action="store", type="string",  default=None, dest="template", help="Template PDB file", metavar="{file}")
    parser.add_option("--dummy", action="store", type="string",  default="/tmp", dest="dummy_dir", help="Dummy directory (default is /tmp)", metavar="{directory}")
    parser.add_option("-c","--clean",default=False, action="store_true", dest="clean", help="Clean mode to restart the PDB numbering (default = False)")
    parser.add_option("-v","--verbose",default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")



    (options, args) = parser.parse_args()
    return options

def main():
    """
    Run the threading-to-triads conversion workflow.

    Workflow:
        1. Parse command-line options and prepare temporary directories.
        2. Load threading data and resolve structural context from PDB/template input.
        3. Build threaded triads consistent with threaded protein and DNA mappings.
        4. Write the resulting triads file.
        5. Remove temporary files.

    Args:
        None.

    Returns:
        None. The threaded triads file is written to disk.
    """

    # Arguments & Options #
    options  = parse_options()
    dummy_dir=options.dummy_dir
    pdb_dir= options.pdb_dir
    template = options.template
    verbose = options.verbose
    clean = options.clean
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    threading_file = options.threading_file
    triad_obj, pdb_obj, x3dna_obj = threading_triads(threading_file, pdb_dir,  template , verbose , dummy_dir,clean)
    # create the threaded object and get the threading file #
    triad_obj.write(options.output)
    shutil.rmtree(dummy_dir)


if __name__ == "__main__":
    main()
