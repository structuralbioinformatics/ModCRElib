import os, sys, re
import configparser
import optparse
import shutil
import subprocess
import socket
import time
from Bio.PDB import PDBIO, Superimposer, PDBParser

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


# Get python path #
python = os.path.join(config.get("Paths", "python_path"), "python")


# Imports jbonet's module #
from SBILib.structure import PDB,Chain
from SBILib.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBILib.structure.residue import ResidueOfNucleotide
from SBILib.structure.atom import AtomOfNucleotide


#alphabet for chains
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWZabcdefghijklmnopqrstuvwxyz1234567890"
alphabet2= []

#for x in alphabet:
#    alphabet2.append(x)
for x in alphabet:
    alphabet2.append(x)
for x in alphabet:
    for y in alphabet:
       alphabet2.append(x+y)
for x in alphabet:
    for y in alphabet:
       for z in alphabet:
           alphabet2.append(x+y+z)


############
# FUNCTIONS
############

def parse_options():
    """
    Parse command-line options for DNA remodeling/extension workflows.

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for input/output and modeling setup.

    """

    parser = optparse.OptionParser("python modify_dna.py -i input_file  [  -o output_file -n EXTENSION --dummy=dummy_dir --verbose --conformation A|B|C|D|Z  --extend ] ")

    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file is a PDB file or a folder of PDB files", metavar="{filename}")
    parser.add_option("-o", default="dna.pdb", action="store", type="string", dest="output_file", help="Output file (default = dna.pdb) or new folder (PDB files named as in the input) ", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode. If not selected the dummy directory will be removed (default = False)", metavar="{boolean}")
    parser.add_option("-e", "--extend", default=False, action="store_true", dest="extend", help="Extend or substitute the DNA. If not selected the DNA will be remodelled, otherwise it will extend (default = False)", metavar="{boolean}")
    parser.add_option("-n",default=10,action="store", type="int",dest="extension", help="Number of base-pairs to extend the length of te DNA at 5' and 3' (default is 10)", metavar="{integer}")
    parser.add_option("-N",default="A",action="store", type="string",dest="nucleotide", help="Dummy Nucleotide to extend the sequence (default is A)", metavar="{string}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("--conformation", default="B", action="store", type="string", dest="conformation", help="Type of conformation to extend or remodel the DNA (only accept A, B, C, D, and Z) (default is B)", metavar="{string}")
    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def dna_superimposition_obj(pdb_obj_static, pdb_obj_mobile, extension,  dummy_dir="/tmp", verbose=False):
    """
        Superimpose a mobile TF-DNA complex onto a static DNA scaffold.
    
        Args:
            pdb_obj_static (PDB): Static DNA-rich structure used as alignment target.
            pdb_obj_mobile (PDB): Mobile protein-DNA structure to superimpose.
            extension (int): Number of basepairs added to each DNA end.
            dummy_dir (str, optional): Temporary directory root.
            verbose (bool, optional): Enable verbose logging.
    
        Returns:
            PDB: Merged PDB object containing static DNA and superimposed proteins.
    
    Args:
        pdb_obj_static (Any): Value used by this routine.
        pdb_obj_mobile (Any): Value used by this routine.
        extension (Any): Value used by this routine.
        dummy_dir (Any): Directory path used by this operation.
        verbose (Any): Boolean flag controlling routine behavior."""

    # Initialize #

    #Temporary folder
    tmp = os.path.join(dummy_dir, str(os.getpid()))
    if not os.path.exists(tmp): os.makedirs(tmp)
    dummy_file = os.path.join(tmp,"superimposition.pdb")

    #Define variables for alignment
    alignment_A = set()
    alignment_B = set()
    backbone_atoms = set(["P", "O1P", "O2P", "O3P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 
                               "OP1", "OP2", "OP3", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*"])
    a_chains=[]
    for c in pdb_obj_static.chains:
     identifier = c.chain
     if identifier not in a_chains and c.chaintype=="N": 
        a_chains.append(identifier)
    b_chains=[]
    for c in pdb_obj_mobile.chains:
     identifier = c.chain
     if identifier not in b_chains and c.chaintype=="N": 
        b_chains.append(identifier)
    sequence_A  = pdb_obj_static.get_chain_by_id(a_chains[0]).gapped_nucleotide_sequence()
    sequence_B  = pdb_obj_mobile.get_chain_by_id(b_chains[0]).gapped_nucleotide_sequence()
    for i in range(len(sequence_B)):
        alignment_B.add((b_chains[0],i+1))
        alignment_A.add((a_chains[0],i+extension))
    #Make pdb files of DNA and TFs
    pdb_file_A=os.path.join(tmp,"static.pdb")
    pdb_file_B=os.path.join(tmp,"mobile.pdb")
    pdb_obj_static.write(pdb_file_A)
    pdb_obj_mobile.write(pdb_file_B)

    # Start the parser #
    pdb_parser = PDBParser(QUIET=True) 
    # Load the PDB files #
    if verbose: print(("Read PDB %s"%(pdb_file_A)))
    static_pdb = pdb_parser.get_structure("static", pdb_file_A)
    if verbose: print(("Read PDB %s"%(pdb_file_B)))
    mobile_pdb = pdb_parser.get_structure("mobile", pdb_file_B)
    # Use the first model in each PDB #
    static_model = static_pdb[0]
    mobile_model = mobile_pdb[0]
    # Get the residues to be aligned #
    static_residues = []
    mobile_residues = []
    # For each chain... #
    for chain in static_model:
        # For each residue... #
        for residue in chain:
            # Initialize #
            residue_num = residue.get_id()
            # If residue and chain in alignment... #
            if (chain.get_id(), residue_num[1]) in alignment_A:
                static_residues.append(residue)
    # For each chain... #
    for chain in mobile_model:
        # For each residue... #
        for residue in chain:
            # Initialize #
            residue_num = residue.get_id()
            # If residue and chain in alignment... #
            if (chain.get_id(), residue_num[1]) in alignment_B:
                mobile_residues.append(residue)

    # Get the atoms to be aligned #
    static_atoms = []
    mobile_atoms = []
    # For each residue... #
    for i in range(len(static_residues)):
        static_residue = static_residues[i]
        mobile_residue = mobile_residues[i]
        # For each backbone atom... #
        for backbone_atom in backbone_atoms:
            # Initialize #
            atom_in_static_residue = False
            atom_in_mobile_residue = False
            # For each atom... #
            for atom in static_residue:
                if atom.name == backbone_atom:
                    atom_in_static_residue = True
                    break
            # For each atom... #
            for atom in mobile_residue:
                if atom.name == backbone_atom:
                    atom_in_mobile_residue = True
                    break
            # If both residues have backbone atom... #
            if atom_in_static_residue and atom_in_mobile_residue:
                static_atoms.append(static_residue[backbone_atom])
                mobile_atoms.append(mobile_residue[backbone_atom])
    # Initiate the superimposer #
    if verbose: print(("Superimposing in %s"%(dummy_file)))
    superimposer = Superimposer() 
    superimposer.set_atoms(static_atoms, mobile_atoms)
    superimposer.apply(mobile_model.get_atoms()) 
    # Save the superposed structure #
    io = PDBIO()
    io.set_structure(mobile_pdb)
    io.save(dummy_file)

    # Get PDB object #
    pdb_superimposed_obj = PDB(dummy_file)
    pdb_obj=PDB()
    number_of_chains = 0
    for chain in pdb_obj_static.chains:
        if chain.chaintype == "N": 
           chain_id = alphabet2[number_of_chains]
           number_of_chains = number_of_chains + 1
           chain.chain=chain_id
           pdb_obj.add_chain(chain)
    for chain in pdb_superimposed_obj.chains:
        if chain.chaintype == "P": 
           chain_id = alphabet2[number_of_chains]
           number_of_chains = number_of_chains + 1
           chain.chain=chain_id
           pdb_obj.add_chain(chain)

    #Clean temporary folder
    if os.path.exists(tmp) : shutil.rmtree(tmp)

    return pdb_obj

def merge_dna(a,b):
   """
      Merge two DNA fragments into canonical `A`/`B` nucleotide chains.
   
      Args:
          a (PDB): First DNA fragment object.
          b (PDB): Second DNA fragment object.
   
      Returns:
          PDB: Combined DNA PDB object with renumbered residues.
   
   Args:
       a (Any): Value used by this routine.
       b (Any): Value used by this routine."""
   a_chains=[]
   for c in a.chains:
    identifier = c.chain
    if identifier not in a_chains: 
       a_chains.append(identifier)

   b_chains=[]
   for c in b.chains:
    identifier = c.chain
    if identifier not in b_chains: 
       b_chains.append(identifier)

   A=ChainOfNucleotide("Complex","A")
   B=ChainOfNucleotide("Complex","B")

   number = 1
   chain_id = a_chains[0]
   c= a.get_chain_by_id(chain_id)

   for nn in c.nucleotides:
    if nn != c.last_nucleotide:
       nn.number = number
       number    = number + 1
       A.add_residue(nn)

   chain_id = b_chains[0]
   c= b.get_chain_by_id(chain_id)

   for nn in c.nucleotides:
       nn.number = number
       number    = number + 1
       A.add_residue(nn)

   chain_id = b_chains[1]
   c= b.get_chain_by_id(chain_id)

   for nn in c.nucleotides:
       nn.number = number
       number    = number + 1
       B.add_residue(nn)


   chain_id = a_chains[1]
   c= a.get_chain_by_id(chain_id)

   for nn in c.nucleotides:
    if nn != c.first_nucleotide:
       nn.number = number
       number    = number + 1
       B.add_residue(nn)



   p=PDB()
   p.add_chain(A)
   p.add_chain(B)

   return p

def remodel_dna_obj(dna_str,sequence,extension, nucleotide_extend,dummy_dir="/tmp", verbose=False):
    """
        Build a remodeled DNA-only structure with symmetric sequence extension.
    
        Args:
            dna_str (str): DNA conformation code (`A`, `B`, `C`, `D`, `Z`).
            sequence (str): Core DNA sequence to preserve.
            extension (int): Number of nucleotides to add at both ends.
            nucleotide_extend (str): Filler nucleotide used in extensions.
            dummy_dir (str, optional): Temporary directory root.
            verbose (bool, optional): Enable verbose logging.
    
        Returns:
            PDB: Remodeled DNA structure as PDB object.
    
    Args:
        dna_str (Any): DNA identifier, sequence, or DNA-related data.
        sequence (Any): Sequence/alignment content used in this step.
        extension (Any): Value used by this routine.
        nucleotide_extend (Any): Value used by this routine.
        dummy_dir (Any): Directory path used by this operation.
        verbose (Any): Boolean flag controlling routine behavior."""

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        os.environ['X3DNA'] = x3dna_path[:-4]
        if verbose: print("Build_dna_structure ")
        #if verbose: print("X3DNA =  " + x3dna_path[:-4])
        # Get current working directory #
        cwd = os.getcwd()
        # Create tmp directory #
        tmp = os.path.join(dummy_dir, str(os.getpid()))
        if not os.path.exists(tmp): os.makedirs(tmp)
        # Change directory #
        os.chdir(tmp)
        # Define a dummy pdb file #
        dna_file=os.path.abspath("dummy.pdb")
        # Define the full DNA sequence
        modeling_seq= extension*nucleotide_extend+sequence+extension*nucleotide_extend

        # Exec process #
        if dna_str == "B":
            if verbose: print(("X3DNA ("+tmp+")  command: " + os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file+" > model_dna.log")

        elif dna_str == "A":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -a " + dna_file+" > model_dna.log"))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-a", dna_file)
    
        elif dna_str == "C":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -c " + dna_file+" > model_dna.log"))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-c", dna_file)

        elif dna_str == "D":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -d " + dna_file+" > model_dna.log"))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-d", dna_file)

        elif dna_str == "Z":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -z " + dna_file+" > model_dna.log"))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-z", dna_file)

        else:
            #by default use B dna straigth line
            if verbose: print(("X3DNA ("+tmp+")  command: " + os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file+" > model_dna.log")
        # Return to original directory #
        os.chdir(cwd)

        #Get DNA PDB obj and return
        pdb_obj=PDB(dna_file)
        # Erase tmp directory #
        if os.path.exists(tmp) : shutil.rmtree(tmp)

    except Exception as e:
        raise ValueError("Could not exec X3DNA REMODELLING with error %s" % (e))


    return pdb_obj


def extended_dna_obj(dna_str, extension, nucleotide_extend, nucleotide, position,  dummy_dir="/tmp", verbose=False):
    """
        Build one 5' or 3' DNA extension fragment and orient it to a reference frame.
    
        Args:
            dna_str (str): DNA conformation code (`A`, `B`, `C`, `D`, `Z`).
            extension (int): Extension length in nucleotides.
            nucleotide_extend (str): Filler nucleotide for extension segment.
            nucleotide (str): Terminal nucleotide used at extension edge.
            position (str): Extension side (`5'` or `3'`).
            dummy_dir (str, optional): Temporary directory root.
            verbose (bool, optional): Enable verbose logging.
    
        Returns:
            PDB: Oriented extension fragment as PDB object.
    
    Args:
        dna_str (Any): DNA identifier, sequence, or DNA-related data.
        extension (Any): Value used by this routine.
        nucleotide_extend (Any): Value used by this routine.
        nucleotide (Any): Value used by this routine.
        position (Any): Value used by this routine.
        dummy_dir (Any): Directory path used by this operation.
        verbose (Any): Boolean flag controlling routine behavior."""

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        os.environ['X3DNA'] = x3dna_path[:-4]
        if verbose: print("Build_dna_structure ")
        #if verbose: print("X3DNA =  " + x3dna_path[:-4])
        if    position == "5'":
              label = "5_"
              offset=extension
              modeling_seq= (extension-1)*nucleotide_extend+nucleotide
        elif  position == "3'":
              label = "3_"
              offset=1
              modeling_seq= nucleotide+(extension-1)*nucleotide_extend
        else:
              label = "5_"
              offset=extension
              modeling_seq= (extension-1)*nucleotide_extend+nucleotide


        # Get current working directory #
        cwd = os.getcwd()
        # Create tmp directory #
        tmp = os.path.join(dummy_dir, str(os.getpid()))
        if not os.path.exists(tmp): os.makedirs(tmp)
        # Change directory #
        os.chdir(tmp)
        # Define a dummy pdb file #
        dna_file=os.path.abspath("dummy.pdb")

        # Exec process #

        if dna_str == "B":
            if verbose: print(("X3DNA ("+tmp+")  command: " + os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file+" > model_dna.log")

        elif dna_str == "A":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -a " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-a", dna_file+" > model_dna.log")
    
        elif dna_str == "C":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -c " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-c", dna_file+" > model_dna.log")

        elif dna_str == "D":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -d " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-d", dna_file+" > model_dna.log")

        elif dna_str == "Z":
            if verbose: print((os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -z " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-z", dna_file+" > model_dna.log")

        else:
            #by default use B dna straigth line
            if verbose: print(("X3DNA ("+tmp+")  command: " + os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file))
            os.system(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file+" > model_dna.log")

        
        # Return to original directory #
        os.chdir(cwd)
        #Get DNA PDB obj and return
        pdb_obj=get_reference_frame(dna_file,offset, label ,dummy_dir,verbose)
        # Erase tmp directory #
        if os.path.exists(tmp) : shutil.rmtree(tmp)

    except Exception as e:
        raise ValueError("Could not exec X3DNA EXTENSION with error %s" % (e))


    return pdb_obj



def get_reference_frame(pdb_file,offset=1, label=None, dummy_dir="/tmp",verbose=False):
    """
        Build a DNA reference frame using X3DNA `find_pair` and `frame_mol`.
    
        Args:
            pdb_file (str): Input PDB file path.
            offset (int, optional): Basepair offset used by `frame_mol`.
            label (str, optional): Prefix for generated reference-frame files.
            dummy_dir (str, optional): Temporary directory root.
            verbose (bool, optional): Enable verbose logging.
    
        Returns:
            PDB: Re-oriented PDB object in the selected reference frame.
    
    Args:
        pdb_file (Any): Path to the input/output file.
        offset (Any): Value used by this routine.
        label (Any): Name/identifier used by this routine.
        dummy_dir (Any): Directory path used by this operation.
        verbose (Any): Boolean flag controlling routine behavior."""
    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        os.environ['X3DNA'] = x3dna_path[:-4]
        # Get current working directory #
        cwd = os.getcwd()
        # Create tmp directory #
        tmp =  os.path.join(dummy_dir, str(os.getpid()))
        if not os.path.exists(tmp): os.makedirs(tmp)
        # Change directory #
        os.chdir(tmp)
        #Output
        if label is not None: 
           reference_frame = label +"_reference_frame"
        else:
           reference_frame = "reference_frame"
        # Exec process #
        if verbose: print(("X3DNA ("+tmp+") command  "+os.path.join(x3dna_path, "find_pair")+" "+pdb_file+" 3dna.out"))
        os.system("\cp %s . "%(pdb_file))
        os.system("%s %s %s > find_pair.log "%(os.path.join(x3dna_path, "find_pair"), os.path.basename(pdb_file), "3dna.out"))
        if verbose: print(("X3DNA ("+tmp+")  command  "+os.path.join(x3dna_path, "frame_mol")+" "+"%d"%(-offset)+" ref_frames.dat "+pdb_file+" "+reference_frame+".pdb"))
        os.system("%s %s %s %s %s > frame_mol.log "%(os.path.join(x3dna_path, "frame_mol"),"%d"%(-offset),"ref_frames.dat",os.path.basename(pdb_file),reference_frame+".pdb"))
        # Get FRAME PDB object #
        pdb_obj = PDB(reference_frame+".pdb")
        # Return to original directory #
        os.chdir(cwd)
        # Erase tmp directory #
        if os.path.exists(tmp) : shutil.rmtree(tmp)
    except Exception as e:
        raise ValueError("Could not exec X3DNA REFERENCE FRAME for %s with error %s" % (pdb_file,e))

    return pdb_obj


def modify_dna(pdb_file, output_file, dna_str, extension, nucleotide, extend, dummy_dir, verbose):
    """
        Remodel or extend DNA in a complex and superimpose original proteins.
    
        Args:
            pdb_file (str): Input PDB file path.
            output_file (str): Output PDB path.
            dna_str (str): DNA conformation code (`A`, `B`, `C`, `D`, `Z`).
            extension (int): Number of nucleotides/basepairs to extend.
            nucleotide (str): Filler nucleotide used for extension/remodeling.
            extend (bool): If True, extend DNA; otherwise remodel in-place.
            dummy_dir (str): Temporary directory root.
            verbose (bool): Enable verbose logging.
    
        Returns:
            None. The modified structure is written to `output_file`.
    
    Args:
        pdb_file (Any): Path to the input/output file.
        output_file (Any): Path to the input/output file.
        dna_str (Any): DNA identifier, sequence, or DNA-related data.
        extension (Any): Value used by this routine.
        nucleotide (Any): Value used by this routine.
        extend (Any): Value used by this routine.
        dummy_dir (Any): Directory path used by this operation.
        verbose (Any): Boolean flag controlling routine behavior."""

    #Get DNA chains and DNA sequence forward
    pdb_obj=PDB(pdb_file)
    dna_chains=[]
    for chain in pdb_obj.chains:
      if chain.chaintype == "N":
       dna_chains.append(chain.chain)
    sequence=pdb_obj.get_chain_by_id(dna_chains[0]).gapped_nucleotide_sequence()

    

    #Make the DNA extended or remodelled
    if not extend:
       try:
          pdb_dna_obj= remodel_dna_obj(dna_str, sequence, extension, nucleotide, tmp, verbose)
       except Exception as e:
          print(("Error %s"%e))
          raise ValueError("Error remodeling DNA %s"%e)

    else:
       # Run extension of DNA
       try:
          a=extended_dna_obj(dna_str, extension, nucleotide,sequence[0] ,"5'", tmp, verbose)
          b=get_reference_frame(pdb_file,1,"model",tmp, verbose)
          p=merge_dna(a,b)
          intermediate=os.path.abspath(os.path.join(tmp,"extended_5_model.pdb"))
          p.write(intermediate,force=True)
          c=get_reference_frame(intermediate,len(sequence)+extension-1,"extended_5_model",tmp,verbose)
          d=extended_dna_obj(dna_str, extension, nucleotide,sequence[-1] ,"3'", tmp, verbose)
          pdb_dna_obj=merge_dna(c,d)
       except Exception as e:
          print(("Error on file %s is %s"%(os.path.basename(pdb_file),e)))
          raise ValueError("Error on file %s is %s"%(os.path.basename(pdb_file),e))

    #Superimpose the TF
    if verbose: print(("Superimposition of %s to dummy "%(pdb_file)))
    result_obj=dna_superimposition_obj(pdb_dna_obj,pdb_obj,extension,tmp,verbose)


    #Write DNA output file
    try:
       if not output_file.endswith("pdb"): output_file=output_file+".pdb"
       if verbose: print(("Write output PDB %s "%(output_file)))
       result_obj.write(output_file,force=True)
    except:
       print(("Error on writing file %s is %s"%(os.path.basename(output_file),e)))
       raise ValueError("Error on writing file %s is %s"%(os.path.basename(output_file),e))

    # Erase tmp directory #
    if os.path.exists(tmp) : shutil.rmtree(tmp)



def main():
    """
    Run the command-line DNA modification workflow.

    Workflow:
        1. Parse command-line options.
        2. Resolve input as one PDB file or a folder of PDB files.
        3. Remodel/extend DNA and superimpose proteins for each input.
        4. Write output PDB files.

    Args:
        None.

    Returns:
        None. Modified PDB files are written to disk.
    """

    # Arguments & Options #
    options = parse_options()
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    tmp = dummy_dir
    if not os.path.exists(tmp): os.makedirs(tmp)
    if options.input_file  is not None:
     input_file = options.input_file
     if not input_file.startswith("/"): input_file = os.path.abspath(options.input_file)
    output      = options.output_file
    dna_str     = options.conformation
    extension   = options.extension
    nucleotide  = options.nucleotide
    if nucleotide not in list("ACGTU"):  nucleotide="A"
    extend      = options.extend
    verbose     = options.verbose


    if os.path.isfile(input_file):
       pdb_file=input_file
       output_file=output
       if verbose: print(("Modify %s to %s"%(pdb_file,output_file)))
       try:
        modify_dna(pdb_file, output_file, dna_str, extension, nucleotide, extend, dummy_dir, verbose)
       except Exception as e:
        print(("Error %s"%e))

    if os.path.isdir(input_file):
       if not os.path.exists(output): os.makedirs(output)  
       pdb_files = [x for x in os.listdir(input_file) if x.endswith(".pdb")]
       print(pdb_files)
       for x in pdb_files:
           pdb_file    = os.path.join(input_file,x)
           output_file = os.path.join(output,x)
           if verbose: print(("Modify %s to %s"%(pdb_file,output_file)))
           try:
             modify_dna(pdb_file, output_file, dna_str, extension, nucleotide, extend, dummy_dir, verbose)
           except Exception as e:
             print(("Error %s"%e))

    print("Done")


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
