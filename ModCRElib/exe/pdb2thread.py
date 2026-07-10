import os, sys, re
import configparser
import optparse
import shutil
import subprocess
import gzip
import traceback
import configparser
import argparse
import shutil
import pandas as pd

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

# Import my functions #
from ModCRElib.structure.contacts import contacts
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.threading import threader


# Imports jbonet's module #
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases,dna_complementary
from SBILib.structure import PDB

#-------------#
# Functions   #
#-------------#
      
     
def parse_user_arguments(*args, **kwds):
    """
        Parse command-line arguments for converting PDB complexes to thread files.
    
        How to run:
            python pdb2thread.py -i INPUT_PDB [-o OUTPUT_PREFIX -s DNA_SEQUENCE --dummy DUMMY_DIR -v]
    
        Example:
            python pdb2thread.py -i complex.pdb -o thread_out -s ACGTACGT -v
    
        The parser configures:
            - Input complex PDB (`-i`) and optional output naming (`-o`).
            - Optional DNA sequence override for interface threading (`-s`).
            - Runtime options (`--dummy`, `-v`).
    
        Args:
            *args: Unused positional passthrough.
            **kwds: Unused keyword passthrough.
    
        Returns:
            argparse.Namespace: Parsed command-line options.
    
    Args:
        *args (Any): Value used by this routine.
        **kwds (Any): Value used by this routine."""
    parser = argparse.ArgumentParser(
        usage="python pdb2thread.py -i INPUT_PDB [-o OUTPUT_PREFIX -s DNA_SEQUENCE --dummy DUMMY_DIR -v]"
    )
    parser.add_argument("--dummy", default="/tmp/", action="store",  dest="dummy_dir", help="Dummy directory (default = /tmp/)") 
    parser.add_argument('-i', '--input_file', dest = 'input_file', action = 'store',
                        help = 'PDB file of TF with DNA')
    parser.add_argument('-o', '--output_name', dest = 'output_name', action = 'store', default = None,
                        help = 'Output name to label the THREADING file')
    parser.add_argument("-s", action="store", dest="dna_sequence", default = None,
                        help="DNA sequence (it must have the same length as in PDB file)")
    parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
 
    options = parser.parse_args()
    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def pdb2thread(pdb_obj,name,x3dna_obj,output_dir, dna = None,verbose=False):
    """
        Build `Threaded` objects from a PDB protein-DNA complex.
    
        Args:
            pdb_obj (PDB): Input PDB complex object.
            name (str): Base name used for threading/template files.
            x3dna_obj (X3DNA): Parsed X3DNA annotation object.
            output_dir (str): Output directory for template/thread files.
            dna (str, optional): DNA sequence override for interface threading.
            verbose (bool, optional): Enable verbose logs.
    
        Returns:
            list: Generated `Threaded` objects (typically one per protein chain).
    
    Args:
        pdb_obj (Any): PDB object or structure file input.
        name (Any): Name/identifier used by this routine.
        x3dna_obj (Any): DNA identifier, sequence, or DNA-related data.
        output_dir (Any): Directory path used by this operation.
        dna (Any): DNA identifier, sequence, or DNA-related data.
        verbose (Any): Boolean flag controlling routine behavior."""

    protein_chains=set()
    dna_chains=set()
    paired={}
    for chain_id in pdb_obj.chain_identifiers:
        chain=pdb_obj.get_chain_by_id(chain_id)
        if chain.chaintype=="P": protein_chains.add(chain_id)
        if chain.chaintype=="N": dna_chains.add(chain_id)
    # Get contacts object #
    if verbose: sys.stdout.write("Run PDB2THREAD ...\n")
    if verbose: sys.stdout.write("\t-- Get contacts ...\n")
    contact = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
    for contact_obj in contact._contacts:
        for prot_chain_id in protein_chains:
            for dna_chain_id in dna_chains:
                if (prot_chain_id in contact_obj._A_chain  and dna_chain_id in  contact_obj._B_chain):
                   paired.setdefault(prot_chain_id,(contact_obj._B_chain[0],contact_obj._B_chain[1]))
    thread_objs  = []
    for chain_id,(dna_id,dna_cmpl) in paired.items():
        if verbose: print(("\t-- Threaded object dummy %s, with chains %s (protein) and  %s (dna) named as %s ...."%(pdb_obj.id,chain_id,dna_id,name)))
        threading_pdb_name   = name
        threading_pdb_chain  = chain_id
        threading_dna_helix  = 1
        threading_identity   = 100.0
        threading_coverage   = 100.0
        threading_query_ali  = pdb_obj.get_chain_by_id(chain_id).gapped_protein_sequence
        if verbose: print(("\t-- Sequence Protein: %s"%(threading_query_ali)))
        threading_hit_ali    = pdb_obj.get_chain_by_id(chain_id).gapped_protein_sequence
        threading_dna        ={}
        threading_dnae       ={}
        kmer =  pdb_obj.get_chain_by_id(dna_id).gapped_nucleotide_sequence()
        if verbose: print(("\t-- Sequence DNA: %s (template pdb) vs %s (input interface) "%(kmer,dna)))
        if dna is not None:
           if len(dna) != len(kmer):
              sys.stdout.write("ERROR/WARNING: DNA sequence does not match the protein-DNA interface. Try standard DNA option for automated modifcation if you are using \"get_best_binding\" if the program has not written the thread or option \"-s %s\" if you use \"pdb2thread\" (please replace \"N\" by any nucleotide {A, C, G, T} or your choice or what you expect from the comparison).\n" % ("N" * len(kmer)))
              #raise ValueError("DNA sequence does not match the protein-DNA interface. Try \"-s %s\" (please replace \"N\" by any nucleotide {A, C, G, T} or your choice)." % ("N" * len(kmer)))
           else:
              kmer=dna
        threading_dna.setdefault(kmer,"0")
        threading_dnae.setdefault(kmer,"0")
        threading_protein={}
        positions=[int(n) for n in pdb_obj.get_chain_by_id(chain_id).protein_idx.split(";")]
        for n in range(len(positions)):
            threading_protein.setdefault((chain_id,positions[n]),threading_query_ali[n])
        #Create specific template
        template_pdb = os.path.join(output_dir,threading_pdb_name+"_"+threading_pdb_chain+".pdb")
        if not os.path.exists(template_pdb):
          template= PDB()
          chain   = pdb_obj.get_chain_by_id(chain_id)
          dna_fwd = pdb_obj.get_chain_by_id(dna_id)
          dna_rev = pdb_obj.get_chain_by_id(dna_cmpl)
          template.add_chain(chain)
          template.add_chain(dna_fwd)
          template.add_chain(dna_rev)
          template.write(template_pdb)
        #Create threading object
        threaded_obj = threader.Threaded(None)
        threaded_obj.set_pdb_name(threading_pdb_name)
        threaded_obj.set_pdb_chain(threading_pdb_chain)
        threaded_obj.set_dna_helix(threading_dna_helix)
        threaded_obj.set_identity(threading_identity)
        threaded_obj.set_coverage(threading_coverage)
        threaded_obj.set_protein(threading_protein)
        threaded_obj.set_dna(threading_dna)
        threaded_obj.set_dna_fixed(threading_dnae)
        threaded_obj.set_query_ali(threading_query_ali)
        threaded_obj.set_pdb_ali(threading_hit_ali)
        thread_objs.append(threaded_obj)

    return thread_objs

#-------------#
# Main        #
#-------------#

def main():
    """
    Run the PDB-to-threading conversion workflow.

    Workflow:
        1. Parse command-line options and prepare temporary directories.
        2. Load the input protein-DNA complex and derive X3DNA annotations.
        3. Build per-chain `Threaded` objects from contacts and interface DNA.
        4. Write template PDB files and threading output files.
        5. Clean temporary files and finish execution.

    Args:
        None.

    Returns:
        None. Template and threading files are written to disk.
    """

    # Arguments & Options #
    options      = parse_user_arguments()
    dummy_dir    = options.dummy_dir
    pdb_file     = options.input_file
    output_name  = options.output_name
    verbose      = options.verbose
    dna_sequence = options.dna_sequence
    if not os.path.exists(dummy_dir):  os.mkdir(dummy_dir)

    # Get PDB object #
    if verbose: print(("Open PDB file: %s"%(os.path.basename(pdb_file))))
    pdb_obj    = PDB(pdb_file)
    name       = pdb_obj.id
    output_dir = os.path.dirname(pdb_file)
    # Get X3DNA object #
    if verbose: sys.stdout.write("\t-- Get DNA object ...\n")
    x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath(pdb_file), dummy_dir)
    if verbose: sys.stdout.write("\t-- Thread the PDB file ...\n")
    try:
       thread_objs = pdb2thread(pdb_obj,name,x3dna_obj,output_dir, dna_sequence,verbose)
    except Exception as e:
       print(("Error: %s"%e))
       exit()
    """
    for threaded_obj in thread_objs:
        name_pdb = threaded_obj.get_pdb_name()
        chain_id = threaded_obj.get_pdb_chain()
        if output_name is not None:
           output_file = os.path.join(output_name+"_"+name_pdb+"_"+chain_id)
        else:
           output_file = os.path.join(name_pdb+"_"+chain_id)
        if dna_sequence is not None:
           output_file = output_file +"_"+dna_sequence
        output_file = output_file +".txt"
        if verbose: sys.stdout.write("\t-- Write threading file %s\n"%(output_file))
        threaded_obj.write(output_file)
    """
    # Write PDB templates and threads with the best DNA sequence
    if verbose: sys.stdout.write("Write threading files\n")
    try:
      input_dir=os.path.dirname(pdb_file)
      if output_name is not None:
       output_dir=os.path.dirname(output_name)
       output_file=os.path.basename(output_name)
      else:
       output_dir=os.path.abspath("./")
       output_file=thread_objs[0].get_pdb_name()
      threading_pdb_name = output_file
      output_complex = os.path.join(output_dir,threading_pdb_name+"_")
      thread_complex = []
      n=0
      for thr_obj in thread_objs:
        n=n+1
        threading_pdb_name  = thr_obj.get_pdb_name()
        threading_pdb_chain = thr_obj.get_pdb_chain()
        template_pdb        = os.path.join(input_dir,threading_pdb_name+"_"+threading_pdb_chain+".pdb")
        if verbose: sys.stdout.write("\t-- Write template file %s\n"%(os.path.basename(template_pdb)))
        if not os.path.exists(os.path.join(output_dir,os.path.basename(template_pdb))): 
               shutil.copy(template_pdb,os.path.join(output_dir,os.path.basename(template_pdb)))
        output_file = os.path.join(output_dir,threading_pdb_name+"_"+threading_pdb_chain)
        output_complex = output_complex + threading_pdb_chain
        kmers=thr_obj.get_kmers()
        for kmr in kmers.keys():
            dna_seq = kmr
        output_file = output_file+"_"+dna_seq+".txt"
        if os.path.exists(output_file):
            if verbose: print(("\t-- Reuse threading file %s\n"%(output_file)))
            thr_obj.write(output_file)
        else:
            if verbose: sys.stdout.write("\t-- Write threading file %s\n"%(output_file))
            thr_obj.write(output_file)
        thread_complex.append(os.path.abspath(output_file))
      if n>1:
        output_complex = output_complex+"_"+dna_seq+".txt"
        if verbose: sys.stdout.write("\n\t-- Write thread complex  %s\n"%(output_complex))
        fo=open(output_complex,"w")
        for output_file in thread_complex:
          fo.write("%s\n"%output_file)
        fo.close()
    except Exception as e:
      print(("Error: %s"%e))
      exit()



    # Clean files #
    try:
       shutil.rmtree(dummy_dir)
    except:
       print("Remaining DUMMY files, not cleanded")

    if verbose: print("Done")


if __name__ == "__main__":
    main()
