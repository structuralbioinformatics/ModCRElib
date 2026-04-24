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

# Import my modules #
from ModCRElib.structure.protein import model_protein
from ModCRElib.structure.dna import x3dna
from ModCRElib.structure.contacts import contacts,interface
from ModCRElib.sequence import homologs as HOMO



# Import jbonet's module #
from SBILib.structure import PDB


class Threaded(object):
    """
    Store, parse, and write threading assignments for protein-DNA templates.

    Object features:
        - Template metadata (`_pdb_name`, `_pdb_chain`, `_helix`).
        - Alignment/quality metadata (`_identity`, `_coverage`, `_pdb_ali`,
          `_query_ali`).
        - Threaded protein map keyed by `(chain, residue_number)`.
        - Threaded DNA maps by k-mer and binding start (`_dna`, `_dna_fixed`).
        - Expanded basepair-to-nucleotide lookup tables for both DNA maps.

    Example serialized threading object:
        # PDB name    : 1ABC
        # PDB chain   : A
        # DNA helix   : 1
        # % identity  : 87.5
        # % coverage  : 93.2
        # PDB align   : MKT--LAV
        # Query align : MRSQQLAV
        >protein
        A;10;M
        A;11;R
        A;12;S
        A;13;L
        A;14;A
        A;15;V
        >dna
        ACGTAC;0
        >dna_fixed
        ACGTAC;0
        //
    """

    def __init__(self, threading_file=None):
        """
                Initialize a threading container and optionally parse a file.
        
                Args:
                    threading_file (str, optional): Path to a threading file.
        
                Returns:
                    None.
        
        Args:
            threading_file (Any): Path to the input/output file."""
        if threading_file is not None: 
           self._file = threading_file
           self._check_parsing()
        else:
           self._file = None
        self._pdb_name = None
        self._pdb_chain = None
        self._helix = None
        self._identity = None
        self._coverage = None
        self._pdb_ali = None
        self._query_ali = None
        self._protein = None
        self._dna = None 
        self._dna_fixed = None 
        self._basepairs = None
        self._basepairs_fixed = None

    def _check_parsing(self):
        """Initialize internal maps and parse file when needed."""
        if self._file is not None:
            self._protein = {}
            self._dna = {}
            self._dna_fixed = {}
            self._basepairs = {}
            self._basepairs_fixed = {}
            self._parse_file()

    def _get_file(self):
        """Return source threading file path."""
        return self._file

    def _parse_file(self):
        """Parse threading file blocks into protein/DNA assignment maps."""
        if os.path.exists(self._get_file()):
            threaded = None
            if self._get_file().endswith(".gz"):
                fd = gzip.open(self._get_file(), "rt")
            else:
                fd = open(self._get_file(), "rt")
            for line in fd:
                # Get helix #
                m = re.search("^# PDB name    : (\S+)", line.strip())
                if m:
                    self._pdb_name = m.group(1)
                m = re.search("^# PDB chain   : (\S)", line.strip())
                if m:
                    self._pdb_chain = m.group(1)
                m = re.search("^# DNA helix   : (\d+)", line.strip())
                if m:
                    self._helix = m.group(1)
                m = re.search("^# %% identity  : (\S+)", line.strip())
                if m:
                    self._identity = float(m.group(1))
                m = re.search("^# %% coverage  : (\S+)", line.strip())
                if m:
                    self._coverage = float(m.group(1))
                m = re.search("^# PDB align   : (\S+)", line.strip())
                if m:
                    self._pdb_ali = m.group(1)
                m = re.search("^# Query align : (\S+)", line.strip())
                if m:
                    self._query_ali = m.group(1)
                m = re.search("^>(\w+)", line.strip())
                if m:
                    threaded = m.group(1)
                if threaded != None:
                    # Skip incorrect lines #
                    if line.startswith("#"): continue
                    if line.startswith(">"): continue
                    if line.startswith("//"): continue
                    word = line.strip().split(";")
                    if threaded == "protein":
                        self._protein.setdefault((word[0], int(word[1])), word[2])
                    elif threaded == "dna":
                        self._dna.setdefault(word[0], word[1])
                        # Get basepairs #
                        basepairs = list(range(int(word[1]), int(word[1]) + len(word[0])))
                        for nucleotide in word[0]:
                            self._basepairs.setdefault((word[0], basepairs.pop(0)), nucleotide)
                    elif threaded == "dna_fixed":
                        self._dna_fixed.setdefault(word[0], word[1])
                        # Get basepairs #
                        basepairs = list(range(int(word[1]), int(word[1]) + len(word[0])))
                        for nucleotide in word[0]:
                            self._basepairs_fixed.setdefault((word[0], basepairs.pop(0)), nucleotide)
            fd.close()
        else:
            raise ValueError("Could not open threading file %s" % self._get_file())

    def get_pdb_name(self):
        """Return template PDB identifier."""
        if self._pdb_name is None: self._check_parsing()
        return self._pdb_name

    def set_pdb_name(self,pdb_name):
        """Set template PDB identifier.
        
        Args:
            pdb_name (Any): Value used by this routine."""
        self._pdb_name=pdb_name

    def get_pdb_chain(self):
        """Return threaded template protein chain identifier."""
        if self._pdb_chain is None: self._check_parsing()
        return self._pdb_chain

    def set_pdb_chain(self,pdb_chain):
        """Set threaded template protein chain identifier.
        
        Args:
            pdb_chain (Any): PDB structure identifier and chain selector."""
        self._pdb_chain=pdb_chain
    
    def get_dna_helix(self):
        """Return DNA helix identifier associated with the thread."""
        if self._helix is None: self._check_parsing()
        return self._helix

    def set_dna_helix(self,dna_helix):
        """Set DNA helix identifier for this thread.
        
        Args:
            dna_helix (Any): DNA identifier, sequence, or DNA-related data."""
        self._helix=dna_helix

    def get_score(self):
        """Return combined score as identity × coverage."""
        if self._identity is None or self._coverage is None: self._check_parsing()
        return self._identity * self._coverage

    def set_identity(self,identity):
        """Set identity metric.
        
        Args:
            identity (Any): Value used by this routine."""
        self._identity=identity

    def set_coverage(self,coverage):
        """Set coverage metric.
        
        Args:
            coverage (Any): Value used by this routine."""
        self._coverage=coverage

    def set_pdb_ali(self,pdb_ali):
        """Set template alignment string.
        
        Args:
            pdb_ali (Any): Value used by this routine."""
        self._pdb_ali=pdb_ali

    def set_query_ali(self,query_ali):
        """Set query alignment string.
        
        Args:
            query_ali (Any): Value used by this routine."""
        self._query_ali=query_ali

    def set_protein(self,protein):
        """Set threaded protein map keyed by `(chain, residue_number)`.
        
        Args:
            protein (Any): Protein identifier, object, or protein-related data."""
        self._protein = protein

    def get_protein(self):
        """Return threaded protein map."""
        if self._protein is None: self._check_parsing()
        return self._protein

    def set_dna(self,dna):
        """Set k-mer DNA assignments and rebuild basepair lookup map.
        
        Args:
            dna (Any): DNA identifier, sequence, or DNA-related data."""
        self._basepairs = {}
        self._dna = dna
        for seq,binding in dna.items():
          basepairs = list(range(int(binding), int(binding) + len(seq)))
          for nucleotide in seq:
            self._basepairs.setdefault((seq, basepairs.pop(0)), nucleotide)

    def set_dna_fixed(self,dna_fixed):
        """Set corrected/fixed DNA assignments and rebuild basepair map.
        
        Args:
            dna_fixed (Any): DNA identifier, sequence, or DNA-related data."""
        self._basepairs_fixed = {}
        self._dna_fixed = dna_fixed
        for seq,binding in dna_fixed.items():
          basepairs_fixed = list(range(int(binding), int(binding) + len(seq)))
          for nucleotide in seq:
            self._basepairs_fixed.setdefault((seq, basepairs_fixed.pop(0)), nucleotide)

    def get_kmers(self):
        """Return threaded k-mer assignments."""
        if self._dna is None or self._dna =={}: self._check_parsing()
        return self._dna

    def get_kmers_fixed(self):
        """Return fixed threaded k-mer assignments."""
        if self._dna_fixed is None or self._dna_fixed=={}:self._check_parsing()
        return self._dna_fixed

    def has_chain(self, chain):
        """Return whether input chain matches threaded chain.
        
        Args:
            chain (Any): Chain identifier."""
        if self._pdb_chain is None: self._check_parsing()
        return chain == self._pdb_chain    

    def has_aminoacid(self, chain, residue_num):
        """Return whether a threaded amino acid exists for chain/position.
        
        Args:
            chain (Any): Chain identifier.
            residue_num (Any): Numeric identifier/index used by this routine."""
        if self._protein is None: self._check_parsing()
        return (chain, residue_num) in iter(self._protein.keys())

    def get_threaded_protein_sequence(self, chain):
        """Return threaded protein sequence for one chain.
        
        Args:
            chain (Any): Chain identifier."""
        if self._protein is None: self._check_parsing()
        if self.has_chain(chain):
            return "".join([self._protein[(key)] for key in sorted(self._protein, key=lambda x: x[-1]) if chain == key[0]])
        return None

    def get_threaded_aminoacid(self, chain, residue_num):
        """Return threaded amino acid at one chain/residue position.
        
        Args:
            chain (Any): Chain identifier.
            residue_num (Any): Numeric identifier/index used by this routine."""
        if self._protein is None: self._check_parsing()
        if self.has_aminoacid(chain, residue_num):
            return self._protein[(chain, residue_num)]
        return None

    def has_basepair(self, kmer, basepair):
        """Return whether a basepair mapping exists in threaded DNA.
        
        Args:
            kmer (Any): Value used by this routine.
            basepair (Any): Value used by this routine."""
        if self._basepairs is None: self._check_parsing()
        return (kmer, basepair) in self._basepairs

    def has_basepair_fixed(self, kmer, basepair):
        """Return whether a basepair mapping exists in fixed threaded DNA.
        
        Args:
            kmer (Any): Value used by this routine.
            basepair (Any): Value used by this routine."""
        if self._basepairs_fixed is None: self._check_parsing()
        return (kmer, basepair) in self._basepairs_fixed

    def get_threaded_nucleotide(self, kmer, basepair):
        """Return threaded nucleotide for one k-mer/basepair coordinate.
        
        Args:
            kmer (Any): Value used by this routine.
            basepair (Any): Value used by this routine."""
        if self._basepairs is None: self._check_parsing()
        if self.has_basepair(kmer, basepair):
            return self._basepairs[(kmer, basepair)]
        return None

    def get_threaded_nucleotide_fixed(self, kmer, basepair):
        """Return fixed threaded nucleotide for one k-mer/basepair coordinate.
        
        Args:
            kmer (Any): Value used by this routine.
            basepair (Any): Value used by this routine."""
        if self._basepairs_fixed is None: self._check_parsing()
        if self.has_basepair_fixed(kmer, basepair):
            return self._basepairs_fixed[(kmer, basepair)]
        return None
   
    def has_kmer(self, kmer):
        """Return whether a k-mer exists in threaded DNA map.
        
        Args:
            kmer (Any): Value used by this routine."""
        if self._dna is None or self._dna=={}: self._check_parsing()
        return kmer in iter(self._dna.keys())
  
    def has_kmer_fixed(self, kmer):
        """Return whether a k-mer exists in fixed DNA map.
        
        Args:
            kmer (Any): Value used by this routine."""
        if self._dna_fixed is None or self._dna_fixed=={}: self._check_parsing()
        return kmer in iter(self._dna_fixed.keys())


    def get_kmer(self, kmer):
        """Return binding start for a threaded k-mer.
        
        Args:
            kmer (Any): Value used by this routine."""
        if self._dna is None or self._dna=={}: self._check_parsing()
        if self.has_kmer(kmer):
            return self._dna[kmer]
        return None

    def get_kmer_fixed(self, kmer):
        """Return binding start for a fixed threaded k-mer.
        
        Args:
            kmer (Any): Value used by this routine."""
        if self._dna_fixed is None or self._dna_fixed=={}: self._check_parsing()
        if self.has_kmer(kmer):
            return self._dna_fixed[kmer]
        return None


    #def get_kmer_escore(self, kmer):
    #    if self.has_kmer(kmer):
    #        return self._dna[kmer][-1]
    #    return None

    def get_kmer_interface(self, kmer):
        """Return interval string `start-end` covered by a threaded k-mer.
        
        Args:
            kmer (Any): Value used by this routine."""
        if self._dna is None or self._dna=={}: self._check_parsing()
        if self.has_kmer(kmer):
            return self._dna[kmer] + "-" + str(int(self._dna[kmer]) + len(kmer) - 1)
        return None

    def get_kmer_fixed_interface(self, kmer):
        """Return interval string `start-end` for a fixed threaded k-mer.
        
        Args:
            kmer (Any): Value used by this routine."""
        if self._dna_fixed is None or self._dna_fixed=={}: self._check_parsing()
        if self.has_kmer_fixed(kmer):
            return self._dna_fixed[kmer] + "-" + str(int(self._dna_fixed[kmer]) + len(kmer) - 1)
        return None


    def get_alignment(self, alignment="pdb"):
        """Return stored template or query alignment string.
        
        Args:
            alignment (Any): Sequence/alignment content used in this step."""
        if alignment != "pdb" and alignment != "query":
            alignment = "pdb"
        if alignment == "pdb":
            if self._pdb_ali is None: self._check_parsing()
            return self._pdb_ali
        if alignment == "query":
            if self._query_ali is None: self._check_parsing()
            return self._query_ali
        return None

    def write(self, file_name):
        """
                Write threading data to canonical thread file format.
        
                Args:
                    file_name (str): Output threading file path.
        
                Returns:
                    None.
        
        Args:
            file_name (Any): Value used by this routine."""
        if self._pdb_name is None or self._protein is None or self._dna is None: self._check_parsing()
        functions.write(file_name,"# PDB name    : %s" % self._pdb_name)
        functions.write(file_name,"# PDB chain   : %s" % self._pdb_chain)
        if self._identity != None: functions.write(file_name,"# %% identity  : %s" % str(self._identity) )
        if self._coverage != None: functions.write(file_name,"# %% coverage  : %s" % str(self._coverage) )
        functions.write(file_name,"# DNA helix   : %s" % str(self._helix))
        functions.write(file_name,"# Query align : %s" % self._query_ali)
        functions.write(file_name,"# PDB align   : %s" % self._pdb_ali)
        functions.write(file_name,">protein")
        #for (name,number),aminoacid in self._protein.iteritems():
        for key in sorted(self._protein, key=lambda x: x[-1]):
            name=key[0]
            number=key[1]
            aminoacid=self._protein.get(key)
            functions.write(file_name,"%s;%s;%s" % (name,number,aminoacid))
        functions.write(file_name,"//")
        functions.write(file_name,">dna")
        for dna_sequence,binding in self._dna.items():
            functions.write(file_name,"%s;%s" % (dna_sequence,binding))
        functions.write(file_name,"//")
        functions.write(file_name,">dna_fixed")
        for dna_sequence,binding in self._dna_fixed.items():
            functions.write(file_name,"%s;%s" % (dna_sequence,binding))
        functions.write(file_name,"//")

  

def parse_options():
    """
    Parse command-line options for homolog-to-template threading generation.

    How to run:
        python threader.py -i HOMOLOG_FILE --pdb TEMPLATE_PDB --chain CHAIN_ID --dna DNA_SEQ
            [--filter --specie SPECIE --code CODE --dummy DUMMY_DIR -o OUTPUT_DIR -v]

    Example:
        python threader.py -i homologs.out --pdb model.pdb --chain A --dna ACGTACGT --specie 9606 -v

    The parser configures:
        - Homolog input source (`-i`) and optional ortholog filtering
          (`--filter`, `--specie`, `--code`).
        - Template context (`--pdb`, `--chain`) and threaded DNA sequence (`--dna`).
        - Runtime/output locations (`--dummy`, `--out`) and verbosity (`-v`).

    Args:
        None.

    Returns:
        optparse.Values: Parsed CLI options for input resources and output.
    """

    parser = optparse.OptionParser(
        "python threader.py -i HOMOLOG_FILE --pdb TEMPLATE_PDB --chain CHAIN_ID --dna DNA_SEQ "
        "[--filter --specie SPECIE --code CODE --dummy DUMMY_DIR -o OUTPUT_DIR -v]"
    )

    parser.add_option("-i", action="store", type="string", dest="input_file", default=None, help="Input blast/hmm file", metavar="INPUT_FILE")
    parser.add_option("--filter", default=False, action="store_true", dest="filter_hits", help="Filter twilight zone hits and  homologs that do not cover 100% of interface (default = False)", metavar="{boolean}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_file", help="Input PDB file (whole complex with DNA)", metavar="PDB_FILE")
    parser.add_option("--chain", action="store", default=None, type="string", dest="chain", help="Selected chain of PDB", metavar="CHAIN")
    parser.add_option("--dna", action="store", default=None, type="string", dest="dna_seq", help="DNA sequence to be thread", metavar="DNA_SEQ")
    parser.add_option("--dummy", action="store", type="string", dest="dummy_dir", default="/tmp", help="Dummy directory (default is /tmp)", metavar="DUMMY_DIR")
    parser.add_option("-o","--out", action="store", default="aux_files", type="string", dest="output_dir", help="Output directory (default is 'aux_files')", metavar="OUTPUT_FILE")
    parser.add_option("--specie", action="store", type="string", dest="specie",  default=None, help="Specie to obtain specific orthologs (i.e. taxon/code/common_name as 9606/HUMAN/'Homo sapiens')", metavar="SPECIE")
    parser.add_option("--code", action="store", type="string", dest="code",  default=None, help="Use only a single sequence with the specific code", metavar="CODE")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("-d", default="basepairs", action="store", type="string", dest="distance_type", help="Distance type (i.e. \"basepairs\", \"dinucleotides\" or \"mindist\"; default = dinucleotides)", metavar="{string}")

    (options, args) = parser.parse_args()

    if options.input_file is None  or options.pdb_file is None  or options.chain is None or options.dna_seq is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options



def main():
    """
    Run the homolog-to-template threading workflow.

    Workflow:
        1. Parse command-line options and prepare output/dummy directories.
        2. Load the template protein-DNA structure and derive DNA/interface context.
        3. Parse homolog/ortholog hits and map aligned residues to template indices.
        4. Build and write per-hit threading files (protein and DNA assignments).
        5. Remove temporary files.

    Args:
        None.

    Returns:
        None. Threading files are written to the output directory.
    """

    # Arguments & Options #
    options   = parse_options()

    # Initialize
    output_dir=options.output_dir
    dummy_dir=options.dummy_dir
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    if not os.path.exists(dummy_dir):  os.mkdir(dummy_dir)
    pdb_chain = options.chain
    code      = options.code


    # Get PDB object #
    dna_sequence={}
    pdb_all   = PDB(os.path.abspath(options.pdb_file))
    pdb_name  = pdb_all.id.strip()
    chain_obj = pdb_all.get_chain_by_id(pdb_chain)
    pdb_obj   = PDB()
    pdb_obj.add_chain(chain_obj)
    for chain_id in pdb_all.chain_identifiers:
        chain=pdb_all.get_chain_by_id(chain_id)
        if chain.chaintype=="N": 
           pdb_obj.add_chain(chain)
           dna_sequence.setdefault(chain_id,chain.nucleotide_sequence())
    dummy_pdb_file=pdb_name+"_"+pdb_chain+".pdb"
    if options.verbose:sys.stdout.write("\t-- Read PDB file with DNA %s \n"%dummy_pdb_file)
    dummy_file=os.path.join(os.path.abspath(options.dummy_dir),dummy_pdb_file )
    if os.path.exists(dummy_file):
       os.remove(dummy_file)
    pdb_obj.write(dummy_file)

    # Get X3DNA object #
    if options.verbose:sys.stdout.write("\t-- Get DNA\n")
    x3dna_obj = x3dna.get_x3dna_obj(dummy_file, os.path.abspath(options.dummy_dir))

    # Get contacts object #
    if options.verbose:sys.stdout.write("\t-- Get contacts\n")
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj, "pdi", options.distance_type, os.path.abspath(options.dummy_dir))

    # Get helix #
    if options.verbose:sys.stdout.write("\t-- Get DNA helices\n")
    dna_helices={}
    for helix in x3dna_obj.get_dna_helices():
        # Get helix dinucleotides #
        dinucleotides = x3dna_obj.get_helix_dinucleotides(helix)
        # For each contact... #
        for contact_obj in contacts_obj.get_contacts():
            if contact_obj._A_chain == pdb_chain:
               dna_helices.setdefault(helix,set()).add(tuple(contact_obj._B_chain))

    # Get interface object #
    if options.verbose:sys.stdout.write("\t-- Get interface\n")
    interface_obj = interface.get_interface_obj(pdb_obj, x3dna_obj, contacts_obj, os.path.abspath(options.dummy_dir))

    # Get sequence/PDB correlation #
    sequence_to_crystal, crystal_to_sequence = model_protein.get_sequence_to_crystal_correlations(pdb_obj, pdb_chain)
    # get homologs #
    homologs_obj  = HOMO.Homologs(os.path.abspath(options.input_file))
    ortholog_hits = homologs_obj.get_orthologs(options.specie)

    # Run 
    for dna_helix in dna_helices.keys():
        #print dna_sequence[list(dna_helices[dna_helix])[0][0]]
        dna_sequence_interface  = dna_sequence[list(dna_helices[dna_helix])[0][0]][interface_obj.get_start():interface_obj.get_interface_length()+interface_obj.get_start()]
        #print dna_sequence_interface
        for ortholog in ortholog_hits.get_hits((0,None)):
          # Add homolog #
           mm = re.search("^(tr|sp)\|(\S+)\|\S+\_(\S+)", ortholog.sequenceID)
           if mm:
              ortholog_name = mm.group(2)
           else:
              if ortholog.sequenceID.startswith("tr") or ortholog.sequenceID.startswith("sp"):
                 ortholog_name=ortholog.sequenceID.strip().split()[0].split("|")[1]
              else:
                 ortholog_name=ortholog.sequenceID.strip().split()[0].split("|")[0]
           if code is not None:
              if ortholog_name != code: 
                 continue
              else:
                 if options.verbose:sys.stdout.write("\t-- Found %s \n"%(ortholog_name))
           threading_file = os.path.join(os.path.abspath(output_dir),ortholog_name + "_" + dna_helix + ".txt")
           threaded_ortholog = Threaded()
           if options.verbose and os.path.exists(threading_file):sys.stdout.write("\t-- Use threading file %s\n"%threading_file)
           if not os.path.exists(threading_file):
                    # Initialize #
                    try:
                     template_alignment = ortholog.sequences[0].sequence.upper()
                     ortholog_alignment = ortholog.sequences[1].sequence.upper()
                     coverage = 100 * float (ortholog.get_min_coverage_of_full_sequence_segment())
                     identity = ortholog.identities_pec
                    except:
                     if options.verbose: sys.stdout.write("\t\t--Skip %s\n"%ortholog.sequenceID)
                     continue
                    # Get the correct aminoacid number
                    first = pdb_obj.get_chain_by_id(pdb_chain).gapped_protein_sequence.find(template_alignment.replace("-", ""))
                    if first<0:first=0
                    positions = list(range(first, len(template_alignment.replace("-", "")) + first))
                    # Thread ortholog sequence
                    protein = {}
                    for position in range(len(template_alignment)):
                      # Skip if hit alignment position is gapped #
                      if template_alignment[position] != "-":
                        sequence_position = positions.pop(0)
                        if sequence_position in sequence_to_crystal:
                           aminoacid_number = sequence_to_crystal[sequence_position]
                           protein.setdefault((pdb_chain,aminoacid_number), ortholog_alignment[position])
                    # Thread DNA
                    dna = {}
                    dnae = {}
                    if len(options.dna_seq) != interface_obj.get_interface_length():
                       if  len(options.dna_seq) < interface_obj.get_interface_length():
                           dna_seq     = options.dna_seq + "N"*(interface_obj.get_interface_length()-len(options.dna_seq) )
                           dna_seq_new = options.dna_seq + dna_sequence_interface[len(options.dna_seq):]
                           if options.verbose: sys.stdout.write("Check %s DNA sequence is too short: %s -> %s\n"%(ortholog_name + "_" + dna_helix + ".txt",dna_seq,dna_seq_new))
                       if  len(options.dna_seq) > interface_obj.get_interface_length():
                           dna_seq     = options.dna_seq
                           dna_seq_new = options.dna_seq[0:interface_obj.get_interface_length()]
                    else:
                           dna_seq     = options.dna_seq
                           dna_seq_new = options.dna_seq
                    dna.setdefault(dna_seq.upper(),interface_obj.get_start())
                    dnae.setdefault(dna_seq_new.upper(),interface_obj.get_start())
                    # Create threaded object #
                    threaded_ortholog.set_pdb_name(pdb_name)
                    threaded_ortholog.set_pdb_chain(pdb_chain)
                    threaded_ortholog.set_dna_helix(dna_helix)
                    threaded_ortholog.set_identity(identity)
                    threaded_ortholog.set_coverage(coverage)
                    threaded_ortholog.set_protein(protein)
                    threaded_ortholog.set_dna(dna)
                    threaded_ortholog.set_dna_fixed(dnae)
                    threaded_ortholog.set_query_ali(ortholog_alignment)
                    threaded_ortholog.set_pdb_ali(template_alignment)

                    # Create threaded file #
                    if options.verbose:sys.stdout.write("\t-- Write threading file %s\n"%threading_file)
                    threaded_ortholog.write(threading_file)
                
    # Clean files #
    os.remove(dummy_file)


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
    
