#!/usr/bin/python

"""
Heuristic TF selector based on PDB, UniProt, and GO-term evidence.

This builder script classifies candidate PDB chains as TF/non-TF, extending an
existing TF list with newly inferred TF entries and a retry list for failures.
"""

#===================================================
#DESCRIPTION
#
#
#===================================================

#---------  Import     ----------------------------#

import urllib.request, urllib.parse, urllib.error
import os, sys, re
from optparse import OptionParser
import configparser

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


#---------  Functions  ----------------------------#

#===================================================
#
# General utility functions
#
#===================================================

def readFile( fileName ):

    f = open( fileName )
    lines = f.readlines( )

    for i in range( 0, len( lines ) ):

        lines[i] = lines[i].rstrip( )

    return lines

def readURL( url ):

    k = urllib.request.urlopen( url )
    lines = k.read( ).decode('utf-8').split("\n")
    k.close( )
    return lines

def getPDBEntry( pdbid ):

    k = readURL( "https://www.rcsb.org/pdb/files/"+pdbid+".pdb" )
    return k

def getChains( pdbEntry ):

    chains = []

    i = 0
    while( i < len( pdbEntry ) and pdbEntry[i][:6] != "SEQRES" ):
        i += 1
    while( i < len( pdbEntry ) and pdbEntry[i][:6] == "SEQRES" ):
        if( pdbEntry[i][7:11] == "  1 " ):
            chain = pdbEntry[i][11:12]
            if( pdbEntry[i][19:20] != " " ):#DNA residues are a space
                chains.append( chain )
        i += 1

    return chains

def getGOTermsFile( pdbid, chain, pdbFile ):

    uniprotId = getMappingFromPDB( chain, pdbFile )

    if( uniprotId == None ):

        uniprotId = getMappingFromPDBSWS( pdbid, chain )

    if( uniprotId != None ):

        k = getGOTermsFromUniProt( uniprotId )

        return k

    else:

        return []

def getOnlyGOTerms( fileLines ):
    terms = []
    for line in fileLines:
        k = line.find( "GO:" )
        if( k != -1 ):
            j = k
            done = False
            while( j < len( line ) and line[j:j+1] != '"' ):
                j += 1
            terms.append( line[k:j] )
    return terms

def getMappingFromPDB( chain, pdbFile ):

    uniprotId = None

    i = 0
    found = False
    while( i < len( pdbFile ) and not found ):

        line = pdbFile[i]

        if( line[0:5] == "DBREF" ):
            if( line[26:33].strip( ) == "UNP" ):
                if( line[12:13] == chain.upper( ) ):
                    uniprotId = line[33:42].strip()
                    found = True
        i += 1

    return uniprotId

def getMappingFromPDBSWS( pdbid, chain ):

    lines = readURL( "http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?qtype=pdb&id="+
                 pdbid+"&chain="+chain )

    uniprotId = None

    i = 0
    found = False

    while( i < len( lines ) and not found ):

        if( len( lines[i] ) >= 4 and lines[i][0:4] == "<tr>" ):
            myLine = lines[ i + 1 ]

            k = 0 # '>' counter
            j = 0 # Substring index
            while( j < len( myLine ) and k < 8 ):
                if( myLine[j:j+1] == ">" ):
                    k += 1
                j += 1

                m = j #other substring index
                while( m < len( myLine ) and myLine[m:m+1] != "<" ):
                    m += 1
                if( myLine[j:m] != "" ):#then PDBSWS has an id
                    uniprotId = myLine[j:m]

            found = True
        i += 1

    return uniprotId


def getGOTermsFromRest( pdbid, chain ):
    import json

    returnValue = None

    k = readURL( "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/"+
                 pdbid.upper( )+"/"+chain.upper( ) )
    annotations=[]
    if( len( k ) > 0 ):
        x=json.loads(".".join(k))
        if "rcsb_polymer_instance_annotation" in x:
           for z in x["rcsb_polymer_instance_annotation"]:
               if "name" in z:
                   annotations.append(z["name"])
               
        #if( k[1] != "<goTerms />" ):
        if len(annotations)>0:
            #returnValue = k
            returnValue = annotations

    return returnValue

def getGOTermsFromUniProt( uniprotId ):

    returnVal = readURL( "http://www.uniprot.org/uniprot/"+uniprotId+".xml" )
    if( len( returnVal ) > 1 ):
        return returnVal#OTHER HANDLING OF ERRONEOUS MAPPING
    else:
        return []

def getStrongNegativeGOTerms( ):

    files_path = config.get("Paths", "files_path")
    NonTF_molecular_function = os.path.join(files_path, config.get("Paths", "nTF_GOMF"))
    k = readFile( NonTF_molecular_function )
    m = []
    for line in k:
        m.append( line.split(" ")[0] )

    return m

def getStrongPositiveGOTerms( ):

    files_path = config.get("Paths", "files_path")
    TF_molecular_function_w = os.path.join(files_path, config.get("Paths", "TF_GOMF"))
    k = readFile( TF_molecular_function_w )
    strongPos = []

    for line in k:

        p = line.split(" ")
        if len(p)<2: continue
        #if ( p[2].split("\t")[0] == "T" ):
        if ( p[2][0:1] == "T" ):
            strongPos.append( p[0] )

    return strongPos

def getWeakPositiveGOTerms( ):

    files_path = config.get("Paths", "files_path")
    TF_molecular_function_w = os.path.join(files_path, config.get("Paths", "TF_GOMF"))
    k = readFile( TF_molecular_function_w )
    weakPos = []

    for line in k:

        p = line.split(" ")
        if len(p)<2: continue

        #if( p[2].split("\t")[0] == "F" ):
        if( p[2][0:1] == "F" ):
            weakPos.append( p[0] )

    return weakPos
def getPositiveBiologicalProcessGOTerms( ):

    files_path = config.get("Paths", "files_path")
    TF_biological_process_w = os.path.join(files_path, config.get("Paths", "TF_GOBP"))
    k = readFile( TF_biological_process_w )

    goTerms = []
    for line in k:
        goTerms.append( line.split(" ")[0] )

    return goTerms
def getNegativeBiologicalProcessGOTerms( ):

    files_path = config.get("Paths", "files_path")
    NonTF_biological_process = os.path.join(files_path, config.get("Paths", "nTF_GOBP"))
    k = readFile( NonTF_biological_process )
    goTerms = []

    for line in k:
        goTerms.append( line.split(" ")[0] )

    return goTerms

def hasIntersection( listA, listB):

    exists = False

    i = 0
    while( i < len( listA ) and not exists ):
        j = 0
        while( j < len( listB ) and not exists ):
            if( listA[i] == listB[j] ):
                exists = True
            j += 1
        i += 1

    return exists


def isTFFromPDBKeywords( pdb ):

    files_path = config.get("Paths", "files_path")
    try:
     negative_keywords = os.path.join(files_path, config.get("Paths", "negKW"))
     positive_keywords = os.path.join(files_path, config.get("Paths", "posKW"))
     negative_keywds = readFile( negative_keywords )
     positive_keywds = readFile( positive_keywords )
    except Exception as e:
     print("Failed isTFFromPDBKeywords, not found files")
     print(e)
     raise(e)
    negative_keywds_dic = {}
    positive_keywds_dic = {}

    tFound = False

    for keywd in negative_keywds:
        negative_keywds_dic[keywd.upper( )] = False

    for keywd in positive_keywds:
        positive_keywds_dic[keywd.upper( )] = False

    select_lines = []#Not all lines are searched

    for line in pdb:
        if( ( len( line ) >= 5 and line[0:5] == "TITLE" ) or ( len( line ) >= 6 and line[0:6] == "KEYWDS" ) ):
            select_lines.append( line )

    i = 0
    done = False

    #CHECK PDB
    while( i < len( select_lines ) and not done ):
        line = select_lines[i].upper( )
        print("CHECK PDB ",line)
        for keywd in negative_keywds_dic:
            #print("CHECK NEGATIVE KEYWORD ",keywd)
            if( line.upper().find( keywd.upper() ) != -1 ):
                print("Found negative KEY",keywd.upper())
                negative_keywds_dic[keywd] = True
                done = True
        for keywd in positive_keywds_dic:
            #print("CHECK POSITIVE KEYWORD ",keywd)
            if( line.upper().find( keywd.upper() ) != -1 ):
                print("Found positive KEY",keywd.upper())
                positive_keywds_dic[keywd] = True
        i += 1

        if ( line.upper().find( "TRANSCRIPTION" ) != -1 ) :
           print("Found positive KEY TRANSCRIPTION")
           tFound = True

    return keywordDecision( positive_keywds_dic, negative_keywds_dic, tFound )

def isTFFromUNPKeywords( uniprotLines ):
    #CHECK UNIPROT
    try:
     files_path = config.get("Paths", "files_path")
     negative_keywords = os.path.join(files_path, config.get("Paths", "negKW"))
     positive_keywords = os.path.join(files_path, config.get("Paths", "posKW"))
     negative_keywds = readFile( negative_keywords )
     positive_keywds = readFile( positive_keywords )
    except Exception as e:
     print("Failed isTFFromUNPKeywords, not found files")
     print(e)
     raise(e)

    negative_keywds_dic = {}
    positive_keywds_dic = {}

    tFound = False

    for keywd in negative_keywds:
        negative_keywds_dic[keywd.upper( )] = False

    for keywd in positive_keywds:
        positive_keywds_dic[keywd.upper( )] = False

    keyword_lines = []

    for line in uniprotLines:
        if( line.find( "<keyword" ) != -1 ):
            keyword_lines.append( line.upper( ) )

    i = 0
    done = False

    while( i < len( keyword_lines ) and not done ):

        line = keyword_lines[i]
        print("CHECK UNIPROT ",line)
        for keywd in negative_keywds_dic:
            if( line.upper().find( keywd.upper() ) != -1 ):
                negative_keywds_dic[keywd] = True
                print("Found negative KEY",keywd.upper())
                done = True
        for keywd in positive_keywds_dic:
            if( line.upper().find( keywd.upper() ) != -1 ):
                positive_keywds_dic[keywd] = True
                print("Found positive KEY",keywd.upper())

        if( line.find( "TRANSCRIPTION" ) != -1 ):
                print("Found positive KEY TRANSCRIPTION")
                tFound = True

        i += 1

    return keywordDecision( positive_keywds_dic, negative_keywds_dic, tFound )

def keywordDecision( positive_keywds_dic, negative_keywds_dic, tFound ):

    isTF = None

    for keywd in negative_keywds_dic:
        if( negative_keywds_dic[keywd] ):
            print("NEGATIVE DECISION ON ",keywd.upper())
            isTF = False

    if( isTF == None ):
        """
        If no negative keywords were found, then
        we will consider a positive keyword to
        indicate this protein is likely a transcription factor
        """
        for keywd in positive_keywds_dic:
            if( positive_keywds_dic[keywd] and tFound ):
                print("POSITIVE DECISION ON ",keywd.upper())
                isTF = True

    return isTF
#===================================================
#
#  Decision functions
#
#  In each case returns: False  !TF
#                        True    TF
#
#  A. Molecular function GO Terms
#     1) Strong negative indicator
#     2) Strong positive indicator
#     3) Weak positive
#  B. Biological process
#     1) Positive or negative
#  C. Keywords
#     1) From PDB
#     2) From Uniprot
#
#===================================================

def isTF( pdbid, chain, pdbLines ):

  try:
    #PDB obsolete information
    #goTermsFile = getGOTermsFromRest( pdbid, chain )
    #Rest is only GO terms, not the entire UniProt entry
    goTermsFile = None

    if( goTermsFile != None ):
        hasNoKeywords = True
    else:
        goTermsFile = getGOTermsFile( pdbid, chain, pdbLines )
        hasNoKeywords = False

    onlyGOTerms = getOnlyGOTerms( goTermsFile )

    print("GO TERMS",onlyGOTerms)
    try:
     molNeg = getStrongNegativeGOTerms( )
     molSPos = getStrongPositiveGOTerms( )
     molWPos = getWeakPositiveGOTerms( )
     posBioProc = getPositiveBiologicalProcessGOTerms( )
     negBioProc = getNegativeBiologicalProcessGOTerms( )
    except Exception as e:
     print("Failed isTF, not found files")
     print(e)
     raise(e)

    print("POSITIVE TERMS",molSPos,molWPos,posBioProc)
    print("NEGATIVE TERMS",molNeg,negBioProc)
    isTF = False
    if ( not hasIntersection( molNeg, onlyGOTerms ) ):
        if( hasIntersection( molSPos, onlyGOTerms ) ):
            print("FOUND POSITIVE FUNCTIONS",set(molSPos).intersection(set(onlyGOTerms)))
            isTF = True
        else:
            if( hasIntersection( molWPos, onlyGOTerms ) ):
                print("FOUND WEAK POSITIVE FUNCTIONS",set(molWPos).intersection(set(onlyGOTerms)))
                if( hasIntersection( posBioProc, onlyGOTerms ) ):
                    print("FOUND POSITIVE PROCESSES",set(posBioProc).intersection(set(onlyGOTerms)))
                    isTF = True
                elif( hasIntersection( negBioProc, onlyGOTerms ) ):
                    print("FOUND NEGATIVE PROCESSES",set(negBioProc).intersection(set(onlyGOTerms)))
                    isTF = False
                else:
                    #WAS OUT ONE TAB IN ( IF !TF )
                    #if( hasNoKeywords ):
                      isTF1 = isTFFromPDBKeywords( pdbLines )
                    #else:
                      isTF2 = isTFFromUNPKeywords( goTermsFile )
                      #True -> check UniProt keywords
                      if isTF1 or isTF2: isTF = True 
            elif( onlyGOTerms != [] ):
                print("NONE OF THE GO TERMS IS ACCEPTED")
                isTF = False
            else:#No UniProt information available
                isTF = isTFFromPDBKeywords( pdbLines )
                #False -> check PDB
    else:
        print("FOUND NEGATIVE FUNCTIONS",set(molNeg).intersection(set(onlyGOTerms)))

    print("isTF ?",isTF)

    if( isTF == None ):
        return ( "Non-TF ( manually verify )" )
    elif( isTF ):
        return ( "TF" )
    else:
        return ( "Non-TF" )
  except Exception as e:
     print("Failed isTF ",e)
     raise(e)


#=========================================================
#
# Read options to run as standalone
# 
#=========================================================

def parse_options():
 """
 Parse CLI options for standalone TF selection runs.

 How to run:
  ``python TFinderSelect.py -i pdb_list.csv -t known_tfs.txt -o tf_out.txt``

 Returns:
  optparse.Values: Namespace with input/output/retry file paths.
 """
 parser = OptionParser( )
 parser.add_option( "-i", "--file", dest="fileName", default=None, help="Supply a file listing PDB IDs w/ or w/o chain IDs" )
 parser.add_option( "-t", "--TFs", dest="fileTF", default=None, help="Supply a file listing PDB IDs w/ or w/o chain IDs of known TFs" )
 parser.add_option( "-r", "--Retry", dest="fileRetry", default=None, help="Supply a file listing PDB IDs w/ or w/o chain IDs of PDBs that failed to be tested" )
 parser.add_option( "-o", "--output_tfindit", dest="fileOut", default=None, help="Output file with format 'CODE,CHAIN' of Transcription Factors" )

 ( options, args ) = parser.parse_args( )

 
 return options

#=========================================================
#
#  Main program to select TF
#  pdbs is a list of tuples: pdb ID, Chain for all pdb entries
#  tfs is the same format list of controlled TFs
#  result returned is a tuple of 2 sets: new list of Tfs not described before in TF list
#  plus a list of the same format with internet connection failures to retry
#
#=========================================================

def fileExist (file):               #Checks if a file exists AND is a file
 if file is not None: return os.path.exists(file) and os.path.isfile(file)
 else: return False



def TFinder(pdbs,tfs):
  """
  Evaluate candidate PDB chains and collect newly inferred TFs.

  Args:
   pdbs (set): ``(pdb_id, chain_or_none)`` candidates.
   tfs (set): Existing ``(pdb_id, chain_or_none, family)`` TF entries.

  Returns:
   tuple: ``(tf_new, retry)`` sets with newly classified TFs and failed cases.
  """
  tf_new=set()
  retry=set()
  tfs_chain={}
  for pdb,chain,family_set in tfs:
   print("PDB",pdb,chain, " in families ",family_set)
   tfs_chain.setdefault(pdb.lower(),set()).add(chain)
  for pdbid,chain in pdbs:
    print("Check %s %s"%(pdbid.lower(),chain))
    if pdbid.lower() in tfs_chain: 
      print("Found in TFs")
      continue
    try:
          print("get PDB",pdbid)
          entry = getPDBEntry( pdbid )
          if( chain == None ):
            chains = getChains( entry )
            print("Chains",pdbid,chains)
            for chain in chains:
               try:
                if (isTF( pdbid.lower(), chain, entry )=="TF"):
                   print("Add TF ",pdbid.lower(), chain)
                   tf_new.add( (pdbid.lower(),chain.upper(),"Unknown") )
                else:
                   print("It is not TF",pdbid.lower(), chain)
               except Exception as e:
                   print("Retry later, failed ",e)
                   retry.add( (pdbid.lower(),chain.upper()) )
                   continue
          else:
            try:
              print("Test as TF ",pdbid.lower(),chain)
              if (isTF( pdbid.lower(), chain , entry )=="TF"):
                  print("Add TF ",pdbid.lower(), chain)
                  tf_new.add( (pdbid.lower(),chain.upper(),"Unknown") )
              else:
                  print("It is not TF",pdbid.lower(), chain)
            except Exception as e:
              print("Retry later, failed ",e)
              retry.add( (pdbid.lower(),chain.upper()) )
              continue
    except Exception as e:
        print("Failed connection, error ",e)
        if( chain == None ):retry.add( (pdbid.lower(),None))
        else:retry.add( (pdbid.lower(),chain.upper()) )
        continue
    
  return tuple([tf_new,retry])


def main():

    # Step 1) Parse CLI options and load the candidate PDB list.
    options = parse_options()
 
    pdbs=set()
    if not fileExist(options.fileName):
     print("Input PDB-list is missing\n")
     sys.exit(10)
    f = open( options.fileName ,'r')
    lines = f.readlines()
    for line in lines:
       data=line.strip().split(",")
       if len(data)>1:
          pdbs.add((data[0],data[1]))
       else:
          pdbs.add((data[0],None))
    f.close()

    # Step 2) Load the existing TF reference list (if provided).
    tfs=set()
    if fileExist(options.fileTF):
     f = open( options.fileTF ,'r')
     lines = f.readlines()
     for line in lines:
       data=line.strip().split(";")
       if len(data)>2:
          tfs.add((data[0],data[1],data[2]))
       elif len(data)==2:
          tfs.add((data[0],data[1],"Unknown"))
       else:
          tfs.add((data[0],None,"Unknown"))
     f.close()
    else:
     print("Input TF-list is missing\n")

  
    # Step 3) Run TF classification and capture retry-able failures.
    retry=set()
    tf_new=set()
    (tf_new,retry)=TFinder(pdbs,tfs)

    # Step 4) Write retry queue and new TF predictions.
    if (options.fileRetry is not None):
     for pdb_name,pdb_chain in retry:
       functions.write(options.fileRetry, "%s;%s" % (pdb_name, pdb_chain))  
  
    for pdb_name,pdb_chain,pdb_family in tf_new:
       functions.write(options.fileOut, "%s;%s;%s" % (pdb_name, pdb_chain,pdb_family))

    

if  __name__ == "__main__":
    main()


