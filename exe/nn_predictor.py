"""
Predict nearby PWM candidates for a TF sequence by BLAST nearest neighbors.

This script searches an input TF sequence against a database of TF sequences,
keeps hits above a user-defined identity threshold, and copies the associated
PWM files into an output directory together with a weight table.
"""

import os,sys
import time
import optparse
import configparser
import shutil

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
from ModCRElib.sequence import blast


#-----------#
#  Options  #
#-----------#

def options():
    """
    Parse the command line options for nearest-neighbor PWM prediction.

    Returns:
        optparse.Values: Parsed CLI options describing the query sequence,
        sequence database, PWM directory, output directory, and identity threshold.
    """

    parser = optparse.OptionParser('python nn_predictor.py -i input_file -d database_path -m pwm_path [-o output_path -p idp_threshold --dummy dummy_dir -v]')
    
    parser.add_option('-i', action = 'store', type = 'string', dest='input_file', help = 'TF sequence file in fasta format.', metavar='{filename}')
    parser.add_option('-d', action = 'store', type = 'string', dest='database_file', help = 'Complete path to the database of TF sequences.', metavar='{filename}')
    parser.add_option('-m', action = 'store', dest = 'pwm_path', help = 'Complete path to the directory with the PWMs of the selected database', metavar='{directory}')
    parser.add_option('-o', action = 'store', default = './output/', type = 'string', dest='output_directory', help = 'Path to the output directory (default = ./output). \n If it does not exist, will be created.', metavar = '{directory}')
    parser.add_option('-p', default = 0, action = 'store', dest = 'idp', help = 'Define the identity percentage threshold (default = 0)')
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode. If not selected the dummy directory will be removed (default = False)", metavar="{boolean}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
 
    (options, args) = parser.parse_args()

    if options.input_file is None or options.database_file is None or  options.pwm_path is None:
        parser.error('missing arguments: type option \"-h"\ for help')
    
    return options


def main():
    """
    Run the command-line nearest-neighbor PWM prediction workflow.

    Workflow:
        1. BLAST the input TF sequence against the provided TF database.
        2. Filter hits by percent identity.
        3. Copy the PWM files associated with the accepted hits.
        4. Write ``PWM_ID_weigths.txt`` with copied PWM paths and weights.

    Returns:
        None. Output files are written into the selected output directory.
    """
    start = time.time()
    
    options = options()
    inp = options.input_file
    db = options.database_file
    idp_threshold = float(options.idp)
    pwm_path = options.pwm_path 
    verbose = options.verbose
    dummy_dir = options.dummy_dir
    out = os.path.abspath(options.output_directory)
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    if not os.path.exists(out): os.makedirs(out)

    blast_object = blast.get_blast_obj(db, inp, dummy_dir)
   
    result={} 
    for hit in blast_object.get_hits():
        ID = hit.sequenceID
        identity = hit.identities
        length = hit.align_length
        idp = identity/length * 100
        if idp >= idp_threshold:
            meme_file = os.path.join(pwm_path,ID + '.MEME')
            if not os.path.exists(meme_file): meme_file = os.path.join(pwm_path,ID + '.meme')
            if not os.path.exists(meme_file): meme_file = os.path.join(pwm_path,ID + '.memes')
            if not os.path.exists(meme_file):
               if verbose: print("Skip file not found %s"%meme_file)
            else:
               outfile = os.path.join(out,os.path.basename(meme_file))
               shutil.copy(meme_file,outfile)
               result.setdefault(outfile,idp)
    end = time.time()
    output = open(os.path.join( out, 'PWM_ID_weigths.txt') , 'w')
    for outfile,idp in sorted(result.items(),key=lambda pair: pair[1],reverse=True):
               output.write("%s\t%10.5f\n"%( outfile,idp))
    output.close()
    if not verbose: shutil.rmtree(dummy_dir)
    if verbose: print('Done, execution time:', end - start)


#-----------------------#
#  Blast TF to DNA PWM  #
#-----------------------#

if __name__ == '__main__':
    main()
