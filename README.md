
#     ModCRElib: ModCRE to be deployed as a downloadable Package      

## Installation instructions:
1) download this repository 
2) download the pbm and pdb database folders
    a) wget http://aleph.upf.edu/modcrefiles/pbm.tgz (60 Gb)
    b) wget http://aleph.upf.edu/modcrefiles/pdb.tgz (57 Gb)
3) place these folders into the downloaded reopository
4) decompress them
    a) tar -xvzf pbm.tgz (110 Gb)
    b) tar -xvzf pdb.tgz (88 Gb)
5) Replace the empty pwm folder in pbm with the provided pmw folder in PWMdatabase
6) Ensure all dependencies are installed
7) Change paths in ModCRElib/configure/config.ini to point to the correct location for your installations. 


## Dependencies:

### Python
----
Bio              → from Biopython (`pip install biopython`)<br/>
SBILib           → comes with package but can be installed seperately (`pip install SBILib` see https://github.com/structuralbioinformatics/SBILib for dependencies)<br/>
bs4              → BeautifulSoup (`pip install beautifulsoup4`)<br/>
bottle           → lightweight web framework (`pip install bottle`)<br/>
ihm              → Integrative Modeling library (`pip install ihm`)<br/>
matplotlib       → plotting library (`pip install matplotlib`)<br/>
numpy            → numerical computing (`pip install numpy`)<br/>
pandas           → data analysis (`pip install pandas`)<br/>
plotly           → interactive plotting (`pip install plotly`)<br/>
scipy            → scientific computing (`pip install scipy`)<br/>
seaborn          → statistical plotting (`pip install seaborn`)<br/>
sklearn          → scikit-learn (`pip install scikit-learn`)<br/>

### Enviornmental
---
BLAST+ (run on 2.12.0)<br/>
CD-HIT (run on 4.8.1)<br/>
Clustal-Omega (run on 1.2.4)<br/>
ClustalW2 (run on 2.1)<br/>
EMBOSS (run on 6.6.0)<br/>
Ghostscript (run on 9.53.3)<br/>
HMMER (run on 3.3.2)<br/>
MEME (run on 5.1.1)<br/>
Modeller (run on 10.3)<br/>
Python (run on 3.8.6)<br/>
TMalign (run on Version 30012025)<br/>
dssp (run on 3.0.0)<br/>
x3dna (run on 2.5.0)<br/>

----------
## Using ModCRElib
----------


## Example 1 
Modelling a Transcription Factor (TF) from an amino acid sequence
#### ------
1) We need a single or multi fasta file containing a.a. sequence of protein to be modelled (example of AHR_Example.fa is given)
2) modelling.sh contains the command to model the protein with the following edits necessary:<br/>
    pdb="/home/pgohl/ModCRE_Package/pdb" Needs to be replaced with the location of the downloaded pdb folder on your machine<br/>
    -i uniput protein in fasta (could be multi fasta)<br/>
    -o output directory (models) <br/>
#### ------
We now have a folder containing the modelled Transcription factor in pdb format.<br/>
You can view the file in chimera or online at https://www.rcsb.org/3d-view by uploading the file.<br/>
#### ------

## Example 2
Predict TF binding specificity
#### ------
3) bin/renumberModels.py is a script to ensure that atom and amino acid numbers of models from modelling.sh is continuous (not allways needed)<br/>
    arg 1 = folder containing the pdbs to be renumbered<br/>
    arg 2 = output folder<br/>
    (python bin/renumberModels.py models remodels)<br/>
4) pwm.sh predicts the binding specificity<br/>
    pdb="/home/pgohl/ModCRE_Package/pdb" Needs to be replaced with the location of the downloaded pdb folder on your machine<br/>
    pbm="/home/pgohl/ModCRE_Package/pbm" Needs to be replaced with the location of the downloaded pbm folder on your machine<br/>
#### ------
We now have a predicted binding specificity of the transcription factor in pwm and meme format in the output folder specified in the previous command.<br/> 
Remember that ModCRE is designed to use PWMs in aggregate to predict binding sites and individual PWM predictions may be more or less accurate.<br/> 
#### ------

## Example 3
Scan dna sequence for binding sites
#### ------
5) pwm/make_scan_ready.py is file that generates a database file of the predicted pwms from the previous step that are in the correct format for scanning.<br/>
    -line 4 references "database.txt" (the generated database file from pwm.sh) and the name of the new file to be used for scanning<br/>
6) scan_sequence.sh is the script to run a scan with.
    -i The dna sequence file in fasta format that is to be scanned<br/>
    -o the name of the output folder<br/>
    -s the species identifier (if restricting scan to TFs of a given specie)<br/>
    -ft the fimo threshold value to be used in designating hits<br/>
    --db the location of the database file to be used (the defaul option to be used is the larger pwm folder provided)<br/>
    -c Cluster complexes into connected binding sites<br/>
#### ------
There are several ways to observe the results. The TFs binding the sequence and their binding sites can be retrieved from the file orthologs.json.<br/>
Here we can view the name of the TF experienceing the hit, the start and end index for binding along the DNA sequence and the various orthologs binding.<br/>
Alternatively, scan_sequence.sh (with parameters used above) will generate thread files of the TF and its orthologs binding the DNA sequence in the folder aux_files.<br/>
#### ------

## Example 4
Generate a model of a TF attached to a predicted binding site along a full length of DNA
#### ------
We can view the TFs binding to the full scanned DNA sequence as predicted in the previous step. For this we will need to process the files a little through.
#### ------
7)  modelling.sh contains the option to use thread files as an input instead of the previously used multi fasta files<br/>
    -i a file containng a list of the location of the thread files to be used (Threads_list.txt is provided)<br/>
    -t indicates that threads are used instead of aa sequences<br/>
    -o the output folder location <br/>
8) We copy the desired models (from the output of the previous step) into a folder containg all the binary interactions that we would like to use in the modelled complex<br/>
9) rename_complex_input.py is a python script that prepares scanning output file names for complex builder (eg. python rename_complex_input.py BinaryInteractions/A9YTQ3.5nj8_1A.18-29.pdb ---> BinaryInteractions/A9YTQ3.5nj8_1A.18-29:1:243_TF.pdb)<br/>
    -the name of the file must follow the following format:<br/>
        {UniprotAccession}.{PDBID}_{Chain}.{index of binding start}-{index of binding end}:{model start index}:{model end index}_{a label}.pdb<br/>
10) BuildComplex.sh is the script that will build a complex based on binary interaction files contained in a given folder<br/>
    a) /soft/system/software/x3dna/2.3/bin/fiber will generate a pdb file for a given DNA sequence<br/>
        -seq the dna sequence to be used (in the case of modelling the scanning results use the same sequence)<br/>
        -b the dna conformation to be applied followed by output location <br/>
    b) exe/complexbuilder.py will build the complex<br/>
        -d the folder containing the binary interaction files<br/>
        -o the output folder location <br/>
#### ------
Now the modelled complex can be view in the output folder (Complex/fragment_1-100/dna__1-100_aa.pdb).<br/> 
#### ------

## Example 5
Generate thread files from a modelled TF 
#### ------
11) get_best_bindings_threads.sh produces thread files for use in modelling and retrieving scores.<br/> 
    pdb="/home/pgohl/ModCRE_Package/pdb" Needs to be replaced with the location of the downloaded pdb folder on your machine<br/>
    pbm="/home/pgohl/ModCRE_Package/pbm" Needs to be replaced with the location of the downloaded pbm folder on your machine<br/>
    -i (the pdb of the transcription factor to be used)<br/>
    -o (output directory)  <br/>
    --seq fasta sequence containing DNA sequence to be bound<br/>
    --dna nucleatide sequence of the binding site<br/> 
    --pwm meme file for the transcription factor<br/>
#### ------
The thread file can be used to score the binding of a TF along any DNA sequence that matches the binding site length. <br/>
previous steps need not be repeated for other DNA sequences, simply create a thread file for the relevant substituted sequence by<br/>
replaceing the DNA sequences at the bottom of the threads folder:<br/>
\>dna<br/>
CAGCTGGCTGTG;0<br/>
//<br/>
\>dna_fixed<br/>
CAGCTGGCTGTG;0<br/>
//<br/>
#### ------

## Example 6
Generate a scoring profile of a TF-DNA interaction
#### ------
12) get_best_score.sh produces a scoring profile for the transcription factor on target dna sequence.<br/>
    -i (the input is a thread file as produced by get_best_bindings_threads)<br/>
    -o (output directory)  <br/>
#### ------
Finally we have a score profile (statistical potentials) for the TF binding along the tested DNA binding site.
#### ------

## Example 7
Generate a scoring profile plot for a TF along a DNA sequence
#### ------
13) get_score_profiles.sh will generate the raw scores that will be used to plot the scoring <br/>
    -i a file containing the folders containing the models to be used<br/>
    -d a fasta format dna sequence to profile the TF binding against<br/>
    -o Output name for tables and plots (default is the name of the input FOLDER)<br/>
#### ------
The output will be stored in each folder provided by the input file. A folder will have been generated and named <br/>
profilerinput.txt_profiling.34_272 by default (profilerinput.txt being the input file, profiling the folder from that file).<br/>
Within can be found individual model scores (pickle files) as well as the mean tables (csv files)<br/>
#### ------
14) bin/plotprofile.py will generate a plot from the mean score file<br/>
    arg 1 = path to the mean table to be plotted <br/>
    arg 2 = location the plot is to be stored at<br/>
    arg 3 = column (score type) to be plotted (options will be printed out if provided isn't in file)<br/>
    (python bin/plotprofile.py profiling/profilerinput.txt_profiling.34_272/profilerinput.txt_profiling.Profiletest_1.mean.csv profiling/energy_plot.out_normal_s3dc_dd.png normal_s3dc_dd)<br/>

There is support for running jobs in parallel. In order to do this the relevant information in the config.ini Cluster field must be filled in.<br/>
After that simply run jobs with the parallel parameter<br/>

## Example 8
Aggregate PWM clusters for a TF
#### ------
15) get_json.sh will generate input (json files) for subsequent steps<br/>
    exe/get_json.py takes positional arguments:<br/>
        home path<br/>
        Fasta format file of the TFs to be run<br/>
        file containing TF codes (uniprot) corresponding to TFs in Fasta file<br/>
        Table of family labels for TF codes (provided in files)<br/>
        Table of nearest neighbors (provided in files is one with cases of 30-100% sequence id)<br/>
        folder of the pwms generated with modcre (output of pwm.sh)<br/>
        output folder name <br/>
        uniprot label indicator (uniprot)<br/>
16) get_aggregates.sh<br/>
    input={name of the output from get_json.sh}<br/>
    modcre={folder of the modcre predicted pwms (pwm.sh output)}<br/>
    pvalue={P-value threshold of TOMTOM similarity between two PWMs (default 0.005)} <br/>
    threshold={Distance threshold of the agglomerative clustering (default is 0.01)}<br/>
    length={Length of the binding site of the output PWMs}<br/>
    ModCRElib/msa/aggregate_pwms.py takes some arguments that the user may wish to change<br/>
        --jaspar location of jaspar pwm database<br/>
        --cisbp location of cisbp pwm database<br/>
        --hocomoco location of hocomoco pwm database<br/>
        -o output folder<br/>
        --info logfile to write to <br/>
        --dummy dummy folder to use <br/>
#### ------
The output will be saved as a folder for each uniprot id provided (in the example given P35869). <br/>
Each of these folders will contain a single folder for each of the clusters detected. <br/>
Within the cluster folder will be found the memes for all of the cluster member PWMs as well as a mean PWM for the cluster. <br/>
#### ------

## Other available functionalities that may be of interest

The exe folder contains other executable programs that, while not the main focus of ModCRElib, may still be usefull. <br/>
- TFinderSelect.py<br/>
    Retrieves info from a pdb or uniprot entry to decide if a protein is a TF<br/>
- build_dna.py<br/>
    Turns a DNA string into a pdb file<br/>
- clean.py<br/>
    Cleans a PDB file<br/>
- contacts.py<br/>
    Gets the contacts calculated from a pdb file. <br/>
- dimers.py<br/>
    Check if input is a component of a dimer and retrieve monomer ids and contacts<br/>
- homologs.py<br/>
    given an input of a blast or hmm file will create a file of homologs<br/>
- interface.py<br/>
    Will search for the interface of interaction within a pdb<br/>
- merge_pwms.py<br/>
    Combine predicted pwms into a single averaged pwm<br/>
- mmcif_to_pdb.py<br/>
    convert mmcif format to pdb<br/>
- mmcifs_to_pdbs.py
    convert a folder of mmcif format files to pdb<br/>
- model_IMP.py<br/>
    Use IMP to model macro-complexes<br/>
- nearest_neighbour.py<br/>
    calculates the closest similar sequences (nearest neighbour) of each TF and compares their PWMs using TOMTOM. It also compares the PWMs modelled for each TF with the dataset of PWMs. Boxplots are built to compare both success using different conditions and scores derived from TOMTOM results<br/>
- pbm.py and pdb.py<br/>
    create the pbm and pdb directories that was provided at http://aleph.upf.edu/modcrefiles/<br/>
- pdb2thread.py<br/>
    converts pdb to a thread file<br/>

