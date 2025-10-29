#!/bin/sh
##$BATCH -l h_vmem=200G
#$BATCH --mem=80G
#purge previous loads
module purge
#start loading the path
#module load modulepath
#load the modules
module load Python/3.8.6-GCCcore-10.2.0
module load ClustalW2/2.1-foss-2020b
module load Clustal-Omega/1.2.4-foss-2020b
#module load dssp/1.0 # /usr/bin/dssp
module load EMBOSS/6.6.0-foss-2020b
module load Modeller
module load MEME/5.1.1-GCCcore-10.2.0-Python-3.8.6
module load BLAST+
module load HMMER/3.3.2-foss-2020b
module load TMalign
#module load hbplus/3.2
module load x3dna
module load CD-HIT/4.8.1-GCC-10.2.0
#module load Rosetta/2016.20
#module load imp/2.21.0-GCC-10.2.0
module load Ghostscript/9.53.3-GCCcore-10.2.0

cd /home/pgohl/ModCRE_Package 

#python exe/get_json.py -home_path -Fasta_TF_Sequence -Uniprot_codes_TF -Family_Labels_TFs -Nearest_Neighbor_file -modcre_generated_pwm_folder -output_Folder_name -id_type
python exe/get_json.py /home/pgohl/ModCRE_Package tf_sequences.fa TF_codes_40-50.txt files/TF_accession_family.csv files/TF_nearest_neighbor_30-100.csv pwm JSON uniprot


