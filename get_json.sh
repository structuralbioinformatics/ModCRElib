#!/bin/sh
source env_setup.sh

#cd ./ 

#python exe/get_json.py -home_path_for_scripts -Fasta_TF_Sequence -Uniprot_codes_TF -Family_Labels_TFs -Nearest_Neighbor_file -modcre_generated_pwm_folder -output_Folder_name -id_type
python exe/get_json.py ./ tf_sequences.fa TF_codes_40-50.txt files/TF_accession_family.csv files/TF_nearest_neighbor_30-100.csv pwm JSON uniprot


