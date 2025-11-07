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

#cd ./

input="JSON"
modcre="pwm"


pvalue="0.005"
threshold=$(echo "$pvalue * 4"| bc -l)
length="50"

python ModCRElib/msa/aggregate_pwms.py --complete 0.95 -i ${input} --pvalue=$pvalue --threshold=$threshold -l $length --jaspar=ExternalPWMs/jaspar_pwms --cisbp=ExternalPWMs/CisBP_pwms --hocomoco=ExternalPWMs/hocomoco_pwms --modcre=${modcre}  -o Aggregates/${input}_MAGG_PV${pvalue}_T${threshold}_L${length}_reference --verbose --info=Aggregates/info_${input}_MAGG_PV${pvalue}_T${threshold}_L${length}.log --reference 0.99 --trim --dummy=Aggregates/dummy_${input}_PV${pvalue}_T${threshold}_L${length}

pvalue="0.005"
threshold=$(echo "$pvalue * 10"| bc -l)
length="50"

python ModCRElib/msa/aggregate_pwms.py --complete 0.95 -i ${input} --pvalue=$pvalue --threshold=$threshold -l $length --jaspar=ExternalPWMs/jaspar_pwms --cisbp=ExternalPWMs/CisBP_pwms --hocomoco=ExternalPWMs/hocomoco_pwms --modcre=${modcre}  -o Aggregates/${input}_MAGG_PV${pvalue}_T${threshold}_L${length}_reference --verbose --info=Aggregates/info_${input}_MAGG_PV${pvalue}_T${threshold}_L${length}.log --reference 0.99 --trim --dummy=Aggregates/dummy_${input}_PV${pvalue}_T${threshold}_L${length}




