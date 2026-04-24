#!/bin/sh
source env_setup.sh

input="JSON"
modcre="pwm"


pvalue="0.005"
threshold=$(echo "$pvalue * 4"| bc -l)
length="50"

python ModCRElib/msa/aggregate_pwms.py --complete 0.95 -i ${input} --pvalue=$pvalue --threshold=$threshold -l $length --jaspar=ExternalPWMs/jaspar_pwms --cisbp=ExternalPWMs/CisBP_pwms --hocomoco=ExternalPWMs/hocomoco_pwms --modcre=${modcre}  -o Aggregates/${input}_MAGG_PV${pvalue}_T${threshold}_L${length}_reference --verbose --info=Aggregates/info_${input}_MAGG_PV${pvalue}_T${threshold}_L${length}.log --reference 0.99 --trim --dummy=Aggregates/dummy_${input}_PV${pvalue}_T${threshold}_L${length}

#pvalue="0.005"
#threshold=$(echo "$pvalue * 10"| bc -l)
#length="50"

#python ModCRElib/msa/aggregate_pwms.py --complete 0.95 -i ${input} --pvalue=$pvalue --threshold=$threshold -l $length --jaspar=ExternalPWMs/jaspar_pwms --cisbp=ExternalPWMs/CisBP_pwms --hocomoco=ExternalPWMs/hocomoco_pwms --modcre=${modcre}  -o Aggregates/${input}_MAGG_PV${pvalue}_T${threshold}_L${length}_reference --verbose --info=Aggregates/info_${input}_MAGG_PV${pvalue}_T${threshold}_L${length}.log --reference 0.99 --trim --dummy=Aggregates/dummy_${input}_PV${pvalue}_T${threshold}_L${length}




