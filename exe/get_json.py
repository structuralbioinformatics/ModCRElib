import sys,os
import pandas as pd
import json


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


# Append scripts path to python path #
sys.path.append(scripts_path)





home      = sys.argv[1]
fasta     = sys.argv[2]
nonredun  = sys.argv[3]
families  = sys.argv[4]
neighbors = sys.argv[5]
modcre    = sys.argv[6]
output    = sys.argv[7]
idformat  = sys.argv[8]


scripts_path = os.path.join(home,"scripts")
config_path  = os.path.join(scripts_path,"ModCRElib","configure")
sys.path.append(scripts_path)

databases_of_pwms={ 
                   "jaspar":os.path.join(home,"jaspar/jaspar_2024/PWMS"),
                   "cisbp":os.path.join(home,"pbm/CisBP_2019/pwms"),
                   "hocomoco":os.path.join(home,"hocomoco/pwms")
                  }

# Imports my functions #
from ModCRElib.beans import functions
from ModCRElib.msa import pwm_pbm as PWM


if not os.path.exists(output):
   os.makedirs(output)

ff=open(fasta,"r")
sequences=set()
for line in ff:
 if line.startswith(">"):
    header="|".join(line.split("|")[0:2])
    sequences.add(header.split("|")[1])
ff.close()


ff=open(nonredun,"r")
selected_acc=set()
for line in ff:
    selected_acc.add(line.rstrip())
ff.close()

set_pwms={}
for db_pwm in databases_of_pwms:
    set_pwms[db_pwm]=set()

acc_fam={}
fam_acc={}
acc_nn={}
tf_fam = pd.read_csv(families)
tf_nn  = pd.read_csv(neighbors)

#Removing this failsafe database check
#Patrick
#for db_pwm in databases_of_pwms:
#    for file in os.listdir(databases_of_pwms[db_pwm]):
#        set_pwms[db_pwm].add(file)

for tf in tf_fam.values:
    acc_fam.setdefault(tf[0],set(tf[1].split(",")))
    for x in tf[1].split(","):
        fam_acc.setdefault(x,set()).add(tf[0])
for tf in tf_nn.values:
    if len(tf[0].split("|"))>1:
       code=tf[0].split("|")[1]
    else:
       code=tf[0].split("|")[0]
    print("Getting nearest neighbor PWMs of %s"%code)
    for nn in tf[1].split(","):
        x=nn.split("|")
        db=x[0]
        seq=x[1]
        pid=x[-1]
        database=None
        if db=="js":database="jaspar"
        if db=="ho":database="hocomoco"
        if db=="bp":database="cisbp"
        if database is None: 
            print("Errror: missing database ",db)
            continue
        for pwm in x[2:-1]:
            #if pwm+".meme" in set_pwms[database]:
            #print("\t Add PWM %s homologs of %s (DB %s) ID: %.2f (Reference %s)"%(pwm,seq,database,100*float(pid),tf[0]))
            data={"database":database, "code":seq, "pwm":pwm+".meme", "similarity":float(pid)}
            acc_nn.setdefault(code,[]).append(data)
acc_modcre={}
print("Reading ModCRE database")
for file in os.listdir(modcre):
    if file.endswith("meme"):
       try:
         msa_obj = PWM.nMSA(os.path.join(modcre,file),file.rstrip(".meme"),"meme")
         binding_site_length = msa_obj.get_binding_site_length()
         if binding_site_length <= 0: 
            print("Skip %s (error: %s)"%(file,"missing binding fragment with 0 length"))
            continue
         if idformat == "uniprot":
            code = "_".join([ x for x in file.split(":")[0].split("_")[2:] ]).split("_")[0]
         else:
            code = "_".join([ x for x in file.split(":")[0].split("_")[2:] ])
         print("\t-Get %s <- %s"%(code,file))
         acc_modcre.setdefault(code,set()).add(file)
         """
         if len(file.split("_"))>2:
            print("\t-Get %s <- %s"%(file.split("_")[2],file))
            acc_modcre.setdefault(file.split("_")[2],set()).add(file)
         else:
            print("\t-Get %s <- %s"%(file.split(".")[0].split("_")[0].upper(),file))
            acc_modcre.setdefault(file.split(".")[0].split("_")[0].upper(),set()).add(file)
         """
       except Exception as e:
         print("Skip %s (error: %s)"%(file,e))

print("Building JSON folder")
for sele in sequences:   
  if  sele not in acc_modcre:
      print("\t- Skip %s (not selected for training/testing) "%sele)
      continue
  if sele not in selected_acc:
      print("\t- Skip %s (not selected for training/testing altough it has PWMs predicted) "%sele)
      continue
  print("\t- make JSON %s"%(sele))
  data=[]
  if sele in acc_nn:
    print("\t\t- Use NN: %s"%( str( [ x["code"]+"|"+x["database"] for x in acc_nn[sele] ] ) ))
    for nn in acc_nn[sele]:
        data.append(nn)
  if sele in acc_modcre:
    print("\t\t- Use ModCRE templates: %s"%( str( [ "_".join(motif.split(".")[-2].split("_")[-3:]) for motif in  acc_modcre[sele]] ) ))
    for motif in acc_modcre[sele]:
        code = "_".join(motif.split(".")[-2].split("_")[-3:])
        mdl_dict={"database":"modcre", "code":code, "pwm":motif, "similarity":0.0}
        data.append(mdl_dict)
  if len(data)>0:
     json_file=open(os.path.join(output,"%s"%(sele+".json")),"w")
     json.dump(data,json_file,indent=4)    

