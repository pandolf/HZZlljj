#! /usr/bin/env python
import sys
import os
import time

#######################################
### usage  python resubmit_job PDname SDname recoType jetAlgo jobNumber
#######################################
if len(sys.argv) != 6:
    print "usage python resubmit_job.py PDname SDname recoType jetAlgo jobNumber"
    print "example : resubmit_job.py QCD_Spring10 Pt80 pf akt5 10"
    print "example : resubmit_job.py EG Run2010A-PromptReco-v4 calo kt6 5"
    sys.exit(1)
PDname = sys.argv[1]
SDname = sys.argv[2]
recoType = sys.argv[3]
jetAlgo = sys.argv[4]
jobNumber = int(sys.argv[5])
dataset_name = PDname + "_" + SDname
pwd = os.environ['PWD']
outputname = dataset_name+"_"+recoType+jetAlgo+"/src/submit_"+str(jobNumber)+".src"
os.system("bsub -q 8nh -o "+pwd+"/"+dataset_name+"_"+recoType+jetAlgo+"/log/"+dataset_name+"_"+str(jobNumber)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset_name+"_"+str(jobNumber))
time.sleep(3.5) #to allow multiple callings
