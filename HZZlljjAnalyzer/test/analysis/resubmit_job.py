#! /usr/bin/env python
import sys
import os
import time

#######################################
### usage  python resubmit_job PDname SDname recoType jetAlgo jobNumber
#######################################
if len(sys.argv) != 3:
    print "usage python resubmit_job.py dataset jobNumber"
    sys.exit(1)
dataset_name = sys.argv[1]
jobNumber = int(sys.argv[2])
pwd = os.environ['PWD']
dir = "HZZlljj_"+dataset_name;
outputname = dir+"/src/submit_"+str(jobNumber)+".src"
os.system("bsub -q 8nh -o "+pwd+"/"+dir+"/log/"+dataset_name+"_"+str(jobNumber)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset_name+"_"+str(jobNumber))
#os.system("bsub -q 2nd -o "+pwd+"/"+dir+"/log/"+dataset_name+"_"+str(jobNumber)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset_name+"_"+str(jobNumber))
time.sleep(3.5) #to allow multiple callings
