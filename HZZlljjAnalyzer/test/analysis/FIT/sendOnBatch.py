#! /usr/bin/env python
import os
import sys
import time
import re

if (len(sys.argv) != 3) :
    print "usage sendOnBatch.py dataset nJobs"
    sys.exit(1)
dataset = sys.argv[1]
njobs = int(sys.argv[2])

queue="8nh"

diskoutputdir = "/cmsrm/pc22_2/pandolf/FitToys_"+dataset
diskoutputmain = diskoutputdir

dir = "/tmp/pandolf/FitToys_" + dataset
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm22 mkdir -p "+diskoutputmain)


pwd = os.environ['PWD']

application = "generateToysForParErrors"

ijob=0

toys_per_job = 5



while ( ijob < njobs ):

    # prepare the script to run
    outputname = dir+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
    outputfile.write('cd /afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_3_patch5/ ; eval `scramv1 runtime -sh` ; cd -\n')
    #outputfile.write('source /afs/cern.ch/sw/lcg/external/gcc/4.3/x86_64-slc5/setup.sh\n')
    #outputfile.write('source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.27.04/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh\n')
    outputfile.write('cp '+pwd+'/HZZlljjRM*.root $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+" "+dataset+" "+str(toys_per_job)+" "+str(ijob)+"\n")
    outputfile.write('rm HZZlljjRM*.root\n')
    outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm22:'+diskoutputmain+'/{}\n') 
    #outputfile.write('cp *.root '+diskoutputmain2+'\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
    ijob = ijob+1
    continue
