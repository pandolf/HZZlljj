#! /usr/bin/env python
import os
import sys
import re

if ( (len(sys.argv) != 2) and (len(sys.argv) != 3) ) :
    print "usage compute_observedLimit_batch.py [dataset] [nbtags=\"ALL\"]"
    sys.exit(1)
dataset = sys.argv[1]

nbtags = "ALL"
if len(sys.argv) == 3:
    nbtags = sys.argv[2]


queue="8nh"

pwd = os.environ['PWD']

datacards_dir = "datacards_" + dataset


massesFile = open('masses_tmp.txt', 'r')


for line in massesFile:

  mass =  line.strip('\n')
  massDir = datacards_dir + "/" + str(mass)
  scriptName = massDir + "/batchScript.src"
  diskoutputmain = '/cmsrm/pc18/pandolf/CMSSW_4_2_8/src/HZZlljj/HZZlljjAnalyzer/test/analysis/FIT/' + massDir
  scriptFile = open(scriptName,'w')
  scriptFile.write('#!/bin/bash\n')
  scriptFile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
  scriptFile.write('export LANGUAGE="C"\n')
  scriptFile.write('export LC_ALL="C"\n')
  scriptFile.write('cd /afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_3_patch5/ ; eval `scramv1 runtime -sh` ; cd -\n')
  scriptFile.write('cp '+pwd+"/"+massDir+'/model*.root $WORKDIR\n')
  scriptFile.write('cd $WORKDIR\n')
  scriptFile.write('echo "Computing upper limit for mass: ' + str(mass) + '\n')
  scriptFile.write('combine model.root -M MarkovChainMC -m ' + str(mass) + ' -H ProfileLikelihood -U >& log.txt\n')
  scriptFile.write('ls log*.txt | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm22:'+diskoutputmain+'/{}\n') 
  scriptFile.close
  os.system("echo bsub -q "+queue+" -o "+pwd+"/log.log source "+pwd+"/"+scriptName)
  os.system("bsub -q "+queue+" -o "+pwd+"/log.log source "+pwd+"/"+scriptName)
  continue



