#! /usr/bin/env python
import os
import sys
import time
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if (len(sys.argv) != 3) and (len(sys.argv) != 4) and (len(sys.argv) != 5):
    print "usage sendOnBatch.py dataset njobs analyzerType=\"HZZlljj\" flags=\"\""
    sys.exit(1)
dataset = sys.argv[1]
inputlist = "files_HZZlljj_"+dataset+".txt"
#settingfile = "config/RSZZsettings.txt"
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = "8nh"
#ijobmax = 40
ijobmax = int(sys.argv[2])
#application = "VecbosApp"
analyzerType = "HZZlljj"
if len(sys.argv) == 4:
    analyzerType = sys.argv[3]
flags = ""
if len(sys.argv) == 5:
    flags = sys.argv[4]
application = "do2ndLevel_"+analyzerType
if flags=="400":
    application = "do2ndLevel_TMVA_400"
if flags=="500":
    application = "do2ndLevel_TMVA_500"
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################
diskoutputdir = "/cmsrm/pc21_2/pandolf/MC/"+dataset
#outputmain = castordir
diskoutputmain = diskoutputdir
# prepare job to write on the cmst3 cluster disks
################################################
dir = analyzerType + "_" + dataset
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")
#outputroot = outputmain+"/root/"
#if castordir != "none": 
#    os.system("rfmkdir -p "+outputmain)
#    os.system("rfmkdir -p "+outputroot)
#    os.system("rfchmod 777 "+outputmain)
#    os.system("rfchmod 777 "+outputroot)
#else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
filesperjob = numfiles/ijobmax
extrafiles  = numfiles%ijobmax
input = open(inputlist)
######################################

for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = pwd+"/"+dir+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    # if it is a normal job get filesperjob lines
    if ijob != (ijobmax-1):
        for line in range(filesperjob):
            ntpfile = input.readline() 
            inputfile.write(ntpfile)
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        while ntpfile != '':
            inputfile.write(ntpfile)
            ntpfile = input.readline()
            continue
    inputfile.close()

    # prepare the script to run
    outputname = dir+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    #    outputfile.write('cd '+pwd)
    #outputfile.write('cp '+pwd+'/Cert_132440-140399_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt $WORKDIR\n')
    #outputfile.write('cp '+pwd+'/lumi_by_LS_132440_140401.csv $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    #outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" _"+str(ijob)+"\n")
    outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+str(ijob)+"\n")
    outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm21:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
    ijob = ijob+1
    time.sleep(4.)
    continue
