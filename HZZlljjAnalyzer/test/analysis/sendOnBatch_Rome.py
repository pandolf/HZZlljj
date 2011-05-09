#! /usr/bin/env python
import os
import sys
import time
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if (len(sys.argv) != 3) and (len(sys.argv) != 4) and (len(sys.argv) != 5):
    print "usage sendOnBatch.py dataset filesPerJob analyzerType=\"HZZlljj\" flags=\"\""
    sys.exit(1)
dataset = sys.argv[1]
inputlist = "files_"+dataset+".txt"
#settingfile = "config/RSZZsettings.txt"
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = "cmslong"
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
#diskoutputdir = "/cmsrm/pc21_2/pandolf/MC/"+dataset
diskoutputdir = "/cmshome/pandolf/STORE/MC/Spring11/"+dataset
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
    os.system("mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0

#copy the configuration in the actual run directory
#os.system("cp -r config "+dataset_name)

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+dir+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(ijobmax,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)


    inputfile.close()

    # prepare the script to run
    outputname = dir+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('cd /cmshome/pandolf/CMSSW_4_1_5/src ; eval `scramv1 runtime -sh` ; cd -\n')
    #    outputfile.write('cd '+pwd)
    #outputfile.write('cd $WORKDIR\n')
    outputfile.write("mkdir -p "+dir+"/res/job_"+str(ijob)+"\n")
    outputfile.write("cp Cert* "+dir+"/res/job_"+str(ijob)+"\n")
    #outputfile.write("cp "+application+" "+dir+"/res/job_"+str(ijob)+"\n")
    outputfile.write("cd "+dir+"/res/job_"+str(ijob)+"\n")
    outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+str(ijob)+"\n")
    outputfile.write('mv *.root '+diskoutputmain+'\n')
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
    ijob = ijob+1

    continue
