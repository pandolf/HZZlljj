#! /bin/bash

DIR=$1

cd $DIR

combineCards.py CMS_hzz2l2q_ee0b=hzz2l2q_ee0b.${DIR}.txt CMS_hzz2l2q_ee1b=hzz2l2q_ee1b.${DIR}.txt CMS_hzz2l2q_ee2b=hzz2l2q_ee2b.${DIR}.txt CMS_hzz2l2q_mm0b=hzz2l2q_mm0b.${DIR}.txt CMS_hzz2l2q_mm1b=hzz2l2q_mm1b.${DIR}.txt CMS_hzz2l2q_mm2b=hzz2l2q_mm2b.${DIR}.txt > CMS_hzz2l2q_${DIR}_6channels.txt

cd -