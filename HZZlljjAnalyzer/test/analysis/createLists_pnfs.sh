#!/bin/bash -f
# usage: ./createLists.sh /pnfs/roma1.infn.it/data/cms/store/user/rahatlou/MC/41xv3DA/

maindir=$1

echo "scanning $maindir"

ls -l $maindir | awk '{print $9}' > datasets.txt
ls -l $maindir | awk '{print "'"$maindir"'" "/" $9}' | xargs -i echo "ls " {} " | grep -v \" 0 \" | awk '{print \"{}/\" \$9}'" > commands.txt 

for i in `ls $maindir`
do 
rm -f files_$i.txt
for j in `ls $maindir/$i`
do
echo dcap://cmsrm-se01.roma1.infn.it$maindir/$i/$j >> files_$i.txt
done
done

