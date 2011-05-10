#!/bin/bash 
# usage: ./makeFileLists.sh [dirname]

maindir=$1
if [ $# -eq 1 ]
then
analyzerType="HZZlljj"
else
analyzerType=$2
fi

echo "scanning $maindir"

for i in `ls $maindir`
do
rm -f files_${analyzerType}_2ndLevel_$i.txt
for j in `ls $maindir/$i/${analyzerType}*`
do
echo $j >> files_${analyzerType}_2ndLevel_$i.txt
#echo $j
#echo $maindir/$i/$j 
done
done

