setenv dataset $1
setenv data_mc $2

cd UpperLimitToys_${dataset}_fit${data_mc}

foreach i ( `ls -1 | grep "[0-9][0-9][0-9]" | grep -v r ` )
cd $i/res
foreach j (`ls -1 outputToy*.tgz`)
tar -xvf $j
end
hadd -f mergedToys.root outputToy/higgsCombineTest.MarkovChainMC.*.root
cd -
end
cd ..
