foreach i ( `ls -1 | grep "[0-9][0-9][0-9]_Run2011A_FULL_4" | grep -v r ` )
cd $i/res
foreach j (`ls -1 outputToy*.tgz`)
tar -xvf $j
end
hadd -f mergedToys.root outputToy/higgsCombineTest.MarkovChainMC.*.root
cd -
end
