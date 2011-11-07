foreach i ( `ls -1 | grep "[0-9][0-9][0-9]_HR11" | grep -v r ` )
cd $i/res
foreach j (`ls -1 outputToy*.tgz`)
tar -xvf $j
end
hadd -f mergedToys.root outputToy/higgsCombineTest.MarkovChainMC.*.root
cd -
end
