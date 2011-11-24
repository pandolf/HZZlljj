setenv dataset $1
setenv data_mc $2
setenv nbtags $3

if ( $# > 2 ) then
  cd UpperLimitToys_${dataset}_${nbtags}btag_fit${data_mc}
else
  cd UpperLimitToys_${dataset}_fit${data_mc}
endif

foreach i ( `ls -1 | grep "[0-9][0-9][0-9]" | grep -v r ` )
if ( $# > 2 ) then
  echo "Merging mass: ${i} (${nbtags} btag)" 
else
  echo "Merging mass: ${i}" 
endif
cd $i/res
foreach j (`ls -1 outputToy*.tgz`)
tar -xvf $j
end
hadd -f mergedToys.root outputToy/higgsCombineTest.MarkovChainMC.*.root
cd -
end
cd ..
