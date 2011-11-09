setenv dataset $1
setenv nbtags $2

cd datacards_${dataset}



foreach i ( `ls -1 | grep "[0-9][0-9][0-9]" | grep -v r ` )
  cd $i
  echo "Computing observed limit for mass: $i"
  if ( $# > 1 ) then
    combine model_${nbtags}btag.root -M MarkovChainMC -m $i -H ProfileLikelihood -U >&! log_${nbtags}btag.txt
  else
    combine model.root -M MarkovChainMC -m $i -H ProfileLikelihood -U >&! log.txt
  endif
  cd -
end
cd ..
