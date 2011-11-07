setenv putype $1

cd datacards_${putype}

foreach i ( `ls -1 | grep "[0-9][0-9][0-9]" | grep -v r ` )
  cd $i
  echo "Computing observed limit for mass: $i"
  combine model.root -M MarkovChainMC -m $i -H ProfileLikelihood -U >&! log.txt
  cd -
end
cd ..
