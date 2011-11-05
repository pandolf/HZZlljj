setenv dataset $1

cd datacards_${dataset}

pwd

foreach i ( `ls -1 | grep "[0-9][0-9][0-9]" | grep -v r ` )
cd $i
echo "Combining datacards for mass: $i"
combineCards.py CMS_hzz2l2q_ee0b=hzz2l2q_ee0b.${i}.txt CMS_hzz2l2q_ee1b=hzz2l2q_ee1b.${i}.txt CMS_hzz2l2q_ee2b=hzz2l2q_ee2b.${i}.txt CMS_hzz2l2q_mm0b=hzz2l2q_mm0b.${i}.txt CMS_hzz2l2q_mm1b=hzz2l2q_mm1b.${i}.txt CMS_hzz2l2q_mm2b=hzz2l2q_mm2b.${i}.txt >! CMS_hzz2l2q_${i}_6channels.txt
text2workspace.py CMS_hzz2l2q_${i}_6channels.txt -b -m $i -o model.root
cd ..
end

cd ..
