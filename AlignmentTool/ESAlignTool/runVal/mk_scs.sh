#!/bin/bash
#rfio://castorcms/?svcClass=cmscafuser&path=

awk '{ print "\x27rfio://castorcms/?svcClass=cmscafuser&path=/castor/cern.ch/cms/store/caf/user/chiyi/3_11_1_p3/RunSkim_atCAF/" $9 "\x27" }' ./filelist_inp.txt > ./test.log

NL=0
NLMAX=`less test.log | wc -l`

i=1
cat test.log | while read line
do
 cat AlignIter0_p0_Upper > AlignIter0_p"$i"_UM
 echo $line >> AlignIter0_p"$i"_UM
 cat AlignIter0_p0_Middle >> AlignIter0_p"$i"_UM

 i=$(($i+1))

done

for((i=1;i<=1;i++))
do
 rm command_iter$i
 for((j=1;j<=$NLMAX;j++))
 do
  sed 's/ITERFIRST/'$i'/;s/ITERSECOND/'$j'/' ./run_iter0_p0.sh > ./run_iter"$i"_p$j.sh

  echo bsub -q cmscaf1nh -o "$PWD"/run_iter"$i"_p"$j".log \""$PWD"/./run_iter"$i"_p$j.sh\" >> command_iter$i
 done

 echo \#bsub -q cmscaf1nh -o "$PWD"/run_iter"$i".log \""$PWD"/./run_iter"$i".sh\" >> command_iter$i
done


#for((i=1;i<=1;i++))
#do
# sed 's/ITERFIRST/'$i'/;' ./run_iter0.sh > ./run_iter"$i".sh
#
#done

chmod +x *.sh
