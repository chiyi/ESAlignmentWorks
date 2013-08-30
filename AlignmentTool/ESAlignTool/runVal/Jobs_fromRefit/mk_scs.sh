#!/bin/bash
#rfio://castorcms/?svcClass=cmscafuser&path=

ALLPARTN=`ls /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames_fromRefit/FileNames_PART_*.txt | wc -l`

for((i=1;i<=$ALLPARTN;i++))
do
 cat AlignIter0_p0_Upper > AlignIter0_p"$i"_UM
 cat /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames_fromRefit/FileNames_PART_"$i".txt >> AlignIter0_p"$i"_UM
 cat AlignIter0_p0_Middle >> AlignIter0_p"$i"_UM

done


#cp file.txt_Sol1 file.txt
#cp file2.txt_Sol1 file2.txt

#for((i=2;i<=7;i++))
#for((i=3;i<=6;i++))
for((i=1;i<=1;i++))
do
 rm command_iter$i
 for((j=1;j<=$ALLPARTN;j++))
 do
  sed 's/ITERFIRST/'$i'/;s/ITERSECOND/'$j'/' ./run_iter0_p0.sh > ./run_iter"$i"_p$j.sh
  
  echo bsub -q cmscaf1nh -o "$PWD"/run_iter"$i"_p"$j".log \""$PWD"/./run_iter"$i"_p$j.sh\" >> command_iter$i
 done

 #for((j=1;j<=$ALLPARTBN;j++))
 #do
 #  sed 's/ITERFIRST/'$i'/;s/ITERSECOND/'$j'/' ./run_iter0_pB0.sh > ./run_iter"$i"_pB$j.sh
   
 # #if(( $j==433 || $j==592 || $j==697 || $j==717 || $j==780 )); then
 # if(( $j==433 )); then
 #  echo $j skip
 #  echo "#"bsub -q 1nh -o "$PWD"/run_iter"$i"_pB"$j".log \""$PWD"/./run_iter"$i"_pB$j.sh\" >> command_iter$i
 # else
 #  echo bsub -q 1nh -o "$PWD"/run_iter"$i"_pB"$j".log \""$PWD"/./run_iter"$i"_pB$j.sh\" >> command_iter$i
 # fi
 #done

 echo \#bsub -q 1nh -o "$PWD"/run_iter"$i".log \""$PWD"/./run_iter"$i".sh\" >> command_iter$i
done


#for((i=4;i<=7;i++))
#do
# sed 's/ITERFIRST/'$i'/;' ./run_iter0.sh > ./run_iter"$i".sh
#
#done
#
chmod +x *.sh
