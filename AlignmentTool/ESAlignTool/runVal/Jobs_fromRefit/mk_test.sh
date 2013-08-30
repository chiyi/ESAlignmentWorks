#!/bin/bash
#rfio://castorcms/?svcClass=cmscafuser&path=

ALLPARTAN=`ls /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames_fromRefit/FileNames_PARTA_*.txt | wc -l`
ALLPARTBN=`ls /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames_fromRefit/FileNames_PARTB_*.txt | wc -l`

for((i=1;i<=1;i++))
do
 cat AlignIter0_p0_Upper > AlignIter0_pA"$i"_UM
 cat /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames_fromRefit/FileNames_PARTA_"$i".txt >> AlignIter0_pA"$i"_UM
 cat AlignIter0_p0_Middle >> AlignIter0_pA"$i"_UM

done

#cp file.txt_Sol1 file.txt
#cp file2.txt_Sol1 file2.txt

#for((i=2;i<=7;i++))
#for((i=3;i<=6;i++))
for((i=1;i<=1;i++))
do
 for((j=1;j<=1;j++))
 do
  sed 's/ITERFIRST/'$i'/;s/ITERSECOND/'$j'/' ./run_iter0_pA0.sh > ./run_iter"$i"_pA$j.sh
  
 done

done


chmod +x *.sh
