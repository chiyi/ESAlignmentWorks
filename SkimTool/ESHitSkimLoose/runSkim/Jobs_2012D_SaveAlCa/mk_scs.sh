#!/bin/bash
#rfio://castorcms/?svcClass=cmscafuser&path=

#ALLPARTAN=`ls /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames/FileNames_PARTA_*.txt | wc -l`
#ALLPARTBN=`ls /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames/FileNames_PARTB_*.txt | wc -l`

ALLPARTN=`ls /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames_2012D/FileNames_PART_*.txt | wc -l`


for((i=1;i<=$ALLPARTN;i++))
do
 cat getSkimLoose_Upper > getSkimLoose_p"$i".py
 #cat /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames/FileNames_PARTA_"$i".txt >> getSkimLoose_p"$i".py
 cat /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames_2012D/FileNames_PART_"$i".txt >> getSkimLoose_p"$i".py
 cat getSkimLoose_Lower >> getSkimLoose_p"$i".py

done
#
#for((i=1;i<=$ALLPARTBN;i++))
#do
# cat getSkimLoose_Upper > getSkimLoose_pB"$i".py
# cat /afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_2_patch4/src/AlignmentTool/ESAlignTool/runAlign/GetFileNames/FileNames_PARTB_"$i".txt >> getSkimLoose_pB"$i".py
# cat getSkimLoose_Lower >> getSkimLoose_pB"$i".py
#
#done

 rm command

 for((j=1;j<=$ALLPARTN;j++))
 do
  sed 's/ITERFIRST/'$i'/;s/NUM/'$j'/g' ./runSkim_p0.sh > ./runSkim_p$j.sh

  echo bsub -q 1nd -o "$PWD"/runSkim_p"$j".log \""$PWD"/./runSkim_p$j.sh\" >> command
 done


chmod +x *.sh
