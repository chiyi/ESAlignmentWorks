#!/bin/bash

rm file.txt
rm file2.txt

for((i=1;i<=11;i=i+1))
do
 rm outputText.txt
 for((j=1;j<=2;j=j+1))
 do
  sed 's/ITERNUM/'$i'/;' ./AlignIter0_p"$j"_Upper > ./AlignIter"$i"_p"$j".py
  cat file.txt >> AlignIter"$i"_p"$j".py
  cat AlignIter0_p"$j"_Lower >> AlignIter"$i"_p"$j".py
  cmsRun AlignIter"$i"_p"$j".py > AlignIter"$i"_p"$j".log 2>&1
  cp AlignmentFile.root AlignmentFile_iter"$i"_p"$j".root
  root -l -q -b Read_elements_wRotation.C+
  echo p$j
  date
 done
 rm CalIter.C
 cat CalIter_wRotation_Upper > CalIter_wRotation.C
 cat outputText.txt >> CalIter_wRotation.C
 echo iterN="$i"';' >> CalIter_wRotation.C
 cat CalIter_wRotation_Lower >> CalIter_wRotation.C
 root -l -q -b CalIter_wRotation.C+
 cp CalIter_wRotation.C CalIter_wRotation$i.C
 echo iter$i
done

