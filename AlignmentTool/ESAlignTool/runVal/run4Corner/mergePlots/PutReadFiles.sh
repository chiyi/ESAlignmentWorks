#!/bin/bash

  cp ReadFiles_Upper ReadFiles_4Corner.C

for((INDNUM=1;$INDNUM<=20;INDNUM++))
do
  #READFILENAME=`echo $line | awk '{print $5}'`
  #SIZE=`echo $line | awk '{print $2}'`

  echo " ind="$INDNUM";" >> ReadFiles_4Corner.C

  READFILENAME="/store/group/alca_ecalcalib/ESalignment/5_3_6/run4CornerPlots/ESpF_4Corner_iter11_P"$INDNUM".root"
  echo " f_pF[ind-1]=TFile::Open(\"root://eoscms//eos/cms"$READFILENAME"\");" >> ReadFiles_4Corner.C
  cat ReadFiles_Middle_pF >> ReadFiles_4Corner.C

  READFILENAME="/store/group/alca_ecalcalib/ESalignment/5_3_6/run4CornerPlots/ESpR_4Corner_iter11_P"$INDNUM".root"
  echo " f_pR[ind-1]=TFile::Open(\"root://eoscms//eos/cms"$READFILENAME"\");" >> ReadFiles_4Corner.C
  cat ReadFiles_Middle_pR >> ReadFiles_4Corner.C

  READFILENAME="/store/group/alca_ecalcalib/ESalignment/5_3_6/run4CornerPlots/ESmF_4Corner_iter11_P"$INDNUM".root"
  echo " f_mF[ind-1]=TFile::Open(\"root://eoscms//eos/cms"$READFILENAME"\");" >> ReadFiles_4Corner.C
  cat ReadFiles_Middle_mF >> ReadFiles_4Corner.C

  READFILENAME="/store/group/alca_ecalcalib/ESalignment/5_3_6/run4CornerPlots/ESmR_4Corner_iter11_P"$INDNUM".root"
  echo " f_mR[ind-1]=TFile::Open(\"root://eoscms//eos/cms"$READFILENAME"\");" >> ReadFiles_4Corner.C
  cat ReadFiles_Middle_mR >> ReadFiles_4Corner.C

done

  cat ReadFiles_Lower >> ReadFiles_4Corner.C
 
