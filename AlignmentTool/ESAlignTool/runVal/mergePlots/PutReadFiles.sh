#!/bin/bash

  cp ReadFiles_Upper ReadFiles_Run2012D.C

INDNUM=0
for FILE in list_rootfiles.txt
do
 while read line
 do
  READFILENAME=`echo $line | awk '{print $5}'`
  SIZE=`echo $line | awk '{print $2}'`

  if(($SIZE>13117)); then
   INDNUM=$(($INDNUM+1))
   echo " ind="$INDNUM";" >> ReadFiles_Run2012D.C
   echo " f[ind-1]=TFile::Open(\"root://eoscms//eos/cms"$READFILENAME"\");" >> ReadFiles_Run2012D.C
   cat ReadFiles_Middle >> ReadFiles_Run2012D.C
  fi
  #echo " 'rfio:"$READFILENAME"'," >> FileNames_PARTA_$IPART.txt
  #echo $SIZE

 done < $FILE
done

  cat ReadFiles_Lower >> ReadFiles_Run2012D.C
 
