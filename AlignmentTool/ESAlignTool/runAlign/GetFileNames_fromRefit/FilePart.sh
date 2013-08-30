#!/bin/bash

for((iPART=1;$iPART<=22;iPART++))
do
 line=`less ChkNameNSize.log | grep ESSkim_wRefit_p"$iPART".root`

 READFILENAME=`echo $line | awk '{print $5}'`
 N=`echo $READFILENAME | grep ESSkim | wc -l`

 #echo " 'rfio:"$READFILENAME"'," > FileNames_PART_$iPART.txt
 if(($N==0)); then
  echo $iPART
  echo "" > FileNames_PART_$iPART.txt
 else
  #echo " 'rfio:"$READFILENAME"'," > FileNames_PART_$iPART.txt
  echo " 'root://eoscms//eos/cms"$READFILENAME"'," > FileNames_PART_$iPART.txt
 fi

done

#IPART=1
#TOTALSIZE=0
#for FILE in ChkNameNSize_B.log
#do
# while read line
# do
#  READFILENAME=`echo $line | awk '{print $9}'`
#  SIZE=`echo $line | awk '{print $5}'`
#
#  TOTALSIZE=$(( $TOTALSIZE + $SIZE ))
#  if(($TOTALSIZE>16000000000)); then
#   IPART=$(($IPART+1))
#   TOTALSIZE=$SIZE
#  fi
#  echo " 'rfio:"$READFILENAME"'," >> FileNames_PART_$IPART.txt
#  echo $TOTALSIZE
#
#
# done < $FILE
#done
