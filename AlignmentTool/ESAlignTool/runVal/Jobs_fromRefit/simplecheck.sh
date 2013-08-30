#!/bin/bash

for((i=1;$i<=691;i++))
#for((i=1;$i<=969;i++))
do
 #N=`less run_iter2_pA"$i".log | grep Done | wc -l`
 #N=`less run_iter7_pB"$i".log | grep Done | wc -l`
 #N=`less OutRootFiles.txt | grep pB"$i".root | wc -l`
 N=`less OutRootFiles.txt | grep pA"$i".root | wc -l`
 if(($N!=1));then
  echo $i
 fi
done
