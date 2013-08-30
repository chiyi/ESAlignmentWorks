#!/bin/bash

for((i=101;$i<=1000;i++))
do
 N=`less OutRootFiles_4Corner.txt | grep ESmF_4Corner_iter11_P"$i".root | wc -l`
 if(($N!=1));then
  echo $i
 fi
done
