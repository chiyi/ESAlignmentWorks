#!/bin/bash
#rfio://castorcms/?svcClass=cmscafuser&path=
awk '{ print "\x27root://eoscms//eos/cms" $1 "?svcClass=default\x27" }' ./filelist_inp.txt > ./test.log

NL=0
NLMAX=`less test.log | wc -l`

i=1
cat test.log | while read line
do
  NL=$(($NL+1))
  #rm getSkimLoose_0_"$i".py
  cat getSkimLoose_Upper > getSkimLoose_0_"$i".py
  echo $line >> getSkimLoose_0_"$i".py
  cat getSkimLoose_Lower >> getSkimLoose_0_"$i".py

  sed 's/PNUM/P'$i'/' ./getSkimLoose_0_"$i".py > getSkimLoose_P"$i".py
  rm getSkimLoose_0_"$i".py

  sed 's/PNUM/P'$i'/g' ./runSkim_p0.sh > ./runSkim_P"$i".sh
  chmod +x runSkim_P"$i".sh

  echo bsub -q cmscaf1nh -o "$PWD"/runSkim_P"$i".log \""$PWD"/runSkim_P"$i".sh\" >> command

  i=$(($i+1))
done


