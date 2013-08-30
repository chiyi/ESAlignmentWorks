#!/bin/bash

rm command
INDNUM=0
for FILE in OutRootFiles.txt
do
 while read line
 do
  READFILENAME=`echo $line | awk '{print $5}'`
  READFILENAME=`echo "$READFILENAME" | sed 's:[]\[\^\$\.\*\/]:\\\\&:g'`
  SIZE=`echo $line | awk '{print $2}'`

  if(( $SIZE > "13117" )); then
   INDNUM=$(($INDNUM+1))
   sed 's/INPUTFILE/'$READFILENAME'/g;s/NUM/'$INDNUM'/g' runPlot_Aligned_P0.sh > runPlot_Aligned_P$INDNUM.sh

   echo bsub -q cmscaf1nh -o "$PWD"/runPlot_Aligned_P"$INDNUM".log \""$PWD"/runPlot_Aligned_P$INDNUM.sh\" >> command

  fi

 done < $FILE
done

chmod +x *.sh
