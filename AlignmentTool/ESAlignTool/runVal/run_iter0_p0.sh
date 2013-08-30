#!/bin/bash
. /afs/cern.ch/cms/caf/setup.sh
export VO_CMS_SW_DIR=/afs/cern.ch/project/gd/apps/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
#alias lcgui='source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh'
#alias crabenv='source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh'

MYScriptPATH="/afs/cern.ch/user/c/chiyi/CMSSW_3_11_1_patch3/src/AlignmentTool/ESAlignTool/runVal"
MYVersionPATH="/afs/cern.ch/user/c/chiyi/CMSSW_3_11_1_patch3/src"
MYOutputPATH="/store/caf/user/chiyi/3_11_1_p3/runValLog_atCAF"


cd $MYVersionPATH
eval `scramv1 ru -sh`
cd -
#cp $MYScriptPATH/ESAlign.* ./
#cp $MYScriptPATH/Read_elements_wRotation.C ./

#rm file.txt
#rm file2.txt
touch $MYScriptPATH/file.txt

i=ITERFIRST
 j=ITERSECOND

  date >> $MYScriptPATH/process.log
  #echo Begin_iter"$i"_p"$j" >> $MYScriptPATH/process.log

  sed 's/ITERNUM/'$i'/;' $MYScriptPATH/AlignIter0_p"$j"_UM > ./AlignIter"$i"_p"$j".py
  cat $MYScriptPATH/file.txt >> AlignIter"$i"_p"$j".py
  cat $MYScriptPATH/Para.txt >> AlignIter"$i"_p"$j".py
  cat $MYScriptPATH/AlignIter0_p0_Lower >> AlignIter"$i"_p"$j".py
  cmsRun AlignIter"$i"_p"$j".py > AlignIter"$i"_p"$j".log 2>&1
  #cp AlignIter"$i"_p"$j".log $MYScriptPATH/
  cmsStageOut AlignIter"$i"_p"$j".log $MYOutputPATH/AlignIter"$i"_p"$j".log >> run_message.log
  rm AlignIter"$i"_p"$j".log
  cmsStageOut AlignmentFile.root $MYOutputPATH/AlignmentFile_iter"$i"_p"$j".root >> run_message.log
  #root -l -q -b Read_elements_wRotation.C+
  #cp outputText.txt $MYScriptPATH/outputText_iter"$i"_p"$j".txt
  echo p$j
  rm AlignmentFile.root
  date

  date >> $MYScriptPATH/process.log
  echo Done_iter"$i"_p"$j" >> $MYScriptPATH/process.log

#  if (( $i!="1" && $i!="11" )); then
#   rm AlignmentFile_iter"$i"_p"$j".root
#  fi

 #rm CalIter_wRotation.C
 #cat $MYScriptPATH/CalIter_wRotation_Upper > CalIter_wRotation.C
 #cat outputText.txt >> CalIter_wRotation.C
 #echo iterN="$i"';' >> CalIter_wRotation.C
 #cat $MYScriptPATH/CalIter_wRotation_Lower >> CalIter_wRotation.C
 #root -l -q -b CalIter_wRotation.C+
 #cmsStageOut CalIter_wRotation.C $MYOutputPATH/CalIter"$i"_wRotation.C  >> run_message.log

 #cmsStageOut file.txt $MYOutputPATH/file.txt  >> run_message.log
 #cmsStageOut file2.txt $MYOutputPATH/file2.txt  >> run_message.log

 rm outputText.txt
 #rm file.txt
 #rm file2.txt
 #rm CalIter_wRotation*.*
 #rm Read_elements_wRotation*.*
 #rm ESAlign*.*
 #rm AlignmentFile.root
 #cp run_message.log $MYScriptPATH/
 rm run_message.log

