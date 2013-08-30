#!/bin/bash
. /afs/cern.ch/cms/caf/setup.sh
export VO_CMS_SW_DIR=/afs/cern.ch/project/gd/apps/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
#alias lcgui='source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh'
#alias crabenv='source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh'

MYScriptPATH="/afs/cern.ch/user/c/chiyi/CMSSW_3_11_1_patch3/src/AlignmentTool/ESAlignTool/runAlign"
MYVersionPATH="/afs/cern.ch/user/c/chiyi/CMSSW_3_11_1_patch3/src"
MYOutputPATH="/store/caf/user/chiyi/3_11_1_p3/runLog_atCAF"


cd $MYVersionPATH
eval `scramv1 ru -sh`
cd -
cp $MYScriptPATH/ESAlign.* ./
cp $MYScriptPATH/Read_elements_wRotation.C ./

touch $MYScriptPATH/file.txt
touch $MYScriptPATH/file2.txt
cp $MYScriptPATH/file.txt ./
cp $MYScriptPATH/file2.txt ./
touch outputText.txt

i=ITERFIRST

 for((j=1;j<=5;j++))
 do
  cat $MYScriptPATH/outputText_iter"$i"_p"$j".txt >> outputText.txt

  #date >> $MYScriptPATH/process.log
 done

 rm CalIter_wRotation.C
 cat $MYScriptPATH/CalIter_wRotation_Upper > CalIter_wRotation.C
 cat outputText.txt >> CalIter_wRotation.C
 echo iterN="$i"';' >> CalIter_wRotation.C
 cat $MYScriptPATH/CalIter_wRotation_Lower >> CalIter_wRotation.C
 root -l -q -b CalIter_wRotation.C+
 cmsStageOut CalIter_wRotation.C $MYOutputPATH/CalIter"$i"_wRotation.C  >> run_message.log

 #cmsStageOut file.txt $MYOutputPATH/file.txt  >> run_message.log
 #cmsStageOut file2.txt $MYOutputPATH/file2.txt  >> run_message.log
 cp file.txt $MYScriptPATH/
 cp file2.txt $MYScriptPATH/

 rm outputText.txt
 rm file.txt
 rm file2.txt
 rm CalIter_wRotation*.*
 rm Read_elements_wRotation*.*
 rm ESAlign*.*
 #cp run_message.log $MYScriptPATH/
 rm run_message.log

