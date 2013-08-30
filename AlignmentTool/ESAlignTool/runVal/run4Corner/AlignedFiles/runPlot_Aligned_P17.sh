#/bin/bash

. /afs/cern.ch/cms/caf/setup.sh
export VO_CMS_SW_DIR=/afs/cern.ch/project/gd/apps/cms
source $VO_CMS_SW_DIR/cmsset_default.sh

MYScriptPATH="/afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src/AlignmentTool/ESAlignTool/runVal/run4Corner/AlignedFiles"
MYVersionPATH="/afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src"
OutputPATH="/store/group/alca_ecalcalib/ESalignment/5_3_6/run4CornerPlots"

cd $MYVersionPATH
eval `scram ru -sh`
cd -

cp $MYScriptPATH/tdrstyle.C ./
cp $MYScriptPATH/ESAlign.* ./
cp $MYScriptPATH/runPlot.C ./

cmsStage /store/group/alca_ecalcalib/ESalignment/5_3_6/Run2012D_Validation/AlignmentFile_iter1_p6.root ./

root -l -q -b runPlot.C+ > runPlot.log 2>&1
cmsStage ESmF_4Corner_iter11.root $OutputPATH/ESmF_4Corner_iter11_P17.root
cmsStage ESmR_4Corner_iter11.root $OutputPATH/ESmR_4Corner_iter11_P17.root
cmsStage ESpF_4Corner_iter11.root $OutputPATH/ESpF_4Corner_iter11_P17.root 
cmsStage ESpR_4Corner_iter11.root $OutputPATH/ESpR_4Corner_iter11_P17.root 

rm $MYScriptPATH/runPlot_Aligned_P17.log
rm *

