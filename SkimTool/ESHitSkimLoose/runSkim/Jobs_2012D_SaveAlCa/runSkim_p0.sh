#!/bin/bash
. /afs/cern.ch/cms/caf/setup.sh
#export VO_CMS_SW_DIR=/afs/cern.ch/cms/sw
export VO_CMS_SW_DIR=/afs/cern.ch/project/gd/apps/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
#alias lcgui='source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh'
#alias crabenv='source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh'

MYScriptPATH="/afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src/SkimTool/ESHitSkimLoose/runSkim/Jobs_2012D_SaveAlCa"
MYVersionPATH="/afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src"
#MYOutputPATH="/store/user/chiyi/5_3_62012_SkimRefit"
MYOutputPATH="/store/group/alca_ecalcalib/ESalignment/5_3_6/Run2012D_SkimRefit"
#MYOutputPATH="/store/caf/user/pchen/chiyi/5_3_62012_SkimRefit"


cd $MYVersionPATH
eval `scram ru -sh`
cd -

#cp $MYScriptPATH/Cert_190456-196531_8TeV_29Jun2012ReReco_Collisions12_JSON.txt ./


cmsRun $MYScriptPATH/getSkimLoose_pNUM.py > tmp.log 2>&1
cmsStage ESSkim_wRefit.root $MYOutputPATH/ESSkim_wRefit_pNUM.root 
rm ESSkim_wRefit.root
rm tmp.log
