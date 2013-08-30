#!/bin/bash
. /afs/cern.ch/cms/caf/setup.sh
#export VO_CMS_SW_DIR=/afs/cern.ch/cms/sw
export VO_CMS_SW_DIR=/afs/cern.ch/project/gd/apps/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
#alias lcgui='source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh'
#alias crabenv='source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh'

MYScriptPATH="/afs/cern.ch/user/c/chiyi/CMSSW_3_11_1_patch3/src/SkimTool/RunSkim_atCaf/jobs"
MYVersionPATH="/afs/cern.ch/user/c/chiyi/CMSSW_3_11_1_patch3/src"
MYOutputPATH="/store/caf/user/chiyi/3_11_1_p3/RunSkim_atCAF"


cd $MYVersionPATH
/afs/cern.ch/user/c/chiyi/CMSSW_3_8_5/src
eval `scramv1 ru -sh`
cd -

cmsRun $MYScriptPATH/getSkimLoose_PNUM.py
cmsStageOut EG_WZpar1_ESSkim_PNUM.root $MYOutputPATH/EG_WZpar1_ESSkim_PNUM.root
rm EG_WZpar1_ESSkim_PNUM.root

