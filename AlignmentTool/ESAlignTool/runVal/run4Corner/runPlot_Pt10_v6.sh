#/bin/bash

. /afs/cern.ch/cms/caf/setup.sh
export VO_CMS_SW_DIR=/afs/cern.ch/project/gd/apps/cms
source $VO_CMS_SW_DIR/cmsset_default.sh

MYScriptPATH="/afs/cern.ch/user/c/chiyi/CMSSW_4_4_0/src/Resolution/run4Corner"
MYVersionPATH="/afs/cern.ch/user/c/chiyi/CMSSW_4_4_0/src"

cd $MYVersionPATH
eval `scramv1 ru -sh`
cd -

cp $MYScriptPATH/tdrstyle.C ./
cp $MYScriptPATH/ESAlign.* ./
cp $MYScriptPATH/runPlot_Pt10_v6.C ./
root -l -q -b runPlot_Pt10_v6.C+ > runPlot_Pt10_v6.log 2>&1
rfcp ESmF_4Corner_iter1.root /castor/cern.ch/user/c/chiyi/4_4_0/ESAlign_runVal/ValPlots/ESmF_4Corner_iter1_Pt10_v6.root
rfcp ESmR_4Corner_iter1.root /castor/cern.ch/user/c/chiyi/4_4_0/ESAlign_runVal/ValPlots/ESmR_4Corner_iter1_Pt10_v6.root
rfcp ESpF_4Corner_iter1.root /castor/cern.ch/user/c/chiyi/4_4_0/ESAlign_runVal/ValPlots/ESpF_4Corner_iter1_Pt10_v6.root
rfcp ESpR_4Corner_iter1.root /castor/cern.ch/user/c/chiyi/4_4_0/ESAlign_runVal/ValPlots/ESpR_4Corner_iter1_Pt10_v6.root
cp runPlot_Pt10_v6.log $MYScriptPATH

rm *

