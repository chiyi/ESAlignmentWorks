#!/bin/bash

. /afs/cern.ch/cms/caf/setup.sh
export VO_CMS_SW_DIR=/afs/cern.ch/project/gd/apps/cms
source $VO_CMS_SW_DIR/cmsset_default.sh

MYScriptPATH="/afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src/AlignmentTool/ESAlignTool/runVal/run4Corner/mergePlots"
MYVersionPATH="/afs/cern.ch/user/c/chiyi/scratch0/ESAlignment/CMSSW_5_3_6/src"

cd $MYVersionPATH
eval `scram runtime -sh`
cd -

cp $MYScriptPATH/ReadFiles_Run2012D.C ./
root -l -q -b ReadFiles_Run2012D.C+ > ReadFiles_Run2012D.log 2>&1
cp ReadFiles_Run2012D.log $MYScriptPATH
cmsStage AddedPlots.root /store/caf/user/chiyi/ESAlignment/2012D/run4Corner
rm *

