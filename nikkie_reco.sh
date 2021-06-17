#!/bin/bash

##
## Check these settings below in order to adjust to the actual environment

runNumber=$1
CMSSW_DIR="/home/chromie/nikkie_tracking_analysis/test_comp/CMSSW_10_1_2/src/Geometry/PixelTelescope/test/"
CMSSW_RECO_CFG="PixelTelescope_BeamData_RECO_tracking.py"

RAW_DATA_RECO="/eos/cms/store/group/dpg_tracker_upgrade/BeamTestTelescope/2018-TELESCOPE-COMMISSIONING/raw/run${runNumber}"
OUTFILE="/eos/cms/store/group/dpg_tracker_upgrade/BeamTestTelescope/2019-NIKKIES_ANALYSIS-DONT_TOUCH/reco_run_${runNumber}.root"
OUT_DQM="/eos/cms/store/group/dpg_tracker_upgrade/BeamTestTelescope/2019-NIKKIES_ANALYSIS-DONT_TOUCH/dqm_run_${runNumber}.root"

cd $RAW_DATA_RECO

if [ -f *complete ]; then
  echo "File found!"
  echo "Nikkie needs to remove the run_complete file :)"
  #rm *complete
  exit 1  
fi

source /cvmfs/cms.cern.ch/cmsset_default.sh

pwd
cd $CMSSW_DIR
pwd
eval `scramv1 runtime -sh`

echo "cmsRun $CMSSW_DIR/$CMSSW_RECO_CFG inputDir=$RAW_DATA_RECO outputFileName=$OUTFILE outputDQMFileName=$OUT_DQM"
cmsRun $CMSSW_DIR/$CMSSW_RECO_CFG inputDir=$RAW_DATA_RECO outputFileName=$OUTFILE outputDQMFileName=$OUT_DQM 2>&1 > ${RAW_DATA_RECO}/blah.root.log
