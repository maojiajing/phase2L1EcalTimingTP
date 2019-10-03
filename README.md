# phase2L1EcalTimingTP
Using the primary recipe from here:
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_10_6_1_patch2
```
cmsrel CMSSW_10_6_1_patch2
cd CMSSW_10_6_1_patch2/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline phase2-l1t-integration-CMSSW_10_6_1_patch2
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v2.22.1-CMSSW_10_6_1_patch2



git cms-addpkg L1Trigger/L1TCommon
git cms-addpkg DataFormats/EcalDigi
git cms-addpkg SimCalorimetry/EcalEBTrigPrimAlgos
git cms-addpkg SimCalorimetry/EcalEBTrigPrimProducers

cd L1Trigger
# to just clone the repo... though it is better if you fork it and clone from your own area!
git clone https://github.com/maojiajing/phase2L1EcalTimingTP.git

cd ../

USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable" scram b -j 8

cd L1Trigger/phase2L1EcalTimingTP/test

# reprocess L1 emulation on phase 2 TTbar or QCD sample
cmsRun reprocess_test_9_3_7.py

#Output can be found in
# TTbar
#'file:/afs/cern.ch/user/j/jmao/work/public/releases/L1_Trigger/MC_Production/FollowCIEMAT/b/CMSSW_10_6_1_patch2/src/step2_2ev_reprocess_slim_jet.root'

#QCD
#'file:/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/Phase2/Phase2L1/EcalTimimg/step2_reprocess_slim_qcd.root'

# test ECAL TP

cmsRun test-EcalTPG-jet.py
```

