# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: repr --processName=REPR --python_filename=reprocess_test_10_5_0_pre1.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 2 --era Phase2 --eventcontent FEVTDEBUGHLT --filein root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step2_2ev_reprocess_slim.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('REPR2',eras.Phase2C4_trigger)
#process = cms.Process('REPR',eras.Phase2C4_trigger)
#process = cms.Process('REPR',eras.Phase2C4_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #'/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_TuneCP5_pythia8/MINIAODSIM/NoPU_103X_upgrade2023_realistic_v2-v2/90000/F063E7F2-A71C-E14C-BF02-665560C6B3CC.root'),
                                #'root://xrootd-cms.infn.it//store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_TuneCP5_pythia8/MINIAODSIM/NoPU_103X_upgrade2023_realistic_v2-v2/90000/F063E7F2-A71C-E14C-BF02-665560C6B3CC.root'),

#                                'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_TuneCP5_pythia8/MINIAODSIM/NoPU_103X_upgrade2023_realistic_v2-v2/90000/F063E7F2-A71C-E14C-BF02-665560C6B3CC.root'),
				  'file:/afs/cern.ch/user/j/jmao/work/public/releases/L1_Trigger/MC_Production/FollowCIEMAT/CMSSW_10_6_0_pre4/src/L1Trigger/phase2L1EcalTimingTP/test/testEvent.root'),
                            
                            secondaryFileNames = cms.untracked.vstring(
                                #'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/30B662A6-E4FD-644A-B124-29107C17ECA1.root',
                                #'root://xrootd-cms.infn.it//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/30B662A6-E4FD-644A-B124-29107C17ECA1.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/30B662A6-E4FD-644A-B124-29107C17ECA1.root',
                            )
                        )


process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('repr nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step2_2ev_reprocess_slim.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1TFullEvent.root"),
    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.out+process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)


process.source.secondaryFileNames = cms.untracked.vstring(
 #"root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext5-v1/40000/FF240C40-5C36-AC43-B89A-E578BB58E32F.root")
				  'file:/afs/cern.ch/user/j/jmao/work/public/releases/L1_Trigger/MC_Production/FollowCIEMAT/CMSSW_10_6_0_pre4/src/L1Trigger/phase2L1EcalTimingTP/test/testEvent.root')
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:10")


#process.source.fileNames = cms.untracked.vstring(
 #'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_pythia8/MINIAODSIM/PU200_pilot_103X_upgrade2023_realistic_v2_ext5-v1/40000/CE417376-6FCB-D74C-8FF9-3D0F214AF6B3.root')


# Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet
process = L1TrackTriggerTracklet(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
