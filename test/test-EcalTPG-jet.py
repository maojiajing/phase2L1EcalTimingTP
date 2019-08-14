# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: repr --processName=REPR --python_filename=reprocess_test_10_5_0_pre1.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 2 --era Phase2 --eventcontent FEVTDEBUGHLT --filein root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step2_2ev_reprocess_slim.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('REPR2',eras.Phase2C4_trigger)
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
    input = cms.untracked.int32(2)
    #input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# Neutrino samples
# "/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/NeutrinoGun_E_10GeV/MINIAODSIM/PU200_103X_upgrade2023_realistic_v2-v1/280000/260E6EA0-EDDF-5A46-8058-7D14AD4D103A.root"),
                               #'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_TuneCP5_pythia8/MINIAODSIM/NoPU_103X_upgrade2023_realistic_v2-v2/90000/F063E7F2-A71C-E14C-BF02-665560C6B3CC.root'),
#				  'file:/afs/cern.ch/user/j/jmao/work/public/releases/L1_Trigger/MC_Production/FollowCIEMAT/CMSSW_10_6_0_pre4/src/L1Trigger/phase2L1EcalTimingTP/test/testEvent.root'),
				#'/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_103X_upgrade2023_realistic_v2-v1/80000/EC78B837-81E3-054E-8A11-05762CDA1980.root'),
#				  'file:/afs/cern.ch/user/j/jmao/work/public/releases/L1_Trigger/MC_Production/FollowCIEMAT/CMSSW_10_6_0_pre4/src/L1Trigger/phase2L1EcalTimingTP/test/testEvent.root'),
				  'file:/afs/cern.ch/user/j/jmao/work/public/releases/L1_Trigger/MC_Production/FollowCIEMAT/CMSSW_10_6_1_patch2/src/step2_2ev_reprocess_slim_jet.root'),

   secondaryFileNames = cms.untracked.vstring(
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/9CDB6CBE-E063-9049-A794-244A09CA53A0.root'
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/354CC7AD-BD24-7249-A26A-93F07D8194FD.root',
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/FF015593-7F87-9142-8363-BCEAAA2B0620.root',
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/E862C343-6377-B74E-B70A-1B7E7EABCB7C.root',
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/E2914C42-0997-7546-93AF-157315C33BFD.root',
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/30B662A6-E4FD-644A-B124-29107C17ECA1.root',
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/E5DA8801-DFA7-9946-8B72-E3F884890188.root',
#       'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/8E9FB60D-C574-8645-916D-F8211D38968D.root',
)
# Neutrino sample
#'root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/NeutrinoGun_E_10GeV/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/280000/CF7695FC-46FE-F649-9A7E-47ABE65D3861.root'),
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:3445-1:3455",
#                                                                    "1:1849-1:1861", "1:2449-1:2462",
#                                                                    "1:1862-1:1872", "1:2463-1:2472", "1:3639-1:3644", "1:3648",
#                                                                    "1:3433-1:3444", "1:3456", "1:3625:3638"
#                                                                    )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('repr nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

#process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
#        filterName = cms.untracked.string('')
#    ),
#    fileName = cms.untracked.string('file:step2_2ev_reprocess_slim.root'),
#    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
#    splitLevel = cms.untracked.int32(0)
#)

# Additional output definition


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

# Path and EndPath definitions
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

############################################################
# L1 Tau object
############################################################

#process.load("L1Trigger.Phase2L1Taus.L1PFTauProducer_cff")
#process.L1PFTauProducer.L1PFObjects = cms.InputTag("l1pfCandidates","PF")
#process.L1PFTauProducer.L1Neutrals = cms.InputTag("l1pfCandidates")
#process.L1PFTauProducer.L1Clusters = cms.InputTag("l1pfCandidates","PF")
#process.L1PFTauProducer.min_pi0pt = cms.double(3)
#process.L1PFTaus = cms.Path(process.L1PFTauProducer)

# L1 EcalTiming Analyzer
process.load("L1Trigger.phase2L1EcalTimingTP.phase2L1EcalTimingAnalyzer_cfi")

process.L1EcalTimingAnalyzer = cms.EDAnalyzer('phase2L1EcalTimingAnalyzer',
                                       ak4GenJets = cms.InputTag("ak4GenJets", "", "HLT"),
                                       ak4GenJetsNoNu = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
                                       ak8GenJets = cms.InputTag("ak8GenJets", "", "HLT"),
                                       ak8GenJetsNoNu = cms.InputTag("ak8GenJetsNoNu", "", "HLT"),
                                       genParticles = cms.InputTag("genParticles", "", "HLT"),
                                       genParticles_t0 = cms.InputTag("genParticles", "t0", "HLT"),
                                       ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","REPR"),
                                       )

process.analyzer = cms.Path(process.L1EcalTimingAnalyzer)

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer-dyll-4FEVT-jet.root"), 
   closeFileFast = cms.untracked.bool(True)
)


# Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
#process.schedule = cms.Schedule(process.L1simulation_step,process.L1PFTaus,process.analyzer,process.endjob_step)
process.schedule = cms.Schedule(process.analyzer,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

#from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet
#process = L1TrackTriggerTracklet(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
