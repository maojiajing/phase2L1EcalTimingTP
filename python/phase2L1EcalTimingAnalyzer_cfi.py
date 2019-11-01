import FWCore.ParameterSet.Config as cms

L1EcalTimingAnalyzer = cms.EDAnalyzer('phase2L1EcalTimingAnalyzer',
                               genParticles     = cms.InputTag("genParticles", "", "HLT"),
                               genParticles_t0     = cms.InputTag("genParticles", "t0", "HLT"),
                               ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","REPR"),
)
 
