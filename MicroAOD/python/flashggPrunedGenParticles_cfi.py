import FWCore.ParameterSet.Config as cms

flashggPrunedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
        "drop  *  ",              # this is the default
        "keep++ abs(pdgId) = 6",  #save t descendants
        "keep++ pdgId = 23",      #save Z descendants
        "keep++ abs(pdgId) = 24", #save W descendants
        "keep++ pdgId = 25"       #save H descendants
        "keep status = 3",
        "keep status = 22",
        "keep status = 23",
        "drop abs(pdgId) > 25"
        )
    )
