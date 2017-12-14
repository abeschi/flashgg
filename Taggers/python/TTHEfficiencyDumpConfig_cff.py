#-----------J. Tao from IHEP-Beijing--------------
import FWCore.ParameterSet.Config as cms

from globalVariables_cff import globalVariables

TTHHadronicEfficiencyDumpConfig = cms.PSet(
    className  = cms.untracked.string("CutBasedTTHHadronicEfficiencyDumper"),
    src = cms.InputTag("flashggTTHHadronicEfficiencyProducer"),
    #src = cms.InputTag(""),
    generatorInfo = cms.InputTag("generator"),
    processId = cms.string("tth"),
    maxCandPerEvent = cms.int32(1),
    lumiWeight = cms.double(1.0),
    classifierCfg = cms.PSet(categories=cms.VPSet()),
    categories = cms.VPSet(),

    workspaceName = cms.untracked.string("cms_hgg_$SQRTS"),
    nameTemplate = cms.untracked.string("$PROCESS_$SQRTS_$LABEL_$SUBCAT"),
    
    dumpHistos = cms.untracked.bool(True),
    dumpWorkspace = cms.untracked.bool(False),
    dumpTrees = cms.untracked.bool(False),
    
    quietRooFit = cms.untracked.bool(False),
    dumpGlobalVariables = cms.untracked.bool(True),
    globalVariables = globalVariables

)


TTHSemiLeptonicEfficiencyDumpConfig = cms.PSet(
    className  = cms.untracked.string("CutBasedTTHSemiLeptonicEfficiencyDumper"),
    src = cms.InputTag("flashggTTHSemiLeptonicEfficiencyProducer"),
    #src = cms.InputTag(""),
    generatorInfo = cms.InputTag("generator"),
    processId = cms.string("tth"),
    maxCandPerEvent = cms.int32(1),
    lumiWeight = cms.double(1.0),
    classifierCfg = cms.PSet(categories=cms.VPSet()),
    categories = cms.VPSet(),

    workspaceName = cms.untracked.string("cms_hgg_$SQRTS"),
    nameTemplate = cms.untracked.string("$PROCESS_$SQRTS_$LABEL_$SUBCAT"),
    
    dumpHistos = cms.untracked.bool(True),
    dumpWorkspace = cms.untracked.bool(False),
    dumpTrees = cms.untracked.bool(False),
    
    quietRooFit = cms.untracked.bool(False),
    dumpGlobalVariables = cms.untracked.bool(True),
    globalVariables = globalVariables

)


TTHFullyLeptonicEfficiencyDumpConfig = cms.PSet(
    className  = cms.untracked.string("CutBasedTTHFullyLeptonicEfficiencyDumper"),
    src = cms.InputTag("flashggTTHFullyLeptonicEfficiencyProducer"),
    #src = cms.InputTag(""),
    generatorInfo = cms.InputTag("generator"),
    processId = cms.string("tth"),
    maxCandPerEvent = cms.int32(1),
    lumiWeight = cms.double(1.0),
    classifierCfg = cms.PSet(categories=cms.VPSet()),
    categories = cms.VPSet(),

    workspaceName = cms.untracked.string("cms_hgg_$SQRTS"),
    nameTemplate = cms.untracked.string("$PROCESS_$SQRTS_$LABEL_$SUBCAT"),
    
    dumpHistos = cms.untracked.bool(True),
    dumpWorkspace = cms.untracked.bool(False),
    dumpTrees = cms.untracked.bool(False),
    
    quietRooFit = cms.untracked.bool(False),
    dumpGlobalVariables = cms.untracked.bool(True),
    globalVariables = globalVariables

)
