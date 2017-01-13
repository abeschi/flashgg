
#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

# maxEvents is the max number of events processed of each file, not globally
inputFiles = "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/160707_152047/0000/myMicroAODOutputFile_1.root"
outputFile = "TTHHadronic.root" 

## I/O SETUP ##
process = cms.Process("TTHHadronicDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(inputFiles))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile))
from flashgg.MetaData.JobConfig import customize
customize.parse()

# load syst producer
from flashgg.Systematics.SystematicsCustomize import *
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")

useEGMTools(process)

## if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processId == 'Data' or customize.processId == 'data':
    print 'data' 
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
    customizePhotonSystematicsForMC(process)
    ##syst (2D) : smearings with EGMTool
    vpset2D   = process.flashggDiPhotonSystematics.SystMethods2D
    newvpset2D = cms.VPSet()
    for pset in vpset2D:
        pset.NSigmas = cms.PSet( firstVar = cms.vint32(), secondVar = cms.vint32() ) # only central value, no up/down syst shifts (2D case)
        if ( pset.Label.value().count("MCSmear") or pset.Label.value().count("SigmaEOverESmearing")):
            pset.ApplyCentralValue = cms.bool(True)
            newvpset2D+= [pset]
    process.flashggDiPhotonSystematics.SystMethods2D = newvpset2D



# flashgg tag sequence (for dipho MVA) and jet collections
from flashgg.Systematics.SystematicsCustomize import *
process.load("flashgg/Taggers/flashggTagSequence_cfi")
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
process.flashggTagSequence.remove(process.flashggUpdatedIdMVADiPhotons) # Needs to be run before systematics
massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggUpdatedIdMVADiPhotons"),cms.InputTag("flashggDiPhotonSystematics"))


# Require standard diphoton trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
    process.dataRequirements += process.hltHighLevel
    process.dataRequirements += process.eeBadScFilter
 
#process.load("flashgg.MicroAOD.flashggMets_cfi")
process.load("flashgg/Taggers/flashggTagTester_cfi")

from flashgg.Taggers.tagsDumpers_cfi import createTagDumper
import flashgg.Taggers.dumperConfigTools as cfgTools

process.TTHHadronicDumper = createTagDumper("TTHHadronicTag")
#process.load("flashgg.Taggers.tthDumper_cfi")
#process.flashggMuMuGamma.PhotonTag=cms.InputTag('flashggUpdatedIdMVAPhotons')
process.TTHHadronicDumper.dumpTrees = True
process.TTHHadronicDumper.dumpHistos = False
process.TTHHadronicDumper.dumpWorkspace = False

#from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
#massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggDiPhotons"),cms.InputTag("flashggPreselectedDiPhotons"))

import flashgg.Taggers.ttHTagVariables as var
all_variables = var.hadronic_variables + var.dipho_variables + var.truth_variables


cfgTools.addCategories(process.TTHHadronicDumper,
                       ## categories definition  
			[	("all","1",0)
			],
			variables = all_variables,
			histograms = []
                     )


process.TTHHadronicDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"

from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents" ,-1)    # max-number of events
customize.setDefault("targetLumi",1e+3) # define integrated lumi
customize(process)


process.p1 = cms.Path(
    process.dataRequirements*
    process.flashggUpdatedIdMVADiPhotons*
    process.flashggDiPhotonSystematics*
    process.flashggTagSequence*
    process.TTHHadronicDumper
    )

print process.p1
