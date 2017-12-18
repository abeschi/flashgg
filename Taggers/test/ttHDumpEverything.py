#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet
import os

# maxEvents is the max number of events processed of each file, not globally
inputFiles = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/ttHToGG_M125_13TeV_powheg_pythia8_v2/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_093929/0000/myMicroAODOutputFile_1.root"
#inputFiles = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_095729/0000/myMicroAODOutputFile_10.root"
#inputFiles = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170121_195858/0000/myMicroAODOutputFile_99.root"
#inputFiles = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170121_195858/0000/myMicroAODOutputFile_1.root"
#inputFiles = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/ReMiniAOD-03Feb2017-2_5_4/2_5_1/DoubleEG/ReMiniAOD-03Feb2017-2_5_4-2_5_1-v0-Run2016C-03Feb2017-v1/170310_111520/0000/myMicroAODOutputFile_7.root"
outputFile = "output.root" 

dropVBFInNonGold = False


## I/O SETUP ##
process = cms.Process("DumpEverything")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif os.environ["CMSSW_VERSION"].count("CMSSW_7_4"):
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' 
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
else:
    raise Exception,"Could not find a sensible CMSSW_VERSION for default globaltag"
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1)

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(inputFiles))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile))
from flashgg.MetaData.JobConfig import customize
customize.parse()


# SYSTEMATICS SECTION
from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process)
if dropVBFInNonGold:
    process.flashggVBFTag.SetArbitraryNonGoldMC = True
    process.flashggVBFTag.DropNonGoldData = True
modifyTagSequenceForSystematics(process, jetSystematicsInputTags)


systlabels = [""]
phosystlabels = []
jetsystlabels = []
metsystlabels = []
elesystlabels = []
musystlabels = []



# Or use the official tool instead
useEGMTools(process)

if customize.processId == "Data":
    customizeSystematicsForData(process)
else:
    customizeSystematicsForBackground(process) # only central corrections, no syst. shifts

cloneTagSequenceForEachSystematic(process, systlabels, phosystlabels, metsystlabels, jetsystlabels, jetSystematicsInputTags)


print 'syst 1D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods])
print 'syst 2D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods2D])


process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
    process.dataRequirements += process.hltHighLevel
    process.dataRequirements += process.eeBadScFilter
 

#import flashgg.Taggers.ttHEfficiencyVariables as var
#if customize.processId == "Data":
#	variables = ''
#else:
#	variables = var.dipho_variables + var.hadronic_variables + var.truth_photon_variables + var.efficiency_variables




process.load("flashgg.Taggers.TTHEfficiencyDumper_cfi")
import flashgg.Taggers.dumperConfigTools as cfgTools


process.TTHHadronicEfficiencyDumper.dumpTrees = True
process.TTHHadronicEfficiencyDumper.dumpWorkspace = False
process.TTHHadronicEfficiencyDumper.quietRooFit = True

import flashgg.Taggers.ttHEfficiencyVariables as var
MyVarsH =  var.truth_photon_variables + var.efficiency_hadronic_variables + var.dipho_variables + var.hadronic_variables
MyVarsL =  var.truth_photon_variables + var.efficiency_leptonic_variables + var.dipho_variables + var.leptonic_variables + var.leptonic_variables_all


# split tree, histogram and datasets by process
process.TTHHadronicEfficiencyDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"

cfgTools.addCategories(process.TTHHadronicEfficiencyDumper,
                       [("all","1>0",0),
                        ],
                        variables = MyVarsH,
                       histograms=[]
                       )

process.TTHSemiLeptonicEfficiencyDumper.dumpTrees = True
process.TTHSemiLeptonicEfficiencyDumper.dumpWorkspace = False
process.TTHSemiLeptonicEfficiencyDumper.quietRooFit = True



# split tree, histogram and datasets by process
process.TTHSemiLeptonicEfficiencyDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"

cfgTools.addCategories(process.TTHSemiLeptonicEfficiencyDumper,
                       [("all","1>0",0),
                        ],
                        variables = MyVarsL,
                       histograms=[]
                       )


process.TTHFullyLeptonicEfficiencyDumper.dumpTrees = True
process.TTHFullyLeptonicEfficiencyDumper.dumpWorkspace = False
process.TTHFullyLeptonicEfficiencyDumper.quietRooFit = True


# split tree, histogram and datasets by process
process.TTHFullyLeptonicEfficiencyDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"

cfgTools.addCategories(process.TTHFullyLeptonicEfficiencyDumper,
                       [("all","1>0",0),
                        ],
                        variables = MyVarsL,
                       histograms=[]
                       )

#process.flashggTagSequence.remove(flashggUpdatedIdMVADiPhotons)


process.p1 = cms.Path(process.dataRequirements*
                     process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
		     process.flashggMetSystematics*
                     process.flashggMuonSystematics*process.flashggElectronSystematics*
                     (process.flashggUnpackedJets*process.jetSystematicsSequence)*
                     process.flashggTagSequence*
                     process.TTHHadronicEfficiencyDumper*
                     process.TTHSemiLeptonicEfficiencyDumper*
                     process.TTHFullyLeptonicEfficiencyDumper)
 

#printSystematicInfo(process)

## set default options if needed
customize.setDefault("maxEvents", -1)
customize.setDefault("targetLumi",1e+3)
## call the customization
customize(process)


print process.p1 

