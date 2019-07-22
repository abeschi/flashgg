#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet
from flashgg.MetaData.MetaConditionsReader import *
import os

# maxEvents is the max number of events processed of each file, not globally
inputFiles = "file:/afs/cern.ch/user/m/malberti/public/xABeschi/Zgamma/myMicroAODOutputFile_ggh.root"
outputFile = "output.root" 



## I/O SETUP ##
process = cms.Process("ZGammaDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(inputFiles))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile))

from flashgg.MetaData.JobConfig import customize
customize.parse()
customize.metaConditions = MetaConditionsReader(customize.metaConditions)
if customize.processId == "Data":
    process.GlobalTag.globaltag = str(customize.metaConditions['globalTags']['data'])
else:
    process.GlobalTag.globaltag = str(customize.metaConditions['globalTags']['MC'])

# SYSTEMATICS SECTION
'''
from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process)
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

systlabels = [""]
phosystlabels = []
jetsystlabels = []
metsystlabels = []
elesystlabels = []
musystlabels = []

# Or use the official tool instead
useEGMTools(process)

#customize.processId = "Data"

if customize.processId == "Data":
    customizeSystematicsForData(process)
else:
    customizeSystematicsForBackground(process) # only central corrections, no syst. shifts

cloneTagSequenceForEachSystematic(process, systlabels, phosystlabels, metsystlabels, jetsystlabels, jetSystematicsInputTags)


print 'syst 1D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods])
print 'syst 2D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods2D])
'''

'''
# Require standard diphoton trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ))
'''
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
    process.dataRequirements += process.hltHighLevel
    process.dataRequirements += process.eeBadScFilter

from flashgg.Taggers.tagsDumpers_cfi import createTagDumper
import flashgg.Taggers.dumperConfigTools as cfgTools

process.ZGammaToEEUntaggedTagDumper = createTagDumper("ZGammaToEEUntaggedTag")

process.ZGammaToEEUntaggedTagDumper.dumpTrees = True
process.ZGammaToEEUntaggedTagDumper.dumpHistos = False
process.ZGammaToEEUntaggedTagDumper.dumpWorkspace = False

import flashgg.Taggers.ZGammaVariables as var
variablesEle_ = var.variables_ele



cfgTools.addCategories(process.ZGammaToEEUntaggedTagDumper,
                       ## categories definition  
			[	("all","1",0)
			],
			variables = variablesEle_,
			histograms = []
                     )

process.ZGammaToMuMuUntaggedTagDumper = createTagDumper("ZGammaToMuMuUntaggedTag")

process.ZGammaToMuMuUntaggedTagDumper.dumpTrees = True
process.ZGammaToMuMuUntaggedTagDumper.dumpHistos = False
process.ZGammaToMuMuUntaggedTagDumper.dumpWorkspace = False

import flashgg.Taggers.ZGammaVariables as var
variablesMu_ = var.variables_mu



cfgTools.addCategories(process.ZGammaToMuMuUntaggedTagDumper,
                       ## categories definition  
			[	("all","1",0)
			],
			variables = variablesMu_,
			histograms = []
                     )

process.load("flashgg/Taggers/flashggZGammaTagSequence_cfi")

#process.flashggZGammaToMuMuUntaggedTag.debug = cms.bool(True)

customize.setDefault("maxEvents" , 1000)    # max-number of events
customize.setDefault("targetLumi",1e+3) # define integrated lumi
customize(process)


process.p1 = cms.Path(process.dataRequirements*
                     process.flashggZGammaTagSequence*
		     process.ZGammaToMuMuUntaggedTagDumper*
                     process.ZGammaToEEUntaggedTagDumper)

print process.p1









