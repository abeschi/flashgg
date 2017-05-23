
#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

# maxEvents is the max number of events processed of each file, not globally
inputFiles = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/170212_190415/0000/myMicroAODOutputFile_293.root"
outputFile = "MuMuGamma.root" 

## I/O SETUP ##
process = cms.Process("mumugammaDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(inputFiles))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile))

from flashgg.MetaData.JobConfig import customize
customize.parse()

process.load("flashgg.Taggers.flashggPhotonWithUpdatedIdMVAProducer_cfi")


# load syst producer
from flashgg.Systematics.SystematicsCustomize import *
process.load("flashgg.Systematics.PhotonSystematics_cfi")



for isyst in [ process.MCScaleHighR9EB, process.MCScaleLowR9EB, process.MCScaleHighR9EE, process.MCScaleLowR9EE ]:
    process.flashggPhotonSystematics.SystMethods.remove(isyst)

    # add EGM scales
for isyst in [ process.MCScaleHighR9EB_EGM, process.MCScaleLowR9EB_EGM, process.MCScaleHighR9EE_EGM, process.MCScaleLowR9EE_EGM ]:
    process.flashggPhotonSystematics.SystMethods.insert(0, isyst)

    # remove old smearings
for isyst in [ process.MCSmearHighR9EE, process.MCSmearLowR9EE, process.MCSmearHighR9EB, process.MCSmearLowR9EB, process.SigmaEOverESmearing ]:
    process.flashggPhotonSystematics.SystMethods.remove(isyst)

    # add EGM smearings (2D)
process.flashggPhotonSystematics.SystMethods2D.extend([
       process.MCSmearHighR9EE_EGM,
       process.MCSmearLowR9EE_EGM,
       process.MCSmearHighR9EB_EGM,
       process.MCSmearLowR9EB_EGM,
       process.SigmaEOverESmearing_EGM])

## if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processId == 'Data' or customize.processId == 'data':
    print 'data' 
    photonScaleBinsData = getattr(process,'photonScaleBinsData',None)
    if hasattr(process,'photonScaleBinsData'):
        print photonScaleBinsData, process.photonScaleBinsData
    process.flashggPhotonSystematics.SystMethods = customizeVPSetForData(process.flashggPhotonSystematics.SystMethods, photonScaleBinsData)
    process.flashggPhotonSystematics.SystMethods2D = customizeVPSetForData(process.flashggPhotonSystematics.SystMethods2D, photonScaleBinsData)
else:
    print 'mc'

    ##syst (2D) : smearings with EGMTool

    photonSmearBins = getattr(process,'photonSmearBins',None)
    photonScaleUncertBins = getattr(process,'photonScaleUncertBins',None)
    for pset in process.flashggPhotonSystematics.SystMethods:
        if photonSmearBins and pset.Label.value().startswith("MCSmear"):
            pset.BinList = photonSmearBins
        elif photonScaleUncertBins and pset.Label.value().count("Scale"):
            pset.BinList = photonScaleUncertBins

    vpset2D   = process.flashggPhotonSystematics.SystMethods2D
    newvpset2D = cms.VPSet()
    for pset in vpset2D:
        pset.NSigmas = cms.PSet( firstVar = cms.vint32(), secondVar = cms.vint32() ) # only central value, no up/down syst shifts (2D case)
        if ( pset.Label.value().count("MCSmear") or pset.Label.value().count("SigmaEOverESmearing")):
            pset.ApplyCentralValue = cms.bool(True)
            newvpset2D+= [pset]
    process.flashggPhotonSystematics.SystMethods2D = newvpset2D


#re-run mumugamma producer
process.load("flashgg.MicroAOD.flashggMuMuGamma_cfi")
process.flashggMuMuGamma.PhotonTag=cms.InputTag('flashggUpdatedIdMVAPhotons')

process.load("flashgg.Taggers.mumugammaDumper_cfi") ##  import mumugammaDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools








# Require trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_DoubleMu*"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
#process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data" or customize.processId == "data":
    process.dataRequirements += process.hltHighLevel
   # process.dataRequirements += process.eeBadScFilter
 





process.load("flashgg.Taggers.mumugammaDumper_cfi") ##  import mumugammaDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools



process.mumugammaDumper.dumpTrees = True
process.mumugammaDumper.dumpWorkspace = False
process.mumugammaDumper.quietRooFit = True


# split tree, histogram and datasets by process
process.mumugammaDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"

cfgTools.addCategory(process.mumugammaDumper,
                     "Reject",
                      " !Is2012FSRZMMG ",
                       -1 ## if nSubcat is -1 do not store anythings
                     )


cfgTools.addCategories(process.mumugammaDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [("all","1>0",0),
                        ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       ## if different variables wanted for different categories, can add categorie one by one with cfgTools.addCategory
                       variables=["mass_mmg            :=mass",
                                  "mass_mumu           :=DiMuPtr.mass", 
                                  "pt_mumu             :=DiMuPtr.pt", 
                                  "leadMuonP           :=DiMuPtr.leadingMuon.p",
                                  "leadMuonPt          :=DiMuPtr.leadingMuon.pt",
                                  "leadMuonEta         :=DiMuPtr.leadingMuon.eta",
                                  "leadMuonPhi         :=DiMuPtr.leadingMuon.phi",
                                  "leadMuonCharge      :=DiMuPtr.leadingMuon.charge",
                                  "leadMuonNtrk        :=DiMuPtr.leadingMuon.innerTrack().hitPattern().trackerLayersWithMeasurement()", ## for rochester corrections
                                  "subleadMuonP        :=DiMuPtr.subleadingMuon.p",
                                  "subleadMuonPt       :=DiMuPtr.subleadingMuon.pt",
                                  "subleadMuonEta      :=DiMuPtr.subleadingMuon.eta",
                                  "subleadMuonPhi      :=DiMuPtr.subleadingMuon.phi",
                                  "subleadMuonCharge   :=DiMuPtr.subleadingMuon.charge",
                                  "subleadMuonNtrk     :=DiMuPtr.subleadingMuon.innerTrack().hitPattern().trackerLayersWithMeasurement()", ## for rochester corrections
                                  "photonE             :=MMG_Photon.energy",
                                  "photonPt            :=MMG_Photon.pt",
                                  "photonScEta         :=MMG_Photon.superCluster.eta",
                                  "photonScPhi         :=MMG_Photon.superCluster.phi",
                                  "photonEta           :=MMG_Photon.eta",
                                  "photonPhi           :=MMG_Photon.phi",
                                  "photonR9            :=MMG_Photon.full5x5_r9",
                                  "photonS4            :=MMG_Photon.s4",
                                  "photonEtaWidth      :=MMG_Photon.superCluster().etaWidth()"
                                  ],
                       ## histograms to be plotted. 
                       ## the variables need to be defined first
                       histograms=[
        ]
                       )


process.mumugammaDumper.nameTemplate = "tree"


customize.setDefault("maxEvents" , 100)    # max-number of events
customize.setDefault("targetLumi",1e+3) # define integrated lumi
customize(process)


process.p1 = cms.Path(
    	process.dataRequirements*
    	process.flashggUpdatedIdMVAPhotons*
	process.flashggPhotonSystematics*
    	process.flashggMuMuGamma*
    	process.mumugammaDumper
    	)

print process.p1
