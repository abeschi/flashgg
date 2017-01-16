#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
          "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_142720/0000/myMicroAODOutputFile_1.root"
        )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test_mmg.root")
)

from flashgg.MetaData.JobConfig import customize
customize.parse()

print customize.processId
if customize.processId == 'data':
    print 'DATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' 
else :
    print 'MCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' 

#update photon shower shapes and photon id
process.load("flashgg.Taggers.flashggPhotonWithUpdatedIdMVAProducer_cfi")

#preselection on single photon???
#process.load("flashgg.Taggers.flashggPreselectedDiPhotons")
#process.flashggPreselectedDiPhotons.src = "flashggUpdatedIdMVAPhotons"

#re-run mumugamma producer
process.load("flashgg.MicroAOD.flashggMuMuGamma_cfi")
process.flashggMuMuGamma.PhotonTag=cms.InputTag('flashggUpdatedIdMVAPhotons')

process.load("flashgg.Taggers.mumugammaDumper_cfi") ##  import mumugammaDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools

process.mumugammaDumper.dumpTrees = True
process.mumugammaDumper.dumpWorkspace = False
process.mumugammaDumper.quietRooFit = True


# split tree, histogram and datasets by process
process.mumugammaDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"

## define categories and associated objects to dump
cfgTools.addCategory(process.mumugammaDumper,
                     "Reject",
                      " !Is2012FSRZMMG ",
                       -1 ## if nSubcat is -1 do not store anythings
                     )

# interestng categories 
cfgTools.addCategories(process.mumugammaDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [("EB","abs(MMG_Photon.superCluster.eta)<1.5",0), ##
                        ("EE","abs(MMG_Photon.superCluster.eta)>1.5",0),##("EE","1",0), ## evereything elese is EB+EE
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

process.dataRequirements = cms.Sequence()

if customize.processId == 'data':
    from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
            #DoubleEG
            "HLT_IsoMu15_v8",
            "HLT_DoubleMu7_v3"
            #SingleEG
            "HLT_Mu24_v3",
            #"HLT_Ele27_WPLoose_Gsf_v*" # 7_6_X
            ##"HLT_Ele27_WPLoose_Gsf_v*",
            ##"HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
            ##"HLT_Ele22_eta2p1_WPLoose_Gsf_v*"
            ) )

if customize.processId == 'data':

	process.p1 = cms.Path(
   		process.hltHighLevel*
    	process.dataRequirements*
    	process.flashggUpdatedIdMVAPhotons*
    	process.flashggMuMuGamma*
    	process.mumugammaDumper
    	)

else :
	process.p1 = cms.Path(
    	process.dataRequirements*
    	process.flashggUpdatedIdMVAPhotons*
    	process.flashggMuMuGamma*
    	process.mumugammaDumper
    	)



customize.setDefault("maxEvents",10000)
customize(process)
