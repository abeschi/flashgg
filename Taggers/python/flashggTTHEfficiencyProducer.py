import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections

bDiscriminator74X = cms.vdouble(0.605,0.890)
bDiscriminator76X = cms.vdouble(0.460,0.800,0.935)
bDiscriminator80XReReco = cms.vdouble(0.5426,0.8484,0.9535)

UnpackedJetCollectionVInputTag = cms.VInputTag()
for i in range(0,maxJetCollections):
    UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))

HTXSInputTags = cms.PSet(stage0cat = cms.InputTag("rivetProducerHTXS","stage0cat"),
                         stage1cat = cms.InputTag("rivetProducerHTXS","stage1cat"),
                         njets     = cms.InputTag("rivetProducerHTXS","njets"),
                         pTH       = cms.InputTag("rivetProducerHTXS","pTH"),
                         pTV       = cms.InputTag("rivetProducerHTXS","pTV"))



flashggTTHHadronicEfficiencyProducer = cms.EDProducer("FlashggTTHHadronicEfficiencyProducer",
                                       DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                       SystLabel=cms.string(""),
                                       MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( 'flashggPrunedGenParticles' ),  
                                       GenJetTag=cms.InputTag( 'slimmedGenJets' ),  
                                       MetTag=cms.InputTag( 'flashggMets' ),  
                                       HLTTag=cms.InputTag( 'TriggerResults', '', 'HLT' ),  
                                       rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
				       tthMVAweightfile = cms.FileInPath("flashgg/Taggers/data/TMVA_tth_hadronic_preMoriond_v0.weights.xml"),
                                       MVAMethod = cms.string("BDT"),
                                       #MVAMethod = cms.string("BDTG_ttH_vs_DiPho_Shepra_newSamplesTest_withCentralObjWeight") 
                                       inputTagJets= UnpackedJetCollectionVInputTag, 
                                       jetPtThreshold = cms.double(0.),
                                       jetEtaThreshold = cms.double(10),
                                       matchingGenPhotons = cms.double(0.1),
                                       bDiscriminator = bDiscriminator80XReReco, #bDiscriminator76X
                                       bTag = cms.string(flashggBTag),
				       useTTHHadronicMVA =  cms.bool(True),
                                       dRJetPhoLeadCut =  cms.double(0.4),
                                       dRJetPhoSubleadCut = cms.double(0.4),                          
                                       leptonPtThreshold = cms.double(20),
                                       muonEtaThreshold = cms.double(2.4), 
                                       muPFIsoSumRelThreshold = cms.double(0.25),
				       muMiniIsoSumRelThreshold = cms.double(0.06),	 
                                       electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                       useStdLeptonID = cms.bool(False),
                                       useElectronMVARecipe = cms.bool(False),
                                       useElectronLooseID = cms.bool(True),
                                       HTXSTags     = HTXSInputTags                                     
				       )

flashggTTHSemiLeptonicEfficiencyProducer = cms.EDProducer("FlashggTTHSemiLeptonicEfficiencyProducer",
                                       DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                       SystLabel=cms.string(""),
                                       MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( 'flashggPrunedGenParticles' ),  
                                       GenJetTag=cms.InputTag( 'slimmedGenJets' ),  
                                       MetTag=cms.InputTag( 'flashggMets' ),  
                                       HLTTag=cms.InputTag( 'TriggerResults', '', 'HLT' ),  
                                       rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                       inputTagJets= UnpackedJetCollectionVInputTag, 
                                       jetPtThreshold = cms.double(25),
                                       jetEtaThreshold = cms.double(2.4),
                                       matchingGenPhotons = cms.double(0.1),
                                       bDiscriminator = bDiscriminator80XReReco, #bDiscriminator76X
                                       bTag = cms.string(flashggBTag),
                                       dRJetPhoLeadCut =  cms.double(0.4),
                                       dRJetPhoSubleadCut = cms.double(0.4),                          
                                       leptonPtThreshold = cms.double(20),
                                       muonEtaThreshold = cms.double(2.4), 
                                       muPFIsoSumRelThreshold = cms.double(0.25),
				       muMiniIsoSumRelThreshold = cms.double(0.06),	 
                                       electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                       useStdLeptonID = cms.bool(False),
                                       useElectronMVARecipe = cms.bool(False),
                                       useElectronLooseID = cms.bool(True),
                                       HTXSTags     = HTXSInputTags                                     
				       )



flashggTTHFullyLeptonicEfficiencyProducer = cms.EDProducer("FlashggTTHFullyLeptonicEfficiencyProducer",
                                       DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                       SystLabel=cms.string(""),
                                       MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( 'flashggPrunedGenParticles' ),  
                                       GenJetTag=cms.InputTag( 'slimmedGenJets' ),  
                                       MetTag=cms.InputTag( 'flashggMets' ),  
                                       HLTTag=cms.InputTag( 'TriggerResults', '', 'HLT' ),  
                                       rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                       inputTagJets= UnpackedJetCollectionVInputTag, 
                                       jetPtThreshold = cms.double(25),
                                       jetEtaThreshold = cms.double(2.4),
                                       matchingGenPhotons = cms.double(0.1),
                                       bDiscriminator = bDiscriminator80XReReco, #bDiscriminator76X
                                       bTag = cms.string(flashggBTag),
                                       dRJetPhoLeadCut =  cms.double(0.4),
                                       dRJetPhoSubleadCut = cms.double(0.4),                          
                                       leptonPtThreshold = cms.double(20),
                                       muonEtaThreshold = cms.double(2.4), 
                                       muPFIsoSumRelThreshold = cms.double(0.25),
				       muMiniIsoSumRelThreshold = cms.double(0.06),	 
                                       electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                       useStdLeptonID = cms.bool(False),
                                       useElectronMVARecipe = cms.bool(False),
                                       useElectronLooseID = cms.bool(True),
                                       HTXSTags     = HTXSInputTags                                     
				       )





