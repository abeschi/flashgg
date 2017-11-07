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
                                       leadPhoOverMassThreshold = cms.double(0.),
                                       leadPhoPtThreshold = cms.double(0),  
                                       leadPhoUseVariableThreshold =  cms.bool(True),
                                       subleadPhoOverMassThreshold = cms.double(0.),
                                       subleadPhoPtThreshold = cms.double(0),
                                       subleadPhoUseVariableThreshold =  cms.bool(True),
                                       MVAThreshold = cms.double(0.),
                                       PhoMVAThreshold = cms.double(-2),
                                       inputTagJets= UnpackedJetCollectionVInputTag, 
                                       jetPtThreshold = cms.double(0.),
                                       jetEtaThreshold = cms.double(10),
                                       bDiscriminator = bDiscriminator80XReReco, #bDiscriminator76X
                                       bTag = cms.string(flashggBTag),
                                       jetsNumberThreshold = cms.int32(0),
                                       bjetsNumberThreshold = cms.int32(0),
				       bjetsLooseNumberThreshold = cms.int32(0),
				       useTTHHadronicMVA =  cms.bool(True),
				       leadPhoOverMassTTHHMVAThreshold = cms.double(0.),
				       MVATTHHMVAThreshold = cms.double(0.),
				       jetsNumberTTHHMVAThreshold = cms.int32(0),
                                       bjetsNumberTTHHMVAThreshold = cms.int32(0),
                                       bjetsLooseNumberTTHHMVAThreshold = cms.int32(0),  
				       tthHadMVAThreshold = cms.double(0.),
                                       dRJetPhoLeadCut =  cms.double(0.),
                                       dRJetPhoSubleadCut = cms.double(0.),                          
                                       leptonPtThreshold = cms.double(20),
                                       muonEtaThreshold = cms.double(2.4), 
                                       muPFIsoSumRelThreshold = cms.double(0.25),
				       muMiniIsoSumRelThreshold = cms.double(0.06),	 
                                       electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                       nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                       nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                       electronIsoThreshold = cms.double(0.15),
				       elMiniIsoEBThreshold = cms.double(0.045),
                                       elMiniIsoEEThreshold = cms.double(0.08),
                                       electronNumOfHitsThreshold = cms.double(1),
                                       TransverseImpactParam = cms.double(0.02),
                                       LongitudinalImpactParam = cms.double(0.2),
				       TransverseImpactParamEB = cms.double(0.0261),
                                       LongitudinalImpactParamEB = cms.double(0.41),
                                       TransverseImpactParamEE = cms.double(0.118),
                                       LongitudinalImpactParamEE = cms.double(0.822),
                                       useStdLeptonID = cms.bool(False),
                                       useElectronMVARecipe = cms.bool(False),
                                       useElectronLooseID = cms.bool(True),
                                       HTXSTags     = HTXSInputTags                                     
				       )
