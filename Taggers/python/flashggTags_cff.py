import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, flashggDeepCSV, maxJetCollections

bDiscriminator74X = cms.vdouble(0.605,0.890)
bDiscriminator76X = cms.vdouble(0.460,0.800,0.935)
bDiscriminator80XReReco = cms.vdouble(0.5426,0.8484,0.9535)
bDiscriminator94X= cms.vdouble(0.1522,0.4941,0.8001)

flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
                                     JetsTag = cms.InputTag("flashggFinalJets"),
                                     NCollections = cms.uint32(maxJetCollections)
                                     )

UnpackedJetCollectionVInputTag = cms.VInputTag()
for i in range(0,maxJetCollections):
    UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))

HTXSInputTags = cms.PSet(stage0cat = cms.InputTag("rivetProducerHTXS","stage0cat"),
                         stage1cat = cms.InputTag("rivetProducerHTXS","stage1cat"),
                         njets     = cms.InputTag("rivetProducerHTXS","njets"),
                         pTH       = cms.InputTag("rivetProducerHTXS","pTH"),
                         pTV       = cms.InputTag("rivetProducerHTXS","pTV"))

flashggUntagged = cms.EDProducer("FlashggUntaggedTagProducer",
#                                 DiPhotonTag=cms.InputTag('flashggDiPhotons'),
                                 DiPhotonTag    = cms.InputTag('flashggPreselectedDiPhotons'),
                                 SystLabel      = cms.string(""),
                                 MVAResultTag   = cms.InputTag('flashggDiPhotonMVA'),
                                 GenParticleTag = cms.InputTag( "flashggPrunedGenParticles" ),
                                 Boundaries     = cms.vdouble(-0.405,0.204,0.564,0.864), #,1.000),
                                 RequireScaledPtCuts = cms.bool(True),
                                 HTXSTags     = HTXSInputTags
)

flashggSigmaMoMpToMTag = cms.EDProducer("FlashggSigmaMpTTagProducer",
#                                 DiPhotonTag=cms.InputTag('flashggDiPhotons'),
                                 DiPhotonTag    = cms.InputTag('flashggPreselectedDiPhotons'),
                                 SystLabel      = cms.string(""),
                                 MVAResultTag   = cms.InputTag('flashggDiPhotonMVA'),
                                 GenParticleTag = cms.InputTag( "flashggPrunedGenParticles" ),
                                 BoundariesSigmaMoM  = cms.vdouble(0.,0.00764,0.0109,0.0288), #boundaries have to be provided including lowest and highest
#                                 BoundariespToM      = cms.vdouble(0.,1.02,1.83,10.0), #,1.000), #boundaries have to be provided including lowest and highest
                                 RequireScaledPtCuts = cms.bool(True)
)



flashggTTHHadronicTag = cms.EDProducer("FlashggTTHHadronicTagProducer",
                                       DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                       SystLabel=cms.string(""),
                                       MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( 'flashggPrunedGenParticles' ),  
                                       rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
				       tthMVAweightfile = cms.FileInPath("flashgg/Taggers/data/TMVA_tth_hadronic_preMoriond_v0.weights.xml"),
                                       MVAMethod = cms.string("BDT"),
                                       #MVAMethod = cms.string("BDTG_ttH_vs_DiPho_Shepra_newSamplesTest_withCentralObjWeight") 
                                       leadPhoOverMassThreshold = cms.double(0.5),
                                       leadPhoPtThreshold = cms.double(20),  
                                       leadPhoUseVariableThreshold =  cms.bool(True),
                                       subleadPhoOverMassThreshold = cms.double(0.25),
                                       subleadPhoPtThreshold = cms.double(20),
                                       subleadPhoUseVariableThreshold =  cms.bool(True),
                                       MVAThreshold = cms.double(0.6),
                                       PhoMVAThreshold = cms.double(-0.9),
                                       inputTagJets= UnpackedJetCollectionVInputTag, 
                                       jetPtThreshold = cms.double(25.),
                                       jetEtaThreshold = cms.double(2.4),
#                                       bDiscriminator = bDiscriminator80XReReco, #bDiscriminator76X
#                                       bTag = cms.string(flashggBTag),
                                       bDiscriminator = bDiscriminator94X,
                                       bTag = cms.string(flashggDeepCSV),
                                       jetsNumberThreshold = cms.int32(5),
                                       bjetsNumberThreshold = cms.int32(1),
				       bjetsLooseNumberThreshold = cms.int32(0),
				       useTTHHadronicMVA =  cms.bool(True),
				       leadPhoOverMassTTHHMVAThreshold = cms.double(0.33333),
				       MVATTHHMVAThreshold = cms.double(0.4),
				       jetsNumberTTHHMVAThreshold = cms.int32(3),
                                       bjetsNumberTTHHMVAThreshold = cms.int32(0),
                                       bjetsLooseNumberTTHHMVAThreshold = cms.int32(1),  
				       tthHadMVAThreshold = cms.double(0.75),
                                       dRJetPhoLeadCut =  cms.double(0.4),
                                       dRJetPhoSubleadCut = cms.double(0.4),                          
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

flashggVBFTag = cms.EDProducer("FlashggVBFTagProducer",
                               DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                               SystLabel=cms.string(""),
                               MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                               VBFDiPhoDiJetMVAResultTag=cms.InputTag('flashggVBFDiPhoDiJetMVA'),
                               VBFMVAResultTag=cms.InputTag('flashggVBFMVA'),
                               GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
                               GenJetTag = cms.InputTag("slimmedGenJets"),
                               #Boundaries=cms.vdouble(0.21,0.6,0.81)
                               #  for the moment we have two categories VBF-0 and VBF-1: to be changed when the diphoton MVA is ready 
                               #Boundaries=cms.vdouble(0.5819, 0.9449)
                               #Boundaries=cms.vdouble(0.62, 0.94),
                               #Boundaries=cms.vdouble(0.634, 0.919),

                               # optimisation of Moriond 17 : 3 VBF categories
                               # These boundaries are recalculated after fixing
                               # the problem with the shape of the BDT output
#                               Boundaries=cms.vdouble(0.215,  0.532,  0.865),
#                               Boundaries=cms.vdouble(0.66633615,  0.89334188,  0.95919197),
#                               Boundaries=cms.vdouble(0.55889473,  0.9087378 ,  0.97044208),
                               Boundaries = cms.vdouble( 0.553, 0.902, 0.957 ),
                               SetArbitraryNonGoldMC = cms.bool(False),
                               DropNonGoldData = cms.bool(False),
                               RequireVBFPreselection = cms.bool(True),
                               VBFPreselLeadPtMin = cms.double(40.),
                               VBFPreselSubleadPtMin = cms.double(30.),
                               VBFPreselPhoIDMVAMin = cms.double(-0.2),
                               GetQCDWeights = cms.bool(False),
                               HTXSTags     = HTXSInputTags
                               )


flashggVHEtTag = cms.EDProducer("FlashggVHEtTagProducer",
                                RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                PATfilters = cms.InputTag('TriggerResults::PAT'),
                                FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                SystLabel=cms.string(""),
                                GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
                                MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                METTag=cms.InputTag('flashggMets'),   
                                useVertex0only=cms.bool(True),
                                leadPhoOverMassThreshold = cms.double(0.375),
                                subleadPhoOverMassThreshold = cms.double(0.25),
                                metPtThreshold = cms.double(70),
                                dPhiDiphotonMetThreshold = cms.double(2.1),
                                diphoMVAThreshold= cms.double(-1.0),
                                phoIdMVAThreshold= cms.double(-0.9),
                                 HTXSTags     = HTXSInputTags
                                #Boundaries=cms.vdouble(0.21,0.6,0.81)
)

flashggTTHLeptonicTag = cms.EDProducer("FlashggTTHLeptonicTagProducer",
                                       DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                       SystLabel=cms.string(""),
                                       MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                       inputTagJets= UnpackedJetCollectionVInputTag,
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
				       MetTag=cms.InputTag( 'flashggMets' ), 
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
				       rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                       leadPhoOverMassThreshold = cms.double(0.25),
                                       subleadPhoOverMassThreshold = cms.double(0.25),
                                       MVAThreshold = cms.double(-2.),
                                       PhoMVAThreshold = cms.double(-0.9), 
                                       jetsNumberThreshold = cms.double(3.),
                                       bjetsNumberThreshold = cms.double(1.),
                                       jetPtThreshold = cms.double(25.), 
                                       jetEtaThreshold= cms.double(2.4),
                                       deltaRJetLeadPhoThreshold = cms.double(0.4),
                                       deltaRJetSubLeadPhoThreshold = cms.double(0.4),
				       deltaRJetLepton = cms.double(0.4),
				       leadingJetPtThreshold = cms.double(50),
                                       bDiscriminator = bDiscriminator94X, #bDiscriminator76X,
                                       bTag = cms.string(flashggDeepCSV),
				       MuonEtaCut = cms.double(2.4),
				       MuonPtCut = cms.double(10),
				       MuonMiniIsoCut = cms.double(0.06),
				       MuonPhotonDrCut = cms.double(0.3),
				       EleEtaCuts = cms.vdouble(1.4442,1.566,2.5),
				       ElePtCut = cms.double(15),
				       ElePhotonDrCut = cms.double(0.3),
				       ElePhotonZMassCut = cms.double(5),
				       DeltaRTrkEle = cms.double(0.35),
                                       HTXSTags     = HTXSInputTags
)

flashggTTHDiLeptonTag = cms.EDProducer("FlashggTTHDiLeptonTagProducer",
                                       DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                       SystLabel=cms.string(""),
                                       MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                       inputTagJets= UnpackedJetCollectionVInputTag,
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
				       MetTag=cms.InputTag( 'flashggMets' ), 
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
				       rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                       leadPhoOverMassThreshold = cms.double(0.25),
                                       subleadPhoOverMassThreshold = cms.double(0.25),
                                       MVAThreshold = cms.double(-2.),
                                       PhoMVAThreshold = cms.double(-0.9), 
                                       jetsNumberThreshold = cms.double(0.),
                                       bjetsNumberThreshold = cms.double(0.),
				       bjetsLooseNumberThreshold = cms.double(0),
                                       jetPtThreshold = cms.double(25.), 
                                       jetEtaThreshold= cms.double(2.4),
                                       deltaRJetLeadPhoThreshold = cms.double(0.4),
                                       deltaRJetSubLeadPhoThreshold = cms.double(0.4),
				       deltaRJetLepton = cms.double(0.4),
				       leadingJetPtThreshold = cms.double(50),
                                       bDiscriminator = bDiscriminator94X, #bDiscriminator76X,
                                       bTag = cms.string(flashggDeepCSV),
				       MuonEtaCut = cms.double(2.4),
				       MuonPtCut = cms.double(10),
				       MuonPhotonDrCut = cms.double(0.2),
				       MuonsDrCut = cms.double(0.1),
				       MuonZMassCut = cms.double(5),
				       MuonsLeadingLeptonPtCut = cms.double(20),
				       EleEtaCuts = cms.vdouble(1.4442,1.566,2.5),
				       ElePtCut = cms.double(10),
				       ElePhotonDrCut = cms.double(0.2),
				       ElesDrCut = cms.double(0.2),
				       ElePhotonZMassCut = cms.double(5),
				       ElesZMassCut = cms.double(5),
				       ElesLeadingLeptonPtCut = cms.double(20),
				       EleMuDrCut = cms.double(0.1),
				       MixedLeadingLeptonPtCut = cms.double(15),
                                       HTXSTags     = HTXSInputTags
)

flashggVHLooseTag = cms.EDProducer("FlashggVHLooseTagProducer",
                                   DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                   SystLabel=cms.string(""),
                                   RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                   PATfilters = cms.InputTag('TriggerResults::PAT'),
                                   FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                   #JetTag=cms.InputTag('flashggSelectedJets'),
                                   inputTagJets= UnpackedJetCollectionVInputTag,
                                   ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                   MuonTag=cms.InputTag('flashggSelectedMuons'),
                                   VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                   METTag=cms.InputTag('flashggMets'),   
                                   useVertex0only=cms.bool(True),
                                   GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
				   rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                   leptonPtThreshold = cms.double(20),
                                   muonEtaThreshold = cms.double(2.4),
                                   leadPhoOverMassThreshold = cms.double(0.375),
                                   subleadPhoOverMassThreshold = cms.double(0.25),
                                   MVAThreshold = cms.double(0.187), #taken from cat2 boundary
                                   deltaRMuonPhoThreshold = cms.double(0.5),
                                   jetsNumberThreshold = cms.double(3.),
                                   jetPtThreshold = cms.double(20.),
                                   jetEtaThreshold= cms.double(2.4),
                                   deltaRPhoLeadJet = cms.double(0.5),
                                   deltaRPhoSubLeadJet = cms.double(0.5),
                                   muPFIsoSumRelThreshold = cms.double(0.25), 
                                   deltaRJetMuonThreshold = cms.double(0.5),
                                   PuIDCutoffThreshold = cms.double(0.8),
                                   PhoMVAThreshold = cms.double(-0.9),
                                   METThreshold = cms.double(45.),
                                   ElectronPtThreshold = cms.double(20.),
                                   DeltaRTrkElec = cms.double(.4),
                                   TransverseImpactParam = cms.double(0.02),
                                   LongitudinalImpactParam = cms.double(0.2),
                                   deltaRPhoElectronThreshold = cms.double(1.),
                                   deltaMassElectronZThreshold = cms.double(10.),
                                   electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                   nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                   nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                   electronIsoThreshold = cms.double(0.15),
                                   electronNumOfHitsThreshold = cms.double(1),
                                   useElectronMVARecipe = cms.bool(False),
                                   useElectronLooseID = cms.bool(True),
                                   HTXSTags     = HTXSInputTags
				    )

flashggVHTightTag = cms.EDProducer("FlashggVHTightTagProducer",
                                   DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                   SystLabel=cms.string(""),
                                   RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                   PATfilters = cms.InputTag('TriggerResults::PAT'),
                                   FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                   #JetTag=cms.InputTag('flashggSelectedJets'),
                                   inputTagJets= UnpackedJetCollectionVInputTag,
                                   ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                   MuonTag=cms.InputTag('flashggSelectedMuons'),
                                   VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                   METTag=cms.InputTag('flashggMets'),   
                                   useVertex0only=cms.bool(True),
                                   GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
				   rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                   leptonPtThreshold = cms.double(20),
                                   muonEtaThreshold = cms.double(2.4),
                                   leadPhoOverMassThreshold = cms.double(0.375),
                                   subleadPhoOverMassThreshold = cms.double(0.25),
                                   MVAThreshold = cms.double(0.187), #taken from cat 3 boundary
                                   deltaRMuonPhoThreshold = cms.double(1),
                                   jetsNumberThreshold = cms.double(3.),
                                   jetPtThreshold = cms.double(20.),
                                   jetEtaThreshold= cms.double(2.4),
                                   muPFIsoSumRelThreshold = cms.double(0.25), 
                                   PuIDCutoffThreshold = cms.double(0.8),
                                   PhoMVAThreshold = cms.double(-0.9),
                                   METThreshold = cms.double(45.),
                                   deltaRJetMuonThreshold = cms.double(0.5),
                                   deltaRJetElectronThreshold = cms.double(0.5),
                                   invMassLepLowThreshold = cms.double(70.),
                                   invMassLepHighThreshold = cms.double(110.),
                                   numberOfLowPtMuonsThreshold = cms.double(2.),
                                   numberOfHighPtMuonsThreshold = cms.double(1.),
                                   leptonLowPtThreshold = cms.double(10.),
                                   deltaRLowPtMuonPhoThreshold = cms.double(0.5),
                                   deltaRPhoLeadJet = cms.double(0.5),
                                   deltaRPhoSubLeadJet = cms.double(0.5),
                                   ElectronPtThreshold = cms.double(20.),
                                   DeltaRTrkElec = cms.double(0.4),
                                   TransverseImpactParam = cms.double(0.02),
                                   LongitudinalImpactParam = cms.double(0.2),
                                   deltaRPhoElectronThreshold = cms.double(1.),
                                   deltaMassElectronZThreshold = cms.double(10.),
                                   electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                   nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                   nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                   electronIsoThreshold = cms.double(0.15),
                                   electronNumOfHitsThreshold = cms.double(1),
                                   useElectronMVARecipe = cms.bool(False),
                                   useElectronLooseID = cms.bool(True),
                                   HTXSTags     = HTXSInputTags
)

flashggVHMetTag = cms.EDProducer("FlashggVHMetTagProducer",
                                 RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                 PATfilters = cms.InputTag('TriggerResults::PAT'),
                                 FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                 inputTagJets= UnpackedJetCollectionVInputTag,
                                 DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                 SystLabel=cms.string(""),
                                 GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
                                 MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                 METTag=cms.InputTag('flashggMets'),
                                 useVertex0only=cms.bool(False),
                                 leadPhoOverMassThreshold = cms.double(0.375),
                                 subleadPhoOverMassThreshold = cms.double(0.25),
                                 metPtThreshold = cms.double(85),
                                 dPhiDiphotonMetThreshold = cms.double(2.4),
                                 dPhiJetMetThreshold = cms.double(999.),
                                 jetPtThreshold = cms.double(50.),
                                 jetEtaThreshold= cms.double(2.4),
                                 deltaRPhoLeadJet = cms.double(0.5),
                                 deltaRPhoSubLeadJet = cms.double(0.5),
                                 diphoMVAThreshold= cms.double(0.6),
                                 phoIdMVAThreshold= cms.double(-0.9),
                                 HTXSTags     = HTXSInputTags
                                 #Boundaries=cms.vdouble(0.21,0.6,0.81)                                                                            
)

flashggZHLeptonicTag = cms.EDProducer("FlashggZHLeptonicTagProducer",
                                   DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                   SystLabel=cms.string(""),
                                   inputTagJets= UnpackedJetCollectionVInputTag,
                                   ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                   MuonTag=cms.InputTag('flashggSelectedMuons'),
                                   VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                   METTag=cms.InputTag('flashggMets'),
                                   useVertex0only=cms.bool(False),
                                   GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
				   rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                   leptonPtThreshold = cms.double(20),
                                   muonEtaThreshold = cms.double(2.4),
                                   leadPhoOverMassThreshold = cms.double(0.375),
                                   subleadPhoOverMassThreshold = cms.double(0.25),
                                   MVAThreshold = cms.double(-0.405),
                                   deltaRMuonPhoThreshold = cms.double(1),
                                   muPFIsoSumRelThreshold = cms.double(0.25),
                                   PhoMVAThreshold = cms.double(-0.9),
                                   invMassLepLowThreshold = cms.double(70.),
                                   invMassLepHighThreshold = cms.double(110.),
                                   deltaRLowPtMuonPhoThreshold = cms.double(0.5),
                                   ElectronPtThreshold = cms.double(20.),
                                   DeltaRTrkElec = cms.double(0.4),
                                   TransverseImpactParam = cms.double(0.02),
                                   LongitudinalImpactParam = cms.double(0.2),
                                   deltaRPhoElectronThreshold = cms.double(1.),
                                   deltaMassElectronZThreshold = cms.double(10.),
                                   electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                   nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                   nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                   electronIsoThreshold = cms.double(0.15),
                                   electronNumOfHitsThreshold = cms.double(1),
                                   useElectronMVARecipe = cms.bool(False),
                                   useElectronLooseID = cms.bool(True),
                                      HTXSTags     = HTXSInputTags
)

flashggWHLeptonicTag = cms.EDProducer("FlashggWHLeptonicTagProducer",
                                   DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                   SystLabel=cms.string(""),
                                   RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                   PATfilters = cms.InputTag('TriggerResults::PAT'),
                                   FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                   inputTagJets= UnpackedJetCollectionVInputTag,
                                   ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                   MuonTag=cms.InputTag('flashggSelectedMuons'),
                                   VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                   METTag=cms.InputTag('flashggMets'),
                                   useVertex0only=cms.bool(False),
                                   GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
				   rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                   leptonPtThreshold = cms.double(20),
                                   muonEtaThreshold = cms.double(2.4),
                                   leadPhoOverMassThreshold = cms.double(0.375),
                                   subleadPhoOverMassThreshold = cms.double(0.25),
                                   MVAThreshold = cms.double(0.0),                                                     
                                   deltaRMuonPhoThreshold = cms.double(0.5),
                                   jetsNumberThreshold = cms.double(3.),
                                   jetPtThreshold = cms.double(20.),
                                   jetEtaThreshold= cms.double(2.4),
                                   deltaRPhoLeadJet = cms.double(0.4),
                                   deltaRPhoSubLeadJet = cms.double(0.4),
                                   muPFIsoSumRelThreshold = cms.double(0.25),
                                   deltaRJetMuonThreshold = cms.double(0.4),
                                   PuIDCutoffThreshold = cms.double(0.8),
                                   PhoMVAThreshold = cms.double(-0.9),
                                   METThreshold = cms.double(45.),
                                   DeltaRTrkElec = cms.double(.4),
                                   TransverseImpactParam = cms.double(0.02),
                                   LongitudinalImpactParam = cms.double(0.2),
                                   deltaRPhoElectronThreshold = cms.double(1.),
                                   deltaMassElectronZThreshold = cms.double(10.),
                                   electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                   nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                   nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                   electronIsoThreshold = cms.double(0.15),
                                   electronNumOfHitsThreshold = cms.double(1),
                                   useElectronMVARecipe = cms.bool(False),
                                   useElectronLooseID = cms.bool(True),
                                      HTXSTags     = HTXSInputTags
                                    )
flashggVHLeptonicLooseTag = cms.EDProducer("FlashggVHLeptonicLooseTagProducer",
                                   DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                   SystLabel=cms.string(""),
                                   RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                   PATfilters = cms.InputTag('TriggerResults::PAT'),
                                   FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                   inputTagJets= UnpackedJetCollectionVInputTag,
                                   ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                   MuonTag=cms.InputTag('flashggSelectedMuons'),
                                   VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                   METTag=cms.InputTag('flashggMets'),
                                   useVertex0only=cms.bool(False),
                                   GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
                                   rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
				   leptonPtThreshold = cms.double(20),
                                   muonEtaThreshold = cms.double(2.4),
                                   leadPhoOverMassThreshold = cms.double(0.375),
                                   subleadPhoOverMassThreshold = cms.double(0.25),
                                   MVAThreshold = cms.double(0.0),
                                   deltaRMuonPhoThreshold = cms.double(1),
                                   jetsNumberThreshold = cms.double(3.),
                                   jetPtThreshold = cms.double(20.),
                                   jetEtaThreshold= cms.double(2.4),
                                   muPFIsoSumRelThreshold = cms.double(0.25),
                                   PuIDCutoffThreshold = cms.double(0.8),
                                   PhoMVAThreshold = cms.double(-0.9),
                                   METThreshold = cms.double(45.),
                                   deltaRJetMuonThreshold = cms.double(0.4),
                                   deltaRJetElectronThreshold = cms.double(0.4),
                                   invMassLepLowThreshold = cms.double(70.),
                                   invMassLepHighThreshold = cms.double(110.),
                                   deltaRPhoLeadJet = cms.double(0.4),
                                   deltaRPhoSubLeadJet = cms.double(0.4),
                                   DeltaRTrkElec = cms.double(0.4),
                                   TransverseImpactParam = cms.double(0.02),
                                   LongitudinalImpactParam = cms.double(0.2),
                                   deltaRPhoElectronThreshold = cms.double(1.),
                                   deltaMassElectronZThreshold = cms.double(10.),
                                   electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                   nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                   nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                   electronIsoThreshold = cms.double(0.15),
                                   electronNumOfHitsThreshold = cms.double(1),
                                   useElectronMVARecipe = cms.bool(False),
                                   useElectronLooseID = cms.bool(True),
                                           HTXSTags     = HTXSInputTags
)



flashggVHHadronicTag = cms.EDProducer("FlashggVHHadronicTagProducer",
                                      DiPhotonTag = cms.InputTag('flashggPreselectedDiPhotons'),
                                      SystLabel=cms.string(""),
                                      MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                      #JetTag = cms.InputTag('flashggSelectedJets'),
                                      inputTagJets= UnpackedJetCollectionVInputTag,
                                      GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
                                      leadPhoOverMassThreshold = cms.double(0.5),
                                      subleadPhoOverMassThreshold = cms.double(0.25),
                                      diphoMVAThreshold = cms.double(0.6),
                                      jetsNumberThreshold = cms.double(2.),
                                      jetPtThreshold = cms.double(40.),
                                      jetEtaThreshold= cms.double(2.4),
                                      dRJetToPhoLThreshold = cms.double(0.4),
                                      dRJetToPhoSThreshold = cms.double(0.4),
                                      dijetMassLowThreshold = cms.double(60.),
                                      dijetMassHighThreshold = cms.double(120.),
                                      cosThetaStarThreshold = cms.double(0.5),
                                      phoIdMVAThreshold = cms.double(-0.9),
                                      HTXSTags     = HTXSInputTags
)

# Tag is for jet studies only, not in default sequence
flashggZPlusJetTag = cms.EDProducer("FlashggZPlusJetTagProducer",
                                    DiPhotonTag    = cms.InputTag('flashggPreselectedDiPhotons'),
                                    SystLabel      = cms.string(""),
                                    MVAResultTag   = cms.InputTag('flashggDiPhotonMVA'),
                                    inputTagJets= UnpackedJetCollectionVInputTag,
                                    GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
                                    GenJetTag = cms.InputTag("slimmedGenJets")
                                    )

