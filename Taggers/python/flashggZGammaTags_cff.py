
import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, flashggDeepCSV, maxJetCollections
bDiscriminator94X= cms.vdouble(0.1522,0.4941,0.8001)

flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
                                     JetsTag = cms.InputTag("flashggFinalJets"),
                                     NCollections = cms.uint32(maxJetCollections)
                                     )

UnpackedJetCollectionVInputTag = cms.VInputTag()
for i in range(0,maxJetCollections):
    UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))

#Possible values: Veto, Loose, Medium, Tight, MVALoose, MVAMedium, MVATight, MVANoIsoLoose, MVANoIsoMedium, MVANoIsoTight
EleFromZId = cms.string("MVAMedium") 
#Possible values: Loose, Medium, Tight
MuFromZId = cms.string("Medium") 

flashggZGammaToEEUntaggedTag = cms.EDProducer("FlashggZGammaToEEUntaggedTagProducer",
                                       SystLabel=cms.string(""),
                                       inputTagJets= UnpackedJetCollectionVInputTag,
                                       EleEleGammaTag=cms.InputTag('flashggEleEleGamma'),
                                       PhotonTag=cms.InputTag('flashggRandomizedPhotons'),
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
				       MetTag=cms.InputTag( 'flashggMets' ), 
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
#                                       MVAweightfile = cms.FileInPath("flashgg/Taggers/data/TMVAClassification_BDT_training_v2.json.weights.xml"),
				       debug = cms.bool(False),
				       useEleEleGammaCandidate = cms.bool(True),
				       chooseByZMass = cms.bool(False),
				       preselectedElePt = cms.double(7),
				       preselectedEleEta = cms.double(2.5),
				       preselectedPhoPt = cms.double(10),
				       MinLeadElePt = cms.double(25),
				       MinSubleadElePt = cms.double(15),
				       MaxEleEta = cms.double(2.5),
				       minDiEleMass = cms.double(50.),
				       minLeptDR = cms.double(0.4),
                                       EleIdWP = EleFromZId,
				       minEEGMass = cms.double(100.),
				       maxEEGMass = cms.double(180.),
				       minEEGPlusDiEleMass = cms.double(185.),
				       minPhoId = cms.double(0.0),
				       minPtGammaOverMass = cms.double(0.14),
				       minLeptPhoDR = cms.double(0.4),
				       mvaBoundaries = cms.vdouble(-2., 0.3, 0.7) #events lower than the first number are thrown away
)

flashggZGammaToMuMuUntaggedTag = cms.EDProducer("FlashggZGammaToMuMuUntaggedTagProducer",
                                       SystLabel=cms.string(""),
                                       inputTagJets= UnpackedJetCollectionVInputTag,
                                       MuMuGammaTag=cms.InputTag('flashggMuMuGamma'),
                                       PhotonTag=cms.InputTag('flashggRandomizedPhotons'),
                                       ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                       MuonTag=cms.InputTag('flashggSelectedMuons'),
				       MetTag=cms.InputTag( 'flashggMets' ), 
                                       VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
#                                       MVAweightfile = cms.FileInPath("flashgg/Taggers/data/TMVAClassification_BDT_training_v2.json.weights.xml"),
				       debug = cms.bool(False),
				       useMuMuGammaCandidate = cms.bool(True),
				       chooseByZMass = cms.bool(False),
				       preselectedMuonPt = cms.double(7),
				       preselectedMuonEta = cms.double(2.5),
				       preselectedPhoPt = cms.double(10),
				       MinLeadMuonPt = cms.double(20),
				       MinSubleadMuonPt = cms.double(10),
				       MaxMuonEta = cms.double(2.5),
				       MuonIso = cms.double(0.4),
				       minDiMuonMass = cms.double(50.),
				       minLeptDR = cms.double(0.4),
                                       MuIdWP = MuFromZId,
				       minMMGMass = cms.double(100.),
				       maxMMGMass = cms.double(180.),
				       minMMGPlusDiMuonMass = cms.double(185.),
				       minPhoId = cms.double(0.0),
				       minPtGammaOverMass = cms.double(0.14),
				       minLeptPhoDR = cms.double(0.4),
				       mvaBoundaries = cms.vdouble(-2., 0.3, 0.7) #events lower than the first number are thrown away
)


