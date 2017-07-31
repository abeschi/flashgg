#-----------J. Tao from IHEP-Beijing--------------

import FWCore.ParameterSet.Config as cms

flashggMuMuGammaRandomizedPhotons2 = cms.EDProducer('FlashggMuMuGammaRandomizedPhotonProducer',
                                  DiMuonTag=cms.InputTag('flashggDiMuons'),
                                  PhotonTag=cms.InputTag('flashggPreselectedPhotons'),
                                  VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                  ##Parameters
				  leadingMuonMinPt=cms.double(10.5),
				  subleadingMuonMinPt=cms.double(10.5),
				  MuonMinSumPt=cms.double(30.),
				  leadingMuonIsoOverPt=cms.double(0.2),
				  subleadingMuonIsoOverPt=cms.double(0.2),
				  MinDimuonMass=cms.double(35.),
				  MinMuMuGammaMass=cms.double(60.),
				  MaxMuMuGammaMass=cms.double(120.),
				  MaxMuMuPlusMuMuGammaMass=cms.double(180.),
				  MaxDeltaR=cms.double(0.8),
				  MinClosestMuonPt=cms.double(21.),                                        
                                  minPhotonPT=cms.double(20.),
                                )
