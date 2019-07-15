#-----------J. Tao from IHEP-Beijing--------------

import FWCore.ParameterSet.Config as cms

flashggMuMuGamma = cms.EDProducer('FlashggMuMuGammaProducer',
                                  DiMuonTag=cms.InputTag('flashggDiMuons'),
                                  #PhotonTag=cms.InputTag('flashggPhotons'),
                                  PhotonTag=cms.InputTag('flashggRandomizedPhotons'),
                                  VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                  ##Parameters                                                
                                  minPhotonPT=cms.double(10.)
                                  )
