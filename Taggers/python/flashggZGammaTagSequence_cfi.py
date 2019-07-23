import FWCore.ParameterSet.Config as cms
from flashgg.Taggers.flashggZGammaTags_cff import *
from flashgg.Taggers.flashggZGammaTagSorter_cfi import flashggZGammaTagSorter
from flashgg.MicroAOD.flashggEleEleGamma_cfi import flashggEleEleGamma
from flashgg.MicroAOD.flashggMuMuGamma_cfi import flashggMuMuGamma
from flashgg.MicroAOD.flashggDiMuons_cfi import flashggDiMuons
from flashgg.MicroAOD.flashggDiElectrons_cfi import flashggDiElectrons
#At some point introduce singlePhoton and DiElectrons and DiMuons corrections
#VBF MVA could be included too

    
flashggZGammaTagSequence = cms.Sequence( flashggDiMuons*
					 flashggDiElectrons*
					 flashggMuMuGamma*
					 flashggEleEleGamma*					 
					 ( flashggZGammaToMuMuUntaggedTag
                                         + flashggZGammaToEEUntaggedTag
                                             )
                                             * flashggZGammaTagSorter
                                           )


