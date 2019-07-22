import FWCore.ParameterSet.Config as cms
from flashgg.Taggers.flashggZGammaTags_cff import *
from flashgg.Taggers.flashggZGammaTagSorter_cfi import flashggZGammaTagSorter
from flashgg.MicroAOD.flashggEleEleGamma_cfi import flashggEleEleGamma
from flashgg.MicroAOD.flashggMuMuGamma_cfi import flashggMuMuGamma
#At some point introduce singlePhoton and DiElectrons and DiMuons corrections
#VBF MVA could be included too

    
flashggZGammaTagSequence = cms.Sequence( flashggEleEleGamma*					 
					 flashggMuMuGamma*
					 ( flashggZGammaToEEUntaggedTag
                                         + flashggZGammaToMuMuUntaggedTag
                                             )
                                             * flashggZGammaTagSorter
                                           )


