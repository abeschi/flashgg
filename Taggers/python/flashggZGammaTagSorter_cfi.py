import FWCore.ParameterSet.Config as cms
from flashgg.Taggers.flashggZGammaTags_cff import *

flashggZGammaTagSorter = cms.EDProducer('FlashggZGammaTagSorter',
                                  # Top of list is highest priority
                                  # Optionally can add category ranges if priority depends on category number
                                  TagPriorityRanges = cms.VPSet(
				        cms.PSet(TagName = cms.InputTag('flashggZGammaToMuMuUntaggedTag')), 
				        cms.PSet(TagName = cms.InputTag('flashggZGammaToEEUntaggedTag')), 
					),
                                  MinObjectWeightException = cms.double(0.0),
                                  MaxObjectWeightException = cms.double(10.),
                                  MinObjectWeightWarning = cms.double(0.4),
                                  MaxObjectWeightWarning = cms.double(2.5),
                                  Debug = cms.untracked.bool(False),
                                  CreateNoTag = cms.bool(False),  # Placeholder for tracking rejected events
                                  )

