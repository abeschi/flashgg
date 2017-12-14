#-----------J. Tao from IHEP-Beijing--------------

import FWCore.ParameterSet.Config as cms

from flashgg.Taggers.TTHEfficiencyDumpConfig_cff import TTHHadronicEfficiencyDumpConfig, TTHSemiLeptonicEfficiencyDumpConfig, TTHFullyLeptonicEfficiencyDumpConfig

TTHHadronicEfficiencyDumper = cms.EDAnalyzer('CutBasedTTHHadronicEfficiencyDumper',
                                **TTHHadronicEfficiencyDumpConfig.parameters_()
                                )

TTHSemiLeptonicEfficiencyDumper = cms.EDAnalyzer('CutBasedTTHSemiLeptonicEfficiencyDumper',
                                **TTHSemiLeptonicEfficiencyDumpConfig.parameters_()
                                )

TTHFullyLeptonicEfficiencyDumper = cms.EDAnalyzer('CutBasedTTHFullyLeptonicEfficiencyDumper',
                                **TTHFullyLeptonicEfficiencyDumpConfig.parameters_()
                                )

