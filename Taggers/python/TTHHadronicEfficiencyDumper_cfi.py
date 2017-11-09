#-----------J. Tao from IHEP-Beijing--------------

import FWCore.ParameterSet.Config as cms

from flashgg.Taggers.TTHHadronicEfficiencyDumpConfig_cff import TTHHadronicEfficiencyDumpConfig

TTHHadronicEfficiencyDumper = cms.EDAnalyzer('CutBasedTTHHadronicEfficiencyDumper',
                                **TTHHadronicEfficiencyDumpConfig.parameters_()
                                )


