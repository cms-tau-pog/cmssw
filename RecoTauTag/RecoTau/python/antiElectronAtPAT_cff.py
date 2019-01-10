##########################################################
# Define new anti-e tau-IDs to be embedded into patTaus
#
# M. Bluj, NCBJ, Poland
# September 2018
##########################################################

import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras

### Define PoolDBESSource with new payloads
from RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi import *
loadRecoTauTagMVAsFromPrepDB_mva6 = loadRecoTauTagMVAsFromPrepDB.clone(
    connect = 'sqlite_file:'+os.getenv('CMSSW_BASE')+'/src/RecoTauTag/TrainingFiles/data/AntiETauId/RecoTauTag_MVAs_2018Oct01.db',
    toGet   = cms.VPSet()
)
antiElectronDiscrMVA6_categories = {
     '0' : "gbr_NoEleMatch_woGwoGSF_BL",
     '2' : "gbr_NoEleMatch_wGwoGSF_BL",
     '5' : "gbr_woGwGSF_BL",
     '7' : "gbr_wGwGSF_BL",
     '8' : "gbr_NoEleMatch_woGwoGSF_EC",
    '10' : "gbr_NoEleMatch_wGwoGSF_EC",
    '13' : "gbr_woGwGSF_EC",
    '15' : "gbr_wGwGSF_EC"
}
#antiElectronDiscrMVA6_WPs = [ "Eff99", "Eff96", "Eff91", "Eff85", "Eff79" ]
antiElectronDiscrMVA6_WPs = [ "eff98", "eff90", "eff80", "eff70", "eff60" ]
antiElectronDiscrMVA6_version = "v2"
for category, gbrForestName in antiElectronDiscrMVA6_categories.items():
    loadRecoTauTagMVAsFromPrepDB_mva6.toGet.append(
        cms.PSet(
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string("RecoTauTag_antiElectronMVA6%s_%s" % (antiElectronDiscrMVA6_version, gbrForestName)),
            label = cms.untracked.string("RecoTauTag_antiElectronMVA6%s_%s" % (antiElectronDiscrMVA6_version, gbrForestName))
        )
    )
    for WP in antiElectronDiscrMVA6_WPs:
        loadRecoTauTagMVAsFromPrepDB_mva6.toGet.append(
            cms.PSet(
                record = cms.string('PhysicsTGraphPayloadRcd'),
                tag = cms.string("RecoTauTag_antiElectronMVA6%s_%s_WP%s" % (antiElectronDiscrMVA6_version, gbrForestName, WP)),
                label = cms.untracked.string("RecoTauTag_antiElectronMVA6%s_%s_WP%s" % (antiElectronDiscrMVA6_version, gbrForestName, WP))
            )
        )

### Define new discriminants
## Raw
from RecoTauTag.RecoTau.PATTauDiscriminationAgainstElectronMVA6_cfi import patTauDiscriminationAgainstElectronMVA6
from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
patTauDiscriminationByElectronRejectionMVA6v2Raw = patTauDiscriminationAgainstElectronMVA6.clone(
    Prediscriminants = noPrediscriminants, #already selected for MiniAOD
    mvaName_NoEleMatch_wGwoGSF_BL = 'RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_BL',
    mvaName_NoEleMatch_wGwoGSF_EC = 'RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_EC',
    mvaName_NoEleMatch_woGwoGSF_BL = 'RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_BL',
    mvaName_NoEleMatch_woGwoGSF_EC = 'RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_EC',
    mvaName_wGwGSF_BL = 'RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_BL',
    mvaName_wGwGSF_EC = 'RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_EC',
    mvaName_woGwGSF_BL = 'RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_BL',
    mvaName_woGwGSF_EC = 'RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_EC'
)
## WPs 
from RecoTauTag.RecoTau.PATTauDiscriminantCutMultiplexer_cfi import patTauDiscriminantCutMultiplexer
# VLoose
patTauDiscriminationByVLooseElectronRejectionMVA6v2 = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = patTauDiscriminationByElectronRejectionMVA6v2Raw.PATTauProducer,
    Prediscriminants = patTauDiscriminationByElectronRejectionMVA6v2Raw.Prediscriminants,
    toMultiplex = cms.InputTag("patTauDiscriminationByElectronRejectionMVA6v2Raw"),
    key = cms.InputTag("patTauDiscriminationByElectronRejectionMVA6v2Raw","category"),
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_BL_WPeff98'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(2),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_BL_WPeff98'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(5),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_BL_WPeff98'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(7),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_BL_WPeff98'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(8),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_EC_WPeff98'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(10),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_EC_WPeff98'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(13),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_EC_WPeff98'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(15),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_EC_WPeff98'),
            variable = cms.string('pt')
        )
    )
)
# Loose
patTauDiscriminationByLooseElectronRejectionMVA6v2 = patTauDiscriminationByVLooseElectronRejectionMVA6v2.clone(
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_BL_WPeff90'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(2),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_BL_WPeff90'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(5),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_BL_WPeff90'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(7),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_BL_WPeff90'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(8),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_EC_WPeff90'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(10),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_EC_WPeff90'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(13),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_EC_WPeff90'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(15),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_EC_WPeff90'),
            variable = cms.string('pt')
        )
    )    
)
# Medium
patTauDiscriminationByMediumElectronRejectionMVA6v2 = patTauDiscriminationByVLooseElectronRejectionMVA6v2.clone(
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_BL_WPeff80'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(2),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_BL_WPeff80'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(5),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_BL_WPeff80'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(7),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_BL_WPeff80'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(8),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_EC_WPeff80'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(10),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_EC_WPeff80'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(13),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_EC_WPeff80'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(15),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_EC_WPeff80'),
            variable = cms.string('pt')
        )
    )    
)
# Tight
patTauDiscriminationByTightElectronRejectionMVA6v2 = patTauDiscriminationByVLooseElectronRejectionMVA6v2.clone(
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_BL_WPeff70'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(2),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_BL_WPeff70'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(5),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_BL_WPeff70'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(7),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_BL_WPeff70'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(8),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_EC_WPeff70'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(10),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_EC_WPeff70'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(13),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_EC_WPeff70'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(15),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_EC_WPeff70'),
            variable = cms.string('pt')
        )
    )    
)
# VTight
patTauDiscriminationByVTightElectronRejectionMVA6v2 = patTauDiscriminationByVLooseElectronRejectionMVA6v2.clone(
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_BL_WPeff60'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(2),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_BL_WPeff60'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(5),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_BL_WPeff60'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(7),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_BL_WPeff60'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(8),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_woGwoGSF_EC_WPeff60'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(10),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_NoEleMatch_wGwoGSF_EC_WPeff60'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(13),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_woGwGSF_EC_WPeff60'),
            variable = cms.string('pt')
        ), 
        cms.PSet(
            category = cms.uint32(15),
            cut = cms.string('RecoTauTag_antiElectronMVA6v2_gbr_wGwGSF_EC_WPeff60'),
            variable = cms.string('pt')
        )
    )    
)
### Put all this stuff to a sequence
patTauDiscriminationByElectronRejectionTask = cms.Task(
    patTauDiscriminationByElectronRejectionMVA6v2Raw,
    patTauDiscriminationByVLooseElectronRejectionMVA6v2,
    patTauDiscriminationByLooseElectronRejectionMVA6v2,
    patTauDiscriminationByMediumElectronRejectionMVA6v2,
    patTauDiscriminationByTightElectronRejectionMVA6v2,
    patTauDiscriminationByVTightElectronRejectionMVA6v2
)
patTauDiscriminationByElectronRejectionSeq = cms.Sequence(patTauDiscriminationByElectronRejectionTask)

againstElectronTauIDSources = cms.PSet(
    againstElectronMVA6Raw2018 = cms.InputTag("patTauDiscriminationByElectronRejectionMVA6v2Raw"),
    againstElectronMVA6category2018 = cms.InputTag("patTauDiscriminationByElectronRejectionMVA6v2Raw","category"),
    againstElectronVLooseMVA62018 = cms.InputTag("patTauDiscriminationByVLooseElectronRejectionMVA6v2"),
    againstElectronLooseMVA62018 = cms.InputTag("patTauDiscriminationByLooseElectronRejectionMVA6v2"),
    againstElectronMediumMVA62018 = cms.InputTag("patTauDiscriminationByMediumElectronRejectionMVA6v2"),
    againstElectronTightMVA62018 = cms.InputTag("patTauDiscriminationByTightElectronRejectionMVA6v2"),
    againstElectronVTightMVA62018 = cms.InputTag("patTauDiscriminationByVTightElectronRejectionMVA6v2")
)
