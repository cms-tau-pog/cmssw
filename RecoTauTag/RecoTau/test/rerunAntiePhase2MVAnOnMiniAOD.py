##########################################################
# (Re)run anti-e tau-IDs payloads for phase-2
#
# M. Bluj, NCBJ, Poland
# January 2020
##########################################################

import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("rerunAntiePhase2MVAnOnMiniAOD")

### Geometry and Detector Conditions (needed for a few tau reco steps)
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2023_realistic_v2')
#process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v5')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

### Source files
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:patMiniAOD_standard.root'
        #'/store/relval/CMSSW_10_6_0/RelValTTbar_14TeV/MINIAODSIM/106X_upgrade2023_realistic_v2_2023D41noPU-v1/10000/63866134-439C-FF4F-B4DF-7398A8787500.root'
        'root://xrootd.unl.edu//store/mc/PhaseIISpr18AODMiniAOD/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/10000/2A39D829-8C65-E811-A420-0025905A60A6.root'
    )
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('rerunTauID_miniAOD_phase2.root')
)

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
### Load PoolDBESSource with mapping to payloads
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')

# MVA phase2
antiElectronDiscrMVA_phase2_categories = {
     '0' : "gbr_NoEleMatch_woGwoGSF_BL",
     '2' : "gbr_NoEleMatch_wGwoGSF_BL",
     '5' : "gbr_woGwGSF_BL",
     '7' : "gbr_wGwGSF_BL",
     '8' : "gbr_NoEleMatch_woGwoGSF_FWEC",
     '9' : "gbr_NoEleMatch_woGwoGSF_VFWEC",
    '10' : "gbr_NoEleMatch_wGwoGSF_FWEC",
    '11' : "gbr_NoEleMatch_wGwoGSF_VFWEC",
    '13' : "gbr_woGwGSF_FWEC",
    '14' : "gbr_woGwGSF_VFWEC",
    '15' : "gbr_wGwGSF_FWEC",
    '16' : "gbr_wGwGSF_VFWEC"
}
antiElectronDiscrMVA_phase2_WPs = [ "Eff98", "Eff90", "Eff80", "Eff70", "Eff60" ]
antiElectronDiscrMVA_phase2_version = "v1"

###
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import *
from RecoTauTag.RecoTau.PATTauDiscriminationAgainstElectronMVA6_cfi import *

process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw = patTauDiscriminationAgainstElectronMVA6.clone(
    PATTauProducer = cms.InputTag('slimmedTaus'),
    Prediscriminants = noPrediscriminants,
    #Prediscriminants = requireLeadTrack,
    vetoEcalCracks = cms.bool(False),
    returnMVA = cms.bool(True),
    method = cms.string("BDTG"),
    mvaName_NoEleMatch_woGwoGSF_BL = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_NoEleMatch_woGwoGSF_BL"),
    mvaName_NoEleMatch_wGwoGSF_BL = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_NoEleMatch_wGwoGSF_BL"),
    mvaName_woGwGSF_BL = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_woGwGSF_BL"),
    mvaName_wGwGSF_BL = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_wGwGSF_BL"),
    mvaName_NoEleMatch_woGwoGSF_EC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_NoEleMatch_woGwoGSF_FWEC"),
    mvaName_NoEleMatch_wGwoGSF_EC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_NoEleMatch_wGwoGSF_FWEC"),
    mvaName_woGwGSF_EC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_woGwGSF_FWEC"),
    mvaName_wGwGSF_EC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_wGwGSF_FWEC"),
    mvaName_NoEleMatch_woGwoGSF_VFEC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_NoEleMatch_woGwoGSF_VFWEC"),
    mvaName_NoEleMatch_wGwoGSF_VFEC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_NoEleMatch_wGwoGSF_VFWEC"),
    mvaName_woGwGSF_VFEC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_woGwGSF_VFWEC"),
    mvaName_wGwGSF_VFEC = cms.string("RecoTauTag_antiElectronPhase2MVA6v1_gbr_wGwGSF_VFWEC"),  
    minMVANoEleMatchWOgWOgsfBL = cms.double(0.0),
    minMVANoEleMatchWgWOgsfBL  = cms.double(0.0),
    minMVAWOgWgsfBL            = cms.double(0.0),
    minMVAWgWgsfBL             = cms.double(0.0),
    minMVANoEleMatchWOgWOgsfEC = cms.double(0.0),
    minMVANoEleMatchWgWOgsfEC  = cms.double(0.0),
    minMVAWOgWgsfEC            = cms.double(0.0),
    minMVAWgWgsfEC             = cms.double(0.0),
    minMVANoEleMatchWOgWOgsfVFEC = cms.double(0.0),
    minMVANoEleMatchWgWOgsfVFEC  = cms.double(0.0),
    minMVAWOgWgsfVFEC		 = cms.double(0.0),
    minMVAWgWgsfVFEC		 = cms.double(0.0),
    srcElectrons = cms.InputTag('slimmedElectronsFromMultiCl'),
    #usePhiAtEcalEntranceExtrapolation = cms.bool(True)
)
## WPs 
from RecoTauTag.RecoTau.PATTauDiscriminantCutMultiplexer_cfi import patTauDiscriminantCutMultiplexer
# VLoose
WP='Eff98'
process.patTauDiscriminationByVLooseElectronRejectionPhase2MVA6v1 = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.PATTauProducer,
    Prediscriminants = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.Prediscriminants,
    toMultiplex = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw"),
    key = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw","category"),
    mapping = cms.VPSet()
)
for category, gbrForestName in antiElectronDiscrMVA_phase2_categories.items():
    pset = cms.PSet(
        category = cms.uint32(int(category)),
        cut = cms.string('RecoTauTag_antiElectronPhase2MVA6%s_%s_WP%s' % (antiElectronDiscrMVA_phase2_version, gbrForestName, WP)),
        variable = cms.string('pt')
    )
    process.patTauDiscriminationByVLooseElectronRejectionPhase2MVA6v1.mapping.append(pset)
# Loose
WP='Eff90'
process.patTauDiscriminationByLooseElectronRejectionPhase2MVA6v1 = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.PATTauProducer,
    Prediscriminants = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.Prediscriminants,
    toMultiplex = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw"),
    key = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw","category"),
    mapping = cms.VPSet()
)
for category, gbrForestName in antiElectronDiscrMVA_phase2_categories.items():
    pset = cms.PSet(
        category = cms.uint32(int(category)),
        cut = cms.string('RecoTauTag_antiElectronPhase2MVA6%s_%s_WP%s' % (antiElectronDiscrMVA_phase2_version, gbrForestName, WP)),
        variable = cms.string('pt')
    )
    process.patTauDiscriminationByLooseElectronRejectionPhase2MVA6v1.mapping.append(pset)
# Medium
WP='Eff80'
process.patTauDiscriminationByMediumElectronRejectionPhase2MVA6v1 = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.PATTauProducer,
    Prediscriminants = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.Prediscriminants,
    toMultiplex = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw"),
    key = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw","category"),
    mapping = cms.VPSet()
)
for category, gbrForestName in antiElectronDiscrMVA_phase2_categories.items():
    pset = cms.PSet(
        category = cms.uint32(int(category)),
        cut = cms.string('RecoTauTag_antiElectronPhase2MVA6%s_%s_WP%s' % (antiElectronDiscrMVA_phase2_version, gbrForestName, WP)),
        variable = cms.string('pt')
    )
    process.patTauDiscriminationByMediumElectronRejectionPhase2MVA6v1.mapping.append(pset)
# Tight
WP='Eff70'
process.patTauDiscriminationByTightElectronRejectionPhase2MVA6v1 = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.PATTauProducer,
    Prediscriminants = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.Prediscriminants,
    toMultiplex = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw"),
    key = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw","category"),
    mapping = cms.VPSet()
)
for category, gbrForestName in antiElectronDiscrMVA_phase2_categories.items():
    pset = cms.PSet(
        category = cms.uint32(int(category)),
        cut = cms.string('RecoTauTag_antiElectronPhase2MVA6%s_%s_WP%s' % (antiElectronDiscrMVA_phase2_version, gbrForestName, WP)),
        variable = cms.string('pt')
    )
    process.patTauDiscriminationByTightElectronRejectionPhase2MVA6v1.mapping.append(pset)
# Tight
WP='Eff60'
process.patTauDiscriminationByVTightElectronRejectionPhase2MVA6v1 = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.PATTauProducer,
    Prediscriminants = process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw.Prediscriminants,
    toMultiplex = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw"),
    key = cms.InputTag("patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw","category"),
    mapping = cms.VPSet()
)
for category, gbrForestName in antiElectronDiscrMVA_phase2_categories.items():
    pset = cms.PSet(
        category = cms.uint32(int(category)),
        cut = cms.string('RecoTauTag_antiElectronPhase2MVA6%s_%s_WP%s' % (antiElectronDiscrMVA_phase2_version, gbrForestName, WP)),
        variable = cms.string('pt')
    )
    process.patTauDiscriminationByVTightElectronRejectionPhase2MVA6v1.mapping.append(pset)

process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Seq = cms.Sequence(
    process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw
    * process.patTauDiscriminationByVLooseElectronRejectionPhase2MVA6v1
    * process.patTauDiscriminationByLooseElectronRejectionPhase2MVA6v1
    * process.patTauDiscriminationByMediumElectronRejectionPhase2MVA6v1
    * process.patTauDiscriminationByTightElectronRejectionPhase2MVA6v1
    * process.patTauDiscriminationByVTightElectronRejectionPhase2MVA6v1
)

process.rerunAntiePhase2MVAnOnMiniAOD = cms.EDAnalyzer('rerunAntiePhase2MVAnOnMiniAOD'
)

process.rerunAntiePhase2MVAnOnMiniAOD.verbosity = cms.int32(0)
process.rerunAntiePhase2MVAnOnMiniAOD.additionalCollectionsAvailable = cms.bool(True)

# embed new id's into tau
embedID = cms.EDProducer("PATTauIDEmbedder",
   src = cms.InputTag('slimmedTaus'),
   tauIDSources = cms.PSet(
       againstElectronPhase2MVA6v1Raw = cms.InputTag('patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw'),
       againstElectronPhase2MVA6v1Category = cms.InputTag('patTauDiscriminationByElectronRejectionPhase2MVA6v1Raw:category'),
       againstElectronPhase2MVA6v1VLoose = cms.InputTag('patTauDiscriminationByVLooseElectronRejectionPhase2MVA6v1'),
       againstElectronPhase2MVA6v1Loose = cms.InputTag('patTauDiscriminationByLooseElectronRejectionPhase2MVA6v1'),
       againstElectronPhase2MVA6v1Medium = cms.InputTag('patTauDiscriminationByMediumElectronRejectionPhase2MVA6v1'),
       againstElectronPhase2MVA6v1Tight = cms.InputTag('patTauDiscriminationByTightElectronRejectionPhase2MVA6v1'),
       againstElectronPhase2MVA6v1VTight = cms.InputTag('patTauDiscriminationByVTightElectronRejectionPhase2MVA6v1')
   ),
)
setattr(process, "newTauIDsEmbedded", embedID)

## added for mvaIsolation on miniAOD testing
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple_newTauIDs.root'),
                               ## save only events passing the full path
                               #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               ## save PAT output; you need a '*' to unpack the list of commands
                               ## 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', "keep *_newTauIDsEmbedded_*_*")
                               )

process.p = cms.Path(
    process.patTauDiscriminationByElectronRejectionPhase2MVA6v1Seq
    * process.newTauIDsEmbedded
    * process.rerunAntiePhase2MVAnOnMiniAOD
)

process.outpath = cms.EndPath(process.out)
