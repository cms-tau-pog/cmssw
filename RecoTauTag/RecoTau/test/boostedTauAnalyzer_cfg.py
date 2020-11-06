import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17')


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISummer16MiniAODv3/ZprimeToZhToZhadhtata_narrow_M-2000_13TeV-madgraph/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/20000/E8530115-67F1-E911-A75F-0242AC130002.root'
    )
)
#################################################################
# A new boostedTau collection is made here and the overlap is removed
########################################################################################
from RecoTauTag.Configuration.boostedHPSPFTaus_cff import ca8PFJetsCHSprunedForBoostedTaus
process.ca8PFJetsCHSprunedForBoostedTausPAT = ca8PFJetsCHSprunedForBoostedTaus.clone(
    src=cms.InputTag("packedPFCandidates"),
    jetCollInstanceName = cms.string('subJetsForSeedingBoostedTausPAT')
)
cleanedBoostedTau = cms.EDProducer("PATBoostedTauCleaner",
   src = cms.InputTag('slimmedTausBoosted'),
   vtxLabel= cms.InputTag('offlineSlimmedPrimaryVertices'),
   ca8JetSrc = cms.InputTag('ca8PFJetsCHSprunedForBoostedTausPAT','subJetsForSeedingBoostedTausPAT'),
   removeOverLap = cms.bool(True),
   )
setattr(process, "cleanedSlimmedTausBoosted", cleanedBoostedTau)

########################################################################################
# Rerun boostedTau Id
########################################################################################
from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')

updatedBoostedTauName = "slimmedBoostedTausNewID" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runBoostedTauIdMVA as tauIdConfig
boostedTauIdEmbedder = tauIdConfig.BoostedTauIDEmbedder(
    process, cms, debug = False,
    updatedTauName = updatedBoostedTauName,
    PATTauProducer = cms.InputTag('cleanedSlimmedTausBoosted'),
    srcChargedIsoPtSum = cms.string('chargedIsoPtSumNoOverLap'),
    srcNeutralIsoPtSum = cms.string('neutralIsoPtSumNoOverLap'),
    toKeep = [
        #"2017v2", "dR0p32017v2", "newDM2017v2", #classic MVAIso tau-Ids
        #"deepTau2017v1", #deepTau Tau-Ids
        #"DPFTau_2016_v0", #D[eep]PF[low] Tau-Id
        "2017v2","deepTau2017v1","againstEle2018"
    ])
boostedTauIdEmbedder.runTauID()


########################################################################################
# Run ntuple Maker
########################################################################################

process.demo = cms.EDAnalyzer('BoostedTauAnalyzer',
    boostedTauSrc = cms.InputTag("slimmedTausBoosted"),
    cleanedBoostedTauSrc = cms.InputTag("slimmedBoostedTausNewID"),
 )

########################################################################################
# Output name
########################################################################################
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ntuple.root')
)
########################################################################################
# path
########################################################################################
#print process.dumpPython()
process.p = cms.Path(
    process.ca8PFJetsCHSprunedForBoostedTausPAT *
    getattr(process, "cleanedSlimmedTausBoosted") *
    process.rerunMvaIsolationBoostSequence *
    getattr(process, updatedBoostedTauName) *
    process.demo
)
