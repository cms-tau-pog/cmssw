import FWCore.ParameterSet.Config as cms

process = cms.Process("lowPtTauCollection")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
#      'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_miniaod.root'
#        'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part0.root'
'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part0.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part1.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part2.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part3.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part4.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part5.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part6.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part7.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part8.root',
#'file:/afs/cern.ch/user/k/klau/myWorkspace/public/ForYuta/UpsilonToTauTau_3prong/UpsilonToTauTau_3prong_miniaod_part9.root'


    )
)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff') 
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') 
from Configuration.AlCa.GlobalTag import GlobalTag 
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v18')

#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v12')


process.lowPtTaus = cms.EDProducer('lowPtTaus',
                                   vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   packedpfcandidates = cms.InputTag('packedPFCandidates'),
                                   genparticles = cms.InputTag('prunedGenParticles'),
                                   abscharge = cms.int32(3),
                                   isMC = cms.bool(True)
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:/tmp/ytakahas/myOutputFile.root'),
#    outputCommands = cms.untracked.vstring('drop *',
#      "keep *_generalTracks_*_*",
#      "keep *_globalMuons_*_*",
#       "keep *_MuonTrackPoints_*_*",
#      "keep *_TrackTrackPoints_*_*")
       
)

  
process.p = cms.Path(process.lowPtTaus)

process.e = cms.EndPath(process.out)
