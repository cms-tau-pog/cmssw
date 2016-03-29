import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.PATTauDiscriminantCutMultiplexer_cfi import *

discriminationByIsolationMVArun2v1raw = cms.EDProducer("PATTauDiscriminationByMVAIsolationRun2",

    # tau collection to discriminate
    PATTauProducer = cms.InputTag('pfTauProducer'),

    # Require leading pion ensures that:
    #  1) these is at least one track above threshold (0.5 GeV) in the signal cone
    #  2) a track OR a pi-zero in the signal cone has pT > 5 GeV
    Prediscriminants = requireLeadTrack,
    loadMVAfromDB = cms.bool(True),
    mvaName = cms.string("tauIdMVAnewDMwLT"),
    mvaOpt = cms.string("newDMwLT"),
    
    srcChargedIsoPtSum = cms.string('chargedIsoPtSum'),
    srcNeutralIsoPtSum = cms.string('neutralIsoPtSum'),
    srcPUcorrPtSum = cms.string('puCorrPtSum'),
    srcPhotonPtSumOutsideSignalCone = cms.string('photonPtSumOutsideSignalCone'),
    srcFootprintCorrection = cms.string('footprintCorrection') 
)

discriminationByIsolationMVArun2v1VLoose = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = cms.InputTag('pfTauProducer'),    
    Prediscriminants = requireLeadTrack,
    toMultiplex = cms.InputTag('discriminationByIsolationMVArun2v1raw'),
    key = cms.InputTag('discriminationByIsolationMVArun2v1raw:category'),
    loadMVAfromDB = cms.bool(True),
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0),
            cut = cms.string("newDMwLTEff80"),
            variable = cms.string("pt"),
        )
    )
)
discriminationByIsolationMVArun2v1Loose = discriminationByIsolationMVArun2v1VLoose.clone()
discriminationByIsolationMVArun2v1Loose.mapping[0].cut = cms.string("newDMwLTEff70")
discriminationByIsolationMVArun2v1Medium = discriminationByIsolationMVArun2v1VLoose.clone()
discriminationByIsolationMVArun2v1Medium.mapping[0].cut = cms.string("newDMwLTEff60")
discriminationByIsolationMVArun2v1Tight = discriminationByIsolationMVArun2v1VLoose.clone()
discriminationByIsolationMVArun2v1Tight.mapping[0].cut = cms.string("newDMwLTEff50")
discriminationByIsolationMVArun2v1VTight = discriminationByIsolationMVArun2v1VLoose.clone()
discriminationByIsolationMVArun2v1VTight.mapping[0].cut = cms.string("newDMwLTEff40")

mvaIsolation2SeqRun2 = cms.Sequence(
   discriminationByIsolationMVArun2v1raw
   + discriminationByIsolationMVArun2v1VLoose
   + discriminationByIsolationMVArun2v1Loose
   + discriminationByIsolationMVArun2v1Medium
   + discriminationByIsolationMVArun2v1Tight
   + discriminationByIsolationMVArun2v1VTight
)