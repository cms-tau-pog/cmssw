import FWCore.ParameterSet.Config as cms

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
from PhysicsTools.PatAlgos.PATTauDiscriminantCutMultiplexer_cfi import *

discriminationByIsolationMVArun2v1raw = cms.EDProducer("PATTauDiscriminationByMVAIsolationRun2",

    # tau collection to discriminate
    PATTauProducer = cms.InputTag('pfTauProducer'),
    Prediscriminants = noPrediscriminants,
    loadMVAfromDB = cms.bool(True),
    mvaName = cms.string("tauIdMVAnewDMwLT"),
    mvaOpt = cms.string("newDMwLT"),
    requireDecayMode = cms.bool(True),
    
    srcChargedIsoPtSum = cms.string('chargedIsoPtSum'),
    srcNeutralIsoPtSum = cms.string('neutralIsoPtSum'),
    srcPUcorrPtSum = cms.string('puCorrPtSum'),
    srcPhotonPtSumOutsideSignalCone = cms.string('photonPtSumOutsideSignalCone'),
    srcFootprintCorrection = cms.string('footprintCorrection') 
)

discriminationByIsolationMVArun2v1VLoose = patTauDiscriminantCutMultiplexer.clone(
    PATTauProducer = cms.InputTag('pfTauProducer'),    
    Prediscriminants = noPrediscriminants,
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
