import FWCore.ParameterSet.Config as cms
from Validation.RecoTau.TauTruthGraphValidator_cfi import TauTruthGraphValidator as _TauTruthGraphValidator

# HLT
hltTauTruthGraphValidation = _TauTruthGraphValidator.clone(
    folder = cms.string("HLT/Tau/TruthGraphTauValidation"),
    hpsTauToken = cms.InputTag("hltHpsPFTauProducer"),
    graphToken = cms.InputTag("truthLogicalGraphProducer"),
    ecalMapToken = cms.InputTag("hltTruthBranchPFClusterEcalAssociators", "hltParticleFlowClusterECALRecoToTruthFixed"),
    hcalMapToken = cms.InputTag("hltTruthBranchPFClusterHcalAssociators","hltParticleFlowClusterHCALRecoToTruthFixed"),
    pfCandidatesToken = cms.InputTag("hltParticleFlowTmp"),
    pfCandidatesBarrelToken = cms.InputTag("hltParticleFlowTmpBarrel"),
    ticlCandidatesToken = cms.InputTag("hltTiclCandidate"),
    hgcalTracksterMapToken = cms.InputTag("hltTruthBranchTracksterAssociators", "hltTiclCandidateRecoToTruthFixed"),
    # Keep clusters collection tags aligned with the HLT associator labels 
    # as in hltTruthBranchPFClusterEcalAssociators and hltTruthBranchPFClusterHcalAssociators.
    ecalClustersToken = cms.InputTag("hltParticleFlowClusterECAL"),
    hcalClustersToken = cms.InputTag("hltParticleFlowClusterHCAL"),
    hgcalTrackstersToken = cms.InputTag("hltTiclCandidate"),
)