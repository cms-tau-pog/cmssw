import FWCore.ParameterSet.Config as cms
from Validation.RecoTau.TauTruthGraphValidator_cfi import TauTruthGraphValidator as _TauTruthGraphValidator

# HLT
hltTauTruthGraphValidation = _TauTruthGraphValidator.clone(
    folder = cms.string("HLT/Tau/TruthGraphTauValidation"),
    hpsTauToken = cms.InputTag("hltHpsPFTauProducer"),
    graphToken = cms.InputTag("truthLogicalGraphProducer"),
    ecalMapToken = cms.InputTag("hltTruthBranchPFClusterEcalAssociators", "hltParticleFlowClusterECALUnseededRecoToTruthFixed"),
    hcalMapToken = cms.InputTag("hltTruthBranchPFClusterHcalAssociators","hltParticleFlowClusterHCALRecoToTruthFixed"),
    # Keep clusters collection tags aligned with the HLT associator labels 
    # as in hltTruthBranchPFClusterEcalAssociators and hltTruthBranchPFClusterHcalAssociators.
    ecalClustersToken = cms.InputTag("hltParticleFlowClusterECALUnseeded"),
    hcalClustersToken = cms.InputTag("hltParticleFlowClusterHCAL"),
)