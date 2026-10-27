import FWCore.ParameterSet.Config as cms
from Validation.RecoTau.TauTruthGraphValidator_cfi import TauTruthGraphValidator as _TauTruthGraphValidator

# HLT
hltTauTruthGraphValidation = _TauTruthGraphValidator.clone(
    folder = cms.string("HLT/Tau/TruthGraphTauValidation"),
)

