/*
 * PATSlimmedTauDiscriminationByIsolation.cc
 *
 *  Created on: Dec 8, 2014
 *      Author: nehrkorn
 */

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminationByIsolationT.h"

typedef RecoTauDiscriminationByIsolationT<pat::TauRef, pat::PackedCandidateCollection, pat::PackedCandidatePtr, PATTauDiscriminationProducerBase> PATSlimmedTauDiscriminationByIsolation;

DEFINE_FWK_MODULE(PATSlimmedTauDiscriminationByIsolation);
