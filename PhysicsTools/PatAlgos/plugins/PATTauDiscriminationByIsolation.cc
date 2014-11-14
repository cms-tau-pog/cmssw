/*
 * PATTauDiscriminationByIsolation.cc
 *
 *  Created on: Oct 7, 2014
 *      Author: nehrkorn
 */

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminationByIsolationT.h"

typedef RecoTauDiscriminationByIsolationT<pat::TauRef, pat::PackedCandidateCollection, edm::Ptr<pat::PackedCandidate>, PATTauDiscriminationProducerBase> PATTauDiscriminationByIsolation;

DEFINE_FWK_MODULE(PATTauDiscriminationByIsolation);
