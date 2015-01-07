/*
 * PATTauDiscriminationByIsolation.cc
 *
 *  Created on: Oct 7, 2014
 *      Author: nehrkorn
 */

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminationByIsolationT.h"

typedef RecoTauDiscriminationByIsolationT<pat::TauRef, reco::PFCandidateCollection, reco::PFCandidatePtr, PATTauDiscriminationProducerBase> PATTauDiscriminationByIsolation;

DEFINE_FWK_MODULE(PATTauDiscriminationByIsolation);
