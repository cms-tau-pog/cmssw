/*
 * PFTauDiscriminationByIsolation.cc
 *
 *  Created on: Oct 7, 2014
 *      Author: nehrkorn
 */

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminationByIsolationT.h"

typedef RecoTauDiscriminationByIsolationT<reco::PFTauRef,reco::PFCandidateCollection,reco::PFCandidatePtr, PFTauDiscriminationProducerBase> PFTauDiscriminationByIsolation;

DEFINE_FWK_MODULE(PFTauDiscriminationByIsolation);
