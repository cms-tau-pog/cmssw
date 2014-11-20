/*
 * PFTauDiscriminationByLeadingObjectPtCut.cc
 *
 *  Created on: Nov 20, 2014
 *      Author: nehrkorn
 */

#include "DataFormats/TauReco/interface/PFTau.h"
#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminationByLeadingObjectPtCutT.h"

typedef RecoTauDiscriminationByLeadingObjectPtCutT<reco::PFTauRef, PFTauDiscriminationProducerBase> PFTauDiscriminationByLeadingObjectPtCut;

DEFINE_FWK_MODULE(PFTauDiscriminationByLeadingObjectPtCut);
