/*
 * PATTauDiscriminationByLeadingObjectPtCut.cc
 *
 *  Created on: Nov 20, 2014
 *      Author: nehrkorn
 */

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminationByLeadingObjectPtCutT.h"

typedef RecoTauDiscriminationByLeadingObjectPtCutT<pat::TauRef, PATTauDiscriminationProducerBase> PATTauDiscriminationByLeadingObjectPtCut;

DEFINE_FWK_MODULE(PATTauDiscriminationByLeadingObjectPtCut);
