/*
 * RecoTauDiscriminationByLeadingObjectPtCutT.h
 *
 *  Created on: Nov 20, 2014
 *      Author: nehrkorn
 */

#ifndef RECOTAUDISCRIMINATIONBYLEADINGOBJECTPTCUTT_H_
#define RECOTAUDISCRIMINATIONBYLEADINGOBJECTPTCUTT_H_

#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"

using namespace reco;

template<typename Ttau, typename Tdiscr>
class RecoTauDiscriminationByLeadingObjectPtCutT : public Tdiscr
{
public:
	explicit RecoTauDiscriminationByLeadingObjectPtCutT(const edm::ParameterSet& iConfig):Tdiscr(iConfig){
		chargedOnly_     = iConfig.getParameter<bool>("UseOnlyChargedHadrons");
		minPtLeadObject_ = iConfig.getParameter<double>("MinPtLeadingObject");
	}
	~RecoTauDiscriminationByLeadingObjectPtCutT(){}
	double discriminate(const Ttau& tau) override;
private:
	bool chargedOnly_;
	double minPtLeadObject_;
};

template<typename Ttau, typename Tdiscr>
double RecoTauDiscriminationByLeadingObjectPtCutT<Ttau, Tdiscr>::discriminate(const Ttau& tau)
{
	double leadObjectPt = -1.;
	if( chargedOnly_ )
	{
		// consider only charged hadrons.  note that the leadPFChargedHadrCand is the highest pt
		// charged signal cone object above the quality cut level (typically 0.5 GeV).
	    if( tau->leadPFChargedHadrCand().isNonnull() )
	    {
	    	leadObjectPt = tau->leadPFChargedHadrCand()->pt();
	    }
	}
	else
	{
		// If using the 'leading pion' option, require that:
		//   1) at least one charged hadron exists above threshold (thePFTauRef->leadPFChargedHadrCand().isNonnull())
		//   2) the lead PFCand exists.  In the case that the highest pt charged hadron is above the PFRecoTauProducer threshold
		//      (typically 5 GeV), the leadPFCand and the leadPFChargedHadrCand are the same object.  If the leadPFChargedHadrCand
		//      is below 5GeV, but there exists a neutral PF particle > 5 GeV, it is set to be the leadPFCand
		if( tau->leadPFCand().isNonnull() && tau->leadPFChargedHadrCand().isNonnull() )
		{
			leadObjectPt = tau->leadPFCand()->pt();
		}
	}

	return ( leadObjectPt > minPtLeadObject_ ? 1. : 0. );
}

#endif /* RECOTAUDISCRIMINATIONBYLEADINGOBJECTPTCUTT_H_ */
