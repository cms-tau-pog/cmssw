/*
 * RecoTauDiscriminationByIsolationT.h
 *
 *  Created on: Oct 1, 2014
 *      Author: nehrkorn
 */

#ifndef RECOTAUDISCRIMINATIONBYISOLATIONT_H_
#define RECOTAUDISCRIMINATIONBYISOLATIONT_H_

#include <functional>
#include <boost/foreach.hpp>
#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"
#include "RecoTauTag/RecoTau/interface/ConeTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TMath.h"
#include "TFormula.h"

using namespace reco;
using namespace std;

namespace candidateExtractor_namespace{
	template<typename Ttau, typename TcandPtr>
	class candidateExtractorT
	{
	public:
		candidateExtractorT(){}
		TcandPtr leadChargedHadr(const Ttau& tau) const { return tau->leadPFChargedHadrCand(); }
		TcandPtr leadNeutral(const Ttau& tau) const	{ return tau->leadPFNeutralCand(); }
		TcandPtr lead(const Ttau& tau) const { return tau->leadPFCand(); }
		std::vector<TcandPtr> isoChargedHadr(const Ttau& tau) const { return tau->isolationPFChargedHadrCands(); }
		std::vector<TcandPtr> isoNeutralHadr(const Ttau& tau) const { return tau->isolationPFNeutrHadrCands(); }
		std::vector<TcandPtr> isoGamma(const Ttau& tau) const { return tau->isolationPFGammaCands(); }
		std::vector<TcandPtr> signalChargedHadr(const Ttau& tau) const { return tau->signalPFChargedHadrCands(); }
	};

	template<>
	class candidateExtractorT<pat::TauRef, pat::PackedCandidatePtr>
	{
	public:
		candidateExtractorT(){}
		pat::PackedCandidatePtr leadChargedHadr(const pat::TauRef& tau) const { return pat::PackedCandidatePtr(tau->leadChargedHadrCand().id(),dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get()),tau->leadChargedHadrCand().key()); }
		pat::PackedCandidatePtr leadNeutral(const pat::TauRef& tau) const { return pat::PackedCandidatePtr(tau->leadNeutralCand().id(),dynamic_cast<const pat::PackedCandidate*>(tau->leadNeutralCand().get()),tau->leadNeutralCand().key()); }
		pat::PackedCandidatePtr lead(const pat::TauRef& tau) const { return pat::PackedCandidatePtr(tau->leadCand().id(),dynamic_cast<const pat::PackedCandidate*>(tau->leadCand().get()),tau->leadCand().key()); }
		std::vector<pat::PackedCandidatePtr> isoChargedHadr(const pat::TauRef& tau) const {
			std::vector<pat::PackedCandidatePtr> ret;
			ret.resize(tau->isolationChargedHadrCands().size());
			if(!tau->isPFTau()){
				for(size_t iICH = 0; iICH < tau->isolationChargedHadrCands().size(); ++iICH){
					ret.at(iICH) = pat::PackedCandidatePtr(tau->isolationChargedHadrCands()[iICH].id(),dynamic_cast<const pat::PackedCandidate*>(tau->isolationChargedHadrCands()[iICH].get()),tau->isolationChargedHadrCands()[iICH].key());
				}
			}
			return ret;
		}
		std::vector<pat::PackedCandidatePtr> isoNeutralHadr(const pat::TauRef& tau) const {
			std::vector<pat::PackedCandidatePtr> ret;
			ret.resize(tau->isolationNeutrHadrCands().size());
			if(!tau->isPFTau()){
				for(size_t iICH = 0; iICH < tau->isolationNeutrHadrCands().size(); ++iICH){
					ret.at(iICH) = pat::PackedCandidatePtr(tau->isolationNeutrHadrCands()[iICH].id(),dynamic_cast<const pat::PackedCandidate*>(tau->isolationNeutrHadrCands()[iICH].get()),tau->isolationNeutrHadrCands()[iICH].key());
				}
			}
			return ret;
		}
		std::vector<pat::PackedCandidatePtr> isoGamma(const pat::TauRef& tau) const {
			std::vector<pat::PackedCandidatePtr> ret;
			ret.resize(tau->isolationGammaCands().size());
			if(!tau->isPFTau()){
				for(size_t iICH = 0; iICH < tau->isolationGammaCands().size(); ++iICH){
					ret.at(iICH) = pat::PackedCandidatePtr(tau->isolationGammaCands()[iICH].id(),dynamic_cast<const pat::PackedCandidate*>(tau->isolationGammaCands()[iICH].get()),tau->isolationGammaCands()[iICH].key());
				}
			}
			return ret;
		}
		std::vector<pat::PackedCandidatePtr> signalChargedHadr(const pat::TauRef& tau) const {
			std::vector<pat::PackedCandidatePtr> ret;
			ret.resize(tau->signalChargedHadrCands().size());
			if(!tau->isPFTau()){
				for(size_t iICH = 0; iICH < tau->signalChargedHadrCands().size(); ++iICH){
					ret.at(iICH) = pat::PackedCandidatePtr(tau->signalChargedHadrCands()[iICH].id(),dynamic_cast<const pat::PackedCandidate*>(tau->signalChargedHadrCands()[iICH].get()),tau->signalChargedHadrCands()[iICH].key());
				}
			}
			return ret;
		}
	};
};

template<typename Ttau, typename TcandColl, typename TcandPtr, typename Tdiscr>
class RecoTauDiscriminationByIsolationT : public Tdiscr
{
 public:
  explicit RecoTauDiscriminationByIsolationT(const edm::ParameterSet& pset)
    : Tdiscr(pset),
      moduleLabel_(pset.getParameter<std::string>("@module_label")),
      qualityCutsPSet_(pset.getParameter<edm::ParameterSet>("qualityCuts"))
  {
    includeTracks_ = pset.getParameter<bool>(
      "ApplyDiscriminationByTrackerIsolation");
    includeGammas_ = pset.getParameter<bool>(
      "ApplyDiscriminationByECALIsolation");

    applyOccupancyCut_ = pset.getParameter<bool>("applyOccupancyCut");
    maximumOccupancy_ = pset.getParameter<uint32_t>("maximumOccupancy");

    applySumPtCut_ = pset.getParameter<bool>("applySumPtCut");
    maximumSumPt_ = pset.getParameter<double>("maximumSumPtCut");

    applyRelativeSumPtCut_ = pset.getParameter<bool>(
      "applyRelativeSumPtCut");
    maximumRelativeSumPt_ = pset.getParameter<double>(
      "relativeSumPtCut");

    storeRawOccupancy_ = pset.exists("storeRawOccupancy") ?
      pset.getParameter<bool>("storeRawOccupancy") : false;
    storeRawSumPt_ = pset.exists("storeRawSumPt") ?
      pset.getParameter<bool>("storeRawSumPt") : false;
    storeRawPUsumPt_ = pset.exists("storeRawPUsumPt") ?
      pset.getParameter<bool>("storeRawPUsumPt") : false;

    // Sanity check on requested options.  We can't apply cuts and store the
    // raw output at the same time
    if ( applySumPtCut_ || applyOccupancyCut_ || applyRelativeSumPtCut_ ) {
      if ( storeRawSumPt_ || storeRawOccupancy_ || storeRawPUsumPt_ ) {
	throw cms::Exception("BadIsoConfig")
	  << "A 'store raw' and a 'apply cut' option have been set to true "
	  << "simultaneously.  These options are mutually exclusive.";
      }
    }

    // Can only store one type
    int numStoreOptions = 0;
    if ( storeRawSumPt_     ) ++numStoreOptions;
    if ( storeRawOccupancy_ ) ++numStoreOptions;
    if ( storeRawPUsumPt_   ) ++numStoreOptions;
    if ( numStoreOptions > 1 ) {
      throw cms::Exception("BadIsoConfig")
	<< "Both 'store sum pt' and 'store occupancy' options are set."
	<< " These options are mutually exclusive.";
    }

    if ( pset.exists("customOuterCone") ) {
      customIsoCone_ = pset.getParameter<double>("customOuterCone");
    } else {
      customIsoCone_ = -1;
    }

    // Get the quality cuts specific to the isolation region
    edm::ParameterSet isolationQCuts = qualityCutsPSet_.getParameterSet(
      "isolationQualityCuts");

    qcuts_.reset(new tau::RecoTauQualityCuts(isolationQCuts));

    vertexAssociator_.reset(
			    new tau::RecoTauVertexAssociator(qualityCutsPSet_,edm::EDProducer::consumesCollector()));

    applyDeltaBeta_ = pset.exists("applyDeltaBetaCorrection") ?
      pset.getParameter<bool>("applyDeltaBetaCorrection") : false;

    if ( applyDeltaBeta_ ) {
      // Factorize the isolation QCuts into those that are used to
      // select PU and those that are not.
      std::pair<edm::ParameterSet, edm::ParameterSet> puFactorizedIsoQCuts =
	reco::tau::factorizePUQCuts(isolationQCuts);

      // Determine the pt threshold for the PU tracks
      // First check if the user specifies explicitly the cut.
      if ( pset.exists("deltaBetaPUTrackPtCutOverride") ) {
	puFactorizedIsoQCuts.second.addParameter<double>(
	  "minTrackPt",
	  pset.getParameter<double>("deltaBetaPUTrackPtCutOverride"));
      } else {
	// Secondly take it from the minGammaEt
	puFactorizedIsoQCuts.second.addParameter<double>(
          "minTrackPt",
	  isolationQCuts.getParameter<double>("minGammaEt"));
      }

      pileupQcutsPUTrackSelection_.reset(new tau::RecoTauQualityCuts(
        puFactorizedIsoQCuts.first));

      pileupQcutsGeneralQCuts_.reset(new tau::RecoTauQualityCuts(
        puFactorizedIsoQCuts.second));

      pfCandSrc_ = pset.getParameter<edm::InputTag>("particleFlowSrc");
      pfCand_token=edm::EDProducer::consumes<TcandColl>(pfCandSrc_);
      vertexSrc_ = pset.getParameter<edm::InputTag>("vertexSrc");
      vertex_token=edm::EDProducer::consumes<reco::VertexCollection>(vertexSrc_);
      deltaBetaCollectionCone_ = pset.getParameter<double>(
        "isoConeSizeForDeltaBeta");
      std::string deltaBetaFactorFormula =
	pset.getParameter<string>("deltaBetaFactor");
      deltaBetaFormula_.reset(
        new TFormula("DB_corr", deltaBetaFactorFormula.c_str()));
    }

    applyRhoCorrection_ = pset.exists("applyRhoCorrection") ?
      pset.getParameter<bool>("applyRhoCorrection") : false;
    if ( applyRhoCorrection_ ) {
      rhoProducer_ = pset.getParameter<edm::InputTag>("rhoProducer");
      rho_token=edm::EDProducer::consumes<double>(rhoProducer_);
      rhoConeSize_ = pset.getParameter<double>("rhoConeSize");
      rhoUEOffsetCorrection_ =
	pset.getParameter<double>("rhoUEOffsetCorrection");
    }

    verbosity_ = ( pset.exists("verbosity") ) ?
      pset.getParameter<int>("verbosity") : 0;
  }

  ~RecoTauDiscriminationByIsolationT(){}

  void beginEvent(const edm::Event& evt, const edm::EventSetup& evtSetup) override;
  double discriminate(const Ttau& Tau) override;

 private:
  std::string moduleLabel_;

  edm::ParameterSet qualityCutsPSet_;
  std::auto_ptr<tau::RecoTauQualityCuts> qcuts_;

  // Inverted QCut which selects tracks with bad DZ/trackWeight
  std::auto_ptr<tau::RecoTauQualityCuts> pileupQcutsPUTrackSelection_;
  std::auto_ptr<tau::RecoTauQualityCuts> pileupQcutsGeneralQCuts_;

  std::auto_ptr<tau::RecoTauVertexAssociator> vertexAssociator_;

  bool includeTracks_;
  bool includeGammas_;
  bool applyOccupancyCut_;
  uint32_t maximumOccupancy_;
  bool applySumPtCut_;
  double maximumSumPt_;
  bool applyRelativeSumPtCut_;
  double maximumRelativeSumPt_;
  double customIsoCone_;

  // Options to store the raw value in the discriminator instead of boolean pass/fail flag
  bool storeRawOccupancy_;
  bool storeRawSumPt_;
  bool storeRawPUsumPt_;

  /* **********************************************************************
     **** Pileup Subtraction Parameters ***********************************
     **********************************************************************/

  // Delta Beta correction
  bool applyDeltaBeta_;
  edm::InputTag pfCandSrc_;
  edm::EDGetTokenT<TcandColl> pfCand_token;
  // Keep track of how many vertices are in the event
  edm::InputTag vertexSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vertex_token;
  std::vector<TcandPtr> chargedPFCandidatesInEvent_;
  std::vector<TcandPtr> isoCharged_;
  std::vector<TcandPtr> isoNeutral_;
  std::vector<TcandPtr> isoPU_;

  // Size of cone used to collect PU tracks
  double deltaBetaCollectionCone_;
  std::auto_ptr<TFormula> deltaBetaFormula_;
  double deltaBetaFactorThisEvent_;

  // Rho correction
  bool applyRhoCorrection_;
  edm::InputTag rhoProducer_;
  edm::EDGetTokenT<double> rho_token;
  double rhoConeSize_;
  double rhoUEOffsetCorrection_;
  double rhoCorrectionThisEvent_;
  double rhoThisEvent_;

  // Flag to enable/disable debug output
  int verbosity_;
};

template<typename Ttau, typename TcandColl, typename TcandPtr, typename Tdiscr>
void RecoTauDiscriminationByIsolationT<Ttau, TcandColl, TcandPtr, Tdiscr>::beginEvent(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  // NB: The use of the PV in this context is necessitated by its use in
  // applying quality cuts to the different objects in the isolation cone
  // The vertex associator contains the logic to select the appropriate vertex
  // We need to pass it the event so it can load the vertices.
  vertexAssociator_->setEvent(event);

  // If we are applying the delta beta correction, we need to get the PF
  // candidates from the event so we can find the PU tracks.
  if ( applyDeltaBeta_ ) {
    // Collect all the PF pile up tracks
    edm::Handle<TcandColl> pfCandidates;
    event.getByToken(pfCand_token, pfCandidates);
    chargedPFCandidatesInEvent_.clear();
    chargedPFCandidatesInEvent_.reserve(pfCandidates->size());
    size_t numPFCandidates = pfCandidates->size();
    for ( size_t i = 0; i < numPFCandidates; ++i ) {
      TcandPtr pfCandidate(pfCandidates, i);
      if ( pfCandidate->charge() != 0 ) {
        chargedPFCandidatesInEvent_.push_back(pfCandidate);
      }
    }
    // Count all the vertices in the event, to parameterize the DB
    // correction factor
    edm::Handle<reco::VertexCollection> vertices;
    event.getByToken(vertex_token, vertices);
    size_t nVtxThisEvent = vertices->size();
    deltaBetaFactorThisEvent_ = deltaBetaFormula_->Eval(nVtxThisEvent);
  }

  if ( applyRhoCorrection_ ) {
    edm::Handle<double> rhoHandle_;
    event.getByToken(rho_token, rhoHandle_);
    rhoThisEvent_ = (*rhoHandle_ - rhoUEOffsetCorrection_)*
      (3.14159)*rhoConeSize_*rhoConeSize_;
  }
}

template<typename Ttau, typename TcandColl, typename TcandPtr, typename Tdiscr>
double
RecoTauDiscriminationByIsolationT<Ttau, TcandColl, TcandPtr, Tdiscr>::discriminate(const Ttau& pfTau)
{
  if ( verbosity_ ) {
    std::cout << "<PFRecoTauDiscriminationByIsolation::discriminate (moduleLabel = " << moduleLabel_ <<")>:" << std::endl;
    std::cout << " tau: Pt = " << pfTau->pt() << ", eta = " << pfTau->eta() << ", phi = " << pfTau->phi() << std::endl;
  }

  candidateExtractor_namespace::candidateExtractorT<Ttau, TcandPtr> candidates;

  // collect the objects we are working with (ie tracks, tracks+gammas, etc)
  isoCharged_.clear();
  isoCharged_.reserve(candidates.isoChargedHadr(pfTau).size());
  isoNeutral_.clear();
  isoNeutral_.reserve(candidates.isoNeutralHadr(pfTau).size());
  isoPU_.clear();
  isoPU_.reserve(chargedPFCandidatesInEvent_.size());

  // Get the primary vertex associated to this tau
  reco::VertexRef pv = vertexAssociator_->associatedVertex(*pfTau);
  // Let the quality cuts know which the vertex to use when applying selections
  // on dz, etc.
  if ( verbosity_ ) {
    if ( pv.isNonnull() ) {
      std::cout << "pv: x = " << pv->position().x() << ", y = " << pv->position().y() << ", z = " << pv->position().z() << std::endl;
    } else {
      std::cout << "pv: N/A" << std::endl;
    }
    if ( candidates.leadChargedHadr(pfTau).isNonnull() ) {
      std::cout << "leadPFChargedHadron:"
    	<< " Pt = " << candidates.leadChargedHadr(pfTau)->pt() << ","
    	<< " eta = " << candidates.leadChargedHadr(pfTau)->eta() << ","
    	<< " phi = " << candidates.leadChargedHadr(pfTau)->phi() << std::endl;
    } else {
      std::cout << "leadPFChargedHadron: N/A" << std::endl;
    }
  }
  // CV: isolation is not well defined in case primary vertex or leading charged hadron do not exist
  if ( !(pv.isNonnull() && candidates.leadChargedHadr(pfTau).isNonnull()) ) return 0.;

  qcuts_->setPV(pv);
  qcuts_->setLeadTrack(candidates.leadChargedHadr(pfTau));

  if ( applyDeltaBeta_ ) {
    pileupQcutsGeneralQCuts_->setPV(pv);
    pileupQcutsGeneralQCuts_->setLeadTrack(candidates.leadChargedHadr(pfTau));
    pileupQcutsPUTrackSelection_->setPV(pv);
    pileupQcutsPUTrackSelection_->setLeadTrack(candidates.leadChargedHadr(pfTau));
  }

  // Load the tracks if they are being used.
  if ( includeTracks_ ) {
	BOOST_FOREACH( const TcandPtr& cand, (const std::vector<TcandPtr>&)candidates.isoChargedHadr(pfTau) ) {
      if ( qcuts_->filterCandRef(cand) ) {
        isoCharged_.push_back(cand);
      }
    }
  }
  if ( includeGammas_ ) {
    BOOST_FOREACH( const TcandPtr& cand, (const std::vector<TcandPtr>&)candidates.isoGamma(pfTau) ) {
      if ( qcuts_->filterCandRef(cand) ) {
        isoNeutral_.push_back(cand);
      }
    }
  }

  typedef reco::tau::cone::DeltaRPtrFilter<TcandPtr> DRFilter;

  // If desired, get PU tracks.
  if ( applyDeltaBeta_ ) {
    // First select by inverted the DZ/track weight cuts. True = invert
    //if ( verbosity_ ) {
    //  std::cout << "Initial PFCands: " << chargedPFCandidatesInEvent_.size() << std::endl;
    //}

    std::vector<TcandPtr> allPU =
      pileupQcutsPUTrackSelection_->filterCandRefs(
          chargedPFCandidatesInEvent_, true);
    //if ( verbosity_ ) {
    //  std::cout << "After track cuts: " << allPU.size() << std::endl;
    //}

    // Now apply the rest of the cuts, like pt, and TIP, tracker hits, etc
    std::vector<TcandPtr> cleanPU =
      pileupQcutsGeneralQCuts_->filterCandRefs(allPU);
    //if ( verbosity_ ) {
    //  std::cout << "After cleaning cuts: " << cleanPU.size() << std::endl;
    //}

    // Only select PU tracks inside the isolation cone.
    DRFilter deltaBetaFilter(pfTau->p4(), 0, deltaBetaCollectionCone_);
    BOOST_FOREACH(const TcandPtr& cand, cleanPU) {
      if ( deltaBetaFilter(cand) ) {
        isoPU_.push_back(cand);
      }
    }
    //if ( verbosity_ ) {
    //  std::cout << "After cone cuts: " << isoPU.size() << std::endl;
    //}
  }

  // Check if we want a custom iso cone
  if ( customIsoCone_ >= 0. ) {
    //std::cout << "<PFRecoTauDiscriminationByIsolation::discriminate (moduleLabel = " << moduleLabel_ <<")>:" << std::endl;
    //std::cout << " customIsoCone = " << customIsoCone_ << std::endl;
    DRFilter filter(pfTau->p4(), 0, customIsoCone_);
    std::vector<TcandPtr> isoCharged_filter;
    std::vector<TcandPtr> isoNeutral_filter;
    // Remove all the objects not in our iso cone
    BOOST_FOREACH( const TcandPtr& isoObject, isoCharged_ ) {
      if ( filter(isoObject) ) isoCharged_filter.push_back(isoObject);
    }
    BOOST_FOREACH( const TcandPtr& isoObject, isoNeutral_ ) {
      if ( filter(isoObject) ) isoNeutral_filter.push_back(isoObject);
    }
    isoCharged_ = isoCharged_filter;
    isoNeutral_ = isoNeutral_filter;
  }

  bool failsOccupancyCut     = false;
  bool failsSumPtCut         = false;
  bool failsRelativeSumPtCut = false;

//--- nObjects requirement
  int neutrals = isoNeutral_.size();

  if ( applyDeltaBeta_ ) {
    neutrals -= TMath::Nint(deltaBetaFactorThisEvent_*isoPU_.size());
  }
  if ( neutrals < 0 ) {
    neutrals = 0;
  }

  size_t nOccupants = isoCharged_.size() + neutrals;

  failsOccupancyCut = ( nOccupants > maximumOccupancy_ );

  double totalPt = 0.0;
  double puPt = 0.0;
//--- Sum PT requirement
  if ( applySumPtCut_ || applyRelativeSumPtCut_ || storeRawSumPt_ || storeRawPUsumPt_ ) {
    double chargedPt = 0.0;
    double neutralPt = 0.0;
    BOOST_FOREACH ( const TcandPtr& isoObject, isoCharged_ ) {
      chargedPt += isoObject->pt();
    }
    BOOST_FOREACH ( const TcandPtr& isoObject, isoNeutral_ ) {
      neutralPt += isoObject->pt();
    }
    BOOST_FOREACH ( const TcandPtr& isoObject, isoPU_ ) {
      puPt += isoObject->pt();
    }
    if ( verbosity_ ) {
      std::cout << "chargedPt = " << chargedPt << std::endl;
      std::cout << "neutralPt = " << neutralPt << std::endl;
      std::cout << "puPt = " << puPt << " (delta-beta corr. = " << (deltaBetaFactorThisEvent_*puPt) << ")" << std::endl;
    }

    if ( applyDeltaBeta_ ) {
      neutralPt -= deltaBetaFactorThisEvent_*puPt;
    }

    if ( applyRhoCorrection_ ) {
      neutralPt -= rhoThisEvent_;
    }

    if ( neutralPt < 0.0 ) {
      neutralPt = 0.0;
    }

    totalPt = chargedPt + neutralPt;
    if ( verbosity_ ) {
      std::cout << "totalPt = " << totalPt << " (cut = " << maximumSumPt_ << ")" << std::endl;
    }

    failsSumPtCut = (totalPt > maximumSumPt_);

//--- Relative Sum PT requirement
    failsRelativeSumPtCut = (totalPt > (pfTau->pt()*maximumRelativeSumPt_));
  }

  bool fails = (applyOccupancyCut_ && failsOccupancyCut) ||
    (applySumPtCut_ && failsSumPtCut) ||
    (applyRelativeSumPtCut_ && failsRelativeSumPtCut);

  // We did error checking in the constructor, so this is safe.
  if ( storeRawSumPt_ ) {
    return totalPt;
  } else if ( storeRawPUsumPt_ ) {
    if ( applyDeltaBeta_ ) return puPt;
    else if ( applyRhoCorrection_ ) return rhoThisEvent_;
    else return 0.;
  } else if ( storeRawOccupancy_ ) {
    return nOccupants;
  } else {
    return (fails ? 0. : 1.);
  }
}

#endif /* RECOTAUDISCRIMINATIONBYISOLATIONT_H_ */
