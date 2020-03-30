/* class PFRecoTauDiscriminationAgainstElectronMVA6
 * created : Nov 2 2015,
 * revised : ,
 * Authorss : Fabio Colombo (KIT)
 */

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>

#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "RecoTauTag/RecoTau/interface/AntiElectronIDMVA6.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoEgamma/EgammaTools/interface/HGCalEgammaIDHelper.h"

#include <iostream>
#include <sstream>
#include <fstream>

using namespace reco;

class PFRecoTauDiscriminationAgainstElectronMVA6 : public PFTauDiscriminationProducerBase {
public:
  explicit PFRecoTauDiscriminationAgainstElectronMVA6(const edm::ParameterSet& cfg)
      : PFTauDiscriminationProducerBase(cfg), mva_(), category_output_() {
    mva_ = std::make_unique<AntiElectronIDMVA6>(cfg);

    srcGsfElectrons_ = cfg.getParameter<edm::InputTag>("srcGsfElectrons");
    GsfElectrons_token = consumes<reco::GsfElectronCollection>(srcGsfElectrons_);
    vetoEcalCracks_ = cfg.getParameter<bool>("vetoEcalCracks");

    verbosity_ = cfg.getParameter<int>("verbosity");

    // add category index
    produces<PFTauDiscriminator>("category");
    
    eIDHelper_.reset(new HGCalEgammaIDHelper(cfg, consumesCollector()));
    radius_ = cfg.getParameter<double>("pcaRadius");
  }

  void beginEvent(const edm::Event&, const edm::EventSetup&) override;

  double discriminate(const PFTauRef&) const override;

  void endEvent(edm::Event&) override;

  ~PFRecoTauDiscriminationAgainstElectronMVA6() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool isInEcalCrack(double) const;

  std::string moduleLabel_;
  std::unique_ptr<AntiElectronIDMVA6> mva_;

  edm::InputTag srcGsfElectrons_;
  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectrons_token;
  edm::Handle<reco::GsfElectronCollection> gsfElectrons_;
  edm::Handle<TauCollection> taus_;

  std::unique_ptr<PFTauDiscriminator> category_output_;

  bool vetoEcalCracks_;

  int verbosity_;
  
  std::unique_ptr<HGCalEgammaIDHelper> eIDHelper_;
  float radius_;
};

void PFRecoTauDiscriminationAgainstElectronMVA6::beginEvent(const edm::Event& evt, const edm::EventSetup& es) {
  mva_->beginEvent(evt, es);

  evt.getByToken(Tau_token, taus_);
  category_output_.reset(new PFTauDiscriminator(TauRefProd(taus_)));

  evt.getByToken(GsfElectrons_token, gsfElectrons_);
  eIDHelper_->eventInit(evt,es);
}

double PFRecoTauDiscriminationAgainstElectronMVA6::discriminate(const PFTauRef& thePFTauRef) const {
  double mvaValue = 1.;
  double category = -1.;
  bool isGsfElectronMatched = false;

  float deltaRDummy = 9.9;

  const float ECALBarrelEndcapEtaBorder = 1.479;
  const float ECALEndcapEtaBorder = 2.4;
  
  float tauEtaAtEcalEntrance = -99.;
  float sumEtaTimesEnergy = 0.;
  float sumEnergy = 0.;
  for (const auto& pfCandidate : thePFTauRef->signalPFCands()) {
    sumEtaTimesEnergy += (pfCandidate->positionAtECALEntrance().eta() * pfCandidate->energy());
    sumEnergy += pfCandidate->energy();
  }
  if (sumEnergy > 0.) {
    tauEtaAtEcalEntrance = sumEtaTimesEnergy / sumEnergy;
  }

  float leadChargedPFCandEtaAtEcalEntrance = -99.;
  float leadChargedPFCandPt = -99.;
  for (const auto& pfCandidate : thePFTauRef->signalPFCands()) {
    const reco::Track* track = nullptr;
    if (pfCandidate->trackRef().isNonnull())
      track = pfCandidate->trackRef().get();
    else if (pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->innerTrack().isNonnull())
      track = pfCandidate->muonRef()->innerTrack().get();
    else if (pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->globalTrack().isNonnull())
      track = pfCandidate->muonRef()->globalTrack().get();
    else if (pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->outerTrack().isNonnull())
      track = pfCandidate->muonRef()->outerTrack().get();
    else if (pfCandidate->gsfTrackRef().isNonnull())
      track = pfCandidate->gsfTrackRef().get();
    if (track) {
      if (track->pt() > leadChargedPFCandPt) {
        leadChargedPFCandEtaAtEcalEntrance = pfCandidate->positionAtECALEntrance().eta();
        leadChargedPFCandPt = track->pt();
      }
    }
  }

  if ((*thePFTauRef).leadChargedHadrCand().isNonnull()) {
    int numSignalGammaCandsInSigCone = 0;
    const std::vector<reco::CandidatePtr>& signalGammaCands = thePFTauRef->signalGammaCands();

    for (const auto& pfGamma : signalGammaCands) {
      double dR = deltaR(pfGamma->p4(), thePFTauRef->leadChargedHadrCand()->p4());
      double signalrad = std::max(0.05, std::min(0.10, 3.0 / std::max(0.0, thePFTauRef->pt())));

      // pfGammas inside the tau signal cone
      if (dR < signalrad) {
        numSignalGammaCandsInSigCone += 1;
      }
    }
    

    // loop over the electrons
   
    const reco::GsfElectron* matchedGsfElectron = nullptr;

    for (const auto& theGsfElectron : *gsfElectrons_) {
    
       if (theGsfElectron.pt() > 10.) {  // CV: only take electrons above some minimal energy/Pt into account...
        double deltaREleTau = deltaR(theGsfElectron.p4(), thePFTauRef->p4());
        deltaRDummy = deltaREleTau;
	
	if ( deltaREleTau < 0.2 && 
	    (matchedGsfElectron == nullptr || (matchedGsfElectron != nullptr && theGsfElectron.pt() > (*matchedGsfElectron).pt())) ) {
	    matchedGsfElectron = &theGsfElectron;
	}
       }	  
     }
     
     if ( matchedGsfElectron != nullptr ) {
	  
	  //HGCAL stuff
	  float HGCAL_sigmaUU =      0.;
          float HGCAL_sigmaVV =      0.;
          float HGCAL_sigmaEE =      0.;
          float HGCAL_sigmaPP =      0.; 
          float HGCAL_nLayers =      0.;  
          float HGCAL_firstLayer =   0.;
          float HGCAL_lastLayer =    0.; 
          float HGCAL_layerEfrac10 = 0.;  
          float HGCAL_layerEfrac90 = 0.;
          float HGCAL_ecEnergyEE =   0.;
          float HGCAL_ecEnergyFH =   0.;
          float HGCAL_measuredDepth = 0.;	  
          float HGCAL_expectedDepth = 0.;	  
          float HGCAL_expectedSigma = 0.;	  
	  float HGCAL_depthCompatibility = 0.;
	  
	  if(!matchedGsfElectron->isEB()) {
	       eIDHelper_->computeHGCAL(*matchedGsfElectron, radius_);
	       hgcal::LongDeps ld(eIDHelper_->energyPerLayer(radius_, true));
               float HGCAL_measuredDepth, HGCAL_expectedDepth, HGCAL_expectedSigma;
               HGCAL_depthCompatibility = eIDHelper_->clusterDepthCompatibility(ld, HGCAL_measuredDepth, HGCAL_expectedDepth, HGCAL_expectedSigma);
	   
	       HGCAL_sigmaUU = eIDHelper_->sigmaUU(); 
	       HGCAL_sigmaVV = eIDHelper_->sigmaVV(); 
	       HGCAL_sigmaEE = eIDHelper_->sigmaEE();
	       HGCAL_sigmaPP = eIDHelper_->sigmaPP();
	       HGCAL_nLayers = ld.nLayers();
               HGCAL_firstLayer = ld.lastLayer();
    	       HGCAL_lastLayer = ld.lastLayer();
    	       HGCAL_layerEfrac10 = ld.layerEfrac10();
    	       HGCAL_layerEfrac90 = ld.layerEfrac90();
    	       HGCAL_ecEnergyEE = ld.energyEE();
    	       HGCAL_ecEnergyFH = ld.energyFH();
	  }
          
	  double mva_match = mva_->MVAValue(*thePFTauRef, *matchedGsfElectron, 
	                                    HGCAL_sigmaUU,
					    HGCAL_sigmaVV,
					    HGCAL_sigmaEE,
					    HGCAL_sigmaPP,
					    HGCAL_nLayers,
                                            HGCAL_firstLayer,
    					    HGCAL_lastLayer,
    					    HGCAL_layerEfrac10,
    					    HGCAL_layerEfrac90,
    					    HGCAL_ecEnergyEE,
    					    HGCAL_ecEnergyFH,
    					    HGCAL_measuredDepth,
    					    HGCAL_expectedDepth,
    					    HGCAL_expectedSigma,
					    HGCAL_depthCompatibility);

          const reco::PFCandidatePtr& lpfch = thePFTauRef->leadPFChargedHadrCand();
          bool hasGsfTrack = false;
          if (lpfch.isNonnull() && abs(lpfch->pdgId()) == 11) {
            //hasGsfTrack = lpfch->gsfTrackRef().isNonnull();
          hasGsfTrack = true;
	  }
          if (!hasGsfTrack)
            hasGsfTrack = (*matchedGsfElectron).gsfTrack().isNonnull();

          //// Veto taus that go to Ecal crack
          if (vetoEcalCracks_ &&
              (isInEcalCrack(tauEtaAtEcalEntrance) || isInEcalCrack(leadChargedPFCandEtaAtEcalEntrance))) {
            // add category index
            category_output_->setValue(tauIndex_, category);
            // return MVA output value
            return -99;
          }
          //// Veto taus that go to Ecal crack

          if (std::abs(tauEtaAtEcalEntrance) < ECALBarrelEndcapEtaBorder) {  // Barrel
            if (numSignalGammaCandsInSigCone == 0 && hasGsfTrack) {
              category = 5.;
            } else if (numSignalGammaCandsInSigCone >= 1 && hasGsfTrack) {
              category = 7.;
            }
          } else if ( std::abs(tauEtaAtEcalEntrance) > ECALBarrelEndcapEtaBorder && 
	              std::abs(tauEtaAtEcalEntrance) < ECALEndcapEtaBorder) {  // Endcap
            if (numSignalGammaCandsInSigCone == 0 && hasGsfTrack) {
              category = 13.;
            } else if (numSignalGammaCandsInSigCone >= 1 && hasGsfTrack) {
              category = 15.;
            }
          }  else if ( std::abs(tauEtaAtEcalEntrance) > ECALEndcapEtaBorder) { 
	    if (numSignalGammaCandsInSigCone == 0 && hasGsfTrack) {
              category = 14.;
            } else if (numSignalGammaCandsInSigCone >= 1 && hasGsfTrack) {
              category = 16.;
            }   
          }
	  
          mvaValue = std::min(mvaValue, mva_match);
          isGsfElectronMatched = true;
     }      // end of loop over electrons

    if (!isGsfElectronMatched) {
      mvaValue = mva_->MVAValue(*thePFTauRef);
      const reco::PFCandidatePtr& lpfch = thePFTauRef->leadPFChargedHadrCand();
      bool hasGsfTrack = false;
      if (lpfch.isNonnull()) {
        hasGsfTrack = lpfch->gsfTrackRef().isNonnull();
      }

      //// Veto taus that go to Ecal crack
      if (vetoEcalCracks_ &&
          (isInEcalCrack(tauEtaAtEcalEntrance) || isInEcalCrack(leadChargedPFCandEtaAtEcalEntrance))) {
        // add category index
        category_output_->setValue(tauIndex_, category);
        // return MVA output value
        return -99;
      }
      //// Veto taus that go to Ecal crack

      if (std::abs(tauEtaAtEcalEntrance) < ECALBarrelEndcapEtaBorder) {  // Barrel
        if (numSignalGammaCandsInSigCone == 0 && !hasGsfTrack) {
          category = 0.;
        } else if (numSignalGammaCandsInSigCone >= 1 && !hasGsfTrack) {
          category = 2.;
        }
      } else if (std::abs(tauEtaAtEcalEntrance) > ECALBarrelEndcapEtaBorder &&
                 std::abs(tauEtaAtEcalEntrance) < ECALEndcapEtaBorder){  // Endcap
        if (numSignalGammaCandsInSigCone == 0 && !hasGsfTrack) {
          category = 8.;
        } else if (numSignalGammaCandsInSigCone >= 1 && !hasGsfTrack) {
          category = 10.;
        }
      }  else if (std::abs(tauEtaAtEcalEntrance) > ECALEndcapEtaBorder) { //VFEndcap
        if (numSignalGammaCandsInSigCone == 0 && !hasGsfTrack) {
          category = 9.;
        } else if (numSignalGammaCandsInSigCone >= 1 && !hasGsfTrack) {
          category = 11.;
        }
      }
    }
  }

  if (verbosity_) {
    edm::LogPrint("PFTauAgainstEleMVA6") << "<PFRecoTauDiscriminationAgainstElectronMVA6::discriminate>:";
    edm::LogPrint("PFTauAgainstEleMVA6") << " tau: Pt = " << thePFTauRef->pt() << ", eta = " << thePFTauRef->eta()
                                         << ", phi = " << thePFTauRef->phi();
    edm::LogPrint("PFTauAgainstEleMVA6") << " deltaREleTau = " << deltaRDummy
                                         << ", isGsfElectronMatched = " << isGsfElectronMatched;
    edm::LogPrint("PFTauAgainstEleMVA6") << " #Prongs = " << thePFTauRef->signalChargedHadrCands().size();
    edm::LogPrint("PFTauAgainstEleMVA6") << " MVA = " << mvaValue << ", category = " << category;
  }

  // add category index
  category_output_->setValue(tauIndex_, category);
  // return MVA output value
  return mvaValue;
}

void PFRecoTauDiscriminationAgainstElectronMVA6::endEvent(edm::Event& evt) {
  // add all category indices to event
  evt.put(std::move(category_output_), "category");
}

bool PFRecoTauDiscriminationAgainstElectronMVA6::isInEcalCrack(double eta) const {
  double absEta = fabs(eta);
  return (absEta > 1.460 && absEta < 1.558);
}

void PFRecoTauDiscriminationAgainstElectronMVA6::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // pfRecoTauDiscriminationAgainstElectronMVA6
  edm::ParameterSetDescription desc;
  desc.add<double>("minMVANoEleMatchWOgWOgsfBL", 0.0);
  desc.add<edm::InputTag>("PFTauProducer", edm::InputTag("pfTauProducer"));
  desc.add<double>("minMVANoEleMatchWgWOgsfBL", 0.0);
  desc.add<std::string>("mvaName_wGwGSF_EC", "gbr_wGwGSF_EC");
  desc.add<double>("minMVAWgWgsfBL", 0.0);
  desc.add<std::string>("mvaName_woGwGSF_EC", "gbr_woGwGSF_EC");
  desc.add<double>("minMVAWOgWgsfEC", 0.0);
  desc.add<std::string>("mvaName_wGwGSF_BL", "gbr_wGwGSF_BL");
  desc.add<std::string>("mvaName_woGwGSF_BL", "gbr_woGwGSF_BL");
  desc.add<bool>("returnMVA", true);
  desc.add<bool>("loadMVAfromDB", true);
  desc.addOptional<edm::FileInPath>("inputFileName");
  {
    edm::ParameterSetDescription pset_Prediscriminants;
    pset_Prediscriminants.add<std::string>("BooleanOperator", "and");
    {
      edm::ParameterSetDescription psd1;
      psd1.add<double>("cut");
      psd1.add<edm::InputTag>("Producer");
      pset_Prediscriminants.addOptional<edm::ParameterSetDescription>("leadTrack", psd1);
    }
    {
      // encountered this at
      // RecoTauTag/Configuration/python/HPSPFTaus_cff.py
      edm::ParameterSetDescription psd1;
      psd1.add<double>("cut");
      psd1.add<edm::InputTag>("Producer");
      pset_Prediscriminants.addOptional<edm::ParameterSetDescription>("decayMode", psd1);
    }
    desc.add<edm::ParameterSetDescription>("Prediscriminants", pset_Prediscriminants);
  }
  desc.add<std::string>("mvaName_NoEleMatch_woGwoGSF_BL", "gbr_NoEleMatch_woGwoGSF_BL");
  desc.add<bool>("vetoEcalCracks", true);
  desc.add<bool>("usePhiAtEcalEntranceExtrapolation", false);
  desc.add<std::string>("mvaName_NoEleMatch_wGwoGSF_BL", "gbr_NoEleMatch_wGwoGSF_BL");
  desc.add<double>("minMVANoEleMatchWOgWOgsfEC", 0.0);
  desc.add<double>("minMVAWOgWgsfBL", 0.0);
  desc.add<double>("minMVAWgWgsfEC", 0.0);
  desc.add<int>("verbosity", 0);
  desc.add<std::string>("mvaName_NoEleMatch_wGwoGSF_EC", "gbr_NoEleMatch_wGwoGSF_EC");
  desc.add<std::string>("method", "BDTG");
  desc.add<edm::InputTag>("srcGsfElectrons", edm::InputTag("gedGsfElectrons"));
  desc.add<std::string>("mvaName_NoEleMatch_woGwoGSF_EC", "gbr_NoEleMatch_woGwoGSF_EC");
  desc.add<double>("minMVANoEleMatchWgWOgsfEC", 0.0);
  desc.add<std::string>("mvaName_wGwGSF_VFEC", "gbr_wGwGSF_VFEC");
  desc.add<std::string>("mvaName_woGwGSF_VFEC", "gbr_woGwGSF_VFEC");
  desc.add<double>("minMVAWOgWgsfVFEC", 0.0);
  desc.add<double>("minMVAWgWgsfVFEC", 0.0);
  desc.add<std::string>("mvaName_NoEleMatch_wGwoGSF_VFEC", "gbr_NoEleMatch_wGwoGSF_VFEC");
  desc.add<std::string>("mvaName_NoEleMatch_woGwoGSF_VFEC", "gbr_NoEleMatch_woGwoGSF_VFEC");
  desc.add<double>("minMVANoEleMatchWgWOgsfVFEC", 0.0);
  desc.add<double>("minMVANoEleMatchWOgWOgsfVFEC", 0.0);
  //for HGCalEgammaIDHelper
  desc.add<double>("pcaRadius", 3.0);
  desc.add<std::vector<double>>("dEdXWeights")->setComment("This must be copie from dEdX_weights in RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi");  
  desc.add<unsigned int>("isoNRings", 5);
  desc.add<double>("isoDeltaR", 0.15);
  desc.add<double>("isoDeltaRmin", 0.0);
  desc.add<edm::InputTag>("EERecHits", edm::InputTag("HGCalRecHit","HGCEERecHits"));
  desc.add<edm::InputTag>("FHRecHits", edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
  desc.add<edm::InputTag>("BHRecHits", edm::InputTag("HGCalRecHit","HGCHEBRecHits"));
  descriptions.add("pfRecoTauDiscriminationAgainstElectronMVA6", desc);
}

DEFINE_FWK_MODULE(PFRecoTauDiscriminationAgainstElectronMVA6);
