//%%%%%%%%%%%%%%%%%%%%%%%%
// A workaround for fixing the cross-cleaning bug of boostedTau reco at miniAOD level
// Abddollah Mohammadi, UW-Madison
//%%%%%%%%%%%%%%%%%%%%%%%
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/PatCandidates/interface/TauPFSpecific.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

class PATBoostedTauCleaner : public edm::stream::EDProducer<> {
public:
  explicit PATBoostedTauCleaner(const edm::ParameterSet &);
  ~PATBoostedTauCleaner() override{};

  void produce(edm::Event &, const edm::EventSetup &) override;

private:
  //--- configuration parameters
  edm::EDGetTokenT<pat::TauCollection> srcToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<edm::View<reco::Jet>> jetsCA8Token_;
  bool removeOverLap_;
};

PATBoostedTauCleaner::PATBoostedTauCleaner(const edm::ParameterSet &cfg) {
  srcToken_ = consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("src"));
  vtxToken_ = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vtxLabel"));
  jetsCA8Token_ = consumes<edm::View<reco::Jet>>(cfg.getParameter<edm::InputTag>("ca8JetSrc"));
  removeOverLap_ = cfg.getParameter<bool>("removeOverLap");
  produces<std::vector<pat::Tau>>();
}

void PATBoostedTauCleaner::produce(edm::Event &evt, const edm::EventSetup &es) {
  edm::Handle<pat::TauCollection> inputTaus;
  evt.getByToken(srcToken_, inputTaus);

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(vtxToken_, vertices);

  edm::Handle<edm::View<reco::Jet>> ca8jetHandle;
  evt.getByToken(jetsCA8Token_, ca8jetHandle);

  auto out = std::make_unique<std::vector<pat::Tau>>();
  out->reserve(inputTaus->size());

  for (std::vector<pat::Tau>::const_iterator it = inputTaus->begin(), ed = inputTaus->end(); it != ed; ++it) {
    out->push_back(*it);
    pat::Tau &tau = out->back();

    float chargedPtIsoSum = it->tauID("chargedIsoPtSum");
    float chargedPtIsoSum03 = it->tauID("chargedIsoPtSumdR03");
    float neutralPtIsoSum = it->tauID("neutralIsoPtSum");
    float neutralPtIsoSum03 = it->tauID("neutralIsoPtSumdR03");

    reco::CandidatePtrVector isolationChHPtrs, isolationNHPtrs, isolationGammaPtrs;

    // Tau vertex
    const pat::PackedCandidate *packedLeadChCand =
        dynamic_cast<const pat::PackedCandidate *>(tau.leadChargedHadrCand().get());

    float minDz = 99;
    int tauVertexIdx = 0;
    int idx = 0;
    for (const auto &vertex : *vertices) {
      float dz = std::abs(packedLeadChCand->dz(vertex.position()));
      if (dz < minDz) {
        minDz = dz;
        tauVertexIdx = idx;
      }
      idx++;
    }

    //############################################################################
    // filling a collection of IsoCandidates which either overlappes with other tau's signal candidates or other subJet contitients
    //############################################################################
    if (removeOverLap_) {
      reco::CandidatePtrVector overLappedIsoCand;
      overLappedIsoCand.clear();

      chargedPtIsoSum = 0;
      chargedPtIsoSum03 = 0;
      neutralPtIsoSum = 0;
      neutralPtIsoSum03 = 0;

      for (const reco::CandidatePtr &isoCand1 : tau.isolationCands()) {
        // Check iso candidate does not overlap with other signal candidates
        for (std::vector<pat::Tau>::const_iterator it2 = inputTaus->begin(), ed2 = inputTaus->end(); it2 != ed2;
             ++it2) {
          if (it2 == it)
            continue;

          if (deltaR2(it2->p4(), tau.p4()) > 1.0)
            continue;

          for (const reco::CandidatePtr &sigCand2 : it2->signalCands()) {
            //if (deltaR2(isoCand1->p4(), sigCand2->p4()) < 1e-8)
            if (isoCand1 == sigCand2) //compare edm::Ptrs rather than dR overlap
	      overLappedIsoCand.push_back(isoCand1);
          }
        }

        // Run on CA8 jets
        for (edm::View<reco::Jet>::const_iterator iJet = ca8jetHandle->begin(); iJet != ca8jetHandle->end(); ++iJet) {
          if (deltaR2(iJet->p4(), tau.p4()) > 1.0)
            continue;

          // Find the subjet that seeds taus : closest subjet to tau
          float dR2Closest = 1000;
	  edm::View<reco::Jet>::const_iterator iTauSeedSubJet = ca8jetHandle->end();
          for (edm::View<reco::Jet>::const_iterator iJet = ca8jetHandle->begin(); iJet != ca8jetHandle->end(); ++iJet) {
            if (iJet->pt() < 14 || std::abs(iJet->eta()) > 2.4)
              continue;
            float dR2 = deltaR2(iJet->p4(), tau.p4());
            if (dR2 < dR2Closest) {
              dR2Closest = dR2;
              iTauSeedSubJet = iJet;
            }
          }

          for (edm::View<reco::Jet>::const_iterator jJet = ca8jetHandle->begin(); jJet != ca8jetHandle->end(); ++jJet) {
            if (deltaR2(jJet->p4(), tau.p4()) > 1.0)
              continue;
            if (jJet == iTauSeedSubJet)
              continue;

            for (size_t id = 0; id < jJet->getJetConstituents().size(); id++) {
              const reco::CandidatePtr daughter = jJet->getJetConstituents().at(id);
              //if (deltaR2(isoCand1->p4(), daughter->p4()) < 1e-8)
	      if (isoCand1 == daughter)
                overLappedIsoCand.push_back(isoCand1);
            }
          }
        }
      }  // end of filling the overLappedIsoCand collection
      //############################################################################
      // Setting IsolationChargedHadrCands and recalculating ChargedIsoPtSums
      //############################################################################

      for (const auto &charged : tau.isolationChargedHadrCands()) {
        bool hasOverLap = false;
        for (const reco::CandidatePtr &overLapCand : overLappedIsoCand) {
          //if (deltaR2(charged->p4(), overLapCand->p4()) < 1e-8) {
          if (charged == overLapCand) {
            hasOverLap = true;
            break;
          }
        }
        if (!hasOverLap) {
          isolationChHPtrs.push_back(charged);

          //q-cuts my selection
          if (charged->pt() <= 0.5)
            continue;
          const reco::Track *track = charged->bestTrack();
          if (track == nullptr)
            continue;
          if (std::abs(track->dxy((*vertices)[tauVertexIdx].position())) >= 0.03)
            continue;
          if (track->normalizedChi2() >= 100)
            continue;
          if (track->numberOfValidHits() < 3)
            continue;
          double dz = std::abs(track->dz((*vertices)[tauVertexIdx].position()));
          double dR2 = deltaR2(charged->p4(), tau.p4());
          if (dz < 0.2) {  //from tau vertex
            //iso cone
            if (dR2 < 0.25)  //0.5^2
              chargedPtIsoSum += charged->pt();
            if (dR2 < 0.09)  //0.3^2
              chargedPtIsoSum03 += charged->pt();
          }

        }  //check has overlap
      }    //loop over all isoCandidates
      tau.setIsolationChargedHadrCands(isolationChHPtrs);

      //############################################################################
      // Setting IsolationNeutralHadrCands
      //############################################################################

      for (const reco::CandidatePtr &neutral : tau.isolationNeutrHadrCands()) {
        bool hasOverLap = false;
        for (const reco::CandidatePtr &overLapCand : overLappedIsoCand) {
          //if (deltaR2(neutral->p4(), overLapCand->p4()) < 1e-8) {
          if (neutral == overLapCand) {
            hasOverLap = true;
            break;
          }
        }
        if (!hasOverLap)
          isolationNHPtrs.push_back(neutral);
      }
      tau.setIsolationNeutralHadrCands(isolationNHPtrs);

      //############################################################################
      // Setting IsolationGammaCands and recalculating neutralIsoPtSum
      //############################################################################

      for (const reco::CandidatePtr &gamma : tau.isolationGammaCands()) {
        bool hasOverLap = false;
        for (const reco::CandidatePtr &overLapCand : overLappedIsoCand) {
          //if (deltaR2(gamma->p4(), overLapCand->p4()) < 1e-8) {
          if (gamma == overLapCand) {
            hasOverLap = true;
            break;
          }
        }
        if (!hasOverLap) {
          isolationGammaPtrs.push_back(gamma);

          // Fill neutralPtIsoSum03 and neutralPtIsoSum
          //q-cuts
          if (gamma->pt() <= 1.)
            continue;
          //iso cone
          double dR2 = deltaR2(gamma->p4(), tau.p4());
          if (dR2 < 0.25)  //0.5^2
            neutralPtIsoSum += gamma->pt();
          if (dR2 < 0.09)  //0.3^2
            neutralPtIsoSum03 += gamma->pt();
        }
      }
      tau.setIsolationGammaCands(isolationGammaPtrs);

      //############################################################################
    }  // check if overLap removal is needed
    //############################################################################

    size_t newTauIds = 4;
    size_t nTauIds = tau.tauIDs().size();
    std::vector<pat::Tau::IdPair> tauIds(nTauIds + newTauIds);

    for (size_t q = 0; q < nTauIds; ++q) {
      tauIds[q] = tau.tauIDs().at(q);
    }
    size_t q = nTauIds;
    tauIds[q].first = "chargedIsoPtSumNoOverLap";
    tauIds[q].second = chargedPtIsoSum;
    q = q + 1;

    tauIds[q].first = "chargedIsoPtSumdR03NoOverLap";
    tauIds[q].second = chargedPtIsoSum03;
    q = q + 1;

    tauIds[q].first = "neutralIsoPtSumNoOverLap";
    tauIds[q].second = neutralPtIsoSum;
    q = q + 1;

    tauIds[q].first = "neutralIsoPtSumdR03NoOverLap";
    tauIds[q].second = neutralPtIsoSum03;
    q = q + 1;

    tau.setTauIDs(tauIds);

  }  // iterate over tau is finished

  //This is for Adding new tau Id
  evt.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATBoostedTauCleaner);
