#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iterator>
#include <array>

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

namespace {
  using BranchMap = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;
}

class TauTruthGraphValidator : public DQMEDAnalyzer {
public:
  explicit TauTruthGraphValidator(const edm::ParameterSet&);
  ~TauTruthGraphValidator() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  double checkForAncestralTau(const edm::Handle<truth::Graph> graph, const auto& associations);

  edm::EDGetTokenT<truth::Graph> graphToken_;
  edm::EDGetTokenT<reco::PFTauCollection> hpsTauToken_;
  std::string folder_;

  edm::EDGetTokenT<reco::PFClusterCollection> ecalClustersToken_;
  edm::EDGetTokenT<BranchMap> ecalMapToken_;
  edm::EDGetTokenT<reco::PFClusterCollection> hcalClustersToken_;
  edm::EDGetTokenT<BranchMap> hcalMapToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesBarrelToken_;
  edm::EDGetTokenT<std::vector<TICLCandidate>> ticlCandidatesToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> hgcalTrackstersToken_;
  edm::EDGetTokenT<BranchMap> hgcalTracksterMapToken_;

  edm::EDGetTokenT<BranchMap> ecalTruthToRecoMapToken_;
  edm::EDGetTokenT<BranchMap> hcalTruthToRecoMapToken_;
  edm::EDGetTokenT<BranchMap> hgcalTracksterTruthToRecoMapToken_;

  struct CaloObjectHistograms {
    MonitorElement* pt = nullptr;
    MonitorElement* ptFake = nullptr;

    MonitorElement* eta = nullptr;
    MonitorElement* etaFake = nullptr;

    MonitorElement* energy = nullptr;
    MonitorElement* energyFake = nullptr;
  };

  CaloObjectHistograms ecalClusterHistograms_;
  CaloObjectHistograms hcalClusterHistograms_;
  CaloObjectHistograms hgcalClusterHistograms_;
};

TauTruthGraphValidator::TauTruthGraphValidator(const edm::ParameterSet& iConfig)
  : folder_(iConfig.getParameter<std::string>("folder")) {
  graphToken_ = consumes<truth::Graph>(iConfig.getParameter<edm::InputTag>("graphToken"));
  hpsTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("hpsTauToken"));
  ecalClustersToken_ = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("ecalClustersToken"));
  ecalMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("ecalMapToken"));
  hcalClustersToken_ = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("hcalClustersToken"));
  hcalMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("hcalMapToken"));
  pfCandidatesToken_ = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesToken"));
  pfCandidatesBarrelToken_ = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesBarrelToken"));
  ticlCandidatesToken_ = consumes<std::vector<TICLCandidate>>(iConfig.getParameter<edm::InputTag>("ticlCandidatesToken"));
  hgcalTrackstersToken_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("hgcalTrackstersToken"));
  hgcalTracksterMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("hgcalTracksterMapToken"));
  ecalTruthToRecoMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("ecalTruthToRecoMapToken"));
  hcalTruthToRecoMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("hcalTruthToRecoMapToken"));
  hgcalTracksterTruthToRecoMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("hgcalTracksterTruthToRecoMapToken"));
}

void TauTruthGraphValidator::bookHistograms(DQMStore::IBooker& ibook, edm::Run const&, edm::EventSetup const&) {
  ibook.setCurrentFolder(folder_);
  auto bookClusterHistograms = [&](const std::string& subfolder, CaloObjectHistograms& histograms) {
    ibook.setCurrentFolder(folder_ + "/" + subfolder);
    histograms.pt = ibook.book1D("caloObject_pt", "Tau signal calorimeter objects; pT [GeV]; objects", 100, 0., 200.);
    histograms.ptFake = ibook.book1D("caloObjectFake_pt", "Tau signal calorimeter objects; pT [GeV]; objects", 100, 0., 200.);
    
    histograms.eta = ibook.book1D("caloObject_eta", "Tau signal calorimeter objects; eta; objects", 100, -5., 5.);
    histograms.etaFake = ibook.book1D("caloObjectFake_eta", "Tau signal calorimeter objects; eta; objects", 100, -5., 5.);
    
    histograms.energy = ibook.book1D("caloObject_energy", "Tau signal calorimeter objects; energy [GeV]; objects", 100, 0., 500.);
    histograms.energyFake = ibook.book1D("caloObjectFake_energy", "Tau signal calorimeter objects; energy [GeV]; objects", 100, 0., 500.);
  };
  bookClusterHistograms("ECAL", ecalClusterHistograms_);
  bookClusterHistograms("HCAL", hcalClusterHistograms_);
  bookClusterHistograms("HGCAL", hgcalClusterHistograms_);
}

double TauTruthGraphValidator::checkForAncestralTau(const edm::Handle<truth::Graph> graph, const auto& associations) {
  double score = 2.0;  // defualt score to return if no tau is found
  for (auto const& assoc : associations) {
    const auto particleId = assoc.index();
    if (particleId >= graph->nParticles()) {
      continue;
    }
    auto const& particle = graph->particles()[particleId];
    if (std::abs(particle.pdgId) == 15) {
      score = assoc.score();
      return score;
    }
  }    
  return score;
};

void TauTruthGraphValidator::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  
  // ——————————————————————————————————————————————————————————————————————————
  // Common inputs
  // ——————————————————————————————————————————————————————————————————————————
  
  edm::Handle<truth::Graph> graph;
  iEvent.getByToken(graphToken_, graph);
  edm::Handle<reco::PFTauCollection> hpsTaus;
  iEvent.getByToken(hpsTauToken_, hpsTaus);
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidatesToken_, pfCandidates);
  edm::Handle<reco::PFCandidateCollection> pfCandidatesBarrel;
  iEvent.getByToken(pfCandidatesBarrelToken_, pfCandidatesBarrel);
  edm::Handle<std::vector<TICLCandidate>> ticlCandidates;
  iEvent.getByToken(ticlCandidatesToken_, ticlCandidates);
  edm::Handle<reco::PFClusterCollection> ecalClusters;
  iEvent.getByToken(ecalClustersToken_, ecalClusters);
  edm::Handle<reco::PFClusterCollection> hcalClusters;
  iEvent.getByToken(hcalClustersToken_, hcalClusters);
  edm::Handle<std::vector<ticl::Trackster>> hgcalTracksters;
  iEvent.getByToken(hgcalTrackstersToken_, hgcalTracksters);

  if (!graph.isValid()) {
    std::cout << "############# graph invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "graph invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hpsTaus.isValid()) {
    std::cout << "############# hpsTaus invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hpsTaus invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!pfCandidates.isValid() || !pfCandidatesBarrel.isValid() || !ticlCandidates.isValid()) {
    std::cout << "############# TICL-to-PF input invalid, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "TICL-to-PF input invalid, skipping event " << iEvent.id();
    return;
  }
  if (!ecalClusters.isValid()) {
    std::cout << "############# ecalClusters invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "ecalClusters invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hcalClusters.isValid()) {
    std::cout << "############# hcalClusters invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hcalClusters invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hgcalTracksters.isValid()) {
    std::cout << "############# hgcalTracksters invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hgcalTracksters invalid collection, skipping event " << iEvent.id();
    return;
  }
  
  // ——————————————————————————————————————————————————————————————————————————
  // RecoToTruth 
  // ——————————————————————————————————————————————————————————————————————————
  
  edm::Handle<BranchMap> ecalMap;
  iEvent.getByToken(ecalMapToken_, ecalMap);
  edm::Handle<BranchMap> hcalMap;
  iEvent.getByToken(hcalMapToken_, hcalMap);
  edm::Handle<BranchMap> hgcalTracksterMap;
  iEvent.getByToken(hgcalTracksterMapToken_, hgcalTracksterMap);

  if (!ecalMap.isValid()) {
    std::cout << "############# ecalMap invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "ecalMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hcalMap.isValid()) {
    std::cout << "############# hcalMap invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hcalMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hgcalTracksterMap.isValid()) {
    std::cout << "############# hgcalTracksterMap invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hgcalTracksterMap invalid collection, skipping event " << iEvent.id();
    return;
  }

  auto fillCaloObjectHistograms = [&](double energy,
                                                const auto& position,
                                                CaloObjectHistograms& histograms) {
    const double rho = std::sqrt(position.x() * position.x() + position.y() * position.y());
    const double magnitude = std::sqrt(rho * rho + position.z() * position.z());
    const double pt = magnitude > 0. ? energy * rho / magnitude : 0.;
    histograms.pt->Fill(pt);
    histograms.eta->Fill(position.eta());
    histograms.energy->Fill(energy);
  };

  auto fillCaloObjectFakeHistograms = [&](double energy,
                                                    const auto& position,
                                                    CaloObjectHistograms& histograms,
                                                    double threshold,
                                                    const auto& associations) {
    double tauScore = checkForAncestralTau(graph, associations);
    if (tauScore > threshold) {
      std::cout << "############# Found tau with score = " << tauScore << std::endl;
      const double rho = std::sqrt(position.x() * position.x() + position.y() * position.y());
      const double magnitude = std::sqrt(rho * rho + position.z() * position.z());
      const double pt = magnitude > 0. ? energy * rho / magnitude : 0.;
      histograms.ptFake->Fill(pt);
      histograms.etaFake->Fill(position.eta());
      histograms.energyFake->Fill(energy);
    }
  };

  double tauScoreThreshold = 0.1;  // anything above this is considered a fake

  for (const auto& hpsTau : *hpsTaus) {
    std::cout << "############# hpsTau.pt() = " << hpsTau.pt() << " hpsTau.eta() = " << hpsTau.eta()  << std::endl;
    const auto& pFCands = hpsTau.signalPFCands();
    std::cout << "############# Number of signal PFCandidates in tau = " << pFCands.size() << std::endl;
    for (const auto& pfCand : pFCands) {
      if (pfCand.isNonnull() && pfCand.id() == pfCandidates.id() && pfCand.key() >= pfCandidatesBarrel->size()) {
        const auto ticlCandidateIndex = pfCand.key() - pfCandidatesBarrel->size();
        if (ticlCandidateIndex < ticlCandidates->size()) {
          for (const auto& trackster : ticlCandidates->at(ticlCandidateIndex).tracksters()) {
            if (!trackster.isNonnull() || trackster.id() != hgcalTracksters.id() ||
                trackster.key() >= hgcalTracksters->size()) {
              continue;
            }
            std::cout << "############# HGCAL trackster.raw_energy() = " << trackster->raw_energy()
                      << " trackster.eta() = " << trackster->barycenter().eta() << std::endl;
            fillCaloObjectHistograms(trackster->raw_energy(),
                                     trackster->barycenter(),
                                     hgcalClusterHistograms_);
            fillCaloObjectFakeHistograms(trackster->raw_energy(),
                                         trackster->barycenter(),
                                         hgcalClusterHistograms_,
                                         tauScoreThreshold,
                                         (*hgcalTracksterMap)[trackster.key()]);
          }
        }
      }
      else {
        for (auto const& [block, index] : pfCand->elementsInBlocks()) {
          auto const& element = block->elements()[index];
          auto const cluster = element.clusterRef();
          if (cluster.isNull()) {
            continue;
          }
          if (cluster.id() == ecalClusters.id()) {
            std::cout << "############# ECAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
            fillCaloObjectHistograms(cluster->energy(),
                                     cluster->position(),
                                     ecalClusterHistograms_);
            fillCaloObjectFakeHistograms(cluster->energy(),
                                         cluster->position(),
                                         ecalClusterHistograms_,
                                         tauScoreThreshold,
                                         (*ecalMap)[cluster.key()]);
          }
          else if (cluster.id() == hcalClusters.id()) {
            std::cout << "############# HCAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
            fillCaloObjectHistograms(cluster->energy(),
                                     cluster->position(),
                                     hcalClusterHistograms_);
            fillCaloObjectFakeHistograms(cluster->energy(),
                                         cluster->position(),
                                         hcalClusterHistograms_,
                                         tauScoreThreshold,
                                         (*hcalMap)[cluster.key()]);  
          }
        }
      }
    }
  }

  // ——————————————————————————————————————————————————————————————————————————
  // TruthToReco
  // ——————————————————————————————————————————————————————————————————————————

  edm::Handle<BranchMap> ecalTruthToRecoMap;
  iEvent.getByToken(ecalTruthToRecoMapToken_, ecalTruthToRecoMap);
  edm::Handle<BranchMap> hcalTruthToRecoMap;
  iEvent.getByToken(hcalTruthToRecoMapToken_, hcalTruthToRecoMap);
  edm::Handle<BranchMap> hgcalTracksterTruthToRecoMap;
  iEvent.getByToken(hgcalTracksterTruthToRecoMapToken_, hgcalTracksterTruthToRecoMap);

  if (!ecalTruthToRecoMap.isValid()) {
    std::cout << "############# ecalTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "ecalTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hcalTruthToRecoMap.isValid()) {
    std::cout << "############# hcalTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hcalTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hgcalTracksterTruthToRecoMap.isValid()) {
    std::cout << "############# hgcalTracksterTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hgcalTracksterTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    return;
  }

  for (const auto& truthTau : graph->particleViews()) {
    if (std::abs(truthTau.pdgId()) != 15) continue;
    std::cout << "############# Truth tau pt = " << truthTau.momentum().pt()
              << " eta = " << truthTau.momentum().eta() << std::endl;
    
    auto const ecalClusterMatches = (*ecalTruthToRecoMap)[truthTau.id()];
    if (!ecalClusterMatches.empty()) {
      std::cout << "############# Truth tau has " << ecalClusterMatches.size() << " associated ECAL clusters" << std::endl;
    }

    auto const hcalClusterMatches = (*hcalTruthToRecoMap)[truthTau.id()];
    if (!hcalClusterMatches.empty()) {
      std::cout << "############# Truth tau has " << hcalClusterMatches.size() << " associated HCAL clusters" << std::endl;
    }

    auto const hgcalTracksterMatches = (*hgcalTracksterTruthToRecoMap)[truthTau.id()];
    if (!hgcalTracksterMatches.empty()) {
      std::cout << "############# Truth tau has " << hgcalTracksterMatches.size() << " associated HGCAL tracksters" << std::endl;
    }
  }
}

void TauTruthGraphValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("hpsTauToken", edm::InputTag("hltHpsPFTauProducer"));
  desc.add<std::string>("folder", "HLT/Tau/TruthGraphTauValidation");
  desc.add<edm::InputTag>("graphToken", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("ecalClustersToken", edm::InputTag("hltParticleFlowClusterECALUnseeded"));
  desc.add<edm::InputTag>("ecalMapToken", edm::InputTag("hltTruthBranchPFClusterEcalAssociators", "hltParticleFlowClusterECALRecoToTruthFixed"));
  desc.add<edm::InputTag>("hcalClustersToken", edm::InputTag("hltParticleFlowClusterHCAL"));
  desc.add<edm::InputTag>("hcalMapToken", edm::InputTag("hltTruthBranchPFClusterHcalAssociators","hltParticleFlowClusterHCALRecoToTruthFixed"));
  desc.add<edm::InputTag>("pfCandidatesToken", edm::InputTag("hltParticleFlowTmp"));
  desc.add<edm::InputTag>("pfCandidatesBarrelToken", edm::InputTag("hltParticleFlowTmpBarrel"));
  desc.add<edm::InputTag>("ticlCandidatesToken", edm::InputTag("hltTiclCandidate"));
  desc.add<edm::InputTag>("hgcalTrackstersToken", edm::InputTag("hltTiclCandidate"));
  desc.add<edm::InputTag>("hgcalTracksterMapToken", edm::InputTag("hltTruthBranchTracksterAssociators", "hltTiclTrackstersCLUE3DHighRecoToTruthFixed"));
  desc.add<edm::InputTag>("ecalTruthToRecoMapToken", edm::InputTag("hltTruthBranchPFClusterEcalAssociators", "hltParticleFlowClusterECALTruthToReco"));
  desc.add<edm::InputTag>("hcalTruthToRecoMapToken", edm::InputTag("hltTruthBranchPFClusterHcalAssociators", "hltParticleFlowClusterHCALTruthToReco"));
  desc.add<edm::InputTag>("hgcalTracksterTruthToRecoMapToken", edm::InputTag("hltTruthBranchTracksterAssociators", "hltTiclCandidateTruthToReco"));
  descriptions.add("TauTruthGraphValidator", desc);
}

DEFINE_FWK_MODULE(TauTruthGraphValidator);