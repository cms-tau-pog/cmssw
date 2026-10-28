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

  struct CaloObjectHistograms {
    MonitorElement* pt = nullptr;
    MonitorElement* eta = nullptr;
    MonitorElement* energy = nullptr;
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
}

void TauTruthGraphValidator::bookHistograms(DQMStore::IBooker& ibook, edm::Run const&, edm::EventSetup const&) {
  ibook.setCurrentFolder(folder_);
  auto bookClusterHistograms = [&](const std::string& subfolder, CaloObjectHistograms& histograms) {
    ibook.setCurrentFolder(folder_ + "/" + subfolder);
    histograms.pt = ibook.book1D("caloObject_pt", "Tau signal calorimeter objects; pT [GeV]; objects", 100, 0., 200.);
    histograms.eta = ibook.book1D("caloObject_eta", "Tau signal calorimeter objects; eta; objects", 100, -5., 5.);
    histograms.energy = ibook.book1D("caloObject_energy", "Tau signal calorimeter objects; energy [GeV]; objects", 100, 0., 500.);
  };
  bookClusterHistograms("ECAL", ecalClusterHistograms_);
  bookClusterHistograms("HCAL", hcalClusterHistograms_);
  bookClusterHistograms("HGCAL", hgcalClusterHistograms_);
}

void TauTruthGraphValidator::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<reco::PFTauCollection> hpsTaus;
  iEvent.getByToken(hpsTauToken_, hpsTaus);
  edm::Handle<truth::Graph> graph;
  iEvent.getByToken(graphToken_, graph);
  edm::Handle<reco::PFClusterCollection> ecalClusters;
  iEvent.getByToken(ecalClustersToken_, ecalClusters);
  edm::Handle<BranchMap> ecalMap;
  iEvent.getByToken(ecalMapToken_, ecalMap);
  edm::Handle<reco::PFClusterCollection> hcalClusters;
  iEvent.getByToken(hcalClustersToken_, hcalClusters);
  edm::Handle<BranchMap> hcalMap;
  iEvent.getByToken(hcalMapToken_, hcalMap);
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidatesToken_, pfCandidates);
  edm::Handle<reco::PFCandidateCollection> pfCandidatesBarrel;
  iEvent.getByToken(pfCandidatesBarrelToken_, pfCandidatesBarrel);
  edm::Handle<std::vector<TICLCandidate>> ticlCandidates;
  iEvent.getByToken(ticlCandidatesToken_, ticlCandidates);
  edm::Handle<std::vector<ticl::Trackster>> hgcalTracksters;
  iEvent.getByToken(hgcalTrackstersToken_, hgcalTracksters);
  edm::Handle<BranchMap> hgcalTracksterMap;
  iEvent.getByToken(hgcalTracksterMapToken_, hgcalTracksterMap);

  if (!hpsTaus.isValid()) {
    std::cout << "############# hpsTaus invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hpsTaus invalid collection, skipping event " << iEvent.id();
    return;
  }

  if (!graph.isValid()) {
    std::cout << "############# graph invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "graph invalid collection, skipping event " << iEvent.id();
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
  if (!pfCandidates.isValid() || !pfCandidatesBarrel.isValid() || !ticlCandidates.isValid()) {
    std::cout << "############# TICL-to-PF input invalid, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "TICL-to-PF input invalid, skipping event " << iEvent.id();
    return;
  }

  if (!hgcalTracksterMap.isValid()) {
    std::cout << "############# hgcalTracksterMap invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hgcalTracksterMap invalid collection, skipping event " << iEvent.id();
    return;
  }

  auto fillCaloObjectHistograms = [&](double energy, const auto& position, CaloObjectHistograms& histograms) {
    const double rho = std::sqrt(position.x() * position.x() + position.y() * position.y());
    const double magnitude = std::sqrt(rho * rho + position.z() * position.z());
    const double pt = magnitude > 0. ? energy * rho / magnitude : 0.;
    histograms.pt->Fill(pt);
    histograms.eta->Fill(position.eta());
    histograms.energy->Fill(energy);
  };

  auto printMatchedTruthParticles = [&](const char* detector, std::size_t recoIndex, const auto& associations) {
    for (auto const& assoc : associations) {
      const auto particleId = assoc.index();
      if (particleId >= graph->nParticles()) {
        continue;
      }
      auto const& particle = graph->particles()[particleId];
      auto const& momentum = particle.momentum;
      std::cout << "############# " << detector << " recoIndex=" << recoIndex
                << " matchedTruthParticleId=" << particleId << " pdgId=" << particle.pdgId
                << " pt=" << momentum.pt() << " eta=" << momentum.eta() << " phi=" << momentum.phi()
                << " energy=" << momentum.E() << " sharedEnergy=" << assoc.sharedEnergy()
                << " score=" << assoc.score() << std::endl;
    }
  };

  for (const auto& hpsTau : *hpsTaus) {
    std::cout << "############# hpsTau.pt() = " << hpsTau.pt() << " hpsTau.eta() = " << hpsTau.eta()  << std::endl;
    const auto PFCands = hpsTau.signalPFCands();
    std::cout << "############# Number of signal PFCandidates in tau = " << PFCands.size() << std::endl;
    for (const auto& pfCand : PFCands) {
      if (pfCand.isNonnull() && pfCand.id() == pfCandidates.id() && pfCand.key() >= pfCandidatesBarrel->size()) {
        const auto ticlCandidateIndex = pfCand.key() - pfCandidatesBarrel->size();
        if (ticlCandidateIndex < ticlCandidates->size()) {
          for (const auto& trackster : ticlCandidates->at(ticlCandidateIndex).tracksters()) {
            if (!trackster.isNonnull() || trackster.id() != hgcalTracksters.id() ||
                trackster.key() >= hgcalTracksters->size()) {
              continue;
            }
            const auto tracksterIndex = trackster.key();
            fillCaloObjectHistograms(trackster->raw_energy(), trackster->barycenter(), hgcalClusterHistograms_);
            std::cout << "############# HGCAL trackster.raw_energy() = " << trackster->raw_energy()
                      << " trackster.eta() = " << trackster->barycenter().eta() << std::endl;
            printMatchedTruthParticles("HGCAL trackster", tracksterIndex, (*hgcalTracksterMap)[tracksterIndex]);
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
            fillCaloObjectHistograms(cluster->energy(), cluster->position(), ecalClusterHistograms_);
            std::cout << "############# ECAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
            printMatchedTruthParticles("ECAL", cluster.key(), (*ecalMap)[cluster.key()]);
          }
          else if (cluster.id() == hcalClusters.id()) {
            fillCaloObjectHistograms(cluster->energy(), cluster->position(), hcalClusterHistograms_);
            std::cout << "############# HCAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
            printMatchedTruthParticles("HCAL", cluster.key(), (*hcalMap)[cluster.key()]);
          }
        }
      }
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
  descriptions.add("TauTruthGraphValidator", desc);
}

DEFINE_FWK_MODULE(TauTruthGraphValidator);