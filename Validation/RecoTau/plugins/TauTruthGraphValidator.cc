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
};

TauTruthGraphValidator::TauTruthGraphValidator(const edm::ParameterSet& iConfig)
  : folder_(iConfig.getParameter<std::string>("folder")) {
  graphToken_ = consumes<truth::Graph>(iConfig.getParameter<edm::InputTag>("graphToken"));
  hpsTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("hpsTauToken"));
  ecalClustersToken_ = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("ecalClustersToken"));
  ecalMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("ecalMapToken"));
  hcalClustersToken_ = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("hcalClustersToken"));
  hcalMapToken_ = consumes<BranchMap>(iConfig.getParameter<edm::InputTag>("hcalMapToken"));
}

void TauTruthGraphValidator::bookHistograms(DQMStore::IBooker& ibook, edm::Run const&, edm::EventSetup const&) {
  ibook.setCurrentFolder(folder_);
  ibook.book1D("my_histo", "; pT [GeV]; entries", 60, 0., 120.);
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

  for (const auto& hpsTau : *hpsTaus) {
    std::cout << "############# hpsTau.pt() = " << hpsTau.pt() << "hpsTau.eta() = " << hpsTau.eta()  << std::endl;
    const auto PFCands = hpsTau.signalPFCands();
    for (const auto& pfCand : PFCands) {
      for (auto const& [block, index] : pfCand->elementsInBlocks()) {
        auto const& element = block->elements()[index];
        auto const cluster = element.clusterRef();
        if (cluster.isNull()) {
          continue;
        }
        if (cluster.id() == ecalClusters.id()) {
          std::cout << "############# ECAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
          for (auto const& assoc : (*ecalMap)[cluster.key()]) {
            std::cout << "############# index " << assoc.index()
              << " sharedEnergy = " << assoc.sharedEnergy()
              << " score = " << assoc.score() << std::endl;
          }
        }
        else if (cluster.id() == hcalClusters.id()) {
          std::cout << "############# HCAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
          for (auto const& assoc : (*hcalMap)[cluster.key()]) {
            std::cout << "############# index " << assoc.index()
              << " sharedEnergy = " << assoc.sharedEnergy()
              << " score = " << assoc.score() << std::endl;
          }
        }
        else {
          std::cout << "############# HGCAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
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
  descriptions.add("TauTruthGraphValidator", desc);
}

DEFINE_FWK_MODULE(TauTruthGraphValidator);