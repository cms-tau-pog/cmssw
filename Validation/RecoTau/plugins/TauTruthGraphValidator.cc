#include <string>
#include <vector>
#include <cmath>
#include <iostream>
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

#include "PhysicsTools/TruthInfo/interface/TauTruthUtils.h"
#include "PhysicsTools/TruthInfo/interface/Branch.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "PhysicsTools/TruthInfo/interface/TruthLevels.h"
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
  int getTauDecayMode(const auto& daughters);

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

  struct TruthObjectHistograms {
    MonitorElement* dm = nullptr;
    MonitorElement* dmMatched = nullptr;

    MonitorElement* pt = nullptr;
    MonitorElement* ptMatched = nullptr;

    MonitorElement* eta = nullptr;
    MonitorElement* etaMatched = nullptr;

    MonitorElement* energy = nullptr;
    MonitorElement* energyMatched = nullptr;
  };

  CaloObjectHistograms ecalClusterHistograms_;
  CaloObjectHistograms hcalClusterHistograms_;
  CaloObjectHistograms hgcalClusterHistograms_;
  TruthObjectHistograms truthTauHistograms_;

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
  auto bookClusterHistograms = [&](const std::string& subfolder, CaloObjectHistograms& histograms) {
    ibook.setCurrentFolder(folder_ + "/" + subfolder);
   
    histograms.pt = ibook.book1D("caloObject_pt", "Tau signal calorimeter objects; pT [GeV]; objects", 50, 0., 30.);
    histograms.ptFake = ibook.book1D("caloObjectFake_pt", "Tau signal calorimeter objects; pT [GeV]; objects", 50, 0., 30.);
    
    histograms.eta = ibook.book1D("caloObject_eta", "Tau signal calorimeter objects; eta; objects", 50, -4., 4.);
    histograms.etaFake = ibook.book1D("caloObjectFake_eta", "Tau signal calorimeter objects; eta; objects", 50, -4., 4.);

    histograms.energy = ibook.book1D("caloObject_energy", "Tau signal calorimeter objects; energy [GeV]; objects", 50, 0., 40.);
    histograms.energyFake = ibook.book1D("caloObjectFake_energy", "Tau signal calorimeter objects; energy [GeV]; objects", 50, 0., 40.);
  };

  auto bookTruthHistograms = [&](TruthObjectHistograms& histograms) {
    ibook.setCurrentFolder(folder_ + "/TruthTau");
    
    histograms.dm = ibook.book1D("truthTau_decayMode", "Truth tau decay mode; decay mode; taus", 13, -0.5, 12.5);
    histograms.dmMatched = ibook.book1D("truthTauMatched_decayMode", "Truth tau decay mode; decay mode; taus", 13, -0.5, 12.5);
    
    histograms.pt = ibook.book1D("truthProng_pt", "Truth prong pT; pT [GeV]; prongs", 50, 0., 40.);
    histograms.ptMatched = ibook.book1D("truthProngMatched_pt", "Truth prong pT; pT [GeV]; prongs", 50, 0., 40.);
    
    histograms.eta = ibook.book1D("truthProng_eta", "Truth prong eta; eta; prongs", 50, -5., 5.);
    histograms.etaMatched = ibook.book1D("truthProngMatched_eta", "Truth prong eta; eta; prongs", 50, -5., 5.);
    
    histograms.energy = ibook.book1D("truthProng_energy", "Truth prong energy; energy [GeV]; prongs", 50, 0., 150.);
    histograms.energyMatched = ibook.book1D("truthProngMatched_energy", "Truth prong energy; energy [GeV]; prongs", 50, 0., 150.);
  };

  bookClusterHistograms("ECAL", ecalClusterHistograms_);
  bookClusterHistograms("HCAL", hcalClusterHistograms_);
  bookClusterHistograms("HGCAL", hgcalClusterHistograms_);
  bookTruthHistograms(truthTauHistograms_);
}

double TauTruthGraphValidator::checkForAncestralTau(const edm::Handle<truth::Graph> graph, const auto& associations) {
  double score = 2.0;  // default score to return if no tau is found
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

int TauTruthGraphValidator::getTauDecayMode(const auto& daughters) {
  int nProngs = 0;
  int nPi0s = 0;
  int nOthers = 0;
  for (const auto& daughter : daughters) {
    if (!daughter.hasGen()) {
      std::cout << "Truth tau has non-gen daughter with pdgId = " << daughter.pdgId() << std::endl;
      continue;
    }
    else if (daughter.pdgId() == 22) {
      continue;
    }  
    else if (std::abs(daughter.pdgId()) == 211 || std::abs(daughter.pdgId()) == 321) {
      nProngs++;
    }
    else if (daughter.pdgId() == 111) {
      nPi0s++;
    }
    else {
      nOthers++;
      std::cout << "Truth tau has unexpected daughter with pdgId = " << daughter.pdgId() << std::endl;
    }
  }
  int decayMode = 12;
  if (nOthers > 0) decayMode = 12;  // other
  else if (nProngs == 1 and nPi0s == 0) decayMode = 0;  // 1-prong, no pi0
  else if (nProngs == 1 and nPi0s == 1) decayMode = 1;  // 1-prong, 1 pi0
  else if (nProngs == 1 and nPi0s == 2) decayMode = 2;  // 1-prong, 2 pi0
  else if (nProngs == 3 and nPi0s == 0)decayMode = 10; // 3-prong, no pi0
  else if (nProngs == 3 and nPi0s == 1) decayMode = 11; // 3-prong, 1 pi0
  else {
    std::cout << "Truth tau has unexpected decay, daughter pdgIds = ";
    for (const auto& daughter : daughters) {
      std::cout << daughter.pdgId() << " ";
    }
    std::cout << std::endl;
  }
  return decayMode;
}

void TauTruthGraphValidator::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  
  std::cout << " ——————————————————————————————————————————————————————————————————————————" << std::endl;
  std::cout << "Analyzing " << iEvent.id() << std::endl;

  // ——————————————————————————————————————————————————————————————————————————
  // Common
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

  bool recoTausPresent = true;

  if (!graph.isValid()) {
    std::cout << "############# graph invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "graph invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hpsTaus.isValid()) {
    std::cout << "############# hpsTaus invalid collection " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "hpsTaus invalid collection " << iEvent.id();
    recoTausPresent = false;
  }
  if (!pfCandidates.isValid() || !pfCandidatesBarrel.isValid() || !ticlCandidates.isValid()) {
    std::cout << "############# TICL-to-PF input invalid, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "TICL-to-PF input invalid, skipping event " << iEvent.id();
    return;
  }
  if (!ecalClusters.isValid()) {
    std::cout << "############# ecalClusters invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "ecalClusters invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hcalClusters.isValid()) {
    std::cout << "############# hcalClusters invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "hcalClusters invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hgcalTracksters.isValid()) {
    std::cout << "############# hgcalTracksters invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "hgcalTracksters invalid collection, skipping event " << iEvent.id();
    return;
  }
  
  std::vector<std::vector<size_t>> pfCandToHpsTau(pfCandidates->size(), std::vector<size_t>{});
  std::vector<std::vector<size_t>> tracksterToPfCand(hgcalTracksters->size(), std::vector<size_t>{});
  std::vector<std::vector<size_t>> ecalClusterToPfCand(ecalClusters->size(), std::vector<size_t>{});
  std::vector<std::vector<size_t>> hcalClusterToPfCand(hcalClusters->size(), std::vector<size_t>{});

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
    std::cout << "############# ecalMap invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "ecalMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hcalMap.isValid()) {
    std::cout << "############# hcalMap invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "hcalMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hgcalTracksterMap.isValid()) {
    std::cout << "############# hgcalTracksterMap invalid collection, skipping event " << iEvent.id() << std::endl;
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

  auto fillCaloObjectFakeHistograms = [&](double energy, const auto& position, CaloObjectHistograms& histograms, double threshold, const auto& associations) {
    double tauScore = checkForAncestralTau(graph, associations);
    if (tauScore > threshold) {
      std::cout << "Association score " << tauScore << " is above threshold "
                << threshold << ", marking as fake" << std::endl;
      const double rho = std::sqrt(position.x() * position.x() + position.y() * position.y());
      const double magnitude = std::sqrt(rho * rho + position.z() * position.z());
      const double pt = magnitude > 0. ? energy * rho / magnitude : 0.;
      histograms.ptFake->Fill(pt);
      histograms.etaFake->Fill(position.eta());
      histograms.energyFake->Fill(energy);
    }
    else {
      std::cout << "Association score " << tauScore << " is below threshold "
                << threshold << ", marking as matched" << std::endl;
    }
  };

  const double recoToTruthThreshold = 0.1;  // anything above this is considered a fake

  if (recoTausPresent) {
    for (size_t hpsTauIndex = 0; hpsTauIndex < hpsTaus->size(); ++hpsTauIndex) {
      const auto& hpsTau = (*hpsTaus)[hpsTauIndex];
      std::cout << "############# hpsTau.pt() = " << hpsTau.pt()
                << " hpsTau.eta() = " << hpsTau.eta()  
                << " hpsTau.phi() = " << hpsTau.phi() << std::endl;
      const auto& pfCandsFromHpsTau = hpsTau.signalPFCands();
      std::cout << "Number of signal PFCandidates in tau = " << pfCandsFromHpsTau.size() << std::endl;
      for (const auto& pfCand : pfCandsFromHpsTau) {
        pfCandToHpsTau[pfCand.key()].push_back(hpsTauIndex);
        if (pfCand.isNonnull() && pfCand.id() == pfCandidates.id() && pfCand.key() >= pfCandidatesBarrel->size()) {
          const auto ticlCandidateIndex = pfCand.key() - pfCandidatesBarrel->size();
          if (ticlCandidateIndex < ticlCandidates->size()) {
            for (const auto& trackster : ticlCandidates->at(ticlCandidateIndex).tracksters()) {
              if (!trackster.isNonnull() || trackster.id() != hgcalTracksters.id() ||
                  trackster.key() >= hgcalTracksters->size()) {
                continue;
              }
              tracksterToPfCand[trackster.key()].push_back(pfCand.key());
              std::cout << "HGCAL trackster.raw_energy() = " << trackster->raw_energy()
                        << " trackster.eta() = " << trackster->barycenter().eta() << std::endl;
              fillCaloObjectHistograms(trackster->raw_energy(),
                                      trackster->barycenter(),
                                      hgcalClusterHistograms_);
              fillCaloObjectFakeHistograms(trackster->raw_energy(),
                                          trackster->barycenter(),
                                          hgcalClusterHistograms_,
                                          recoToTruthThreshold,
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
              ecalClusterToPfCand[cluster.key()].push_back(pfCand.key());
              std::cout << "ECAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
              fillCaloObjectHistograms(cluster->energy(),
                                      cluster->position(),
                                      ecalClusterHistograms_);
              fillCaloObjectFakeHistograms(cluster->energy(),
                                          cluster->position(),
                                          ecalClusterHistograms_,
                                          recoToTruthThreshold,
                                          (*ecalMap)[cluster.key()]);
            }
            else if (cluster.id() == hcalClusters.id()) {
              hcalClusterToPfCand[cluster.key()].push_back(pfCand.key());
              std::cout << "HCAL pfCand.pt() = " << pfCand->pt() << " pfCand.eta() = " << pfCand->eta() << std::endl;
              fillCaloObjectHistograms(cluster->energy(),
                                      cluster->position(),
                                      hcalClusterHistograms_);
              fillCaloObjectFakeHistograms(cluster->energy(),
                                          cluster->position(),
                                          hcalClusterHistograms_,
                                          recoToTruthThreshold,
                                          (*hcalMap)[cluster.key()]);  
            }
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
    std::cout << "************* ecalTruthToRecoMap invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "ecalTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hcalTruthToRecoMap.isValid()) {
    std::cout << "************* hcalTruthToRecoMap invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "hcalTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    return;
  }
  if (!hgcalTracksterTruthToRecoMap.isValid()) {
    std::cout << "************* hgcalTracksterTruthToRecoMap invalid collection, skipping event " << iEvent.id() << std::endl;
    edm::LogWarning("TauTruthGraphValidator") << "hgcalTracksterTruthToRecoMap invalid collection, skipping event " << iEvent.id();
    return;
  }

  auto checkIfMatched = [&](const auto& matches, const auto& mapToPf, const double& threshold) {
    if (matches.empty()) return false;
    for (const auto& match : matches) {
      if (match.score() > threshold) continue;  // Skip matches above the threshold
      const auto& pfCandsFromMatch = mapToPf[match.index()];
      if (!pfCandsFromMatch.empty()) {
        for (const auto& pfCand : pfCandsFromMatch) {
          const auto& hpsTausFromPfCand = pfCandToHpsTau[pfCand];
          if (!hpsTausFromPfCand.empty()) {
            std::cout << "Truth tau matched to a cluster with score " << match.score()
                      << " which traces to a HPS tau" << std::endl;
            return true;  // Found a match to a reconstructed tau
          }
        }
      }
    }
    return false;  // No match found
  };

  const double truthToRecoThreshold = 1.0;  // anything below this is considered matched

  for (const uint32_t id : truth::levelAntichain(*graph, truth::levelFromName("tauVisibleHadronic"))) {
    const auto& truthTau = graph->particle(id);

    std::cout << "************* Truth tau pt = " << truthTau.momentum().pt()
              << " eta = " << truthTau.momentum().eta()
              << " phi = " << truthTau.momentum().phi() << std::endl;

    const truth::Branch tauDecayBranch(&truthTau, truth::ClosureSpec::untilPdgId(TauTruthUtils::kTauDecayHadrons));
    const auto daughters = tauDecayBranch.closureLeaves();
    const int decayMode = getTauDecayMode(daughters);
    bool truthTauMatchedToHps = false;

    for (const auto& prong : daughters) {
      if (std::abs(prong.pdgId()) != 211 && std::abs(prong.pdgId()) != 321) continue;
      const auto& prongMomentum = prong.momentum();
      truthTauHistograms_.pt->Fill(prongMomentum.pt());
      truthTauHistograms_.eta->Fill(prongMomentum.eta());
      truthTauHistograms_.energy->Fill(prongMomentum.energy());
      
      const auto ecalClusterMatches = (*ecalTruthToRecoMap)[prong.id()];
      const auto hcalClusterMatches = (*hcalTruthToRecoMap)[prong.id()];
      const auto hgcalTracksterMatches = (*hgcalTracksterTruthToRecoMap)[prong.id()];

      const bool ecalMatchedToHps = checkIfMatched(ecalClusterMatches, ecalClusterToPfCand, truthToRecoThreshold);
      const bool hcalMatchedToHps = checkIfMatched(hcalClusterMatches, hcalClusterToPfCand, truthToRecoThreshold);
      const bool hgcalMatchedToHps = checkIfMatched(hgcalTracksterMatches, tracksterToPfCand, truthToRecoThreshold);

      if (ecalMatchedToHps || hcalMatchedToHps || hgcalMatchedToHps) {
        truthTauHistograms_.ptMatched->Fill(prongMomentum.pt());
        truthTauHistograms_.etaMatched->Fill(prongMomentum.eta());
        truthTauHistograms_.energyMatched->Fill(prongMomentum.energy());
        truthTauMatchedToHps = true;
      }
    }
    truthTauHistograms_.dm->Fill(decayMode);
    if (truthTauMatchedToHps) {
      truthTauHistograms_.dmMatched->Fill(decayMode);
    }
    else {
      std::cout << "Truth tau with DM = " << decayMode << " and daughters: ";
      for (const auto& daughter : daughters) {
        std::cout << daughter.pdgId() << " ";
      }
      std::cout << "not matched to any reconstructed tau" << std::endl;
    }
  }
}

void TauTruthGraphValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("hpsTauToken", edm::InputTag("hltHpsPFTauProducer"));
  desc.add<std::string>("folder", "HLT/Tau/TruthGraphTauValidation");
  desc.add<edm::InputTag>("graphToken", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("ecalClustersToken", edm::InputTag("hltParticleFlowClusterECAL"));
  desc.add<edm::InputTag>("ecalMapToken", edm::InputTag("hltTruthBranchPFClusterEcalAssociators", "hltParticleFlowClusterECALRecoToTruthFixed"));
  desc.add<edm::InputTag>("hcalClustersToken", edm::InputTag("hltParticleFlowClusterHCAL"));
  desc.add<edm::InputTag>("hcalMapToken", edm::InputTag("hltTruthBranchPFClusterHcalAssociators","hltParticleFlowClusterHCALRecoToTruthFixed"));
  desc.add<edm::InputTag>("pfCandidatesToken", edm::InputTag("hltParticleFlowTmp"));
  desc.add<edm::InputTag>("pfCandidatesBarrelToken", edm::InputTag("hltParticleFlowTmpBarrel"));
  desc.add<edm::InputTag>("ticlCandidatesToken", edm::InputTag("hltTiclCandidate"));
  desc.add<edm::InputTag>("hgcalTrackstersToken", edm::InputTag("hltTiclCandidate"));
  desc.add<edm::InputTag>("hgcalTracksterMapToken", edm::InputTag("hltTruthBranchTracksterAssociators", "hltTiclCandidateRecoToTruthFixed"));
  desc.add<edm::InputTag>("ecalTruthToRecoMapToken", edm::InputTag("hltTruthBranchPFClusterEcalAssociators", "hltParticleFlowClusterECALTruthToReco"));
  desc.add<edm::InputTag>("hcalTruthToRecoMapToken", edm::InputTag("hltTruthBranchPFClusterHcalAssociators", "hltParticleFlowClusterHCALTruthToReco"));
  desc.add<edm::InputTag>("hgcalTracksterTruthToRecoMapToken", edm::InputTag("hltTruthBranchTracksterAssociators", "hltTiclCandidateTruthToReco"));
  descriptions.add("TauTruthGraphValidator", desc);
}

DEFINE_FWK_MODULE(TauTruthGraphValidator);
