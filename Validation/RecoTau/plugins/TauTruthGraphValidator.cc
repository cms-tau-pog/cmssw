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
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "SimDataFormats/CaloAnalysis/interface/SimTauCPLink.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

class TauTruthGraphValidator : public DQMEDAnalyzer {
public:
  explicit TauTruthGraphValidator(const edm::ParameterSet&);
  ~TauTruthGraphValidator() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<reco::PFTauCollection> hpsTauToken_;
  std::string folder_;
};

TauTruthGraphValidator::TauTruthGraphValidator(const edm::ParameterSet& iConfig)
  : folder_(iConfig.getParameter<std::string>("folder")) {
  hpsTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("hpsTauToken"));
}

void TauTruthGraphValidator::bookHistograms(DQMStore::IBooker& ibook, edm::Run const&, edm::EventSetup const&) {
  ibook.setCurrentFolder(folder_);
  ibook.book1D("my_histo", "; pT [GeV]; entries", 60, 0., 120.);
}

void TauTruthGraphValidator::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<reco::PFTauCollection> hpsTaus;
  iEvent.getByToken(hpsTauToken_, hpsTaus);

  if (!hpsTaus.isValid()) {
    std::cout << "hpsTaus invalid collection, skipping event " << iEvent.id();
    edm::LogWarning("TauTruthGraphValidator") << "hpsTaus invalid collection, skipping event " << iEvent.id();
    return;
  }

  for (const auto& hpsTau : *hpsTaus) {
    std::cout << "hpsTau.pt() = " << hpsTau.pt() << std::endl;
  }

}

void TauTruthGraphValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("hpsTauToken", edm::InputTag("hltHpsPFTauProducer"));
  desc.add<std::string>("folder", "RecoTauV/TruthGraphTauValidator");
  descriptions.add("TauTruthGraphValidator", desc);
}

DEFINE_FWK_MODULE(TauTruthGraphValidator);