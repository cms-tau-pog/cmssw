// -*- C++ -*-
//
// Package:    BoostTau/BoostedTauAnalyzer
// Class:      BoostedTauAnalyzer
//
/**\class BoostedTauAnalyzer BoostedTauAnalyzer.cc BoostTau/BoostedTauAnalyzer/plugins/BoostedTauAnalyzer.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Abdollah Mohammadi
//         Created:  Thu, 21 May 2020 17:33:45 GMT
//
//

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class BoostedTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit BoostedTauAnalyzer(const edm::ParameterSet&);
  //      ~BoostedTauAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void fillBoostedTaus(const edm::Event&);

  void branchesBoostedTaus(TTree*);

  // ----------member data ---------------------------

  edm::EDGetTokenT<std::vector<pat::Tau> > boostedTauCollection_;
  edm::EDGetTokenT<std::vector<pat::Tau> > boostedTauCollectionNoOverLap_;

  TTree* tree_;
};

//
// constructors and destructor
//
BoostedTauAnalyzer::BoostedTauAnalyzer(const edm::ParameterSet& iConfig)
    : boostedTauCollection_(consumes<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("boostedTauSrc"))),
      boostedTauCollectionNoOverLap_(
          consumes<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("cleanedBoostedTauSrc"))) {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("EventTree", "Event data");
  branchesBoostedTaus(tree_);
}

// ------------ method called for each event  ------------
void BoostedTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  fillBoostedTaus(iEvent);
  tree_->Fill();
}

//Tau Id & Isolation

vector<float> boostedTauByIsolationMVArun2v1DBoldDMwLTraw__;
vector<float> boostedTauByIsolationMVArun2v1DBoldDMwLTrawNew__;

vector<bool> boostedTauByLooseIsolationMVArun2v1DBoldDMwLTNew__;
vector<bool> boostedTauByLooseIsolationMVArun2v1DBoldDMwLT__;

vector<float> boostedTauByDeepTau2017v1VSjetraw__;
vector<float> boostedTauByDeepTau2017v1VSeraw__;
vector<float> boostedTauByDeepTau2017v1VSmuraw__;

void BoostedTauAnalyzer::branchesBoostedTaus(TTree* tree) {
  tree->Branch("boostedTauByIsolationMVArun2v1DBoldDMwLTraw", &boostedTauByIsolationMVArun2v1DBoldDMwLTraw__);
  tree->Branch("boostedTauByIsolationMVArun2v1DBoldDMwLTrawNew", &boostedTauByIsolationMVArun2v1DBoldDMwLTrawNew__);

  tree->Branch("boostedTauByLooseIsolationMVArun2v1DBoldDMwLT", &boostedTauByLooseIsolationMVArun2v1DBoldDMwLT__);
  tree->Branch("boostedTauByLooseIsolationMVArun2v1DBoldDMwLTNew", &boostedTauByLooseIsolationMVArun2v1DBoldDMwLTNew__);

  tree->Branch("boostedTauByDeepTau2017v1VSjetraw", &boostedTauByDeepTau2017v1VSjetraw__);
  tree->Branch("boostedTauByDeepTau2017v1VSeraw", &boostedTauByDeepTau2017v1VSeraw__);
  tree->Branch("boostedTauByDeepTau2017v1VSmuraw", &boostedTauByDeepTau2017v1VSmuraw__);
}

void BoostedTauAnalyzer::fillBoostedTaus(const edm::Event& e) {
  // Tau Id & Isolation

  boostedTauByIsolationMVArun2v1DBoldDMwLTraw__.clear();
  boostedTauByIsolationMVArun2v1DBoldDMwLTrawNew__.clear();

  boostedTauByLooseIsolationMVArun2v1DBoldDMwLT__.clear();
  boostedTauByLooseIsolationMVArun2v1DBoldDMwLTNew__.clear();

  boostedTauByDeepTau2017v1VSjetraw__.clear();
  boostedTauByDeepTau2017v1VSeraw__.clear();
  boostedTauByDeepTau2017v1VSmuraw__.clear();

  edm::Handle<vector<pat::Tau> > boostedTauHandleNoOverLap;
  e.getByToken(boostedTauCollectionNoOverLap_, boostedTauHandleNoOverLap);

  if (!boostedTauHandleNoOverLap.isValid()) {
    edm::LogWarning("BoostedTauAnalyzer") << "no pat::Tau in event boostedTauHandleNoOverLap";
    return;
  }

  for (vector<pat::Tau>::const_iterator itau = boostedTauHandleNoOverLap->begin();
       itau != boostedTauHandleNoOverLap->end();
       ++itau) {
    boostedTauByIsolationMVArun2v1DBoldDMwLTraw__.push_back(itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw"));
    boostedTauByLooseIsolationMVArun2v1DBoldDMwLT__.push_back(itau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT"));

    boostedTauByIsolationMVArun2v1DBoldDMwLTrawNew__.push_back(itau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
    boostedTauByLooseIsolationMVArun2v1DBoldDMwLTNew__.push_back(
        itau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));

    boostedTauByDeepTau2017v1VSjetraw__.push_back(itau->tauID("byDeepTau2017v1VSjetraw"));
    boostedTauByDeepTau2017v1VSeraw__.push_back(itau->tauID("byDeepTau2017v1VSeraw"));
    boostedTauByDeepTau2017v1VSmuraw__.push_back(itau->tauID("byDeepTau2017v1VSmuraw"));

  }  // loop over tau candidates
}

// ------------ method called once each job just before starting event loop  ------------
void BoostedTauAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void BoostedTauAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BoostedTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BoostedTauAnalyzer);
