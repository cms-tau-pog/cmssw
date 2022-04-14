#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterAssociation.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
/* 
 * class PFRecoTauDiscriminationByIPCut
 */

using namespace reco;

class PFRecoTauDiscriminationByIPCut : public PFTauDiscriminationProducerBase {
public:
  explicit PFRecoTauDiscriminationByIPCut(const edm::ParameterSet& iConfig)
      : PFTauDiscriminationProducerBase(iConfig),
        m_TausTIPToken(consumes<PFTauTIPAssociationByRef>(iConfig.getParameter<edm::InputTag>("TausIP"))),
        m_tauTIPSelectorString(iConfig.getParameter<std::string>("Cut")),
        m_tauTIPSelector(m_tauTIPSelectorString)
       {
  }
  ~PFRecoTauDiscriminationByIPCut() override {}
  void beginEvent(const edm::Event&, const edm::EventSetup&) override;
  double discriminate(const PFTauRef& pfTau) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  typedef edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef> >       PFTauTIPAssociationByRef;
  edm::EDGetTokenT<PFTauTIPAssociationByRef> m_TausTIPToken;
  edm::Handle<PFTauTIPAssociationByRef> TausTIP;

  const std::string m_tauTIPSelectorString;
  const StringCutObjectSelector<reco::PFTauTransverseImpactParameter> m_tauTIPSelector;
};

double PFRecoTauDiscriminationByIPCut::discriminate(const PFTauRef& thePFTauRef) const {

  return m_tauTIPSelector(*(*TausTIP)[thePFTauRef]) ? 1. : 0.;
}

void PFRecoTauDiscriminationByIPCut::beginEvent(const edm::Event& event, const edm::EventSetup& eventSetup) {
  event.getByToken(m_TausTIPToken, TausTIP);
}

void PFRecoTauDiscriminationByIPCut::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("TausIP", edm::InputTag("hltTauIPCollection"));
  desc.add<std::string>("Cut", "abs(dxy) > -999.");
  {
    edm::ParameterSetDescription psd0;
    psd0.add<std::string>("BooleanOperator", "and");
    desc.add<edm::ParameterSetDescription>("Prediscriminants", psd0);
  }
  desc.add<edm::InputTag>("PFTauProducer", edm::InputTag("pfRecoTauProducer"));
  descriptions.add("pfRecoTauDiscriminationByIPCut", desc);
}

DEFINE_FWK_MODULE(PFRecoTauDiscriminationByIPCut);
