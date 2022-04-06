#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "RecoTauTag/RecoTau/interface/DeepTauBase.h"
#include "RecoTauTag/RecoTau/interface/DeepTauHelper.h"

namespace deep_tau {
  constexpr int NumberOfOutputs = 4;
}

using namespace deep_tau_2017;
using namespace deeptau_helper;

class DeepTauIdSonicProducer : public TritonEDProducer<> {
public:
  explicit DeepTauIdSonicProducer(edm::ParameterSet const& cfg)
      : TritonEDProducer<>(cfg),
        tausToken_(consumes<TauCollection>(cfg.getParameter<edm::InputTag>("taus"))),
        pfcandToken_(consumes<CandidateCollection>(cfg.getParameter<edm::InputTag>("pfcands"))),
        vtxToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
        electrons_token_(consumes<std::vector<pat::Electron>>(cfg.getParameter<edm::InputTag>("electrons"))),
        muons_token_(consumes<std::vector<pat::Muon>>(cfg.getParameter<edm::InputTag>("muons"))),
        rho_token_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
        basicTauDiscriminators_inputToken_(consumes<reco::TauDiscriminatorContainer>(
            cfg.getUntrackedParameter<edm::InputTag>("basicTauDiscriminators"))),
        basicTauDiscriminatorsdR03_inputToken_(consumes<reco::TauDiscriminatorContainer>(
            cfg.getUntrackedParameter<edm::InputTag>("basicTauDiscriminatorsdR03"))),
        pfTauTransverseImpactParameters_token_(
            consumes<edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>>(
                cfg.getParameter<edm::InputTag>("pfTauTransverseImpactParameters"))),
        debug_level(cfg.getParameter<int>("debug_level")),
        disable_dxy_pca_(cfg.getParameter<bool>("disable_dxy_pca")),
        disable_hcalFraction_workaround_(cfg.getParameter<bool>("disable_hcalFraction_workaround")),
        disable_CellIndex_workaround_(cfg.getParameter<bool>("disable_CellIndex_workaround")),
        doSplitVersion_(cfg.getParameter<bool>("doSplitVersion")),
        outputdiscs_(GetOutputs()) {
    for (const auto& output_desc : outputdiscs_) {
      produces<TauDiscriminator>(output_desc.first);
      const auto& cut_list = cfg.getParameter<std::vector<std::string>>(output_desc.first + "WP");
      for (const std::string& cut_str : cut_list) {
        workingPoints_[output_desc.first].push_back(std::make_unique<Cutter>(cut_str));
      }
    }

    // prediscriminant operator
    // require the tau to pass the following prediscriminants
    const edm::ParameterSet& prediscriminantConfig = cfg.getParameter<edm::ParameterSet>("Prediscriminants");

    // determine boolean operator used on the prediscriminants
    std::string pdBoolOperator = prediscriminantConfig.getParameter<std::string>("BooleanOperator");
    // convert string to lowercase
    transform(pdBoolOperator.begin(), pdBoolOperator.end(), pdBoolOperator.begin(), ::tolower);

    if (pdBoolOperator == "and") {
      andPrediscriminants_ = 0x1;  //use chars instead of bools so we can do a bitwise trick later
    } else if (pdBoolOperator == "or") {
      andPrediscriminants_ = 0x0;
    } else {
      throw cms::Exception("TauDiscriminationProducerBase")
          << "PrediscriminantBooleanOperator defined incorrectly, options are: AND,OR";
    }

    // get the list of prediscriminants
    std::vector<std::string> prediscriminantsNames =
        prediscriminantConfig.getParameterNamesForType<edm::ParameterSet>();

    for (auto const& iDisc : prediscriminantsNames) {
      const edm::ParameterSet& iPredisc = prediscriminantConfig.getParameter<edm::ParameterSet>(iDisc);
      const edm::InputTag& label = iPredisc.getParameter<edm::InputTag>("Producer");
      double cut = iPredisc.getParameter<double>("cut");

      PATTauDiscInfo thisDiscriminator;
      thisDiscriminator.label = label;
      thisDiscriminator.cut = cut;
      thisDiscriminator.disc_token = consumes<pat::PATTauDiscriminator>(label);
      patPrediscriminants_.push_back(thisDiscriminator);
    }
  }

  using TauDiscriminator = deep_tau::DeepTauBase::TauDiscriminator;
  using TauCollection = deep_tau::DeepTauBase::TauCollection;
  using CandidateCollection = deep_tau::DeepTauBase::CandidateCollection;
  using TauRef = deep_tau::DeepTauBase::TauRef;
  using TauRefProd = deep_tau::DeepTauBase::TauRefProd;
  using ElectronCollection = deep_tau::DeepTauBase::ElectronCollection;
  using MuonCollection = deep_tau::DeepTauBase::MuonCollection;
  using Cutter = deep_tau::TauWPThreshold;
  using CutterPtr = deep_tau::DeepTauBase::CutterPtr;
  using WPList = deep_tau::DeepTauBase::WPList;
  using BasicDiscriminator = deep_tau::DeepTauBase::BasicDiscriminator;
  using PATTauDiscInfo = deep_tau::DeepTauBase::TauDiscInfo<pat::PATTauDiscriminator>;
  using OutputCollection = deep_tau::DeepTauBase::OutputCollection;

  void acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) override;
  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) override;
  void createOutputs(edm::Event& event, const std::vector<std::vector<float>>& pred, edm::Handle<TauCollection> taus);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  float scale_and_rm_outlier(float val, float scale);

  // select boolean operation on prediscriminants (and = 0x01, or = 0x00)
  uint8_t andPrediscriminants_;
  std::vector<PATTauDiscInfo> patPrediscriminants_;

private:
  edm::EDGetTokenT<TauCollection> tausToken_;
  edm::EDGetTokenT<CandidateCollection> pfcandToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electrons_token_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muons_token_;
  edm::EDGetTokenT<double> rho_token_;
  edm::EDGetTokenT<reco::TauDiscriminatorContainer> basicTauDiscriminators_inputToken_;
  edm::EDGetTokenT<reco::TauDiscriminatorContainer> basicTauDiscriminatorsdR03_inputToken_;
  edm::EDGetTokenT<edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>>
      pfTauTransverseImpactParameters_token_;
  std::string input_layer_, output_layer_;
  const int debug_level;
  const bool disable_dxy_pca_;
  const bool disable_hcalFraction_workaround_;
  const bool disable_CellIndex_workaround_;
  const bool doSplitVersion_;

  OutputCollection outputdiscs_;
  std::map<std::string, WPList> workingPoints_;

  std::vector<size_t> tau_indices_;

  static constexpr float pi = M_PI;
  std::map<BasicDiscriminator, size_t> basicDiscrIndexMap_;
  std::map<BasicDiscriminator, size_t> basicDiscrdR03IndexMap_;

  template <typename CandidateCastType, typename TauCastType>
  void getPredictionsV2(TauCollection::const_reference& tau,
                        const size_t tau_index,
                        const edm::RefToBase<reco::BaseTau> tau_ref,
                        const std::vector<pat::Electron>* electrons,
                        const std::vector<pat::Muon>* muons,
                        const edm::View<reco::Candidate>& pfCands,
                        const reco::Vertex& pv,
                        double rho,
                        const TauFunc& tau_funcs,
                        std::vector<float>* p_tauBlockInputs,
                        std::vector<float>* p_egammaInnerBlockInputs,
                        std::vector<float>* p_muonInnerBlockInputs,
                        std::vector<float>* p_hadronInnerBlockInputs,
                        std::vector<float>* p_egammaOuterBlockInputs,
                        std::vector<float>* p_muonOuterBlockInputs,
                        std::vector<float>* p_hadronOuterBlockInputs,
                        std::vector<int64_t>* p_innerGridposInputs,
                        std::vector<int64_t>* p_outerGridposInputs);

  template <typename CandidateCastType, typename TauCastType>
  void createConvFeatures(const TauCastType& tau,
                          const size_t tau_index,
                          const edm::RefToBase<reco::BaseTau> tau_ref,
                          const reco::Vertex& pv,
                          double rho,
                          const std::vector<pat::Electron>* electrons,
                          const std::vector<pat::Muon>* muons,
                          const edm::View<reco::Candidate>& pfCands,
                          const CellGrid& grid,
                          const TauFunc& tau_funcs,
                          bool is_inner,
                          std::vector<float>* p_egammaBlockInputs,
                          std::vector<float>* p_muonBlockInputs,
                          std::vector<float>* p_hadronBlockInputs,
                          std::vector<int64_t>* p_GridposInputs);
};

void DeepTauIdSonicProducer::acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) {
  edm::Handle<TauCollection> taus;
  iEvent.getByToken(tausToken_, taus);

  edm::ProductID tauProductID = taus.id();

  // load prediscriminators
  size_t nPrediscriminants = patPrediscriminants_.size();
  for (size_t iDisc = 0; iDisc < nPrediscriminants; ++iDisc) {
    edm::ProductID discKeyId;
    patPrediscriminants_[iDisc].fill(iEvent);
    discKeyId = patPrediscriminants_[iDisc].handle->keyProduct().id();

    // Check to make sure the product is correct for the discriminator.
    // If not, throw a more informative exception.
    if (tauProductID != discKeyId) {
      throw cms::Exception("MisconfiguredPrediscriminant")
          << "The tau collection has product ID: " << tauProductID
          << " but the pre-discriminator is keyed with product ID: " << discKeyId << std::endl;
    }
  }

  const reco::TauDiscriminatorContainer basicTauDiscriminators_default;
  const reco::TauDiscriminatorContainer basicTauDiscriminatorsdR03_default;
  const edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>
      pfTauTransverseImpactParameters_default;

  const std::vector<pat::Electron>* electron_collection;
  const std::vector<pat::Muon>* muon_collection;
  const reco::TauDiscriminatorContainer* basicTauDiscriminators;
  const reco::TauDiscriminatorContainer* basicTauDiscriminatorsdR03;
  const edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>*
      pfTauTransverseImpactParameters;

  electron_collection = &iEvent.get(electrons_token_);
  muon_collection = &iEvent.get(muons_token_);
  pfTauTransverseImpactParameters = &pfTauTransverseImpactParameters_default;
  basicTauDiscriminators = &basicTauDiscriminators_default;
  basicTauDiscriminatorsdR03 = &basicTauDiscriminatorsdR03_default;

  TauFunc tauIDs = {basicTauDiscriminators,
                    basicTauDiscriminatorsdR03,
                    pfTauTransverseImpactParameters,
                    basicDiscrIndexMap_,
                    basicDiscrdR03IndexMap_};

  edm::Handle<edm::View<reco::Candidate>> pfCands;
  iEvent.getByToken(pfcandToken_, pfCands);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  edm::Handle<double> rho;
  iEvent.getByToken(rho_token_, rho);

  // vector to store the indices for the taus passing the selections
  tau_indices_.clear();

  for (size_t tau_index = 0; tau_index < taus->size(); ++tau_index) {
    const edm::RefToBase<reco::BaseTau> tauRef = taus->refAt(tau_index);
    bool passesPrediscriminants =
        tauIDs.passPrediscriminants<std::vector<PATTauDiscInfo>>(patPrediscriminants_, andPrediscriminants_, tauRef);
    if (!passesPrediscriminants)
      continue;

    // tau index that passes the selection
    tau_indices_.push_back(tau_index);
  }

  if (tau_indices_.empty()) {
    // no tau passing the requirement
    // no need to run acquire and inference
    client_->setBatchSize(0);
    return;
  }

  int n_taus = tau_indices_.size();
  if (doSplitVersion_) {
    // always set the batch size to 1, since the 'batch' for
    // deeptau is different with the traditional ones
    client_->setBatchSize(1);
  } else {
    // for the regular non-split deep-tau model, set the
    // batch size to the number of taus per event
    client_->setBatchSize(n_taus);
  }

  auto& input_tauBlock = iInput.at("input_tau");
  auto& input_innerEgammaBlock = iInput.at("input_inner_egamma");
  auto& input_outerEgammaBlock = iInput.at("input_outer_egamma");
  auto& input_innerMuonBlock = iInput.at("input_inner_muon");
  auto& input_outerMuonBlock = iInput.at("input_outer_muon");
  auto& input_innerHadronBlock = iInput.at("input_inner_hadrons");
  auto& input_outerHadronBlock = iInput.at("input_outer_hadrons");

  TritonInputData* input_innerGridposBlock = nullptr;
  TritonInputData* input_outerGridposBlock = nullptr;
  TritonInputContainer<int64_t> data_innerGridposBlock = nullptr;
  TritonInputContainer<int64_t> data_outerGridposBlock = nullptr;
  std::vector<int64_t>* p_vdata_innerGridposBlock = nullptr;
  std::vector<int64_t>* p_vdata_outerGridposBlock = nullptr;

  if (doSplitVersion_) {
    // for inner and outer grids
    // usually less than 10 inner grids and 50 outer grids per tau
    // set these numbers temporarily for vector reservation
    int n_inner_cells = 10 * n_taus;
    int n_outer_cells = 50 * n_taus;
    input_tauBlock.setShape(0, n_taus);
    input_innerEgammaBlock.setShape(0, n_inner_cells);
    input_outerEgammaBlock.setShape(0, n_outer_cells);
    input_innerMuonBlock.setShape(0, n_inner_cells);
    input_outerMuonBlock.setShape(0, n_outer_cells);
    input_innerHadronBlock.setShape(0, n_inner_cells);
    input_outerHadronBlock.setShape(0, n_outer_cells);

    // coordinates of the inner grids: n_inner_cells x 3 (i_tau, j_eta, k_phi)
    input_innerGridposBlock = &(iInput.at("input_inner_pos"));
    input_innerGridposBlock->setShape(0, n_inner_cells);
    data_innerGridposBlock = input_innerGridposBlock->allocate<int64_t>();
    p_vdata_innerGridposBlock = &((*data_innerGridposBlock)[0]);

    // coordinates of the outer grids: n_outer_cells x 3 (i_tau, j_eta, k_phi)
    input_outerGridposBlock = &(iInput.at("input_outer_pos"));
    input_outerGridposBlock->setShape(0, n_outer_cells);
    data_outerGridposBlock = input_outerGridposBlock->allocate<int64_t>();
    p_vdata_outerGridposBlock = &((*data_outerGridposBlock)[0]);
  }

  auto data_tauBlock = input_tauBlock.allocate<float>();
  auto* p_vdata_tauBlock = &((*data_tauBlock)[0]);

  auto data_innerEgammaBlock = input_innerEgammaBlock.allocate<float>();
  auto* p_vdata_innerEgammaBlock = &((*data_innerEgammaBlock)[0]);

  auto data_outerEgammaBlock = input_outerEgammaBlock.allocate<float>();
  auto* p_vdata_outerEgammaBlock = &((*data_outerEgammaBlock)[0]);

  auto data_innerMuonBlock = input_innerMuonBlock.allocate<float>();
  auto* p_vdata_innerMuonBlock = &((*data_innerMuonBlock)[0]);

  auto data_outerMuonBlock = input_outerMuonBlock.allocate<float>();
  auto* p_vdata_outerMuonBlock = &((*data_outerMuonBlock)[0]);

  auto data_innerHadronBlock = input_innerHadronBlock.allocate<float>();
  auto* p_vdata_innerHadronBlock = &((*data_innerHadronBlock)[0]);

  auto data_outerHadronBlock = input_outerHadronBlock.allocate<float>();
  auto* p_vdata_outerHadronBlock = &((*data_outerHadronBlock)[0]);

  for (unsigned itau_passed = 0; itau_passed < tau_indices_.size(); ++itau_passed) {
    if (!doSplitVersion_) {
      // for the non-split version,
      // batch size is not one, needs to go to corresponding itau
      p_vdata_tauBlock = &((*data_tauBlock)[itau_passed]);

      p_vdata_innerEgammaBlock = &((*data_innerEgammaBlock)[itau_passed]);
      p_vdata_innerMuonBlock = &((*data_innerMuonBlock)[itau_passed]);
      p_vdata_innerHadronBlock = &((*data_innerHadronBlock)[itau_passed]);

      p_vdata_outerEgammaBlock = &((*data_outerEgammaBlock)[itau_passed]);
      p_vdata_outerMuonBlock = &((*data_outerMuonBlock)[itau_passed]);
      p_vdata_outerHadronBlock = &((*data_outerHadronBlock)[itau_passed]);
    }

    int tau_index = tau_indices_[itau_passed];
    const edm::RefToBase<reco::BaseTau> tauRef = taus->refAt(tau_index);
    getPredictionsV2<pat::PackedCandidate, pat::Tau>(taus->at(tau_index),
                                                     tau_index,
                                                     tauRef,
                                                     electron_collection,
                                                     muon_collection,
                                                     *pfCands,
                                                     vertices->at(0),
                                                     *rho,
                                                     tauIDs,
                                                     p_vdata_tauBlock,
                                                     p_vdata_innerEgammaBlock,
                                                     p_vdata_innerMuonBlock,
                                                     p_vdata_innerHadronBlock,
                                                     p_vdata_outerEgammaBlock,
                                                     p_vdata_outerMuonBlock,
                                                     p_vdata_outerHadronBlock,
                                                     p_vdata_innerGridposBlock,
                                                     p_vdata_outerGridposBlock);
  }

  if (doSplitVersion_) {
    // insert one more set of zeros to calculate the 'ZeroOutputTensor'
    // i.e., the output from inner/outer network when the input is zero
    // this tensor will be used to pad the core network for the cells without any particle
    p_vdata_innerEgammaBlock->insert(
        p_vdata_innerEgammaBlock->end(), dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs, 0.);
    p_vdata_innerMuonBlock->insert(
        p_vdata_innerMuonBlock->end(), dnn_inputs_2017_v2::MuonBlockInputs::NumberOfInputs, 0.);
    p_vdata_innerHadronBlock->insert(
        p_vdata_innerHadronBlock->end(), dnn_inputs_2017_v2::HadronBlockInputs::NumberOfInputs, 0.);

    p_vdata_outerEgammaBlock->insert(
        p_vdata_outerEgammaBlock->end(), dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs, 0.);
    p_vdata_outerMuonBlock->insert(
        p_vdata_outerMuonBlock->end(), dnn_inputs_2017_v2::MuonBlockInputs::NumberOfInputs, 0.);
    p_vdata_outerHadronBlock->insert(
        p_vdata_outerHadronBlock->end(), dnn_inputs_2017_v2::HadronBlockInputs::NumberOfInputs, 0.);

    // the actual number of inner cells in the event + 1
    // Note the last element of the Egamma, Muon, and Hadron Block is the zero-paddled vector
    // for retriving outputs from the inner network when the inputs are zero, which will be
    // used to paddle the inputs for the core network
    int n_inner_cells = (p_vdata_innerEgammaBlock->size() / dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs);
    input_innerEgammaBlock.setShape(0, n_inner_cells);
    input_innerMuonBlock.setShape(0, n_inner_cells);
    input_innerHadronBlock.setShape(0, n_inner_cells);

    // the actual number of outer cells in the event + 1
    // Note the last element of the Egamma, Muon, and Hadron Block is the zero-paddled vector
    // for retriving outputs from the outer network when the inputs are zero, which will be
    // used to paddle the inputs for the core network
    int n_outer_cells = (p_vdata_outerEgammaBlock->size() / dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs);
    input_outerEgammaBlock.setShape(0, n_outer_cells);
    input_outerMuonBlock.setShape(0, n_outer_cells);
    input_outerHadronBlock.setShape(0, n_outer_cells);

    // grid coordinates (i-th tau, j-th eta, k-th phi) of the inner and outer cells
    // The last element from the inner and outer network is zero-paddled vector
    // subtract it when setting the Gridpos shape
    input_innerGridposBlock->setShape(0, n_inner_cells - 1);
    input_innerGridposBlock->toServer(data_innerGridposBlock);
    input_outerGridposBlock->setShape(0, n_outer_cells - 1);
    input_outerGridposBlock->toServer(data_outerGridposBlock);
  }

  // tau
  input_tauBlock.toServer(data_tauBlock);

  input_innerEgammaBlock.toServer(data_innerEgammaBlock);
  input_innerMuonBlock.toServer(data_innerMuonBlock);
  input_innerHadronBlock.toServer(data_innerHadronBlock);

  input_outerEgammaBlock.toServer(data_outerEgammaBlock);
  input_outerMuonBlock.toServer(data_outerMuonBlock);
  input_outerHadronBlock.toServer(data_outerHadronBlock);
}

void DeepTauIdSonicProducer::produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) {
  if (tau_indices_.empty()) {
    edm::LogInfo("DeepTauIdSonicProducer") << "no tau sent to the server; skip this event in produce";
    return;
  }
  edm::Handle<TauCollection> taus;
  iEvent.getByToken(tausToken_, taus);
  const auto& output_tauval = iOutput.at("main_output/Softmax");
  const auto& outputs_tauval = output_tauval.fromServer<float>();

  // fill the taus passing the selections with the results from produce,
  //  and the taus failing the selections with 2 and -1.0
  //  (based on https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_0_pre2/RecoTauTag/RecoTau/plugins/DeepTauId.cc#L1541-L1543 )
  std::vector<std::vector<float>> pred_all(taus->size(), std::vector<float>(deep_tau::NumberOfOutputs, 0.));
  for (unsigned itau = 0; itau < taus->size(); ++itau) {
    for (int k = 0; k < deep_tau::NumberOfOutputs; ++k) {
      pred_all[itau][k] = (k == 2) ? -1.f : 2.f;
    }
  }
  for (unsigned itau_passed = 0; itau_passed < tau_indices_.size(); ++itau_passed) {
    int tau_index = tau_indices_[itau_passed];
    if (doSplitVersion_) {
      // the current mode always runs with batchSize of 1
      std::copy(outputs_tauval[0].begin() + deep_tau::NumberOfOutputs * itau_passed,
                outputs_tauval[0].begin() + deep_tau::NumberOfOutputs * (itau_passed + 1),
                pred_all[tau_index].begin());
    } else {
      std::copy(outputs_tauval[itau_passed].begin(), outputs_tauval[itau_passed].end(), pred_all[tau_index].begin());
    }

    if (debug_level >= 2) {
      for (int i = 0; i < 4; ++i) {
        std::cout << "tau index " << itau_passed << " k " << i << " pred " << pred_all[tau_index][i] << std::endl;
      }
    }
  }

  createOutputs(iEvent, pred_all, taus);
}

void DeepTauIdSonicProducer::createOutputs(edm::Event& event,
                                           const std::vector<std::vector<float>>& pred,
                                           edm::Handle<TauCollection> taus) {
  for (const auto& output_desc : outputdiscs_) {
    const WPList* working_points = nullptr;
    if (workingPoints_.find(output_desc.first) != workingPoints_.end()) {
      working_points = &workingPoints_.at(output_desc.first);
    }
    auto result = output_desc.second.get_value(taus, pred, working_points, false);
    event.put(std::move(result), output_desc.first);
  }
}

template <typename CandidateCastType, typename TauCastType>
void DeepTauIdSonicProducer::getPredictionsV2(TauCollection::const_reference& tau,
                                              const size_t tau_index,
                                              const edm::RefToBase<reco::BaseTau> tau_ref,
                                              const std::vector<pat::Electron>* electrons,
                                              const std::vector<pat::Muon>* muons,
                                              const edm::View<reco::Candidate>& pfCands,
                                              const reco::Vertex& pv,
                                              double rho,
                                              const TauFunc& tau_funcs,
                                              std::vector<float>* p_tauBlockInputs,
                                              std::vector<float>* p_egammaInnerBlockInputs,
                                              std::vector<float>* p_muonInnerBlockInputs,
                                              std::vector<float>* p_hadronInnerBlockInputs,
                                              std::vector<float>* p_egammaOuterBlockInputs,
                                              std::vector<float>* p_muonOuterBlockInputs,
                                              std::vector<float>* p_hadronOuterBlockInputs,
                                              std::vector<int64_t>* p_innerGridposInputs,
                                              std::vector<int64_t>* p_outerGridposInputs) {
  CellGrid inner_grid(dnn_inputs_2017_v2::number_of_inner_cell,
                      dnn_inputs_2017_v2::number_of_inner_cell,
                      0.02,
                      0.02,
                      disable_CellIndex_workaround_);
  CellGrid outer_grid(dnn_inputs_2017_v2::number_of_outer_cell,
                      dnn_inputs_2017_v2::number_of_outer_cell,
                      0.05,
                      0.05,
                      disable_CellIndex_workaround_);
  auto tau_casted = dynamic_cast<const TauCastType&>(tau);
  // fill in the inner and outer grids for electrons, muons, and pfCands
  fillGrids(tau_casted, *electrons, inner_grid, outer_grid);
  fillGrids(tau_casted, *muons, inner_grid, outer_grid);
  fillGrids(tau_casted, pfCands, inner_grid, outer_grid);

  unsigned tauBlockSize = p_tauBlockInputs->size();
  p_tauBlockInputs->insert(p_tauBlockInputs->end(), dnn_inputs_2017_v2::TauBlockInputs::NumberOfInputs, 0.);
  std::vector<float>::iterator tauIter = p_tauBlockInputs->begin() + tauBlockSize;

  createTauBlockInputs<CandidateCastType>(
      tau_casted, tau_index, tau_ref, pv, rho, tau_funcs, tauIter, disable_dxy_pca_);
  using namespace dnn_inputs_2017_v2;

  // egamma, muon, and hadron inner and outer inputs for the grids
  createConvFeatures<CandidateCastType>(tau_casted,
                                        tau_index,
                                        tau_ref,
                                        pv,
                                        rho,
                                        electrons,
                                        muons,
                                        pfCands,
                                        inner_grid,
                                        tau_funcs,
                                        true,
                                        p_egammaInnerBlockInputs,
                                        p_muonInnerBlockInputs,
                                        p_hadronInnerBlockInputs,
                                        p_innerGridposInputs);
  createConvFeatures<CandidateCastType>(tau_casted,
                                        tau_index,
                                        tau_ref,
                                        pv,
                                        rho,
                                        electrons,
                                        muons,
                                        pfCands,
                                        outer_grid,
                                        tau_funcs,
                                        false,
                                        p_egammaOuterBlockInputs,
                                        p_muonOuterBlockInputs,
                                        p_hadronOuterBlockInputs,
                                        p_outerGridposInputs);
}

template <typename CandidateCastType, typename TauCastType>
void DeepTauIdSonicProducer::createConvFeatures(const TauCastType& tau,
                                                const size_t tau_index,
                                                const edm::RefToBase<reco::BaseTau> tau_ref,
                                                const reco::Vertex& pv,
                                                double rho,
                                                const std::vector<pat::Electron>* electrons,
                                                const std::vector<pat::Muon>* muons,
                                                const edm::View<reco::Candidate>& pfCands,
                                                const CellGrid& grid,
                                                const TauFunc& tau_funcs,
                                                bool is_inner,
                                                std::vector<float>* p_egammaBlockInputs,
                                                std::vector<float>* p_muonBlockInputs,
                                                std::vector<float>* p_hadronBlockInputs,
                                                std::vector<int64_t>* p_GridposInputs) {
  // fill in the block inputs with zeros
  int n_cells = 0;
  if (doSplitVersion_) {
    n_cells = grid.num_valid_cells();
  } else {
    n_cells = is_inner ? (dnn_inputs_2017_v2::number_of_inner_cell * dnn_inputs_2017_v2::number_of_inner_cell)
                       : (dnn_inputs_2017_v2::number_of_outer_cell * dnn_inputs_2017_v2::number_of_outer_cell);
  }

  unsigned egammaBlockSize = p_egammaBlockInputs->size();
  p_egammaBlockInputs->insert(
      p_egammaBlockInputs->end(), n_cells * dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs, 0.);
  std::vector<float>::iterator egammaIter = p_egammaBlockInputs->begin() + egammaBlockSize;

  unsigned muonBlockSize = p_muonBlockInputs->size();
  p_muonBlockInputs->insert(
      p_muonBlockInputs->end(), n_cells * dnn_inputs_2017_v2::MuonBlockInputs::NumberOfInputs, 0.);
  std::vector<float>::iterator muonIter = p_muonBlockInputs->begin() + muonBlockSize;

  unsigned hadronBlockSize = p_hadronBlockInputs->size();
  p_hadronBlockInputs->insert(
      p_hadronBlockInputs->end(), n_cells * dnn_inputs_2017_v2::HadronBlockInputs::NumberOfInputs, 0.);
  std::vector<float>::iterator hadronIter = p_hadronBlockInputs->begin() + hadronBlockSize;

  unsigned idx = 0;
  for (int eta = -grid.maxEtaIndex(); eta <= grid.maxEtaIndex(); ++eta) {
    for (int phi = -grid.maxPhiIndex(); phi <= grid.maxPhiIndex(); ++phi) {
      if (debug_level >= 2) {
        std::cout << "processing ( eta = " << eta << ", phi = " << phi << " )" << std::endl;
      }
      const CellIndex cell_index{eta, phi};
      const int eta_index = grid.getEtaTensorIndex(cell_index);
      const int phi_index = grid.getPhiTensorIndex(cell_index);

      const auto cell_iter = grid.find(cell_index);
      if (cell_iter != grid.end()) {
        if (debug_level >= 2) {
          std::cout << " creating inputs for ( eta = " << eta << ", phi = " << phi << " ): idx = " << idx << std::endl;
        }
        const Cell& cell = cell_iter->second;
        createEgammaBlockInputs<CandidateCastType>(
            idx, tau, tau_index, tau_ref, pv, rho, electrons, pfCands, cell, tau_funcs, is_inner, egammaIter);
        createMuonBlockInputs<CandidateCastType>(
            idx, tau, tau_index, tau_ref, pv, rho, muons, pfCands, cell, tau_funcs, is_inner, muonIter);
        createHadronsBlockInputs<CandidateCastType>(idx,
                                                    tau,
                                                    tau_index,
                                                    tau_ref,
                                                    pv,
                                                    rho,
                                                    pfCands,
                                                    cell,
                                                    tau_funcs,
                                                    is_inner,
                                                    hadronIter,
                                                    disable_hcalFraction_workaround_);

        if (doSplitVersion_) {
          p_GridposInputs->push_back(tau_index);
          p_GridposInputs->push_back(eta_index);
          p_GridposInputs->push_back(phi_index);
        }
        idx += 1;
      } else {
        if (debug_level >= 2) {
          std::cout << " skipping creation of inputs, because ( eta = " << eta << ", phi = " << phi
                    << " ) is not in the grid !!" << std::endl;
        }
        if (!doSplitVersion_) {
          // for the nonSplitVersion, we need to fill in the zeros for the input tensors
          idx += 1;
        }
      }
    }
  }
}

void DeepTauIdSonicProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("taus", edm::InputTag("slimmedTaus"));
  desc.add<edm::InputTag>("pfcands", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoAll"));
  desc.add<int>("debug_level", 0);
  desc.add<bool>("disable_dxy_pca", false);
  desc.add<bool>("disable_hcalFraction_workaround", false);
  desc.add<bool>("disable_CellIndex_workaround", false);
  desc.add<bool>("doSplitVersion", true);

  desc.add<std::vector<std::string>>("VSeWP", {"-1."});
  desc.add<std::vector<std::string>>("VSmuWP", {"-1."});
  desc.add<std::vector<std::string>>("VSjetWP", {"-1."});

  desc.addUntracked<edm::InputTag>("basicTauDiscriminators", edm::InputTag("basicTauDiscriminators"));
  desc.addUntracked<edm::InputTag>("basicTauDiscriminatorsdR03", edm::InputTag("basicTauDiscriminatorsdR03"));
  desc.add<edm::InputTag>("pfTauTransverseImpactParameters", edm::InputTag("hpsPFTauTransverseImpactParameters"));

  {
    edm::ParameterSetDescription pset_Prediscriminants;
    pset_Prediscriminants.add<std::string>("BooleanOperator", "and");
    {
      edm::ParameterSetDescription psd1;
      psd1.add<double>("cut");
      psd1.add<edm::InputTag>("Producer");
      pset_Prediscriminants.addOptional<edm::ParameterSetDescription>("decayMode", psd1);
    }
    desc.add<edm::ParameterSetDescription>("Prediscriminants", pset_Prediscriminants);
  }

  descriptions.add("DeepTauIdSonicProducer", desc);
}

DEFINE_FWK_MODULE(DeepTauIdSonicProducer);
