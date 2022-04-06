/*
 * \class DeepTauId
 *
 * Tau identification using Deep NN.
 *
 * \author Konstantin Androsov, INFN Pisa
 *         Christian Veelken, Tallinn
 */

#include "RecoTauTag/RecoTau/interface/DeepTauBase.h"
#include "RecoTauTag/RecoTau/interface/DeepTauHelper.h"

namespace deep_tau {
  constexpr int NumberOfOutputs = 4;
}

using namespace deep_tau_2017;
using namespace deeptau_helper;

using bd = deep_tau::DeepTauBase::BasicDiscriminator;
const std::map<bd, std::string> deep_tau::DeepTauBase::stringFromDiscriminator_{
    {bd::ChargedIsoPtSum, "ChargedIsoPtSum"},
    {bd::NeutralIsoPtSum, "NeutralIsoPtSum"},
    {bd::NeutralIsoPtSumWeight, "NeutralIsoPtSumWeight"},
    {bd::FootprintCorrection, "TauFootprintCorrection"},
    {bd::PhotonPtSumOutsideSignalCone, "PhotonPtSumOutsideSignalCone"},
    {bd::PUcorrPtSum, "PUcorrPtSum"}};
const std::vector<bd> deep_tau::DeepTauBase::requiredBasicDiscriminators_ = {bd::ChargedIsoPtSum,
                                                                             bd::NeutralIsoPtSum,
                                                                             bd::NeutralIsoPtSumWeight,
                                                                             bd::PhotonPtSumOutsideSignalCone,
                                                                             bd::PUcorrPtSum};
const std::vector<bd> deep_tau::DeepTauBase::requiredBasicDiscriminatorsdR03_ = {bd::ChargedIsoPtSum,
                                                                                 bd::NeutralIsoPtSum,
                                                                                 bd::NeutralIsoPtSumWeight,
                                                                                 bd::PhotonPtSumOutsideSignalCone,
                                                                                 bd::FootprintCorrection};

class DeepTauId : public deep_tau::DeepTauBase {
public:
  const std::map<BasicDiscriminator, size_t> matchDiscriminatorIndices(
      edm::Event& event,
      edm::EDGetTokenT<reco::TauDiscriminatorContainer> discriminatorContainerToken,
      std::vector<BasicDiscriminator> requiredDiscr) {
    std::map<std::string, size_t> discrIndexMapStr;
    auto const aHandle = event.getHandle(discriminatorContainerToken);
    auto const aProv = aHandle.provenance();
    if (aProv == nullptr)
      aHandle.whyFailed()->raise();
    const auto& psetsFromProvenance = edm::parameterSet(aProv->stable(), event.processHistory());
    auto const idlist = psetsFromProvenance.getParameter<std::vector<edm::ParameterSet>>("IDdefinitions");
    for (size_t j = 0; j < idlist.size(); ++j) {
      std::string idname = idlist[j].getParameter<std::string>("IDname");
      if (discrIndexMapStr.count(idname)) {
        throw cms::Exception("DeepTauId")
            << "basic discriminator " << idname << " appears more than once in the input.";
      }
      discrIndexMapStr[idname] = j;
    }

    //translate to a map of <BasicDiscriminator, index> and check if all discriminators are present
    std::map<BasicDiscriminator, size_t> discrIndexMap;
    for (size_t i = 0; i < requiredDiscr.size(); i++) {
      if (discrIndexMapStr.find(stringFromDiscriminator_.at(requiredDiscr[i])) == discrIndexMapStr.end())
        throw cms::Exception("DeepTauId") << "Basic Discriminator " << stringFromDiscriminator_.at(requiredDiscr[i])
                                          << " was not provided in the config file.";
      else
        discrIndexMap[requiredDiscr[i]] = discrIndexMapStr[stringFromDiscriminator_.at(requiredDiscr[i])];
    }
    return discrIndexMap;
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
    desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
    desc.add<edm::InputTag>("taus", edm::InputTag("slimmedTaus"));
    desc.add<edm::InputTag>("pfcands", edm::InputTag("packedPFCandidates"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
    desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoAll"));
    desc.add<std::vector<std::string>>("graph_file",
                                       {"RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v2p6_e6.pb"});
    desc.add<bool>("mem_mapped", false);
    desc.add<unsigned>("version", 2);
    desc.add<int>("debug_level", 0);
    desc.add<bool>("disable_dxy_pca", false);
    desc.add<bool>("disable_hcalFraction_workaround", false);
    desc.add<bool>("disable_CellIndex_workaround", false);
    desc.add<bool>("save_inputs", false);
    desc.add<bool>("is_online", false);

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

    descriptions.add("DeepTau", desc);
  }

public:
  explicit DeepTauId(const edm::ParameterSet& cfg, const deep_tau::DeepTauCache* cache)
      : DeepTauBase(cfg, GetOutputs(), cache),
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
        version_(cfg.getParameter<unsigned>("version")),
        debug_level(cfg.getParameter<int>("debug_level")),
        disable_dxy_pca_(cfg.getParameter<bool>("disable_dxy_pca")),
        disable_hcalFraction_workaround_(cfg.getParameter<bool>("disable_hcalFraction_workaround")),
        disable_CellIndex_workaround_(cfg.getParameter<bool>("disable_CellIndex_workaround")),
        save_inputs_(cfg.getParameter<bool>("save_inputs")),
        json_file_(nullptr),
        file_counter_(0) {
    if (version_ == 1) {
      input_layer_ = cache_->getGraph().node(0).name();
      output_layer_ = cache_->getGraph().node(cache_->getGraph().node_size() - 1).name();
      const auto& shape = cache_->getGraph().node(0).attr().at("shape").shape();
      if (shape.dim(1).size() != dnn_inputs_2017v1::NumberOfInputs)
        throw cms::Exception("DeepTauId")
            << "number of inputs does not match the expected inputs for the given version";
    } else if (version_ == 2) {
      tauBlockTensor_ = std::make_unique<tensorflow::Tensor>(
          tensorflow::DT_FLOAT, tensorflow::TensorShape{1, dnn_inputs_2017_v2::TauBlockInputs::NumberOfInputs});
      for (size_t n = 0; n < 2; ++n) {
        const bool is_inner = n == 0;
        const auto n_cells =
            is_inner ? dnn_inputs_2017_v2::number_of_inner_cell : dnn_inputs_2017_v2::number_of_outer_cell;
        eGammaTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
            tensorflow::DT_FLOAT,
            tensorflow::TensorShape{1, 1, 1, dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs});
        muonTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
            tensorflow::DT_FLOAT,
            tensorflow::TensorShape{1, 1, 1, dnn_inputs_2017_v2::MuonBlockInputs::NumberOfInputs});
        hadronsTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
            tensorflow::DT_FLOAT,
            tensorflow::TensorShape{1, 1, 1, dnn_inputs_2017_v2::HadronBlockInputs::NumberOfInputs});
        convTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
            tensorflow::DT_FLOAT,
            tensorflow::TensorShape{1, n_cells, n_cells, dnn_inputs_2017_v2::number_of_conv_features});
        zeroOutputTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
            tensorflow::DT_FLOAT, tensorflow::TensorShape{1, 1, 1, dnn_inputs_2017_v2::number_of_conv_features});

        eGammaTensor_[is_inner]->flat<float>().setZero();
        muonTensor_[is_inner]->flat<float>().setZero();
        hadronsTensor_[is_inner]->flat<float>().setZero();

        setCellConvFeatures(*zeroOutputTensor_[is_inner], getPartialPredictions(is_inner), 0, 0, 0);
      }
    } else {
      throw cms::Exception("DeepTauId") << "version " << version_ << " is not supported.";
    }
  }

  static std::unique_ptr<deep_tau::DeepTauCache> initializeGlobalCache(const edm::ParameterSet& cfg) {
    return DeepTauBase::initializeGlobalCache(cfg);
  }

  static void globalEndJob(const deep_tau::DeepTauCache* cache_) { return DeepTauBase::globalEndJob(cache_); }

private:
  inline void checkInputs(const tensorflow::Tensor& inputs,
                          const std::string& block_name,
                          int n_inputs,
                          const CellGrid* grid = nullptr) const {
    if (debug_level >= 1) {
      std::cout << "<checkInputs>: block_name = " << block_name << std::endl;
      if (block_name == "input_tau") {
        for (int input_index = 0; input_index < n_inputs; ++input_index) {
          float input = inputs.matrix<float>()(0, input_index);
          if (edm::isNotFinite(input)) {
            throw cms::Exception("DeepTauId")
                << "in the " << block_name
                << ", input is not finite, i.e. infinite or NaN, for input_index = " << input_index;
          }
          if (debug_level >= 2) {
            std::cout << block_name << "[var = " << input_index << "] = " << std::setprecision(5) << std::fixed << input
                      << std::endl;
          }
        }
      } else {
        assert(grid);
        int n_eta, n_phi;
        if (block_name.find("input_inner") != std::string::npos) {
          n_eta = 5;
          n_phi = 5;
        } else if (block_name.find("input_outer") != std::string::npos) {
          n_eta = 10;
          n_phi = 10;
        } else
          assert(0);
        int eta_phi_index = 0;
        for (int eta = -n_eta; eta <= n_eta; ++eta) {
          for (int phi = -n_phi; phi <= n_phi; ++phi) {
            const CellIndex cell_index{eta, phi};
            const auto cell_iter = grid->find(cell_index);
            if (cell_iter != grid->end()) {
              for (int input_index = 0; input_index < n_inputs; ++input_index) {
                float input = inputs.tensor<float, 4>()(eta_phi_index, 0, 0, input_index);
                if (edm::isNotFinite(input)) {
                  throw cms::Exception("DeepTauId")
                      << "in the " << block_name << ", input is not finite, i.e. infinite or NaN, for eta = " << eta
                      << ", phi = " << phi << ", input_index = " << input_index;
                }
                if (debug_level >= 2) {
                  std::cout << block_name << "[eta = " << eta << "][phi = " << phi << "][var = " << input_index
                            << "] = " << std::setprecision(5) << std::fixed << input << std::endl;
                }
              }
              eta_phi_index += 1;
            }
          }
        }
      }
    }
  }

  inline void saveInputs(const tensorflow::Tensor& inputs,
                         const std::string& block_name,
                         int n_inputs,
                         const CellGrid* grid = nullptr) {
    if (debug_level >= 1) {
      std::cout << "<saveInputs>: block_name = " << block_name << std::endl;
    }
    if (!is_first_block_)
      (*json_file_) << ", ";
    (*json_file_) << "\"" << block_name << "\": [";
    if (block_name == "input_tau") {
      for (int input_index = 0; input_index < n_inputs; ++input_index) {
        float input = inputs.matrix<float>()(0, input_index);
        if (input_index != 0)
          (*json_file_) << ", ";
        (*json_file_) << input;
      }
    } else {
      assert(grid);
      int n_eta, n_phi;
      if (block_name.find("input_inner") != std::string::npos) {
        n_eta = 5;
        n_phi = 5;
      } else if (block_name.find("input_outer") != std::string::npos) {
        n_eta = 10;
        n_phi = 10;
      } else
        assert(0);
      int eta_phi_index = 0;
      for (int eta = -n_eta; eta <= n_eta; ++eta) {
        if (eta != -n_eta)
          (*json_file_) << ", ";
        (*json_file_) << "[";
        for (int phi = -n_phi; phi <= n_phi; ++phi) {
          if (phi != -n_phi)
            (*json_file_) << ", ";
          (*json_file_) << "[";
          const CellIndex cell_index{eta, phi};
          const auto cell_iter = grid->find(cell_index);
          for (int input_index = 0; input_index < n_inputs; ++input_index) {
            float input = 0.;
            if (cell_iter != grid->end()) {
              input = inputs.tensor<float, 4>()(eta_phi_index, 0, 0, input_index);
            }
            if (input_index != 0)
              (*json_file_) << ", ";
            (*json_file_) << input;
          }
          if (cell_iter != grid->end()) {
            eta_phi_index += 1;
          }
          (*json_file_) << "]";
        }
        (*json_file_) << "]";
      }
    }
    (*json_file_) << "]";
    is_first_block_ = false;
  }

private:
  tensorflow::Tensor getPredictions(edm::Event& event, edm::Handle<TauCollection> taus) override {
    // Empty dummy vectors
    const std::vector<pat::Electron> electron_collection_default;
    const std::vector<pat::Muon> muon_collection_default;
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

    if (!is_online_) {
      electron_collection = &event.get(electrons_token_);
      muon_collection = &event.get(muons_token_);
      pfTauTransverseImpactParameters = &pfTauTransverseImpactParameters_default;
      basicTauDiscriminators = &basicTauDiscriminators_default;
      basicTauDiscriminatorsdR03 = &basicTauDiscriminatorsdR03_default;
    } else {
      electron_collection = &electron_collection_default;
      muon_collection = &muon_collection_default;
      pfTauTransverseImpactParameters = &event.get(pfTauTransverseImpactParameters_token_);
      basicTauDiscriminators = &event.get(basicTauDiscriminators_inputToken_);
      basicTauDiscriminatorsdR03 = &event.get(basicTauDiscriminatorsdR03_inputToken_);

      // Get indices for discriminators
      if (!discrIndicesMapped_) {
        basicDiscrIndexMap_ =
            matchDiscriminatorIndices(event, basicTauDiscriminators_inputToken_, requiredBasicDiscriminators_);
        basicDiscrdR03IndexMap_ =
            matchDiscriminatorIndices(event, basicTauDiscriminatorsdR03_inputToken_, requiredBasicDiscriminatorsdR03_);
        discrIndicesMapped_ = true;
      }
    }

    TauFunc tauIDs = {basicTauDiscriminators,
                      basicTauDiscriminatorsdR03,
                      pfTauTransverseImpactParameters,
                      basicDiscrIndexMap_,
                      basicDiscrdR03IndexMap_};

    edm::Handle<edm::View<reco::Candidate>> pfCands;
    event.getByToken(pfcandToken_, pfCands);

    edm::Handle<reco::VertexCollection> vertices;
    event.getByToken(vtxToken_, vertices);

    edm::Handle<double> rho;
    event.getByToken(rho_token_, rho);

    tensorflow::Tensor predictions(tensorflow::DT_FLOAT, {static_cast<int>(taus->size()), deep_tau::NumberOfOutputs});

    for (size_t tau_index = 0; tau_index < taus->size(); ++tau_index) {
      const edm::RefToBase<reco::BaseTau> tauRef = taus->refAt(tau_index);

      std::vector<tensorflow::Tensor> pred_vector;

      bool passesPrediscriminants;
      if (is_online_) {
        passesPrediscriminants = tauIDs.passPrediscriminants<std::vector<TauDiscInfo<reco::PFTauDiscriminator>>>(
            recoPrediscriminants_, andPrediscriminants_, tauRef);
      } else {
        passesPrediscriminants = tauIDs.passPrediscriminants<std::vector<TauDiscInfo<pat::PATTauDiscriminator>>>(
            patPrediscriminants_, andPrediscriminants_, tauRef);
      }

      if (passesPrediscriminants) {
        if (version_ == 1) {
          if (is_online_)
            getPredictionsV1<reco::PFCandidate, reco::PFTau>(
                taus->at(tau_index), tau_index, tauRef, electron_collection, muon_collection, pred_vector, tauIDs);
          else
            getPredictionsV1<pat::PackedCandidate, pat::Tau>(
                taus->at(tau_index), tau_index, tauRef, electron_collection, muon_collection, pred_vector, tauIDs);
        } else if (version_ == 2) {
          if (is_online_) {
            getPredictionsV2<reco::PFCandidate, reco::PFTau>(taus->at(tau_index),
                                                             tau_index,
                                                             tauRef,
                                                             electron_collection,
                                                             muon_collection,
                                                             *pfCands,
                                                             vertices->at(0),
                                                             *rho,
                                                             pred_vector,
                                                             tauIDs);
          } else
            getPredictionsV2<pat::PackedCandidate, pat::Tau>(taus->at(tau_index),
                                                             tau_index,
                                                             tauRef,
                                                             electron_collection,
                                                             muon_collection,
                                                             *pfCands,
                                                             vertices->at(0),
                                                             *rho,
                                                             pred_vector,
                                                             tauIDs);
        } else {
          throw cms::Exception("DeepTauId") << "version " << version_ << " is not supported.";
        }

        for (int k = 0; k < deep_tau::NumberOfOutputs; ++k) {
          const float pred = pred_vector[0].flat<float>()(k);
          if (!(pred >= 0 && pred <= 1))
            throw cms::Exception("DeepTauId")
                << "invalid prediction = " << pred << " for tau_index = " << tau_index << ", pred_index = " << k;
          predictions.matrix<float>()(tau_index, k) = pred;
        }
      } else {
        // This else statement was added as a part of the DeepTau@HLT development. It does not affect the current state
        // of offline DeepTauId code as there the preselection is not used (it was added in the DeepTau@HLT). It returns
        // default values for deepTau score if the preselection failed. Before this statement the values given for this tau
        // were random. k == 2 corresponds to the tau score and all other k values to e, mu and jets. By defining in this way
        // the final score is -1.
        for (int k = 0; k < deep_tau::NumberOfOutputs; ++k) {
          predictions.matrix<float>()(tau_index, k) = (k == 2) ? -1.f : 2.f;
        }
      }
    }
    return predictions;
  }

  template <typename CandidateCastType, typename TauCastType>
  void getPredictionsV1(TauCollection::const_reference& tau,
                        const size_t tau_index,
                        const edm::RefToBase<reco::BaseTau> tau_ref,
                        const std::vector<pat::Electron>* electrons,
                        const std::vector<pat::Muon>* muons,
                        std::vector<tensorflow::Tensor>& pred_vector,
                        TauFunc tau_funcs) {
    const tensorflow::Tensor& inputs = createInputsV1<dnn_inputs_2017v1, const CandidateCastType>(
        dynamic_cast<const TauCastType&>(tau), tau_index, tau_ref, electrons, muons, tau_funcs);
    tensorflow::run(&(cache_->getSession()), {{input_layer_, inputs}}, {output_layer_}, &pred_vector);
  }

  template <typename CandidateCastType, typename TauCastType>
  void getPredictionsV2(TauCollection::const_reference& tau,
                        const size_t tau_index,
                        const edm::RefToBase<reco::BaseTau> tau_ref,
                        const std::vector<pat::Electron>* electrons,
                        const std::vector<pat::Muon>* muons,
                        const edm::View<reco::Candidate>& pfCands,
                        const reco::Vertex& pv,
                        double rho,
                        std::vector<tensorflow::Tensor>& pred_vector,
                        const TauFunc& tau_funcs) {
    if (debug_level >= 2) {
      std::cout << "<DeepTauId::getPredictionsV2 (moduleLabel = " << moduleDescription().moduleLabel()
                << ")>:" << std::endl;
      std::cout << " tau: pT = " << tau.pt() << ", eta = " << tau.eta() << ", phi = " << tau.phi() << std::endl;
    }
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
    fillGrids(dynamic_cast<const TauCastType&>(tau), *electrons, inner_grid, outer_grid);
    fillGrids(dynamic_cast<const TauCastType&>(tau), *muons, inner_grid, outer_grid);
    fillGrids(dynamic_cast<const TauCastType&>(tau), pfCands, inner_grid, outer_grid);

    tauBlockTensor_->flat<float>().setZero();
    createTauBlockInputs<CandidateCastType>(dynamic_cast<const TauCastType&>(tau),
                                            tau_index,
                                            tau_ref,
                                            pv,
                                            rho,
                                            tau_funcs,
                                            *tauBlockTensor_,
                                            disable_dxy_pca_);
    using namespace dnn_inputs_2017_v2;
    checkInputs(*tauBlockTensor_, "input_tau", TauBlockInputs::NumberOfInputs);
    createConvFeatures<CandidateCastType>(dynamic_cast<const TauCastType&>(tau),
                                          tau_index,
                                          tau_ref,
                                          pv,
                                          rho,
                                          electrons,
                                          muons,
                                          pfCands,
                                          inner_grid,
                                          tau_funcs,
                                          true);
    checkInputs(*eGammaTensor_[true], "input_inner_egamma", EgammaBlockInputs::NumberOfInputs, &inner_grid);
    checkInputs(*muonTensor_[true], "input_inner_muon", MuonBlockInputs::NumberOfInputs, &inner_grid);
    checkInputs(*hadronsTensor_[true], "input_inner_hadrons", HadronBlockInputs::NumberOfInputs, &inner_grid);
    createConvFeatures<CandidateCastType>(dynamic_cast<const TauCastType&>(tau),
                                          tau_index,
                                          tau_ref,
                                          pv,
                                          rho,
                                          electrons,
                                          muons,
                                          pfCands,
                                          outer_grid,
                                          tau_funcs,
                                          false);
    checkInputs(*eGammaTensor_[false], "input_outer_egamma", EgammaBlockInputs::NumberOfInputs, &outer_grid);
    checkInputs(*muonTensor_[false], "input_outer_muon", MuonBlockInputs::NumberOfInputs, &outer_grid);
    checkInputs(*hadronsTensor_[false], "input_outer_hadrons", HadronBlockInputs::NumberOfInputs, &outer_grid);

    if (save_inputs_) {
      std::string json_file_name = "DeepTauId_" + std::to_string(file_counter_) + ".json";
      json_file_ = new std::ofstream(json_file_name.data());
      is_first_block_ = true;
      (*json_file_) << "{";
      saveInputs(*tauBlockTensor_, "input_tau", dnn_inputs_2017_v2::TauBlockInputs::NumberOfInputs);
      saveInputs(*eGammaTensor_[true],
                 "input_inner_egamma",
                 dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs,
                 &inner_grid);
      saveInputs(
          *muonTensor_[true], "input_inner_muon", dnn_inputs_2017_v2::MuonBlockInputs::NumberOfInputs, &inner_grid);
      saveInputs(*hadronsTensor_[true],
                 "input_inner_hadrons",
                 dnn_inputs_2017_v2::HadronBlockInputs::NumberOfInputs,
                 &inner_grid);
      saveInputs(*eGammaTensor_[false],
                 "input_outer_egamma",
                 dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs,
                 &outer_grid);
      saveInputs(
          *muonTensor_[false], "input_outer_muon", dnn_inputs_2017_v2::MuonBlockInputs::NumberOfInputs, &outer_grid);
      saveInputs(*hadronsTensor_[false],
                 "input_outer_hadrons",
                 dnn_inputs_2017_v2::HadronBlockInputs::NumberOfInputs,
                 &outer_grid);
      (*json_file_) << "}";
      delete json_file_;
      ++file_counter_;
    }

    tensorflow::run(&(cache_->getSession("core")),
                    {{"input_tau", *tauBlockTensor_},
                     {"input_inner", *convTensor_.at(true)},
                     {"input_outer", *convTensor_.at(false)}},
                    {"main_output/Softmax"},
                    &pred_vector);
    if (debug_level >= 1) {
      std::cout << "output = { ";
      for (int idx = 0; idx < deep_tau::NumberOfOutputs; ++idx) {
        if (idx > 0)
          std::cout << ", ";
        std::string label;
        if (idx == 0)
          label = "e";
        else if (idx == 1)
          label = "mu";
        else if (idx == 2)
          label = "tau";
        else if (idx == 3)
          label = "jet";
        else
          assert(0);
        std::cout << label << " = " << pred_vector[0].flat<float>()(idx);
      }
      std::cout << " }" << std::endl;
    }
  }

  tensorflow::Tensor getPartialPredictions(bool is_inner) {
    std::vector<tensorflow::Tensor> pred_vector;
    if (is_inner) {
      tensorflow::run(&(cache_->getSession("inner")),
                      {
                          {"input_inner_egamma", *eGammaTensor_.at(is_inner)},
                          {"input_inner_muon", *muonTensor_.at(is_inner)},
                          {"input_inner_hadrons", *hadronsTensor_.at(is_inner)},
                      },
                      {"inner_all_dropout_4/Identity"},
                      &pred_vector);
    } else {
      tensorflow::run(&(cache_->getSession("outer")),
                      {
                          {"input_outer_egamma", *eGammaTensor_.at(is_inner)},
                          {"input_outer_muon", *muonTensor_.at(is_inner)},
                          {"input_outer_hadrons", *hadronsTensor_.at(is_inner)},
                      },
                      {"outer_all_dropout_4/Identity"},
                      &pred_vector);
    }
    return pred_vector.at(0);
  }

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
                          bool is_inner) {
    if (debug_level >= 2) {
      std::cout << "<DeepTauId::createConvFeatures (is_inner = " << is_inner << ")>:" << std::endl;
    }
    tensorflow::Tensor& convTensor = *convTensor_.at(is_inner);
    eGammaTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
        tensorflow::DT_FLOAT,
        tensorflow::TensorShape{
            (long long int)grid.num_valid_cells(), 1, 1, dnn_inputs_2017_v2::EgammaBlockInputs::NumberOfInputs});
    muonTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
        tensorflow::DT_FLOAT,
        tensorflow::TensorShape{
            (long long int)grid.num_valid_cells(), 1, 1, dnn_inputs_2017_v2::MuonBlockInputs::NumberOfInputs});
    hadronsTensor_[is_inner] = std::make_unique<tensorflow::Tensor>(
        tensorflow::DT_FLOAT,
        tensorflow::TensorShape{
            (long long int)grid.num_valid_cells(), 1, 1, dnn_inputs_2017_v2::HadronBlockInputs::NumberOfInputs});

    eGammaTensor_[is_inner]->flat<float>().setZero();
    muonTensor_[is_inner]->flat<float>().setZero();
    hadronsTensor_[is_inner]->flat<float>().setZero();

    unsigned idx = 0;
    for (int eta = -grid.maxEtaIndex(); eta <= grid.maxEtaIndex(); ++eta) {
      for (int phi = -grid.maxPhiIndex(); phi <= grid.maxPhiIndex(); ++phi) {
        if (debug_level >= 2) {
          std::cout << "processing ( eta = " << eta << ", phi = " << phi << " )" << std::endl;
        }
        const CellIndex cell_index{eta, phi};
        const auto cell_iter = grid.find(cell_index);
        if (cell_iter != grid.end()) {
          if (debug_level >= 2) {
            std::cout << " creating inputs for ( eta = " << eta << ", phi = " << phi << " ): idx = " << idx
                      << std::endl;
          }
          const Cell& cell = cell_iter->second;
          createEgammaBlockInputs<CandidateCastType>(idx,
                                                     tau,
                                                     tau_index,
                                                     tau_ref,
                                                     pv,
                                                     rho,
                                                     electrons,
                                                     pfCands,
                                                     cell,
                                                     tau_funcs,
                                                     is_inner,
                                                     *eGammaTensor_.at(is_inner));
          createMuonBlockInputs<CandidateCastType>(idx,
                                                   tau,
                                                   tau_index,
                                                   tau_ref,
                                                   pv,
                                                   rho,
                                                   muons,
                                                   pfCands,
                                                   cell,
                                                   tau_funcs,
                                                   is_inner,
                                                   *muonTensor_.at(is_inner));
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
                                                      *hadronsTensor_.at(is_inner),
                                                      disable_hcalFraction_workaround_);
          idx += 1;
        } else {
          if (debug_level >= 2) {
            std::cout << " skipping creation of inputs, because ( eta = " << eta << ", phi = " << phi
                      << " ) is not in the grid !!" << std::endl;
          }
        }
      }
    }

    const auto predTensor = getPartialPredictions(is_inner);
    idx = 0;
    for (int eta = -grid.maxEtaIndex(); eta <= grid.maxEtaIndex(); ++eta) {
      for (int phi = -grid.maxPhiIndex(); phi <= grid.maxPhiIndex(); ++phi) {
        const CellIndex cell_index{eta, phi};
        const int eta_index = grid.getEtaTensorIndex(cell_index);
        const int phi_index = grid.getPhiTensorIndex(cell_index);

        const auto cell_iter = grid.find(cell_index);
        if (cell_iter != grid.end()) {
          setCellConvFeatures(convTensor, predTensor, idx, eta_index, phi_index);
          idx += 1;
        } else {
          setCellConvFeatures(convTensor, *zeroOutputTensor_[is_inner], 0, eta_index, phi_index);
        }
      }
    }
  }

  void setCellConvFeatures(tensorflow::Tensor& convTensor,
                           const tensorflow::Tensor& features,
                           unsigned batch_idx,
                           int eta_index,
                           int phi_index) {
    for (int n = 0; n < dnn_inputs_2017_v2::number_of_conv_features; ++n) {
      convTensor.tensor<float, 4>()(0, eta_index, phi_index, n) = features.tensor<float, 4>()(batch_idx, 0, 0, n);
    }
  }

  template <typename dnn, typename CandidateCastType, typename TauCastType>
  tensorflow::Tensor createInputsV1(const TauCastType& tau,
                                    const size_t tau_index,
                                    const edm::RefToBase<reco::BaseTau> tau_ref,
                                    const std::vector<pat::Electron>* electrons,
                                    const std::vector<pat::Muon>* muons,
                                    TauFunc tau_funcs) const {
    static constexpr bool check_all_set = false;
    static constexpr float default_value_for_set_check = -42;

    tensorflow::Tensor inputs(tensorflow::DT_FLOAT, {1, dnn_inputs_2017v1::NumberOfInputs});
    const auto& get = [&](int var_index) -> float& { return inputs.matrix<float>()(0, var_index); };
    auto leadChargedHadrCand = dynamic_cast<const CandidateCastType*>(tau.leadChargedHadrCand().get());

    if (check_all_set) {
      for (int var_index = 0; var_index < dnn::NumberOfInputs; ++var_index) {
        get(var_index) = default_value_for_set_check;
      }
    }

    get(dnn::pt) = tau.p4().pt();
    get(dnn::eta) = tau.p4().eta();
    get(dnn::mass) = tau.p4().mass();
    get(dnn::decayMode) = tau.decayMode();
    get(dnn::chargedIsoPtSum) = tau_funcs.getChargedIsoPtSum(tau, tau_ref);
    get(dnn::neutralIsoPtSum) = tau_funcs.getNeutralIsoPtSum(tau, tau_ref);
    get(dnn::neutralIsoPtSumWeight) = tau_funcs.getNeutralIsoPtSumWeight(tau, tau_ref);
    get(dnn::photonPtSumOutsideSignalCone) = tau_funcs.getPhotonPtSumOutsideSignalCone(tau, tau_ref);
    get(dnn::puCorrPtSum) = tau_funcs.getPuCorrPtSum(tau, tau_ref);
    get(dnn::dxy) = tau_funcs.getdxy(tau, tau_index);
    get(dnn::dxy_sig) = tau_funcs.getdxySig(tau, tau_index);
    get(dnn::dz) = leadChargedHadrCand ? candFunc::getTauDz(*leadChargedHadrCand) : default_value;
    get(dnn::ip3d) = tau_funcs.getip3d(tau, tau_index);
    get(dnn::ip3d_sig) = tau_funcs.getip3dSig(tau, tau_index);
    get(dnn::hasSecondaryVertex) = tau_funcs.getHasSecondaryVertex(tau, tau_index);
    get(dnn::flightLength_r) = tau_funcs.getFlightLength(tau, tau_index).R();
    get(dnn::flightLength_dEta) = dEta(tau_funcs.getFlightLength(tau, tau_index), tau.p4());
    get(dnn::flightLength_dPhi) = dPhi(tau_funcs.getFlightLength(tau, tau_index), tau.p4());
    get(dnn::flightLength_sig) = tau_funcs.getFlightLengthSig(tau, tau_index);
    get(dnn::leadChargedHadrCand_pt) = leadChargedHadrCand ? leadChargedHadrCand->p4().Pt() : default_value;
    get(dnn::leadChargedHadrCand_dEta) =
        leadChargedHadrCand ? dEta(leadChargedHadrCand->p4(), tau.p4()) : default_value;
    get(dnn::leadChargedHadrCand_dPhi) =
        leadChargedHadrCand ? dPhi(leadChargedHadrCand->p4(), tau.p4()) : default_value;
    get(dnn::leadChargedHadrCand_mass) = leadChargedHadrCand ? leadChargedHadrCand->p4().mass() : default_value;
    get(dnn::pt_weighted_deta_strip) = reco::tau::pt_weighted_deta_strip(tau, tau.decayMode());
    get(dnn::pt_weighted_dphi_strip) = reco::tau::pt_weighted_dphi_strip(tau, tau.decayMode());
    get(dnn::pt_weighted_dr_signal) = reco::tau::pt_weighted_dr_signal(tau, tau.decayMode());
    get(dnn::pt_weighted_dr_iso) = reco::tau::pt_weighted_dr_iso(tau, tau.decayMode());
    get(dnn::leadingTrackNormChi2) = tau_funcs.getLeadingTrackNormChi2(tau);
    get(dnn::e_ratio) = reco::tau::eratio(tau);
    get(dnn::gj_angle_diff) = calculateGottfriedJacksonAngleDifference(tau, tau_index, tau_funcs);
    get(dnn::n_photons) = reco::tau::n_photons_total(tau);
    get(dnn::emFraction) = tau_funcs.getEmFraction(tau);
    get(dnn::has_gsf_track) = leadChargedHadrCand && std::abs(leadChargedHadrCand->pdgId()) == 11;
    get(dnn::inside_ecal_crack) = isInEcalCrack(tau.p4().Eta());
    auto gsf_ele = findMatchedElectron(tau, electrons, 0.3);
    get(dnn::gsf_ele_matched) = gsf_ele != nullptr;
    get(dnn::gsf_ele_pt) = gsf_ele != nullptr ? gsf_ele->p4().Pt() : default_value;
    get(dnn::gsf_ele_dEta) = gsf_ele != nullptr ? dEta(gsf_ele->p4(), tau.p4()) : default_value;
    get(dnn::gsf_ele_dPhi) = gsf_ele != nullptr ? dPhi(gsf_ele->p4(), tau.p4()) : default_value;
    get(dnn::gsf_ele_mass) = gsf_ele != nullptr ? gsf_ele->p4().mass() : default_value;
    calculateElectronClusterVars(gsf_ele, get(dnn::gsf_ele_Ee), get(dnn::gsf_ele_Egamma));
    get(dnn::gsf_ele_Pin) = gsf_ele != nullptr ? gsf_ele->trackMomentumAtVtx().R() : default_value;
    get(dnn::gsf_ele_Pout) = gsf_ele != nullptr ? gsf_ele->trackMomentumOut().R() : default_value;
    get(dnn::gsf_ele_EtotOverPin) = get(dnn::gsf_ele_Pin) > 0
                                        ? (get(dnn::gsf_ele_Ee) + get(dnn::gsf_ele_Egamma)) / get(dnn::gsf_ele_Pin)
                                        : default_value;
    get(dnn::gsf_ele_Eecal) = gsf_ele != nullptr ? gsf_ele->ecalEnergy() : default_value;
    get(dnn::gsf_ele_dEta_SeedClusterTrackAtCalo) =
        gsf_ele != nullptr ? gsf_ele->deltaEtaSeedClusterTrackAtCalo() : default_value;
    get(dnn::gsf_ele_dPhi_SeedClusterTrackAtCalo) =
        gsf_ele != nullptr ? gsf_ele->deltaPhiSeedClusterTrackAtCalo() : default_value;
    get(dnn::gsf_ele_mvaIn_sigmaEtaEta) = gsf_ele != nullptr ? gsf_ele->mvaInput().sigmaEtaEta : default_value;
    get(dnn::gsf_ele_mvaIn_hadEnergy) = gsf_ele != nullptr ? gsf_ele->mvaInput().hadEnergy : default_value;
    get(dnn::gsf_ele_mvaIn_deltaEta) = gsf_ele != nullptr ? gsf_ele->mvaInput().deltaEta : default_value;

    get(dnn::gsf_ele_Chi2NormGSF) = default_value;
    get(dnn::gsf_ele_GSFNumHits) = default_value;
    get(dnn::gsf_ele_GSFTrackResol) = default_value;
    get(dnn::gsf_ele_GSFTracklnPt) = default_value;
    if (gsf_ele != nullptr && gsf_ele->gsfTrack().isNonnull()) {
      get(dnn::gsf_ele_Chi2NormGSF) = gsf_ele->gsfTrack()->normalizedChi2();
      get(dnn::gsf_ele_GSFNumHits) = gsf_ele->gsfTrack()->numberOfValidHits();
      if (gsf_ele->gsfTrack()->pt() > 0) {
        get(dnn::gsf_ele_GSFTrackResol) = gsf_ele->gsfTrack()->ptError() / gsf_ele->gsfTrack()->pt();
        get(dnn::gsf_ele_GSFTracklnPt) = std::log10(gsf_ele->gsfTrack()->pt());
      }
    }

    get(dnn::gsf_ele_Chi2NormKF) = default_value;
    get(dnn::gsf_ele_KFNumHits) = default_value;
    if (gsf_ele != nullptr && gsf_ele->closestCtfTrackRef().isNonnull()) {
      get(dnn::gsf_ele_Chi2NormKF) = gsf_ele->closestCtfTrackRef()->normalizedChi2();
      get(dnn::gsf_ele_KFNumHits) = gsf_ele->closestCtfTrackRef()->numberOfValidHits();
    }
    get(dnn::leadChargedCand_etaAtEcalEntrance) = tau_funcs.getEtaAtEcalEntrance(tau);
    get(dnn::leadChargedCand_pt) = leadChargedHadrCand->pt();

    get(dnn::leadChargedHadrCand_HoP) = default_value;
    get(dnn::leadChargedHadrCand_EoP) = default_value;
    if (leadChargedHadrCand->pt() > 0) {
      get(dnn::leadChargedHadrCand_HoP) = tau_funcs.getEcalEnergyLeadingChargedHadr(tau) / leadChargedHadrCand->pt();
      get(dnn::leadChargedHadrCand_EoP) = tau_funcs.getHcalEnergyLeadingChargedHadr(tau) / leadChargedHadrCand->pt();
    }

    MuonHitMatchV1 muon_hit_match;
    if (tau.leadPFChargedHadrCand().isNonnull() && tau.leadPFChargedHadrCand()->muonRef().isNonnull())
      muon_hit_match.addMatchedMuon(*tau.leadPFChargedHadrCand()->muonRef(), tau);

    auto matched_muons = muon_hit_match.findMatchedMuons(tau, muons, 0.3, 5);
    for (auto muon : matched_muons)
      muon_hit_match.addMatchedMuon(*muon, tau);
    muon_hit_match.fillTensor<dnn>(get, tau, default_value);

    LorentzVectorXYZ signalChargedHadrCands_sumIn, signalChargedHadrCands_sumOut;
    processSignalPFComponents(tau,
                              tau.signalChargedHadrCands(),
                              signalChargedHadrCands_sumIn,
                              signalChargedHadrCands_sumOut,
                              get(dnn::signalChargedHadrCands_sum_innerSigCone_pt),
                              get(dnn::signalChargedHadrCands_sum_innerSigCone_dEta),
                              get(dnn::signalChargedHadrCands_sum_innerSigCone_dPhi),
                              get(dnn::signalChargedHadrCands_sum_innerSigCone_mass),
                              get(dnn::signalChargedHadrCands_sum_outerSigCone_pt),
                              get(dnn::signalChargedHadrCands_sum_outerSigCone_dEta),
                              get(dnn::signalChargedHadrCands_sum_outerSigCone_dPhi),
                              get(dnn::signalChargedHadrCands_sum_outerSigCone_mass),
                              get(dnn::signalChargedHadrCands_nTotal_innerSigCone),
                              get(dnn::signalChargedHadrCands_nTotal_outerSigCone));

    LorentzVectorXYZ signalNeutrHadrCands_sumIn, signalNeutrHadrCands_sumOut;
    processSignalPFComponents(tau,
                              tau.signalNeutrHadrCands(),
                              signalNeutrHadrCands_sumIn,
                              signalNeutrHadrCands_sumOut,
                              get(dnn::signalNeutrHadrCands_sum_innerSigCone_pt),
                              get(dnn::signalNeutrHadrCands_sum_innerSigCone_dEta),
                              get(dnn::signalNeutrHadrCands_sum_innerSigCone_dPhi),
                              get(dnn::signalNeutrHadrCands_sum_innerSigCone_mass),
                              get(dnn::signalNeutrHadrCands_sum_outerSigCone_pt),
                              get(dnn::signalNeutrHadrCands_sum_outerSigCone_dEta),
                              get(dnn::signalNeutrHadrCands_sum_outerSigCone_dPhi),
                              get(dnn::signalNeutrHadrCands_sum_outerSigCone_mass),
                              get(dnn::signalNeutrHadrCands_nTotal_innerSigCone),
                              get(dnn::signalNeutrHadrCands_nTotal_outerSigCone));

    LorentzVectorXYZ signalGammaCands_sumIn, signalGammaCands_sumOut;
    processSignalPFComponents(tau,
                              tau.signalGammaCands(),
                              signalGammaCands_sumIn,
                              signalGammaCands_sumOut,
                              get(dnn::signalGammaCands_sum_innerSigCone_pt),
                              get(dnn::signalGammaCands_sum_innerSigCone_dEta),
                              get(dnn::signalGammaCands_sum_innerSigCone_dPhi),
                              get(dnn::signalGammaCands_sum_innerSigCone_mass),
                              get(dnn::signalGammaCands_sum_outerSigCone_pt),
                              get(dnn::signalGammaCands_sum_outerSigCone_dEta),
                              get(dnn::signalGammaCands_sum_outerSigCone_dPhi),
                              get(dnn::signalGammaCands_sum_outerSigCone_mass),
                              get(dnn::signalGammaCands_nTotal_innerSigCone),
                              get(dnn::signalGammaCands_nTotal_outerSigCone));

    LorentzVectorXYZ isolationChargedHadrCands_sum;
    processIsolationPFComponents(tau,
                                 tau.isolationChargedHadrCands(),
                                 isolationChargedHadrCands_sum,
                                 get(dnn::isolationChargedHadrCands_sum_pt),
                                 get(dnn::isolationChargedHadrCands_sum_dEta),
                                 get(dnn::isolationChargedHadrCands_sum_dPhi),
                                 get(dnn::isolationChargedHadrCands_sum_mass),
                                 get(dnn::isolationChargedHadrCands_nTotal));

    LorentzVectorXYZ isolationNeutrHadrCands_sum;
    processIsolationPFComponents(tau,
                                 tau.isolationNeutrHadrCands(),
                                 isolationNeutrHadrCands_sum,
                                 get(dnn::isolationNeutrHadrCands_sum_pt),
                                 get(dnn::isolationNeutrHadrCands_sum_dEta),
                                 get(dnn::isolationNeutrHadrCands_sum_dPhi),
                                 get(dnn::isolationNeutrHadrCands_sum_mass),
                                 get(dnn::isolationNeutrHadrCands_nTotal));

    LorentzVectorXYZ isolationGammaCands_sum;
    processIsolationPFComponents(tau,
                                 tau.isolationGammaCands(),
                                 isolationGammaCands_sum,
                                 get(dnn::isolationGammaCands_sum_pt),
                                 get(dnn::isolationGammaCands_sum_dEta),
                                 get(dnn::isolationGammaCands_sum_dPhi),
                                 get(dnn::isolationGammaCands_sum_mass),
                                 get(dnn::isolationGammaCands_nTotal));

    get(dnn::tau_visMass_innerSigCone) = (signalGammaCands_sumIn + signalChargedHadrCands_sumIn).mass();

    if (check_all_set) {
      for (int var_index = 0; var_index < dnn::NumberOfInputs; ++var_index) {
        if (get(var_index) == default_value_for_set_check)
          throw cms::Exception("DeepTauId: variable with index = ") << var_index << " is not set.";
      }
    }

    return inputs;
  }

  static void calculateElectronClusterVars(const pat::Electron* ele, float& elecEe, float& elecEgamma) {
    if (ele) {
      elecEe = elecEgamma = 0;
      auto superCluster = ele->superCluster();
      if (superCluster.isNonnull() && superCluster.isAvailable() && superCluster->clusters().isNonnull() &&
          superCluster->clusters().isAvailable()) {
        for (auto iter = superCluster->clustersBegin(); iter != superCluster->clustersEnd(); ++iter) {
          const double energy = (*iter)->energy();
          if (iter == superCluster->clustersBegin())
            elecEe += energy;
          else
            elecEgamma += energy;
        }
      }
    } else {
      elecEe = elecEgamma = default_value;
    }
  }

  template <typename CandidateCollection, typename TauCastType>
  static void processSignalPFComponents(const TauCastType& tau,
                                        const CandidateCollection& candidates,
                                        LorentzVectorXYZ& p4_inner,
                                        LorentzVectorXYZ& p4_outer,
                                        float& pt_inner,
                                        float& dEta_inner,
                                        float& dPhi_inner,
                                        float& m_inner,
                                        float& pt_outer,
                                        float& dEta_outer,
                                        float& dPhi_outer,
                                        float& m_outer,
                                        float& n_inner,
                                        float& n_outer) {
    p4_inner = LorentzVectorXYZ(0, 0, 0, 0);
    p4_outer = LorentzVectorXYZ(0, 0, 0, 0);
    n_inner = 0;
    n_outer = 0;

    const double innerSigCone_radius = getInnerSignalConeRadius(tau.pt());
    for (const auto& cand : candidates) {
      const double dR = reco::deltaR(cand->p4(), tau.leadChargedHadrCand()->p4());
      const bool isInside_innerSigCone = dR < innerSigCone_radius;
      if (isInside_innerSigCone) {
        p4_inner += cand->p4();
        ++n_inner;
      } else {
        p4_outer += cand->p4();
        ++n_outer;
      }
    }

    pt_inner = n_inner != 0 ? p4_inner.Pt() : default_value;
    dEta_inner = n_inner != 0 ? dEta(p4_inner, tau.p4()) : default_value;
    dPhi_inner = n_inner != 0 ? dPhi(p4_inner, tau.p4()) : default_value;
    m_inner = n_inner != 0 ? p4_inner.mass() : default_value;

    pt_outer = n_outer != 0 ? p4_outer.Pt() : default_value;
    dEta_outer = n_outer != 0 ? dEta(p4_outer, tau.p4()) : default_value;
    dPhi_outer = n_outer != 0 ? dPhi(p4_outer, tau.p4()) : default_value;
    m_outer = n_outer != 0 ? p4_outer.mass() : default_value;
  }

  template <typename CandidateCollection, typename TauCastType>
  static void processIsolationPFComponents(const TauCastType& tau,
                                           const CandidateCollection& candidates,
                                           LorentzVectorXYZ& p4,
                                           float& pt,
                                           float& d_eta,
                                           float& d_phi,
                                           float& m,
                                           float& n) {
    p4 = LorentzVectorXYZ(0, 0, 0, 0);
    n = 0;

    for (const auto& cand : candidates) {
      p4 += cand->p4();
      ++n;
    }

    pt = n != 0 ? p4.Pt() : default_value;
    d_eta = n != 0 ? dEta(p4, tau.p4()) : default_value;
    d_phi = n != 0 ? dPhi(p4, tau.p4()) : default_value;
    m = n != 0 ? p4.mass() : default_value;
  }

  template <typename TauCastType>
  static const pat::Electron* findMatchedElectron(const TauCastType& tau,
                                                  const std::vector<pat::Electron>* electrons,
                                                  double deltaR) {
    const double dR2 = deltaR * deltaR;
    const pat::Electron* matched_ele = nullptr;
    for (const auto& ele : *electrons) {
      if (reco::deltaR2(tau.p4(), ele.p4()) < dR2 && (!matched_ele || matched_ele->pt() < ele.pt())) {
        matched_ele = &ele;
      }
    }
    return matched_ele;
  }

private:
  edm::EDGetTokenT<std::vector<pat::Electron>> electrons_token_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muons_token_;
  edm::EDGetTokenT<double> rho_token_;
  edm::EDGetTokenT<reco::TauDiscriminatorContainer> basicTauDiscriminators_inputToken_;
  edm::EDGetTokenT<reco::TauDiscriminatorContainer> basicTauDiscriminatorsdR03_inputToken_;
  edm::EDGetTokenT<edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>>
      pfTauTransverseImpactParameters_token_;
  std::string input_layer_, output_layer_;
  const unsigned version_;
  const int debug_level;
  const bool disable_dxy_pca_;
  const bool disable_hcalFraction_workaround_;
  const bool disable_CellIndex_workaround_;
  std::unique_ptr<tensorflow::Tensor> tauBlockTensor_;
  std::array<std::unique_ptr<tensorflow::Tensor>, 2> eGammaTensor_, muonTensor_, hadronsTensor_, convTensor_,
      zeroOutputTensor_;
  const bool save_inputs_;
  std::ofstream* json_file_;
  bool is_first_block_;
  int file_counter_;

  //boolean to check if discriminator indices are already mapped
  bool discrIndicesMapped_ = false;
  std::map<BasicDiscriminator, size_t> basicDiscrIndexMap_;
  std::map<BasicDiscriminator, size_t> basicDiscrdR03IndexMap_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DeepTauId);
