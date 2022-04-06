#include "RecoTauTag/RecoTau/interface/DeepTauHelper.h"

namespace deep_tau_2017 {

  float getTauID(const pat::Tau& tau, const std::string& tauID, float default_value, bool assert_input) {
    static tbb::concurrent_unordered_set<std::string> isFirstWarning;
    if (tau.isTauIDAvailable(tauID)) {
      return tau.tauID(tauID);
    } else {
      if (assert_input) {
        throw cms::Exception("DeepTauId")
            << "Exception in <getTauID>: No tauID '" << tauID << "' available in pat::Tau given as function argument.";
      }
      if (isFirstWarning.insert(tauID).second) {
        edm::LogWarning("DeepTauID") << "Warning in <getTauID>: No tauID '" << tauID
                                     << "' available in pat::Tau given as function argument."
                                     << " Using default_value = " << default_value << " instead." << std::endl;
      }
      return default_value;
    }
  }

  void MuonHitMatchV1::addMatchedMuon(const pat::Muon& muon, reco::BaseTau const& tau) {
    static constexpr int n_stations = 4;

    ++n_muons;
    const double dR2 = reco::deltaR2(tau.p4(), muon.p4());
    if (!best_matched_muon || dR2 < deltaR2_best_match) {
      best_matched_muon = &muon;
      deltaR2_best_match = dR2;
    }

    for (const auto& segment : muon.matches()) {
      if (segment.segmentMatches.empty())
        continue;
      if (n_matches.count(segment.detector()))
        ++n_matches.at(segment.detector()).at(segment.station() - 1);
    }

    if (muon.outerTrack().isNonnull()) {
      const auto& hit_pattern = muon.outerTrack()->hitPattern();
      for (int hit_index = 0; hit_index < hit_pattern.numberOfAllHits(reco::HitPattern::TRACK_HITS); ++hit_index) {
        auto hit_id = hit_pattern.getHitPattern(reco::HitPattern::TRACK_HITS, hit_index);
        if (hit_id == 0)
          break;
        if (hit_pattern.muonHitFilter(hit_id) && (hit_pattern.getHitType(hit_id) == TrackingRecHit::valid ||
                                                  hit_pattern.getHitType(hit_id == TrackingRecHit::bad))) {
          const int station = hit_pattern.getMuonStation(hit_id) - 1;
          if (station > 0 && station < n_stations) {
            std::vector<UInt_t>* muon_n_hits = nullptr;
            if (hit_pattern.muonDTHitFilter(hit_id))
              muon_n_hits = &n_hits.at(MuonSubdetId::DT);
            else if (hit_pattern.muonCSCHitFilter(hit_id))
              muon_n_hits = &n_hits.at(MuonSubdetId::CSC);
            else if (hit_pattern.muonRPCHitFilter(hit_id))
              muon_n_hits = &n_hits.at(MuonSubdetId::RPC);

            if (muon_n_hits)
              ++muon_n_hits->at(station);
          }
        }
      }
    }
  }

  unsigned MuonHitMatchV1::countMuonStationsWithMatches(size_t first_station, size_t last_station) const {
    static const std::map<int, std::vector<bool>> masks = {
        {MuonSubdetId::DT, {false, false, false, false}},
        {MuonSubdetId::CSC, {true, false, false, false}},
        {MuonSubdetId::RPC, {false, false, false, false}},
    };
    unsigned cnt = 0;
    for (unsigned n = first_station; n <= last_station; ++n) {
      for (const auto& match : n_matches) {
        if (!masks.at(match.first).at(n) && match.second.at(n) > 0)
          ++cnt;
      }
    }
    return cnt;
  }

  unsigned MuonHitMatchV1::countMuonStationsWithHits(size_t first_station, size_t last_station) const {
    static const std::map<int, std::vector<bool>> masks = {
        {MuonSubdetId::DT, {false, false, false, false}},
        {MuonSubdetId::CSC, {false, false, false, false}},
        {MuonSubdetId::RPC, {false, false, false, false}},
    };

    unsigned cnt = 0;
    for (unsigned n = first_station; n <= last_station; ++n) {
      for (const auto& hit : n_hits) {
        if (!masks.at(hit.first).at(n) && hit.second.at(n) > 0)
          ++cnt;
      }
    }
    return cnt;
  }

  template <>
  CellObjectType GetCellObjectType(const pat::Electron&) {
    return CellObjectType::Electron;
  }

  template <>
  CellObjectType GetCellObjectType(const pat::Muon&) {
    return CellObjectType::Muon;
  }

  template <>
  CellObjectType GetCellObjectType(reco::Candidate const& cand) {
    static const std::map<int, CellObjectType> obj_types = {{11, CellObjectType::PfCand_electron},
                                                            {13, CellObjectType::PfCand_muon},
                                                            {22, CellObjectType::PfCand_gamma},
                                                            {130, CellObjectType::PfCand_neutralHadron},
                                                            {211, CellObjectType::PfCand_chargedHadron}};

    auto iter = obj_types.find(std::abs(cand.pdgId()));
    if (iter == obj_types.end())
      return CellObjectType::Other;
    return iter->second;
  }

  bool CellGrid::tryGetCellIndex(double deltaEta, double deltaPhi, CellIndex& cellIndex) const {
    const auto getCellIndex = [this](double x, double maxX, double size, int& index) {
      const double absX = std::abs(x);
      if (absX > maxX)
        return false;
      double absIndex;
      if (disable_CellIndex_workaround_) {
        // CV: use consistent definition for CellIndex
        //     in DeepTauId.cc code and new DeepTau trainings
        absIndex = std::floor(absX / size + 0.5);
      } else {
        // CV: backwards compatibility with DeepTau training v2p1 used during Run 2
        absIndex = std::floor(std::abs(absX / size - 0.5));
      }
      index = static_cast<int>(std::copysign(absIndex, x));
      return true;
    };

    return getCellIndex(deltaEta, maxDeltaEta(), cellSizeEta, cellIndex.eta) &&
           getCellIndex(deltaPhi, maxDeltaPhi(), cellSizePhi, cellIndex.phi);
  }

}  // namespace deep_tau_2017

namespace deeptau_helper {
  const deep_tau::DeepTauBase::OutputCollection& GetOutputs() {
    static constexpr size_t e_index = 0, mu_index = 1, tau_index = 2, jet_index = 3;
    static const deep_tau::DeepTauBase::OutputCollection outputs_ = {
        {"VSe", deep_tau::DeepTauBase::Output({tau_index}, {e_index, tau_index})},
        {"VSmu", deep_tau::DeepTauBase::Output({tau_index}, {mu_index, tau_index})},
        {"VSjet", deep_tau::DeepTauBase::Output({tau_index}, {jet_index, tau_index})},
    };
    return outputs_;
  }

  bool isAbove(double value, double min) { return std::isnormal(value) && value > min; }

  bool calculateElectronClusterVarsV2(const pat::Electron& ele,
                                      float& cc_ele_energy,
                                      float& cc_gamma_energy,
                                      int& cc_n_gamma) {
    cc_ele_energy = cc_gamma_energy = 0;
    cc_n_gamma = 0;
    const auto& superCluster = ele.superCluster();
    if (superCluster.isNonnull() && superCluster.isAvailable() && superCluster->clusters().isNonnull() &&
        superCluster->clusters().isAvailable()) {
      for (auto iter = superCluster->clustersBegin(); iter != superCluster->clustersEnd(); ++iter) {
        const float energy = static_cast<float>((*iter)->energy());
        if (iter == superCluster->clustersBegin())
          cc_ele_energy += energy;
        else {
          cc_gamma_energy += energy;
          ++cc_n_gamma;
        }
      }
      return true;
    } else
      return false;
  }

  double getInnerSignalConeRadius(double pt) {
    static constexpr double min_pt = 30., min_radius = 0.05, cone_opening_coef = 3.;
    // This is equivalent of the original formula (std::max(std::min(0.1, 3.0/pt), 0.05)
    return std::max(cone_opening_coef / std::max(pt, min_pt), min_radius);
  }

  bool isInEcalCrack(double eta) {
    const double abs_eta = std::abs(eta);
    return abs_eta > 1.46 && abs_eta < 1.558;
  }
}  // namespace deeptau_helper
