#ifndef TauStudies_TruthAnalysis_plugins_TauTruthUtils_h
#define TauStudies_TruthAnalysis_plugins_TauTruthUtils_h

#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <cmath> // for std::sqrt
#include <algorithm> // for std::max
#include <vector>
#include <map>

#include "PhysicsTools/TruthInfo/interface/Branch.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/Particle.h"
#include "SimDataFormats/TruthInfo/interface/Vertex.h"

namespace TauTruthUtils {

/// Common hadrons seen in tau decays
/// Note: do not count intermediate resonance
///       rho0 (113), rho+(213), a0 (9000.11), eta (221), eta' (331), omega (223), ...
/// https://pdg.lbl.gov/2026/reviews/rpp2026-rev-monte-carlo-numbering.pdf#page=6
static inline const std::vector<int32_t> kTauDecayHadrons = {
    22, // photon
    111, 211, -211, // pions
    130, // neutral kaons
    321, -321, // charged kaons
    //221, // eta: stop here, so its photons are not followed into Geant4
    //22, // to capture rare hadrons+photon decay: 15 -> 16 22 111 -211 without secondaries
};

/// Inclsuive tau decay products
static inline const std::vector<int32_t> kTauDecayProducts = [] {
    std::vector<int32_t> v = kTauDecayHadrons;
    const std::vector<int32_t> kOther = { 11, -11, 12, -12, 13, -13, 14, -14 };
    v.insert(v.end(), kOther.begin(), kOther.end());
    return v;
}();

}  // namespace TauTruthUtils

#endif
