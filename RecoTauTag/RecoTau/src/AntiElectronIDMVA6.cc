#include "RecoTauTag/RecoTau/interface/AntiElectronIDMVA6.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include <TMath.h>
#include <TFile.h>
#include <array>


namespace HGCal_helpers {

class coordinates {
 public:
  coordinates() : x(0), y(0), z(0), eta(100), phi(0) {}
  float x, y, z, eta, phi;
  inline math::XYZTLorentzVectorD toVector() { return math::XYZTLorentzVectorD(x, y, z, 0); }
};

class simpleTrackPropagator {
 public:
  simpleTrackPropagator(MagneticField const *f)
      : field_(f), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {
    ROOT::Math::SMatrixIdentity id;
    AlgebraicSymMatrix55 C(id);
    C *= 0.001;
    err_ = CurvilinearTrajectoryError(C);
  }
  void setPropagationTargetZ(const float &z);

  bool propagate(const double px, const double py, const double pz, const double x, const double y,
                 const double z, const float charge, coordinates &coords) const;

  bool propagate(const math::XYZTLorentzVectorD &momentum, const math::XYZTLorentzVectorD &position,
                 const float charge, coordinates &coords) const;

 private:
  simpleTrackPropagator() : field_(0), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {}
  const RKPropagatorInS &RKProp() const { return prod_.propagator; }
  Plane::PlanePointer targetPlaneForward_, targetPlaneBackward_;
  MagneticField const *field_;
  CurvilinearTrajectoryError err_;
  defaultRKPropagator::Product prod_;
  float absz_target_;
};

void simpleTrackPropagator::setPropagationTargetZ(const float &z) {
  targetPlaneForward_ = Plane::build(Plane::PositionType(0, 0, std::abs(z)), Plane::RotationType());
  targetPlaneBackward_ =
      Plane::build(Plane::PositionType(0, 0, -std::abs(z)), Plane::RotationType());
  absz_target_ = std::abs(z);
}

bool simpleTrackPropagator::propagate(const double px, const double py, const double pz,
                                      const double x, const double y, const double z,
                                      const float charge, coordinates &output) const {
  output = coordinates();

  typedef TrajectoryStateOnSurface TSOS;
  GlobalPoint startingPosition(x, y, z);
  GlobalVector startingMomentum(px, py, pz);
  Plane::PlanePointer startingPlane =
      Plane::build(Plane::PositionType(x, y, z), Plane::RotationType());
  TSOS startingStateP(
      GlobalTrajectoryParameters(startingPosition, startingMomentum, charge, field_), err_,
      *startingPlane);

  TSOS trackStateP;
  
  if (pz > 0) {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneForward_);
  } else {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneBackward_);
  }
  if (trackStateP.isValid()) {
    output.x = trackStateP.globalPosition().x();
    output.y = trackStateP.globalPosition().y();
    output.z = trackStateP.globalPosition().z();
    output.phi = trackStateP.globalPosition().phi();
    output.eta = trackStateP.globalPosition().eta();
    return true;
  }
  return false;
}

bool simpleTrackPropagator::propagate(const math::XYZTLorentzVectorD &momentum,
                                      const math::XYZTLorentzVectorD &position, const float charge,
                                      coordinates &output) const {
  return propagate(momentum.px(), momentum.py(), momentum.pz(), position.x(), position.y(),
                   position.z(), charge, output);
}

}
AntiElectronIDMVA6::AntiElectronIDMVA6(const edm::ParameterSet& cfg)
  : isInitialized_(false),
    mva_NoEleMatch_woGwoGSF_BL_(nullptr),
    mva_NoEleMatch_wGwoGSF_BL_(nullptr),
    mva_woGwGSF_BL_(nullptr),
    mva_wGwGSF_BL_(nullptr),
    mva_NoEleMatch_woGwoGSF_EC_(nullptr),
    mva_NoEleMatch_wGwoGSF_EC_(nullptr),
    mva_woGwGSF_EC_(nullptr),
    mva_wGwGSF_EC_(nullptr),
    mva_NoEleMatch_woGwoGSF_VFEC_(nullptr),
    mva_NoEleMatch_wGwoGSF_VFEC_(nullptr),
    mva_woGwGSF_VFEC_(nullptr),
    mva_wGwGSF_VFEC_(nullptr)	 
{
  loadMVAfromDB_ = cfg.exists("loadMVAfromDB") ? cfg.getParameter<bool>("loadMVAfromDB"): false;
  if ( !loadMVAfromDB_ ) {
    if(cfg.exists("inputFileName")){
      inputFileName_ = cfg.getParameter<edm::FileInPath>("inputFileName");
    }else throw cms::Exception("MVA input not defined") << "Requested to load tau MVA input from ROOT file but no file provided in cfg file";
    
  }

  mvaName_NoEleMatch_woGwoGSF_BL_ = cfg.getParameter<std::string>("mvaName_NoEleMatch_woGwoGSF_BL");
  mvaName_NoEleMatch_wGwoGSF_BL_ = cfg.getParameter<std::string>("mvaName_NoEleMatch_wGwoGSF_BL");
  mvaName_woGwGSF_BL_ = cfg.getParameter<std::string>("mvaName_woGwGSF_BL");
  mvaName_wGwGSF_BL_ = cfg.getParameter<std::string>("mvaName_wGwGSF_BL");
  mvaName_NoEleMatch_woGwoGSF_EC_ = cfg.getParameter<std::string>("mvaName_NoEleMatch_woGwoGSF_EC");
  mvaName_NoEleMatch_wGwoGSF_EC_ = cfg.getParameter<std::string>("mvaName_NoEleMatch_wGwoGSF_EC");
  mvaName_woGwGSF_EC_ = cfg.getParameter<std::string>("mvaName_woGwGSF_EC");
  mvaName_wGwGSF_EC_ = cfg.getParameter<std::string>("mvaName_wGwGSF_EC");
  mvaName_NoEleMatch_woGwoGSF_VFEC_ = cfg.getParameter<std::string>("mvaName_NoEleMatch_woGwoGSF_VFEC");
  mvaName_NoEleMatch_wGwoGSF_VFEC_ = cfg.getParameter<std::string>("mvaName_NoEleMatch_wGwoGSF_VFEC");
  mvaName_woGwGSF_VFEC_ = cfg.getParameter<std::string>("mvaName_woGwGSF_VFEC");
  mvaName_wGwGSF_VFEC_ = cfg.getParameter<std::string>("mvaName_wGwGSF_VFEC");
  
  
  //usePhiAtEcalEntranceExtrapolation_ = cfg.getParameter<bool>("usePhiAtEcalEntranceExtrapolation");

  Var_NoEleMatch_woGwoGSF_Barrel_ = new Float_t[9];
  Var_NoEleMatch_wGwoGSF_Barrel_ = new Float_t[17];
  Var_woGwGSF_Barrel_ = new Float_t[27];
  Var_wGwGSF_Barrel_ = new Float_t[36];
  Var_NoEleMatch_woGwoGSF_Endcap_ = new Float_t[6];
  Var_NoEleMatch_wGwoGSF_Endcap_ = new Float_t[14];
  Var_woGwGSF_Endcap_ = new Float_t[31];
  Var_wGwGSF_Endcap_ = new Float_t[38];
  Var_NoEleMatch_woGwoGSF_VFEndcap_ = new Float_t[6];
  Var_NoEleMatch_wGwoGSF_VFEndcap_ = new Float_t[14];
  Var_woGwGSF_VFEndcap_ = new Float_t[32];
  Var_wGwGSF_VFEndcap_ = new Float_t[40];

  bField_ = 0;
  verbosity_ = 0;
}

AntiElectronIDMVA6::~AntiElectronIDMVA6()
{
  delete[] Var_NoEleMatch_woGwoGSF_Barrel_;
  delete[] Var_NoEleMatch_wGwoGSF_Barrel_;
  delete[] Var_woGwGSF_Barrel_;
  delete[] Var_wGwGSF_Barrel_;
  delete[] Var_NoEleMatch_woGwoGSF_Endcap_;
  delete[] Var_NoEleMatch_wGwoGSF_Endcap_;
  delete[] Var_woGwGSF_Endcap_;
  delete[] Var_wGwGSF_Endcap_;
  delete[] Var_NoEleMatch_woGwoGSF_VFEndcap_;
  delete[] Var_NoEleMatch_wGwoGSF_VFEndcap_;
  delete[] Var_woGwGSF_VFEndcap_;
  delete[] Var_wGwGSF_VFEndcap_;

  if (!loadMVAfromDB_) {
    delete mva_NoEleMatch_woGwoGSF_BL_;
    delete mva_NoEleMatch_wGwoGSF_BL_;
    delete mva_woGwGSF_BL_;
    delete mva_wGwGSF_BL_;
    delete mva_NoEleMatch_woGwoGSF_EC_;
    delete mva_NoEleMatch_wGwoGSF_EC_;
    delete mva_woGwGSF_EC_;
    delete mva_wGwGSF_EC_;
    delete mva_NoEleMatch_woGwoGSF_VFEC_;
    delete mva_NoEleMatch_wGwoGSF_VFEC_;
    delete mva_woGwGSF_VFEC_;
    delete mva_wGwGSF_VFEC_;   
  }


  for ( std::vector<TFile*>::iterator it = inputFilesToDelete_.begin();
	it != inputFilesToDelete_.end(); ++it ) {
    delete (*it);
  }
}

namespace
{
  const GBRForest* loadMVAfromFile(TFile* inputFile, const std::string& mvaName)
  {
    const GBRForest* mva = (GBRForest*)inputFile->Get(mvaName.data());
    if ( !mva )
      throw cms::Exception("PFRecoTauDiscriminationAgainstElectronMVA6::loadMVA")
        << " Failed to load MVA = " << mvaName.data() << " from file " << " !!\n";

    return mva;
  }

  const GBRForest* loadMVAfromDB(const edm::EventSetup& es, const std::string& mvaName)
  {
    edm::ESHandle<GBRForest> mva;
    es.get<GBRWrapperRcd>().get(mvaName, mva);
    return mva.product();
  }
}

void AntiElectronIDMVA6::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( !isInitialized_ ) {
    if ( loadMVAfromDB_ ) {
      mva_NoEleMatch_woGwoGSF_BL_ = loadMVAfromDB(es, mvaName_NoEleMatch_woGwoGSF_BL_);
      mva_NoEleMatch_wGwoGSF_BL_ = loadMVAfromDB(es, mvaName_NoEleMatch_wGwoGSF_BL_);
      mva_woGwGSF_BL_ = loadMVAfromDB(es, mvaName_woGwGSF_BL_);
      mva_wGwGSF_BL_ = loadMVAfromDB(es, mvaName_wGwGSF_BL_);
      mva_NoEleMatch_woGwoGSF_EC_ = loadMVAfromDB(es, mvaName_NoEleMatch_woGwoGSF_EC_);
      mva_NoEleMatch_wGwoGSF_EC_ = loadMVAfromDB(es, mvaName_NoEleMatch_wGwoGSF_EC_);
      mva_woGwGSF_EC_ = loadMVAfromDB(es, mvaName_woGwGSF_EC_);
      mva_wGwGSF_EC_ = loadMVAfromDB(es, mvaName_wGwGSF_EC_);
      mva_NoEleMatch_woGwoGSF_VFEC_ = loadMVAfromDB(es, mvaName_NoEleMatch_woGwoGSF_VFEC_);
      mva_NoEleMatch_wGwoGSF_VFEC_ = loadMVAfromDB(es, mvaName_NoEleMatch_wGwoGSF_VFEC_);
      mva_woGwGSF_VFEC_ = loadMVAfromDB(es, mvaName_woGwGSF_VFEC_);
      mva_wGwGSF_VFEC_ = loadMVAfromDB(es, mvaName_wGwGSF_VFEC_);
    } else {
          if ( inputFileName_.location() == edm::FileInPath::Unknown ) throw cms::Exception("PFRecoTauDiscriminationAgainstElectronMVA6::loadMVA")
          << " Failed to find File = " << inputFileName_ << " !!\n";
          TFile* inputFile = new TFile(inputFileName_.fullPath().data());

      mva_NoEleMatch_woGwoGSF_BL_ = loadMVAfromFile(inputFile, mvaName_NoEleMatch_woGwoGSF_BL_);
      mva_NoEleMatch_wGwoGSF_BL_ = loadMVAfromFile(inputFile, mvaName_NoEleMatch_wGwoGSF_BL_);
      mva_woGwGSF_BL_ = loadMVAfromFile(inputFile, mvaName_woGwGSF_BL_);
      mva_wGwGSF_BL_ = loadMVAfromFile(inputFile, mvaName_wGwGSF_BL_);
      mva_NoEleMatch_woGwoGSF_EC_ = loadMVAfromFile(inputFile, mvaName_NoEleMatch_woGwoGSF_EC_);
      mva_NoEleMatch_wGwoGSF_EC_ = loadMVAfromFile(inputFile, mvaName_NoEleMatch_wGwoGSF_EC_);
      mva_woGwGSF_EC_ = loadMVAfromFile(inputFile, mvaName_woGwGSF_EC_);
      mva_wGwGSF_EC_ = loadMVAfromFile(inputFile, mvaName_wGwGSF_EC_);
      mva_NoEleMatch_woGwoGSF_VFEC_ = loadMVAfromFile(inputFile, mvaName_NoEleMatch_woGwoGSF_VFEC_);
      mva_NoEleMatch_wGwoGSF_VFEC_ = loadMVAfromFile(inputFile, mvaName_NoEleMatch_wGwoGSF_VFEC_);
      mva_woGwGSF_VFEC_ = loadMVAfromFile(inputFile, mvaName_woGwGSF_VFEC_);
      mva_wGwGSF_VFEC_ = loadMVAfromFile(inputFile, mvaName_wGwGSF_VFEC_);
      inputFilesToDelete_.push_back(inputFile);  
    }
    isInitialized_ = true;
  }

  edm::ESHandle<MagneticField> pSetup;
  es.get<IdealMagneticFieldRecord>().get(pSetup);
  bField_ = pSetup->inTesla(GlobalPoint(0,0,0)).z();
  
  //HGCAL
  aField_ = &(*pSetup);
  recHitTools_.getEventSetup(es);
  
  for (unsigned ilayer = 1; ilayer <= 52; ++ilayer) {
    const GlobalPoint pos = recHitTools_.getPositionLayer(ilayer);
    layerPositions_.push_back(pos.z());
  }

}

double AntiElectronIDMVA6::MVAValue(Float_t TauPt,                                         
                                    Float_t TauEtaAtEcalEntrance,
                                    Float_t TauPhi,
                                    Float_t TauLeadChargedPFCandPt,
                                    Float_t TauLeadChargedPFCandEtaAtEcalEntrance,
                                    Float_t TauEmFraction,
                                    Float_t TauLeadPFChargedHadrHoP,
                                    Float_t TauLeadPFChargedHadrEoP,
                                    Float_t TauVisMassIn,
                                    Float_t TaudCrackEta,
                                    Float_t TauHasGsf,
                                    Int_t TauSignalPFGammaCandsIn,
                                    Int_t TauSignalPFGammaCandsOut,
                                    const std::vector<Float_t>& GammasdEtaInSigCone,
                                    const std::vector<Float_t>& GammasdPhiInSigCone,
                                    const std::vector<Float_t>& GammasPtInSigCone,
                                    const std::vector<Float_t>& GammasdEtaOutSigCone,
                                    const std::vector<Float_t>& GammasdPhiOutSigCone,
                                    const std::vector<Float_t>& GammasPtOutSigCone,
                                    Float_t ElecEta,
                                    Float_t ElecPhi,
                                    Float_t ElecChi2NormGSF,
                                    Float_t ElecChi2NormKF,
                                    Float_t ElecGSFNumHits,
                                    Float_t ElecKFNumHits,
                                    Float_t ElecGSFTrackResol,
                                    Float_t ElecGSFTracklnPt,
                                    Float_t ElecPin,
                                    Float_t ElecPout,
                                    Float_t ElecEecal,
                                    Float_t ElecDeltaEta,
                                    Float_t ElecDeltaPhi,
                                    Float_t ElecMvaInSigmaEtaEta,
                                    Float_t ElecMvaInHadEnergy,
                                    Float_t ElecMvaInDeltaEta,
				    Float_t HGCAL_sigmaUU,
     	 			    Float_t HGCAL_sigmaVV,
     	 			    Float_t HGCAL_sigmaEE,
     	 			    Float_t HGCAL_sigmaPP,
     	 			    Float_t HGCAL_nLayers,
				    Float_t HGCAL_firstLayer,
     	 			    Float_t HGCAL_lastLayer,
     	 			    Float_t HGCAL_layerEfrac10,
     	 			    Float_t HGCAL_layerEfrac90,
     	 			    Float_t HGCAL_ecEnergyEE,
     	 			    Float_t HGCAL_ecEnergyFH,
     	 			    Float_t HGCAL_measuredDepth,
     	 			    Float_t HGCAL_expectedDepth,
     	 			    Float_t HGCAL_expectedSigma,
     	 			    Float_t HGCAL_depthCompatibility,
     	 			    Float_t ElecESeedClusterOverPout,	  
     	 			    Float_t ElecSuperClusterEtaWidth,
     	 			    Float_t ElecSuperClusterPhiWidth,      	 		
     	 			    Float_t ElecSigmaIEtaIEta5x5,
     	 			    Float_t ElecSigmaIPhiIPhi5x5,
     	 			    Float_t ElecShowerCircularity,
     	 			    Float_t ElecR9
				    ) 
{ 

  double sumPt  = 0.;
  double dEta2  = 0.;
  double dPhi2  = 0.;
  double sumPt2 = 0.;
  for ( unsigned int i = 0 ; i < GammasPtInSigCone.size() ; ++i ) {
    double pt_i  = GammasPtInSigCone[i];
    double phi_i = GammasdPhiInSigCone[i];
    if ( GammasdPhiInSigCone[i] > M_PI ) phi_i = GammasdPhiInSigCone[i] - 2*M_PI;
    else if ( GammasdPhiInSigCone[i] < -M_PI ) phi_i = GammasdPhiInSigCone[i] + 2*M_PI;
    double eta_i = GammasdEtaInSigCone[i];
    sumPt  +=  pt_i;
    sumPt2 += (pt_i*pt_i);
    dEta2  += (pt_i*eta_i*eta_i);
    dPhi2  += (pt_i*phi_i*phi_i);
  }
  Float_t TauGammaEnFracIn = -99.;
  if ( TauPt > 0. ) {
    TauGammaEnFracIn = sumPt/TauPt;
  }
  if ( sumPt > 0. ) {
    dEta2 /= sumPt;
    dPhi2 /= sumPt;
  }
  Float_t TauGammaEtaMomIn = std::sqrt(dEta2)*std::sqrt(TauGammaEnFracIn)*TauPt;
  Float_t TauGammaPhiMomIn = std::sqrt(dPhi2)*std::sqrt(TauGammaEnFracIn)*TauPt;

  sumPt  = 0.;
  dEta2  = 0.;
  dPhi2  = 0.;
  sumPt2 = 0.;
  for ( unsigned int i = 0 ; i < GammasPtOutSigCone.size() ; ++i ) {
    double pt_i  = GammasPtOutSigCone[i];
    double phi_i = GammasdPhiOutSigCone[i];
    if ( GammasdPhiOutSigCone[i] > M_PI ) phi_i = GammasdPhiOutSigCone[i] - 2*M_PI;
    else if ( GammasdPhiOutSigCone[i] < -M_PI ) phi_i = GammasdPhiOutSigCone[i] + 2*M_PI;
    double eta_i = GammasdEtaOutSigCone[i];
    sumPt  +=  pt_i;
    sumPt2 += (pt_i*pt_i);
    dEta2  += (pt_i*eta_i*eta_i);
    dPhi2  += (pt_i*phi_i*phi_i);
  }
  Float_t TauGammaEnFracOut = sumPt/TauPt;
  if ( sumPt > 0. ) {
    dEta2 /= sumPt;
    dPhi2 /= sumPt;
  }
  Float_t TauGammaEtaMomOut = std::sqrt(dEta2)*std::sqrt(TauGammaEnFracOut)*TauPt;
  Float_t TauGammaPhiMomOut = std::sqrt(dPhi2)*std::sqrt(TauGammaEnFracOut)*TauPt;
  
  return MVAValue(TauPt,
                  TauEtaAtEcalEntrance,
                  TauPhi,
                  TauLeadChargedPFCandPt,
                  TauLeadChargedPFCandEtaAtEcalEntrance,
                  TauEmFraction,
                  TauLeadPFChargedHadrHoP,
                  TauLeadPFChargedHadrEoP,
                  TauVisMassIn,
                  TaudCrackEta,
                  TauHasGsf,
                  TauSignalPFGammaCandsIn,
                  TauSignalPFGammaCandsOut,
                  TauGammaEtaMomIn,
                  TauGammaEtaMomOut,
                  TauGammaPhiMomIn,
                  TauGammaPhiMomOut,
                  TauGammaEnFracIn,
                  TauGammaEnFracOut,
                  ElecEta,
                  ElecPhi,
                  ElecChi2NormGSF,
                  ElecChi2NormKF,
                  ElecGSFNumHits,
                  ElecKFNumHits,
                  ElecGSFTrackResol,
                  ElecGSFTracklnPt,
                  ElecPin,
                  ElecPout,
                  ElecEecal,
                  ElecDeltaEta,
                  ElecDeltaPhi,
                  ElecMvaInSigmaEtaEta,
                  ElecMvaInHadEnergy,
                  ElecMvaInDeltaEta,
		  HGCAL_sigmaUU,
     	 	  HGCAL_sigmaVV,
     	 	  HGCAL_sigmaEE,
     	 	  HGCAL_sigmaPP,
     	 	  HGCAL_nLayers,
     	 	  HGCAL_firstLayer,
     	 	  HGCAL_lastLayer,
     	 	  HGCAL_layerEfrac10,
     	 	  HGCAL_layerEfrac90,
     	 	  HGCAL_ecEnergyEE,
     	 	  HGCAL_ecEnergyFH,
     	 	  HGCAL_measuredDepth,
     	 	  HGCAL_expectedDepth,
     	 	  HGCAL_expectedSigma,
     	 	  HGCAL_depthCompatibility,
     	 	  ElecESeedClusterOverPout,	
     	 	  ElecSuperClusterEtaWidth,				    
     	 	  ElecSuperClusterPhiWidth, 
     	 	  ElecSigmaIEtaIEta5x5,
     	 	  ElecSigmaIPhiIPhi5x5,
     	 	  ElecShowerCircularity,
     	 	  ElecR9
		  );
}

double AntiElectronIDMVA6::MVAValue(Float_t TauPt,                                         
                                    Float_t TauEtaAtEcalEntrance,
                                    Float_t TauPhi,
                                    Float_t TauLeadChargedPFCandPt,
                                    Float_t TauLeadChargedPFCandEtaAtEcalEntrance,
                                    Float_t TauEmFraction,
                                    Float_t TauLeadPFChargedHadrHoP,
                                    Float_t TauLeadPFChargedHadrEoP,
                                    Float_t TauVisMassIn,
                                    Float_t TaudCrackEta,
                                    Float_t TauHasGsf,
                                    Int_t TauSignalPFGammaCandsIn,
                                    Int_t TauSignalPFGammaCandsOut,
                                    Float_t TauGammaEtaMomIn,
                                    Float_t TauGammaEtaMomOut,
                                    Float_t TauGammaPhiMomIn,
                                    Float_t TauGammaPhiMomOut,
                                    Float_t TauGammaEnFracIn,
                                    Float_t TauGammaEnFracOut,
                                    Float_t ElecEta,
                                    Float_t ElecPhi,
                                    Float_t ElecChi2NormGSF,
                                    Float_t ElecChi2NormKF,
                                    Float_t ElecGSFNumHits,
                                    Float_t ElecKFNumHits,
                                    Float_t ElecGSFTrackResol,
                                    Float_t ElecGSFTracklnPt,
                                    Float_t ElecPin,
                                    Float_t ElecPout,
                                    Float_t ElecEecal,
                                    Float_t ElecDeltaEta,
                                    Float_t ElecDeltaPhi,
                                    Float_t ElecMvaInSigmaEtaEta,
                                    Float_t ElecMvaInHadEnergy,
                                    Float_t ElecMvaInDeltaEta,    
     	 			    Float_t HGCAL_sigmaUU,
     	 			    Float_t HGCAL_sigmaVV,
     	 			    Float_t HGCAL_sigmaEE,
     	 			    Float_t HGCAL_sigmaPP,
     	 			    Float_t HGCAL_nLayers,
				    Float_t HGCAL_firstLayer,
     	 			    Float_t HGCAL_lastLayer,
     	 			    Float_t HGCAL_layerEfrac10,
     	 			    Float_t HGCAL_layerEfrac90,
     	 			    Float_t HGCAL_ecEnergyEE,
     	 			    Float_t HGCAL_ecEnergyFH,
     	 			    Float_t HGCAL_measuredDepth,
     	 			    Float_t HGCAL_expectedDepth,
     	 			    Float_t HGCAL_expectedSigma,
     	 			    Float_t HGCAL_depthCompatibility,
     	 			    Float_t ElecESeedClusterOverPout,	  
     	 			    Float_t ElecSuperClusterEtaWidth,
     	 			    Float_t ElecSuperClusterPhiWidth, 
     	 			    Float_t ElecSigmaIEtaIEta5x5,
     	 			    Float_t ElecSigmaIPhiIPhi5x5,
     	 			    Float_t ElecShowerCircularity,
     	 			    Float_t ElecR9
				    )     	
{
if (!isInitialized_) {
    throw cms::Exception("ClassNotInitialized") << " AntiElectronMVA not properly initialized !!\n";
  }

  double mvaValue = -99.;

  const float ECALBarrelEndcapEtaBorder = 1.479;
  const float ECALBarrelVFEndcapEtaBorder = 2.4;
  const float HGCALBorder = 3.0;

				     
  float ElecDeltaPinPoutOverPin = std::abs(ElecPin - ElecPout) / ElecPin;
  float ElecEecalOverPout = (ElecEecal / ElecPout);
  float ElecNumHitsDiffOverSum = (ElecGSFNumHits - ElecKFNumHits) / (ElecGSFNumHits + ElecKFNumHits);
  
                                     				     
   if (TauSignalPFGammaCandsIn == 0 && TauHasGsf < 0.5) {
    if (std::abs(TauEtaAtEcalEntrance) < ECALBarrelEndcapEtaBorder) {
      Var_NoEleMatch_woGwoGSF_Barrel_[0] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_NoEleMatch_woGwoGSF_Barrel_[1] = std::log(std::max(float(1.), TauPt));
      Var_NoEleMatch_woGwoGSF_Barrel_[2] = TauEmFraction;
      Var_NoEleMatch_woGwoGSF_Barrel_[3] = TauLeadPFChargedHadrHoP;
      Var_NoEleMatch_woGwoGSF_Barrel_[4] = TauLeadPFChargedHadrEoP;
      Var_NoEleMatch_woGwoGSF_Barrel_[5] = TauVisMassIn;
      Var_NoEleMatch_woGwoGSF_Barrel_[6] = TaudCrackEta;
      Var_NoEleMatch_woGwoGSF_Barrel_[7] = TauEtaAtEcalEntrance;
      Var_NoEleMatch_woGwoGSF_Barrel_[8] = TauLeadChargedPFCandEtaAtEcalEntrance;
      mvaValue = mva_NoEleMatch_woGwoGSF_BL_->GetClassifier(Var_NoEleMatch_woGwoGSF_Barrel_);      
        
    } else if (std::abs(TauEtaAtEcalEntrance) > ECALBarrelEndcapEtaBorder
            && std::abs(TauEtaAtEcalEntrance) < ECALBarrelVFEndcapEtaBorder) {
      Var_NoEleMatch_woGwoGSF_Endcap_[0] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_NoEleMatch_woGwoGSF_Endcap_[1] = std::log(std::max(float(1.), TauPt));
      Var_NoEleMatch_woGwoGSF_Endcap_[2] = TauVisMassIn;
      Var_NoEleMatch_woGwoGSF_Endcap_[3] = TaudCrackEta;
      Var_NoEleMatch_woGwoGSF_Endcap_[4] = TauEtaAtEcalEntrance;
      Var_NoEleMatch_woGwoGSF_Endcap_[5] = TauLeadChargedPFCandEtaAtEcalEntrance;
      mvaValue = mva_NoEleMatch_woGwoGSF_EC_->GetClassifier(Var_NoEleMatch_woGwoGSF_Endcap_);
    } else if (std::abs(TauEtaAtEcalEntrance) > ECALBarrelVFEndcapEtaBorder
            && std::abs(TauEtaAtEcalEntrance) < HGCALBorder){
      Var_NoEleMatch_woGwoGSF_VFEndcap_[0] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_NoEleMatch_woGwoGSF_VFEndcap_[1] = std::log(std::max(float(1.), TauPt));
      Var_NoEleMatch_woGwoGSF_VFEndcap_[2] = TauVisMassIn;
      Var_NoEleMatch_woGwoGSF_VFEndcap_[3] = TaudCrackEta;
      Var_NoEleMatch_woGwoGSF_VFEndcap_[4] = TauEtaAtEcalEntrance;
      Var_NoEleMatch_woGwoGSF_VFEndcap_[5] = TauLeadChargedPFCandEtaAtEcalEntrance;
      mvaValue = mva_NoEleMatch_woGwoGSF_VFEC_->GetClassifier(Var_NoEleMatch_woGwoGSF_VFEndcap_);          
    } 
  } 
  
  else if (TauSignalPFGammaCandsIn > 0 && TauHasGsf < 0.5) {
    if (std::abs(TauEtaAtEcalEntrance) < ECALBarrelEndcapEtaBorder) {
      Var_NoEleMatch_wGwoGSF_Barrel_[0] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_NoEleMatch_wGwoGSF_Barrel_[1] = std::log(std::max(float(1.), TauPt));
      Var_NoEleMatch_wGwoGSF_Barrel_[2] = TauEmFraction;
      Var_NoEleMatch_wGwoGSF_Barrel_[3] = TauSignalPFGammaCandsIn;
      Var_NoEleMatch_wGwoGSF_Barrel_[4] = TauSignalPFGammaCandsOut;
      Var_NoEleMatch_wGwoGSF_Barrel_[5] = TauLeadPFChargedHadrHoP;
      Var_NoEleMatch_wGwoGSF_Barrel_[6] = TauLeadPFChargedHadrEoP;
      Var_NoEleMatch_wGwoGSF_Barrel_[7] = TauVisMassIn;
      Var_NoEleMatch_wGwoGSF_Barrel_[8] = TauGammaEtaMomIn;
      Var_NoEleMatch_wGwoGSF_Barrel_[9] = TauGammaEtaMomOut;
      Var_NoEleMatch_wGwoGSF_Barrel_[10] = TauGammaPhiMomIn;
      Var_NoEleMatch_wGwoGSF_Barrel_[11] = TauGammaPhiMomOut;
      Var_NoEleMatch_wGwoGSF_Barrel_[12] = TauGammaEnFracIn;
      Var_NoEleMatch_wGwoGSF_Barrel_[13] = TauGammaEnFracOut;
      Var_NoEleMatch_wGwoGSF_Barrel_[14] = TaudCrackEta;
      Var_NoEleMatch_wGwoGSF_Barrel_[15] = TauEtaAtEcalEntrance;
      Var_NoEleMatch_wGwoGSF_Barrel_[16] = TauLeadChargedPFCandEtaAtEcalEntrance;     
      mvaValue = mva_NoEleMatch_wGwoGSF_BL_->GetClassifier(Var_NoEleMatch_wGwoGSF_Barrel_);
    } else if (std::abs(TauEtaAtEcalEntrance) > ECALBarrelEndcapEtaBorder &&                       
               std::abs(TauEtaAtEcalEntrance) < ECALBarrelVFEndcapEtaBorder) {
      Var_NoEleMatch_wGwoGSF_Endcap_[0] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_NoEleMatch_wGwoGSF_Endcap_[1] = std::log(std::max(float(1.), TauPt));
      Var_NoEleMatch_wGwoGSF_Endcap_[2] = TauSignalPFGammaCandsIn;
      Var_NoEleMatch_wGwoGSF_Endcap_[3] = TauSignalPFGammaCandsOut;
      Var_NoEleMatch_wGwoGSF_Endcap_[4] = TauVisMassIn;
      Var_NoEleMatch_wGwoGSF_Endcap_[5] = TauGammaEtaMomIn;
      Var_NoEleMatch_wGwoGSF_Endcap_[6] = TauGammaEtaMomOut;
      Var_NoEleMatch_wGwoGSF_Endcap_[7] = TauGammaPhiMomIn;
      Var_NoEleMatch_wGwoGSF_Endcap_[8] = TauGammaPhiMomOut;
      Var_NoEleMatch_wGwoGSF_Endcap_[9] = TauGammaEnFracIn;
      Var_NoEleMatch_wGwoGSF_Endcap_[10] = TauGammaEnFracOut;
      Var_NoEleMatch_wGwoGSF_Endcap_[11] = TaudCrackEta;
      Var_NoEleMatch_wGwoGSF_Endcap_[12] = TauEtaAtEcalEntrance;
      Var_NoEleMatch_wGwoGSF_Endcap_[13] = TauLeadChargedPFCandEtaAtEcalEntrance;
      mvaValue = mva_NoEleMatch_wGwoGSF_EC_->GetClassifier(Var_NoEleMatch_wGwoGSF_Endcap_);
    } else if (std::abs(TauEtaAtEcalEntrance) > ECALBarrelVFEndcapEtaBorder 
            && std::abs(TauEtaAtEcalEntrance) < HGCALBorder){
      Var_NoEleMatch_wGwoGSF_VFEndcap_[0] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_NoEleMatch_wGwoGSF_VFEndcap_[1] = std::log(std::max(float(1.), TauPt));
      Var_NoEleMatch_wGwoGSF_VFEndcap_[2] = TauSignalPFGammaCandsIn;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[3] = TauSignalPFGammaCandsOut;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[4] = TauVisMassIn;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[5] = TauGammaEtaMomIn;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[6] = TauGammaEtaMomOut;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[7] = TauGammaPhiMomIn;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[8] = TauGammaPhiMomOut;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[9] = TauGammaEnFracIn;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[10] = TauGammaEnFracOut;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[11] = TaudCrackEta;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[12] = TauEtaAtEcalEntrance;
      Var_NoEleMatch_wGwoGSF_VFEndcap_[13] = TauLeadChargedPFCandEtaAtEcalEntrance;   
      mvaValue = mva_NoEleMatch_wGwoGSF_VFEC_->GetClassifier(Var_NoEleMatch_wGwoGSF_VFEndcap_);
    }  
  } 
  
  else if (TauSignalPFGammaCandsIn == 0 && TauHasGsf > 0.5) {
  
    if (std::abs(TauEtaAtEcalEntrance) < ECALBarrelEndcapEtaBorder) {
      Var_woGwGSF_Barrel_[0] = std::log(std::max(float(0.1), ElecChi2NormGSF));
      Var_woGwGSF_Barrel_[1] = ElecGSFNumHits;
      Var_woGwGSF_Barrel_[2] = std::log(std::max(float(0.1), ElecGSFTrackResol));
      Var_woGwGSF_Barrel_[3] = ElecGSFTracklnPt;
      Var_woGwGSF_Barrel_[4] = ElecNumHitsDiffOverSum;
      Var_woGwGSF_Barrel_[5] = std::log(std::max(float(0.1), ElecChi2NormKF));
      Var_woGwGSF_Barrel_[6] = std::min(ElecDeltaPinPoutOverPin, float(1.));
      Var_woGwGSF_Barrel_[7] = std::min(ElecEecalOverPout, float(20.));
      Var_woGwGSF_Barrel_[8] = std::min(ElecMvaInSigmaEtaEta, float(0.01));
      Var_woGwGSF_Barrel_[9] = std::min(ElecMvaInHadEnergy, float(20.));
      Var_woGwGSF_Barrel_[10] = std::min(ElecMvaInDeltaEta, float(0.1));
      Var_woGwGSF_Barrel_[11] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_woGwGSF_Barrel_[12] = std::log(std::max(float(1.), TauPt));     
      Var_woGwGSF_Barrel_[13] = TauEmFraction;
      Var_woGwGSF_Barrel_[14] = TauLeadPFChargedHadrHoP;
      Var_woGwGSF_Barrel_[15] = TauLeadPFChargedHadrEoP;
      Var_woGwGSF_Barrel_[16] = TauVisMassIn;
      Var_woGwGSF_Barrel_[17] = TaudCrackEta;
      Var_woGwGSF_Barrel_[18] = TauEtaAtEcalEntrance;
      Var_woGwGSF_Barrel_[19] = TauLeadChargedPFCandEtaAtEcalEntrance;   
      Var_woGwGSF_Barrel_[20] = ElecDeltaEta;
      Var_woGwGSF_Barrel_[21] = ElecDeltaPhi;
      Var_woGwGSF_Barrel_[22] = ElecSigmaIEtaIEta5x5;
      Var_woGwGSF_Barrel_[23] = ElecShowerCircularity;
      Var_woGwGSF_Barrel_[24] = ElecR9;
      Var_woGwGSF_Barrel_[25] = ElecSuperClusterEtaWidth;
      Var_woGwGSF_Barrel_[26] = ElecSuperClusterPhiWidth;
      mvaValue = mva_woGwGSF_BL_->GetClassifier(Var_woGwGSF_Barrel_);
    } else if (std::abs(TauEtaAtEcalEntrance) > ECALBarrelEndcapEtaBorder 
            && std::abs(TauEtaAtEcalEntrance) < ECALBarrelVFEndcapEtaBorder ) {
      Var_woGwGSF_Endcap_[0] = std::log(std::max(float(0.1), ElecChi2NormGSF));
      Var_woGwGSF_Endcap_[1] = ElecGSFNumHits;
      Var_woGwGSF_Endcap_[2] = std::log(std::max(float(0.1), ElecGSFTrackResol));
      Var_woGwGSF_Endcap_[3] = ElecGSFTracklnPt;
      Var_woGwGSF_Endcap_[4] = ElecNumHitsDiffOverSum;
      Var_woGwGSF_Endcap_[5] = std::log(std::max(float(0.1), ElecChi2NormKF));
      Var_woGwGSF_Endcap_[6] = ElecEecal;
      Var_woGwGSF_Endcap_[7] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_woGwGSF_Endcap_[8] = std::log(std::max(float(1.), TauPt));
      Var_woGwGSF_Endcap_[9] = TauVisMassIn;
      Var_woGwGSF_Endcap_[10] = TaudCrackEta;
      Var_woGwGSF_Endcap_[11] = TauEtaAtEcalEntrance;
      Var_woGwGSF_Endcap_[12] = TauLeadChargedPFCandEtaAtEcalEntrance;
      Var_woGwGSF_Endcap_[13] = HGCAL_sigmaUU;
      Var_woGwGSF_Endcap_[14] = HGCAL_sigmaVV;
      Var_woGwGSF_Endcap_[15] = HGCAL_sigmaEE;
      Var_woGwGSF_Endcap_[16] = HGCAL_sigmaPP;
      Var_woGwGSF_Endcap_[17] = HGCAL_nLayers;
      Var_woGwGSF_Endcap_[18] = HGCAL_lastLayer;
      Var_woGwGSF_Endcap_[19] = HGCAL_layerEfrac10;
      Var_woGwGSF_Endcap_[20] = HGCAL_layerEfrac90;
      Var_woGwGSF_Endcap_[21] = HGCAL_ecEnergyEE;
      Var_woGwGSF_Endcap_[22] = HGCAL_ecEnergyFH;
      Var_woGwGSF_Endcap_[23] = HGCAL_measuredDepth;
      Var_woGwGSF_Endcap_[24] = HGCAL_expectedDepth;
      Var_woGwGSF_Endcap_[25] = HGCAL_depthCompatibility;
      Var_woGwGSF_Endcap_[26] = ElecDeltaEta;
      Var_woGwGSF_Endcap_[27] = ElecDeltaPhi;
      Var_woGwGSF_Endcap_[28] = ElecESeedClusterOverPout;
      Var_woGwGSF_Endcap_[29] = ElecSuperClusterEtaWidth;
      Var_woGwGSF_Endcap_[30] = ElecSuperClusterPhiWidth;
      mvaValue = mva_woGwGSF_EC_->GetClassifier(Var_woGwGSF_Endcap_);
    } else  if (std::abs(TauEtaAtEcalEntrance) > ECALBarrelVFEndcapEtaBorder
            && std::abs(TauEtaAtEcalEntrance) < HGCALBorder){
      Var_woGwGSF_VFEndcap_[0] = std::log(std::max(float(0.1), ElecChi2NormGSF));
      Var_woGwGSF_VFEndcap_[1] = ElecGSFNumHits;
      Var_woGwGSF_VFEndcap_[2] = std::log(std::max(float(0.1), ElecGSFTrackResol));
      Var_woGwGSF_VFEndcap_[3] = ElecGSFTracklnPt;
      Var_woGwGSF_VFEndcap_[4] = ElecNumHitsDiffOverSum;
      Var_woGwGSF_VFEndcap_[5] = std::log(std::max(float(0.1), ElecChi2NormKF));
      Var_woGwGSF_VFEndcap_[6] = ElecEecal;
      Var_woGwGSF_VFEndcap_[7] = TMath::Min(float(2.), TauLeadChargedPFCandPt/TMath::Max(float(1.), TauPt));
      Var_woGwGSF_VFEndcap_[8] = std::log(TMath::Max(float(1.), TauPt));
      Var_woGwGSF_VFEndcap_[9] = TauVisMassIn;						   
      Var_woGwGSF_VFEndcap_[10] = TaudCrackEta;
      Var_woGwGSF_VFEndcap_[11] = TauEtaAtEcalEntrance;
      Var_woGwGSF_VFEndcap_[12] = TauLeadChargedPFCandEtaAtEcalEntrance;
      Var_woGwGSF_VFEndcap_[13] = HGCAL_sigmaUU;
      Var_woGwGSF_VFEndcap_[14] = HGCAL_sigmaVV;
      Var_woGwGSF_VFEndcap_[15] = HGCAL_sigmaEE;
      Var_woGwGSF_VFEndcap_[16] = HGCAL_sigmaPP;
      Var_woGwGSF_VFEndcap_[17] = HGCAL_nLayers;
      Var_woGwGSF_VFEndcap_[18] = HGCAL_lastLayer;
      Var_woGwGSF_VFEndcap_[19] = HGCAL_layerEfrac10;
      Var_woGwGSF_VFEndcap_[20] = HGCAL_layerEfrac90;
      Var_woGwGSF_VFEndcap_[21] = HGCAL_ecEnergyEE;
      Var_woGwGSF_VFEndcap_[22] = HGCAL_ecEnergyFH;
      Var_woGwGSF_VFEndcap_[23] = HGCAL_measuredDepth;
      Var_woGwGSF_VFEndcap_[24] = HGCAL_expectedDepth;
      Var_woGwGSF_VFEndcap_[25] = HGCAL_expectedSigma;
      Var_woGwGSF_VFEndcap_[26] = HGCAL_depthCompatibility;
      Var_woGwGSF_VFEndcap_[27] = ElecDeltaEta;
      Var_woGwGSF_VFEndcap_[28] = ElecDeltaPhi;
      Var_woGwGSF_VFEndcap_[29] = ElecESeedClusterOverPout;
      Var_woGwGSF_VFEndcap_[30] = ElecSuperClusterEtaWidth;
      Var_woGwGSF_VFEndcap_[31] = ElecSuperClusterPhiWidth;
      mvaValue = mva_woGwGSF_VFEC_->GetClassifier(Var_woGwGSF_VFEndcap_); //here
     }			 
  } 
  else if (TauSignalPFGammaCandsIn > 0 && TauHasGsf > 0.5) {
    if (std::abs(TauEtaAtEcalEntrance) < ECALBarrelEndcapEtaBorder) {
      Var_wGwGSF_Barrel_[0] = std::log(std::max(float(0.1), ElecChi2NormGSF));
      Var_wGwGSF_Barrel_[1] = ElecGSFNumHits;
      Var_wGwGSF_Barrel_[2] = std::log(std::max(float(0.1), ElecGSFTrackResol));
      Var_wGwGSF_Barrel_[3] = ElecGSFTracklnPt;
      Var_wGwGSF_Barrel_[4] = ElecNumHitsDiffOverSum;
      Var_wGwGSF_Barrel_[5] = std::log(std::max(float(0.1), ElecChi2NormKF));
      Var_wGwGSF_Barrel_[6] = std::min(ElecDeltaPinPoutOverPin, float(1.));
      Var_wGwGSF_Barrel_[7] = std::min(ElecEecalOverPout, float(20.));
      Var_wGwGSF_Barrel_[8] = std::min(ElecMvaInSigmaEtaEta, float(0.01));
      Var_wGwGSF_Barrel_[9] = std::min(ElecMvaInHadEnergy, float(20.));
      Var_wGwGSF_Barrel_[10] = std::min(ElecMvaInDeltaEta, float(0.1));
      Var_wGwGSF_Barrel_[11] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_wGwGSF_Barrel_[12] = std::log(std::max(float(1.), TauPt));
      Var_wGwGSF_Barrel_[13] = TauEmFraction;
      Var_wGwGSF_Barrel_[14] = TauSignalPFGammaCandsIn;
      Var_wGwGSF_Barrel_[15] = TauSignalPFGammaCandsOut;
      Var_wGwGSF_Barrel_[16] = TauLeadPFChargedHadrHoP;
      Var_wGwGSF_Barrel_[17] = TauLeadPFChargedHadrEoP;
      Var_wGwGSF_Barrel_[18] = TauVisMassIn;
      Var_wGwGSF_Barrel_[19] = TauGammaEtaMomIn;
      Var_wGwGSF_Barrel_[20] = TauGammaEtaMomOut;
      Var_wGwGSF_Barrel_[21] = TauGammaPhiMomIn;
      Var_wGwGSF_Barrel_[22] = TauGammaPhiMomOut;
      Var_wGwGSF_Barrel_[23] = TauGammaEnFracIn;
      Var_wGwGSF_Barrel_[24] = TauGammaEnFracOut;
      Var_wGwGSF_Barrel_[25] = TaudCrackEta;
      Var_wGwGSF_Barrel_[26] = TauEtaAtEcalEntrance;
      Var_wGwGSF_Barrel_[27] = TauLeadChargedPFCandEtaAtEcalEntrance;
      Var_wGwGSF_Barrel_[28] = ElecDeltaEta;
      Var_wGwGSF_Barrel_[29] = ElecDeltaPhi;
      Var_wGwGSF_Barrel_[30] = ElecSigmaIPhiIPhi5x5;
      Var_wGwGSF_Barrel_[31] = ElecSigmaIEtaIEta5x5;
      Var_wGwGSF_Barrel_[32] = ElecShowerCircularity;
      Var_wGwGSF_Barrel_[33] = ElecESeedClusterOverPout;  
      Var_wGwGSF_Barrel_[34] = ElecSuperClusterEtaWidth;  
      Var_wGwGSF_Barrel_[35] = ElecSuperClusterPhiWidth;  
      mvaValue = mva_wGwGSF_BL_->GetClassifier(Var_wGwGSF_Barrel_);
    } else if(std::abs(TauEtaAtEcalEntrance) > ECALBarrelEndcapEtaBorder 
           && std::abs(TauEtaAtEcalEntrance) < ECALBarrelVFEndcapEtaBorder   ) {
      Var_wGwGSF_Endcap_[0] = std::log(std::max(float(0.1), ElecChi2NormGSF));
      Var_wGwGSF_Endcap_[1] = ElecGSFNumHits;
      Var_wGwGSF_Endcap_[2] = std::log(std::max(float(0.1), ElecGSFTrackResol));
      Var_wGwGSF_Endcap_[3] = ElecGSFTracklnPt;
      Var_wGwGSF_Endcap_[4] = ElecNumHitsDiffOverSum;
      Var_wGwGSF_Endcap_[5] = std::log(std::max(float(0.1), ElecChi2NormKF));
      Var_wGwGSF_Endcap_[6] = ElecEecal;
      Var_wGwGSF_Endcap_[7] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_wGwGSF_Endcap_[8] = std::log(std::max(float(1.), TauPt));
      Var_wGwGSF_Endcap_[9] = TauSignalPFGammaCandsIn;
      Var_wGwGSF_Endcap_[10] = TauSignalPFGammaCandsOut;
      Var_wGwGSF_Endcap_[11] = TauVisMassIn;
      Var_wGwGSF_Endcap_[12] = TauGammaEtaMomIn;
      Var_wGwGSF_Endcap_[13] = TauGammaEtaMomOut;
      Var_wGwGSF_Endcap_[14] = TauGammaPhiMomIn;
      Var_wGwGSF_Endcap_[15] = TauGammaPhiMomOut;
      Var_wGwGSF_Endcap_[16] = TauGammaEnFracIn;
      Var_wGwGSF_Endcap_[17] = TauGammaEnFracOut;
      Var_wGwGSF_Endcap_[18] = TaudCrackEta;
      Var_wGwGSF_Endcap_[19] = TauEtaAtEcalEntrance;
      Var_wGwGSF_Endcap_[20] = TauLeadChargedPFCandEtaAtEcalEntrance;
      Var_wGwGSF_Endcap_[21] = HGCAL_sigmaVV;
      Var_wGwGSF_Endcap_[22] = HGCAL_sigmaEE;
      Var_wGwGSF_Endcap_[23] = HGCAL_sigmaPP;
      Var_wGwGSF_Endcap_[24] = HGCAL_nLayers;
      Var_wGwGSF_Endcap_[25] = HGCAL_firstLayer;
      Var_wGwGSF_Endcap_[26] = HGCAL_lastLayer;
      Var_wGwGSF_Endcap_[27] = HGCAL_layerEfrac10;
      Var_wGwGSF_Endcap_[28] = HGCAL_layerEfrac90;
      Var_wGwGSF_Endcap_[29] = HGCAL_ecEnergyEE;
      Var_wGwGSF_Endcap_[30] = HGCAL_ecEnergyFH;
      Var_wGwGSF_Endcap_[31] = HGCAL_measuredDepth;
      Var_wGwGSF_Endcap_[32] = HGCAL_expectedDepth;
      Var_wGwGSF_Endcap_[33] = ElecDeltaEta;
      Var_wGwGSF_Endcap_[34] = ElecDeltaPhi;
      Var_wGwGSF_Endcap_[35] = ElecESeedClusterOverPout;     
      Var_wGwGSF_Endcap_[36] = ElecSuperClusterEtaWidth;
      Var_wGwGSF_Endcap_[37] = ElecSuperClusterPhiWidth; 
      mvaValue = mva_wGwGSF_EC_->GetClassifier(Var_wGwGSF_Endcap_);
    } else if (std::abs(TauEtaAtEcalEntrance) > ECALBarrelVFEndcapEtaBorder
            && std::abs(TauEtaAtEcalEntrance) < HGCALBorder){
      Var_wGwGSF_VFEndcap_[0] = std::log(std::max(float(0.1), ElecChi2NormGSF));
      Var_wGwGSF_VFEndcap_[1] = ElecGSFNumHits;
      Var_wGwGSF_VFEndcap_[2] = std::log(std::max(float(0.1), ElecGSFTrackResol));
      Var_wGwGSF_VFEndcap_[3] = ElecGSFTracklnPt;
      Var_wGwGSF_VFEndcap_[4] = ElecNumHitsDiffOverSum;
      Var_wGwGSF_VFEndcap_[5] = std::log(std::max(float(0.1), ElecChi2NormKF));
      Var_wGwGSF_VFEndcap_[6] = ElecEecal;
      Var_wGwGSF_VFEndcap_[7] = std::min(float(2.), TauLeadChargedPFCandPt / std::max(float(1.), TauPt));
      Var_wGwGSF_VFEndcap_[8] = std::log(std::max(float(1.), TauPt));
      Var_wGwGSF_VFEndcap_[9] = TauSignalPFGammaCandsIn;
      Var_wGwGSF_VFEndcap_[10] = TauSignalPFGammaCandsOut;
      Var_wGwGSF_VFEndcap_[11] = TauVisMassIn;
      Var_wGwGSF_VFEndcap_[12] = TauGammaEtaMomIn;
      Var_wGwGSF_VFEndcap_[13] = TauGammaEtaMomOut;
      Var_wGwGSF_VFEndcap_[14] = TauGammaPhiMomIn;
      Var_wGwGSF_VFEndcap_[15] = TauGammaPhiMomOut;
      Var_wGwGSF_VFEndcap_[16] = TauGammaEnFracIn;
      Var_wGwGSF_VFEndcap_[17] = TauGammaEnFracOut;
      Var_wGwGSF_VFEndcap_[18] = TaudCrackEta;
      Var_wGwGSF_VFEndcap_[19] = TauEtaAtEcalEntrance;
      Var_wGwGSF_VFEndcap_[20] = TauLeadChargedPFCandEtaAtEcalEntrance;
      Var_wGwGSF_VFEndcap_[21] = HGCAL_sigmaUU;
      Var_wGwGSF_VFEndcap_[22] = HGCAL_sigmaVV;
      Var_wGwGSF_VFEndcap_[23] = HGCAL_sigmaEE;
      Var_wGwGSF_VFEndcap_[24] = HGCAL_sigmaPP;
      Var_wGwGSF_VFEndcap_[25] = HGCAL_nLayers;
      Var_wGwGSF_VFEndcap_[26] = HGCAL_lastLayer;
      Var_wGwGSF_VFEndcap_[27] = HGCAL_layerEfrac10;
      Var_wGwGSF_VFEndcap_[28] = HGCAL_layerEfrac90;
      Var_wGwGSF_VFEndcap_[29] = HGCAL_ecEnergyEE;
      Var_wGwGSF_VFEndcap_[30] = HGCAL_ecEnergyFH;
      Var_wGwGSF_VFEndcap_[31] = HGCAL_measuredDepth;
      Var_wGwGSF_VFEndcap_[32] = HGCAL_expectedDepth;
      Var_wGwGSF_VFEndcap_[33] = HGCAL_expectedSigma;
      Var_wGwGSF_VFEndcap_[34] = HGCAL_depthCompatibility;
      Var_wGwGSF_VFEndcap_[35] = ElecDeltaEta;
      Var_wGwGSF_VFEndcap_[36] = ElecDeltaPhi;
      Var_wGwGSF_VFEndcap_[37] = ElecESeedClusterOverPout;     
      Var_wGwGSF_VFEndcap_[38] = ElecSuperClusterEtaWidth;
      Var_wGwGSF_VFEndcap_[39] = ElecSuperClusterPhiWidth; 
      mvaValue = mva_wGwGSF_VFEC_->GetClassifier(Var_wGwGSF_VFEndcap_);
    }
  }

  return mvaValue;
 
}

double AntiElectronIDMVA6::MVAValue(const reco::PFTau& thePFTau, const reco::GsfElectron& theGsfEle, float sigmaUU, float sigmaVV, float sigmaEE, float sigmaPP, float nLayers, float firstLayer, float lastLayer, float layerEfrac10, float layerEfrac90, float ecEnergyEE, float ecEnergyFH, float measuredDepth, float expectedDepth, float expectedSigma, float depthCompatibility)
{
  // === tau variables ===

  float TauEtaAtEcalEntrance = -99.;
   
  float sumEtaTimesEnergy_B = 0.;
  float sumEnergy_B = 0.;
  
  float sumEtaTimesEnergy_F = 0.;
  float sumEnergy_F = 0.;
  
  
  const std::vector<reco::PFCandidatePtr>& signalPFCands = thePFTau.signalPFCands();
  for ( const auto & pfCandidate : signalPFCands ) {    
    
    reco::Candidate const* signalCand = pfCandidate.get();
    float eta_B = pfCandidate->positionAtECALEntrance().eta();
    float eta_F = atECalEntrance(signalCand);
   
    if (eta_B != -99.) sumEtaTimesEnergy_B += eta_B*pfCandidate->energy();
    sumEnergy_B += pfCandidate->energy();
    
    if (eta_F != -99.) sumEtaTimesEnergy_F += eta_F*pfCandidate->energy();
    sumEnergy_F += pfCandidate->energy();
    
  }
  
  if ( fabs(thePFTau.eta())> 1.479 && sumEnergy_F > 0. ) {
    TauEtaAtEcalEntrance = sumEtaTimesEnergy_F/sumEnergy_F;
  }
  else if ( fabs(thePFTau.eta())< 1.479 && sumEnergy_B > 0. ) TauEtaAtEcalEntrance = sumEtaTimesEnergy_B/sumEnergy_B;
    
  
  
  float TauLeadChargedPFCandEtaAtEcalEntrance = -99.;
  float TauLeadChargedPFCandPt = -99.;
  for ( const auto & pfCandidate : signalPFCands ) {
    const reco::Track* track = nullptr;
    if ( pfCandidate->trackRef().isNonnull() ) track = pfCandidate->trackRef().get();
    else if ( pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->innerTrack().isNonnull()  ) track = pfCandidate->muonRef()->innerTrack().get();
    else if ( pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->globalTrack().isNonnull() ) track = pfCandidate->muonRef()->globalTrack().get();
    else if ( pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->outerTrack().isNonnull()  ) track = pfCandidate->muonRef()->outerTrack().get();
    else if ( pfCandidate->gsfTrackRef().isNonnull() ) track = pfCandidate->gsfTrackRef().get();
    if ( track ) {
      if ( track->pt() > TauLeadChargedPFCandPt ) {
	
	reco::Candidate const* signalCand = pfCandidate.get();
	if ( fabs(thePFTau.eta()) > 1.479 ) TauLeadChargedPFCandEtaAtEcalEntrance = atECalEntrance(signalCand);
	else TauLeadChargedPFCandEtaAtEcalEntrance = pfCandidate->positionAtECALEntrance().eta();
	TauLeadChargedPFCandPt = track->pt();
      }
    }
  }

  Float_t TauPt = thePFTau.pt();
  Float_t TauEmFraction = TMath::Max(thePFTau.leadPFChargedHadrCand()->ecalEnergy()/(thePFTau.leadPFChargedHadrCand()->ecalEnergy()+thePFTau.leadPFChargedHadrCand()->hcalEnergy()), float(0.)); 
    
  Float_t TauLeadPFChargedHadrHoP = 0.;
  Float_t TauLeadPFChargedHadrEoP = 0.;
  if ( thePFTau.leadPFChargedHadrCand()->p() > 0. ) {
    TauLeadPFChargedHadrHoP = thePFTau.leadPFChargedHadrCand()->hcalEnergy()/thePFTau.leadPFChargedHadrCand()->p();
    TauLeadPFChargedHadrEoP = thePFTau.leadPFChargedHadrCand()->ecalEnergy()/thePFTau.leadPFChargedHadrCand()->p();
  }

  std::vector<Float_t> GammasdEtaInSigCone;
  std::vector<Float_t> GammasdPhiInSigCone;
  std::vector<Float_t> GammasPtInSigCone;
  std::vector<Float_t> GammasdEtaOutSigCone;
  std::vector<Float_t> GammasdPhiOutSigCone;
  std::vector<Float_t> GammasPtOutSigCone;
  reco::Candidate::LorentzVector pfGammaSum(0,0,0,0);
  reco::Candidate::LorentzVector pfChargedSum(0,0,0,0);
  
  for ( const auto & gamma : thePFTau.signalPFGammaCands() ) {
    float dR = deltaR(gamma->p4(), thePFTau.leadPFChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0/std::max(0.0, thePFTau.pt())));

    // pfGammas inside the tau signal cone
    if (dR < signalrad) {
      if ( thePFTau.leadPFChargedHadrCand().isNonnull() ) {
        GammasdEtaInSigCone.push_back(gamma->eta() - thePFTau.leadPFChargedHadrCand()->eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - thePFTau.leadPFChargedHadrCand()->phi());
      }
      else {
        GammasdEtaInSigCone.push_back(gamma->eta() - thePFTau.eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - thePFTau.phi());
      }
      GammasPtInSigCone.push_back(gamma->pt());
      pfGammaSum += gamma->p4();
    }
    // pfGammas outside the tau signal cone
    else {
      if ( thePFTau.leadPFChargedHadrCand().isNonnull() ) {
        GammasdEtaOutSigCone.push_back(gamma->eta() - thePFTau.leadPFChargedHadrCand()->eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - thePFTau.leadPFChargedHadrCand()->phi());
      } 
      else {
        GammasdEtaOutSigCone.push_back(gamma->eta() - thePFTau.eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - thePFTau.phi());
      }
      GammasPtOutSigCone.push_back(gamma->pt());
    }
  }
  
  for ( const auto & charged : thePFTau.signalPFChargedHadrCands() ) {
    float dR = deltaR(charged->p4(), thePFTau.leadPFChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0/std::max(0.0, thePFTau.pt())));
  
    // charged particles inside the tau signal cone
    if (dR < signalrad) {
        pfChargedSum += charged->p4();
    }
  }
  
  Int_t TauSignalPFGammaCandsIn = GammasPtInSigCone.size();
  Int_t TauSignalPFGammaCandsOut = GammasPtOutSigCone.size();
  Float_t TauVisMassIn = (pfGammaSum + pfChargedSum).mass();

  Float_t TauPhi = thePFTau.phi();

  Float_t TaudCrackEta = dCrackEta(TauEtaAtEcalEntrance);
  Float_t TauHasGsf = 1.;//thePFTau.leadPFChargedHadrCand()->gsfTrackRef().isNonnull();

  
  // === electron variables ===
  
  Float_t ElecEta = theGsfEle.eta();
  Float_t ElecPhi = theGsfEle.phi();
                  			
  // Phase 2 related variables 
  Float_t ElecESeedClusterOverPout = theGsfEle.eSeedClusterOverPout();
  Float_t ElecSuperClusterEtaWidth  = theGsfEle.superCluster()->etaWidth(); 
  Float_t ElecSuperClusterPhiWidth  = theGsfEle.superCluster()->phiWidth();
  Float_t ElecShowerCircularity = 1.-theGsfEle.e1x5()/theGsfEle.e5x5();
  Float_t ElecR9 = theGsfEle.r9();
  Float_t ElecSigmaIEtaIEta5x5 = theGsfEle.full5x5_sigmaIetaIeta();
  Float_t ElecSigmaIPhiIPhi5x5 = theGsfEle.full5x5_sigmaIphiIphi();
  
  Float_t ElecPin = std::sqrt(theGsfEle.trackMomentumAtVtx().Mag2());
  Float_t ElecPout = std::sqrt(theGsfEle.trackMomentumOut().Mag2());  
  Float_t ElecEecal = theGsfEle.ecalEnergy();
  Float_t ElecDeltaEta = theGsfEle.deltaEtaEleClusterTrackAtCalo();
  Float_t ElecDeltaPhi = theGsfEle.deltaPhiEleClusterTrackAtCalo();
  Float_t ElecMvaInSigmaEtaEta = (theGsfEle).mvaInput().sigmaEtaEta;
  Float_t ElecMvaInHadEnergy = (theGsfEle).mvaInput().hadEnergy;
  Float_t ElecMvaInDeltaEta = (theGsfEle).mvaInput().deltaEta;
  
  //Variables related to the GsfTrack
  Float_t ElecChi2NormGSF = -99.;
  Float_t ElecGSFNumHits = -99.;
  Float_t ElecGSFTrackResol = -99.;
  Float_t ElecGSFTracklnPt = -99.;
  if ( theGsfEle.gsfTrack().isNonnull() ) {
    ElecChi2NormGSF = (theGsfEle).gsfTrack()->normalizedChi2();
    ElecGSFNumHits = (theGsfEle).gsfTrack()->numberOfValidHits();
    if ( theGsfEle.gsfTrack()->pt() > 0. ) {
      ElecGSFTrackResol = theGsfEle.gsfTrack()->ptError()/theGsfEle.gsfTrack()->pt();
      ElecGSFTracklnPt = log(theGsfEle.gsfTrack()->pt())*M_LN10;
    }
  }

  //Variables related to the CtfTrack
  Float_t ElecChi2NormKF = -99.;
  Float_t ElecKFNumHits = -99.;
  if ( theGsfEle.closestCtfTrackRef().isNonnull() ) {
    ElecChi2NormKF = (theGsfEle).closestCtfTrackRef()->normalizedChi2();
    ElecKFNumHits = (theGsfEle).closestCtfTrackRef()->numberOfValidHits();
  }


   return MVAValue(TauPt,
   TauEtaAtEcalEntrance,
   TauPhi,
   TauLeadChargedPFCandPt,
   TauLeadChargedPFCandEtaAtEcalEntrance,
   TauEmFraction,
   TauLeadPFChargedHadrHoP,
   TauLeadPFChargedHadrEoP,
   TauVisMassIn,
   TaudCrackEta,
   TauHasGsf,
   TauSignalPFGammaCandsIn,
   TauSignalPFGammaCandsOut,
   GammasdEtaInSigCone, 						     
   GammasdPhiInSigCone, 						   
   GammasPtInSigCone,							   
   GammasdEtaOutSigCone,						   
   GammasdPhiOutSigCone,						   
   GammasPtOutSigCone,
   ElecEta,
   ElecPhi,
   ElecChi2NormGSF,
   ElecChi2NormKF,
   ElecGSFNumHits,
   ElecKFNumHits,
   ElecGSFTrackResol,
   ElecGSFTracklnPt,
   ElecPin,
   ElecPout,
   ElecEecal,
   ElecDeltaEta,
   ElecDeltaPhi,
   ElecMvaInSigmaEtaEta,
   ElecMvaInHadEnergy,
   ElecMvaInDeltaEta,
   sigmaUU,
   sigmaVV,
   sigmaEE,
   sigmaPP,
   nLayers,
   firstLayer,
   lastLayer,
   layerEfrac10,
   layerEfrac90,
   ecEnergyEE,
   ecEnergyFH,
   measuredDepth,
   expectedDepth,
   expectedSigma,
   depthCompatibility,
   ElecESeedClusterOverPout,	 
   ElecSuperClusterEtaWidth,
   ElecSuperClusterPhiWidth, 
   ElecSigmaIEtaIEta5x5,
   ElecSigmaIPhiIPhi5x5,
   ElecShowerCircularity,
   ElecR9
   );	 
}

double AntiElectronIDMVA6::MVAValue(const reco::PFTau& thePFTau)
{
  // === tau variables ===
  float TauEtaAtEcalEntrance = -99.;
  
  float sumEtaTimesEnergy_B = 0.;
  float sumEnergy_B = 0.;
  
  float sumEtaTimesEnergy_F = 0.;
  float sumEnergy_F = 0.;
  
  
  const std::vector<reco::PFCandidatePtr>& signalPFCands = thePFTau.signalPFCands();
  for ( const auto & pfCandidate : signalPFCands ) {    
    
    reco::Candidate const* signalCand = pfCandidate.get();
    float eta_B = pfCandidate->positionAtECALEntrance().eta();
    float eta_F = atECalEntrance(signalCand);
   
    if (eta_B != -99.) sumEtaTimesEnergy_B += eta_B*pfCandidate->energy();
    sumEnergy_B += pfCandidate->energy();
    
    if (eta_F != -99.) sumEtaTimesEnergy_F += eta_F*pfCandidate->energy();
    sumEnergy_F += pfCandidate->energy();
    
  }
  
  if ( fabs(thePFTau.eta())> 1.479 && sumEnergy_F > 0. ) {
    TauEtaAtEcalEntrance = sumEtaTimesEnergy_F/sumEnergy_F;
  }
  else if ( fabs(thePFTau.eta())< 1.479 && sumEnergy_B > 0. ) TauEtaAtEcalEntrance = sumEtaTimesEnergy_B/sumEnergy_B;
    
  
  
  float TauLeadChargedPFCandEtaAtEcalEntrance = -99.;
  float TauLeadChargedPFCandPt = -99.;
  for ( const auto & pfCandidate : signalPFCands ) {
    const reco::Track* track = nullptr;
    if ( pfCandidate->trackRef().isNonnull() ) track = pfCandidate->trackRef().get();
    else if ( pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->innerTrack().isNonnull()  ) track = pfCandidate->muonRef()->innerTrack().get();
    else if ( pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->globalTrack().isNonnull() ) track = pfCandidate->muonRef()->globalTrack().get();
    else if ( pfCandidate->muonRef().isNonnull() && pfCandidate->muonRef()->outerTrack().isNonnull()  ) track = pfCandidate->muonRef()->outerTrack().get();
    else if ( pfCandidate->gsfTrackRef().isNonnull() ) track = pfCandidate->gsfTrackRef().get();
    if ( track ) {
      if ( track->pt() > TauLeadChargedPFCandPt ) {
	 
	 reco::Candidate const* signalCand = pfCandidate.get();
	 
	if ( fabs(thePFTau.eta()) > 1.479 ) TauLeadChargedPFCandEtaAtEcalEntrance = atECalEntrance(signalCand);
	else TauLeadChargedPFCandEtaAtEcalEntrance = pfCandidate->positionAtECALEntrance().eta();
	TauLeadChargedPFCandPt = track->pt();
      }
    }
  }

  Float_t TauPt = thePFTau.pt();
  Float_t TauEmFraction = TMath::Max(thePFTau.leadPFChargedHadrCand()->ecalEnergy()/(thePFTau.leadPFChargedHadrCand()->ecalEnergy()+thePFTau.leadPFChargedHadrCand()->hcalEnergy()), float(0.));
  
  Float_t TauLeadPFChargedHadrHoP = 0.;
  Float_t TauLeadPFChargedHadrEoP = 0.;
  if ( thePFTau.leadPFChargedHadrCand()->p() > 0. ) {
    TauLeadPFChargedHadrHoP = thePFTau.leadPFChargedHadrCand()->hcalEnergy()/thePFTau.leadPFChargedHadrCand()->p();
    TauLeadPFChargedHadrEoP = thePFTau.leadPFChargedHadrCand()->ecalEnergy()/thePFTau.leadPFChargedHadrCand()->p();
  }

  std::vector<Float_t> GammasdEtaInSigCone;
  std::vector<Float_t> GammasdPhiInSigCone;
  std::vector<Float_t> GammasPtInSigCone;
  std::vector<Float_t> GammasdEtaOutSigCone;
  std::vector<Float_t> GammasdPhiOutSigCone;
  std::vector<Float_t> GammasPtOutSigCone;
  reco::Candidate::LorentzVector pfGammaSum(0,0,0,0);
  reco::Candidate::LorentzVector pfChargedSum(0,0,0,0);
  
  for ( const auto & gamma : thePFTau.signalPFGammaCands() ) {
    float dR = deltaR(gamma->p4(), thePFTau.leadPFChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0/std::max(0.0, thePFTau.pt())));

    // pfGammas inside the tau signal cone
    if (dR < signalrad) {
      if ( thePFTau.leadPFChargedHadrCand().isNonnull() ) {
        GammasdEtaInSigCone.push_back(gamma->eta() - thePFTau.leadPFChargedHadrCand()->eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - thePFTau.leadPFChargedHadrCand()->phi());
      }
      else {
        GammasdEtaInSigCone.push_back(gamma->eta() - thePFTau.eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - thePFTau.phi());
      }
      GammasPtInSigCone.push_back(gamma->pt());
      pfGammaSum += gamma->p4();
    }
    // pfGammas outside the tau signal cone
    else {
      if ( thePFTau.leadPFChargedHadrCand().isNonnull() ) {
        GammasdEtaOutSigCone.push_back(gamma->eta() - thePFTau.leadPFChargedHadrCand()->eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - thePFTau.leadPFChargedHadrCand()->phi());
      } 
      else {
        GammasdEtaOutSigCone.push_back(gamma->eta() - thePFTau.eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - thePFTau.phi());
      }
      GammasPtOutSigCone.push_back(gamma->pt());
    }
  }
  
  for ( const auto & charged : thePFTau.signalPFChargedHadrCands() ) {
    float dR = deltaR(charged->p4(), thePFTau.leadPFChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0/std::max(0.0, thePFTau.pt())));
  
    // charged particles inside the tau signal cone
    if (dR < signalrad) {
        pfChargedSum += charged->p4();
    }
  }
  
  Int_t TauSignalPFGammaCandsIn = GammasPtInSigCone.size();
  Int_t TauSignalPFGammaCandsOut = GammasPtOutSigCone.size();
  Float_t TauVisMassIn = (pfGammaSum + pfChargedSum).mass();

  Float_t TauPhi = thePFTau.phi();
  
  Float_t TaudCrackEta = dCrackEta(TauEtaAtEcalEntrance);
  Float_t TauHasGsf = thePFTau.leadPFChargedHadrCand()->gsfTrackRef().isNonnull();

  
  // === electron variables ===
  
  Float_t dummyElecEta = 9.9;
    
  return MVAValue(TauPt,
                  TauEtaAtEcalEntrance,
                  TauPhi,
                  TauLeadChargedPFCandPt,
                  TauLeadChargedPFCandEtaAtEcalEntrance,
                  TauEmFraction,
                  TauLeadPFChargedHadrHoP,
                  TauLeadPFChargedHadrEoP,
                  TauVisMassIn,
                  TaudCrackEta,
                  TauHasGsf,
                  TauSignalPFGammaCandsIn,
                  TauSignalPFGammaCandsOut,
                  GammasdEtaInSigCone,
                  GammasdPhiInSigCone,
                  GammasPtInSigCone,
                  GammasdEtaOutSigCone,
                  GammasdPhiOutSigCone,
                  GammasPtOutSigCone,
                  dummyElecEta,           
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.,
        	  0.
		  );
}

double AntiElectronIDMVA6::MVAValue(const pat::Tau& theTau, const pat::Electron& theEle)
{
  // === tau variables ===
  
  pat::PackedCandidate const* packedLeadTauCand =
       dynamic_cast<pat::PackedCandidate const*>(theTau.leadChargedHadrCand().get());

  float TauLeadChargedPFCandPt = theTau.ptLeadChargedCand();

  Float_t TauPt = theTau.pt();
  Float_t TauEmFraction = TMath::Max(theTau.ecalEnergyLeadChargedHadrCand()/(theTau.hcalEnergyLeadChargedHadrCand()+theTau.ecalEnergyLeadChargedHadrCand()), float(0.));
  
  Float_t TauLeadPFChargedHadrHoP = 0.;
  Float_t TauLeadPFChargedHadrEoP = 0.;
  if ( theTau.leadChargedHadrCand()->p() > 0. ) {
    TauLeadPFChargedHadrHoP = theTau.hcalEnergyLeadChargedHadrCand()/theTau.leadChargedHadrCand()->p();
    TauLeadPFChargedHadrEoP = theTau.ecalEnergyLeadChargedHadrCand()/theTau.leadChargedHadrCand()->p();
  }
  
  std::vector<Float_t> GammasdEtaInSigCone;
  std::vector<Float_t> GammasdPhiInSigCone;
  std::vector<Float_t> GammasPtInSigCone;
  std::vector<Float_t> GammasdEtaOutSigCone;
  std::vector<Float_t> GammasdPhiOutSigCone;
  std::vector<Float_t> GammasPtOutSigCone;
  reco::Candidate::LorentzVector pfGammaSum(0, 0, 0, 0);
  reco::Candidate::LorentzVector pfChargedSum(0, 0, 0, 0);

  const reco::CandidatePtrVector signalGammaCands = theTau.signalGammaCands();
  for (const auto& gamma : signalGammaCands) {
    float dR = deltaR(gamma->p4(), theTau.leadChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0 / std::max(0.0, theTau.pt())));

    // pfGammas inside the tau signal cone
    if (dR < signalrad) {
      if (theTau.leadChargedHadrCand().isNonnull()) {
        GammasdEtaInSigCone.push_back(gamma->eta() - theTau.leadChargedHadrCand()->eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - theTau.leadChargedHadrCand()->phi());
        //A.-C. please check whether this change is safe against future trainings
        //GammasdPhiInSigCone.push_back(deltaPhi((*gamma)->phi(), theTau.leadChargedHadrCand()->phi()));
      } else {
        GammasdEtaInSigCone.push_back(gamma->eta() - theTau.eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - theTau.phi());
        //A.-C. please check whether this change is safe against future trainings
        //GammasdPhiInSigCone.push_back(deltaPhi(gamma->phi(), theTau.phi()));
      }
      GammasPtInSigCone.push_back(gamma->pt());
      pfGammaSum += gamma->p4();
    } 
    // pfGammas outside the tau signal cone
    else {
      if (theTau.leadChargedHadrCand().isNonnull()) {
        GammasdEtaOutSigCone.push_back(gamma->eta() - theTau.leadChargedHadrCand()->eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - theTau.leadChargedHadrCand()->phi());
        //A.-C. please check whether this change is safe against future trainings
        //GammasdPhiOutSigCone.push_back(deltaPhi(gamma->phi(), theTau.leadChargedHadrCand()->phi()));
      } else {
        GammasdEtaOutSigCone.push_back(gamma->eta() - theTau.eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - theTau.phi());
        //A.-C. please chaekc whether this change is safe against future trainings
        //GammasdPhiOutSigCone.push_back(deltaPhi(gamma->phi(), theTau.phi()));
      }
      GammasPtOutSigCone.push_back(gamma->pt());
    }
  }
 
  const reco::CandidatePtrVector signalChargedCands = theTau.signalChargedHadrCands();
  for (const auto& charged : signalChargedCands) {
    float dR = deltaR(charged->p4(), theTau.leadChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0 / std::max(0.0, theTau.pt())));

    // charged particles inside the tau signal cone
    if (dR < signalrad) {
      pfChargedSum += charged->p4();
    }
  }

  Int_t TauSignalPFGammaCandsIn = GammasPtInSigCone.size();
  Int_t TauSignalPFGammaCandsOut = GammasPtOutSigCone.size();
  Float_t TauVisMassIn = (pfGammaSum + pfChargedSum).mass();
  
  Float_t TauPhi = theTau.phi();
  
  const reco::CandidatePtrVector signalCands = theTau.signalCands();
  float sumEtaTimesEnergy = 0.;
  float sumEnergy = 0.;
  //const reco::CandidatePtrVector signalCands = theTau.signalCands();
  for (const auto& signalCandPtr : signalCands) {
    reco::Candidate const* signalCand = signalCandPtr.get();

    float eta = atECalEntrance(signalCand);
    if (eta != -99.) sumEtaTimesEnergy += eta * signalCand->energy();   
    sumEnergy += signalCand->energy();
    
  }
     
   Float_t TauEtaAtEcalEntrance = -99.;
   if ( fabs(theTau.eta())> 1.479 && sumEnergy > 0.) TauEtaAtEcalEntrance = sumEtaTimesEnergy/sumEnergy;
   else TauEtaAtEcalEntrance = theTau.etaAtEcalEntrance();
  
   
   Float_t TauLeadChargedPFCandEtaAtEcalEntrance = -99.;
   if ( fabs(theTau.eta()) > 1.479 ) TauLeadChargedPFCandEtaAtEcalEntrance = atECalEntrance(packedLeadTauCand);    
   else TauLeadChargedPFCandEtaAtEcalEntrance = theTau.etaAtEcalEntranceLeadChargedCand();
     
     
  Float_t TaudCrackEta = dCrackEta(TauEtaAtEcalEntrance); 
  
  Float_t TauHasGsf = 1.;//0;
  //if( abs(packedLeadTauCand->pdgId()) == 11 ) TauHasGsf = 1;
  
  // === electron variables ===
  
  Float_t ElecEta = theEle.eta();
  Float_t ElecPhi = theEle.phi();
  
           
  //Phase2 related variables
  
  Float_t ElecESeedClusterOverPout = (theEle).eSeedClusterOverPout();
  Float_t ElecSuperClusterEtaWidth  = (theEle).superCluster()->etaWidth(); 
  Float_t ElecSuperClusterPhiWidth  = (theEle).superCluster()->phiWidth();
  Float_t ElecShowerCircularity = 1.-(theEle).e1x5()/(theEle).e5x5();
  Float_t ElecR9 = (theEle).r9();
  Float_t ElecSigmaIEtaIEta5x5 = (theEle).full5x5_sigmaIetaIeta();
  Float_t ElecSigmaIPhiIPhi5x5 = (theEle).full5x5_sigmaIphiIphi();
 
  Float_t HGCAL_sigmaUU = (theEle).userFloat("hgcElectronID:sigmaUU");
  Float_t HGCAL_sigmaVV = (theEle).userFloat("hgcElectronID:sigmaVV");
  Float_t HGCAL_sigmaEE = (theEle).userFloat("hgcElectronID:sigmaEE");
  Float_t HGCAL_sigmaPP = (theEle).userFloat("hgcElectronID:sigmaPP");
  Float_t HGCAL_nLayers = (theEle).userFloat("hgcElectronID:nLayers");
  Float_t HGCAL_firstLayer = (theEle).userFloat("hgcElectronID:firstLayer");
  Float_t HGCAL_lastLayer = (theEle).userFloat("hgcElectronID:lastLayer");
  Float_t HGCAL_layerEfrac10 = (theEle).userFloat("hgcElectronID:layerEfrac10");
  Float_t HGCAL_layerEfrac90 = (theEle).userFloat("hgcElectronID:layerEfrac90");
  Float_t HGCAL_ecEnergyEE = (theEle).userFloat("hgcElectronID:ecEnergyEE");
  Float_t HGCAL_ecEnergyFH = (theEle).userFloat("hgcElectronID:ecEnergyFH");
  Float_t HGCAL_measuredDepth = (theEle).userFloat("hgcElectronID:measuredDepth");
  Float_t HGCAL_expectedDepth = (theEle).userFloat("hgcElectronID:expectedDepth");
  Float_t HGCAL_expectedSigma = (theEle).userFloat("hgcElectronID:expectedSigma");
  Float_t HGCAL_depthCompatibility = (theEle).userFloat("hgcElectronID:depthCompatibility");
  
  
  Float_t ElecPin = std::sqrt(theEle.trackMomentumAtVtx().Mag2());
  Float_t ElecPout = std::sqrt(theEle.trackMomentumOut().Mag2());  
  Float_t ElecEecal = theEle.ecalEnergy();
  Float_t ElecDeltaEta = theEle.deltaEtaEleClusterTrackAtCalo();
  Float_t ElecDeltaPhi = theEle.deltaPhiEleClusterTrackAtCalo();
  Float_t ElecMvaInSigmaEtaEta = (theEle).mvaInput().sigmaEtaEta;
  Float_t ElecMvaInHadEnergy = (theEle).mvaInput().hadEnergy;
  Float_t ElecMvaInDeltaEta = (theEle).mvaInput().deltaEta;
  
  //Variables related to the GsfTrack
  Float_t ElecChi2NormGSF = -99.;
  Float_t ElecGSFNumHits = -99.;
  Float_t ElecGSFTrackResol = -99.;
  Float_t ElecGSFTracklnPt = -99.;
  if ( theEle.gsfTrack().isNonnull() ) {
    ElecChi2NormGSF = (theEle).gsfTrack()->normalizedChi2();
    ElecGSFNumHits = (theEle).gsfTrack()->numberOfValidHits();
    if ( theEle.gsfTrack()->pt() > 0. ) {
      ElecGSFTrackResol = theEle.gsfTrack()->ptError()/theEle.gsfTrack()->pt();
      ElecGSFTracklnPt = log(theEle.gsfTrack()->pt())*M_LN10;
    }
  }

  //Variables related to the CtfTrack
  Float_t ElecChi2NormKF = -99.;
  Float_t ElecKFNumHits = -99.;
  if ( theEle.closestCtfTrackRef().isNonnull() ) {
    ElecChi2NormKF = (theEle).closestCtfTrackRef()->normalizedChi2();
    ElecKFNumHits = (theEle).closestCtfTrackRef()->numberOfValidHits();
  }
  
  return MVAValue(TauPt,
                  TauEtaAtEcalEntrance,
                  TauPhi,
                  TauLeadChargedPFCandPt,
                  TauLeadChargedPFCandEtaAtEcalEntrance,
                  TauEmFraction,
                  TauLeadPFChargedHadrHoP,
                  TauLeadPFChargedHadrEoP,
                  TauVisMassIn,
                  TaudCrackEta,
                  TauHasGsf,
                  TauSignalPFGammaCandsIn,
                  TauSignalPFGammaCandsOut,
                  GammasdEtaInSigCone,  						    
                  GammasdPhiInSigCone,  						  
                  GammasPtInSigCone,							  
                  GammasdEtaOutSigCone, 						  
                  GammasdPhiOutSigCone, 						  
                  GammasPtOutSigCone,
                  ElecEta,
                  ElecPhi,
                  ElecChi2NormGSF,
                  ElecChi2NormKF,
                  ElecGSFNumHits,
                  ElecKFNumHits,
                  ElecGSFTrackResol,
                  ElecGSFTracklnPt,
                  ElecPin,
                  ElecPout,
                  ElecEecal,
                  ElecDeltaEta,
                  ElecDeltaPhi,
                  ElecMvaInSigmaEtaEta,
                  ElecMvaInHadEnergy,
                  ElecMvaInDeltaEta,
		  HGCAL_sigmaUU,
     	 	  HGCAL_sigmaVV,
     	 	  HGCAL_sigmaEE,
     	 	  HGCAL_sigmaPP,
     	 	  HGCAL_nLayers,
		  HGCAL_firstLayer,
     	 	  HGCAL_lastLayer,
     	 	  HGCAL_layerEfrac10,
     	 	  HGCAL_layerEfrac90,
     	 	  HGCAL_ecEnergyEE,
     	 	  HGCAL_ecEnergyFH,
     	 	  HGCAL_measuredDepth,
     	 	  HGCAL_expectedDepth,
     	 	  HGCAL_expectedSigma,
     	 	  HGCAL_depthCompatibility,
     	 	  ElecESeedClusterOverPout,	
     	 	  ElecSuperClusterEtaWidth,
     	 	  ElecSuperClusterPhiWidth, 
    	 	  ElecSigmaIEtaIEta5x5,
     	 	  ElecSigmaIPhiIPhi5x5,
     	 	  ElecShowerCircularity,
     	 	  ElecR9
		  );		
}

double AntiElectronIDMVA6::MVAValue(const pat::Tau& theTau)
{
  // === tau variables ===
  
  pat::PackedCandidate const* packedLeadTauCand =
       dynamic_cast<pat::PackedCandidate const*>(theTau.leadChargedHadrCand().get());
  
  float TauLeadChargedPFCandPt = theTau.ptLeadChargedCand();

  Float_t TauPt = theTau.pt();
  Float_t TauEmFraction = TMath::Max(theTau.ecalEnergyLeadChargedHadrCand()/(theTau.hcalEnergyLeadChargedHadrCand()+theTau.ecalEnergyLeadChargedHadrCand()), float(0.));
  
  Float_t TauLeadPFChargedHadrHoP = 0.;
  Float_t TauLeadPFChargedHadrEoP = 0.;
  if ( theTau.leadChargedHadrCand()->p() > 0. ) {
    TauLeadPFChargedHadrHoP = theTau.hcalEnergyLeadChargedHadrCand()/theTau.leadChargedHadrCand()->p();
    TauLeadPFChargedHadrEoP = theTau.ecalEnergyLeadChargedHadrCand()/theTau.leadChargedHadrCand()->p();
  }
 

  std::vector<Float_t> GammasdEtaInSigCone;
  std::vector<Float_t> GammasdPhiInSigCone;
  std::vector<Float_t> GammasPtInSigCone;
  std::vector<Float_t> GammasdEtaOutSigCone;
  std::vector<Float_t> GammasdPhiOutSigCone;
  std::vector<Float_t> GammasPtOutSigCone;
  reco::Candidate::LorentzVector pfGammaSum(0, 0, 0, 0);
  reco::Candidate::LorentzVector pfChargedSum(0, 0, 0, 0);

  const reco::CandidatePtrVector signalGammaCands = theTau.signalGammaCands();
  for (const auto& gamma : signalGammaCands) {
    float dR = deltaR(gamma->p4(), theTau.leadChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0 / std::max(0.0, theTau.pt())));

    // pfGammas inside the tau signal cone
    if (dR < signalrad) {
      if (theTau.leadChargedHadrCand().isNonnull()) {
        GammasdEtaInSigCone.push_back(gamma->eta() - theTau.leadChargedHadrCand()->eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - theTau.leadChargedHadrCand()->phi());
        //A.-C. please check whether this change is safe against future trainings
        //GammasdPhiInSigCone.push_back(deltaPhi((*gamma)->phi(), theTau.leadChargedHadrCand()->phi()));
      } else {
        GammasdEtaInSigCone.push_back(gamma->eta() - theTau.eta());
        GammasdPhiInSigCone.push_back(gamma->phi() - theTau.phi());
        //A.-C. please check whether this change is safe against future trainings
        //GammasdPhiInSigCone.push_back(deltaPhi(gamma->phi(), theTau.phi()));
      }
      GammasPtInSigCone.push_back(gamma->pt());
      pfGammaSum += gamma->p4();
    } 
    // pfGammas outside the tau signal cone
    else {
      if (theTau.leadChargedHadrCand().isNonnull()) {
        GammasdEtaOutSigCone.push_back(gamma->eta() - theTau.leadChargedHadrCand()->eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - theTau.leadChargedHadrCand()->phi());
        //A.-C. please check whether this change is safe against future trainings
        //GammasdPhiOutSigCone.push_back(deltaPhi(gamma->phi(), theTau.leadChargedHadrCand()->phi()));
      } else {
        GammasdEtaOutSigCone.push_back(gamma->eta() - theTau.eta());
        GammasdPhiOutSigCone.push_back(gamma->phi() - theTau.phi());
        //A.-C. please chaekc whether this change is safe against future trainings
        //GammasdPhiOutSigCone.push_back(deltaPhi(gamma->phi(), theTau.phi()));
      }
      GammasPtOutSigCone.push_back(gamma->pt());
    }
  }
 
  const reco::CandidatePtrVector signalChargedCands = theTau.signalChargedHadrCands();
  for (const auto& charged : signalChargedCands) {
    float dR = deltaR(charged->p4(), theTau.leadChargedHadrCand()->p4());
    float signalrad = std::max(0.05, std::min(0.10, 3.0 / std::max(0.0, theTau.pt())));

    // charged particles inside the tau signal cone
    if (dR < signalrad) {
      pfChargedSum += charged->p4();
    }
  }

  Int_t TauSignalPFGammaCandsIn = GammasPtInSigCone.size();
  Int_t TauSignalPFGammaCandsOut = GammasPtOutSigCone.size();
  Float_t TauVisMassIn = (pfGammaSum + pfChargedSum).mass();
  Float_t TauPhi = theTau.phi();
  
  const reco::CandidatePtrVector signalCands = theTau.signalCands();
  //float sumPhiTimesEnergy = 0.;
  float sumEtaTimesEnergy = 0.;
  float sumEnergy = 0.;
  //const reco::CandidatePtrVector signalCands = theTau.signalCands();
  for (const auto& signalCandPtr : signalCands) {
    reco::Candidate const* signalCand = signalCandPtr.get();
 
    float eta = atECalEntrance(signalCand);
    if (eta != -99.) sumEtaTimesEnergy += eta * signalCand->energy(); 
    
    sumEnergy += signalCand->energy();
    
  }
  
   Float_t TauEtaAtEcalEntrance = -99.;
   if ( fabs(theTau.eta())> 1.479 && sumEnergy > 0.) TauEtaAtEcalEntrance = sumEtaTimesEnergy/sumEnergy;
   else TauEtaAtEcalEntrance = theTau.etaAtEcalEntrance();
  
   
    Float_t TauLeadChargedPFCandEtaAtEcalEntrance = -99.;
    if ( fabs(theTau.eta()) > 1.479 ) TauLeadChargedPFCandEtaAtEcalEntrance = atECalEntrance(packedLeadTauCand);   
    else 
     TauLeadChargedPFCandEtaAtEcalEntrance = theTau.etaAtEcalEntranceLeadChargedCand();
     
  Float_t TaudCrackEta = dCrackEta(TauEtaAtEcalEntrance); 
  
  Float_t TauHasGsf = 0.;
  
 
  // === electron variables ===
  Float_t dummyElecEta = 9.9;
 
  return MVAValue(TauPt,                                                 
                  TauEtaAtEcalEntrance, 			         
                  TauPhi,					         
                  TauLeadChargedPFCandPt,			         
                  TauLeadChargedPFCandEtaAtEcalEntrance,	        
                  TauEmFraction,				         
                  TauLeadPFChargedHadrHoP,			         
                  TauLeadPFChargedHadrEoP,			         
                  TauVisMassIn, 				         
                  TaudCrackEta, 				         
                  TauHasGsf,					         
                  TauSignalPFGammaCandsIn,			         
                  TauSignalPFGammaCandsOut,			         
                  GammasdEtaInSigCone,				         
                  GammasdPhiInSigCone,				         
                  GammasPtInSigCone,				         
                  GammasdEtaOutSigCone, 			         
                  GammasdPhiOutSigCone, 			         
                  GammasPtOutSigCone,				         
                  dummyElecEta, 				         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
                  0.,						         
		  0.,							
                  0.,						     	 
                  0.,						     	 
     	 	  0.,						     	 
     	 	  0.,						     	 
     	 	  0.,							
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,						     	 	  
     	 	  0.,   					     	 	  
     	 	  0.,						     	 	  
     	 	  0., 						     	 	  
     	 	  0.
		  );
}

double AntiElectronIDMVA6::minimum(double a, double b)
{
  if ( std::abs(b) < std::abs(a) ) return b;
  else return a;
}

namespace {

  // IN: define locations of the 18 phi-cracks
 std::array<double,18> fill_cPhi() {
   constexpr double pi = M_PI; // 3.14159265358979323846;
   std::array<double,18> cPhi;
   // IN: define locations of the 18 phi-cracks
   cPhi[0] = 2.97025;
   for ( unsigned iCrack = 1; iCrack <= 17; ++iCrack ) {
     cPhi[iCrack] = cPhi[0] - 2.*iCrack*pi/18;
   }
   return cPhi;
 }
     
  const std::array<double,18> cPhi = fill_cPhi();

}

double AntiElectronIDMVA6::dCrackPhi(double phi, double eta)
{
//--- compute the (unsigned) distance to the closest phi-crack in the ECAL barrel  

  constexpr double pi = M_PI; // 3.14159265358979323846;

  // IN: shift of this location if eta < 0
  constexpr double delta_cPhi = 0.00638;

  double retVal = 99.; 

  if ( eta >= -1.47464 && eta <= 1.47464 ) {

    // the location is shifted
    if ( eta < 0. ) phi += delta_cPhi;

    // CV: need to bring-back phi into interval [-pi,+pi]
    if ( phi >  pi ) phi -= 2.*pi;
    if ( phi < -pi ) phi += 2.*pi;

    if ( phi >= -pi && phi <= pi ) {

      // the problem of the extrema:
      if ( phi < cPhi[17] || phi >= cPhi[0] ) {
	if ( phi < 0. ) phi += 2.*pi;
	retVal = minimum(phi - cPhi[0], phi - cPhi[17] - 2.*pi);        	
      } else {
	// between these extrema...
	bool OK = false;
	unsigned iCrack = 16;
	while( !OK ) {
	  if ( phi < cPhi[iCrack] ) {
	    retVal = minimum(phi - cPhi[iCrack + 1], phi - cPhi[iCrack]);
	    OK = true;
	  } else {
	    iCrack -= 1;
	  }
	}
      }
    } else {
      retVal = 0.; // IN: if there is a problem, we assume that we are in a crack
    }
  } else {
    return -99.;       
  }
  
  return std::abs(retVal);
}

double AntiElectronIDMVA6::dCrackEta(double eta)
{
//--- compute the (unsigned) distance to the closest eta-crack in the ECAL barrel
  
  // IN: define locations of the eta-cracks
  double cracks[5] = { 0., 4.44747e-01, 7.92824e-01, 1.14090e+00, 1.47464e+00 };
  
  double retVal = 99.;
  
  for ( int iCrack = 0; iCrack < 5 ; ++iCrack ) {
    double d = minimum(eta - cracks[iCrack], eta + cracks[iCrack]);
    if ( std::abs(d) < std::abs(retVal) ) {
      retVal = d;
    }
  }

  return std::abs(retVal);
}

float AntiElectronIDMVA6::atECalEntrance(const pat::PackedCandidate &part)
{

  float result = -99.;
  
  HGCal_helpers::simpleTrackPropagator toHGCalPropagator(aField_);
  toHGCalPropagator.setPropagationTargetZ(layerPositions_[0]);
  
  
  if (std::abs(part.vertex().z()) >= layerPositions_[0]) return result;
  
  HGCal_helpers::coordinates hitsHGCal;
  
  
  toHGCalPropagator.propagate(part.px(),part.py(),part.pz(),
                              part.vertex().x(),part.vertex().y(),part.vertex().z(),
			      part.charge(),
			      hitsHGCal);
  
  result = hitsHGCal.eta;

  return result;
}

float AntiElectronIDMVA6::atECalEntrance(const reco::Candidate* part) {
  
  float result = -99.;
  
  HGCal_helpers::simpleTrackPropagator toHGCalPropagator(aField_);
  toHGCalPropagator.setPropagationTargetZ(layerPositions_[0]);
  
  
  if (std::abs(part->vertex().z()) >= layerPositions_[0]) return result;
  
  HGCal_helpers::coordinates hitsHGCal;
  
  
  toHGCalPropagator.propagate(part->px(),part->py(),part->pz(),
                              part->vertex().x(),part->vertex().y(),part->vertex().z(),
			      part->charge(),
			      hitsHGCal);
  
  result = hitsHGCal.eta;
  
  return result;
}

