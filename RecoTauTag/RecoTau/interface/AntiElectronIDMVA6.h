//--------------------------------------------------------------------------------------------------
// AntiElectronIDMVA6
//
// Helper Class for applying MVA anti-electron discrimination
//
// Authors: F.Colombo, C.Veelken
//--------------------------------------------------------------------------------------------------

#ifndef RECOTAUTAG_RECOTAU_AntiElectronIDMVA6_H
#define RECOTAUTAG_RECOTAU_AntiElectronIDMVA6_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
//#include "FastSimulation/Particle/interface/RawParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

#include "DataFormats/Common/interface/ValueMap.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeGeometry/interface/MagVolumeOutsideValidity.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


class AntiElectronIDMVA6 
{
  public:

   AntiElectronIDMVA6(const edm::ParameterSet&);
   ~AntiElectronIDMVA6(); 

   void beginEvent(const edm::Event&, const edm::EventSetup&);
   
   double MVAValue(Float_t TauPt,
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
		  );

  double MVAValue(Float_t TauPt,
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
		  );

   // this function can be called for all categories
   double MVAValue(const reco::PFTau& thePFTau, const reco::GsfElectron& theGsfEle, float sigmaUU, float sigmaVV, float sigmaEE, float sigmaPP, float nLayers, float firstLayer, float lastLayer, float layerEfrac10, float layerEfrac90, float ecEnergyEE, float ecEnergyFH, float measuredDepth, float expectedDepth, float expectedSigma, float depthCompatibility);
   // this function can be called for category 1 only !!
   double MVAValue(const reco::PFTau& thePFTau);

   // this function can be called for all categories
   double MVAValue(const pat::Tau& theTau, 
		   const pat::Electron& theEle);
   // this function can be called for category 1 only !!
   double MVAValue(const pat::Tau& theTau);
   // track extrapolation to ECAL entrance (used to re-calculate varibales that might not be available on miniAOD)
  
   float atECalEntrance(const pat::PackedCandidate &part);
   float atECalEntrance(const reco::Candidate* part);

 private:   

   double dCrackEta(double eta);
   double minimum(double a,double b);
   double dCrackPhi(double phi, double eta);

   bool isInitialized_;
   bool loadMVAfromDB_;
   edm::FileInPath inputFileName_;
   
   std::string mvaName_NoEleMatch_woGwoGSF_BL_;
   std::string mvaName_NoEleMatch_wGwoGSF_BL_;
   std::string mvaName_woGwGSF_BL_;
   std::string mvaName_wGwGSF_BL_;
   std::string mvaName_NoEleMatch_woGwoGSF_EC_;
   std::string mvaName_NoEleMatch_wGwoGSF_EC_;
   std::string mvaName_woGwGSF_EC_;
   std::string mvaName_wGwGSF_EC_;
   std::string mvaName_NoEleMatch_woGwoGSF_VFEC_;
   std::string mvaName_NoEleMatch_wGwoGSF_VFEC_;
   std::string mvaName_woGwGSF_VFEC_;
   std::string mvaName_wGwGSF_VFEC_;

   //bool usePhiAtEcalEntranceExtrapolation_;
   
   Float_t* Var_NoEleMatch_woGwoGSF_Barrel_;
   Float_t* Var_NoEleMatch_wGwoGSF_Barrel_;
   Float_t* Var_woGwGSF_Barrel_;
   Float_t* Var_wGwGSF_Barrel_;
   Float_t* Var_NoEleMatch_woGwoGSF_Endcap_;
   Float_t* Var_NoEleMatch_wGwoGSF_Endcap_;
   Float_t* Var_woGwGSF_Endcap_;
   Float_t* Var_wGwGSF_Endcap_;
   Float_t* Var_NoEleMatch_woGwoGSF_VFEndcap_;
   Float_t* Var_NoEleMatch_wGwoGSF_VFEndcap_;
   Float_t* Var_woGwGSF_VFEndcap_;
   Float_t* Var_wGwGSF_VFEndcap_;

   const GBRForest* mva_NoEleMatch_woGwoGSF_BL_;
   const GBRForest* mva_NoEleMatch_wGwoGSF_BL_;
   const GBRForest* mva_woGwGSF_BL_;
   const GBRForest* mva_wGwGSF_BL_;
   const GBRForest* mva_NoEleMatch_woGwoGSF_EC_;
   const GBRForest* mva_NoEleMatch_wGwoGSF_EC_;
   const GBRForest* mva_woGwGSF_EC_;
   const GBRForest* mva_wGwGSF_EC_;
   const GBRForest* mva_NoEleMatch_woGwoGSF_VFEC_;
   const GBRForest* mva_NoEleMatch_wGwoGSF_VFEC_;
   const GBRForest* mva_woGwGSF_VFEC_;
   const GBRForest* mva_wGwGSF_VFEC_;

   std::vector<TFile*> inputFilesToDelete_;

   double bField_;
   int verbosity_;
   
   //HGCAL
   MagneticField const *aField_;
   hgcal::RecHitTools recHitTools_;
   std::vector<float> layerPositions_;

 
};

#endif
