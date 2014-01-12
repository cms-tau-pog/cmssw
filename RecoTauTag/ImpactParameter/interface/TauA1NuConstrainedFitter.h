#ifndef TauA1NuConstrainedFitter_H
#define TauA1NuConstrainedFitter_H

/* From SimpleFits Package
 * Designed an written by
 * author: Ian M. Nugent
 * Humboldt Foundations
 */
#include "RecoTauTag/ImpactParameter/interface/MultiProngTauSolver.h"
#include "RecoTauTag/ImpactParameter/interface/LorentzVectorParticle.h"
#include "RecoTauTag/ImpactParameter/interface/ErrorMatrixPropagator.h"
#include <vector>
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"

class TauA1NuConstrainedFitter : public MultiProngTauSolver{
 public:
  TauA1NuConstrainedFitter(unsigned int ambiguity,LorentzVectorParticle A1,TVector3 PVertex, TMatrixTSym<double> VertexCov);
  virtual ~TauA1NuConstrainedFitter(){};

  enum Pars{tau_phi=0,tau_theta,a1_px,a1_py,a1_pz,a1_m,nu_px,nu_py,nu_pz,npar};
  enum ExpandedPars{a1_vx=9,a1_vy,a1_vz,nexpandedpar};
  enum OrignialPars{norigpar=13};

  std::vector<LorentzVectorParticle> GetReFitDaughters();
  LorentzVectorParticle GetMother();

  bool Fit();
    
 private:
  static TMatrixT<double> ComputeInitalExpPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeNuLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeA1LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> SolveAmbiguityAnalytically(TMatrixT<double> &inpar);
  static TMatrixT<double> SolveAmbiguityAnalyticallywithRot(TMatrixT<double> &inpar);

  void UpdateExpandedPar();
  static void CovertParToObjects(TVectorD &v,TLorentzVector &a1,TLorentzVector &nu,double &phi,double &theta,TVector3 &TauDir);

  TVectorD par_0;
  TVectorD par;
  TMatrixTSym<double> cov_0;
  TMatrixTSym<double> cov;

  TMatrixT<double> exppar;
  TMatrixTSym<double> expcov;
  std::vector<LorentzVectorParticle> particles_;
  unsigned int ambiguity_;

  static  unsigned int static_amb;
};
#endif
