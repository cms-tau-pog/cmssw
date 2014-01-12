#ifndef MultiProngTauSolver_H
#define MultiProngTauSolver_H
/* From SimpleFits Package
 * Designed an written by
 * author: Ian M. Nugent
 * Humboldt Foundations
 */
#include "RecoTauTag/ImpactParameter/interface/MultiProngTauSolver.h"
#include "RecoTauTag/ImpactParameter/interface/LorentzVectorParticle.h"
#include "RecoTauTag/ImpactParameter/interface/ErrorMatrixPropagator.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class  MultiProngTauSolver {
 public:
  enum Ambiguity{zero,minus,plus,NAmbiguity};

  // constructor and Destructor
  MultiProngTauSolver(){};
  virtual ~MultiProngTauSolver(){};
      
  static void quadratic(double &x_plus,double &x_minus,double a, double b, double c, bool &isReal);
  static void AnalyticESolver(TLorentzVector &nu_plus,TLorentzVector &nu_minus,TLorentzVector A1,bool &isReal);
  static void NumericalESolver(TLorentzVector &nu_plus,TLorentzVector &nu_minus,TLorentzVector A1);
  static void SolvebyRotation(TVector3 TauDir,TLorentzVector A1, TLorentzVector &Tau_plus,TLorentzVector &Tau_minus,TLorentzVector &nu_plus,TLorentzVector &nu_minus,bool &isReal,bool rotateback=true);
  static bool SetTauDirectionatThetaGJMax(TLorentzVector a1, double &theta,double &phi,double scale=1.0);  
  static double ThetaGJMax(TLorentzVector a1);
  static LorentzVectorParticle EstimateNu(LorentzVectorParticle &a1,TVector3 pv,int ambiguity,TLorentzVector &tau);

  static TMatrixT<double> RotateToTauFrame(TMatrixT<double> &inpar);
};

#endif
