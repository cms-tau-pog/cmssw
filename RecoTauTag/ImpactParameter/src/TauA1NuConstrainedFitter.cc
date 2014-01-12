/* From SimpleFits Package
 * Designed an written by
 * author: Ian M. Nugent
 * Humboldt Foundations
 */
#include "RecoTauTag/ImpactParameter/interface/TauA1NuConstrainedFitter.h"
#include "RecoTauTag/ImpactParameter/interface/PDGInfo.h"
#include "RecoTauTag/ImpactParameter/interface/LagrangeMultipliersFitter.h"
#include <iostream>

unsigned int TauA1NuConstrainedFitter::static_amb;

TauA1NuConstrainedFitter::TauA1NuConstrainedFitter(unsigned int ambiguity,LorentzVectorParticle A1,TVector3 PVertex, TMatrixTSym<double> VertexCov):
  MultiProngTauSolver(),
  ambiguity_(ambiguity)
{
  TLorentzVector Tau(0,0,0,0);
  //dummy substitution not used later
  LorentzVectorParticle Nu(TMatrixT<double>(LorentzVectorParticle::NLorentzandVertexPar,1),TMatrixTSym<double>(LorentzVectorParticle::NLorentzandVertexPar),PDGInfo::nu_tau,0.0,A1.BField()); 
  particles_.push_back(A1);
  particles_.push_back(Nu);
  
  // setup 13 by 13 matrix
  int size=LorentzVectorParticle::NVertex+particles_.size()*LorentzVectorParticle::NLorentzandVertexPar;
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);

  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;
  inpar(LorentzVectorParticle::vx,0)=PVertex.X();
  inpar(LorentzVectorParticle::vy,0)=PVertex.Y();
  inpar(LorentzVectorParticle::vz,0)=PVertex.Z();  
  for(int i=0; i<LorentzVectorParticle::NVertex;i++){
    for(int j=0; j<LorentzVectorParticle::NVertex;j++)incov(i,j)=VertexCov(i,j);
  }
  int A1offset=LorentzVectorParticle::NVertex;
  int Nuoffset=LorentzVectorParticle::NLorentzandVertexPar+LorentzVectorParticle::NVertex;
  for(int i=0; i<LorentzVectorParticle::NLorentzandVertexPar;i++){
    inpar(i+A1offset,0)=A1.Parameter(i);
    inpar(i+Nuoffset,0)=Nu.Parameter(i)+1.0;// offset by 1 GeV to prevent convergence on first iteration
    for(int j=0; j<LorentzVectorParticle::NLorentzandVertexPar;j++){
      incov(i+A1offset,j+A1offset)=A1.Covariance(i,j);
      incov(i+Nuoffset,j+Nuoffset)=Nu.Covariance(i,j);
    }
  }
  exppar.ResizeTo(nexpandedpar,1);
  exppar=ComputeInitalExpPar(inpar);
  expcov.ResizeTo(nexpandedpar,nexpandedpar);
  expcov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeInitalExpPar,inpar,incov);
  TMatrixT<double> PAR_0(npar,1);
  par_0.ResizeTo(npar);
  cov_0.ResizeTo(npar,npar);
  PAR_0=ComputeExpParToPar(exppar);
  for(int i=0; i<npar;i++)par_0(i)=PAR_0(i,0);
  cov_0=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeExpParToPar,exppar,expcov);

  for(int i=0; i<npar;i++){
    for(int j=0;j<npar;j++){cov_0(i,j)=expcov(i,j);}
  }
  par.ResizeTo(npar);
  par=par_0;
  cov.ResizeTo(npar,npar);
  cov=cov_0;
}

TMatrixT<double> TauA1NuConstrainedFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(nexpandedpar,1);
  int offset=LorentzVectorParticle::NVertex;// for A1
  TVector3 pv(inpar(LorentzVectorParticle::vx,0),inpar(LorentzVectorParticle::vy,0),inpar(LorentzVectorParticle::vz,0));
  TVector3 sv(inpar(LorentzVectorParticle::vx+offset,0),inpar(LorentzVectorParticle::vy+offset,0),inpar(LorentzVectorParticle::vz+offset,0));
  TVector3 TauDir=sv-pv;
  outpar(tau_phi,0)=TauDir.Phi();
  outpar(tau_theta,0)=TauDir.Theta();
  outpar(a1_px,0)=inpar(LorentzVectorParticle::px+offset,0);
  outpar(a1_py,0)=inpar(LorentzVectorParticle::py+offset,0);
  outpar(a1_pz,0)=inpar(LorentzVectorParticle::pz+offset,0);
  outpar(a1_m,0)=inpar(LorentzVectorParticle::m+offset,0);
  outpar(a1_vx,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(a1_vy,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(a1_vz,0)=inpar(LorentzVectorParticle::vz+offset,0);
  offset+=LorentzVectorParticle::NLorentzandVertexPar; // for Nu
  outpar(nu_px,0)=inpar(LorentzVectorParticle::px+offset,0);
  outpar(nu_py,0)=inpar(LorentzVectorParticle::py+offset,0);
  outpar(nu_pz,0)=inpar(LorentzVectorParticle::pz+offset,0);
  return outpar;
}

TMatrixT<double> TauA1NuConstrainedFitter::ComputeExpParToPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npar,1);
  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}

TMatrixT<double> TauA1NuConstrainedFitter::ComputeNuLorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(LorentzVectorParticle::NLorentzandVertexPar,1);
  outpar(LorentzVectorParticle::vx,0)=inpar(a1_vx,0);
  outpar(LorentzVectorParticle::vy,0)=inpar(a1_vy,0);
  outpar(LorentzVectorParticle::vz,0)=inpar(a1_vz,0);
  outpar(LorentzVectorParticle::px,0)=inpar(nu_px,0);
  outpar(LorentzVectorParticle::py,0)=inpar(nu_py,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(nu_pz,0);
  outpar(LorentzVectorParticle::m,0)=0;
  return outpar;
}

TMatrixT<double> TauA1NuConstrainedFitter::ComputeA1LorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(LorentzVectorParticle::NLorentzandVertexPar,1);
  outpar(LorentzVectorParticle::vx,0)=inpar(a1_vx,0);
  outpar(LorentzVectorParticle::vy,0)=inpar(a1_vy,0);
  outpar(LorentzVectorParticle::vz,0)=inpar(a1_vz,0);
  outpar(LorentzVectorParticle::px,0)=inpar(a1_px,0);
  outpar(LorentzVectorParticle::py,0)=inpar(a1_py,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(a1_pz,0);
  outpar(LorentzVectorParticle::m,0)=inpar(a1_m,0);
  return outpar;
}

TMatrixT<double> TauA1NuConstrainedFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixT<double> nupar=ComputeNuLorentzVectorPar(inpar);
  TMatrixT<double> a1par=ComputeA1LorentzVectorPar(inpar);
  for(int i=0;i<LorentzVectorParticle::NLorentzandVertexPar;i++){
    if(i<LorentzVectorParticle::NVertex){outpar(i,0)=a1par(i,0);}
    else{outpar(i,0)=nupar(i,0)+a1par(i,0);}
    //if(i==LorentzVectorParticle::m) outpar(i,0)=PDGInfo::tau_mass();
  }
  double Enu2=pow(nupar(LorentzVectorParticle::px,0),2.0)+pow(nupar(LorentzVectorParticle::py,0),2.0)+pow(nupar(LorentzVectorParticle::pz,0),2.0);
  double Ea12=pow(a1par(LorentzVectorParticle::px,0),2.0)+pow(a1par(LorentzVectorParticle::py,0),2.0)+pow(a1par(LorentzVectorParticle::pz,0),2.0)+pow(a1par(LorentzVectorParticle::m,0),2.0);
  double P2=pow(outpar(LorentzVectorParticle::px,0),2.0)+pow(outpar(LorentzVectorParticle::py,0),2.0)+pow(outpar(LorentzVectorParticle::pz,0),2.0);
  outpar(LorentzVectorParticle::m,0)=sqrt(fabs(pow(sqrt(Enu2)+sqrt(Ea12),2.0)-P2));
  return outpar;
}

void TauA1NuConstrainedFitter::UpdateExpandedPar(){
  // assumes changes to a1 correlation to vertex is small
  if(par.GetNrows()==npar && cov.GetNrows()==npar && exppar.GetNrows()==npar && expcov.GetNrows()==npar) return;
  for(int i=0; i<npar;i++){
    exppar(i,0)=par(i);
    for(int j=0; j<npar;j++){expcov(i,j)=cov(i,j);}
  }
}

std::vector<LorentzVectorParticle> TauA1NuConstrainedFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> a1=ComputeA1LorentzVectorPar(exppar);
  TMatrixTSym<double> a1cov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeA1LorentzVectorPar,exppar,expcov);
  refitParticles.push_back(LorentzVectorParticle(a1,a1cov,particles_.at(0).PDGID(),c,b));
  TMatrixT<double> nu=ComputeNuLorentzVectorPar(exppar);
  TMatrixTSym<double> nucov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeNuLorentzVectorPar,exppar,expcov);
  refitParticles.push_back(LorentzVectorParticle(nu,nucov,PDGInfo::nu_tau,0.0,b));
  return refitParticles;
}

LorentzVectorParticle TauA1NuConstrainedFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeMotherLorentzVectorPar,exppar,expcov);
  LorentzVectorParticle mymother= LorentzVectorParticle(m,mcov,(int)(-1.0*fabs(PDGInfo::tau_minus)*c),c,b);
  return mymother;
}

void TauA1NuConstrainedFitter::CovertParToObjects(TVectorD &v,TLorentzVector &a1,TLorentzVector &nu,double &phi,double &theta,TVector3 &TauDir){
  a1=TLorentzVector(v(a1_px),v(a1_py),v(a1_pz),sqrt(v(a1_m)*v(a1_m)+v(a1_px)*v(a1_px)+v(a1_py)*v(a1_py)+v(a1_pz)*v(a1_pz)));
  nu=TLorentzVector(v(nu_px),v(nu_py),v(nu_pz),sqrt(v(nu_px)*v(nu_px)+v(nu_py)*v(nu_py)+v(nu_pz)*v(nu_pz)));
  phi=v(tau_phi);
  theta=v(tau_theta);
  TauDir.SetMagThetaPhi(1.0,theta,phi);
}

bool TauA1NuConstrainedFitter::Fit(){
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Check if Tau Direction is unphysical and if nessicary set the starting point to Theta_{GJ-Max}
  TLorentzVector a1(par(a1_px),par(a1_py),par(a1_pz),sqrt(par(a1_m)*par(a1_m)+par(a1_px)*par(a1_px)+par(a1_py)*par(a1_py)+par(a1_pz)*par(a1_pz)));
  double phi(par(tau_phi)),theta(par(tau_theta));
  TLorentzVector Tau_plus,Tau_minus,nu_plus,nu_minus;
  TVector3 TauDir(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
  bool isReal;
  SolvebyRotation(TauDir,a1,Tau_plus,Tau_minus,nu_plus,nu_minus,isReal);
  TMatrixT<double>    thepar=LagrangeMultipliersFitter::convertToMatrix(par);
  static_amb=ambiguity_;

  //case 1: is real then solve analytically
  if(isReal && (ambiguity_==plus || ambiguity_==minus)){
    // popogate errors
    TMatrixT<double> par_tmp=TauA1NuConstrainedFitter::SolveAmbiguityAnalytically(thepar);
    cov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::SolveAmbiguityAnalytically,thepar,cov_0);
    for(int i=0; i<npar;i++) par(i)=par_tmp(i,0);
    return true;
  }
  // case 2 is in unphsyical region - rotate and substitue \theta_{GJ} with \theta_{GJ}^{Max} and then solve analytically
  else if(ambiguity_==zero){
    TMatrixT<double> par_tmp=TauA1NuConstrainedFitter::SolveAmbiguityAnalyticallywithRot(thepar);
    cov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::SolveAmbiguityAnalyticallywithRot,thepar,cov_0);
    for(int i=0; i<npar;i++) par(i)=par_tmp(i,0);
    return true;
  }
  return false;
}

TMatrixT<double> TauA1NuConstrainedFitter::SolveAmbiguityAnalytically(TMatrixT<double> &inpar){
  // Solve equation quadratic equation
  TMatrixT<double> outpar(inpar.GetNrows(),1);
  TLorentzVector a1,nu;
  double phi(0),theta(0);
  TVector3 TauDir;
  TVectorT<double> v=LagrangeMultipliersFitter::convertToVector(inpar);
 CovertParToObjects(v,a1,nu,phi,theta,TauDir);
  TLorentzVector a1_d=a1;
  TLorentzVector nu_d=nu;
  TLorentzVector Tau_plus,Tau_minus,nu_plus,nu_minus;
  bool isReal;
  SolvebyRotation(TauDir,a1_d,Tau_plus,Tau_minus,nu_plus,nu_minus,isReal,true);
  if(static_amb==plus)nu=nu_plus;
  else nu=nu_minus;
  for(int i=0; i<outpar.GetNrows();i++){ outpar(i,0)=v(i);}
  outpar(nu_px,0)=nu.Px();  
  outpar(nu_py,0)=nu.Py();
  outpar(nu_pz,0)=nu.Pz(); 
  return outpar;
}

TMatrixT<double> TauA1NuConstrainedFitter::SolveAmbiguityAnalyticallywithRot(TMatrixT<double> &inpar){
  // Rotate and subsitute \theta_{GJ} with \theta_{GJ}^{Max} - assumes uncertianty on thata and phi of the a1 or small compared to the tau direction. 
  TMatrixT<double> outpar(inpar.GetNrows(),1);
  TVectorT<double> v=LagrangeMultipliersFitter::convertToVector(inpar);
  TLorentzVector a1,nu;
  double phi(0),theta(0);
  TVector3 TauDir;
  CovertParToObjects(v,a1,nu,phi,theta,TauDir);
  double theta_a1(a1.Theta()),phi_a1(a1.Phi()),theta_GJMax(ThetaGJMax(a1));
  TauDir.RotateZ(-phi_a1);
  TauDir.RotateY(-theta_a1);
  double phiprime(TauDir.Phi());
  TauDir=TVector3(sin(theta_GJMax)*cos(phiprime),sin(theta_GJMax)*sin(phiprime),cos(theta_GJMax));
  TauDir.RotateY(theta_a1);
  TauDir.RotateZ(phi_a1);
  for(int i=0; i<outpar.GetNrows();i++) outpar(i,0)=v(i);
  outpar(tau_phi,0)=TauDir.Phi();
  outpar(tau_theta,0)=TauDir.Theta();
  return SolveAmbiguityAnalytically(outpar);
}

