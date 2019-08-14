#ifndef DataFormats_lowPtTau_h
#define DataFormats_lowPtTau_h

#include <vector>

// a simple class
struct lowPtTau
{
  //  explicit lowPtTau(int v):value_(v) { }

  //  lowPtTau():value_(0) { }

  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t mass;
  Int_t charge;
  Float_t vprob;
  Float_t vz;
  Float_t dz2pv;
  Float_t dxy2pv;
  Float_t flightSig3D;
  Float_t flightLength3D;
  Float_t flightLengthErr3D;
  Float_t flightSig2D;
  Float_t flightLength2D;
  Float_t flightLengthErr2D;
  Int_t pfid1;
  Int_t pfid2;
  Int_t pfid3;
  Float_t pt1;
  Float_t pt2;
  Float_t pt3;
  Float_t eta1;
  Float_t eta2;
  Float_t eta3;
  Float_t phi1;
  Float_t phi2;
  Float_t phi3;
  Float_t mass1;
  Float_t mass2;
  Float_t mass3;
  Float_t isopt03;
  Float_t isopt04;
  Float_t max_dr_3prong;	
  Bool_t gen_isRight;
  Float_t gen_pt;
  Float_t gen_eta;
  Float_t gen_phi;
  Int_t gen_parent_pdgId;

  bool operator<(const lowPtTau& another) const { 
    return pt > another.pt;
  }

  //  int value_;
};

// this is our new product, it is simply a 
// collection of lowPtTau held in an std::vector
typedef std::vector<lowPtTau> lowPtTauCollection;

#endif
