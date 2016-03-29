#ifndef DataFormats_PatCandidates_PATTauDiscriminator_h
#define DataFormats_PatCandidates_PATTauDiscriminator_h
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include <vector>

namespace pat {
  typedef edm::AssociationVector<pat::TauRefProd,std::vector<float> > PATTauDiscriminatorBase;
  
  class PATTauDiscriminator : public PATTauDiscriminatorBase {
  public:
    PATTauDiscriminator() :
      PATTauDiscriminatorBase()
      { }
    
    PATTauDiscriminator(const pat::TauRefProd & ref) :
      PATTauDiscriminatorBase(ref)
      { }
    
    PATTauDiscriminator(const PATTauDiscriminatorBase &v) :
      PATTauDiscriminatorBase(v)
      { }
  };
  
  typedef PATTauDiscriminator::value_type PATTauDiscriminatorVT;  
  typedef edm::Ref<PATTauDiscriminator> PATTauDiscriminatorRef;  
  typedef edm::RefProd<PATTauDiscriminator> PATTauDiscriminatorRefProd;  
  typedef edm::RefVector<PATTauDiscriminator> PATTauDiscriminatorRefVector; 
}
#endif
