//////////////////////////////////////////////////////////////
///
///Class:		DecayVectors
///Description:
///            Abstract class, derive methods for 
///            Generating LorentzVectors for S->1,2,3...
#pragma once

#include "Particle.h"
#include "DecayModel.h"
#include <vector>

namespace elSpectro{

 
  class DecayVectors {

  public:

    DecayVectors()=default;
    virtual ~DecayVectors()=default;
    DecayVectors(const DecayVectors& other); //need the virtual destructor...so rule of 5
    DecayVectors(DecayVectors&&)=default;
    DecayVectors& operator=(const DecayVectors& other);
    DecayVectors& operator=(DecayVectors&& other) = default;

    ///Function will decay parent to products and return a weight that
    ///can be used to get phase space distribution
    virtual double  Generate(const LorentzVector& parent,
			     const particle_ptrs& products)  = 0;

    
    virtual void PostInit(ReactionInfo* info){};

  protected:
    
    mutable double _weight={1};
 
  private:

    
    ClassDef(elSpectro::DecayVectors,1); //class DecayVectors
 

  };

}
  
