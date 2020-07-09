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
    ///Function will decay parent to products and return a weight that
    ///can be used to get phase space distribution
    virtual double  Generate(const LorentzVector& parent,
			     const particle_ptrs& products)  = 0;
    
  private:

    
    ClassDef(elSpectro::DecayVectors,1); //class DecayVectors
 

  };

}
  
