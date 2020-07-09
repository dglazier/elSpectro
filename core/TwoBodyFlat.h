//////////////////////////////////////////////////////////////
///
///Class:		TwoBodyFlat
///Description:
///            Derive methods for 
///            Generating LorentzVectors for parent->1,2
#pragma once

#include "DecayVectors.h"
#include <Math/VectorUtil.h> //for boosts etc.

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;

  class TwoBodyFlat : public DecayVectors {

 
  public:

    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;
    
  private:

    LorentzVector _a;
    LorentzVector _b;
    
    ClassDef(elSpectro::TwoBodyFlat,1); //class DecayVectors
 

  };

}
  
