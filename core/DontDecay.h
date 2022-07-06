//////////////////////////////////////////////////////////////
///
///Class:		DontDecay
///Description:
///            Derive methods for 
///            Generating LorentzVectors for parent->1,2
#pragma once

#include "DecayVectors.h"
#include "LorentzVector.h"
#include <TRandom.h> //for gRandom
#include <Math/VectorUtil.h> //for boosts etc.

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;
  
  class DontDecay : public DecayVectors {
    
  public:

    DontDecay()=default;
    virtual ~DontDecay()=default;
    DontDecay(const DontDecay& other); //need the virtual destructor...so rule of 5
    DontDecay(DontDecay&&)=default;
    DontDecay& operator=(const DontDecay& other);
    DontDecay& operator=(DontDecay&& other) = default;
 

    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;

  
  protected :

    double W() const noexcept{return _W;}
  
  private:

    LorentzVector _a;
    double _W={0};
 
 
    ClassDef(elSpectro::DontDecay,1); //class DecayVectors
 

  };
 
  
}
  
