//////////////////////////////////////////////////////////////
///
///Class:		TwoBodyFlat
///Description:
///            Derive methods for 
///            Generating LorentzVectors for parent->1,2
#pragma once

#include "DecayVectors.h"
#include <TRandom.h> //for gRandom
#include <Math/VectorUtil.h> //for boosts etc.
#include <Math/RotationY.h>

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;

  class TwoBodyFlat : public DecayVectors {
    
  public:

    TwoBodyFlat()=default;
    virtual ~TwoBodyFlat()=default;
    TwoBodyFlat(const TwoBodyFlat& other); //need the virtual destructor...so rule of 5
    TwoBodyFlat(TwoBodyFlat&&)=default;
    TwoBodyFlat& operator=(const TwoBodyFlat& other);
    TwoBodyFlat& operator=(TwoBodyFlat&& other) = default;
 

    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;

    virtual double RandomCosTh() const noexcept{ return gRandom->Uniform(-1,1); }
    virtual double RandomPhi() const noexcept { return gRandom->Uniform(-TMath::Pi(),TMath::Pi()); }

  protected :

    double W() const noexcept{return _W;}
    void RotateZaxisToCMDirection(const LorentzVector& parent);

  private:

    LorentzVector _a;
    LorentzVector _b;
    double _W={0};
 
    ROOT::Math::RotationY _rotateToZaxis;
    LorentzVector _cachedParent;

    ClassDef(elSpectro::TwoBodyFlat,1); //class DecayVectors
 

  };
  inline void TwoBodyFlat::RotateZaxisToCMDirection(const LorentzVector& parent){
    if(_cachedParent!=parent){ //SetAngle is expensive (sin,cos calls) only call if necessary
      _cachedParent = parent;
      _rotateToZaxis.SetAngle(_cachedParent.Theta());
     }
    _a=_rotateToZaxis*_a;
    _b=_rotateToZaxis*_b;
  }
  
}
  
