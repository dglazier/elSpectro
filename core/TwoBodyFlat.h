//////////////////////////////////////////////////////////////
///
///Class:		TwoBodyFlat
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

    virtual double MyRandomCosTh() const noexcept{ return gRandom->Uniform(-1,1); }
    virtual double RandomCosTh() const noexcept{
      if(_forIntegral==true)
	//integration requires flat distribution
	return gRandom->Uniform(-1,1);
      return MyRandomCosTh();
    }
 

  protected :

    double W() const noexcept{return _W;}
    // void RotateZaxisToCMDirection(const LorentzVector& parent);

  private:

    LorentzVector _a;
    LorentzVector _b;
    double _W={0};
 
 
    ClassDef(elSpectro::TwoBodyFlat,1); //class DecayVectors
 

  };
 
  
}
  
