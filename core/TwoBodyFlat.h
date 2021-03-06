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

    virtual double RandomCosTh() const noexcept{
      return gRandom->Uniform(-1,1);
    }
    
    double Probability() const{return 1./4/TMath::Pi();}


  protected :

    double W() const noexcept{return _W;}
  
  private:

    LorentzVector _a;
    LorentzVector _b;
    double _W={0};
 
 
    ClassDef(elSpectro::TwoBodyFlat,1); //class DecayVectors
 

  };
 
  
}
  
