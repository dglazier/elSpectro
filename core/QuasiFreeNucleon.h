//////////////////////////////////////////////////////////////
///
///Class:		QuasiFreeNucleon
///Description:
///            Derive methods for 
///            Generating LorentzVectors for Nucleus->nucleon + spectator 
#pragma once

#include "DecayVectors.h"
#include "Distribution.h"
#include "DistTF1.h"
#include "LorentzVector.h"
#include <TRandom.h> //for gRandom
#include <Math/VectorUtil.h> //for boosts etc.

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;
  
  class QuasiFreeNucleon : public DecayVectors {
    
  public:

    //default deuteron Fermi distribution, can overwrite when create object 
    QuasiFreeNucleon(Distribution* dist=new DistTF1{TF1("deuteronFermiDist","(x*(0.26**2-0.0456**2)/(x**2+0.0456**2)/(x**2+0.26**2))**2",0.,1)}):_fermiDist{dist}{};
    
     virtual ~QuasiFreeNucleon()=default;
    QuasiFreeNucleon(const QuasiFreeNucleon& other); //need the virtual destructor...so rule of 5
    QuasiFreeNucleon(QuasiFreeNucleon&&)=default;
    QuasiFreeNucleon& operator=(const QuasiFreeNucleon& other);
    QuasiFreeNucleon& operator=(QuasiFreeNucleon&& other) = default;
 

    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;

    double RandomCosTh() const noexcept{
      return gRandom->Uniform(-1,1);
    }


  protected :

 
  private:

    LorentzVector _nucleon;
    LorentzVector _spectator;

    std::unique_ptr<Distribution> _fermiDist;
    
    ClassDef(elSpectro::QuasiFreeNucleon,1); //class DecayVectors
 

  };
 
  
}
  
