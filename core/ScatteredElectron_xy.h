//////////////////////////////////////////////////////////////
///
///Class:		ScatteredElectron_xy
///Description:
///            Generating LorentzVectors for electron and virtual photon
///            From sampled scattered electron x and y
#pragma once

#include "DecayVectors.h"
#include "Distribution.h"
#include <Math/VectorUtil.h> //for boosts etc.
#include <Math/RotationY.h>

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;

  class ScatteredElectron_xy : public DecayVectors {

 
  public:
    
    ScatteredElectron_xy(Distribution* dist);
    
    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;

    Distribution* Dist()const {return _random_xy.get();}

  protected:

    void RotateZaxisToCMDirection(const LorentzVector& parent);
    
  private:

    LorentzVector _scattered;
    LorentzVector _gamma_ion; //residual gamma* + ion system

    dist_uptr _random_xy;
    
    ROOT::Math::RotationY _rotateToZaxis;
    LorentzVector _cachedParent;
    
    ClassDef(elSpectro::ScatteredElectron_xy,1); //class DecayVectors
 

  };

  inline void ScatteredElectron_xy::RotateZaxisToCMDirection(const LorentzVector& parent){
    if(_cachedParent!=parent){ //SetAngle is expensive (sin,cos calls) only call if necessary
      _cachedParent = parent;
      _rotateToZaxis.SetAngle(_cachedParent.Theta());
     }
    _scattered=_rotateToZaxis*_scattered;
  }
  
}
  
