//////////////////////////////////////////////////////////////
///
///Class:		ScatteredElectron_xy
///Description:
///            Generating LorentzVectors for electron and virtual photon
///            From sampled scattered electron x and y
#pragma once

#include "DecayVectors.h"
#include "DistVirtPhotFlux_xy.h"
#include <Math/VectorUtil.h> //for boosts etc.
#include <Math/RotationY.h>
#include <TH1.h>

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;

  class ScatteredElectron_xy : public DecayVectors {

 
  public:
    
    ScatteredElectron_xy(double eb, double mion, double Wmin);
    
    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;

    DistVirtPhotFlux_xy &Dist(){return _random_xy;}


    void SetModel(DecayModel* model){_gStarNmodel=model;}
   protected:

    void RotateZaxisToCMDirection(const LorentzVector& parent);
    
  private:

    DecayModel* _gStarNmodel{nullptr};
    
    DistVirtPhotFlux_xy _random_xy;
 
    LorentzVector _scattered;
    LorentzVector _gamma_ion; //residual gamma* + ion system

    
    ROOT::Math::RotationY _rotateToZaxis;
    LorentzVector _cachedParent;
    
    // TH1D *hist, *histW;
    TH1D histy={"ydist","ydist",1000,0,1};
    TH1D histW={"genWdist","genWdist",1000,0,100};
    
    ClassDefOverride(elSpectro::ScatteredElectron_xy,1); //class DecayVectors
 

  };

  inline void ScatteredElectron_xy::RotateZaxisToCMDirection(const LorentzVector& parent){
    if(_cachedParent!=parent){ //SetAngle is expensive (sin,cos calls) only call if necessary
      _cachedParent = parent;
      _rotateToZaxis.SetAngle(_cachedParent.Theta());
     }
    _scattered=_rotateToZaxis*_scattered;
  }
  
}
  
