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
#include <Math/RotationZYX.h>
#include <TH1.h>
#include <TH2.h>

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;

  class ScatteredElectron_xy : public DecayVectors {

 
  public:
    
    ScatteredElectron_xy(double eb, double mion, double Wmin);
    
    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;
    double GenerateGivenXandY(const LorentzVector& parent,
		    const particle_ptrs& products,double xx, double yy)  ;
    
    DistVirtPhotFlux_xy &Dist(){return _random_xy;}


    void SetModel(DecayModel* model){_gStarNmodel=model;}

    double dsigma() const override{return _random_xy.CurrentValue();} // dxdy=1

    double Probability() const final{return _random_xy.Probability();}

  protected:

    void RotateZaxisToCMDirection(const LorentzVector& parent, LorentzVector& child);
    double CompleteGivenXandY(const LorentzVector& parent, const particle_ptrs& products,double xx, double yy);

  private:

    DecayModel* _gStarNmodel{nullptr};
    
    DistVirtPhotFlux_xy _random_xy;
 
    LorentzVector _scattered;
    LorentzVector _gamma_ion; //residual gamma* + ion system
    LorentzVector _parent_in_elFrame;
    
    ROOT::Math::RotationZYX _sc_rotateToZaxis;
    LorentzVector _cachedParent;
    
    // TH1D *hist, *histW;
    TH1D histy={"ydist","ydist",1000,0,1};
    TH2D histyQ2={"yQ2dist","yQ2dist",200,0,20,100,0,1};
    TH2D histyx={"yxdist","yxdist",100,0,1,100,0,1};
    TH2D histyCosTh={"yCosThdist","yCosThdist",100,-1,1,100,0,1};
    TH1D histW={"genWdist","genWdist",100,5,30};
    
    ClassDefOverride(elSpectro::ScatteredElectron_xy,1); //class DecayVectors
 

  };

  inline void ScatteredElectron_xy::RotateZaxisToCMDirection(const LorentzVector& parent, LorentzVector& child){
    if(_cachedParent!=parent){ //SetAngle is expensive (sin,cos calls) only call if necessary
      _cachedParent = parent;
      //_rotateToZaxis.SetAngle(_cachedParent.Theta());
      _sc_rotateToZaxis.SetComponents(-_cachedParent.Phi(),-_cachedParent.Theta(),0);
     }
    child=_sc_rotateToZaxis*child;
  }
  
}
  
