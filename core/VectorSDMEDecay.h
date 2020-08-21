//////////////////////////////////////////////////////////////
///
///Class:		VectorSDMEDecay
///Description:
///             Calculate intensity based on vector SDME
///             See Schilling and Wolf formalism

#pragma once

#include "DecayModel.h"
#include "DecayingParticle.h"
#include "SDME.h"
#include "PhotonPolarisationVector.h"

namespace elSpectro{

 
  class VectorSDMEDecay : public DecayModel {

  public:
    
    VectorSDMEDecay()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    VectorSDMEDecay( particle_ptrs , const std::vector<int> pdgs );

    // Each model must define its intensity
    double Intensity() const final;
    
    //Keep trying with new DecayVectors until pass
    bool RegenerateOnFail() const  noexcept final {return false;}

    void PostInit(ReactionInfo* info) final;
    
    bool CanUseSDME()const noexcept final{return true;}

  private:
 
    const SDME* _rho ={nullptr};

    LorentzVector* _photon={nullptr};
    LorentzVector* _meson={nullptr};
    LorentzVector* _baryon={nullptr};
    LorentzVector* _child1={nullptr};
    PhotonPolarisationVector* _photonPol={nullptr};
    
    
    ClassDefOverride(elSpectro::VectorSDMEDecay,1); //class VectorSDMEDecay
    
  };//class VectorSDMEDecay

}//namespace elSpectro
