//////////////////////////////////////////////////////////////
///
///Class:		SDMEDecay
///Description:
///             Calculate intensity based on vector SDME
///             See Schilling and Wolf formalism

#pragma once

#include "DecayModel.h"
#include "DecayingParticle.h"
#include "SDME.h"
#include "PhotonPolarisationVector.h"

namespace elSpectro{

 
  class SDMEDecay : public DecayModel {

  public:
    
    SDMEDecay()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    SDMEDecay( particle_ptrs , const std::vector<int> pdgs );

    // Each model must define its intensity
    double Intensity() const override =0;
    
    //Keep trying with new DecayVectors until pass
    bool RegenerateOnFail() const  noexcept final {return true;}

    void PostInit(ReactionInfo* info) override;
    
    bool CanUseSDME()const noexcept final{return true;}
    virtual short Spin() const { return 0;}

    
  protected:
 
    const SDME* _rho ={nullptr};

    LorentzVector* _photon={nullptr};
    LorentzVector* _meson={nullptr};
    LorentzVector* _baryon={nullptr};
    LorentzVector* _child1={nullptr};
    PhotonPolarisationVector* _photonPol={nullptr};
    
   private:
  
    ClassDefOverride(elSpectro::SDMEDecay,1); //class SDMEDecay
    
  };//class SDMEDecay

}//namespace elSpectro
