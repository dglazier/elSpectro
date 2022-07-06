//////////////////////////////////////////////////////////////
///
///Class:		Bremsstrahlung
///Description:
///         Handle bremsstrahlung photon production, neglect brem e-

#pragma once

#include "DecayModel.h"

namespace elSpectro{

 
  class Bremsstrahlung : public DecayModel {

  public:
    
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    //Bremsstrahlung( particle_ptrs , const std::vector<int> pdgs );
    Bremsstrahlung( );

    // Each model must define its intensity
    // Phase space intensity is handled by MassPhaseSpace
    double Intensity() const final{
      return 1.;
    }
    bool RegenerateOnFail() const  noexcept final {return false;}
    ///void SetParent(DecayingParticle* pa);
    //void PostInit(ReactionInfo* info) override;
    
  private:
     ClassDefOverride(elSpectro::Bremsstrahlung,1); //class Bremsstrahlung
    
  };//class Bremsstrahlung



}//namespace elSpectro
