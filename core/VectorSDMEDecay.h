//////////////////////////////////////////////////////////////
///
///Class:		VectorSDMEDecay
///Description:
///             Calculate intensity based on vector SDME
///             See Schilling and Wolf formalism

#pragma once

#include "SDMEDecay.h"

namespace elSpectro{

 
  class VectorSDMEDecay : public SDMEDecay {

  public:
    
    VectorSDMEDecay()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    VectorSDMEDecay( particle_ptrs , const std::vector<int> pdgs );

    // Each model must define its intensity
    double Intensity() const final;
    
    //void PostInit(ReactionInfo* info) final;
    
    short Spin() const final{ return 1;}

  private:
 
   ClassDefOverride(elSpectro::VectorSDMEDecay,1); //class VectorSDMEDecay
    
  };//class VectorSDMEDecay

}//namespace elSpectro
