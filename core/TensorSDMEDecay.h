//////////////////////////////////////////////////////////////
///
///Class:		TensorSDMEDecay
///Description:
///             Calculate intensity based on vector SDME
///             See Schilling and Wolf formalism

#pragma once

#include "SDMEDecay.h"

namespace elSpectro{

 
  class TensorSDMEDecay : public SDMEDecay {

  public:
    
    TensorSDMEDecay()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    TensorSDMEDecay( particle_ptrs , const std::vector<int> pdgs );

    // Each model must define its intensity
    double Intensity() const final;
    
    //void PostInit(ReactionInfo* info) final;
    
     short Spin() const final{ return 2;}

  private:
 
   ClassDefOverride(elSpectro::TensorSDMEDecay,1); //class TensorSDMEDecay
    
  };//class TensorSDMEDecay

}//namespace elSpectro
