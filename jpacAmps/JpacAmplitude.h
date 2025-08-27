#pragma once
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacAmps {
  using namespace jpacPhoto;


  class JpacAmplitude{

  public:
    
    JpacAmplitude(double mass){
      _kine  = new_kinematics( mass );
    }
    
    virtual amplitude Get() = 0;

    const kinematics& Kine() const {return _kine;}
    
  protected :
    kinematics _kine;

    
  };

}
