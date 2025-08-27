#pragma once

#include "Hadron.h"

namespace phoPro{

  class Y4260Meson : public Hadron{

    static constexpr double myPDG = 4431; //Jpsi with additional 1

  public:
    
  Y4260Meson():Hadron(myPDG) {

      SetMass(4.220);
      SetWidth(0.05);
      
      JpsiMeson jpsi;
      
      AddDecay(0.01,{jpsi,Hadron(211),Hadron(-211)}); //Jpsi,pi+,pi-

      RegisterParticle(BreitWignerDistribution());
      
      
    }

    
  };
}
