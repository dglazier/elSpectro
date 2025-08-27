#pragma once

#include "Hadron.h"

namespace phoPro{

  class X3872Meson : public Hadron{

    static constexpr double X3872PDG = 204431;

  public:
    
  X3872Meson():Hadron(X3872PDG) {

      SetMass(3.872);
      SetWidth(0.001);
      
      RhoMeson rho; //Make sure registered
      JpsiMeson jpsi;
      
      AddDecay(0.05,{rho,jpsi});

      RegisterParticle(BreitWignerDistribution());
      
      
    }

    
  };

}
