#pragma once

#include "Hadron.h"

namespace phoPro{

  class RhoMeson : public Hadron{
    
    static constexpr int myPDG=113;
    
  public:
    
  RhoMeson():Hadron(myPDG,"Meson"){
     
      AddDecay(1,{Hadron(211),Hadron(-211)});//100% decay to 2 pions
      
      RegisterParticle(BreitWignerDistribution());

    }
    
  };

  class OmegaMeson : public Hadron{
    
    static constexpr int myPDG=223;
    
  public:
    
  OmegaMeson():Hadron(myPDG,"Meson"){
     
      AddDecay(1,{Hadron(111),Hadron(211),Hadron(-211)});//100% decay to 2 pions
      
      //RegisterParticle(BreitWignerDistribution());

    }
    
  };

  class PhiMeson : public Hadron{
    
    static constexpr int myPDG=333;
    
  public:
    
  PhiMeson():Hadron(myPDG,"Meson"){
     
      AddDecay(1,{Hadron(321),Hadron(-321)});//100% decay to 2 charged kaons
      
      RegisterParticle(BreitWignerDistribution());

    }
    
  };


  class JpsiMeson : public Hadron{

    static constexpr int myPDG=443;

  public:
  JpsiMeson():Hadron(myPDG,"Meson"){

      AddDecay(0.05971,{Hadron(11),Hadron(-11)});//6% to e+e-
      // SetWidth(0.0000926);
      
      RegisterParticle(BreitWignerDistribution()); 
    }
    
  };

  class UpsilonMeson : public Hadron{

    static constexpr int myPDG=553;

  public:
  UpsilonMeson():Hadron(myPDG,"Meson"){

      AddDecay(0.0238,{Hadron(11),Hadron(-11)});//0.0238
      //SetWidth(0.000054);
   
      RegisterParticle(BreitWignerDistribution());
    }
    
  };


  

}
