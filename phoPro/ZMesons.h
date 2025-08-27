#pragma once

#include "Hadron.h"

namespace phoPro{

  class Zc3900Meson : public Hadron{

    static constexpr double myPDG = 104431; //hc(1P) with additional 1

  public:
    
  Zc3900Meson():Hadron(myPDG,"Meson") {

      SetMass(3.8884);
      SetWidth(0.05);
      
      JpsiMeson jpsi;
      
      AddDecay(0.1,{jpsi,Hadron{211}}); //10% to Jpsi,pi+

      RegisterParticle(BreitWignerDistribution());
      
      
    }

    
  };

  class Zb10610Meson : public Hadron{

    static constexpr double myPDG = 105531; //hb(1P) with additional 1

  public:
    
    Zb10610Meson():Hadron(myPDG,"Meson") {

      SetMass(10.6072);
      SetWidth(0.0184);
      
      UpsilonMeson ups;
      
      AddDecay(0.036,{ups,Hadron{211}}); //Upsilon,pi+, 0.036
      //Note BB* is 86% BR

      RegisterParticle(BreitWignerDistribution());
      
      
    }

    
  };

  class Zb10650Meson : public Hadron{

    static constexpr double myPDG = 105532; //hb(1P) with additional 2

  public:
    
    Zb10650Meson():Hadron(myPDG,"Meson") {

      SetMass(10.652);
      SetWidth(0.0115);
      
      UpsilonMeson ups;
      
      AddDecay(0.036,{ups,Hadron{211}}); //Upsilon,pi+, 0.036
 
      RegisterParticle(BreitWignerDistribution());
      
      
    }

    
  };

}
