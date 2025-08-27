///////////////////////////////////////////////////////////
// taken from $JPACPHOTO/scripts/exclusive_XYZ/Z_mesons.cpp

#pragma once

#include "JpacAmplitude.h" //interface

#include "constants.hpp"
#include "kinematics.hpp"
#include "blended.hpp"

#include "analytic/pomeron_exchange.hpp"
#include "covariant/pomeron_exchange.hpp"

namespace jpacAmps {

  // using namespace jpacPhoto;

  //J^PC = 1^+ 1^--
  class Y : public JpacAmplitude {

  public:

    Y(double mass): JpacAmplitude(mass) {
      Kine()->set_meson_JP(VECTOR);
    }
      
  protected:
   // Pomeron trajectory (high energies)
    double alpha_0H = 1.15;
    double alpha_PH = 0.11;

    // Pomeron trajectory (low energies)
    double alpha_0L = 0.94;
    double alpha_PL = 0.36;

    // Slope parameter and normalization (high)
    double bH = 1.01;
    double AH = E * 0.16;

    // Slope parameter and normalization (low)
    double bL = 0.12;
    double AL = E * 0.38;

    // Scale parameters
    double R_Jpsi  = 1;
    double R_Psi2S = 0.55;
    double R_Y     = 0.84;

    
    double low_sLimit = 7.*7;
    double high_sLimit = 20.*20;
 

  };
  class Y4260 : public  Y {
    
    
  public :
    
    Y4260() : Y(M_Y4260) {};
    
    //interface getter
    amplitude Get() override {return Amp_Blend();}
   
    amplitude Amp_High(){
      // Y(4260)
      amplitude   YH = new_amplitude<analytic::pomeron_exchange>(Kine(), "#it{Y}(4260)");
      YH->set_parameters({AH*R_Y, bH, alpha_0H, alpha_PH});

      return YH;
    }
  
    amplitude Amp_Low(){
      // Y(4260)
      amplitude   YL = new_amplitude<covariant::pomeron_exchange>(Kine(), "#it{Y}(4260)");
      YL->set_parameters({AL*R_Y, bL, alpha_0L, alpha_PL});
      
      return YL;
    }
    
    amplitude Amp_Blend(){
      return new_amplitude<blended>(Kine(),std::array<amplitude,2>{Amp_Low(),Amp_High()} ,std::array<double,2>{low_sLimit,high_sLimit} , "Zc blended");
       
    }
    
  };

}
