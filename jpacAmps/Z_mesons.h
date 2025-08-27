///////////////////////////////////////////////////////////
// taken from $JPACPHOTO/scripts/exclusive_XYZ/Z_mesons.cpp

#pragma once

#include "JpacAmplitude.h" //interface

#include "constants.hpp"
#include "kinematics.hpp"
#include "blended.hpp"

#include "analytic/pseudoscalar_exchange.hpp"
#include "covariant/pseudoscalar_exchange.hpp"
#include "regge/pseudoscalar_exchange.hpp"

namespace jpacAmps {

  // using namespace jpacPhoto;

  //J^PC = 1^+ 1^+-
  class Zc : public JpacAmplitude {

  public:

    Zc(double mass): JpacAmplitude(mass) {
      Kine()->set_meson_JP(AXIALVECTOR);
    }
      
  protected:

    // Bottom vertex coupling (pi - nucleon - nucleon)
    double g_piNN = sqrt(2) * sqrt(4*PI*13.81); 

    // Cutoff for exponential form factor
    double lambda_pi = .900;  // MeV 

    // Pion trajectory parameters
    double slope = 0.7; // GeV-2
    double inter = - M2_PION * slope;

    // Zc(3900) couplings 
    double gc_jpsi  = 1.91; // psi coupling before VMD scaling
    double gc_gamma = E * F_JPSI * gc_jpsi / M_JPSI;

 

    double low_sLimit = 15.*15;
    double high_sLimit = 20.*20;
    
  };
  
  class Zc3900 : public  Zc {
    
    
  public :
    
    Zc3900() : Zc(M_ZC3900) {};
    
    //interface getter
    amplitude Get() override {return Amp_Blend();}
   
    amplitude Amp_High(){
      // Zc(3900)
      amplitude   ZcR = new_amplitude<regge::pseudoscalar_exchange>(Kine(), "#it{Z}_{c}(3900)");
      ZcR->set_parameters({gc_gamma, g_piNN, lambda_pi, inter, slope});
      return ZcR;
    }
  
    amplitude Amp_Low(){
      // Zc(3900)
      amplitude   Zc = new_amplitude<analytic::pseudoscalar_exchange>(Kine(),   M_PION, "#it{Z}_{c}(3900)");
      Zc->set_parameters({gc_gamma, g_piNN, lambda_pi});
      return Zc;

    }
    
    amplitude Amp_Blend(){
      return new_amplitude<blended>(Kine(),std::array<amplitude,2>{Amp_Low(),Amp_High()} ,std::array<double,2>{low_sLimit,high_sLimit} , "Zc blended");
       
    }
    
  };
  /////////////////////////////////////////////////////////////////////
 class Zb : public JpacAmplitude {

  public:

    Zb(double mass): JpacAmplitude(mass) {
      Kine()->set_meson_JP(AXIALVECTOR);
    }
      
  protected:

    // Bottom vertex coupling (pi - nucleon - nucleon)
    double g_piNN = sqrt(2) * sqrt(4*PI*13.81); 

    // Cutoff for exponential form factor
    double lambda_pi = .900;  // MeV 

    // Pion trajectory parameters
    double slope = 0.7; // GeV-2
    double inter = - M2_PION * slope;

  
    // Zb(10610) couplings
    double gb_upsilon1S = 0.49, gb_upsilon2S = 3.30, gb_upsilon3S = 9.22;
    double gb_gamma = E * (F_UPSILON1S * gb_upsilon1S / M_UPSILON1S 
                         + F_UPSILON2S * gb_upsilon2S / M_UPSILON2S
                         + F_UPSILON3S * gb_upsilon3S / M_UPSILON3S);  

    // Zb'(10650) couplings
    double gbp_upsilon1S = 0.21, gbp_upsilon2S = 1.47, gbp_upsilon3S = 4.8;
    double gbp_gamma = E * (F_UPSILON1S * gbp_upsilon1S / M_UPSILON1S 
                         +  F_UPSILON2S * gbp_upsilon2S / M_UPSILON2S
                         +  F_UPSILON3S * gbp_upsilon3S / M_UPSILON3S);  


    double low_sLimit = 20.*20; //10GeV above threshold...
    double high_sLimit = 22.*22;
    
  };

 //Zb(10610)
  class Zb10610 : public  Zb {
    
    

  public :
    
    Zb10610() : Zb(M_ZB10610) {};
    
    //interface getter
    amplitude Get() override {return Amp_Blend();}
   
    amplitude Amp_High(){
      // Zb(10610)
      amplitude   ZbR = new_amplitude<regge::pseudoscalar_exchange>(Kine(), "#it{Z}_{b}(10610)");
      ZbR->set_parameters({gb_gamma, g_piNN, lambda_pi, inter, slope});

      return ZbR;
    }
  
    amplitude Amp_Low(){
      // Zb(10610)
      amplitude   Zb = new_amplitude<analytic::pseudoscalar_exchange>(Kine(),   M_PION, "#it{Z}_{b}(10610)");
      Zb->set_parameters({gb_gamma, g_piNN, lambda_pi});

      return Zb;
    }
    
    amplitude Amp_Blend(){
      return new_amplitude<blended>(Kine(),std::array<amplitude,2>{Amp_Low(),Amp_High()} ,std::array<double,2>{low_sLimit,high_sLimit} , "Zb blended");
       
    }
    
  };
  
 //Zb(10650)
  class Zb10650 : public  Zb {
    
    

  public :
    
    Zb10650() : Zb(M_ZB10650) {};
    
    //interface getter
    amplitude Get() override {return Amp_Blend();}
   
    amplitude Amp_High(){
      // Zb'(10650)
      amplitude   ZbpR = new_amplitude<regge::pseudoscalar_exchange>(Kine(), "#it{Z}_{b}'(10650)");
      ZbpR->set_parameters({gbp_gamma, g_piNN, lambda_pi, inter, slope});
      return   ZbpR;
    }
  
    amplitude Amp_Low(){
      // Zb'(10650)
      amplitude   Zbp = new_amplitude<analytic::pseudoscalar_exchange>(Kine(), M_PION, "#it{Z}_{b}'(10650)");
      Zbp->set_parameters({gbp_gamma, g_piNN, lambda_pi});

      return Zbp;
    }
    
    amplitude Amp_Blend(){
      return new_amplitude<blended>(Kine(),std::array<amplitude,2>{Amp_Low(),Amp_High()} ,std::array<double,2>{low_sLimit,high_sLimit} , "Zb blended");
       
    }
    
  };
}
