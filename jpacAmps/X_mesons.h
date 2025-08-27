///////////////////////////////////////////////////////////
// taken from $JPACPHOTO/scripts/exclusive_XYZ/X_mesons.cpp

#pragma once

#include "JpacAmplitude.h" //interface

#include "constants.hpp"
#include "kinematics.hpp"
#include "blended.hpp"

#include "analytic/vector_exchange.hpp"
#include "covariant/vector_exchange.hpp"
#include "regge/vector_exchange.hpp"

//example usage :
//    jpacAmps::X3872 X3872amp;
//    auto jpac_amp = X3872amp.Get();
//    jpac_amp->probability_distribution(54,-4); //s,t

namespace jpacAmps {

 
  //J^PC = 1^++
  class Chic1 : public JpacAmplitude {

  public:

    Chic1(double mass): JpacAmplitude(mass) {
      Kine()->set_meson_JP(AXIALVECTOR);
    }
      
  protected:

    // Nucleon couplings 
    double gV_omega = 16., gT_omega = 0.;
    double gV_rho = 2.4,  gT_rho = 14.6;
    double gV_phi = -6.2, gT_phi = 2.1;
    double gV_psi = 1.6E-3, gT_psi = 0.;
    
    // Photon couplings
    double gChi_omega   = 5.2E-4;
    double gChi_rho     = 9.2E-4;
    double gChi_phi     = 4.2E-4;
    double gChi_psi     = 1.;
    
    // Form factor cutoffs
    double LamOmega = 1.2;
    double LamRho   = 1.4; 
    double gX_omega     = 8.2E-3;
    double gX_rho       = 3.6E-3;
 
    // Reggeon trajectory
    double inter = 0.5;
    double slope = 0.9;


    double low_sLimit = 7.*7;
    double high_sLimit = 20.*20;
    
  };
    

  class X3872 : public  Chic1 {
    
    

  public :
    
    X3872() : Chic1(M_X3872) {};
    
    //interface getter
    amplitude Get() override {return Amp_Blend();}
   
    amplitude Amp_High(){
      
   
      // X(3872) High Energy
      amplitude X_omegaH = new_amplitude<regge::vector_exchange>(Kine(), "#omega exchange");
      X_omegaH->set_parameters({gX_omega, gV_omega, gT_omega, LamOmega, inter, slope});

      amplitude X_rhoH = new_amplitude<regge::vector_exchange>(Kine(), "#rho exchange");
      X_rhoH->set_parameters({gX_rho, gV_rho, gT_rho, LamRho, inter, slope});
    
      // Total is the sum of the above exchanges
      amplitude X_H = X_omegaH + X_rhoH;
      X_H->set_id("#it{X}(3872) High");

      return X_H; //give up ownership
    }
  
    amplitude Amp_Low(){
      amplitude X_omegaL = new_amplitude<analytic::vector_exchange>(Kine(), M_OMEGA, "#omega exchange");
      X_omegaL->set_parameters({gX_omega, gV_omega, gT_omega, LamOmega});

      amplitude X_rhoL = new_amplitude<analytic::vector_exchange>(Kine(), M_RHO, "#rho exchange");
      X_rhoL->set_parameters({gX_rho, gV_rho, gT_rho, LamRho});
    
      // Total is the sum of the above exchanges
      amplitude X_L = X_omegaL + X_rhoL;
      X_L->set_id("#it{X}(3872)");

      return X_L; //give up ownership
    }

    amplitude Amp_Blend(){
      return new_amplitude<blended>(Kine(),std::array<amplitude,2>{Amp_Low(),Amp_High()} ,std::array<double,2>{low_sLimit,high_sLimit} , "X blended");
       
    }
  
  
  private :
  
  };//X class



  class Chic1_3510 : public Chic1 {

  public:
   
    Chic1_3510() : Chic1(M_CHIC1) {}
   
    //interface getter
    amplitude Get() override {return Amp_Blend();}
   
    amplitude Amp_High(){
      // chi_c1
      amplitude ChiC1_omegaH = new_amplitude<regge::vector_exchange>(Kine(), "#omega exchange");
      ChiC1_omegaH->set_parameters({gChi_omega, gV_omega, gT_omega, LamOmega, inter, slope});

      amplitude ChiC1_rhoH = new_amplitude<regge::vector_exchange>(Kine(),  "#rho exchange");
      ChiC1_rhoH->set_parameters({gChi_rho, gV_rho, gT_rho, LamRho, inter, slope});

      amplitude ChiC1_H = ChiC1_omegaH + ChiC1_rhoH;
      ChiC1_H->set_id("#chi_{c1}");

      return ChiC1_H;
    
    }

    amplitude Amp_Low(){

      // chi_c1
      amplitude ChiC1_omegaL = new_amplitude<analytic::vector_exchange>(Kine(), M_OMEGA, "#omega exchange");
      ChiC1_omegaL->set_parameters({gChi_omega, gV_omega, gT_omega, LamOmega});

      amplitude ChiC1_rhoL = new_amplitude<analytic::vector_exchange>(Kine(), M_RHO, "#rho exchange");
      ChiC1_rhoL->set_parameters({gChi_rho, gV_rho, gT_rho, LamRho});

      amplitude ChiC1_phiL = new_amplitude<analytic::vector_exchange>(Kine(), M_PHI, "#phi exchange");
      ChiC1_phiL->set_option(analytic::vector_exchange::kNoFF);
      ChiC1_phiL->set_parameters({gChi_phi, gV_phi, gT_phi});

      amplitude ChiC1_psiL = new_amplitude<analytic::vector_exchange>(Kine(), M_JPSI, "#it{J}/#psi exchange");
      ChiC1_psiL->set_option(analytic::vector_exchange::kNoFF);

      ChiC1_psiL->set_parameters({gChi_psi, gV_psi, gT_psi});
    
      amplitude ChiC1_L = ChiC1_omegaL + ChiC1_rhoL + ChiC1_phiL + ChiC1_psiL;
      ChiC1_L->set_id("#chi_{c1}");

      return ChiC1_L;
    
    }

    amplitude Amp_Blend(){
      return new_amplitude<blended>(Kine(),std::array<amplitude,2>{Amp_Low(),Amp_High()} ,std::array<double,2>{low_sLimit,high_sLimit} , "chi_c1 blended");
       
    }
  
  private:
     
     
  };
   
}
 
