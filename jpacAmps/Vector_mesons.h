///////////////////////////////////////////////////////////
// taken from $JPACPHOTO/scripts/vector_SDMEs/dsigdt.cpp

#pragma once

#include "JpacAmplitude.h" //interface

#include "constants.hpp"

#include "regge/reggeon_exchange.hpp"

namespace jpacAmps {

  // using namespace jpacPhoto;

  //J^PC = 1^+ 1^+-
  class Vector_mesons : public JpacAmplitude {

  public:

  Vector_mesons(double mass): JpacAmplitude(mass) {
      
      Kine()->set_meson_JP(VECTOR);
      rel[0] = {1., 0.,             0.};
      rel[1] = {1., -0.95/sqrt(2.), -0.56};
      rel[2] = {1., -0.83/sqrt(2.),  0.} ;  
      
    }
    
  protected:
    
    //-------------------------------------------------
    // Unnatural exchanges
    
    double bU = 0.; // no form-factor 
    
    // Top couplings
    // Multiplied by masses and half to match normalization
    std::array<double,3> gT_pi  = {0.252*M_RHO*M_RHO/2, 0.696*M_OMEGA*M_OMEGA/2, 0.040*M_PHI*M_PHI/2};
    std::array<double,3> gT_eta = {0.136*M_RHO*M_RHO/2, 0.479*M_OMEGA*M_OMEGA/2, 0.210*M_PHI*M_PHI/2};

    // Bottom couplings
    double g_piNN = 13.26,  g_etaNN = 2.24;
    std::array<double, 2> betaB_pi  = {0., -2.*M_PROTON*g_piNN};  // Only flip allowed
    std::array<double, 2> betaB_eta = {0., -2.*M_PROTON*g_etaNN}; // Only flip allowed

    // Trajectory
    double alpha0U = -0.013, alphaPU = 0.7;
    //-------------------------------------------------
    // Natural exchanges

    // Cutoffs
    double b_pom     = 3.60,  b_f2     = 0.55,  b_a2     = 0.53;

    // Trajectories
    double alp0_pom  = 1.08,  alp0_f2  = 0.5,   alp0_a2  = 0.5;
    double alpP_pom  = 0.2,   alpP_f2  = 0.9,   alpP_a2  = 0.9;

    // Overall normalization
    std::array<double,3> beta_pom = {2.506*2/alpP_pom/PI, 0.739*2/alpP_pom/PI, 0.932*2/alpP_pom/PI};
    std::array<double,3> beta_f2  = {2.476*2/alpP_f2/PI,  0.730*2/alpP_f2/PI,  0.}; 
    std::array<double,3> beta_a2  = {0.370*2/alpP_a2/PI,  1.256*2/alpP_a2/PI,  0.};

    // Relative couplings
    std::array<std::array<double,3>,3> rel;

    // Bottom couplings 
    double kappa_pom = 0.,    kappa_f2 = 0,     kappa_a2 = -8.0;

  };//Vector_mesons

  /////////////////////////////////////////////////////////
  //rho vector meson
  
  class rho : public Vector_mesons {
    
  public:
    
    rho() : Vector_mesons(M_RHO) {}
    
    //interface getter
    amplitude Get() override {return Amp_Total();}
    
    amplitude Amp_Total(){
      amplitude rho_total = Amp_Pomeron() + Amp_f2() + Amp_a2() + Amp_eta() + Amp_pi();
      return rho_total;
    }
    amplitude Amp_Pomeron(){
      amplitude rho_pomeron = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "Pomeron");
      rho_pomeron->set_parameters( {alp0_pom, alpP_pom, b_pom, 
	  beta_pom[index]*rel[0][0], beta_pom[index]*rel[0][1], beta_pom[index]*rel[0][2], 
	  1., kappa_pom} );
      
      return rho_pomeron;
    }
    
    amplitude Amp_f2(){
     amplitude rho_f2 = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "#it{f}_{2}");
     
     rho_f2->set_parameters( { alp0_f2, alpP_f2, b_f2, 
	 beta_f2[index]*rel[1][0], beta_f2[index]*rel[1][1], beta_f2[index]*rel[1][2], 
	 1., kappa_f2} );

     return rho_f2;  
   }
   
   amplitude Amp_a2(){
     amplitude rho_a2 = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "#it{a}_{2}");

     rho_a2->set_parameters( { alp0_a2, alpP_a2, b_a2, 
	 beta_a2[index]*rel[2][0], beta_a2[index]*rel[2][1], beta_a2[index]*rel[2][2], 
	 1., kappa_a2} );

     return rho_a2;
   }
   amplitude Amp_pi(){
     amplitude rho_pi  = new_amplitude<regge::reggeon_exchange>(Kine(), -1, "#pi");
     
     rho_pi->set_parameters( {alpha0U, alphaPU, bU, 
	 gT_pi[index],   sqrt(2)*gT_pi[index],  gT_pi[index],  
	 betaB_pi[0],  betaB_pi[1]});

     return rho_pi;
   }
   amplitude Amp_eta(){
     amplitude rho_eta = new_amplitude<regge::reggeon_exchange>(Kine(), -1, "#eta");

     rho_eta->set_parameters({alpha0U, alphaPU, bU, 
	 gT_eta[index],   sqrt(2)*gT_eta[index],  gT_eta[index],
	 betaB_eta[0], betaB_eta[1]});

     return rho_eta;
   }

  private:
    static constexpr uint index=0;
  };
 class omega : public Vector_mesons {
    
  public:
    
  omega() : Vector_mesons(M_OMEGA) {}
    
    //interface getter
    amplitude Get() override {return Amp_Total();}
    
    amplitude Amp_Total(){
      amplitude rho_total = Amp_Pomeron() + Amp_f2() + Amp_a2() + Amp_eta() + Amp_pi();
      return rho_total;
    }
    amplitude Amp_Pomeron(){
      amplitude rho_pomeron = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "Pomeron");
      rho_pomeron->set_parameters( {alp0_pom, alpP_pom, b_pom, 
	  beta_pom[index]*rel[0][0], beta_pom[index]*rel[0][1], beta_pom[index]*rel[0][2], 
	  1., kappa_pom} );
      
      return rho_pomeron;
    }
    
    amplitude Amp_f2(){
     amplitude rho_f2 = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "#it{f}_{2}");
     
     rho_f2->set_parameters( { alp0_f2, alpP_f2, b_f2, 
	 beta_f2[index]*rel[1][0], beta_f2[index]*rel[1][1], beta_f2[index]*rel[1][2], 
	 1., kappa_f2} );

     return rho_f2;  
   }
   
   amplitude Amp_a2(){
     amplitude rho_a2 = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "#it{a}_{2}");

     rho_a2->set_parameters( { alp0_a2, alpP_a2, b_a2, 
	 beta_a2[index]*rel[2][0], beta_a2[index]*rel[2][1], beta_a2[index]*rel[2][2], 
	 1., kappa_a2} );

     return rho_a2;
   }
   amplitude Amp_pi(){
     amplitude rho_pi  = new_amplitude<regge::reggeon_exchange>(Kine(), -1, "#pi");
     
     rho_pi->set_parameters( {alpha0U, alphaPU, bU, 
	 gT_pi[index],   sqrt(2)*gT_pi[index],  gT_pi[index],  
	 betaB_pi[0],  betaB_pi[1]});

     return rho_pi;
   }
   amplitude Amp_eta(){
     amplitude rho_eta = new_amplitude<regge::reggeon_exchange>(Kine(), -1, "#eta");

     rho_eta->set_parameters({alpha0U, alphaPU, bU, 
	 gT_eta[index],   sqrt(2)*gT_eta[index],  gT_eta[index],
	 betaB_eta[0], betaB_eta[1]});

     return rho_eta;
   }

  private:
    static constexpr uint index=1;
  };

  ///////////////////////////////////////////////////////
  //phi
  class phi : public Vector_mesons {
    
  public:
    
  phi() : Vector_mesons(M_PHI) {}
    
    //interface getter
    amplitude Get() override {return Amp_Total();}
    
    amplitude Amp_Total(){
      amplitude rho_total = Amp_Pomeron() + Amp_f2() + Amp_a2() + Amp_eta() + Amp_pi();
      return rho_total;
    }
    amplitude Amp_Pomeron(){
      amplitude rho_pomeron = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "Pomeron");
      rho_pomeron->set_parameters( {alp0_pom, alpP_pom, b_pom, 
	  beta_pom[index]*rel[0][0], beta_pom[index]*rel[0][1], beta_pom[index]*rel[0][2], 
	  1., kappa_pom} );
      
      return rho_pomeron;
    }
    
    amplitude Amp_f2(){
     amplitude rho_f2 = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "#it{f}_{2}");
     
     rho_f2->set_parameters( { alp0_f2, alpP_f2, b_f2, 
	 beta_f2[index]*rel[1][0], beta_f2[index]*rel[1][1], beta_f2[index]*rel[1][2], 
	 1., kappa_f2} );

     return rho_f2;  
   }
   
   amplitude Amp_a2(){
     amplitude rho_a2 = new_amplitude<regge::reggeon_exchange>(Kine(), +1, "#it{a}_{2}");

     rho_a2->set_parameters( { alp0_a2, alpP_a2, b_a2, 
	 beta_a2[index]*rel[2][0], beta_a2[index]*rel[2][1], beta_a2[index]*rel[2][2], 
	 1., kappa_a2} );

     return rho_a2;
   }
   amplitude Amp_pi(){
     amplitude rho_pi  = new_amplitude<regge::reggeon_exchange>(Kine(), -1, "#pi");
     
     rho_pi->set_parameters( {alpha0U, alphaPU, bU, 
	 gT_pi[index],   sqrt(2)*gT_pi[index],  gT_pi[index],  
	 betaB_pi[0],  betaB_pi[1]});

     return rho_pi;
   }
   amplitude Amp_eta(){
     amplitude rho_eta = new_amplitude<regge::reggeon_exchange>(Kine(), -1, "#eta");

     rho_eta->set_parameters({alpha0U, alphaPU, bU, 
	 gT_eta[index],   sqrt(2)*gT_eta[index],  gT_eta[index],
	 betaB_eta[0], betaB_eta[1]});

     return rho_eta;
   }

  private:
    static constexpr uint index=2;
  };
}
