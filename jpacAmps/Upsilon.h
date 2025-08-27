#pragma once

#include "JpacAmplitude.h" //interface

#include "constants.hpp"
#include "covariant/pomeron_exchange.hpp"
#include "analytic/pomeron_exchange.hpp"
//#include "baryon_resonance.hpp"


namespace jpacAmps {

  // using namespace jpacPhoto;

  //J^PC = 1^+ 1^--
  class Jpsi_P : public JpacAmplitude {

  public:

    Jpsi_P(double mass=3.0969000): JpacAmplitude(mass) {
      Kine()->set_meson_JP(VECTOR);
    }

    amplitude Get()  override {return Pomeron_Exchange();}

    // ----------------------------------------------------------------
    // S - CHANNEL
    /* amplitude P_c4450() const { */
    /*   auto P_c4450 = new_amplitude<baryon_resonance>(Kine(), +1,  3, -1, 4.45, 0.040, "P_{c}(4450)"); */
    /*   P_c4450->set_parameters( {0.01, .7071} ); */
    /*   return P_c4450; */
    /* } */
    /* amplitude P_c4380() const { */
    /*   auto P_c4380 = new_amplitude<baryon_resonance>(Kine(),5, +1, 4.38, 0.205, "P_{c}(4380)" ); */
    /*   P_c4380->set_parameters( {0.01, .7071} ); */
    /*   return P_c4380; */
    /* } */
    // -----------------------------------------------------------------
    // T - CHANNEL
    
    amplitude Pomeron_Exchange()  {
      // Set up pomeron trajectory
      //      auto alpha =new linear_trajectory(+1, 0.941, 0.364, "pomeron");
      // Create amplitude with kinematics and trajectory
      auto  background = new_amplitude<analytic::pomeron_exchange>(Kine(), "Background");
      // normalization, t-slope and linear trajectory
      background->set_parameters({0.379, 0.12, 0.941, 0.364,});
      return background;
    }
    
  protected:
  

  };
}

