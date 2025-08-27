#pragma once

#include "JpacAmplitude.h" //interface

#include "constants.hpp"
#include "covariant/pomeron_exchange.hpp"
#include "analytic/pomeron_exchange.hpp"
#include "analytic/baryon_BW.hpp"


namespace jpacAmps {

  // using namespace jpacPhoto;

  //J^PC = 1^+ 1^--
  class Jpsi_P : public JpacAmplitude {

  public:

    Jpsi_P(double mass=3.0969000): JpacAmplitude(mass) {
      Kine()->set_meson_JP(VECTOR);
    }

    amplitude Get()  override {return Pomeron_Exchange();}
    amplitude GetPentaquarks()  {
      return Pomeron_Exchange() + P_c4312() + P_c4380() + P_c4457();
    }

    //----------------------------------------------------------------
    //S - CHANNEL
    amplitude P_c4312() const {
      // 2% BR coupling to JPSI by VMD with equal photocouplings

      auto P_c4312 = new_amplitude<baryon_BW>(Kine(), 1, std::array<double,2>{4.312, 0.0098}, "P_{c}(4312)");
      P_c4312->set_parameters({0.02, F_JPSI, 0.7071});
      return P_c4312;
    }
   amplitude P_c4380() const {
     auto P_c4440 = new_amplitude<baryon_BW>(Kine(),  -3, std::array<double,2>{4.38, 0.021}, "P_{c}(4440)");
     //      P_c4440->set_parameters({0.01, F_JPSI, 0.7071});
      P_c4440->set_parameters({0.02, F_JPSI, 0.7071});
      return P_c4440;
    }
    amplitude P_c4457() const {
      auto P_c4457 = new_amplitude<baryon_BW>(Kine(), 5, std::array<double,2>{4.457, 0.0064}, "P_{c}(4457)" );
      //      P_c4457->set_parameters({0.005, F_JPSI, 0.7071});
      P_c4457->set_parameters({0.02, F_JPSI, 0.7071});
      return P_c4457;
    }
    // -----------------------------------------------------------------
    // T - CHANNEL
    
    amplitude Pomeron_Exchange()  {
      // Set up pomeron trajectory
      //      auto alpha =new linear_trajectory(+1, 0.941, 0.364, "pomeron");
      // Create amplitude with kinematics and trajectory
      auto  background = new_amplitude<covariant::pomeron_exchange>(Kine(), "Background");
      // normalization, t-slope and linear trajectory
      const double E        = sqrt(4. * PI * ALPHA);
      background->set_parameters({E*0.379, 0.12, 0.941, 0.364});
      return background;
    }
     
  protected:
  

  };
}

