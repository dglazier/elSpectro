//////////////////////////////////////////////////////////////
///
///Class:		JpacTwoBody
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) preconfigured jpacPhoto amplitude
///             2) it decay as a function of s and t JpacDecayst
///
///            Note derived classes should include a constructor to initialise
///            JpacTwoBody( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "TwoBodyProduction.h"
#include "SDME.h"
#include "FunctionsForElectronScattering.h"
#include "core/amplitude.hpp"

namespace elSpectro{

  using jpacAmp_ptr = jpacPhoto::amplitude*;

  class JpacTwoBody : public TwoBodyProduction {

  public:
    
    JpacTwoBody()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    JpacTwoBody( jpacAmp_ptr amp, particle_ptrs parts,
		  const std::vector<int> pdgs  );
    
    bool HasAngularDistribution() override{return true; } //I have an angular distribution


    double MatrixElementsSquared_T() const override {

      //   std::cout<<"JpacTwoBody::MatrixElementsSquared_T "<< GetMeson() <<" "<<get_s()<<" "<<get_t()<<std::endl;
      _amp->_kinematics->set_meson_mass( GetMeson()->Mass() );
      if(get_W()<_amp->_kinematics->Wth()) return 0;
      return _amp->probability_distribution(get_s(),get_t())/4;// Average over initial state helicites;
    }
    
  private:

    jpacAmp_ptr _amp={nullptr}; //I am not the owner

    ClassDefOverride(elSpectro::JpacTwoBody,1); //class JpacTwoBody
    
  };//class JpacTwoBody

}//namespace elSpectro
