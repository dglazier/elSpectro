//////////////////////////////////////////////////////////////
///
///Class:		JpacModelstM2
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) preconfigured jpacPhoto amplitude
///             2) it decay as a function of s and t JpacDecayst
///
///            Note derived classes should include a constructor to initialise
///            JpacModelstM2( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModelstM2.h"
#include "FunctionsForElectronScattering.h"
#include "inclusive/inclusive_production.hpp"

namespace elSpectro{

 
  class JpacModelstM2 : public DecayModelstM2 {

  public:
    
    JpacModelstM2()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    JpacModelstM2( jpacPhoto::inclusive_production* inc, particle_ptrs parts,
		  const std::vector<int> pdgs  );
    
    bool HasAngularDistribution() override{return true; } //I have an angular distribution
    // bool HasMassDistribution() override{return true; } //I have an angular distribution


    double d3sigma() const override {
      _inc->_kinematics->set_meson_mass( GetMeson()->Mass() );
     // _amp->_kinematics->set_Q2( get_Q2() );
      //std::cout<<"me "<<GetMeson()->Mass()<<" Q2 "<< get_Q2()<<" t "<<get_t()<<" s "<<get_s()<<" W "<<get_W()<<" jpac "<<_amp->_kinematics->Wth()<<" VAL "<<_amp->probability_distribution(get_s(),get_t())/4<<std::endl;
      if(get_W()<_amp->_kinematics->Wth()) return 0;
      return _inc->dsigma_dtdM2(get_s(),get_t(),get_M2());
    }
    
    
  private:

    jpacPhoto::inclusive_production* _inc={nullptr}; //I am not the owner

    ClassDefOverride(elSpectro::JpacModelstM2,1); //class JpacModelstM2
    
  };//class JpacModelstM2

}//namespace elSpectro
