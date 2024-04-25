//////////////////////////////////////////////////////////////
///
///Class:		JpacModelst
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) preconfigured jpacPhoto amplitude
///             2) it decay as a function of s and t JpacDecayst
///
///            Note derived classes should include a constructor to initialise
///            JpacModelst( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModelst.h"
#include "SDME.h"
#include "FunctionsForElectronScattering.h"
#include "core/amplitude.hpp"

namespace elSpectro{

  using jpacAmp_ptr = jpacPhoto::amplitude*;

  class JpacModelst : public DecayModelst {

  public:
    
    JpacModelst()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    JpacModelst( jpacAmp_ptr amp, particle_ptrs parts,
		  const std::vector<int> pdgs  );
    
    bool HasAngularDistribution() override{return true; } //I have an angular distribution


    double MatrixElementsSquared_T() const override {
      _amp->_kinematics->set_meson_mass( GetMeson()->Mass() );
     // _amp->_kinematics->set_Q2( get_Q2() );
      //std::cout<<"me "<<GetMeson()->Mass()<<" Q2 "<< get_Q2()<<" t "<<get_t()<<" s "<<get_s()<<" W "<<get_W()<<" jpac "<<_amp->_kinematics->Wth()<<" VAL "<<_amp->probability_distribution(get_s(),get_t())/4<<std::endl;
      if(get_W()<_amp->_kinematics->Wth()) return 0;
      return _amp->probability_distribution(get_s(),get_t())/4;// Average over initial state helicites;
    }
    
    void CalcMesonSDMEs() const  override ;
    void CalcBaryonSDMEs() const   override ;
    
  private:

    jpacAmp_ptr _amp={nullptr}; //I am not the owner

    ClassDefOverride(elSpectro::JpacModelst,1); //class JpacModelst
    
  };//class JpacModelst

}//namespace elSpectro
