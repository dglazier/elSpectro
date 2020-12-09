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
#include "amplitudes/amplitude.hpp"

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

    //double Intensity() const override;

    //void PostInit(ReactionInfo* info) override;

    double MatrixElementsSquared_T() const override {
     _amp->kinematics->set_mX( GetMeson()->Mass() );


     // std::cout<<" MatrixElementsSquared_T s "<<get_s()<<" t "<<get_t()<<" mm "<<GetMeson()->Mass()<<" "<<_amp->probability_distribution(get_s(),get_t())<<std::endl;
      //return -get_t()/10;// + get_s()/10;
     // return get_s()/10;
      return _amp->probability_distribution(get_s(),get_t())/4;// Average over initial state helicites;
    }
    
    void CalcMesonSDMEs() const  override ;
    void CalcBaryonSDMEs() const   override ;
    
  private:

    jpacAmp_ptr _amp={nullptr}; //I am not the owner

    ClassDefOverride(elSpectro::JpacModelst,1); //class JpacModelst
    
  };//class JpacModelst

}//namespace elSpectro
