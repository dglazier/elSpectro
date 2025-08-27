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
#include "src/amplitude.hpp"

namespace elSpectro{

  using jpacAmp_ptr = std::shared_ptr<jpacPhoto::raw_amplitude>;

  class JpacTwoBody : public TwoBodyProduction {

  public:
    
    JpacTwoBody()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    JpacTwoBody(  jpacAmp_ptr amp ,const decaying_objs& decs, const particle_objs& stables);
    
    bool HasAngularDistribution() override{return true; } //I have an angular distribution


    double MatrixElementsSquared_T() const override {
      //set current value of mass
      // _amp->get_kinematics()->set_meson_mass( GetMeson()->PdgMass() );
       _amp->get_kinematics()->set_meson_mass( GetMeson()->Mass() );
       //std::cout<<"me "<<GetMeson()->Mass()<<" t "<<get_t()<<" s "<<get_s()<<" W "<<get_W()<<" jpac "<<_amp->get_kinematics()->Wth()<<" VAL "<<_amp->probability_distribution(get_s(),get_t())/4<<" "<<_amp->differential_xsection(get_s(),get_t())<<std::endl;
      //std::cout<<"JpacTwoBody::MatrixElementsSquared_T "<< GetMeson()->Mass() <<" "<<get_s()<<" "<<get_t()<<" "<<get_cosThCM()<<" "<<get_W()<<std::endl;
      // if(-get_t()<_amp->get_kinematics()->t_min(get_s())) return 0.;
      //if(-get_t()<_amp->get_kinematics()->t_max(get_s())) return 0.;
      // std::cout<<"JpacTwoBody::MatrixElementsSquared_T CHECK "<< _amp->probability_distribution(28,-0.01)<<std::endl;
      //_amp->get_kinematics()->set_meson_mass( GetMeson()->Mass() );
      //if(get_W()<_amp->get_kinematics()->Wth()) return 0;
      // auto res = _amp->probability_distribution(get_s(),get_t())/4;
      //std::cout<<"JpacTwoBody::MatrixElementsSquared_T "<< GetMeson()->Mass() <<" "<<get_s()<<" "<<get_t()<<" "<<get_cosThCM()<<" "<<get_W()<<" "<<std::endl;
     //auto res = _amp->probability_distribution(get_s(),get_t())/4;
     // std::cout<<"JpacTwoBody::MatrixElementsSquared_T "<< GetMeson()->Mass() <<" "<<get_s()<<" "<<get_t()<<" "<<get_cosThCM()<<" "<<get_W()<<" "<<_amp->probability_distribution(get_s(),_amp->get_kinematics()->t_min(get_s()))/4<<std::endl;
      // if(get_W()>2.32) exit(0);

      return _amp->probability_distribution(get_s(),get_t())/4;// Average over initial state helicites;
    }
    double MatrixElementsSquared_T_at_tmin() const override {

      //std::cout<<"JpacTwoBody::MatrixElementsSquared_T "<< GetMeson()->Mass() <<" "<<get_s()<<" "<<get_t()<<" "<<get_cosThCM()<<" "<<get_W()<<std::endl;
      // if(-get_t()<_amp->get_kinematics()->t_min(get_s())) return 0.;
      //if(-get_t()<_amp->get_kinematics()->t_max(get_s())) return 0.;
      
       _amp->get_kinematics()->set_meson_mass( GetMeson()->PdgMass() );
      // _amp->get_kinematics()->set_meson_mass( GetMeson()->Mass() );
       //if(get_W()<_amp->get_kinematics()->Wth()) return 0;
        return _amp->probability_distribution(get_s(),_amp->get_kinematics()->t_min(get_s()))/4;// Average over initial state helicites;
    }

    double IntegratedCrossSection(double W){
      return _amp->integrated_xsection(W*W);
    }
  private:

    jpacAmp_ptr _amp={nullptr}; //shared_ptr

    ClassDefOverride(elSpectro::JpacTwoBody,1); //class JpacTwoBody
    
  };//class JpacTwoBody

}//namespace elSpectro
