#include "TwoBody_stu.h"
#include "FunctionsForKinematics.h"
#include <TRandom.h>
#include <TMath.h>

namespace elSpectro{

  TwoBody_stu::TwoBody_stu(double s,double t,double t_slope,double u,double u_slope):
    _s_strength{s},
    _t_strength{t},_t_slope{t_slope},
    _u_strength{u},_u_slope{u_slope}
  {
    //renormalise total strength =1
    double total_strength=s+t+u;
    _s_strength/=total_strength;
    _t_strength/=total_strength;
    _u_strength/=total_strength;
    if(_u_strength) {
      std::cerr<<"TwoBody_stu::TwoBody_stu u channel not implemented yet..."<<std::endl;
      exit(0);
    }
    
  }

  /////////////////////////////////////////////////////////
  void TwoBody_stu::PostInit(ReactionInfo* info) {
    auto reactionInfo = static_cast<ReactionElectroProd*> (info);
    _p1=reactionInfo->_photon;//gamma
    _p2=reactionInfo->_target;
    _p3=reactionInfo->_meson;
    _p4=reactionInfo->_baryon;
    _CM=reactionInfo->_photoN;
  }
  
}
